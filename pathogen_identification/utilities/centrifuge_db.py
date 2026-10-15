"""
Build a Centrifuge index from INSaFLU reference sources.

Classes:
    CentrifugeDBExporter  - select references from the DB and write input FASTA + conversion table
    CentrifugeDBTaxonomy  - check NCBI nodes.dmp / names.dmp and taxids used by the export
    CentrifugeDBNaming    - index base names (centrifuge_<YYYYMMDD>) and their directories
    CentrifugeDBBuilder   - run centrifuge-build, verify the index, promote tmp -> final
"""

import gzip
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import date
from typing import Callable, Dict, Iterator, List, Optional, Set, Tuple


class CentrifugeDBError(Exception):
    """Raised when the database cannot be built safely."""


def _televir_metadata_constants():
    # imported lazily: constants pull in the full settings / software config
    from constants.constants import Televir_Metadata_Constants

    return Televir_Metadata_Constants


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------


@dataclass
class ExportPlan:
    """Result of querying the DB. No sequence data is read yet."""

    accid_to_taxid: Dict[str, str] = field(default_factory=dict)
    accid_to_file: Dict[str, str] = field(default_factory=dict)
    conflicts: Dict[str, List[str]] = field(default_factory=dict)
    missing_files: List[str] = field(default_factory=list)
    no_taxid: List[str] = field(default_factory=list)

    @property
    def taxids(self) -> Set[str]:
        return set(self.accid_to_taxid.values())


def iter_fasta_records(path: str, wanted: Set[str]) -> Iterator[Tuple[str, str]]:
    """
    Stream (name, sequence) for records whose name is in `wanted`.
    name = first whitespace-delimited token of the header.
    Supports plain and gzipped FASTA.
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        name = None
        keep = False
        chunks: List[str] = []
        for line in handle:
            if line.startswith(">"):
                if keep:
                    yield name, "".join(chunks).upper()
                tokens = line[1:].split()
                name = tokens[0] if tokens else ""
                keep = name in wanted
                chunks = []
            elif keep:
                chunks.append(line.strip())
        if keep:
            yield name, "".join(chunks).upper()


class CentrifugeDBExporter:
    """
    Selects references visible to a user and writes the centrifuge inputs.

    Visible = system files (owner is null) + files owned by `user`.
    Deleted files are skipped. Only ReferenceSourceFileMap.STATUS_GOOD rows are used.
    Headers written to the FASTA are the accession (`ReferenceSource.accid`),
    which is also the first column of the conversion table.
    """

    def __init__(self, user, include_phage: bool = False):
        self.user = user
        self.include_phage = include_phage

    def get_maps(self):
        from django.db.models import Q

        from pathogen_identification.models import ReferenceSourceFileMap

        maps = ReferenceSourceFileMap.objects.filter(
            Q(reference_source_file__owner__isnull=True)
            | Q(reference_source_file__owner=self.user),
            reference_source_file__is_deleted=False,
            status=ReferenceSourceFileMap.STATUS_GOOD,
            reference_source__isnull=False,
            reference_source_file__isnull=False,
        ).select_related(
            "reference_source",
            "reference_source__taxid",
            "reference_source_file",
        )

        if not self.include_phage:
            maps = maps.exclude(reference_source__description__icontains="phage")

        return maps

    @staticmethod
    def normalise_taxid(taxid) -> Optional[str]:
        if taxid is None:
            return None
        taxid = str(taxid).strip()
        if not taxid.isdigit():
            return None
        return str(int(taxid))

    def build_plan(self) -> ExportPlan:
        plan = ExportPlan()
        accid_taxids: Dict[str, Set[str]] = {}

        for ref_map in self.get_maps().order_by("reference_source__accid"):
            accid = (ref_map.reference_source.accid or "").strip()
            if not accid:
                continue

            taxid_obj = ref_map.reference_source.taxid
            taxid = self.normalise_taxid(taxid_obj.taxid if taxid_obj else None)
            if taxid is None:
                plan.no_taxid.append(accid)
                continue

            accid_taxids.setdefault(accid, set()).add(taxid)

            if accid not in plan.accid_to_file:
                plan.accid_to_file[accid] = ref_map.reference_source_file.filepath

        for accid, taxids in accid_taxids.items():
            if len(taxids) > 1:
                plan.conflicts[accid] = sorted(taxids)
            else:
                plan.accid_to_taxid[accid] = next(iter(taxids))

        plan.no_taxid = sorted(set(plan.no_taxid))
        plan.missing_files = sorted(
            {
                path
                for accid, path in plan.accid_to_file.items()
                if accid in plan.accid_to_taxid and not os.path.isfile(path)
            }
        )
        return plan

    def read_sequences(self, plan: ExportPlan) -> Tuple[Dict[str, str], List[str]]:
        """
        Read sequences for every accid in the plan.
        Returns (accid -> sequence, accids not found in their source file).
        """
        by_file: Dict[str, Set[str]] = {}
        for accid, path in plan.accid_to_file.items():
            if accid in plan.accid_to_taxid:
                by_file.setdefault(path, set()).add(accid)

        sequences: Dict[str, str] = {}
        for path, wanted in sorted(by_file.items()):
            if not os.path.isfile(path):
                continue
            for name, seq in iter_fasta_records(path, wanted):
                if name not in sequences and seq:
                    sequences[name] = seq

        not_found = sorted(set(plan.accid_to_taxid) - set(sequences))
        return sequences, not_found

    @staticmethod
    def write(
        plan: ExportPlan,
        sequences: Dict[str, str],
        input_fna: str,
        conversion_tsv: str,
    ) -> int:
        """
        Write input FASTA and conversion table for accids that have a sequence.
        Returns number of records written.
        """
        accids = sorted(a for a in sequences if a in plan.accid_to_taxid)

        with open(input_fna, "w") as fna, open(conversion_tsv, "w") as conv:
            for accid in accids:
                fna.write(f">{accid}\n")
                seq = sequences[accid]
                for i in range(0, len(seq), 80):
                    fna.write(seq[i : i + 80] + "\n")
                conv.write(f"{accid}\t{plan.accid_to_taxid[accid]}\n")

        return len(accids)


# ---------------------------------------------------------------------------
# Taxonomy
# ---------------------------------------------------------------------------


class CentrifugeDBTaxonomy:
    """
    NCBI taxonomy dump checks. Expects nodes.dmp and names.dmp in `taxonomy_dir`.
    Field separator in these files is '\\t|\\t'; splitting on '|' and stripping works for both.
    """

    NODES = "nodes.dmp"
    NAMES = "names.dmp"

    def __init__(self, taxonomy_dir: str):
        self.taxonomy_dir = taxonomy_dir
        self.nodes_path = os.path.join(taxonomy_dir, self.NODES)
        self.names_path = os.path.join(taxonomy_dir, self.NAMES)
        self._node_ids: Optional[Set[str]] = None
        self._scientific_ids: Optional[Set[str]] = None

    def check_files(self):
        missing = [p for p in (self.nodes_path, self.names_path) if not os.path.isfile(p)]
        if missing:
            raise CentrifugeDBError(
                "Taxonomy files not found: {}. Put NCBI taxdump nodes.dmp and names.dmp "
                "there or pass --taxonomy-dir.".format(", ".join(missing))
            )

    @staticmethod
    def _fields(line: str) -> List[str]:
        return [f.strip() for f in line.rstrip("\n").rstrip("|").split("|")]

    def node_ids(self) -> Set[str]:
        if self._node_ids is None:
            self.check_files()
            with open(self.nodes_path) as handle:
                self._node_ids = {
                    self._fields(line)[0] for line in handle if line.strip()
                }
        return self._node_ids

    def scientific_name_ids(self) -> Set[str]:
        if self._scientific_ids is None:
            self.check_files()
            ids = set()
            with open(self.names_path) as handle:
                for line in handle:
                    fields = self._fields(line)
                    if len(fields) >= 4 and fields[3] == "scientific name":
                        ids.add(fields[0])
            self._scientific_ids = ids
        return self._scientific_ids

    def missing_taxids(self, taxids: Set[str]) -> List[str]:
        """Taxids absent from nodes.dmp or lacking a scientific name."""
        nodes = self.node_ids()
        names = self.scientific_name_ids()
        return sorted(t for t in taxids if t not in nodes or t not in names)


# ---------------------------------------------------------------------------
# Naming
# ---------------------------------------------------------------------------


class CentrifugeDBNaming:
    """
    Index names: <index_dir>/centrifuge_<YYYYMMDD>.{1,2,3}.cf
    Temporary build base: centrifuge_<YYYYMMDD>.tmp
    """

    PREFIX = "centrifuge"
    SUBDIR = "centrifuge"

    def __init__(self, base_dir: Optional[str] = None, build_date: Optional[date] = None):
        if base_dir is None:
            db_dir = _televir_metadata_constants().SOURCE["DBDIR_MAIN"]
            base_dir = os.path.join(db_dir, self.SUBDIR)
        self.base_dir = base_dir
        self.build_date = build_date or date.today()

    @property
    def name(self) -> str:
        return f"{self.PREFIX}_{self.build_date:%Y%m%d}"

    @property
    def final_base(self) -> str:
        return os.path.join(self.base_dir, self.name)

    @property
    def tmp_base(self) -> str:
        return os.path.join(self.base_dir, self.name + ".tmp")

    @property
    def work_dir(self) -> str:
        return os.path.join(self.base_dir, "work_" + self.name)

    @staticmethod
    def index_files(base: str) -> List[str]:
        """
        Index files of exactly this base: <base>.<digits>.cf.
        Must not match other bases such as <base>.tmp.1.cf.
        """
        pattern = re.compile(re.escape(os.path.basename(base)) + r"\.\d+\.cf$")
        directory = os.path.dirname(base) or "."
        return sorted(
            os.path.join(directory, name)
            for name in os.listdir(directory)
            if pattern.match(name)
        ) if os.path.isdir(directory) else []

    def will_replace(self) -> bool:
        return len(self.index_files(self.final_base)) > 0


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------


class CentrifugeDBRegistry:
    """
    Registers a built index in the Televir utility registry (`database` table),
    the same table the pipeline reads through UtilityDB.
    """

    CATEGORY = "centrifuge"
    DB_TYPE = "classification"
    SOFTWARE = "centrifuge"

    def __init__(self, registry_dir: Optional[str] = None):
        if registry_dir is None:
            from constants.televir_directories import \
                Televir_Directory_Constants

            registry_dir = Televir_Directory_Constants.docker_app_directory
        self.registry_dir = registry_dir

    @classmethod
    def build_record(cls, naming: "CentrifugeDBNaming", n_accids: int, n_taxids: int) -> dict:
        return {
            "name": naming.name,
            "db_category": cls.CATEGORY,
            "db_name": naming.name,
            "db_type": cls.DB_TYPE,
            "path": naming.final_base,
            "installed": "True",  # stored as text, matching existing rows
            "software": cls.SOFTWARE,
            "date": naming.build_date.isoformat(),
            "version": None,
            "source_url": None,
            "file_mod_date": None,
            "description": f"{n_accids} accessions, {n_taxids} taxids (NCBI taxonomy)",
        }

    def register(self, record: dict) -> None:
        from pathogen_identification.utilities.utilities_televir_dbs import \
            Utility_Repository

        repository = Utility_Repository(db_path=self.registry_dir, install_type="docker")
        repository.upsert_database(record)


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------


class CentrifugeDBBuilder:
    """
    Runs centrifuge-build and centrifuge-inspect. Commands are passed as lists (no shell).
    `runner` is injectable for tests; it must behave like subprocess.run.
    """

    def __init__(self, threads: int = 1, runner: Callable = subprocess.run):
        self.threads = threads
        self.runner = runner

    @staticmethod
    def resolve_binary(name: str) -> str:
        configured = _televir_metadata_constants()().get_software_binary(name)
        if os.path.isfile(configured):
            return configured
        found = shutil.which(name)
        if found:
            return found
        raise CentrifugeDBError(f"Binary '{name}' not found (looked in {configured} and PATH).")

    def build_command(
        self,
        input_fna: str,
        conversion_tsv: str,
        taxonomy: CentrifugeDBTaxonomy,
        tmp_base: str,
    ) -> List[str]:
        return [
            self.resolve_binary("centrifuge-build"),
            "-p",
            str(self.threads),
            "--conversion-table",
            conversion_tsv,
            "--taxonomy-tree",
            taxonomy.nodes_path,
            "--name-table",
            taxonomy.names_path,
            input_fna,
            tmp_base,
        ]

    def build(
        self,
        input_fna: str,
        conversion_tsv: str,
        taxonomy: CentrifugeDBTaxonomy,
        tmp_base: str,
        log_path: str,
    ) -> None:
        cmd = self.build_command(input_fna, conversion_tsv, taxonomy, tmp_base)
        with open(log_path, "w") as log:
            result = self.runner(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)
        if result.returncode != 0:
            raise CentrifugeDBError(
                f"centrifuge-build failed with exit code {result.returncode}. See {log_path}"
            )

    def count_index_sequences(self, index_base: str) -> int:
        cmd = [self.resolve_binary("centrifuge-inspect"), "-n", index_base]
        result = self.runner(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise CentrifugeDBError(f"centrifuge-inspect failed on {index_base}")
        return len([line for line in result.stdout.splitlines() if line.strip()])

    def verify(self, tmp_base: str, expected: int) -> None:
        if not CentrifugeDBNaming.index_files(tmp_base):
            raise CentrifugeDBError(f"No index files written for {tmp_base}")
        found = self.count_index_sequences(tmp_base)
        if found != expected:
            raise CentrifugeDBError(
                f"Index sequence count {found} differs from exported count {expected}"
            )

    @staticmethod
    def promote(tmp_base: str, final_base: str) -> List[str]:
        """Rename tmp index files to the final base. Replaces existing files of the same base."""
        files = CentrifugeDBNaming.index_files(tmp_base)
        if not files:
            raise CentrifugeDBError(f"Nothing to promote for {tmp_base}")

        for old_path in CentrifugeDBNaming.index_files(final_base):
            os.remove(old_path)

        new_paths = []
        for old_path in files:
            suffix = old_path[len(tmp_base):]
            new_path = final_base + suffix
            os.replace(old_path, new_path)
            new_paths.append(new_path)
        return new_paths

    @staticmethod
    def remove_partial(tmp_base: str) -> None:
        for path in CentrifugeDBNaming.index_files(tmp_base):
            os.remove(path)
