"""
Build a Kraken 2 database from INSaFLU reference sources.

Reuses from centrifuge_db: CentrifugeDBExporter (input.fna + conversion table),
CentrifugeDBTaxonomy (nodes/names checks), CentrifugeDBRegistry (registration).

Kraken 2 database = a directory:
    <tmp>/library/input.fna        sequences, headers = accession
    <tmp>/seqid2taxid.map          accession<TAB>taxid (read by k2 build)
    <tmp>/taxonomy/nodes.dmp       NCBI taxonomy
    <tmp>/taxonomy/names.dmp
    <tmp>/taxonomy/merged.dmp      optional, remaps retired taxids
Build output: hash.k2d, opts.k2d, taxo.k2d
"""

import os
import shutil
import subprocess
from datetime import date
from typing import Callable, List, Optional, Set

from pathogen_identification.utilities.centrifuge_db import (
    CentrifugeDBBuilder, CentrifugeDBError, CentrifugeDBExporter,
    CentrifugeDBRegistry, CentrifugeDBTaxonomy)


class Kraken2DBError(CentrifugeDBError):
    """Raised when a Kraken 2 database cannot be built safely."""


KRAKEN2_DEFAULT_KMER_LEN = 35
KRAKEN2_DEFAULT_MINIMIZER_LEN = 31
KRAKEN2_MAX_MINIMIZER_LEN_NUCLEOTIDE = 31
OUTPUT_FILES = ("hash.k2d", "opts.k2d", "taxo.k2d")


def validate_kmer_settings(kmer_len: int, minimizer_len: int) -> None:
    """Kraken 2 rules for nucleotide databases."""
    if kmer_len < 1:
        raise Kraken2DBError("--kmer-len must be >= 1")
    if minimizer_len < 1:
        raise Kraken2DBError("--minimizer-len must be >= 1")
    if minimizer_len > KRAKEN2_MAX_MINIMIZER_LEN_NUCLEOTIDE:
        raise Kraken2DBError(
            f"--minimizer-len must be <= {KRAKEN2_MAX_MINIMIZER_LEN_NUCLEOTIDE} "
            "for nucleotide databases"
        )
    if minimizer_len > kmer_len:
        raise Kraken2DBError("--minimizer-len must be <= --kmer-len")


class Kraken2DBNaming:
    """
    Directory names: <base>/kraken2_<YYYYMMDD>
    Temporary build directory: <base>/kraken2_<YYYYMMDD>.tmp
    """

    PREFIX = "kraken2"
    SUBDIR = "kraken2"

    def __init__(self, base_dir: Optional[str] = None, build_date: Optional[date] = None):
        if base_dir is None:
            from pathogen_identification.utilities.centrifuge_db import \
                _televir_metadata_constants

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

    def will_replace(self) -> bool:
        return os.path.isdir(self.final_base)


class Kraken2DBTaxonomy(CentrifugeDBTaxonomy):
    """
    Adds merged.dmp handling. Taxids retired in NCBI are remapped by k2 build
    when merged.dmp is present; taxids absent with no merge target are errors.
    """

    MERGED = "merged.dmp"

    def __init__(self, taxonomy_dir: str):
        super().__init__(taxonomy_dir)
        self.merged_path = os.path.join(taxonomy_dir, self.MERGED)

    def merged_map(self) -> dict:
        merged = {}
        if os.path.isfile(self.merged_path):
            with open(self.merged_path) as handle:
                for line in handle:
                    fields = self._fields(line)
                    if len(fields) >= 2 and fields[0]:
                        merged[fields[0]] = fields[1]
        return merged

    def classify_missing(self, taxids: Set[str]):
        """
        Returns (unresolvable, remappable).
        unresolvable: absent from nodes/names and not in merged.dmp -> error.
        remappable:   absent but merged.dmp maps it -> k2 build remaps it.
        """
        missing = self.missing_taxids(taxids)
        merged = self.merged_map()
        remappable = [t for t in missing if t in merged]
        unresolvable = [t for t in missing if t not in merged]
        return unresolvable, remappable


class Kraken2DBStaging:
    """Lays out the Kraken 2 input directory."""

    @staticmethod
    def stage(
        plan,
        sequences: dict,
        taxonomy: Kraken2DBTaxonomy,
        tmp_dir: str,
    ) -> int:
        library_dir = os.path.join(tmp_dir, "library")
        taxonomy_dir = os.path.join(tmp_dir, "taxonomy")
        os.makedirs(library_dir, exist_ok=True)
        os.makedirs(taxonomy_dir, exist_ok=True)

        written = CentrifugeDBExporter.write(
            plan,
            sequences,
            os.path.join(library_dir, "input.fna"),
            os.path.join(tmp_dir, "seqid2taxid.map"),
        )

        shutil.copyfile(taxonomy.nodes_path, os.path.join(taxonomy_dir, "nodes.dmp"))
        shutil.copyfile(taxonomy.names_path, os.path.join(taxonomy_dir, "names.dmp"))
        if os.path.isfile(taxonomy.merged_path):
            shutil.copyfile(taxonomy.merged_path, os.path.join(taxonomy_dir, "merged.dmp"))

        return written


class Kraken2DBBuilder:
    """
    Runs k2 build and k2 inspect, and checks dustmasker is on PATH
    (k2 build masks low-complexity sequence with it). `runner` is injectable.
    """

    BINARY = "k2"
    MASKER = "dustmasker"

    def __init__(self, threads: int = 1, runner: Callable = subprocess.run):
        self.threads = threads
        self.runner = runner

    @staticmethod
    def require_masker() -> str:
        found = shutil.which(Kraken2DBBuilder.MASKER)
        if not found:
            raise Kraken2DBError(
                "dustmasker not found on PATH (BLAST+). k2 build needs it for masking."
            )
        return found

    def resolve_binary(self) -> str:
        return CentrifugeDBBuilder.resolve_binary(self.BINARY)

    def build_command(
        self,
        db_dir: str,
        kmer_len: int,
        minimizer_len: int,
        max_db_size: Optional[str] = None,
    ) -> List[str]:
        cmd = [
            self.resolve_binary(),
            "build",
            "--db", db_dir,
            "--threads", str(self.threads),
            "--kmer-len", str(kmer_len),
            "--minimizer-len", str(minimizer_len),
        ]
        if max_db_size:
            cmd += ["--max-db-size", str(max_db_size)]
        return cmd

    def build(
        self,
        db_dir: str,
        kmer_len: int,
        minimizer_len: int,
        log_path: str,
        max_db_size: Optional[str] = None,
    ) -> None:
        cmd = self.build_command(db_dir, kmer_len, minimizer_len, max_db_size)
        with open(log_path, "w") as log:
            result = self.runner(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)
        if result.returncode != 0:
            raise Kraken2DBError(
                f"k2 build failed with exit code {result.returncode}. See {log_path}"
            )

    def verify(self, db_dir: str) -> None:
        missing = [f for f in OUTPUT_FILES if not os.path.isfile(os.path.join(db_dir, f))]
        if missing:
            raise Kraken2DBError(f"Missing database files: {', '.join(missing)}")

        cmd = [self.resolve_binary(), "inspect", "--db", db_dir, "--skip-counts"]
        result = self.runner(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise Kraken2DBError(f"k2 inspect failed on {db_dir}")

    @staticmethod
    def promote(tmp_dir: str, final_dir: str) -> None:
        """Replace final_dir (same-day build) with tmp_dir."""
        if os.path.isdir(final_dir):
            shutil.rmtree(final_dir)
        os.replace(tmp_dir, final_dir)

    @staticmethod
    def remove_partial(tmp_dir: str) -> None:
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)


class Kraken2DBRegistry(CentrifugeDBRegistry):
    """Same table and upsert as centrifuge; category and software = kraken2."""

    CATEGORY = "kraken2"
    DB_TYPE = "classification"
    SOFTWARE = "kraken2"
