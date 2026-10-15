import os
import shutil

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand, CommandError

from constants.constants import Televir_Metadata_Constants
from managing_files.models import ProcessControler
from pathogen_identification.utilities.centrifuge_db import (
    CentrifugeDBError, CentrifugeDBExporter)
from pathogen_identification.utilities.kraken2_db import (
    KRAKEN2_DEFAULT_KMER_LEN, KRAKEN2_DEFAULT_MINIMIZER_LEN, Kraken2DBBuilder,
    Kraken2DBError, Kraken2DBNaming, Kraken2DBRegistry, Kraken2DBStaging,
    Kraken2DBTaxonomy, validate_kmer_settings)
from utils.process_SGE import ProcessSched


class Command(BaseCommand):
    help = (
        "Build a Kraken 2 database from INSaFLU references (system files + the user's files), "
        "using NCBI taxonomy IDs. Output: <ref_db>/kraken2/kraken2_<YYYYMMDD>/"
    )

    def add_arguments(self, parser):
        parser.add_argument("--user_id", type=int, default=1, help="user id (default 1)")
        parser.add_argument("--threads", type=int, default=1, help="k2 build threads")
        parser.add_argument(
            "--kmer-len",
            type=int,
            default=KRAKEN2_DEFAULT_KMER_LEN,
            help=f"k-mer length (default {KRAKEN2_DEFAULT_KMER_LEN})",
        )
        parser.add_argument(
            "--minimizer-len",
            type=int,
            default=KRAKEN2_DEFAULT_MINIMIZER_LEN,
            help=f"minimizer length, <= 31 and <= kmer-len (default {KRAKEN2_DEFAULT_MINIMIZER_LEN})",
        )
        parser.add_argument(
            "--max-db-size",
            type=str,
            default=None,
            help="max hash table size, e.g. 10GiB (k2 downsamples the library to fit)",
        )
        parser.add_argument(
            "--include-phage",
            action="store_true",
            help="keep references whose description contains 'phage' (excluded by default)",
        )
        parser.add_argument(
            "--taxonomy-dir",
            type=str,
            default=None,
            help="dir with NCBI nodes.dmp, names.dmp (and optional merged.dmp) "
                 "(default: <metadata>/taxonomy)",
        )
        parser.add_argument(
            "--out-dir",
            type=str,
            default=None,
            help="database output dir (default: <ref_db>/kraken2)",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="query DB, check taxids and tools; write nothing",
        )
        parser.add_argument(
            "--force",
            action="store_true",
            help="run even if a build for this user is marked as running",
        )

    def handle(self, *args, **options):
        try:
            validate_kmer_settings(options["kmer_len"], options["minimizer_len"])
        except Kraken2DBError as e:
            raise CommandError(str(e))

        user = User.objects.get(pk=options["user_id"])
        process_controler = ProcessControler()
        process_sched = ProcessSched()
        process_name = process_controler.get_name_build_kraken2_db(user.pk)

        naming = Kraken2DBNaming(base_dir=options["out_dir"])
        taxonomy_dir = options["taxonomy_dir"] or os.path.join(
            Televir_Metadata_Constants.SOURCE["METAD"], "taxonomy"
        )

        if options["dry_run"]:
            self._dry_run(user, naming, taxonomy_dir, options)
            return

        running = ProcessControler.objects.filter(
            owner__id=user.pk,
            name=process_name,
            is_running=True,
            is_finished=False,
            is_error=False,
        ).exists()
        if running and not options["force"]:
            print("A kraken2 build is already running. Use --force to override.")
            return

        process_sched.set_process_controlers(user, process_name, 0)
        process_sched.set_process_controler(
            user, process_name, ProcessControler.FLAG_RUNNING
        )

        try:
            self._build(user, naming, taxonomy_dir, options)
        except Exception as e:
            print(f"ERROR: {e}")
            Kraken2DBBuilder.remove_partial(naming.tmp_base)
            process_sched.set_process_controler(
                user, process_name, ProcessControler.FLAG_ERROR
            )
            return

        process_sched.set_process_controler(
            user, process_name, ProcessControler.FLAG_FINISHED
        )
        print(f"Kraken 2 database written: {naming.final_base}")

    # ------------------------------------------------------------------ steps

    @staticmethod
    def _preflight():
        Kraken2DBBuilder.require_masker()
        Kraken2DBBuilder().resolve_binary()

    def _export_checks(self, user, taxonomy_dir, options):
        exporter = CentrifugeDBExporter(user, include_phage=options["include_phage"])
        plan = exporter.build_plan()

        if plan.conflicts:
            lines = [f"  {acc}: taxids {', '.join(t)}" for acc, t in plan.conflicts.items()]
            raise CentrifugeDBError(
                "Accessions mapped to more than one taxid:\n" + "\n".join(lines)
            )

        taxonomy = Kraken2DBTaxonomy(taxonomy_dir)
        taxonomy.check_files()
        unresolvable, remappable = taxonomy.classify_missing(plan.taxids)
        if unresolvable:
            raise CentrifugeDBError(
                "Taxids not in nodes.dmp/names.dmp and not in merged.dmp: "
                + ", ".join(unresolvable)
            )

        if plan.missing_files:
            raise CentrifugeDBError(
                "Reference files missing on disk: " + ", ".join(plan.missing_files)
            )

        if not plan.accid_to_taxid:
            raise CentrifugeDBError("No references selected for this user.")

        return exporter, plan, taxonomy, remappable

    def _dry_run(self, user, naming, taxonomy_dir, options):
        try:
            _, plan, _, remappable = self._export_checks(user, taxonomy_dir, options)
            self._preflight()
        except (CentrifugeDBError, Kraken2DBError) as e:
            print(f"DRY RUN FAILED: {e}")
            return

        if remappable:
            print(f"NOTE: {len(remappable)} retired taxids will be remapped by k2 (merged.dmp)")
        print("DRY RUN OK (no files written, sequences not read)")
        print(f"  database dir:        {naming.final_base}")
        print(f"  accessions:          {len(plan.accid_to_taxid)}")
        print(f"  distinct taxids:     {len(plan.taxids)}")
        print(f"  no taxid skipped:    {len(plan.no_taxid)}")
        print(f"  k / minimizer:       {options['kmer_len']} / {options['minimizer_len']}")
        print(f"  will replace db:     {naming.will_replace()}")
        record = Kraken2DBRegistry.build_record(naming, len(plan.accid_to_taxid), len(plan.taxids))
        print(f"  registry record:     {record}")

    def _build(self, user, naming, taxonomy_dir, options):
        self._preflight()
        exporter, plan, taxonomy, _ = self._export_checks(user, taxonomy_dir, options)

        os.makedirs(naming.base_dir, exist_ok=True)
        tmp_dir = naming.tmp_base
        Kraken2DBBuilder.remove_partial(tmp_dir)

        sequences, not_found = exporter.read_sequences(plan)
        if not_found:
            print(
                f"WARNING: {len(not_found)} accessions had no sequence in their file "
                f"and are skipped (e.g. {', '.join(not_found[:5])})"
            )

        written = Kraken2DBStaging.stage(plan, sequences, taxonomy, tmp_dir)
        if written == 0:
            raise CentrifugeDBError("No sequences were exported.")
        print(f"Exported {written} sequences, {len(plan.taxids)} taxids")

        builder = Kraken2DBBuilder(threads=options["threads"])
        log_path = os.path.join(tmp_dir, f"k2-build-{naming.name}.log")
        builder.build(
            tmp_dir,
            options["kmer_len"],
            options["minimizer_len"],
            log_path,
            max_db_size=options["max_db_size"],
        )
        builder.verify(tmp_dir)

        # intermediate inputs are not needed to classify (as in `kraken2-build --clean`)
        for item in ("library", "taxonomy", "seqid2taxid.map"):
            path = os.path.join(tmp_dir, item)
            if os.path.isdir(path):
                shutil.rmtree(path, ignore_errors=True)
            elif os.path.isfile(path):
                os.remove(path)

        Kraken2DBBuilder.promote(tmp_dir, naming.final_base)

        used_taxids = {plan.accid_to_taxid[a] for a in sequences if a in plan.accid_to_taxid}
        record = Kraken2DBRegistry.build_record(naming, written, len(used_taxids))
        try:
            Kraken2DBRegistry().register(record)
        except Exception as e:
            raise Kraken2DBError(
                f"Database written to {naming.final_base} but registration failed: {e}"
            ) from e
