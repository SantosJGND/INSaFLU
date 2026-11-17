import os
import shutil

from django.contrib.auth.models import User
from django.core.management.base import BaseCommand

from constants.constants import Televir_Metadata_Constants
from managing_files.models import ProcessControler
from pathogen_identification.utilities.centrifuge_db import (
    CentrifugeDBBuilder, CentrifugeDBError, CentrifugeDBExporter,
    CentrifugeDBNaming, CentrifugeDBRegistry, CentrifugeDBTaxonomy)
from utils.process_SGE import ProcessSched


class Command(BaseCommand):
    help = (
        "Build a Centrifuge index from INSaFLU references (system files + the user's files), "
        "using NCBI taxonomy IDs. Output: <ref_db>/centrifuge/centrifuge_<YYYYMMDD>.*.cf"
    )

    def add_arguments(self, parser):
        parser.add_argument("--user_id", type=int, default=1, help="user id (default 1)")
        parser.add_argument("--threads", type=int, default=1, help="centrifuge-build threads")
        parser.add_argument(
            "--include-phage",
            action="store_true",
            help="keep references whose description contains 'phage' (excluded by default)",
        )
        parser.add_argument(
            "--taxonomy-dir",
            type=str,
            default=None,
            help="dir with NCBI nodes.dmp and names.dmp (default: <metadata>/taxonomy)",
        )
        parser.add_argument(
            "--out-dir",
            type=str,
            default=None,
            help="index output dir (default: <ref_db>/centrifuge)",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="query DB and check conflicts/taxids only; write nothing",
        )
        parser.add_argument(
            "--force",
            action="store_true",
            help="run even if a build for this user is marked as running",
        )

    def handle(self, *args, **options):
        user = User.objects.get(pk=options["user_id"])
        process_controler = ProcessControler()
        process_sched = ProcessSched()
        process_name = process_controler.get_name_build_centrifuge_db(user.pk)

        naming = CentrifugeDBNaming(base_dir=options["out_dir"])
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
            print("A centrifuge build is already running. Use --force to override.")
            return

        process_sched.set_process_controlers(user, process_name, 0)
        process_sched.set_process_controler(
            user, process_name, ProcessControler.FLAG_RUNNING
        )

        try:
            self._build(user, naming, taxonomy_dir, options)
        except Exception as e:
            print(f"ERROR: {e}")
            naming_tmp = naming.tmp_base
            CentrifugeDBBuilder.remove_partial(naming_tmp)
            if os.path.isdir(naming.work_dir):
                shutil.rmtree(naming.work_dir, ignore_errors=True)
            process_sched.set_process_controler(
                user, process_name, ProcessControler.FLAG_ERROR
            )
            return

        process_sched.set_process_controler(
            user, process_name, ProcessControler.FLAG_FINISHED
        )
        print(f"Centrifuge index written: {naming.final_base}")

    # ------------------------------------------------------------------ steps

    def _export_checks(self, user, taxonomy_dir, options):
        exporter = CentrifugeDBExporter(user, include_phage=options["include_phage"])
        plan = exporter.build_plan()

        if plan.conflicts:
            lines = [f"  {acc}: taxids {', '.join(t)}" for acc, t in plan.conflicts.items()]
            raise CentrifugeDBError(
                "Accessions mapped to more than one taxid:\n" + "\n".join(lines)
            )

        taxonomy = CentrifugeDBTaxonomy(taxonomy_dir)
        missing = taxonomy.missing_taxids(plan.taxids)
        if missing:
            raise CentrifugeDBError(
                "Taxids not in nodes.dmp/names.dmp (scientific name): " + ", ".join(missing)
            )

        if plan.missing_files:
            raise CentrifugeDBError(
                "Reference files missing on disk: " + ", ".join(plan.missing_files)
            )

        if not plan.accid_to_taxid:
            raise CentrifugeDBError("No references selected for this user.")

        return exporter, plan, taxonomy

    def _dry_run(self, user, naming, taxonomy_dir, options):
        try:
            _, plan, _ = self._export_checks(user, taxonomy_dir, options)
        except CentrifugeDBError as e:
            print(f"DRY RUN FAILED: {e}")
            return

        print("DRY RUN OK (no files written, sequences not read)")
        print(f"  index base:          {naming.final_base}")
        print(f"  accessions:          {len(plan.accid_to_taxid)}")
        print(f"  distinct taxids:     {len(plan.taxids)}")
        print(f"  no taxid skipped:    {len(plan.no_taxid)}")
        print(f"  will replace index:  {naming.will_replace()}")
        record = CentrifugeDBRegistry.build_record(naming, len(plan.accid_to_taxid), len(plan.taxids))
        print(f"  registry record:     {record}")

    def _build(self, user, naming, taxonomy_dir, options):
        exporter, plan, taxonomy = self._export_checks(user, taxonomy_dir, options)

        os.makedirs(naming.base_dir, exist_ok=True)
        os.makedirs(naming.work_dir, exist_ok=True)

        input_fna = os.path.join(naming.work_dir, "input.fna")
        conversion_tsv = os.path.join(naming.work_dir, "conversion.tsv")
        log_path = os.path.join(naming.work_dir, f"centrifuge-build-{naming.name}.log")

        sequences, not_found = exporter.read_sequences(plan)
        if not_found:
            print(
                f"WARNING: {len(not_found)} accessions had no sequence in their file "
                f"and are skipped (e.g. {', '.join(not_found[:5])})"
            )
        written = CentrifugeDBExporter.write(plan, sequences, input_fna, conversion_tsv)
        if written == 0:
            raise CentrifugeDBError("No sequences were exported.")
        print(f"Exported {written} sequences, {len(plan.taxids)} taxids")

        builder = CentrifugeDBBuilder(threads=options["threads"])
        CentrifugeDBBuilder.remove_partial(naming.tmp_base)
        builder.build(input_fna, conversion_tsv, taxonomy, naming.tmp_base, log_path)
        builder.verify(naming.tmp_base, written)
        CentrifugeDBBuilder.promote(naming.tmp_base, naming.final_base)

        # register only after the index is in its final place
        used_taxids = {plan.accid_to_taxid[a] for a in sequences if a in plan.accid_to_taxid}
        record = CentrifugeDBRegistry.build_record(naming, written, len(used_taxids))
        try:
            CentrifugeDBRegistry().register(record)
        except Exception as e:
            raise CentrifugeDBError(
                f"Index written to {naming.final_base} but registration failed: {e}"
            ) from e

        shutil.rmtree(naming.work_dir, ignore_errors=True)
