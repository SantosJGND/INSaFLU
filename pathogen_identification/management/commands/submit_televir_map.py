import argparse
import logging
import os
import shutil

import pandas as pd
from django.core.management.base import BaseCommand

from constants.constants import Televir_Metadata_Constants as Televir_Metadata
from managing_files.models import ProcessControler
from pathogen_identification.constants_settings import MEDIA_ROOT, ConstantsSettings
from pathogen_identification.install_registry import Params_Illumina, Params_Nanopore
from pathogen_identification.models import FinalReport, Projects, RawReference
from pathogen_identification.modules.metadata_handler import RunMetadataHandler
from pathogen_identification.modules.object_classes import (
    Read_class,
    Sample_runClass,
    SoftwareDetail,
)
from pathogen_identification.modules.remap_class import (
    Mapping_Instance,
    Mapping_Manager,
)
from pathogen_identification.utilities.televir_parameters import TelevirParameters
from pathogen_identification.utilities.update_DBs import (
    Update_FinalReport,
    Update_ReferenceMap,
)
from pathogen_identification.utilities.utilities_general import simplify_name_lower
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager
from pathogen_identification.utilities.utilities_views import (
    ReportSorter,
    TelevirParameters,
)
from settings.constants_settings import ConstantsSettings as CS
from utils.process_SGE import ProcessSGE


class RunMain:
    remap_manager: Mapping_Manager
    mapping_instance: Mapping_Instance
    metadata_tool: RunMetadataHandler
    remap_params: TelevirParameters
    ##  metadata
    sift_query: str
    max_remap: int
    taxid_limit: int

    ## directories.
    root: str

    input_reads_dir: str
    filtered_reads_dir: str
    depleted_reads_dir: str

    log_dir: str

    dir_classification: str = f"classification_reports"
    dir_plots: str = f"plots"
    igv_dir: str = f"igv"

    def __init__(self, config: dict, method_args: pd.DataFrame, project: Projects):
        self.sample_name = config["sample_name"]
        self.type = config["type"]
        self.project_name = project.name
        self.username = project.owner.username
        self.prefix = "none"
        self.config = config
        self.taxid = config["taxid"]
        self.accid = config["accid"]
        self.threads = config["threads"]
        self.house_cleaning = False
        self.clean = config["clean"]
        self.project_pk = project.pk

        self.full_report = os.path.join(
            self.config["directories"][CS.PIPELINE_NAME_remapping], "full_report.csv"
        )

        self.logger_level_main = logging.INFO
        self.logger_level_detail = logging.ERROR
        self.logger = logging.getLogger("main {}".format(self.prefix))
        self.logger.setLevel(self.logger_level_main)

        logFormatter = logging.Formatter(
            fmt="{} {} %(levelname)s :%(message)s".format(
                config["sample_name"], self.prefix
            )
        )

        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)
        self.logger.propagate = False

        ######## DIRECTORIES ########

        self.deployment_root_dir = config["deployment_root_dir"]
        self.substructure_dir = config["sub_directory"]
        self.deployment_dir = os.path.join(
            self.deployment_root_dir, self.substructure_dir
        )

        self.media_dir = os.path.join(
            ConstantsSettings.media_directory, self.substructure_dir
        )
        self.static_dir = os.path.join(
            ConstantsSettings.static_directory, self.substructure_dir
        )

        self.media_dir_logdir = os.path.join(
            self.media_dir,
            "logs",
        )

        ####################################

        self.r1 = Read_class(
            config["r1"],
            config["directories"][CS.PIPELINE_NAME_read_quality_analysis],
            config["directories"]["reads_enriched_dir"],
            config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(
                config["bin"], CS.PIPELINE_NAME_read_quality_analysis
            ),
        )

        self.r2 = Read_class(
            config["r2"],
            config["directories"][CS.PIPELINE_NAME_read_quality_analysis],
            config["directories"]["reads_enriched_dir"],
            config["directories"]["reads_depleted_dir"],
            bin=get_bindir_from_binaries(
                config["bin"], CS.PIPELINE_NAME_read_quality_analysis
            ),
        )

        self.sample = Sample_runClass(
            self.r1,
            self.r2,
            self.sample_name,
            self.project_name,
            self.username,
            self.config["technology"],
            self.type,
            0,
            ",".join(
                [os.path.basename(self.r1.current), os.path.basename(self.r2.current)]
            ),
            bin=get_bindir_from_binaries(
                config["bin"], CS.PIPELINE_NAME_read_quality_analysis
            ),
            threads=self.threads,
        )

        ### mapping parameters
        self.min_scaffold_length = config["assembly_contig_min_length"]
        self.minimum_coverage = int(config["minimum_coverage_threshold"])
        self.maximum_coverage = 1000000000

        ### metadata
        self.remap_params = TelevirParameters.get_remap_software(
            self.username, self.project_name
        )

        self.metadata_tool = RunMetadataHandler(
            self.username,
            self.config,
            sift_query=config["sift_query"],
            prefix=self.prefix,
        )

        self.max_remap = self.remap_params.max_accids
        self.taxid_limit = self.remap_params.max_taxids

        ### methods
        self.remapping_method = SoftwareDetail(
            CS.PIPELINE_NAME_remapping,
            method_args,
            config,
            self.prefix,
        )

        ###

        self.media_dir_classification = os.path.join(
            self.media_dir,
            self.dir_classification,
        )

        self.static_dir_plots = os.path.join(
            self.substructure_dir,
            self.dir_plots,
        )

        self.media_dir_igv = os.path.join(
            self.static_dir,
            self.igv_dir,
        )

        os.makedirs(
            self.media_dir_classification,
            exist_ok=True,
        )

        os.makedirs(
            os.path.join(ConstantsSettings.static_directory, self.static_dir_plots),
            exist_ok=True,
        )

        os.makedirs(
            self.media_dir_igv,
            exist_ok=True,
        )

        self.filtered_reads_dir = config["directories"][
            CS.PIPELINE_NAME_read_quality_analysis
        ]
        self.log_dir = config["directories"]["log_dir"]

    def generate_targets(self):
        result_df = pd.DataFrame(columns=["qseqid", "taxid"])
        if self.taxid:
            result_df = pd.DataFrame([{"taxid": self.taxid, "qseqid": ""}])
        elif self.accid:
            result_df = pd.DataFrame([{"qseqid": self.accid, "taxid": ""}])

        self.metadata_tool.match_and_select_targets(
            result_df,
            pd.DataFrame(columns=["qseqid", "taxid"]),
            self.max_remap,
            self.taxid_limit,
        )

    def deploy_REMAPPING(self):
        """ """
        self.remap_manager = Mapping_Manager(
            self.metadata_tool.remap_targets,
            self.sample.r1,
            self.sample.r2,
            self.remapping_method,
            "Dummy",
            self.type,
            self.prefix,
            self.threads,
            get_bindir_from_binaries(self.config["bin"], CS.PIPELINE_NAME_remapping),
            self.logger_level_detail,
            self.house_cleaning,
            remap_params=self.remap_params,
            logdir=self.config["directories"]["log_dir"],
        )
        self.logger.info(
            f"{self.prefix} remapping # targets: {len(self.metadata_tool.remap_targets)}"
        )

        print("moving to : ", self.static_dir_plots)
        print("moving to : ", self.media_dir_igv)

        self.remap_manager.run_mappings_move_clean(
            self.static_dir_plots, self.media_dir_igv
        )
        self.remap_manager.export_reference_fastas_if_failed(self.media_dir_igv)

        self.remap_manager.merge_mapping_reports()
        self.remap_manager.collect_final_report_summary_statistics()

    def run(self):
        self.deploy_REMAPPING()
        print("remap_manager.report")
        self.report = self.remap_manager.report
        self.export_final_reports()

    def export_final_reports(self):
        ### main report
        self.report.to_csv(
            self.full_report,
            index=False,
            sep="\t",
            header=True,
        )


def get_bindir_from_binaries(binaries, key, value: str = ""):
    if value == "":
        try:
            return os.path.join(binaries["ROOT"], binaries[key]["default"], "bin")
        except KeyError:
            return ""
    else:
        try:
            return os.path.join(binaries["ROOT"], binaries[key][value], "bin")
        except KeyError:
            return ""


class Input_Generator:
    method_args: pd.DataFrame

    def __init__(self, reference: RawReference, output_dir: str, threads: int = 4):
        self.utils = Utils_Manager()
        self.reference = reference
        self.install_registry = Televir_Metadata()

        self.dir_branch = os.path.join(
            ConstantsSettings.televir_subdirectory,
            f"{reference.run.project.owner.pk}",
            f"{reference.run.project.pk}",
            f"{reference.run.sample.pk}",
            "remapping",
            f"{reference.pk}",
        )

        self.technology = reference.run.project.technology
        self.threads = threads
        self.prefix = reference.accid
        self.project = reference.run.project.name
        self.clean = False
        self.deployment_root_dir = os.path.join(
            output_dir, "remapping", f"ref_{reference.pk}"
        )
        self.dir = os.path.join(self.deployment_root_dir, self.dir_branch)

        os.makedirs(self.dir, exist_ok=True)

        self.r1_path = reference.run.sample.sample.path_name_1.path
        self.r2_path = (
            reference.run.sample.sample.path_name_2.path
            if reference.run.sample.sample.exist_file_2()
            else ""
        )

        self.taxid = reference.taxid
        self.accid = reference.accid

        if self.technology == "ONT":
            self.params = Params_Nanopore
        else:
            self.params = Params_Illumina

    def input_read_project_path(self, filepath):
        if not os.path.isfile(filepath):
            return ""
        rname = os.path.basename(filepath)

        new_rpath = os.path.join(self.dir, "reads") + "/" + rname
        shutil.copy(filepath, new_rpath)
        return new_rpath

    def generate_method_args(self):
        parameter_leaf = self.reference.run.parameter_set.leaf
        run_df = self.utils.get_leaf_parameters(parameter_leaf)

        self.method_args = run_df[run_df.module == CS.PIPELINE_NAME_remapping]

        if self.method_args.empty:
            raise ValueError(
                f"no remapping parameters found for {self.reference.accid} in leaf {parameter_leaf}"
            )

    def generate_config_file(self):
        self.config = {
            "sample_name": simplify_name_lower(
                os.path.basename(self.r1_path).replace(".fastq.gz", "")
            ),
            "source": self.install_registry.SOURCE,
            "technology": self.technology,
            "deployment_root_dir": self.deployment_root_dir,
            "sub_directory": self.dir_branch,
            "directories": {},
            "threads": self.threads,
            "prefix": self.prefix,
            "project_name": self.project,
            "metadata": self.install_registry.metadata_full_path,
            "bin": self.install_registry.BINARIES,
            "taxid": self.taxid,
            "accid": self.accid,
            "clean": self.clean,
        }

        for dr, g in ConstantsSettings.DIRS.items():
            self.config["directories"][dr] = os.path.join(self.dir, g)
            os.makedirs(self.config["directories"][dr], exist_ok=True)

        self.config["r1"] = self.input_read_project_path(self.r1_path)
        self.config["r2"] = self.input_read_project_path(self.r2_path)
        self.config["type"] = [
            ConstantsSettings.SINGLE_END,
            ConstantsSettings.PAIR_END,
        ][int(os.path.isfile(self.config["r2"]))]

        self.config.update(self.params.CONSTANTS)

    def update_raw_reference_status_mapped(self):
        self.reference.status = RawReference.STATUS_MAPPED
        self.reference.save()

    def update_raw_reference_status_fail(self):
        self.reference.status = RawReference.STATUS_FAIL
        self.reference.save()

    def update_final_report(self, run_class: RunMain):
        run = self.reference.run
        sample = run.sample

        Update_FinalReport(run_class, run, sample)

        for ref_map in run_class.remap_manager.mapped_instances:
            Update_ReferenceMap(ref_map, run, sample)

    def run_reference_overlap_analysis(self):
        run = self.reference.run
        sample = run.sample
        final_report = FinalReport.objects.filter(sample=sample, run=run).order_by(
            "-coverage"
        )
        #
        if sample is None:
            return
        report_layout_params = TelevirParameters.get_report_layout_params(run_pk=run.pk)
        report_sorter = ReportSorter(sample, final_report, report_layout_params)
        report_sorter.sort_reports_save()


class Command(BaseCommand):
    help = "deploy run"

    def add_arguments(self, parser):
        parser.add_argument(
            "--ref_id",
            type=int,
            help="user deploying the run (pk)",
        )

        parser.add_argument(
            "-o",
            "--outdir",
            type=str,
            help="output directory",
        )

    def handle(self, *args, **options):
        ###
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()

        raw_reference_id = int(options["ref_id"])

        reference = RawReference.objects.get(pk=raw_reference_id)

        if reference is None:
            print("reference not found")
            return

        if reference.run is None:
            print("run not found")
            return

        if reference.run.project is None:
            print("project not found")
            return

        project_name = reference.run.project.name
        user = reference.run.project.owner
        project_name = reference.run.project.name

        ######## register map
        process_SGE.set_process_controler(
            user,
            process_controler.get_name_televir_map(reference.pk),
            ProcessControler.FLAG_RUNNING,
        )

        ########

        input_generator = Input_Generator(
            reference, options["outdir"], threads=ConstantsSettings.MAPPING_THREADS
        )

        try:
            input_generator.generate_method_args()
            input_generator.generate_config_file()
            print("config generated")

            run_engine = RunMain(
                input_generator.config,
                input_generator.method_args,
                reference.run.project,
            )
            run_engine.generate_targets()
            run_engine.run()

            input_generator.update_raw_reference_status_mapped()
            input_generator.update_final_report(run_engine)
            print("done")
            input_generator.run_reference_overlap_analysis()

            ######## register map sucess
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_map(reference.pk),
                ProcessControler.FLAG_FINISHED,
            )

            ########

        except Exception as e:
            print(e)
            input_generator.update_raw_reference_status_fail()

            ######## register map sucess
            process_SGE.set_process_controler(
                user,
                process_controler.get_name_televir_map(reference.pk),
                ProcessControler.FLAG_ERROR,
            )

            ########
