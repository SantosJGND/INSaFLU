import datetime
import os
import shutil
import traceback

import pandas as pd
from constants.constants import Constants, FileType, TypePath
from django.contrib.auth.models import User
from django.db import IntegrityError, transaction
from managing_files.models import ProcessControler
from utils.process_SGE import ProcessSGE

from pathogen_identification.constants_settings import ConstantsSettings
from pathogen_identification.install_registry import Deployment_Params
from pathogen_identification.models import (
    ParameterSet,
    PIProject_Sample,
    Processed,
    Projects,
    SoftwareTree,
    SoftwareTreeNode,
    Submitted,
)
from pathogen_identification.modules.run_main import RunMain_class
from pathogen_identification.utilities.update_DBs import Update_Sample_Runs
from pathogen_identification.utilities.utilities_general import simplify_name
from pathogen_identification.utilities.utilities_pipeline import Utils_Manager


class PathogenIdentification_deployment:

    project: str
    prefix: str
    rdir: str
    threads: int = 3
    run_engine: RunMain_class
    params = dict
    run_params_db = pd.DataFrame()
    pk: int = 0
    username: str

    def __init__(
        self,
        pipeline_index: int,
        sample,  # sample name
        project: str = "test",
        prefix: str = "main",
        username: str = "admin",
        technology: str = "ONT",
        pk: int = 0,
        deployment_root_dir: str = "/tmp/insaflu/insaflu_something",
        dir_branch: str = "deployment",
    ) -> None:

        self.pipeline_index = pipeline_index

        self.username = username
        self.project = project
        self.sample = sample
        self.prefix = prefix
        self.deployment_root_dir = deployment_root_dir
        self.dir_branch = dir_branch
        self.dir = os.path.join(self.deployment_root_dir, dir_branch)

        self.prefix = prefix
        self.pk = pk
        self.technology = technology
        self.install_registry = Deployment_Params
        self.parameter_set = ParameterSet.objects.get(pk=pk)
        self.tree_makup = self.parameter_set.leaf.software_tree.global_index

    def input_read_project_path(self, filepath) -> str:
        """copy input reads to project directory and return new path"""

        if not os.path.isfile(filepath):
            return ""

        rname = os.path.basename(filepath)
        new_rpath = os.path.join(self.dir, "reads") + "/" + rname
        shutil.copy(filepath, new_rpath)
        return new_rpath

    def configure(self, r1_path: str, r2_path: str = "") -> bool:
        """generate config dictionary for run_main, and copy input reads to project directory."""
        self.get_constants()
        branch_exists = self.configure_params()

        if not branch_exists:
            return False

        self.generate_config_file()
        self.prep_test_env()

        new_r1_path = self.input_read_project_path(r1_path)
        new_r2_path = self.input_read_project_path(r2_path)

        self.config["sample_name"] = self.sample
        self.config["r1"] = new_r1_path
        self.config["r2"] = new_r2_path
        self.config["type"] = ["SE", "PE"][int(os.path.isfile(self.config["r2"]))]

        return True

    def get_constants(self):
        """set constants for technology"""
        if self.technology in "Illumina/IonTorrent":
            self.constants = ConstantsSettings.CONSTANTS_ILLUMINA
        if self.technology == "ONT":
            self.constants = ConstantsSettings.CONSTANTS_ONT

    def configure_params(self):
        """get pipeline parameters from database"""

        utils = Utils_Manager()

        all_paths = utils.get_all_technology_pipelines(self.technology, self.tree_makup)

        leaf_index = self.pipeline_index

        self.run_params_db = all_paths.get(self.pipeline_index, None)

        if self.run_params_db is None:
            print("Pipeline index not found")
            return False

        return True

    def generate_config_file(self):

        self.config = {
            "project": self.project,
            "source": self.install_registry.SOURCE,
            "technology": self.technology,
            "deployment_root_dir": self.deployment_root_dir,
            "sub_directory": self.dir_branch,
            "directories": {},
            "threads": self.threads,
            "prefix": self.prefix,
            "project_name": self.project,
            "metadata": {
                x: os.path.join(self.install_registry.METADATA["ROOT"], g)
                for x, g in self.install_registry.METADATA.items()
            },
            "technology": self.technology,
            "bin": self.install_registry.BINARIES,
            "actions": {},
        }

        for dr, g in ConstantsSettings.DIRS.items():
            self.config["directories"][dr] = os.path.join(self.dir, g)

        for dr, g in ConstantsSettings.ACTIONS.items():
            self.config["actions"][dr] = g

        self.config.update(self.constants)

    def prep_test_env(self):
        """
        from main directory bearing scripts, params.py and main.sh, create metagenome run directory

        :return:
        """
        os.makedirs(self.dir, exist_ok=True)
        os.makedirs(
            os.path.join(ConstantsSettings.media_directory, self.dir_branch),
            exist_ok=True,
        )
        os.makedirs(
            os.path.join(ConstantsSettings.static_directory, self.dir_branch),
            exist_ok=True,
        )

        for directory in self.config["directories"].values():
            os.makedirs(directory, exist_ok=True)

    def close(self):
        if os.path.exists(self.dir):
            shutil.rmtree(self.dir)

    def run_main_prep(self):

        self.run_engine = RunMain_class(self.config, self.run_params_db, self.username)


class Run_Main_from_Leaf:
    user: User
    file_r1: str
    file_r2: str
    file_sample: str
    technology: str
    project_name: str
    description: str
    date_created: str
    date_modified: str
    pk: int
    date_submitted = datetime.datetime

    container: PathogenIdentification_deployment

    def __init__(
        self,
        user: User,
        input_data: PIProject_Sample,
        project: Projects,
        pipeline_leaf: SoftwareTreeNode,
        pipeline_tree: SoftwareTree,
        odir: str,
    ):
        self.user = user
        self.sample = input_data
        self.project = project
        self.pipeline_leaf = pipeline_leaf
        self.pipeline_tree = pipeline_tree
        prefix = f"{simplify_name(input_data.name)}_{user.pk}_{project.pk}_{pipeline_leaf.pk}"
        self.date_submitted = datetime.datetime.now()

        self.file_r1 = input_data.sample.get_fastq_available(TypePath.MEDIA_ROOT, True)
        if input_data.sample.exist_file_2():
            self.file_r2 = input_data.sample.get_fastq_available(
                TypePath.MEDIA_ROOT, False
            )
        else:
            self.file_r2 = ""

        self.technology = input_data.sample.get_type_technology()
        self.project_name = project.name
        self.date_created = project.creation_date
        self.date_modified = project.last_change_date

        self.deployment_directory_structure = os.path.join(
            ConstantsSettings.televir_subdirectory,
            self.user.username,
            self.project_name,
            simplify_name(os.path.basename(input_data.name)),
            prefix,
        )

        self.unique_id = prefix

        self.deployment_directory = os.path.join(
            odir, self.deployment_directory_structure
        )

        self.parameter_set = self.register_parameter_set()
        self.pk = self.parameter_set.pk

        self.container = PathogenIdentification_deployment(
            pipeline_index=pipeline_leaf.index,
            sample=input_data.name,
            project=self.project_name,
            prefix=prefix,
            username=self.user.username,
            deployment_root_dir=odir,
            technology=self.technology,
            dir_branch=self.deployment_directory_structure,
            pk=self.pk,
        )

        self.is_available = self.check_availability()

    def get_status(self):
        return self.parameter_set.status

    def check_finished(self):
        return self.parameter_set.status == ParameterSet.STATUS_FINISHED

    def check_availability(self):
        return self.parameter_set.status not in [
            ParameterSet.STATUS_RUNNING,
            ParameterSet.STATUS_FINISHED,
        ]

    def get_in_line(self):

        self.parameter_set.status = ParameterSet.STATUS_QUEUED
        self.parameter_set.save()

    def check_submission(self):
        if self.parameter_set.status in [
            ParameterSet.STATUS_RUNNING,
        ]:

            return True

        else:
            return False

    def check_processed(self):
        if self.parameter_set.status in [
            ParameterSet.STATUS_FINISHED,
        ]:

            return True

        else:
            return False

    def register_parameter_set(self):

        try:
            new_run = ParameterSet.objects.get(
                leaf=self.pipeline_leaf,
                sample=self.sample,
                project=self.project,
            )

            return new_run

        except ParameterSet.DoesNotExist:
            new_run = ParameterSet.objects.create(
                leaf=self.pipeline_leaf,
                sample=self.sample,
                project=self.project,
            )

            return new_run

    def configure(self):
        configured = self.container.configure(
            self.file_r1,
            r2_path=self.file_r2,
        )
        return configured

    def set_run_process_running(self):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_project(self.unique_id),
            ProcessControler.FLAG_RUNNING,
        )

    def set_run_process_error(self):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_project(self.unique_id),
            ProcessControler.FLAG_ERROR,
        )

    def set_run_process_finished(self):
        process_controler = ProcessControler()
        process_SGE = ProcessSGE()
        process_SGE.set_process_controler(
            self.user,
            process_controler.get_name_televir_project(self.unique_id),
            ProcessControler.FLAG_FINISHED,
        )

    def Deploy(self):

        try:
            self.container.run_main_prep()
            self.container.run_engine.Run()
            self.container.run_engine.export_sequences()
            self.container.run_engine.Summarize()
            self.container.run_engine.generate_output_data_classes()
            self.container.run_engine.export_logdir()
            return True
        except Exception as e:
            print(traceback.format_exc())
            print(e)
            return False

    def Update_dbs(self):

        db_updated = Update_Sample_Runs(self.container.run_engine, self.parameter_set)

        return db_updated

    def register_submission(self):

        self.set_run_process_running()
        new_run = ParameterSet.objects.get(pk=self.pk)
        new_run.register_subprocess()

    def register_error(self):

        self.set_run_process_error()

        new_run = ParameterSet.objects.get(pk=self.pk)
        new_run.register_error()

    def register_completion(self):

        self.set_run_process_finished()

        new_run = ParameterSet.objects.get(pk=self.pk)
        new_run.register_finished()

    def update_project_change_date(self):
        self.project.last_change_date = datetime.datetime.now()
        self.project.save()

    def Submit(self):

        if not self.check_submission() and not self.check_processed():
            self.register_submission()
            configured = self.configure()

            if configured:
                run_success = self.Deploy()
            else:
                print("Error in configuration")
                self.register_error()
                return

            if run_success:
                update_successful = self.Update_dbs()
                if update_successful:
                    self.register_completion()
                    self.update_project_change_date()

                else:
                    print("Error in updating database")
                    self.register_error()

            else:
                print("Error in run")
                self.register_error()