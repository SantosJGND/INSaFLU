import mimetypes
import os
from datetime import datetime
from typing import Dict, List, Optional

import pandas as pd
from Bio import SeqIO
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.files.temp import NamedTemporaryFile
from django.http import HttpResponse, JsonResponse
from django.template.loader import render_to_string
from django.utils.safestring import mark_safe
from django.utils.translation import ugettext_lazy as _
from django.views.decorators.csrf import csrf_protect
from django.views.decorators.http import require_POST

from constants.constants import Constants, FileExtensions, FileType, TypePath
from constants.meta_key_and_values import MetaKeyAndValue
from constants.software_names import SoftwareNames
from fluwebvirus.settings import BASE_DIR, STATIC_ROOT, STATIC_URL
from managing_files.models import ProcessControler
from managing_files.models import ProjectSample as InsafluProjectSample
from pathogen_identification.constants_settings import \
    ConstantsSettings as PICS
from pathogen_identification.models import (FinalReport, ParameterSet,
                                            PIProject_Sample, Projects,
                                            RawReference, ReferenceMap_Main,
                                            ReferencePanel,
                                            ReferenceSourceFileMap, RunMain,
                                            SoftwareTreeNode, TelefluMapping,
                                            TeleFluProject, TeleFluSample)
from pathogen_identification.tables import ReferenceSourceTable
from pathogen_identification.utilities.reference_utils import (
    check_file_reference_submitted, check_raw_reference_submitted,
    check_user_reference_exists, create_combined_reference)
from pathogen_identification.utilities.televir_bioinf import TelevirBioinf
from pathogen_identification.utilities.televir_parameters import \
    TelevirParameters
from pathogen_identification.utilities.utilities_general import \
    get_services_dir
from pathogen_identification.utilities.utilities_pipeline import (
    SoftwareTreeUtils, Utils_Manager)
from pathogen_identification.utilities.utilities_views import (
    RawReferenceUtils, ReportSorter, SampleReferenceManager,
    set_control_reports)
from pathogen_identification.views import inject__added_references
from settings.constants_settings import ConstantsSettings as CS
from utils.process_SGE import ProcessSGE
from utils.software import Software
from utils.utils import Utils


def simplify_name(name: str):
    return (
        name.replace("_", "_")
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
        .lower()
    )


@login_required
@require_POST
def submit_sample_metagenomics_televir(request):
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False, "no_references": False}

        process_SGE = ProcessSGE()

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))

        user = sample.project.owner
        project = sample.project

        software_utils = SoftwareTreeUtils(user, project, sample=sample)
        runs_to_deploy = software_utils.check_runs_to_submit_metagenomics_sample(sample)
        reference_manager = SampleReferenceManager(sample)
        reference_utils = RawReferenceUtils(sample)

        count_references = reference_utils.query_sample_compound_references_regressive()

        if count_references.exists() is False:
            data["no_references"] = True
            return JsonResponse(data)

        try:
            if len(runs_to_deploy) > 0:
                for sample, leaves_to_deploy in runs_to_deploy.items():
                    for leaf in leaves_to_deploy:
                        metagenomics_run = reference_manager.mapping_run_from_leaf(leaf)
                        _ = process_SGE.set_submit_televir_sample_metagenomics(
                            user=request.user,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                            combined_analysis=True,
                            map_run_pk=metagenomics_run.pk,
                        )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def submit_sample_screening_televir(request):
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))

        user = sample.project.owner
        project = sample.project
        reference_manager = SampleReferenceManager(sample)

        software_utils = SoftwareTreeUtils(user, project, sample=sample)
        runs_to_deploy = software_utils.check_runs_to_submit_screening_sample(sample)

        reference_utils = RawReferenceUtils(sample)

        count_references = reference_utils.query_sample_compound_references_regressive()

        if count_references.exists() is False:
            data["no_references"] = True
            return JsonResponse(data)

        try:
            if len(runs_to_deploy) > 0:
                for sample, leaves_to_deploy in runs_to_deploy.items():
                    for leaf in leaves_to_deploy:
                        screening_run = reference_manager.screening_run_from_leaf(leaf)
                        _ = process_SGE.set_submit_televir_sample_metagenomics(
                            user=request.user,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                            map_run_pk=screening_run.pk,
                        )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


def check_reference_mapped(sample_id, reference: RawReference):
    return RawReference.objects.filter(
        accid=reference.accid,
        run__status__in=[
            RunMain.STATUS_PREP,
            RunMain.STATUS_RUNNING,
            RunMain.STATUS_FINISHED,
        ],
        run__run_type__in=[RunMain.RUN_TYPE_MAP_REQUEST],
        run__sample__pk=sample_id,
        run__parameter_set__status__in=[
            ParameterSet.STATUS_QUEUED,
            ParameterSet.STATUS_RUNNING,
            ParameterSet.STATUS_FINISHED,
        ],
    ).exists()


def deploy_remap(
    sample: PIProject_Sample, project: Projects, reference_id_list: list = []
):
    process_SGE = ProcessSGE()

    data = {"is_ok": True, "is_deployed": False, "is_empty": False, "message": ""}
    user = sample.project.owner
    sample_id = sample.pk
    reference_manager = SampleReferenceManager(sample)

    ### check if all references are already mapped
    # reference_id_list = request.POST.getlist("reference_ids[]")
    added_references = RawReference.objects.filter(
        run__sample__pk=sample_id, run__run_type=RunMain.RUN_TYPE_STORAGE
    )

    if len(reference_id_list) == 0 and len(added_references) == 0:
        data["is_empty"] = True
        data["message"] = "No references to deploy"
        return data

    #### check among reference id list

    references_mapped_dict: Dict[str, list] = {}

    for reference_id in reference_id_list:
        reference = RawReference.objects.get(pk=int(reference_id))
        if check_reference_mapped(sample_id, reference):

            if reference.accid is None:
                continue
            indices = (
                RawReference.objects.filter(
                    accid=reference.accid,
                    run__status__in=[
                        RunMain.STATUS_RUNNING,
                        RunMain.STATUS_FINISHED,
                    ],
                    run__run_type__in=[RunMain.RUN_TYPE_MAP_REQUEST],
                    run__sample__pk=sample_id,
                )
                .values_list("run__parameter_set__leaf__index", flat=True)
                .distinct()
            )

            references_mapped_dict[reference.accid] = list(indices)

    ##### check among added references

    for added_reference in added_references:
        if check_reference_mapped(sample_id, added_reference):

            ref_workflows = (
                RawReference.objects.filter(
                    accid=added_reference.accid,
                    run__status__in=[
                        RunMain.STATUS_RUNNING,
                        RunMain.STATUS_FINISHED,
                    ],
                    run__run_type__in=[RunMain.RUN_TYPE_MAP_REQUEST],
                    run__sample__pk=sample_id,
                )
                .values_list("run__parameter_set__leaf__index", flat=True)
                .distinct()
            )

            if added_reference.accid in references_mapped_dict:
                references_mapped_dict[added_reference.accid].extend(
                    list(ref_workflows)
                )
            else:
                references_mapped_dict[added_reference.accid] = list(ref_workflows)

    software_utils = SoftwareTreeUtils(user, project, sample=sample)
    runs_to_deploy, workflow_deployed_dict = (
        software_utils.check_runs_to_submit_mapping_only(sample)
    )

    if len(runs_to_deploy) == 0:
        return data

    try:

        if len(runs_to_deploy) > 0:
            deployed = 0
            deployed_refs = {
                sample: {
                    leaf: {"deployed": [], "not_deployed": []}
                    for leaf in runs_to_deploy[sample]
                }
                for sample in runs_to_deploy
            }
            for sample, leaves_to_deploy in runs_to_deploy.items():

                reference_utils = RawReferenceUtils(sample)

                count_references = (
                    reference_utils.query_sample_compound_references_regressive()
                )

                if count_references.exists() is False:
                    continue

                for leaf in leaves_to_deploy:

                    references_added = []
                    references_to_add: List[RawReference] = []

                    for reference_id in reference_id_list:
                        reference = RawReference.objects.get(pk=int(reference_id))
                        if reference.accid in references_mapped_dict:
                            if leaf.index in references_mapped_dict[reference.accid]:
                                deployed_refs[sample][leaf]["not_deployed"].append(
                                    reference.accid
                                )
                                continue

                        references_to_add.append(reference)
                        reference.save()

                        references_added.append(reference.accid)

                    for added_reference in added_references:
                        if added_reference.accid in references_mapped_dict:
                            if (
                                leaf.index
                                in references_mapped_dict[added_reference.accid]
                            ):
                                deployed_refs[sample][leaf]["not_deployed"].append(
                                    added_reference.accid
                                )
                                continue
                        added_reference.pk = None

                        references_to_add.append(added_reference)
                        references_added.append(added_reference.accid)

                    if len(references_added) == 0:
                        continue

                    mapping_run = reference_manager.mapping_request_run_from_leaf(leaf)
                    for reference in references_to_add:
                        reference.pk = None
                        reference.run = mapping_run
                        reference.save()

                    taskID = process_SGE.set_submit_televir_sample_metagenomics(
                        user=user,
                        sample_pk=sample.pk,
                        leaf_pk=leaf.pk,
                        mapping_request=True,
                        map_run_pk=mapping_run.pk,
                    )
                    deployed += 1
                    deployed_refs[sample][leaf]["deployed"] = references_added

            data["is_deployed"] = True
            message = f"Sample {sample.name}, deployed {deployed} run{['', 's'][int(deployed != 1)]}. "
            total_runs_deployed = 0

            for sample, leaves_to_deploy in runs_to_deploy.items():
                for leaf in leaves_to_deploy:
                    if len(deployed_refs[sample][leaf]["deployed"]) > 0:
                        message += f"Leaf {leaf.index}, deployed references: {', '.join(deployed_refs[sample][leaf]['deployed'])}. "
                        total_runs_deployed += 1

                    if len(deployed_refs[sample][leaf]["not_deployed"]) > 0:
                        message += f"References: {', '.join(deployed_refs[sample][leaf]['not_deployed'])} already mapped for Leaf {leaf.index}. "

            if total_runs_deployed == 0:
                data["is_deployed"] = False

            data["message"] = message
            data["is_ok"] = True

    except Exception as e:
        print(e)
        data["is_ok"] = False

    return data


@login_required
@require_POST
def submit_sample_mapping_televir(request):
    if request.is_ajax():
        data = {"is_ok": True, "is_deployed": False, "is_empty": False, "message": ""}

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = sample.project
        reference_id_list = request.POST.getlist("reference_ids[]")

        data = deploy_remap(sample, project, reference_id_list)

        return JsonResponse(data)

    return JsonResponse({"is_ok": False})


@login_required
@require_POST
def available_televir_files(request):

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False, "files": []}

        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=user_id)
        files = ReferenceSourceFile.objects.filter(
            owner=user, is_deleted=False
        ).distinct("file")

        data["is_ok"] = True
        data["files"] = {file.pk: file.file for file in files}

        return JsonResponse(data)


@login_required
@require_POST
def submit_sample_mapping_panels(request):
    if request.is_ajax():
        process_SGE = ProcessSGE()
        user = request.user
        data = {
            "is_ok": True,
            "is_deployed": False,
            "is_empty": False,
            "params_empty": False,
            "message": "",
        }

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))

        project = sample.project
        software_utils = SoftwareTreeUtils(user, project, sample=sample)
        runs_to_deploy, _ = software_utils.check_runs_to_submit_mapping_only(sample)

        if len(runs_to_deploy) == 0:
            data["params_empty"] = True
            return JsonResponse(data)

        sample_panels = sample.panels_added

        if len(sample_panels) == 0:
            data["is_empty"] = True
            return JsonResponse(data)

        try:
            for sample, leaves_to_deploy in runs_to_deploy.items():
                for leaf in leaves_to_deploy:
                    for panel in sample_panels:
                        if panel.is_deleted:
                            continue

                        taskID = process_SGE.set_submit_televir_sample_panel_map(
                            user=request.user,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                            mapping_request=True,
                            panel_pk=panel.pk,
                        )
                        data["is_deployed"] = True

        except Exception as e:
            print(e)
            print("error")

            data["is_ok"] = False
        return JsonResponse(data)


@login_required
@require_POST
def submit_samples_mapping_panels(request):
    if request.is_ajax():
        data = {
            "is_ok": True,
            "is_deployed": False,
            "is_empty": False,
            "samples_deployed": 0,
            "message": "",
        }

        process_SGE = ProcessSGE()

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))
        user = request.user

        project_samples = PIProject_Sample.objects.filter(project=project)

        sample_ids = request.POST.getlist("sample_ids[]")
        check_box_all_checked = request.POST.get("check_box_all_checked", False)
        if check_box_all_checked:
            sample_ids = []
        else:
            sample_ids = [int(sample_id) for sample_id in sample_ids]

        if len(sample_ids) > 0:
            project_samples = project_samples.filter(pk__in=sample_ids)

        try:
            samples_submitted = 0
            errors = ""

            for sample in project_samples:
                reference_manager = SampleReferenceManager(sample)

                software_utils = SoftwareTreeUtils(user, project, sample=sample)
                runs_to_deploy, workflow_deployed_dict = (
                    software_utils.check_runs_to_submit_mapping_only(sample)
                )

                if len(runs_to_deploy) == 0:
                    continue

                sample_panels = sample.panels_added

                if len(sample_panels) == 0:
                    data["is_empty"] = True
                    continue

                try:
                    for sample, leaves_to_deploy in runs_to_deploy.items():
                        for leaf in leaves_to_deploy:
                            for panel in sample_panels:
                                if panel.is_deleted:
                                    continue

                                references = RawReference.objects.filter(panel=panel)
                                run_panel_copy = reference_manager.copy_panel(panel)

                                panel_mapping_run = reference_manager.mapping_request_panel_run_from_leaf(
                                    leaf, panel_pk=run_panel_copy.pk
                                )
                                for reference in references:
                                    reference.pk = None
                                    reference.run = panel_mapping_run
                                    reference.panel = run_panel_copy
                                    reference.save()

                                taskID = (
                                    process_SGE.set_submit_televir_sample_metagenomics(
                                        user=request.user,
                                        sample_pk=sample.pk,
                                        leaf_pk=leaf.pk,
                                        mapping_request=True,
                                        map_run_pk=panel_mapping_run.pk,
                                    )
                                )
                                data["is_deployed"] = True
                                samples_submitted += 1

                except Exception as e:
                    print(e)
                    print("error")

                    data["is_ok"] = False
                    errors += f" Error deploying sample {sample.name}"

        except Exception as e:
            print(e)
            print("error")

            data["is_ok"] = False
            data["message"] = errors

        data["message"] = f"Deployed {samples_submitted} samples. {errors}"

        return JsonResponse(data)


@login_required
@require_POST
def submit_project_samples_mapping_televir(request):
    if request.is_ajax():
        data = {
            "is_ok": True,
            "is_deployed": False,
            "is_empty": False,
            "samples_deployed": 0,
            "message": "",
        }

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        project_samples = PIProject_Sample.objects.filter(project=project)

        sample_ids = request.POST.getlist("sample_ids[]")
        sample_ids = [int(sample_id) for sample_id in sample_ids]

        if len(sample_ids) > 0:
            project_samples = project_samples.filter(pk__in=sample_ids)

        try:
            data_dump = []
            for sample in project_samples:
                data = deploy_remap(sample, project, [])
                data_dump.append(data)

            data["is_deployed"] = any([d["is_deployed"] for d in data_dump])
            data["is_ok"] = True
            data["samples_deployed"] = sum([d["is_deployed"] for d in data_dump])

            data["message"] = "\n".join(
                [d["message"] for d in data_dump if d["message"]]
            )

        except Exception as e:
            print(e)
            data["is_ok"] = False

        return JsonResponse(data)


@login_required
@require_POST
def get_all_samples_selected(request):
    """
    get all sample ids for selected samples
    """
    if request.is_ajax():
        data = {"is_ok": False, "sample_ids": []}

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        sample_ids = PIProject_Sample.objects.filter(
            project=project, is_deleted_in_file_system=False
        ).values_list("pk", flat=True)
        sample_ids = [str(sample_id) for sample_id in sample_ids]

        if len(sample_ids) > 0:
            data["sample_ids"] = sample_ids

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def deploy_ProjectPI(request):
    """
    prepare data for deployment of pathogen identification.
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=int(user_id))

        samples = PIProject_Sample.objects.filter(
            project=project, is_deleted_in_file_system=False
        )

        sample_ids = request.POST.getlist("sample_ids[]")
        sample_ids = [int(sample_id) for sample_id in sample_ids]
        check_box_all_checked = request.POST.get("check_box_all_checked", False)
        if check_box_all_checked:
            samples = samples.filter(is_deleted_in_file_system=False)
        elif len(sample_ids) > 0:
            samples = samples.filter(pk__in=sample_ids)

        software_utils = SoftwareTreeUtils(user, project)

        try:
            for sample in samples:

                runs_to_deploy = software_utils.check_runs_to_deploy_sample(sample)

                if len(runs_to_deploy) > 0:
                    for sample, _ in runs_to_deploy.items():
                        taskID = process_SGE.set_submit_televir_sample(
                            user=user,
                            project_pk=project.pk,
                            sample_pk=sample.pk,
                        )

                    data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def deploy_ProjectPI_runs(request):
    """
    prepare data for deployment of pathogen identification.
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=int(user_id))

        software_utils = SoftwareTreeUtils(user, project)
        runs_to_deploy = software_utils.check_runs_to_deploy_project()

        sample_ids = request.POST.getlist("sample_ids[]")
        check_box_all_checked = request.POST.get("check_box_all_checked", False)
        if check_box_all_checked:
            sample_ids = []
        else:
            sample_ids = [int(sample_id) for sample_id in sample_ids]

        try:
            if len(runs_to_deploy) > 0:
                for sample, leaves_to_deploy in runs_to_deploy.items():
                    if len(sample_ids) > 0:
                        if sample.pk not in sample_ids:
                            continue
                    for leaf in leaves_to_deploy:
                        taskID = process_SGE.set_submit_televir_run(
                            user=request.user,
                            project_pk=project.pk,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                        )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def deploy_ProjectPI_combined_runs(request):
    """
    prepare data for deployment of pathogen identification.
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=int(user_id))

        samples = PIProject_Sample.objects.filter(
            project=project, is_deleted_in_file_system=False
        )

        sample_ids = request.POST.getlist("sample_ids[]")
        check_box_all_checked = request.POST.get("check_box_all_checked", False)
        if check_box_all_checked:
            samples = samples.filter(is_deleted_in_file_system=False)
        elif len(sample_ids) > 0:
            samples = samples.filter(
                pk__in=[int(sample_id) for sample_id in sample_ids]
            )

        try:

            samples_deployed = 0

            for sample in samples:

                reference_utils = RawReferenceUtils(sample)

                count_references = (
                    reference_utils.query_sample_compound_references_regressive()
                )

                if count_references.exists() is False:
                    continue

                software_utils = SoftwareTreeUtils(user, project, sample=sample)
                runs_to_deploy = (
                    software_utils.check_runs_to_submit_metagenomics_sample(sample)
                )
                for sample, leaves_to_deploy in runs_to_deploy.items():

                    reference_manager = SampleReferenceManager(sample)
                    for leaf in leaves_to_deploy:
                        metagenomics_run = reference_manager.mapping_run_from_leaf(leaf)

                        taskID = process_SGE.set_submit_televir_sample_metagenomics(
                            user=request.user,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                            combined_analysis=True,
                            map_run_pk=metagenomics_run.pk,
                        )

                samples_deployed += 1

            if samples_deployed > 0:
                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def submit_televir_project_sample_runs(request):
    """
    submit a new sample to televir project
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = Projects.objects.get(id=int(sample.project.pk))

        software_utils = SoftwareTreeUtils(user, project)
        runs_to_deploy = software_utils.check_runs_to_deploy_sample(sample)

        try:
            if len(runs_to_deploy) > 0:
                for sample, leafs_to_deploy in runs_to_deploy.items():
                    for leaf in leafs_to_deploy:
                        taskID = process_SGE.set_submit_televir_run(
                            user=request.user,
                            project_pk=project.pk,
                            sample_pk=sample.pk,
                            leaf_pk=leaf.pk,
                        )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def submit_televir_project_sample(request):
    """
    submit a new sample to televir project
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}
        process_SGE = ProcessSGE()
        user = request.user

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = Projects.objects.get(id=int(sample.project.pk))

        software_utils = SoftwareTreeUtils(user, project=project)
        runs_to_deploy = software_utils.check_runs_to_deploy_sample(sample)

        try:
            if len(runs_to_deploy) > 0:
                for sample, leafs_to_deploy in runs_to_deploy.items():
                    taskID = process_SGE.set_submit_televir_sample(
                        user=request.user,
                        project_pk=project.pk,
                        sample_pk=sample.pk,
                    )

                data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def Project_explify_merge(request):
    """
    submit a new sample to televir project
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user
        utils: Utils = Utils()
        try:
            temp_directory = utils.get_temp_dir()

            project_id = int(request.POST["project_id"])
            project = Projects.objects.get(id=int(project_id))

            rpip_file = request.FILES["rpip_file"]
            upip_file = request.FILES["upip_file"]

            rpip_report_path = os.path.join(
                temp_directory,
                rpip_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )
            upip_report_path = os.path.join(
                temp_directory,
                upip_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )

            with open(rpip_report_path, "wb") as f:
                f.write(rpip_file.file.read())
            with open(upip_report_path, "wb") as f:
                f.write(upip_file.file.read())

        except Exception as e:
            print(e)
            return JsonResponse(data)

        ### check process not running for this project
        process_controler = ProcessControler()
        try:
            ProcessControler.objects.get(
                owner__id=user.pk,
                name=process_controler.get_name_televir_project_merge_explify(
                    project_pk=project.pk,
                ),
                is_running=True,
            )

        except ProcessControler.DoesNotExist:
            all_reports = FinalReport.objects.filter(
                run__project__pk=int(project.pk)
            ).order_by("-coverage")

            televir_reports = pd.DataFrame.from_records(all_reports.values())

            report_path = os.path.join(
                temp_directory, f"televir_project.{project.pk}.tsv"
            )
            televir_reports.to_csv(report_path, sep="\t", index=False)

            _ = process_SGE.set_submit_televir_explify_merge(
                user=request.user,
                project_pk=project.pk,
                rpip_filepath=rpip_report_path,
                upip_filepath=upip_report_path,
                televir_report_filepath=report_path,
                out_dir=temp_directory,
            )

            data["is_deployed"] = True

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def Project_explify_merge_external(request):
    """
    merge explify rpip and upip reports to televir report, all provided by user.
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user
        utils: Utils = Utils()
        try:
            temp_directory = utils.get_temp_dir()

            rpip_file = request.FILES["rpip_file"]
            upip_file = request.FILES["upip_file"]
            project_file = request.FILES["project_file"]

            rpip_report_path = os.path.join(
                temp_directory,
                rpip_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )
            upip_report_path = os.path.join(
                temp_directory,
                upip_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )

            report_path = os.path.join(
                temp_directory,
                project_file.name.replace(" ", "_").replace("(", "_").replace(")", "_"),
            )

            with open(rpip_report_path, "wb") as f:
                f.write(rpip_file.file.read())
            with open(upip_report_path, "wb") as f:
                f.write(upip_file.file.read())
            with open(report_path, "wb") as f:
                f.write(project_file.file.read())

        except Exception as e:
            print(e)
            return JsonResponse(data)

        ### check process not running for this project
        process_controler = ProcessControler()
        try:
            ProcessControler.objects.get(
                owner__id=user.pk,
                name=process_controler.get_name_televir_project_merge_explify_external(
                    user_pk=user.pk,
                ),
                is_running=True,
            )

        except ProcessControler.DoesNotExist:
            taskID = process_SGE.set_submit_televir_explify_merge_external(
                user=request.user,
                rpip_filepath=rpip_report_path,
                upip_filepath=upip_report_path,
                televir_report_filepath=report_path,
                out_dir=temp_directory,
            )

            data["is_deployed"] = True

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def Update_televir_project(request):
    """
    update televir project
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))

        try:
            _ = process_SGE.set_submit_update_televir_project(
                user=request.user,
                project_id=project.pk,
            )

            data["is_deployed"] = True

        except Exception as e:
            print(e)
            data["is_deployed"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def Project_explify_delete_external(request):
    """
    delete external televir report
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user
        utils: Utils = Utils()
        try:
            output_file_merged = os.path.join(
                get_services_dir(user), "merged_televir_explify.tsv"
            )
            os.remove(output_file_merged)
        except Exception as e:
            print(e)
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def kill_televir_project_sample(request):
    """
    kill all processes a sample, set queued to false
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = Projects.objects.get(id=int(sample.project.pk))

        runs_params = ParameterSet.objects.filter(
            sample=sample,
            status__in=[
                ParameterSet.STATUS_RUNNING,
                ParameterSet.STATUS_QUEUED,
            ],
        )

        for single_run_param in runs_params:
            try:  # kill process

                if single_run_param.status == ParameterSet.STATUS_RUNNING:
                    single_run_param.delete_run_data()

                process_SGE.kill_televir_process_controler_runs(
                    user.pk, project.pk, sample.pk, single_run_param.leaf.pk
                )

                single_run_param.status = ParameterSet.STATUS_KILLED
                single_run_param.save()

            except ProcessControler.DoesNotExist as e:
                print(e)
                print("ProcessControler.DoesNotExist")
                pass

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def kill_televir_project_tree_sample(request):
    """
    kill all processes a sample, set queued to false
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(id=int(sample_id))
        project = Projects.objects.get(id=int(sample.project.pk))

        runs_params = ParameterSet.objects.filter(
            sample=sample,
            status__in=[
                ParameterSet.STATUS_RUNNING,
                ParameterSet.STATUS_QUEUED,
            ],
        )
        killed = 0

        for single_run_param in runs_params:
            try:  # kill process

                if single_run_param.status == ParameterSet.STATUS_RUNNING:
                    single_run_param.delete_run_data()

                process_SGE.kill_televir_process_controler_samples(
                    user.pk, project.pk, sample.pk, single_run_param.leaf.pk
                )
                killed += 1

            except ProcessControler.DoesNotExist as e:
                print(e)
                print("ProcessControler.DoesNotExist")
                pass

            single_run_param.status = ParameterSet.STATUS_KILLED
            single_run_param.save()

        if killed > 0:
            data["is_deployed"] = True

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def kill_televir_project_all_sample(request):
    """
    kill all processes a sample, set queued to false
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False, "is_empty": True}

        process_SGE = ProcessSGE()
        user = request.user

        project_id = int(request.POST["project_id"])
        project = Projects.objects.get(id=int(project_id))
        sample_ids = request.POST.getlist("sample_ids[]")
        check_box_all_checked = request.POST.get("check_box_all_checked", False)
        if check_box_all_checked:
            sample_ids = []
        else:
            sample_ids = [int(sample_id) for sample_id in sample_ids]

        samples = PIProject_Sample.objects.filter(project__id=int(project_id))

        if len(sample_ids) > 0:
            samples = samples.filter(pk__in=sample_ids)

        killed = 0

        for sample in samples:

            runs_params = ParameterSet.objects.filter(
                sample=sample,
                status__in=[
                    ParameterSet.STATUS_RUNNING,
                    ParameterSet.STATUS_QUEUED,
                ],
            )
            for single_run_param in runs_params:
                try:  # kill process

                    if single_run_param.status == ParameterSet.STATUS_RUNNING:
                        single_run_param.delete_run_data()

                    process_SGE.kill_televir_process_controler_samples(
                        user.pk, project.pk, sample.pk, single_run_param.leaf.pk
                    )
                    killed += 1

                except ProcessControler.DoesNotExist as e:
                    print(e)
                    print("ProcessControler.DoesNotExist")
                    pass

                single_run_param.status = ParameterSet.STATUS_KILLED
                single_run_param.save()

        if killed > 0:
            data["is_empty"] = False

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def sort_report_projects(request):
    """
    sort report projects
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}
        process_SGE = ProcessSGE()
        samples = PIProject_Sample.objects.filter(
            project__pk=int(request.POST["project_id"])
        )

        project = Projects.objects.get(id=int(request.POST["project_id"]))
        samples = PIProject_Sample.objects.filter(project=project)
        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=project.pk
        )
        try:
            for sample in samples:
                final_reports = FinalReport.objects.filter(sample=sample)

                report_sorter = ReportSorter(
                    sample, final_reports, report_layout_params
                )

                if report_sorter.reports_availble is False:
                    pass
                elif report_sorter.check_analyzed():
                    pass
                else:
                    taskID = process_SGE.set_submit_televir_sort_pisample_reports(
                        user=request.user,
                        pisample_pk=sample.pk,
                    )
                    data["is_deployed"] = True

        except Exception as e:
            print(e)
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def sort_report_sample(request):
    """
    sort report projects
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}
        process_SGE = ProcessSGE()
        sample = PIProject_Sample.objects.get(pk=int(request.POST["sample_id"]))

        project = sample.project
        report_layout_params = TelevirParameters.get_report_layout_params(
            project_pk=project.pk
        )
        try:
            final_reports = FinalReport.objects.filter(sample=sample)

            report_sorter = ReportSorter(sample, final_reports, report_layout_params)

            if report_sorter.reports_availble is False:
                pass
            elif report_sorter.check_analyzed():
                pass
            else:
                taskID = process_SGE.set_submit_televir_sort_pisample_reports(
                    user=request.user,
                    pisample_pk=sample.pk,
                )
                data["is_deployed"] = True

        except Exception as e:
            print(e)
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def teleflu_igv_create(request):
    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        teleflu_project_pk = int(request.POST["pk"])
        teleflu_project = TeleFluProject.objects.get(pk=teleflu_project_pk)
        insaflu_project = teleflu_project.insaflu_project

        ### get reference
        reference = teleflu_project.reference
        if reference is None:
            return JsonResponse(data)

        reference_file = reference.get_reference_fasta(TypePath.MEDIA_ROOT)

        samples = InsafluProjectSample.objects.filter(project=insaflu_project)
        # samples= [sample.sample for sample in samples]
        sample_dict = {}

        ### get sample files
        software_names = SoftwareNames()

        for sample in samples:
            bam_file = sample.get_file_output(
                TypePath.MEDIA_ROOT, FileType.FILE_BAM, software_names.get_snippy_name()
            )
            bam_file_index = sample.get_file_output(
                TypePath.MEDIA_ROOT,
                FileType.FILE_BAM_BAI,
                software_names.get_snippy_name(),
            )
            vcf_file = sample.get_file_output(
                TypePath.MEDIA_ROOT, FileType.FILE_VCF, software_names.get_snippy_name()
            )

            if bam_file and bam_file_index and vcf_file:
                sample_dict[sample.sample.pk] = {
                    "name": sample.sample.name,
                    "bam_file": bam_file,
                    "bam_file_index": bam_file_index,
                    "vcf_file": vcf_file,
                }

        ### merge vcf files
        televir_bioinf = TelevirBioinf()
        vcf_files = [files["vcf_file"] for sample_pk, files in sample_dict.items()]
        group_vcf = teleflu_project.project_vcf
        stacked_html = teleflu_project.project_igv_report_media

        os.makedirs(teleflu_project.project_vcf_directory, exist_ok=True)

        merged_success = televir_bioinf.merge_vcf_files(vcf_files, group_vcf)

        try:

            televir_bioinf.create_igv_report(
                reference_file,
                vcf_file=group_vcf,
                tracks=sample_dict,
                output_html=stacked_html,
            )

        except Exception as e:
            print(e)
            return JsonResponse(data)

        return JsonResponse(data)


@login_required
@require_POST
def create_insaflu_reference_from_raw(request):
    if request.is_ajax():
        data = {"is_ok": False, "exists": False}

        ref_id = int(request.POST["ref_id"])
        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=user_id)
        process_SGE = ProcessSGE()

        try:
            raw_ref = RawReference.objects.get(id=ref_id)

            description = raw_ref.description
            accid = raw_ref.accid

            if check_user_reference_exists(
                description, accid, user_id
            ) or check_raw_reference_submitted(ref_id=ref_id, user_id=user_id):
                data["is_ok"] = True
                data["exists"] = True
                return JsonResponse(data)
            # success = create_reference(ref_id, user_id)
            taskID = process_SGE.set_submit_raw_televir_teleflu_create(user, ref_id)

        except Exception as e:
            print(e)
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def create_insaflu_reference_from_filemap(request):
    if request.is_ajax():
        data = {"is_ok": False, "exists": False}

        ref_id = int(request.POST["ref_id"])
        user_id = int(request.POST["user_id"])
        user = User.objects.get(id=user_id)
        process_SGE = ProcessSGE()

        try:
            reference = ReferenceSourceFileMap.objects.get(id=ref_id)
            accid = reference.reference_source.accid
            description = reference.description

            if check_user_reference_exists(
                description, accid, user_id
            ) or check_file_reference_submitted(ref_id=ref_id, user_id=user_id):
                data["is_ok"] = True
                data["exists"] = True
                return JsonResponse(data)
            # success = create_reference(ref_id, user_id)
            _ = process_SGE.set_submit_file_televir_teleflu_create(user, ref_id)

        except Exception as e:
            print(e)
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def add_references_to_sample(request):
    """
    add references to sample
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False, "is_empty": False}
        sample_id = int(request.POST["sample_id"])
        sample = PIProject_Sample.objects.get(pk=sample_id)

        reference_id_list = request.POST.getlist("reference_ids[]")

        if len(reference_id_list) == 0:
            data["is_empty"] = True
            return JsonResponse(data)

        reference_id_list = [int(x) for x in reference_id_list]

        sample_reference_manager = SampleReferenceManager(sample)

        try:
            ref_sources = ReferenceSourceFileMap.objects.filter(
                pk__in=reference_id_list
            )

            for reference in ref_sources:

                sample_reference_manager.add_reference(reference)

        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)

        data = {"is_ok": True}
        return JsonResponse(data)


def inject_references(references: list, request):
    context = {}
    data = {}

    context["references_table"] = ReferenceSourceTable(references)

    template_table_html = os.path.join(
        BASE_DIR,
        "templates",
        "pathogen_identification/references_table_table_only.html",
    )
    template_table_html = "pathogen_identification/references_table_table_only.html"

    # render tamplate using context
    rendered_table = render_to_string(template_table_html, context, request=request)
    data["my_content"] = rendered_table
    data["references_count"] = len(references)
    data["is_empty"] = len(references) == 0

    return data


@login_required
@csrf_protect
def create_teleflu_project(request):
    """
    create teleflu project
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False, "exists": False, "is_empty": False}

        ref_ids = request.POST.getlist("ref_ids[]")
        sample_ids = request.POST.getlist("sample_ids[]")
        check_box_all_checked = request.POST.get("check_box_all_checked", False)

        def teleflu_project_name_from_refs(ref_ids):

            refs = [RawReference.objects.get(pk=int(x)) for x in ref_ids]
            date_now_str = datetime.now().strftime("%Y%m%d")
            if len(refs) == 0:
                return "teleflu_project"

            if len(refs) == 1:
                return f"televir_project_{refs[0].accid}_{date_now_str}"

            return f"televir_project_multiple_refs_{date_now_str}"

        def teleflu_project_description(ref_ids):
            if len(ref_ids) == 0:
                return "teleflu_project"
            if len(ref_ids) == 1:
                return "single reference project"
            return "multiple references project"

        project_name = teleflu_project_name_from_refs(ref_ids)

        first_ref = RawReference.objects.get(pk=int(ref_ids[0]))

        project = first_ref.run.project
        if check_box_all_checked:
            sample_ids = PIProject_Sample.objects.filter(project=project).values_list(
                "pk", flat=True
            )
        date = datetime.now()

        try:

            metareference = create_combined_reference(ref_ids, project_name)

            if not metareference:
                data["is_error"] = True
                return JsonResponse(data)

            try:
                teleflu_project = TeleFluProject.objects.get(
                    televir_project=project,
                    name=project_name,
                    raw_reference=metareference,
                    is_deleted=False,
                )

                data["is_ok"] = True
                data["exists"] = True
                data["project_id"] = teleflu_project.pk
                data["project_name"] = teleflu_project.name

                return JsonResponse(data)

            except TeleFluProject.DoesNotExist:

                teleflu_project = TeleFluProject(
                    televir_project=project,
                    name=project_name,
                    last_change_date=date,
                    description=teleflu_project_description(ref_ids),
                    raw_reference=metareference,
                )
                teleflu_project.save()

            for sample_id in sample_ids:
                sample = PIProject_Sample.objects.get(pk=int(sample_id))
                TeleFluSample.objects.create(
                    teleflu_project=teleflu_project,
                    televir_sample=sample,
                )

            data["is_ok"] = True
            data["project_id"] = teleflu_project.pk
            data["project_name"] = teleflu_project.name

        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)

        return JsonResponse(data)


@login_required
@csrf_protect
def query_teleflu_projects(request):
    """
    query teleflu_projects
    """
    if request.is_ajax():
        data = {
            "is_ok": False,
            "is_error": False,
            "is_empty": False,
            "teleflu_data": [],
        }

        televir_project_id = int(request.GET["project_id"])

        ## TeleFlu Projects
        try:
            teleflu_projects = TeleFluProject.objects.filter(
                televir_project=televir_project_id, is_deleted=False
            ).order_by("-last_change_date")

            teleflu_data = []
            for tproj in teleflu_projects:

                tproject_data = {
                    "id": tproj.pk,
                    "samples": tproj.nsamples,
                    "ref_description": tproj.raw_reference.description_first,
                    "ref_accid": tproj.raw_reference.accids_str,
                    "ref_taxid": tproj.raw_reference.taxids_str,
                    "insaflu_project": False if tproj.insaflu_project is None else True,
                }

                insaflu_project = tproj.insaflu_project
                if insaflu_project is None:
                    tproject_data["insaflu_project"] = "None"
                else:
                    count_not_finished = InsafluProjectSample.objects.filter(
                        project__id=insaflu_project.id,
                        is_deleted=False,
                        is_error=False,
                        is_finished=False,
                    ).count()

                    if count_not_finished == 0:
                        tproject_data["insaflu_project"] = "Finished"

                    else:
                        tproject_data["insaflu_project"] = "Processing"

                teleflu_data.append(tproject_data)

        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)

        if len(teleflu_data) == 0:
            data["is_empty"] = True
            return JsonResponse(data)

        data["teleflu_projects"] = teleflu_data
        data["is_ok"] = True

        return JsonResponse(data)


@login_required
@csrf_protect
def delete_teleflu_project(request):
    """
    delete teleflu project
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False}

        try:
            teleflu_project_id = int(request.POST["project_id"])
            teleflu_project = TeleFluProject.objects.get(pk=teleflu_project_id)

            teleflu_project.is_deleted = True
            teleflu_project.save()

            data["is_ok"] = True

        except Exception as e:
            print(e)
            data["is_error"] = True

        return JsonResponse(data)


@login_required
@csrf_protect
def create_insaflu_project(request):
    """
    create insaflu project associated with teleflu map project"""
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False, "exists": False}
        teleflu_project_id = int(request.POST["project_id"])
        teleflu_project = TeleFluProject.objects.get(pk=teleflu_project_id)
        process_SGE = ProcessSGE()

        if teleflu_project.insaflu_project is not None:
            data["exists"] = True
            return JsonResponse(data)

        try:
            taskID = process_SGE.set_submit_televir_teleflu_project_create(
                user=request.user,
                project_pk=teleflu_project.pk,
            )
            data["is_ok"] = True
            return JsonResponse(data)
        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)


@csrf_protect
def set_teleflu_check_box_values(request):
    """
    manage check boxes through ajax
    """
    if request.is_ajax():
        data = {"is_ok": False}
        utils = Utils()
        if Constants.GET_CHECK_BOX_SINGLE in request.GET:
            data["is_ok"] = True
            for key in request.session.keys():
                if (
                    key.startswith(Constants.TELEFLU_CHECK_BOX)
                    and len(key.split("_")) == 4
                    and utils.is_integer(key.split("_")[3])
                ):
                    data[key] = request.session[key]
        ## change single status of a check_box_single
        elif Constants.GET_CHANGE_CHECK_BOX_SINGLE in request.GET:
            data["is_ok"] = True
            key_name = "{}_{}".format(
                Constants.TELEFLU_CHECK_BOX, request.GET.get(Constants.CHECK_BOX_VALUE)
            )
            for key in request.session.keys():
                if (
                    key.startswith(Constants.TELEFLU_CHECK_BOX)
                    and len(key.split("_")) == 4
                    and utils.is_integer(key.split("_")[3])
                ):
                    if request.session[key]:
                        data[key] = False
                    if key == key_name:
                        request.session[key] = utils.str2bool(
                            request.GET.get(Constants.GET_CHANGE_CHECK_BOX_SINGLE)
                        )
                    else:
                        request.session[key] = False

        return JsonResponse(data)


@login_required
@require_POST
def add_teleflu_sample(request):
    """add samples to teleflu_project"""

    if request.is_ajax():
        data = {
            "is_ok": False,
            "is_error": False,
            "is_empty": False,
            "not_added": False,
        }

        sample_ids = request.POST.getlist("sample_ids[]")
        ref_id = int(request.POST["teleflu_id"])
        check_box_all_checked = request.POST.get("check_box_all_checked", False)

        teleflu_project = TeleFluProject.objects.get(pk=ref_id)
        if check_box_all_checked:
            sample_ids = PIProject_Sample.objects.filter(
                project=teleflu_project.televir_project
            ).values_list("pk", flat=True)

        if len(sample_ids) == 0:
            data["is_empty"] = True
            return JsonResponse(data)

        added = 0

        for sample_id in sample_ids:

            sample = PIProject_Sample.objects.get(pk=int(sample_id))

            try:
                TeleFluSample.objects.get(
                    teleflu_project=teleflu_project,
                    televir_sample=sample,
                )
            except TeleFluSample.DoesNotExist:

                TeleFluSample.objects.create(
                    teleflu_project=teleflu_project,
                    televir_sample=sample,
                )

            added += 1

        if added == 0:
            data["not_added"] = True
            return JsonResponse(data)
        else:
            data["added"] = True

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@require_POST
def add_teleflu_mapping_workflow(request):
    """
    create mapping workflow for teleflu project
    """

    if request.is_ajax():
        data = {"is_ok": False}
        project_id = int(request.POST["project_id"])
        leaf_id = int(request.POST["leaf_id"])

        project = TeleFluProject.objects.get(pk=project_id)
        leaf = SoftwareTreeNode.objects.get(pk=leaf_id)

        try:
            TelefluMapping.objects.get(
                teleflu_project=project,
                leaf=leaf,
            )

            data["exists"] = True

        except TelefluMapping.DoesNotExist:
            TelefluMapping.objects.create(
                teleflu_project=project,
                leaf=leaf,
            )
            data["is_ok"] = True

        except Exception as e:
            print(e)
            data["is_ok"] = False

        return JsonResponse(data)


############
# TELEFLU


def excise_paths_leaf_last(string_with_paths: str):
    """
    if the string has paths, find / and return the last
    """
    split_space = string_with_paths.split(" ")
    new_string = ""
    if len(split_space) > 1:
        for word in split_space:
            new_word = word
            if "/" in word:
                new_word = word.split("/")[-1]
            new_string += new_word + " "
        return new_string
    else:
        return string_with_paths


def teleflu_node_info(node, params_df, node_pk):
    node_info = {
        "pk": node_pk,
        "node": node,
        "modules": [],
    }

    for pipeline_step in CS.vect_pipeline_televir_workflows_display:
        acronym = [x[0] for x in pipeline_step.split(" ")]
        acronym = "".join(acronym).upper()
        params = params_df[params_df.module == pipeline_step].to_dict("records")
        if params:  # if there are parameters for this module
            software = params[0].get("software")
            software = software.split("_")[0]

            params = params[0].get("value")
            params = excise_paths_leaf_last(params)
            params = f"{software} {params}"
        else:
            params = ""

        node_info["modules"].append(
            {
                "module": pipeline_step,
                "parameters": params,
                "short_name": acronym,
                "available": (
                    "software_on"
                    if pipeline_step in params_df.module.values
                    else "software_off"
                ),
            }
        )

    return node_info


@login_required
def load_teleflu_workflows(request):
    """
    load teleflu workflows
    """
    data = {"is_ok": False, "mapping_workflows": []}

    if request.is_ajax():

        teleflu_project_pk = int(request.GET["project_id"])

        utils_manager = Utils_Manager()

        mappings = TelefluMapping.objects.filter(teleflu_project__pk=teleflu_project_pk)
        teleflu_project = TeleFluProject.objects.get(pk=teleflu_project_pk)
        mapping_workflows = []
        existing_mapping_pks = []

        for mapping in mappings:
            if mapping.leaf is None:
                continue

            params_df = utils_manager.get_leaf_parameters(mapping.leaf)
            node_info = node_info = teleflu_node_info(
                mapping.leaf.index, params_df, mapping.leaf.pk
            )

            samples_mapped = mapping.mapped_samples

            samples_stacked = mapping.stacked_samples_televir
            node_info["running_or_queued"] = mapping.queued_or_running_mappings_exist
            node_info["pk"] = mapping.pk
            node_info["samples_stacked"] = samples_stacked.count()
            node_info["samples_to_stack"] = samples_mapped.exclude(
                pk__in=samples_stacked.values_list("pk", flat=True)
            ).exists()

            sample_summary, mapped_samples, mapped_success = mapping.sample_summary

            mapped_fail = mapped_samples - mapped_success

            node_info["samples_mapped"] = samples_mapped.count()
            node_info["mapped_success"] = mapped_success
            node_info["mapped_fail"] = mapped_fail
            node_info["left_to_map"] = (
                teleflu_project.nsamples - mapped_success - mapped_fail
            ) > 0

            node_info["sample_summary"] = sample_summary
            existing_mapping_pks.append(mapping.leaf.pk)
            node_info["stacked_html_exists"] = os.path.exists(
                mapping.mapping_igv_report
            )
            node_info["stacked_html"] = mapping.mapping_igv_report.replace(
                "/insaflu_web/INSaFLU", ""
            )

            node_info["stacked_variants_vcf"] = mapping.variants_vcf_media_path

            mapping_workflows.append(node_info)

        data["mapping_workflows"] = mapping_workflows
        data["is_ok"] = True
        data["teleflu_project_pk"] = teleflu_project_pk
        data["project_nsamples"] = teleflu_project.nsamples

        return JsonResponse(data)

    return JsonResponse(data)


@login_required
@require_POST
def map_teleflu_workflow_samples(request):

    if request.is_ajax():

        data = {"is_ok": False, "is_error": False, "is_empty": False}
        project_id = int(request.POST["project_id"])
        workflow_id = int(request.POST["workflow_id"])

        teleflu_project = TeleFluProject.objects.get(pk=project_id)
        teleflu_samples = TeleFluSample.objects.filter(teleflu_project=teleflu_project)
        teleflu_mapping = TelefluMapping.objects.get(pk=workflow_id)
        workflow_leaf = teleflu_mapping.leaf
        user = request.user

        mapping = TelefluMapping.objects.get(
            leaf=workflow_leaf, teleflu_project=teleflu_project
        )

        samples_to_map = teleflu_samples.exclude(
            televir_sample__in=mapping.mapped_samples
        )

        references_to_map = teleflu_project.raw_reference.references
        process_SGE = ProcessSGE()

        if len(samples_to_map) == 0:
            data["is_empty"] = True
            return JsonResponse(data)

        deployed = 0
        for sample in samples_to_map:
            reference_manager = SampleReferenceManager(sample.televir_sample)
            mapping_run = reference_manager.mapping_request_run_from_leaf(workflow_leaf)
            for reference in references_to_map:
                reference.pk = None
                reference.run = mapping_run
                reference.save()

            taskID = process_SGE.set_submit_televir_sample_metagenomics(
                user=user,
                sample_pk=sample.televir_sample.pk,
                leaf_pk=workflow_leaf.pk,
                mapping_request=True,
                map_run_pk=mapping_run.pk,
            )
            deployed += 1

        if deployed == 0:
            data["is_error"] = True
            return JsonResponse(data)

        data["is_ok"] = True
        return JsonResponse(data)


@login_required
@csrf_protect
def stack_igv_teleflu_workflow(request):
    """
    create insaflu project associated with teleflu map project"""
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False, "exists": False, "running": False}
        teleflu_project_id = int(request.POST["project_id"])
        mapping_id = int(request.POST["workflow_id"])

        try:
            teleflu_mapping = TelefluMapping.objects.get(
                pk=mapping_id,
            )

        except TelefluMapping.DoesNotExist:
            data["is_error"] = True
            return JsonResponse(data)

        process_SGE = ProcessSGE()
        process_controler = ProcessControler()

        if ProcessControler.objects.filter(
            owner=request.user,
            name=process_controler.get_name_televir_teleflu_igv_stack(
                teleflu_mapping_id=mapping_id,
            ),
            is_running=True,
        ).exists():
            data["running"] = True

            return JsonResponse(data)

        try:
            process_SGE.set_submit_teleflu_map(
                user=request.user,
                leaf_pk=teleflu_mapping.leaf.pk,
                project_pk=teleflu_project_id,
            )
            data["is_ok"] = True

        except Exception as e:
            print(e)
            data["is_error"] = True

        return JsonResponse(data)


@login_required
@require_POST
def add_references_all_samples(request):
    """
    add references to sample
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False, "is_empty": False}
        project_id = int(request.POST["ref_id"])
        project = Projects.objects.get(pk=project_id)
        samples = PIProject_Sample.objects.filter(project=project)

        reference_id_list = request.POST.getlist("reference_ids[]")
        data["empty_content"] = inject_references([], request)["my_content"]
        if len(reference_id_list) == 0:
            data["is_empty"] = True
            return JsonResponse(data)

        reference_id_list = [int(x) for x in reference_id_list]

        references_existing = []

        try:
            ref_sources = ReferenceSourceFileMap.objects.filter(
                pk__in=reference_id_list
            )
            for sample in samples:
                sample_reference_manager = SampleReferenceManager(sample)
                sample_id = sample.pk

                for reference in ref_sources:
                    if RawReference.objects.filter(
                        accid=reference.accid,
                        run__sample__pk=sample_id,
                        run__run_type=RunMain.RUN_TYPE_STORAGE,
                    ).exists():
                        references_existing.append(reference.accid)
                        continue

                    # truncate reference description: max 150 characters
                    sample_reference_manager.add_reference(reference)

        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)

        data["is_ok"] = True

        return JsonResponse(data)


@login_required
@require_POST
def remove_added_reference(request):
    """
    remove added reference
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_error": False}

        reference_id = int(request.POST["reference_id"])
        sample_id = int(request.POST["sample_id"])

        try:
            reference = RawReference.objects.get(pk=reference_id)
            reference.delete()
        except Exception as e:
            print(e)
            data["is_error"] = True
            return JsonResponse(data)

        query_set_added_manual = RawReference.objects.filter(
            run__sample__pk=sample_id, run__run_type=RunMain.RUN_TYPE_STORAGE
        )

        context = inject__added_references(query_set_added_manual, request)
        data["added_references"] = context["my_content"]
        data["is_empty"] = query_set_added_manual.count() == 0

        data = {"is_ok": True}
        return JsonResponse(data)

    return JsonResponse({"is_ok": False})


@login_required
@require_POST
def deploy_televir_map(request):
    """
    prepare data for deployment of pathogen identification.
    """

    if request.is_ajax():
        data = {"is_ok": False, "is_deployed": False}

        process_SGE = ProcessSGE()
        user = request.user

        reference_id = int(request.POST["reference_id"])
        project_id = int(request.POST["project_id"])
        taskID = process_SGE.set_submit_televir_map(
            user, reference_pk=reference_id, project_pk=project_id
        )

        data["is_ok"] = True

        return JsonResponse(data)


#######################################################################
####################################################################### PANELS


def check_metadata_table_clean(metadata_table_file) -> Optional[pd.DataFrame]:
    """
    check metadata table
    """

    metadata_table = metadata_table_file.read().decode("utf-8")
    ### determine line terminator
    lineterminator = "\n"
    if "\r\n" in metadata_table:
        lineterminator = "\r\n"

    metadata_table = metadata_table.split(lineterminator)
    file_extention = os.path.splitext(metadata_table_file.name)[1]

    sep = "\t" if file_extention == ".tsv" else ","
    # check for the line terminator

    metadata_table = [x.split(sep) for x in metadata_table]
    metadata_table = [x for x in metadata_table if len(x) > 1]

    try:
        metadata_table = pd.DataFrame(metadata_table[1:], columns=metadata_table[0])

    except Exception as e:
        print(e)
        return None

    if len(metadata_table) == 0:
        return None
    return metadata_table


def check_table_columns(metadata_table):

    if "Accession ID" not in metadata_table.columns:
        return False

    if "TaxID" not in metadata_table.columns:
        return False

    return True


@login_required
@require_POST
def check_panel_upload_clean(request):
    """
    check if fasta and metadata coherent.
    """

    if request.is_ajax():

        data = {
            "is_ok": False,
            "is_error": False,
            "error_message": "",
            "summary": "",
            "success": "",
            "log": [],
            "pass": False,
        }

        description = request.POST.get("description", "").strip()

        if request.FILES.get("metadata", None) is None:
            data["is_error"] = True
            data["error_message"] = "Metadata file is empty."
            return JsonResponse(data)

        if request.FILES.get("fasta_file", None) is None:
            data["is_error"] = True
            data["error_message"] = "Fasta file is empty."
            return JsonResponse(data)

        data[Constants.SEQUENCES_TO_PASS] = []
        software = Software()
        utils = Utils()

        # name = request.POST.get("name", "").strip()
        reference_metadata_table_file = request.FILES.get("metadata", None)
        reference_fasta_file = request.FILES.get("fasta_file", None)

        try:
            ReferenceSourceFile.objects.get(
                file=reference_fasta_file.name, owner=request.user, is_deleted=False
            )

            data["is_error"] = True
            data["error_message"] = "Fasta file already exists."
            return JsonResponse(data)
        except ReferenceSourceFile.DoesNotExist:
            pass
        except ReferenceSourceFile.MultipleObjectsReturned:
            data["is_error"] = True
            data["error_message"] = "Fasta file already exists."
            return JsonResponse(data)

        if os.path.splitext(reference_metadata_table_file.name)[1] not in [
            ".tsv",
            ".csv",
        ]:
            data["is_error"] = True
            data["error_message"] = "Metadata file must be in tsv or csv format."
            return JsonResponse(data)

        if os.path.splitext(reference_fasta_file.name)[1] not in [".fasta", ".fa"]:
            data["is_error"] = True
            data["error_message"] = (
                "Fasta file must be in fasta (.fa or .fasta) format."
            )
            return JsonResponse(data)

        reference_metadata_table = pd.DataFrame(
            columns=["accid", "taxid", "description"]
        )

        if reference_metadata_table is None:
            data["is_error"] = True
            data["error_message"] = "Metadata table is empty."
            return JsonResponse(data)

        reference_metadata_table = check_metadata_table_clean(
            reference_metadata_table_file
        )

        if check_table_columns(reference_metadata_table) is False:
            data["is_error"] = True
            data["error_message"] = "Metadata table has not the right columns."
            return JsonResponse(data)

        error_message = ""
        data["error_message"] = ""

        if "Description" not in reference_metadata_table.columns:
            description_ref = f"Panel reference"
            if description != "":
                description_ref = description_ref + f" - {description}"
            reference_metadata_table["Description"] = description_ref
            error_message = "Description column not found. Added default description."

        ### testing file names
        ## testing fasta
        some_error_in_files = False
        reference_fasta_temp_file_name = NamedTemporaryFile(
            prefix="flu_fa_", delete=False
        )
        try:
            file_data = reference_fasta_file.read()
            reference_fasta_temp_file_name.write(file_data)
            reference_fasta_temp_file_name.flush()
            reference_fasta_temp_file_name.close()
            software.dos_2_unix(reference_fasta_temp_file_name.name)
            data["temp_file"] = reference_fasta_temp_file_name.name
        except Exception as e:
            some_error_in_files = True
            error_message = "Error in the fasta file"
            data["is_error"] = True
            return JsonResponse(data)

        try:
            number_locus = utils.is_fasta(reference_fasta_temp_file_name.name)
            data["log"].append("Number of sequences in fasta: {}".format(number_locus))

            ## test the max numbers
            if number_locus > Constants.MAX_SEQUENCES_FROM_CONTIGS_FASTA:
                error_message = "Max allow number of contigs in Multi-Fasta: {}".format(
                    Constants.MAX_SEQUENCES_FROM_CONTIGS_FASTA
                )
                some_error_in_files = True

            total_length_fasta = utils.get_total_length_fasta(
                reference_fasta_temp_file_name.name
            )

            data["log"].append(
                "\nTotal length of the sequences in fasta: {}".format(
                    total_length_fasta
                )
            )

            if (
                not some_error_in_files
                and total_length_fasta > PICS.MAX_LENGTH_SEQUENCE_TOTAL_REFERENCE_FASTA
            ):
                some_error_in_files = True
                error_message = (
                    "The max sum length of the sequences in fasta: {}".format(
                        PICS.MAX_LENGTH_SEQUENCE_TOTAL_REFERENCE_FASTA
                    )
                )
                data["is_error"] = True
                data["error_message"] = error_message
                return JsonResponse(data)

            n_seq_name_bigger_than = utils.get_number_seqs_names_bigger_than(
                reference_fasta_temp_file_name.name,
                Constants.MAX_LENGTH_CONTIGS_SEQ_NAME,
                0,
            )

            if not some_error_in_files and n_seq_name_bigger_than > 0:
                some_error_in_files = True
                if n_seq_name_bigger_than == 1:
                    error_message = "There is one sequence name length bigger than {0}. The max. length name is {0}.".format(
                        Constants.MAX_LENGTH_CONTIGS_SEQ_NAME
                    )
                else:
                    error_message = "There are {0} sequences with name length bigger than {1}. The max. length name is {1}.".format(
                        n_seq_name_bigger_than,
                        Constants.MAX_LENGTH_CONTIGS_SEQ_NAME,
                    )

            ## if some errors in the files, fasta or genBank, return
            if some_error_in_files:
                data["is_error"] = True
                data["error_message"] = error_message
                return data

            ### check if there all seq names are present in the database yet
            b_pass = False
            vect_error, vect_fail_seqs, vect_pass_seqs = [], [], []
            dict_names = {}
            retained_rows = []
            already_existing_ids = []
            with open(reference_fasta_temp_file_name.name) as handle_in:
                for record in SeqIO.parse(handle_in, "fasta"):

                    if record.id in reference_metadata_table["Accession ID"].values:
                        dict_names[record.id] = 1
                    if (
                        len(reference_metadata_table)
                    ) > 0 and not record.id in reference_metadata_table[
                        "Accession ID"
                    ].values:
                        vect_fail_seqs.append(record.id)
                        vect_error.append(
                            f"Sequence name '{record.id}' does not have match in the metadata table."
                        )
                        continue

                    ## try to upload
                    referrence_exists = ReferenceSourceFileMap.objects.filter(
                        reference_source__accid__iexact=record.id,
                        reference_source_file__owner=request.user,
                    ).exists()
                    if referrence_exists:
                        already_existing_ids.append(record.id)

                    record_metadata = reference_metadata_table[
                        reference_metadata_table["Accession ID"] == record.id
                    ]

                    retained_rows = retained_rows + [record_metadata]
                    vect_pass_seqs.append(record.id)
                    b_pass = True

            ## if none of them pass throw an error
            error_message = ""
            data["log"].append(
                (
                    f"\n {len(vect_pass_seqs)} sequence(s) name(s) found in metadata table."
                )
            )
            if len(vect_pass_seqs) > 0:
                data["pass"] = True
                retained_rows = pd.concat(retained_rows)
                taxid_count = retained_rows["TaxID"].value_counts()
                accid_count = retained_rows["Accession ID"].value_counts()

                data["log"].append(f"\n {len(taxid_count)} taxid(s).")
                data["log"].append(f"\n {len(accid_count)} accid(s).")

            if not b_pass:
                some_error_in_files = True
                if len(reference_metadata_table) > 0 and len(vect_pass_seqs) == 0:
                    error_message = (
                        "None of these names: '{}' match to the sequences names".format(
                            "', '".join(reference_metadata_table["Accession ID"].values)
                        )
                    )
                    data["error_message"] = error_message
                    data["is_error"] = True
                    return JsonResponse(data)

                for message in vect_error:
                    error_message += message + " "

                    data["log"].append(f"\n{error_message}")
                    data["is_error"] = True
            else:
                ## if empty load all
                data[Constants.SEQUENCES_TO_PASS] = vect_pass_seqs

            if len(already_existing_ids) > 0:
                data["log"].append(
                    f"\n{len(already_existing_ids)} sequences already exist in the database."
                )
                data["log"].append(f"\n{', '.join(already_existing_ids)}")
            ### some sequences names suggested are not present in the file
            vect_fail_seqs = [key for key in dict_names if dict_names[key] == 0]
            if len(vect_fail_seqs) > 0:
                error_message += (
                    "Sequences names '{}' does not have match in the file".format(
                        ", '".join(vect_fail_seqs)
                    )
                )
        except Exception as e:
            print(e)
            some_error_in_files = True
            error_message = "Error in parsing fasta file"
            data["error_message"] = error_message
            data["is_error"] = True
            return JsonResponse(data)

        data["is_ok"] = True
        ## remove temp files
        os.unlink(reference_fasta_temp_file_name.name)
        return JsonResponse(data)


from pathogen_identification.models import ReferenceSourceFile


@login_required
@require_POST
def delete_reference_file(request):
    """
    delete reference panel
    """
    if request.is_ajax():
        data = {"is_ok": False, "is_error": False, "message": ""}
        panel_id = int(request.POST["file_id"])
        try:
            panel = ReferenceSourceFile.objects.get(pk=panel_id)
            panel.is_deleted = True
            panel.save()
            data["is_ok"] = True
            return JsonResponse(data)
        except Exception as e:
            print(e)
            data["is_ok"] = False
            data["message"] = "Error deleting the file"
            return JsonResponse(data)


@login_required
@require_POST
def set_sample_reports_control(request):
    """
    set sample reports control
    """
    if request.is_ajax():
        data = {"is_ok": False}
        data["set_control"] = False
        sample_id = int(request.POST["sample_id"])
        try:
            sample = PIProject_Sample.objects.get(pk=int(sample_id))

            sample_control_flag = False if sample.is_control else True

            sample_reports = FinalReport.objects.filter(sample=sample)
            for report in sample_reports:
                report.control_flag = FinalReport.CONTROL_FLAG_NONE

                report.save()

            sample.is_control = sample_control_flag
            sample.save()

            set_control_reports(sample.project.pk)

            data["is_ok"] = True
            data["set_control"] = sample_control_flag
            return JsonResponse(data)

        except Exception as e:
            print(e)
            data["is_ok"] = False
            return JsonResponse(data)


@csrf_protect
def create_reference_panel(request):
    """
    create a reference panel"""
    if request.is_ajax():
        user = request.user
        name = request.POST.get("name")
        icon = request.POST.get("icon", "")

        try:
            ReferencePanel.objects.get(
                name=name,
                owner=user,
                is_deleted=False,
                icon=icon,
                panel_type=ReferencePanel.PANEL_TYPE_MAIN,
            )
        except ReferencePanel.DoesNotExist:
            ReferencePanel.objects.create(
                name=name,
                owner=user,
                icon=icon,
            )

        return JsonResponse({"is_ok": True})

    return JsonResponse({"is_ok": False})


@csrf_protect
def add_references_to_panel(request):
    """
    add references to panel"""
    if request.is_ajax():
        panel_id = int(request.POST.get("ref_id"))
        panel = ReferencePanel.objects.get(pk=panel_id)

        reference_ids = request.POST.getlist("reference_ids[]")
        already_added = []
        for reference_id in reference_ids:
            if panel.references_count >= PICS.MAX_REFERENCES_PANEL:
                return JsonResponse(
                    {"is_full": True, "max_references": PICS.MAX_REFERENCES_PANEL}
                )
            try:
                reference = ReferenceSourceFileMap.objects.get(pk=int(reference_id))
                description = reference.description
                if RawReference.objects.filter(
                    accid=reference.reference_source.accid, panel=panel
                ).exists():
                    already_added.append(reference.reference_source.accid)
                    continue

                if description is None:
                    description = ""

                if len(description) > 200:
                    description = description[:200]

                _ = RawReference.objects.create(
                    accid=reference.reference_source.accid,
                    taxid=reference.reference_source.taxid,
                    description=description,
                    panel=panel,
                )
            except Exception as e:
                print(e)
                return JsonResponse({"is_ok": False})

        return JsonResponse({"is_ok": True})

    return JsonResponse({"is_ok": False})


@csrf_protect
def add_file_to_panel(request):
    """
    add references to panel"""
    if request.is_ajax():
        panel_id = int(request.POST.get("panel_id"))
        file_id = int(request.POST.get("file_id"))

        panel = ReferencePanel.objects.get(pk=panel_id)
        file = ReferenceSourceFile.objects.get(pk=file_id)
        refs = ReferenceSourceFileMap.objects.filter(
            reference_source_file=file
        ).distinct()

        try:
            for reference in refs:
                if panel.references_count >= PICS.MAX_REFERENCES_PANEL:
                    return JsonResponse(
                        {
                            "is_full": True,
                            "max_references": PICS.MAX_REFERENCES_PANEL,
                            "is_ok": True,
                        }
                    )

                if RawReference.objects.filter(
                    accid=reference.reference_source.accid, panel=panel
                ).exists():
                    continue

                RawReference.objects.create(
                    accid=reference.reference_source.accid,
                    taxid=reference.reference_source.taxid,
                    description=reference.reference_source.description,
                    panel=panel,
                )
        except Exception as e:
            print(e)
            return JsonResponse({"is_ok": False})

        return JsonResponse({"is_ok": True})

    return JsonResponse({"is_ok": False})


@csrf_protect
def get_panels(request):
    if request.is_ajax():
        user = request.user
        panels = ReferencePanel.objects.filter(
            owner=user,
            is_deleted=False,
            project_sample=None,
            panel_type=ReferencePanel.PANEL_TYPE_MAIN,
        ).order_by("-creation_date")
        panel_data = [
            {
                "id": panel.pk,
                "name": panel.name,
                "references_count": panel.references_count,
                "icon": panel.icon,
            }
            for panel in panels
        ]

        data = {
            "is_ok": True,
            "panels": panel_data,
        }

        return JsonResponse(data)

    return JsonResponse({"is_ok": False})


@csrf_protect
def remove_panel_reference(request):
    if request.is_ajax():
        panel_id = int(request.POST.get("panel_id"))
        reference_id = int(request.POST.get("reference_id"))

        reference = RawReference.objects.get(pk=reference_id)

        reference.panel = None
        reference.save()

        return JsonResponse({"is_ok": True})

    return JsonResponse({"is_ok": False})


@csrf_protect
def delete_reference_panel(request):
    """
    delete a panel"""

    if request.is_ajax():

        panel_id = int(request.POST.get("panel_id"))
        panel = ReferencePanel.objects.get(pk=panel_id)
        panel.is_deleted = True
        panel.save()

        return JsonResponse({"is_ok": True})

    return JsonResponse({"is_ok": False})


@csrf_protect
def get_panel_references(request):
    """
    get panel references"""
    if request.is_ajax():
        user = request.user
        name = request.GET.get("name")
        panel_id = request.GET.get("panel_id")

        panel = ReferencePanel.objects.get(pk=panel_id)

        references = RawReference.objects.filter(panel=panel)

        data = {
            "is_ok": True,
            "references": list(
                references.values("id", "taxid", "accid", "description")
            ),
            "panel_name": panel.name,
        }

        return JsonResponse(data)

    return JsonResponse({"is_ok": False})


@csrf_protect
def add_panels_to_sample(request):
    """
    add panels to sample"""
    if request.is_ajax():
        sample_id = int(request.POST.get("sample_id"))
        panel_ids = request.POST.getlist("panel_ids[]")

        sample = PIProject_Sample.objects.get(pk=sample_id)

        for panel_id in panel_ids:
            sample.add_panel(int(panel_id))

        return JsonResponse({"is_ok": True})

    return JsonResponse({"is_ok": False})


@csrf_protect
def add_panels_to_project(request):
    """
    add panels to sample"""
    if request.is_ajax():
        project_id = int(request.POST.get("project_id"))
        panel_ids = request.POST.getlist("panel_ids[]")

        samples = PIProject_Sample.objects.filter(project__pk=project_id)

        for panel_id in panel_ids:
            for sample in samples:
                sample.add_panel(int(panel_id))

        return JsonResponse({"is_ok": True})

    return JsonResponse({"is_ok": False})


@csrf_protect
def remove_sample_panel(request):
    """
    remove sample panel"""

    if request.is_ajax():
        panel_id = int(request.POST.get("panel_id"))
        sample_id = int(request.POST.get("sample_id"))
        sample = PIProject_Sample.objects.get(pk=sample_id)
        sample.remove_panel(panel_id)

        return JsonResponse({"is_ok": True})

    return JsonResponse({"is_ok": False})


@csrf_protect
def get_sample_panels(request):
    """
    get sample panels"""
    if request.is_ajax():
        sample_id = request.GET.get("sample_id")

        sample = PIProject_Sample.objects.get(pk=sample_id)

        panels = sample.panels_added

        panel_data = [
            {
                "id": panel.pk,
                "name": panel.name,
                "references_count": panel.references_count,
                "icon": panel.icon,
            }
            for panel in panels
        ]
        data = {"is_ok": True, "panels": panel_data}

        return JsonResponse(data)

    return JsonResponse({"is_ok": False})


def get_sample_panel_suggestions(request):
    """
    get sample panel updates"""
    if request.is_ajax():
        user = request.user
        sample_id = request.GET.get("sample_id")

        sample = PIProject_Sample.objects.get(pk=sample_id)

        panels_sample = sample.panels_added

        panels_global_names = panels_sample.values_list("pk", flat=True)
        panels_suggest = ReferencePanel.objects.filter(
            project_sample=None,
            is_deleted=False,
            panel_type=ReferencePanel.PANEL_TYPE_MAIN,
            owner__in=[user, None],
        ).exclude(pk__in=panels_global_names)

        panel_data = [
            {
                "id": panel.id,
                "name": panel.name,
                "references_count": panel.references_count,
                "icon": panel.icon,
            }
            for panel in panels_suggest
        ]
        data = {"is_ok": True, "panels": panel_data}

        return JsonResponse(data)

    return JsonResponse({"is_ok": False})


def get_project_panel_suggestions(request):
    """
    get sample panel updates"""
    if request.is_ajax():

        owner = request.user

        panels_suggest = ReferencePanel.objects.filter(
            project_sample=None,
            is_deleted=False,
            panel_type=ReferencePanel.PANEL_TYPE_MAIN,
            owner__in=[owner, None],
        )

        panel_data = [
            {
                "id": panel.id,
                "name": panel.name,
                "references_count": panel.references_count,
                "icon": panel.icon,
            }
            for panel in panels_suggest
        ]
        data = {"is_ok": True, "panels": panel_data}

        return JsonResponse(data)

    return JsonResponse({"is_ok": False})


@csrf_protect
def validate_project_name(request):
    """
    test if exist this project name
    """
    if request.is_ajax():
        project_name = request.GET.get("project_name")

        data = {
            "is_taken": Projects.objects.filter(
                name__iexact=project_name,
                is_deleted=False,
                owner__username=request.user.username,
            ).exists()
        }

        ## check if name has spaces:
        if " " in project_name:
            data["has_spaces"] = True
            data["error_message"] = _("Spaces are not allowed in the project name.")

        ## check if name has special characters:
        if not project_name.replace("_", "").isalnum():
            data["has_special_characters"] = True
            data["error_message"] = _(
                "Special characters are not allowed in the project name."
            )

        if data["is_taken"]:
            data["error_message"] = _("Exists a project with this name.")

        return JsonResponse(data)


@csrf_protect
def IGV_display(request):
    """display python plotly app"""

    if request.is_ajax():
        data = {"is_ok": False}
        if request.method == "GET":
            sample_pk = request.GET.get("sample_pk")
            run_pk = request.GET.get("run_pk")
            reference = request.GET.get("accid")
            unique_id = request.GET.get("unique_id")

            sample = PIProject_Sample.objects.get(pk=int(sample_pk))
            sample_name = sample.name
            run = RunMain.objects.get(pk=int(run_pk))

            ref_map = ReferenceMap_Main.objects.get(
                reference=unique_id, sample=sample, run=run
            )

            def remove_pre_static(path: str) -> str:
                cwd = os.getcwd()
                if path.startswith(cwd):
                    path = path[len(cwd) :]

                path = path.replace(STATIC_ROOT, STATIC_URL)

                return path

            #################################################
            ### bam file
            path_name_bam = remove_pre_static(
                ref_map.bam_file_path,
            )
            path_name_bai = remove_pre_static(
                ref_map.bai_file_path,
            )
            path_name_reference = remove_pre_static(
                ref_map.fasta_file_path,
            )
            path_name_reference_index = remove_pre_static(
                ref_map.fai_file_path,
            )
            path_name_vcf = remove_pre_static(ref_map.vcf)

            data["is_ok"] = True
            data["path_bam"] = mark_safe(request.build_absolute_uri(path_name_bam))

            data["path_reference"] = mark_safe(
                request.build_absolute_uri(path_name_reference)
            )
            data["path_reference_index"] = mark_safe(
                request.build_absolute_uri(path_name_reference_index)
            )
            data["reference_name"] = reference

            #### other files
            data["bam_file_id"] = mark_safe(
                '<strong>Bam file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_bam,
                    os.path.basename(path_name_bam),
                    os.path.basename(path_name_bam),
                )
            )
            data["bai_file_id"] = mark_safe(
                '<strong>Bai file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_bai,
                    os.path.basename(path_name_bai),
                    os.path.basename(path_name_bai),
                )
            )
            data["vcf_file_id"] = mark_safe(
                '<strong>Vcf file:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_vcf,
                    os.path.basename(path_name_vcf),
                    os.path.basename(path_name_vcf),
                )
            )
            data["reference_id"] = mark_safe(
                '<strong>Reference:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_reference,
                    os.path.basename(path_name_reference),
                    os.path.basename(path_name_reference),
                )
            )
            data["reference_index_id"] = mark_safe(
                '<strong>Ref. index:</strong> <a href="{}" download="{}"> {}</a>'.format(
                    path_name_reference_index,
                    os.path.basename(path_name_reference_index),
                    os.path.basename(path_name_reference_index),
                )
            )

            data["static_dir"] = run.static_dir
            data["sample_name"] = sample_name

        return JsonResponse(data)
