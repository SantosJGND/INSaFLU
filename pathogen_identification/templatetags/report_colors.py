import os

from django import template
from django.utils.safestring import mark_safe

from pathogen_identification.models import FinalReport
from pathogen_identification.utilities.televir_parameters import TelevirParameters

register = template.Library()
cell_color = "2, 155, 194"
cell_color = "255, 210, 112"


@register.filter(name="color_code")
def coverage_col(coverage_value):
    if coverage_value is None or coverage_value == "":
        coverage_value = 0
    ncol = f"background-color: rgba({cell_color}, {int(coverage_value)}%);"
    return ncol


@register.filter(name="round_str")
def round_str(value):
    if value is None or value == "":
        value = 0
    value = float(value)
    ncol = round(value, 2)
    return ncol


@register.filter(name="round_str_perc_invert")
def round_str(value):
    if value is None or value == "":
        value = 0
    value = float(value)
    ncol = round(100 - value, 2)
    return ncol


@register.filter(name="round")
def round_num(value):
    if value is None or value == "":
        value = 0
    ncol = round(value, 2)
    return ncol


@register.filter(name="round_to_int")
def round_to_int(value):
    if value is None or value == "":
        value = 0

    value = round(value, 5)
    return value


@register.filter(name="round_to_percent")
def round_to_percent(value):
    if value is None or value == "":
        value = 0

    value = round(value * 100, 2)
    return value


@register.filter("reconvert_string_to_int")
def reconvert_string_to_int(value):
    if value is None or value == "":
        value = 0
    elif "," in value:
        value = value.replace(",", "")

    return int(value)


@register.filter("convert_int_to_str_format")
def convert_int_to_str_format(value):
    if value is None or value == "":
        value = 0

    return f"{int(value):,}"


@register.filter(name="success_code")
def map_success_col(success_count):
    ncol = f"background-color: rgba({cell_color}, {50 * success_count}%);"
    return ncol


@register.simple_tag
def depth_color_error(depth_value, max_value):
    if depth_value and max_value:
        ncol = float(depth_value) * 100 / float(max_value)
    else:
        ncol = 0

    ncol = 100 - ncol

    ncol = f"background-color: rgba({cell_color}, {int(ncol)}%);"
    return ncol


@register.simple_tag
def depth_color(depth_value, max_value):
    if depth_value and max_value:
        ncol = float(depth_value) * 100 / float(max_value)
    else:
        ncol = 0

    ncol = f"background-color: rgba({cell_color}, {int(ncol)}%);"
    return ncol


@register.simple_tag
def depth_color_gaps(depth_value, max_value):
    if depth_value and max_value:
        ncol = float(depth_value) * 100 / float(max_value)
    else:
        ncol = 0

    ncol = 100 - ncol

    ncol = f"background-color: rgba({cell_color}, {int(ncol)}%);"
    return ncol


@register.simple_tag
def depth_color_windows(window_value: str, max_prop):
    if window_value and max_prop:
        if "/" in window_value:
            window_value = window_value.split("/")

            window_value = int(window_value[0]) / int(window_value[1])
        elif window_value == "NA":
            window_value = 0
        ncol = float(window_value) * 100 / float(max_prop)
    else:
        ncol = 0

    if str(ncol) == "nan":
        ncol = 0

    ncol = f"background-color: rgba({cell_color}, {int(ncol)}%);"
    return ncol


@register.simple_tag
def flag_false_positive(depth, depthc, coverage, mapped, windows_covered, project_pk):
    flag_build = TelevirParameters.get_flag_build(project_pk)
    flag_build = flag_build(depth, depthc, coverage, mapped, windows_covered)

    if flag_build.assert_false_positive():
        return "Likely False Positive"

    elif flag_build.assert_vestigial():
        return "Vestigial Mapping"

    return ""


@register.simple_tag
def success_count_color(success):
    counts_dict = {
        "": 0,
        "none": 0,
        "reads": 1,
        "contigs": 1,
        "reads and contigs": 2,
    }

    ncol = f"background-color: rgba({cell_color}, {counts_dict[success] * 50}%);"
    return ncol


@register.simple_tag
def flag_false_positive_color(
    depth, depthc, coverage, mapped, windows_covered, project_pk
):
    flag_build = TelevirParameters.get_flag_build(project_pk=project_pk)

    flag_build = flag_build(depth, depthc, coverage, mapped, windows_covered)

    if flag_build.assert_false_positive():
        return "background-color: rgba(255, 0, 0, 0.5);"

    elif flag_build.assert_vestigial():
        return "background-color: rgba(255, 0, 0, 0.5);"

    return "background-color: rgba(169, 169, 169, 1);"  # Dark gray


@register.simple_tag
def flag_control_color(flag):
    """Set background color to red if flag is True, otherwise set it to dark gray."""
    if flag in [FinalReport.CONTROL_FLAG_PRESENT]:
        return "background-color: rgba(255, 0, 0, 0.5);"
    return "background-color: rgba(169, 169, 169, 1);"  # Dark gray
