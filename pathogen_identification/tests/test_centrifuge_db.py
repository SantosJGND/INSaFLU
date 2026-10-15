"""
Unit tests for pathogen_identification/utilities/centrifuge_db.py.

No Django setup needed: DB querying is not covered here; binaries are faked.
"""

import gzip
import os
import types
from datetime import date

import pytest

from pathogen_identification.utilities.centrifuge_db import (
    CentrifugeDBBuilder, CentrifugeDBError, CentrifugeDBExporter,
    CentrifugeDBNaming, CentrifugeDBTaxonomy, ExportPlan,
    iter_fasta_records)


# ----------------------------------------------------------------- helpers


def write_text(path, text):
    with open(path, "w") as handle:
        handle.write(text)
    return str(path)


def fake_result(returncode=0, stdout=""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr="")


# ------------------------------------------------------------------- FASTA


def test_iter_fasta_records_selects_and_uppercases(tmp_path):
    path = write_text(
        tmp_path / "refs.fa",
        ">NC_1.1 some description\nacgt\nAC\n>NC_2.1\nGGGG\n>NC_3.1\nTTTT\n",
    )
    records = dict(iter_fasta_records(path, {"NC_1.1", "NC_3.1"}))
    assert records == {"NC_1.1": "ACGTAC", "NC_3.1": "TTTT"}


def test_iter_fasta_records_gzip(tmp_path):
    path = str(tmp_path / "refs.fa.gz")
    with gzip.open(path, "wt") as handle:
        handle.write(">A1\nacgt\n>B1\nccc\n")
    assert dict(iter_fasta_records(path, {"B1"})) == {"B1": "CCC"}


# ------------------------------------------------------------- exporter


def test_normalise_taxid():
    assert CentrifugeDBExporter.normalise_taxid("00562") == "562"
    assert CentrifugeDBExporter.normalise_taxid(" 11320 ") == "11320"
    assert CentrifugeDBExporter.normalise_taxid("abc") is None
    assert CentrifugeDBExporter.normalise_taxid("") is None
    assert CentrifugeDBExporter.normalise_taxid(None) is None


def test_exporter_write_outputs_headers_and_conversion(tmp_path):
    plan = ExportPlan(
        accid_to_taxid={"A1.1": "11320", "B1.1": "562", "C1.1": "9606"},
        accid_to_file={"A1.1": "x", "B1.1": "x", "C1.1": "x"},
    )
    sequences = {"A1.1": "A" * 100, "B1.1": "CCGT"}  # C1.1 has no sequence

    fna = str(tmp_path / "input.fna")
    tsv = str(tmp_path / "conversion.tsv")
    written = CentrifugeDBExporter.write(plan, sequences, fna, tsv)

    assert written == 2
    assert open(tsv).read() == "A1.1\t11320\nB1.1\t562\n"

    lines = open(fna).read().splitlines()
    assert lines[0] == ">A1.1"
    assert lines[1] == "A" * 80  # wrapped at 80
    assert lines[2] == "A" * 20
    assert ">B1.1" in lines and "CCGT" in lines
    assert not any("C1.1" in line for line in lines)


def test_read_sequences_reports_not_found(tmp_path):
    ref = write_text(tmp_path / "ref.fa", ">A1.1\nACGT\n")
    plan = ExportPlan(
        accid_to_taxid={"A1.1": "1", "Z9.9": "2"},
        accid_to_file={"A1.1": ref, "Z9.9": ref},
    )
    seqs, not_found = CentrifugeDBExporter(user=None).read_sequences(plan)
    assert seqs == {"A1.1": "ACGT"}
    assert not_found == ["Z9.9"]


# ------------------------------------------------------------ taxonomy


NODES = "1\t|\t1\t|\tno rank\t|\n11320\t|\t1\t|\tspecies\t|\n562\t|\t1\t|\tspecies\t|\n"
NAMES = (
    "1\t|\troot\t|\t\t|\tscientific name\t|\n"
    "11320\t|\tInfluenza A virus\t|\t\t|\tscientific name\t|\n"
    "11320\t|\tflu\t|\t\t|\tcommon name\t|\n"
    "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n"
)


def make_taxonomy(tmp_path, nodes=NODES, names=NAMES):
    write_text(tmp_path / "nodes.dmp", nodes)
    write_text(tmp_path / "names.dmp", names)
    return CentrifugeDBTaxonomy(str(tmp_path))


def test_taxonomy_missing_files_raise(tmp_path):
    tax = CentrifugeDBTaxonomy(str(tmp_path))
    with pytest.raises(CentrifugeDBError, match="nodes.dmp"):
        tax.check_files()


def test_taxonomy_missing_taxids(tmp_path):
    tax = make_taxonomy(tmp_path)
    # 9606 not in nodes; 1 has a scientific name; 11320/562 are fine
    assert tax.missing_taxids({"11320", "562", "9606"}) == ["9606"]


def test_taxonomy_only_scientific_names_count(tmp_path):
    names = "77\t|\tsome common\t|\t\t|\tcommon name\t|\n"
    nodes = "77\t|\t1\t|\tspecies\t|\n"
    tax = make_taxonomy(tmp_path, nodes=nodes, names=names)
    assert tax.missing_taxids({"77"}) == ["77"]


# ------------------------------------------------------------- naming


def test_naming_uses_date_and_tmp_suffix(tmp_path):
    naming = CentrifugeDBNaming(base_dir=str(tmp_path), build_date=date(2026, 10, 14))
    assert naming.name == "centrifuge_20261014"
    assert naming.final_base == os.path.join(str(tmp_path), "centrifuge_20261014")
    assert naming.tmp_base == os.path.join(str(tmp_path), "centrifuge_20261014.tmp")
    assert naming.will_replace() is False


def test_naming_will_replace_when_same_day_index_exists(tmp_path):
    naming = CentrifugeDBNaming(base_dir=str(tmp_path), build_date=date(2026, 10, 14))
    write_text(tmp_path / "centrifuge_20261014.1.cf", "x")
    assert naming.will_replace() is True


def test_naming_different_date_does_not_replace(tmp_path):
    write_text(tmp_path / "centrifuge_20261013.1.cf", "x")
    naming = CentrifugeDBNaming(base_dir=str(tmp_path), build_date=date(2026, 10, 14))
    assert naming.will_replace() is False


# ------------------------------------------------------------- builder


@pytest.fixture
def fake_binaries(monkeypatch):
    monkeypatch.setattr(
        CentrifugeDBBuilder,
        "resolve_binary",
        staticmethod(lambda name: f"/fake/bin/{name}"),
    )


def test_build_command_arguments(tmp_path, fake_binaries):
    tax = make_taxonomy(tmp_path)
    builder = CentrifugeDBBuilder(threads=4)
    cmd = builder.build_command("in.fna", "conv.tsv", tax, "/x/centrifuge_1.tmp")
    assert cmd == [
        "/fake/bin/centrifuge-build",
        "-p", "4",
        "--conversion-table", "conv.tsv",
        "--taxonomy-tree", tax.nodes_path,
        "--name-table", tax.names_path,
        "in.fna",
        "/x/centrifuge_1.tmp",
    ]


def test_build_raises_on_nonzero_exit(tmp_path, fake_binaries):
    tax = make_taxonomy(tmp_path)
    calls = []

    def runner(cmd, **kwargs):
        calls.append(cmd)
        return fake_result(returncode=3)

    builder = CentrifugeDBBuilder(runner=runner)
    with pytest.raises(CentrifugeDBError, match="exit code 3"):
        builder.build("in.fna", "conv.tsv", tax, "/x/t", str(tmp_path / "log"))
    assert calls[0][0] == "/fake/bin/centrifuge-build"


def test_verify_count_mismatch_raises(tmp_path, fake_binaries):
    base = str(tmp_path / "centrifuge_1.tmp")
    write_text(tmp_path / "centrifuge_1.tmp.1.cf", "x")
    builder = CentrifugeDBBuilder(runner=lambda cmd, **kw: fake_result(stdout="A\nB\n"))
    with pytest.raises(CentrifugeDBError, match="differs"):
        builder.verify(base, expected=3)


def test_verify_passes_on_matching_count(tmp_path, fake_binaries):
    base = str(tmp_path / "centrifuge_1.tmp")
    write_text(tmp_path / "centrifuge_1.tmp.1.cf", "x")
    builder = CentrifugeDBBuilder(runner=lambda cmd, **kw: fake_result(stdout="A\nB\n"))
    builder.verify(base, expected=2)


def test_verify_without_index_files_raises(tmp_path, fake_binaries):
    builder = CentrifugeDBBuilder(runner=lambda cmd, **kw: fake_result(stdout=""))
    with pytest.raises(CentrifugeDBError, match="No index files"):
        builder.verify(str(tmp_path / "missing"), expected=0)


# ------------------------------------------------------- promote / cleanup


def test_promote_renames_and_replaces_existing(tmp_path):
    tmp_base = str(tmp_path / "centrifuge_1.tmp")
    final_base = str(tmp_path / "centrifuge_1")
    for ext in ("1.cf", "2.cf", "3.cf"):
        write_text(tmp_path / f"centrifuge_1.tmp.{ext}", "new")
    write_text(tmp_path / "centrifuge_1.1.cf", "old")

    new_paths = CentrifugeDBBuilder.promote(tmp_base, final_base)

    assert sorted(os.path.basename(p) for p in new_paths) == [
        "centrifuge_1.1.cf", "centrifuge_1.2.cf", "centrifuge_1.3.cf",
    ]
    assert open(tmp_path / "centrifuge_1.1.cf").read() == "new"
    assert not any("tmp" in f for f in os.listdir(tmp_path))


def test_promote_without_files_raises(tmp_path):
    with pytest.raises(CentrifugeDBError, match="Nothing to promote"):
        CentrifugeDBBuilder.promote(str(tmp_path / "a.tmp"), str(tmp_path / "a"))


def test_remove_partial_only_removes_tmp(tmp_path):
    write_text(tmp_path / "centrifuge_1.tmp.1.cf", "x")
    write_text(tmp_path / "centrifuge_1.1.cf", "keep")
    CentrifugeDBBuilder.remove_partial(str(tmp_path / "centrifuge_1.tmp"))
    assert not (tmp_path / "centrifuge_1.tmp.1.cf").exists()
    assert (tmp_path / "centrifuge_1.1.cf").exists()


# ------------------------------------------------------------- registry


from pathogen_identification.utilities.centrifuge_db import CentrifugeDBRegistry  # noqa: E402
from pathogen_identification.utilities.utilities_televir_dbs import \
    Utility_Repository  # noqa: E402
from sqlalchemy import text  # noqa: E402


def _rows(repo, name):
    with repo.engine.connect() as conn:
        return conn.execute(
            text("SELECT name, db_category, db_type, software, installed, path FROM database WHERE name = :n"),
            {"n": name},
        ).fetchall()


def test_upsert_inserts_and_is_queryable(tmp_path):
    repo = Utility_Repository(db_path=str(tmp_path / "reg"), install_type="local")
    naming = CentrifugeDBNaming(base_dir=str(tmp_path / "idx"), build_date=date(2026, 10, 14))
    record = CentrifugeDBRegistry.build_record(naming, 2, 2)
    repo.upsert_database(record)

    rows = _rows(repo, "centrifuge_20261014")
    assert len(rows) == 1
    name, category, db_type, software, installed, path = rows[0]
    assert (category, db_type, software) == ("centrifuge", "classification", "centrifuge")
    assert path == naming.final_base
    assert installed == "True"


def test_upsert_same_name_replaces_row(tmp_path):
    repo = Utility_Repository(db_path=str(tmp_path / "reg"), install_type="local")
    naming = CentrifugeDBNaming(base_dir=str(tmp_path / "idx"), build_date=date(2026, 10, 14))
    repo.upsert_database(CentrifugeDBRegistry.build_record(naming, 2, 2))
    repo.upsert_database(CentrifugeDBRegistry.build_record(naming, 5, 3))
    assert len(_rows(repo, "centrifuge_20261014")) == 1


def test_different_date_adds_row(tmp_path):
    repo = Utility_Repository(db_path=str(tmp_path / "reg"), install_type="local")
    for d in (date(2026, 10, 13), date(2026, 10, 14)):
        naming = CentrifugeDBNaming(base_dir=str(tmp_path / "idx"), build_date=d)
        repo.upsert_database(CentrifugeDBRegistry.build_record(naming, 1, 1))
    assert len(_rows(repo, "centrifuge_20261013")) == 1
    assert len(_rows(repo, "centrifuge_20261014")) == 1


def test_upsert_handles_quotes_safely(tmp_path):
    repo = Utility_Repository(db_path=str(tmp_path / "reg"), install_type="local")
    repo.upsert_database({"name": "x'y", "path": "/a'b", "installed": True})
    assert len(_rows(repo, "x'y")) == 1


def test_upsert_requires_name(tmp_path):
    repo = Utility_Repository(db_path=str(tmp_path / "reg"), install_type="local")
    with pytest.raises(ValueError):
        repo.upsert_database({"path": "/a"})


def test_build_record_fields(tmp_path):
    naming = CentrifugeDBNaming(base_dir=str(tmp_path), build_date=date(2026, 10, 14))
    rec = CentrifugeDBRegistry.build_record(naming, 7, 4)
    assert rec["db_category"] == "centrifuge"
    assert rec["db_type"] == "classification"
    assert rec["software"] == "centrifuge"
    assert rec["date"] == "2026-10-14"
    assert "7 accessions" in rec["description"]
