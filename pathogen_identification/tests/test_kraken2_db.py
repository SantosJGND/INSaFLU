"""
Unit tests for pathogen_identification/utilities/kraken2_db.py.
No Django setup; k2 / dustmasker are faked, registry uses a temp SQLite file.
"""

import os
import types
from datetime import date

import pytest
from sqlalchemy import text

from pathogen_identification.utilities.centrifuge_db import ExportPlan
from pathogen_identification.utilities.kraken2_db import (
    Kraken2DBBuilder, Kraken2DBError, Kraken2DBNaming, Kraken2DBRegistry,
    Kraken2DBStaging, Kraken2DBTaxonomy, validate_kmer_settings)
from pathogen_identification.utilities.utilities_televir_dbs import \
    Utility_Repository

NODES = "1\t|\t1\t|\tno rank\t|\n11320\t|\t1\t|\tspecies\t|\n562\t|\t1\t|\tspecies\t|\n"
NAMES = (
    "1\t|\troot\t|\t\t|\tscientific name\t|\n"
    "11320\t|\tInfluenza A virus\t|\t\t|\tscientific name\t|\n"
    "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n"
)


def write_text(path, text_):
    with open(path, "w") as handle:
        handle.write(text_)
    return str(path)


def fake_result(returncode=0, stdout=""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr="")


@pytest.fixture
def fake_k2(monkeypatch):
    monkeypatch.setattr(Kraken2DBBuilder, "resolve_binary", lambda self: "/fake/k2")


def make_taxonomy(tmp_path, nodes=NODES, names=NAMES, merged=None):
    write_text(tmp_path / "nodes.dmp", nodes)
    write_text(tmp_path / "names.dmp", names)
    if merged is not None:
        write_text(tmp_path / "merged.dmp", merged)
    return Kraken2DBTaxonomy(str(tmp_path))


# ------------------------------------------------------------ validation


@pytest.mark.parametrize("k,l", [(35, 31), (31, 31), (21, 11)])
def test_valid_kmer_settings(k, l):
    validate_kmer_settings(k, l)


@pytest.mark.parametrize(
    "k,l,msg",
    [
        (35, 32, "<= 31"),
        (20, 25, "<= --kmer-len"),
        (0, 1, "kmer-len"),
        (35, 0, "minimizer-len"),
    ],
)
def test_invalid_kmer_settings(k, l, msg):
    with pytest.raises(Kraken2DBError, match=msg):
        validate_kmer_settings(k, l)


# ---------------------------------------------------------------- naming


def test_naming_directories(tmp_path):
    naming = Kraken2DBNaming(base_dir=str(tmp_path), build_date=date(2026, 10, 14))
    assert naming.name == "kraken2_20261014"
    assert naming.final_base == os.path.join(str(tmp_path), "kraken2_20261014")
    assert naming.tmp_base == os.path.join(str(tmp_path), "kraken2_20261014.tmp")
    assert naming.will_replace() is False
    os.makedirs(naming.final_base)
    assert naming.will_replace() is True


# -------------------------------------------------------------- taxonomy


def test_merged_classification(tmp_path):
    tax = make_taxonomy(tmp_path, merged="9999\t|\t562\t|\n")
    unresolvable, remappable = tax.classify_missing({"562", "9999", "123456"})
    assert "9999" in remappable
    assert unresolvable == ["123456"]


def test_merged_absent_is_fine(tmp_path):
    tax = make_taxonomy(tmp_path)
    assert tax.merged_map() == {}
    assert tax.classify_missing({"562"}) == ([], [])


# --------------------------------------------------------------- staging


def test_staging_layout(tmp_path):
    src = tmp_path / "tax"
    src.mkdir()
    tax = make_taxonomy(src, merged="7\t|\t562\t|\n")
    plan = ExportPlan(
        accid_to_taxid={"NC_1.1": "11320", "NC_2.1": "562"},
        accid_to_file={"NC_1.1": "x", "NC_2.1": "x"},
    )
    tmp_dir = str(tmp_path / "kraken2_20261014.tmp")
    written = Kraken2DBStaging.stage(
        plan, {"NC_1.1": "ACGT", "NC_2.1": "GGCC"}, tax, tmp_dir
    )

    assert written == 2
    assert open(os.path.join(tmp_dir, "seqid2taxid.map")).read() == "NC_1.1\t11320\nNC_2.1\t562\n"
    fna = open(os.path.join(tmp_dir, "library", "input.fna")).read()
    assert fna == ">NC_1.1\nACGT\n>NC_2.1\nGGCC\n"
    assert os.path.isfile(os.path.join(tmp_dir, "taxonomy", "nodes.dmp"))
    assert os.path.isfile(os.path.join(tmp_dir, "taxonomy", "names.dmp"))
    assert os.path.isfile(os.path.join(tmp_dir, "taxonomy", "merged.dmp"))


# --------------------------------------------------------------- builder


def test_build_command_arguments(fake_k2):
    cmd = Kraken2DBBuilder(threads=4).build_command("/db", 35, 31)
    assert cmd == [
        "/fake/k2", "build", "--db", "/db", "--threads", "4",
        "--kmer-len", "35", "--minimizer-len", "31",
    ]


def test_build_command_max_db_size(fake_k2):
    cmd = Kraken2DBBuilder(threads=2).build_command("/db", 31, 21, max_db_size="10GiB")
    assert cmd[-2:] == ["--max-db-size", "10GiB"]


def test_build_raises_on_nonzero_exit(tmp_path, fake_k2):
    calls = []

    def runner(cmd, **kwargs):
        calls.append(cmd)
        return fake_result(returncode=2)

    with pytest.raises(Kraken2DBError, match="exit code 2"):
        Kraken2DBBuilder(runner=runner).build("/db", 35, 31, str(tmp_path / "log"))
    assert calls[0][1] == "build"


def test_verify_missing_files(tmp_path, fake_k2):
    builder = Kraken2DBBuilder(runner=lambda cmd, **kw: fake_result())
    with pytest.raises(Kraken2DBError, match="hash.k2d"):
        builder.verify(str(tmp_path))


def test_verify_inspect_failure(tmp_path, fake_k2):
    for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        write_text(tmp_path / f, "x")
    builder = Kraken2DBBuilder(runner=lambda cmd, **kw: fake_result(returncode=1))
    with pytest.raises(Kraken2DBError, match="inspect failed"):
        builder.verify(str(tmp_path))


def test_verify_passes(tmp_path, fake_k2):
    for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        write_text(tmp_path / f, "x")
    calls = []

    def runner(cmd, **kw):
        calls.append(cmd)
        return fake_result()

    Kraken2DBBuilder(runner=runner).verify(str(tmp_path))
    assert calls[0] == ["/fake/k2", "inspect", "--db", str(tmp_path), "--skip-counts"]


def test_require_masker(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(Kraken2DBError, match="dustmasker"):
        Kraken2DBBuilder.require_masker()


# ------------------------------------------------------ promote / cleanup


def test_promote_replaces_same_day_directory(tmp_path):
    tmp_dir = str(tmp_path / "kraken2_1.tmp")
    final_dir = str(tmp_path / "kraken2_1")
    os.makedirs(tmp_dir)
    write_text(os.path.join(tmp_dir, "hash.k2d"), "new")
    os.makedirs(final_dir)
    write_text(os.path.join(final_dir, "hash.k2d"), "old")

    Kraken2DBBuilder.promote(tmp_dir, final_dir)

    assert open(os.path.join(final_dir, "hash.k2d")).read() == "new"
    assert not os.path.exists(tmp_dir)


def test_promote_leaves_other_dates(tmp_path):
    other = str(tmp_path / "kraken2_20261013")
    os.makedirs(other)
    tmp_dir = str(tmp_path / "kraken2_20261014.tmp")
    os.makedirs(tmp_dir)
    Kraken2DBBuilder.promote(tmp_dir, str(tmp_path / "kraken2_20261014"))
    assert os.path.isdir(other)


def test_remove_partial(tmp_path):
    tmp_dir = tmp_path / "kraken2_1.tmp"
    tmp_dir.mkdir()
    write_text(tmp_dir / "x", "x")
    Kraken2DBBuilder.remove_partial(str(tmp_dir))
    assert not tmp_dir.exists()


# --------------------------------------------------------------- registry


def _rows(repo, name):
    with repo.engine.connect() as conn:
        return conn.execute(
            text("SELECT name, db_category, db_type, software, installed, path "
                 "FROM database WHERE name = :n"),
            {"n": name},
        ).fetchall()


def test_kraken2_record_fields(tmp_path):
    naming = Kraken2DBNaming(base_dir=str(tmp_path), build_date=date(2026, 10, 14))
    rec = Kraken2DBRegistry.build_record(naming, 7, 4)
    assert rec["name"] == "kraken2_20261014"
    assert rec["db_category"] == "kraken2"
    assert rec["software"] == "kraken2"
    assert rec["db_type"] == "classification"
    assert rec["installed"] == "True"
    assert rec["path"] == naming.final_base


def test_kraken2_registration_roundtrip(tmp_path):
    repo = Utility_Repository(db_path=str(tmp_path / "reg"), install_type="local")
    naming = Kraken2DBNaming(base_dir=str(tmp_path / "db"), build_date=date(2026, 10, 14))
    repo.upsert_database(Kraken2DBRegistry.build_record(naming, 2, 2))
    rows = _rows(repo, "kraken2_20261014")
    assert rows == [("kraken2_20261014", "kraken2", "classification", "kraken2", "True", naming.final_base)]


def test_upsert_scoped_by_category(tmp_path):
    """A row with the same name in another category must not be deleted."""
    repo = Utility_Repository(db_path=str(tmp_path / "reg"), install_type="local")
    repo.upsert_database({"name": "viral", "db_category": "centrifuge", "path": "/a", "installed": "True"})
    repo.upsert_database({"name": "viral", "db_category": "kraken2", "path": "/b", "installed": "True"})
    repo.upsert_database({"name": "viral", "db_category": "kraken2", "path": "/c", "installed": "True"})

    with repo.engine.connect() as conn:
        rows = conn.execute(
            text("SELECT db_category, path FROM database WHERE name = 'viral' ORDER BY db_category")
        ).fetchall()
    assert rows == [("centrifuge", "/a"), ("kraken2", "/c")]


def test_upsert_null_category_scoped(tmp_path):
    repo = Utility_Repository(db_path=str(tmp_path / "reg"), install_type="local")
    repo.upsert_database({"name": "n1", "path": "/a"})
    repo.upsert_database({"name": "n1", "path": "/b"})
    with repo.engine.connect() as conn:
        count = conn.execute(text("SELECT COUNT(*) FROM database WHERE name = 'n1'")).scalar()
    assert count == 1
