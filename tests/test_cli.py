import csv
import io
import json
from importlib.metadata import version

import httpx
import pytest

from structpick.cli import main
from tests.conftest import fixture_handler


def run(client, *argv):
    out = io.StringIO()
    code = main(list(argv), client=client, out=out)
    return code, out.getvalue()


def test_cli_reports_installed_package_version(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])

    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"structpick {version('structpick')}"


def test_list_table(client):
    code, out = run(client, "list", "P69905", "--top", "3")
    assert code == 0
    lines = out.splitlines()
    assert lines[0].startswith("P69905 (HBA1): Hemoglobin subunit alpha")
    assert lines[1].split() == [
        "rank",
        "pdb_id",
        "chain",
        "method",
        "resolution",
        "coverage",
        "residues",
        "ligands",
    ]
    # 2dn2 (1.25 Å, 99%) outranks 1bab (1.5 Å, 100%)
    assert lines[2].split()[1] == "2dn2"
    assert lines[2].split()[-1] == "HEM"
    assert lines[-1].startswith("AlphaFold: AF-P69905-F1 v6")


def test_list_csv_and_filters(client):
    code, out = run(client, "list", "P69905", "--format", "csv", "--max-resolution", "1.6")
    assert code == 0
    rows = list(csv.DictReader(io.StringIO(out)))
    assert {r["pdb_id"] for r in rows} == {"2dn2", "1bab"}


def test_list_json_with_ligand_filter(client):
    code, out = run(client, "list", "P69905", "--format", "json", "--ligand", "HEM")
    assert code == 0
    data = json.loads(out)
    assert data["target"]["accession"] == "P69905"
    assert {s["pdb_id"] for s in data["structures"]} == {
        "2dn2",
        "1bab",
    }  # 2hhb fixture has no ligands
    assert all(s["bound_ligands"] == ["HEM"] for s in data["structures"])


def test_fetch_best_structure(client, tmp_path):
    code, out = run(client, "fetch", "P69905", "-o", str(tmp_path))
    assert code == 0
    assert (tmp_path / "2dn2.cif").read_bytes() == b"data_2DN2.cif"
    assert "chain A" in out


def test_fetch_falls_back_to_alphafold(client, tmp_path):
    code, out = run(client, "fetch", "P69905", "--max-resolution", "0.5", "-o", str(tmp_path))
    assert code == 0
    assert (tmp_path / "AF-P69905-F1-model_v6.cif").exists()
    assert "AlphaFold" in out


def test_unknown_accession_exit_code(client, capsys):
    code, _ = run(client, "list", "NOPE00")
    assert code == 2
    assert "no entry" in capsys.readouterr().err


def test_network_error_exit_code(make_client, capsys):
    def broken(request):
        if request.url.path.startswith("/uniprotkb/"):
            return fixture_handler(request)
        return httpx.Response(500)

    code, _ = run(make_client(broken), "list", "P69905")
    assert code == 3
    assert "request failed" in capsys.readouterr().err
