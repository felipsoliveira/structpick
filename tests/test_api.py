import httpx
import pytest

from structpick.api import TargetNotFoundError


def test_target(client):
    target = client.target("P69905")
    assert target.name == "Hemoglobin subunit alpha"
    assert target.gene == "HBA1"
    assert target.organism == "Homo sapiens"
    assert target.length == 142


def test_unknown_target_raises(client):
    with pytest.raises(TargetNotFoundError):
        client.target("NOPE00")


def test_malformed_accession_raises_not_found(make_client):
    client = make_client(lambda request: httpx.Response(400, json={"messages": ["invalid"]}))
    with pytest.raises(TargetNotFoundError):
        client.target("XXXX99")


def test_structures(client):
    structures = client.structures("P69905")
    assert len(structures) == 6
    first = structures[0]
    assert (first.pdb_id, first.chain, first.resolution, first.coverage) == ("1bab", "A", 1.5, 1.0)


def test_no_structures_is_empty(client):
    assert client.structures("Q00000") == []


def test_ligands_are_batched_once_per_unique_entry(make_client):
    from tests.conftest import fixture_handler

    calls = []

    def counting(request):
        calls.append(request.url.path)
        return fixture_handler(request)

    client = make_client(counting)
    structures = client.with_ligands(client.structures("P69905"))
    ligand_calls = [c for c in calls if "ligand_monomers" in c]
    assert len(ligand_calls) == 1
    assert ligand_calls[0].endswith("/1bab,2dn2,2hhb")
    assert {lig.code for lig in structures[0].bound_ligands} == {"HEM"}


def test_ligand_batches_limit_url_length(make_client):
    calls = []

    def handler(request):
        entry_ids = request.url.path.rsplit("/", 1)[1].split(",")
        calls.append(entry_ids)
        return httpx.Response(200, json={entry_id: [] for entry_id in entry_ids})

    client = make_client(handler)
    entry_ids = [f"{i:04x}" for i in range(501)]

    assert len(client.ligands_many(entry_ids)) == 501
    assert [len(batch) for batch in calls] == [500, 1]


def test_alphafold(client):
    model = client.alphafold("P69905")
    assert model is not None
    assert model.mean_plddt > 90
    assert model.cif_url.endswith(".cif")
    assert client.alphafold("Q00000") is None


def test_retries_on_rate_limit(make_client):
    responses = iter(
        [httpx.Response(429), httpx.Response(503), httpx.Response(200, json={"Q1": []})]
    )
    client = make_client(lambda request: next(responses))
    assert client.structures("Q1") == []


def test_gives_up_after_retries(make_client):
    client = make_client(lambda request: httpx.Response(503))
    with pytest.raises(httpx.HTTPStatusError):
        client.structures("Q1")
