from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path

import httpx
import pytest

from structpick.api import Client

FIXTURES = Path(__file__).parent / "fixtures"


def load(name: str) -> object:
    return json.loads((FIXTURES / name).read_text())


ROUTES: dict[str, str] = {
    "/uniprotkb/P69905": "uniprot_P69905.json",
    "/pdbe/api/mappings/best_structures/P69905": "best_structures_P69905.json",
    "/api/prediction/P69905": "alphafold_P69905.json",
    "/pdbe/api/pdb/entry/ligand_monomers/1bab": "ligands_1bab.json",
    "/pdbe/api/pdb/entry/ligand_monomers/2dn2": "ligands_2dn2.json",
}


def fixture_handler(request: httpx.Request) -> httpx.Response:
    path = request.url.path
    if path in ROUTES:
        return httpx.Response(200, json=load(ROUTES[path]))
    if path.startswith("/pdbe/api/pdb/entry/ligand_monomers/"):
        entry_ids = path.rsplit("/", 1)[1].lower().split(",")
        if len(entry_ids) > 1:
            data: dict[str, object] = {}
            for entry_id in entry_ids:
                route = f"/pdbe/api/pdb/entry/ligand_monomers/{entry_id}"
                if route in ROUTES:
                    payload = load(ROUTES[route])
                    if isinstance(payload, dict):
                        data.update(payload)
                else:
                    data[entry_id] = []
            return httpx.Response(200, json=data)
        return httpx.Response(200, json={entry_ids[0]: []})
    if path.startswith("/download/") or path.startswith("/files/"):
        return httpx.Response(200, content=b"data_" + path.rsplit("/", 1)[1].encode())
    return httpx.Response(404)


@pytest.fixture
def make_client() -> Callable[..., Client]:
    def factory(handler: Callable[[httpx.Request], httpx.Response] = fixture_handler) -> Client:
        return Client(http=httpx.Client(transport=httpx.MockTransport(handler)), backoff=0)

    return factory


@pytest.fixture
def client(make_client: Callable[..., Client]) -> Client:
    return make_client()
