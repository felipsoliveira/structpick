"""Smoke tests against the real services. Run with: pytest -m network"""

import pytest

from structpick.api import Client
from structpick.select import Filters, select_structures

pytestmark = pytest.mark.network


def test_hemoglobin_end_to_end():
    with Client() as client:
        target = client.target("P69905")
        ranked = select_structures(client.structures(target.accession), Filters(max_resolution=2.0))
        assert target.gene == "HBA1"
        assert ranked
        assert ranked[0].resolution is not None and ranked[0].resolution <= 2.0
        assert client.alphafold(target.accession) is not None
