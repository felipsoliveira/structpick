from __future__ import annotations

import time
from collections.abc import Iterable
from dataclasses import replace
from typing import Any

import httpx

from structpick import __version__
from structpick.models import AlphaFoldModel, Ligand, Structure, Target

PDBE = "https://www.ebi.ac.uk/pdbe/api"
ALPHAFOLD = "https://alphafold.ebi.ac.uk/api"
UNIPROT = "https://rest.uniprot.org"

RETRY_STATUS = {429, 500, 502, 503, 504}
MAX_LIGAND_BATCH_SIZE = 500


class TargetNotFoundError(LookupError):
    pass


class Client:
    """Small wrapper around the UniProt, PDBe and AlphaFold DB REST APIs."""

    def __init__(
        self,
        http: httpx.Client | None = None,
        retries: int = 3,
        backoff: float = 1.0,
    ) -> None:
        self._http = http or httpx.Client(
            timeout=httpx.Timeout(30.0),
            headers={"User-Agent": f"structpick/{__version__}"},
            follow_redirects=True,
        )
        self._retries = retries
        self._backoff = backoff

    def __enter__(self) -> Client:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def close(self) -> None:
        self._http.close()

    def _get(
        self,
        url: str,
        params: dict[str, str] | None = None,
        missing: frozenset[int] = frozenset({404}),
    ) -> Any | None:
        """GET JSON. Returns None for `missing` statuses; retries rate limits and server errors."""
        for attempt in range(self._retries + 1):
            try:
                response = self._http.get(url, params=params)
            except httpx.TransportError:
                if attempt == self._retries:
                    raise
            else:
                if response.status_code in missing:
                    return None
                if response.status_code not in RETRY_STATUS or attempt == self._retries:
                    response.raise_for_status()
                    return response.json()
            time.sleep(self._backoff * 2**attempt)
        raise AssertionError("unreachable")

    def target(self, accession: str) -> Target:
        data = self._get(
            f"{UNIPROT}/uniprotkb/{accession}",
            params={
                "fields": "accession,protein_name,gene_names,organism_name,length",
                "format": "json",
            },
            # UniProt answers 400 for strings that aren't valid accessions.
            missing=frozenset({400, 404}),
        )
        if data is None:
            raise TargetNotFoundError(f"UniProt has no entry {accession!r}")
        description: dict[str, Any] = data.get("proteinDescription", {})
        submissions: list[dict[str, Any]] = description.get("submissionNames", [])
        name = (
            description.get("recommendedName", {}).get("fullName", {}).get("value")
            or (submissions[0].get("fullName", {}).get("value") if submissions else None)
            or accession
        )
        genes = data.get("genes") or []
        gene = genes[0].get("geneName", {}).get("value") if genes else None
        return Target(
            accession=data["primaryAccession"],
            name=name,
            gene=gene,
            organism=data.get("organism", {}).get("scientificName", ""),
            length=data.get("sequence", {}).get("length", 0),
        )

    def structures(self, accession: str) -> list[Structure]:
        data = self._get(f"{PDBE}/mappings/best_structures/{accession}")
        if not data:
            return []
        rows: list[dict[str, Any]] = data.get(accession) or next(iter(data.values()), [])
        return [
            Structure(
                pdb_id=row["pdb_id"],
                chain=row["chain_id"],
                method=row.get("experimental_method", ""),
                resolution=row.get("resolution"),
                coverage=float(row.get("coverage", 0.0)),
                unp_start=row.get("unp_start", 0),
                unp_end=row.get("unp_end", 0),
            )
            for row in rows
        ]

    def ligands_many(self, pdb_ids: Iterable[str]) -> dict[str, tuple[Ligand, ...]]:
        """Fetch ligand metadata for multiple PDB entries in bounded batches."""
        entry_ids = list(dict.fromkeys(pdb_id.lower() for pdb_id in pdb_ids))
        if not entry_ids:
            return {}
        result: dict[str, tuple[Ligand, ...]] = {}
        for offset in range(0, len(entry_ids), MAX_LIGAND_BATCH_SIZE):
            batch = entry_ids[offset : offset + MAX_LIGAND_BATCH_SIZE]
            data = self._get(f"{PDBE}/pdb/entry/ligand_monomers/{','.join(batch)}") or {}
            result.update(
                {
                    entry_id: tuple(
                        Ligand(
                            code=row["chem_comp_id"],
                            name=row.get("chem_comp_name", ""),
                            chain=row["chain_id"],
                        )
                        for row in data.get(entry_id, data.get(entry_id.upper(), []))
                    )
                    for entry_id in batch
                }
            )
        return result

    def ligands(self, pdb_id: str) -> tuple[Ligand, ...]:
        """Fetch ligand metadata for a single PDB entry."""
        return self.ligands_many((pdb_id,)).get(pdb_id.lower(), ())

    def with_ligands(self, structures: list[Structure]) -> list[Structure]:
        """Attach ligand metadata, batching unique PDB entries into bounded requests."""
        ligands_by_entry = self.ligands_many(s.pdb_id for s in structures)
        return [replace(s, ligands=ligands_by_entry.get(s.pdb_id.lower(), ())) for s in structures]

    def alphafold(self, accession: str) -> AlphaFoldModel | None:
        data = self._get(f"{ALPHAFOLD}/prediction/{accession}")
        if not data:
            return None
        row = next((r for r in data if r.get("uniprotAccession") == accession), data[0])
        return AlphaFoldModel(
            accession=row["uniprotAccession"],
            version=row.get("latestVersion", 0),
            mean_plddt=row.get("globalMetricValue", 0.0),
            fraction_very_high=row.get("fractionPlddtVeryHigh", 0.0),
            start=row.get("sequenceStart", 0),
            end=row.get("sequenceEnd", 0),
            cif_url=row["cifUrl"],
            pdb_url=row["pdbUrl"],
        )

    def download(self, url: str) -> bytes:
        response = self._http.get(url)
        response.raise_for_status()
        return response.content


def entry_url(pdb_id: str, fmt: str = "cif") -> str:
    return f"https://files.rcsb.org/download/{pdb_id.upper()}.{fmt}"
