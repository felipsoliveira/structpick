from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from structpick.models import ADDITIVES, Structure, is_additive

__all__ = ["ADDITIVES", "Filters", "is_additive", "rank_key", "select_structures"]


METHOD_ALIASES = {
    "xray": "x-ray diffraction",
    "em": "electron microscopy",
    "nmr": "solution nmr",
}


@dataclass(frozen=True)
class Filters:
    max_resolution: float | None = None
    min_coverage: float = 0.0
    method: str | None = None
    ligand: str | None = None
    require_ligand: bool = False


def _method_matches(structure: Structure, method: str) -> bool:
    wanted = METHOD_ALIASES.get(method.lower(), method.lower())
    return structure.method.lower() == wanted


def _keep(structure: Structure, filters: Filters) -> bool:
    if filters.method and not _method_matches(structure, filters.method):
        return False
    if structure.coverage < filters.min_coverage:
        return False
    if filters.max_resolution is not None and (
        structure.resolution is None or structure.resolution > filters.max_resolution
    ):
        return False
    bound = {lig.code.upper() for lig in structure.bound_ligands}
    if filters.require_ligand and not bound:
        return False
    return not (filters.ligand and filters.ligand.upper() not in bound)


def rank_key(structure: Structure) -> tuple[bool, float, float, str, str]:
    """Coverage first (in 10% steps, so 99% and 100% tie), then resolution.

    Entries without a resolution (NMR) go last.
    """
    no_resolution = structure.resolution is None
    resolution = structure.resolution if structure.resolution is not None else 0.0
    return (
        no_resolution,
        -round(structure.coverage, 1),
        resolution,
        structure.pdb_id,
        structure.chain,
    )


def select_structures(
    structures: Iterable[Structure], filters: Filters | None = None
) -> list[Structure]:
    filters = filters or Filters()
    return sorted((s for s in structures if _keep(s, filters)), key=rank_key)
