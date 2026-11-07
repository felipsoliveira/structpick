from __future__ import annotations

from dataclasses import dataclass, field

# Things that show up in crystal structures but aren't what you'd dock against:
# solvents, cryoprotectants, buffers, common ions and detergents.
ADDITIVES = frozenset(
    {
        # solvents / cryoprotectants
        "HOH", "DOD", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "P6G", "MPD", "DMS", "IPA",
        "EOH", "MOH", "ACT", "ACE", "FMT", "BME", "DTT", "TRS", "MES", "EPE", "CIT", "FLC",
        "TAR", "MLI", "IMD", "BTB", "NHE", "CAPS", "PO4", "SO4", "NO3", "SCN", "AZI",
        # ions
        "NA", "K", "LI", "CL", "BR", "IOD", "F", "MG", "CA", "MN", "ZN", "FE", "FE2", "CO",
        "NI", "CU", "CU1", "CD", "HG", "SR", "BA", "CS", "RB", "NH4", "UNX",
        # detergents
        "BOG", "LDA", "SDS", "LMT", "DDM", "OLC",
    }
)  # fmt: skip


def is_additive(code: str) -> bool:
    return code.upper() in ADDITIVES


@dataclass(frozen=True)
class Ligand:
    code: str
    name: str
    chain: str


@dataclass(frozen=True)
class Structure:
    pdb_id: str
    chain: str
    method: str
    resolution: float | None
    coverage: float
    unp_start: int
    unp_end: int
    ligands: tuple[Ligand, ...] = field(default=())

    @property
    def bound_ligands(self) -> tuple[Ligand, ...]:
        """Ligands in this chain that aren't crystallisation additives or ions."""
        return tuple(
            lig for lig in self.ligands if lig.chain == self.chain and not is_additive(lig.code)
        )


@dataclass(frozen=True)
class AlphaFoldModel:
    accession: str
    version: int
    mean_plddt: float
    fraction_very_high: float
    start: int
    end: int
    cif_url: str
    pdb_url: str


@dataclass(frozen=True)
class Target:
    accession: str
    name: str
    gene: str | None
    organism: str
    length: int
