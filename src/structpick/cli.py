from __future__ import annotations

import argparse
import csv
import json
import sys
from collections.abc import Sequence
from dataclasses import asdict, replace
from pathlib import Path
from typing import TextIO

import httpx

from structpick import __version__
from structpick.api import Client, TargetNotFoundError, entry_url
from structpick.models import AlphaFoldModel, Structure, Target
from structpick.select import Filters, select_structures

COLUMNS = ["rank", "pdb_id", "chain", "method", "resolution", "coverage", "residues", "ligands"]


def _row(rank: int, s: Structure) -> dict[str, str]:
    return {
        "rank": str(rank),
        "pdb_id": s.pdb_id,
        "chain": s.chain,
        "method": s.method,
        "resolution": f"{s.resolution:.2f}" if s.resolution is not None else "",
        "coverage": f"{s.coverage:.0%}",
        "residues": f"{s.unp_start}-{s.unp_end}",
        "ligands": ",".join(sorted({lig.code for lig in s.bound_ligands})),
    }


def _print_table(
    target: Target, rows: list[dict[str, str]], af: AlphaFoldModel | None, out: TextIO
) -> None:
    gene = f" ({target.gene})" if target.gene else ""
    print(
        f"{target.accession}{gene}: {target.name} | {target.organism} | {target.length} aa",
        file=out,
    )
    if rows:
        widths = {c: max(len(c), *(len(r[c]) for r in rows)) for c in COLUMNS}
        print("  ".join(c.ljust(widths[c]) for c in COLUMNS).rstrip(), file=out)
        for r in rows:
            print("  ".join(r[c].ljust(widths[c]) for c in COLUMNS).rstrip(), file=out)
    else:
        print("No experimental structure matches the filters.", file=out)
    if af is not None:
        print(
            f"AlphaFold: AF-{af.accession}-F1 v{af.version}, mean pLDDT {af.mean_plddt:.1f}, "
            f"{af.fraction_very_high:.0%} of residues very high confidence",
            file=out,
        )


def _filters(args: argparse.Namespace) -> Filters:
    return Filters(
        max_resolution=args.max_resolution,
        min_coverage=args.min_coverage,
        method=args.method,
        ligand=args.ligand,
        require_ligand=args.holo,
    )


def rank(
    client: Client, accession: str, filters: Filters, top: int | None = None
) -> list[Structure]:
    """Filter on the cheap fields first, then look up ligands only for what's left."""
    ligand_free = replace(filters, ligand=None, require_ligand=False)
    ranked = select_structures(client.structures(accession), ligand_free)
    if filters.ligand or filters.require_ligand:
        ranked = select_structures(client.with_ligands(ranked), filters)
        return ranked[:top]
    return client.with_ligands(ranked[:top])


def cmd_list(client: Client, args: argparse.Namespace, out: TextIO) -> int:
    target = client.target(args.accession)
    ranked = rank(client, target.accession, _filters(args), top=args.top)
    af = client.alphafold(target.accession)
    rows = [_row(i, s) for i, s in enumerate(ranked, 1)]

    if args.format == "json":
        payload = {
            "target": asdict(target),
            "structures": [
                {**asdict(s), "bound_ligands": [lig.code for lig in s.bound_ligands]}
                for s in ranked
            ],
            "alphafold": asdict(af) if af else None,
        }
        json.dump(payload, out, indent=2)
        out.write("\n")
    elif args.format == "csv":
        writer = csv.DictWriter(out, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    else:
        _print_table(target, rows, af, out)
    return 0


def cmd_fetch(client: Client, args: argparse.Namespace, out: TextIO) -> int:
    target = client.target(args.accession)
    ranked = rank(client, target.accession, _filters(args), top=1)
    args.output.mkdir(parents=True, exist_ok=True)

    if ranked and not args.alphafold:
        best = ranked[0]
        path = args.output / f"{best.pdb_id}.cif"
        path.write_bytes(client.download(entry_url(best.pdb_id)))
        print(f"{path}  (PDB {best.pdb_id} chain {best.chain}, {best.method})", file=out)
        return 0

    af = client.alphafold(target.accession)
    if af is None:
        print(f"No structure or AlphaFold model found for {target.accession}.", file=sys.stderr)
        return 1
    path = args.output / f"AF-{af.accession}-F1-model_v{af.version}.cif"
    path.write_bytes(client.download(af.cif_url))
    print(f"{path}  (AlphaFold, mean pLDDT {af.mean_plddt:.1f})", file=out)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="structpick",
        description=(
            "Rank PDB structures of a UniProt target for docking, with AlphaFold as fallback."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("accession", help="UniProt accession, e.g. P69905")
    common.add_argument(
        "--max-resolution", type=float, metavar="Å", help="drop entries worse than this"
    )
    common.add_argument("--min-coverage", type=float, default=0.0, metavar="FRACTION",
                        help="minimum fraction of the UniProt sequence covered (0-1)")  # fmt: skip
    common.add_argument("--method", help="xray, em, nmr (or the full PDBe method name)")
    common.add_argument(
        "--ligand", metavar="CODE", help="only entries with this ligand bound, e.g. HEM"
    )
    common.add_argument(
        "--holo", action="store_true", help="only entries with any non-additive ligand"
    )

    p_list = sub.add_parser("list", parents=[common], help="show ranked structures")
    p_list.add_argument("--top", type=int, default=10, help="how many to show (default 10)")
    p_list.add_argument("--format", choices=["table", "csv", "json"], default="table")
    p_list.set_defaults(func=cmd_list)

    p_fetch = sub.add_parser("fetch", parents=[common], help="download the best structure as mmCIF")
    p_fetch.add_argument("-o", "--output", type=Path, default=Path("."), help="output directory")
    p_fetch.add_argument(
        "--alphafold", action="store_true", help="download the AlphaFold model instead"
    )
    p_fetch.set_defaults(func=cmd_fetch)
    return parser


def main(
    argv: Sequence[str] | None = None, client: Client | None = None, out: TextIO | None = None
) -> int:
    args = build_parser().parse_args(argv)
    out = out or sys.stdout
    try:
        if client is not None:
            return int(args.func(client, args, out))
        with Client() as owned:
            return int(args.func(owned, args, out))
    except TargetNotFoundError as exc:
        print(f"structpick: {exc}", file=sys.stderr)
        return 2
    except httpx.HTTPError as exc:
        print(f"structpick: request failed: {exc}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
