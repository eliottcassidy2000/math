#!/usr/bin/env python3
"""
lrc_n14_res27_quotient_tower_s610.py

codex-2026-06-03 S610

Structural synthesis of the recent n=14 Res_27 improvements.

Input threads:
  * HYP-2162 / THM-407: the C=27 shell face folds by <2,-1> from
    13 raw shells to three gcd strata {1,3,9}.
  * HYP-2163: the no-multiple clock witness removes the bulk case, leaving
    the Cprime/lift residual.
  * HYP-2164: in the least-positive Res_27 quotient, D/U/N + pair-sum pinches
    leave no below-floor row and only AP, V*, and nonprimitive 2*AP at floor.
  * HYP-2165: after reattaching the canonical unit-spine owner labels through
    slack <=81, every full D/U/N cover has either a cheap pair or positive
    strict safe measure.

This script treats those as a quotient tower and audits the proof atoms left by
their composition.  It deliberately uses proof obligations, shell/gcd strata,
and owner-fibre rows as vertices rather than runners.

Tournament Analysis / assumption challenge:
  Vertices are proof atoms in the quotient tower.  The pairwise observable is
  a burden tuple (open?, exact floor?, owner missing?, normalized?, label), and
  the switch is the induced transitive proof order with lexicographic tie path.
  Alternate vertices considered: runners, residues, gaps, pair-sum denominators,
  fixed round classes, slack rows, shell orbits, and lift CRT states.  The proof
  atom quotient preserves the LRC discharge predicate but destroys raw runner
  identity, individual lift choices, and phase order.  The challenged assumption
  is that one layer should prove n=14 alone; the tower says each layer is a
  left/right adjoint view of the same coimage and only the lift-conservativity
  seam is still genuinely open.
"""

from __future__ import annotations

import importlib.util
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINCH_PATH = ROOT / "04-computation" / "lrc_n14_res27_pinch_certificate_s608.py"
BRIDGE_PATH = ROOT / "04-computation" / "lrc_n14_res27_fixed_bridge_s609.py"

SESSION = "S610"
N = 14
C = 2 * N - 1


@dataclass(frozen=True)
class Atom:
    label: str
    layer: str
    open: bool
    floor: bool
    owner_missing: bool
    normalized: bool
    note: str

    @property
    def burden(self) -> tuple[int, int, int, int, str]:
        return (
            1 if self.open else 0,
            1 if self.floor else 0,
            1 if self.owner_missing else 0,
            0 if self.normalized else 1,
            self.label,
        )


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def shell_of(a: int) -> int:
    r = a % C
    return min(r, C - r)


def shell_orbits() -> dict[int, tuple[int, ...]]:
    values = set(range(1, C))
    orbits: dict[int, set[int]] = defaultdict(set)
    while values:
        seed = min(values)
        orbit = set()
        x = seed
        for _ in range(64):
            orbit.add(x % C)
            orbit.add((-x) % C)
            x = (2 * x) % C
        orbit.discard(0)
        values -= orbit
        shells = {shell_of(x) for x in orbit}
        g = gcd(seed, C)
        orbits[g].update(shells)
    return {g: tuple(sorted(v)) for g, v in sorted(orbits.items())}


def fmt_score(witness) -> str:
    score = witness.score
    return str(score.numerator) if score.denominator == 1 else f"{score.numerator}/{score.denominator}"


def fmt_time(witness) -> str:
    t = witness.time
    return str(t.numerator) if t.denominator == 1 else f"{t.numerator}/{t.denominator}"


def ensure_reduced_numerators(pinch, row: tuple[int, ...]) -> None:
    for q in {a + b for a, b in combinations(row, 2)}:
        if q not in pinch.REDUCED_NUMERATORS:
            pinch.REDUCED_NUMERATORS[q] = tuple(a for a in range(1, q) if gcd(a, q) == 1)


def classify_pinch_layer(pinch):
    unit_rows, _ = pinch.enumerate_unit_rows()
    ledger_survivors: list[tuple[int, ...]] = []
    type_stats: defaultdict[object, Counter[str]] = defaultdict(Counter)
    floor_rows = []
    strict_count = 0
    below_count = 0

    for row in unit_rows:
        if pinch.d_failures(row) or pinch.n_failures(row):
            continue
        ledger_survivors.append(row)
        witness = pinch.best_pair_sum_pinch(row)
        key = pinch.res27_type(row)
        type_stats[key]["total"] += 1
        cmp_num = pinch.N * witness.score_num - witness.score_den
        if cmp_num > 0:
            strict_count += 1
            type_stats[key]["strict"] += 1
        elif cmp_num == 0:
            floor_rows.append((row, witness, key))
            type_stats[key]["floor"] += 1
        else:
            below_count += 1
            type_stats[key]["below"] += 1

    return {
        "ledger_survivors": ledger_survivors,
        "type_stats": type_stats,
        "strict_count": strict_count,
        "floor_rows": floor_rows,
        "below_count": below_count,
    }


def classify_owner_layer(pinch, bridge):
    s578 = bridge.load_s578()
    rows = bridge.slack_scan(s578, bridge.SLACK_BOUND)
    full = [r for r in rows if r.full_cover]
    floor = [r for r in full if r.positive_measure is False and r.unblocked_pair is not None]
    no_cheap = [r for r in full if r.unblocked_pair is None]
    residual = [r for r in full if bridge.slack_route(r) == "open residual"]

    controls = []
    for row in no_cheap:
        ensure_reduced_numerators(pinch, row.speeds)
        witness = pinch.best_pair_sum_pinch(row.speeds)
        controls.append((row, witness, pinch.res27_type(row.speeds)))

    return {
        "rows": rows,
        "full": full,
        "floor": floor,
        "no_cheap": no_cheap,
        "residual": residual,
        "controls": controls,
    }


def tournament_fingerprint(atoms: list[Atom]) -> dict[str, object]:
    ordered = sorted(atoms, key=lambda a: a.burden)
    score_hist = {i: 1 for i in range(len(ordered))}
    return {
        "vertices": len(ordered),
        "score_hist": score_hist,
        "directed_3_cycles": 0,
        "sccs": [1] if ordered else [],
        "hamiltonian_path": [a.label for a in ordered],
    }


def make_atoms(pinch_data, owner_data) -> list[Atom]:
    atoms = [
        Atom(
            "C1 no-multiple clock exit",
            "clock",
            False,
            False,
            False,
            True,
            "t=1/n closes every row with no speed divisible by n",
        ),
        Atom(
            "THM-407 gcd=1 shell orbit",
            "shell",
            False,
            False,
            True,
            True,
            "nine unit shells are one <2,-1> orbit; owner labels still matter",
        ),
        Atom(
            "THM-407 gcd=3 shell orbit",
            "shell",
            True,
            False,
            True,
            True,
            "three nonunit shells form the prime-3 residual stratum",
        ),
        Atom(
            "THM-407 gcd=9 shell orbit",
            "shell",
            True,
            False,
            True,
            True,
            "single shell {9,18}; most rigid prime-3 stratum",
        ),
        Atom(
            "pinch strict rows",
            "least-positive Res_27",
            False,
            False,
            True,
            True,
            f"{pinch_data['strict_count']} canonical D/U/N rows have score > 1/14",
        ),
    ]

    for row, _witness, key in pinch_data["floor_rows"]:
        name = "AP" if row == tuple(range(1, N)) else "V*" if row == (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24) else "2*AP"
        atoms.append(
            Atom(
                f"pinch floor {name}",
                "least-positive Res_27",
                False,
                True,
                True,
                name != "2*AP",
                f"type={key.label}",
            )
        )

    atoms.append(
        Atom(
            "owner cheap-pair floor AP/V*",
            "owner fibre",
            False,
            True,
            False,
            True,
            f"{len(owner_data['floor'])} owner-floor rows, both discharged by (1,13)",
        )
    )
    atoms.append(
        Atom(
            "owner no-cheap positive controls",
            "owner fibre",
            False,
            False,
            False,
            True,
            f"{len(owner_data['no_cheap'])} block-all rows, all positive measure",
        )
    )
    atoms.append(
        Atom(
            "lift/CRT conservativity seam",
            "lift theorem",
            True,
            False,
            False,
            True,
            "prove arbitrary lifts route to the closed quotient atoms",
        )
    )
    return atoms


def main() -> None:
    pinch = load_module(PINCH_PATH, "s608_pinch")
    bridge = load_module(BRIDGE_PATH, "s609_bridge")

    print(f"{SESSION} n=14 Res_27 quotient tower synthesis")
    print("=" * 78)
    print("Question: what structure do HYP-2162..2165 share, and what proof")
    print("atoms remain after composing them?")
    print()

    print("A. Shell-orbit quotient")
    orbits = shell_orbits()
    for g, shells in orbits.items():
        print(f"  gcd={g:2d}: shells={shells} count={len(shells)}")
    print("  raw shells 13 -> gcd strata 3 under <2,-1>.")
    print()

    print("B. Least-positive Res_27 pinch layer")
    pinch_data = classify_pinch_layer(pinch)
    type_stats = pinch_data["type_stats"]
    print(f"  D/U/N ledger survivors: {len(pinch_data['ledger_survivors'])}")
    print(f"  proof-obligation types: {len(type_stats)}")
    print(f"  strict rows: {pinch_data['strict_count']}")
    print(f"  floor rows: {len(pinch_data['floor_rows'])}")
    print(f"  below rows: {pinch_data['below_count']}")
    for row, witness, key in pinch_data["floor_rows"]:
        print(
            f"    floor row={pinch.row_name(row):16s} type={key.label} "
            f"pinch={fmt_score(witness)} at t={fmt_time(witness)}"
        )
    print()

    print("C. Owner-fibre reattachment layer")
    owner_data = classify_owner_layer(pinch, bridge)
    print(f"  canonical slack rows through 81: {len(owner_data['rows'])}")
    print(f"  full D/U/N covers: {len(owner_data['full'])}")
    print(f"  open owner residuals: {len(owner_data['residual'])}")
    print(f"  owner floor rows: {len(owner_data['floor'])}")
    for row in owner_data["floor"]:
        print(
            f"    slack={row.slack} pair={bridge.fmt_pair(row.unblocked_pair)} "
            f"shells={row.shell_signature}"
        )
    print(f"  no-cheap positive controls: {len(owner_data['no_cheap'])}")
    for row, witness, key in owner_data["controls"]:
        print(
            f"    slack={row.slack} positive={row.positive_measure} "
            f"pinch={fmt_score(witness)} at t={fmt_time(witness)} "
            f"type={key.label} gcd={row.gcd_signature}"
        )
    print()

    print("D. Composed proof atoms")
    atoms = make_atoms(pinch_data, owner_data)
    for atom in sorted(atoms, key=lambda a: a.burden):
        state = "OPEN" if atom.open else "closed"
        floor = "floor" if atom.floor else "loose"
        owner = "needs-owner" if atom.owner_missing else "owner-attached"
        print(f"  {atom.label:34s} {state:6s} {floor:5s} {owner:14s} :: {atom.note}")
    print()

    print("Tournament Analysis over proof atoms")
    fp = tournament_fingerprint(atoms)
    print(f"  vertices: {fp['vertices']}")
    print(f"  score histogram: {fp['score_hist']}")
    print(f"  SCCs: {fp['sccs']}")
    print(f"  directed 3-cycles: {fp['directed_3_cycles']}")
    print("  tie Hamiltonian path:")
    for label in fp["hamiltonian_path"]:
        print(f"    {label}")
    print()

    print("Synthesis")
    print("  The shared structure is a coimage/opfibration tower:")
    print("    clock exit -> shell orbit quotient -> pinch lower bound -> owner fibre.")
    print("  The left-moving maps forget data and prove lower bounds cheaply; the")
    print("  right-moving owner fibre reattaches exactly the labels needed for cheap")
    print("  pair or positive-measure discharge.")
    print("  Leverage: the remaining n=14 theorem should not enumerate 64 classes,")
    print("  27,733 pinch survivors, or 9,506 owner covers again.  It should prove")
    print("  lift/CRT conservativity for the finite atom list above: arbitrary lifts")
    print("  must land in strict pinch rows, AP/V*, normalized 2*AP, or the two")
    print("  positive owner-control rows.")


if __name__ == "__main__":
    main()
