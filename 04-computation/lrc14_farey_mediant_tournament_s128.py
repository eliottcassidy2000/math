#!/usr/bin/env python3
"""S128: Farey-mediant tournament interface for the LRC14 summit.

The exact arithmetic hook is this:

    M(S) = p/q, threshold = 1/14, excess e = 14*p - q.

Then M(S) > 1/14 iff e > 0, and e = 1 iff p/q is a Farey
neighbor of 1/14.  In that unit-excess case,

    p/q = p/(14*p - 1)

is the iterated Stern-Brocot mediant child above 1/14.  For example
2/27 = mediant(1/14, 1/13) and 3/41 = mediant(1/14, 2/27).

This script checks how that mediant arithmetic interacts with the tournament
quotient from S127.  It reports:

1. The apex winding tournament class at t = 1/14.
2. The exact optimum winding tournament class at t = M(S).
3. A row-level proof-priority tournament whose pairwise observable is escape
   height M(S)-1/14, with the listed row order as the tie Hamiltonian path.

The result is not an LRC14 proof.  It is a precise proof interface: a future
theorem can try to show every non-AP/GW apex-regular survivor must enter a
non-tight Farey-mediant class, hence has M(S) > 1/14.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s124 = load_module("s124_ap_gw", REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py")
s127 = load_module("s127_tournaments", REPO / "04-computation" / "lrc14_tournament_realizability_summit_codex_s127.py")

THR = F(1, 14)
AP = tuple(range(1, 14))
GW = tuple(list(range(1, 12)) + [13, 24])
RESIDUE_LIAR_26 = tuple(list(range(1, 12)) + [13, 26])
NEAR_MISS_36 = tuple(list(range(1, 12)) + [13, 36])
PETAL_8_16 = tuple(sorted((set(AP) - {8}) | {16}))
PETAL_10_20 = tuple(sorted((set(AP) - {10}) | {20}))


@dataclass(frozen=True)
class Row:
    label: str
    speeds: tuple[int, ...]
    family: str


def unique_rows() -> list[Row]:
    rows: list[Row] = [
        Row("AP", AP, "known tight"),
        Row("GW 12->24", GW, "known tight"),
        Row("residue-liar 12->26", RESIDUE_LIAR_26, "q-threshold loose"),
        Row("near-miss 12->36", NEAR_MISS_36, "Farey unit-excess loose"),
        Row("petal 8->16", PETAL_8_16, "minimal petal loose"),
        Row("petal 10->20", PETAL_10_20, "minimal petal loose"),
    ]

    for v in range(7, 14):
        new = 2 * v
        row = tuple(sorted((set(AP) - {v}) | {new}))
        if len(row) == 13 and new not in AP:
            rows.append(Row(f"petal {v}->{new}", row, "all minimal petals"))

    for m in range(1, 9):
        row = tuple(sorted(set(list(range(1, 12)) + [13, 12 * m])))
        if len(row) == 13:
            rows.append(Row(f"12m family m={m}", row, "12m tail"))

    seen: set[tuple[int, ...]] = set()
    out: list[Row] = []
    for row in rows:
        if row.speeds not in seen:
            seen.add(row.speeds)
            out.append(row)
    return out


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def winding_mask(row: tuple[int, ...], t: F) -> int:
    """Tournament mask on speeds sorted increasingly at scale t.

    Arc i->j iff frac((s_i-s_j)*t) lies in (0, 1/2).  Boundary ties at 0 or
    1/2 are oriented by the increasing-speed Hamiltonian path.
    """
    speeds = tuple(sorted(row))
    n = len(speeds)
    mask = 0
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            rel = frac_part(F(speeds[i] - speeds[j]) * t)
            if rel == 0 or rel == F(1, 2):
                orient_forward = True
            else:
                orient_forward = rel < F(1, 2)
            if orient_forward:
                mask |= 1 << bit
            bit += 1
    return mask


def exact_M(row: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    M, pts = s124.M_exact(row)
    return M, tuple(sorted(pts))


def farey_excess(value: F) -> int:
    return 14 * value.numerator - value.denominator


def farey_mediant_label(value: F) -> str:
    p, q = value.numerator, value.denominator
    e = farey_excess(value)
    if value == THR:
        return "tight floor e=0"
    if e == 1:
        if p == 1:
            return "right Farey parent 1/13"
        prev = F(p - 1, 14 * (p - 1) - 1)
        med = F(1 + prev.numerator, 14 + prev.denominator)
        assert med == value
        return f"unit child: med(1/14,{prev})"
    return f"nonunit excess e={e}"


def class_groups(rows: list[Row], masks: list[int], prefix: str) -> tuple[list[int], dict[int, list[str]]]:
    ids = s127.classify_large_masks(masks, 13)
    groups: dict[int, list[str]] = {}
    for row, cid in zip(rows, ids, strict=True):
        groups.setdefault(cid, []).append(row.label)
    print(f"  {prefix} achieved class groups:")
    for cid, labels in sorted(groups.items()):
        print(f"    {prefix}{cid}: {', '.join(labels)}")
    return ids, groups


def row_risk_tournament(rows: list[Row], Ms: list[F]) -> int:
    """Row-level tournament: harder row beats easier row.

    Observable: exact escape height M(S)-1/14.  Smaller escape is harder.
    Tie Hamiltonian path: listed row order.
    """
    n = len(rows)
    mask = 0
    bit = 0
    keys = [(M - THR, idx) for idx, M in enumerate(Ms)]
    for i in range(n):
        for j in range(i + 1, n):
            if keys[i] <= keys[j]:
                mask |= 1 << bit
            bit += 1
    return mask


def print_farey_algebra() -> None:
    print("[Farey algebra around 1/14]")
    print("  threshold floor = 1/14")
    print("  right parent    = 1/13")
    for p in range(2, 7):
        child = F(p, 14 * p - 1)
        prev = F(p - 1, 14 * (p - 1) - 1)
        print(
            f"  child p={p}: {child} = mediant(1/14,{prev}), "
            f"gap={child - THR}, determinant={child.denominator - 14 * child.numerator}"
        )
    print()


def print_row_table(rows: list[Row]) -> None:
    Ms = [exact_M(row.speeds)[0] for row in rows]
    apex_masks = [winding_mask(row.speeds, THR) for row in rows]
    opt_masks = [winding_mask(row.speeds, M) for row, M in zip(rows, Ms, strict=True)]
    apex_ids, _ = class_groups(rows, apex_masks, "A")
    opt_ids, _ = class_groups(rows, opt_masks, "F")
    print()
    print("[Rows: exact M, Farey status, apex class, optimum class]")
    print(
        f"  {'row':24s} {'q':>2s} {'M':>7s} {'e':>3s} "
        f"{'apex':>5s} {'opt':>5s} {'c3':>4s} {'hp':>10s}  Farey status"
    )
    for row, M, aid, oid, opt_mask in zip(rows, Ms, apex_ids, opt_ids, opt_masks, strict=True):
        fp = s127.tournament_fingerprint(opt_mask, 13)
        print(
            f"  {row.label:24s} {s124.q_threshold(row.speeds):2d} {str(M):>7s} "
            f"{farey_excess(M):3d} A{aid:<4d} F{oid:<4d} "
            f"{fp['c3']:4d} {fp['hp']:10d}  {farey_mediant_label(M)}"
        )
    print()

    risk = row_risk_tournament(rows, Ms)
    fp = s127.tournament_fingerprint(risk, len(rows))
    print("[Row-level proof-priority tournament]")
    print("  vertices: listed candidate rows")
    print("  observable: exact escape height M(S)-1/14")
    print("  switch/gauge: smaller escape beats larger escape; ties follow listed order")
    print("  quotient preserves: proof priority by closeness to the LRC14 floor")
    print("  quotient destroys: runner geometry and exact winding class")
    print(
        f"  fingerprint: score_hist={fp['score_hist']} c3={fp['c3']} "
        f"scc={fp['scc']} hp={fp['hp']}"
    )
    print()


def print_proof_readout() -> None:
    print("[Proof readout]")
    print("  Unit excess e=1 is the exact Farey-neighbor condition.")
    print("  The near miss 12->36 has M=3/41, the second iterated mediant child")
    print("  above 1/14: 3/41 = mediant(1/14, 2/27).")
    print("  A tournament-only apex proof cannot work, because apex classes are")
    print("  magnitude-blind.  The plausible tournament proof must be multi-scale:")
    print("    apex class at 1/14 + first Farey-mediant escape class + q-threshold data.")
    print("  Local theorem target:")
    print("    every non-AP/GW q=14 survivor either has nonunit excess, or enters a")
    print("    unit-excess Farey child class outside the tight optimum classes.")
    print("  That would turn the mediant into a certificate M(S)>1/14, not just a")
    print("  numerical near miss.")


def main() -> None:
    print("S128 LRC14 FAREY-MEDIANT TOURNAMENT INTERFACE")
    print("=" * 78)
    print("[Assumption challenge]")
    print("  considered vertex sets: runners, residues, Farey fractions, row families,")
    print("    optimum witness times, and proof obligations.")
    print("  chosen carriers: runner winding tournaments at t=1/14 and t=M(S), plus")
    print("    a row-level proof-priority tournament by escape height.")
    print("  preserved: Farey excess, achieved winding class, and proof ordering.")
    print("  destroyed: continuous neighborhoods away from the chosen scales.")
    print()
    rows = unique_rows()
    print_farey_algebra()
    print_row_table(rows)
    print_proof_readout()


if __name__ == "__main__":
    main()
