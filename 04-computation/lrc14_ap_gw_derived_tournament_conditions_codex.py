#!/usr/bin/env python3
"""Derived AP/GW necessary-condition atlas for LRC14.

This scout extends the older S124/S127 atlases.

The guiding analogy is Tao's "derived multiplicative functions": a derived
object is the infinitesimal coefficient of a rigid multiplicative identity.  In
LRC14 language, the Goddyn-Wong row is treated as a first-order tangent of AP:

    AP + epsilon * (e_10 - e_12),    epsilon^2 = 0.

The script keeps this metaphor exact by converting it into finite residue,
Farey, Jacobsthal, and tournament-isomorphism tests.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
S124_PATH = REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py"
S127_PATH = REPO / "04-computation" / "lrc14_tournament_realizability_summit_codex_s127.py"


def load_module(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


s124 = load_module(S124_PATH, "s124_ap_gw")
s127 = load_module(S127_PATH, "s127_tournament_summit")


N = 14
AP = s124.AP
GW = s124.GW
RESIDUE_LIAR_26 = s124.RESIDUE_LIAR_26
NEAR_MISS_36 = s124.NEAR_MISS_36
EVEN_SLOTS = {
    "A": frozenset((2, 12)),
    "B": frozenset((4, 10)),
    "C": frozenset((6, 8)),
}


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]


def replace_one(v: int, w: int) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - {v}) | {w}))


ROWS = [
    Row("AP", AP),
    Row("petal 8->16", replace_one(8, 16)),
    Row("petal 10->20", replace_one(10, 20)),
    Row("GW 12->24", GW),
    Row("residue-liar 12->26", RESIDUE_LIAR_26),
    Row("Farey child 12->36", NEAR_MISS_36),
    Row("odd petal 13->26", replace_one(13, 26)),
]


def residue_defect(row: tuple[int, ...]) -> tuple[int, ...]:
    c = s124.residue_counts(row)
    ap = s124.residue_counts(AP)
    return tuple(c[i] - ap[i] for i in range(N))


def slot_of_residue(r: int) -> str:
    r %= N
    for name, values in EVEN_SLOTS.items():
        if r in values:
            return name
    return "-"


def tangent_multiplier(row: tuple[int, ...]) -> tuple[str, int | None, int | None]:
    """Classify a row as zero, k*tangent, residue-vertical, or other."""
    if set(row) == set(AP):
        return ("zero", None, None)
    removed = sorted(set(AP) - set(row))
    added = sorted(set(row) - set(AP))
    if len(removed) != 1 or len(added) != 1:
        return ("other", None, None)
    v, w = removed[0], added[0]
    if w % N == v % N:
        return ("residue-vertical", v, None)
    if w % v == 0:
        return (f"{w // v}-tangent", v, w // v)
    return ("affine-slip", v, None)


def farey_index(m: Fraction) -> int:
    """Signed distance from the denominator-14 ray p/14p."""
    return m.denominator - N * m.numerator


def row_class_ids(rows: list[Row]) -> list[int]:
    masks = [s127.mask_from_runner_row(row.speeds) for row in rows]
    return s127.classify_large_masks(masks, 13)


def derived_tangent_or_ap(row: tuple[int, ...]) -> bool:
    kind, v, k = tangent_multiplier(row)
    if kind == "zero":
        return True
    return k == 2 and v in (8, 10, 12) and set(row) == set(replace_one(v, 2 * v))


def jacobsthal_tangent_gate_or_ap(row: tuple[int, ...]) -> bool:
    kind, v, k = tangent_multiplier(row)
    if kind == "zero":
        return True
    return k == 2 and v is not None and s124.gw_gate_passes(v)


def ap_gw_kind(row: tuple[int, ...]) -> bool:
    return (
        s124.q_punctured_cover(row)
        and s124.exact_unit_shell(row)
        and s124.exact_odd_skeleton(row)
        and s124.complement_sum14(row)
        and s124.zsum_max(row)
        and s124.single_even_dipole_or_zero(row)
        and derived_tangent_or_ap(row)
        and jacobsthal_tangent_gate_or_ap(row)
    )


def condition_tournament() -> None:
    bank = s124.bank_single_swaps(300)
    conditions = [
        ("q14_core", s124.q_punctured_cover),
        ("unit_shell", s124.exact_unit_shell),
        ("odd_skeleton", s124.exact_odd_skeleton),
        ("literal_unit_binders", s124.complement_sum14),
        ("zsum_max", s124.zsum_max),
        ("even_dipole_or_zero", s124.single_even_dipole_or_zero),
        ("dual_number_tangent", derived_tangent_or_ap),
        ("jacobsthal_gate", jacobsthal_tangent_gate_or_ap),
        ("ap_gw_kind", ap_gw_kind),
    ]
    pass_sets = [{row for row in bank if pred(row)} for _, pred in conditions]
    adj = [[False] * len(conditions) for _ in conditions]
    incomparable = 0
    flips = 0
    for i, j in combinations(range(len(conditions)), 2):
        a, b = pass_sets[i], pass_sets[j]
        if a < b:
            adj[i][j] = True
        elif b < a:
            adj[j][i] = True
            flips += 1
        elif a == b:
            adj[i][j] = True
        else:
            incomparable += 1
            if (len(a), i) <= (len(b), j):
                adj[i][j] = True
            else:
                adj[j][i] = True
                flips += 1
    scores = [sum(row) for row in adj]
    cyc3 = 0
    for a, b, c in combinations(range(len(conditions)), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cyc3 += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            cyc3 += 1
    print("[Tournament Analysis: derived condition strength]")
    print("  universe: AP single replacements v<=300")
    print("  relation: X -> Y iff X has a smaller pass set; ties by listed order")
    print("  condition pass sizes:")
    for (name, _), ps in zip(conditions, pass_sets, strict=True):
        print(f"    {name:<22} {len(ps):4d}")
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3cycles={cyc3}, incomparable_pairs={incomparable}, flips={flips}")
    print()


def print_row_table() -> None:
    class_ids = row_class_ids(ROWS)
    print("[Row diagnostics: AP, GW, and nearby impostors]")
    print(
        "  name                   class q   M      fidx defect_kind             "
        "tangent             slots  gate tight"
    )
    for row, cid in zip(ROWS, class_ids, strict=True):
        m, _ = s124.M_exact(row.speeds)
        defect_kind, missing, extra = s124.even_dipole_data(row.speeds)
        kind, v, k = tangent_multiplier(row.speeds)
        if v is None:
            slots = "--"
            gate = "-"
        else:
            src = slot_of_residue(v)
            tgt = slot_of_residue((k or 1) * v)
            slots = f"{src}->{tgt}"
            gate = str(s124.gw_gate_passes(v)) if k == 2 else "-"
        tight = "yes" if m == Fraction(1, N) else "no"
        print(
            f"  {row.name:<22} T{cid:<3d} {s124.q_threshold(row.speeds):<3d} "
            f"{str(m):<6} {farey_index(m):>4d} {defect_kind:<23} "
            f"{kind:<19} {slots:<6} {gate:<5} {tight}"
        )
    print()


def print_even_slot_cycle() -> None:
    print("[Even antipodal binder-slot cycle]")
    print("  slots: A={2,12}, B={4,10}, C={6,8}")
    print("  AP has one binder in each slot.  An even tail petal moves one binder.")
    print("  edge  row           M      Farey index  GW no-unit gate")
    for v in (8, 10, 12):
        row = replace_one(v, 2 * v)
        m, _ = s124.M_exact(row)
        src = slot_of_residue(v)
        tgt = slot_of_residue(2 * v)
        print(
            f"  {src}->{tgt}  {v:2d}->{2*v:<2d}       {str(m):<6} "
            f"{farey_index(m):>5d}        {s124.gw_gate_passes(v)}"
        )
    print("  Only A->B, the 12->24 edge, is tight at denominator 14.")
    print()


def print_kuratowski_mediants() -> None:
    print("[Kuratowski/Wagner obstruction mediants]")
    print("  Use vertex-per-edge density v/e, so disjoint union takes mediants.")
    print("  mK5+nK33 gives (5m+6n)/(10m+9n).")
    print("  m n  v/e     status")
    for m in range(0, 4):
        for n in range(0, 4):
            if not (m or n):
                continue
            frac = Fraction(5 * m + 6 * n, 10 * m + 9 * n)
            if (m, n) in ((1, 0), (0, 1)):
                status = "minimal obstruction atom"
            else:
                status = "mixture; nonplanar but not new minimal atom"
            print(f"  {m} {n}  {str(frac):<7} {status}")
    print("  LRC transfer: mediants mix obstruction packets; they do not by themselves")
    print("  create a new irreducible proof atom unless the labelled tournament class is new.")
    print()


def print_necessary_conditions() -> None:
    print("[Convoluted necessary-condition stack for AP/GW-kind tightness]")
    conditions = [
        "NC0 primitive cardinality: |S|=13 and gcd(S)=1.",
        "NC1 q-threshold exactness: q(S)=14, i.e. every q=2..13 is hit but 14 is not.",
        "NC2 unit/apex shell: unit residues and residue 7 occur exactly once.",
        "NC3 literal binders: the unique pairs in residues (1,13), (3,11), (5,9) sum to 14.",
        "NC4 antipodal derivative zero: the total mod-14 antipodal binder count remains 6.",
        "NC5 defect support: residue defect from AP is zero or a rank-one even-shell dipole.",
        "NC6 dual-number tangent: the nonzero dipole is e_{2v mod 14}-e_v from one tail v.",
        "NC7 binder-slot legality: the tangent edge lies in A,B,C with source->target defined by v->2v.",
        "NC8 Jacobsthal/GW gate: [14-v,27-2v] contains no integer coprime to v.",
        "NC9 Farey tight index: for exact M=p/q, q-14p=0; index -1 is a child/near-miss corridor.",
        "NC10 terminal tournament class: the augmented class is the AP class with q=14 or the GW class.",
        "NC11 quotient warning: raw tournament T0 also contains residue-liar 12->26, so q is external data.",
        "NC12 state-lift target: a bad residual must map to an achieved augmented class, not merely a graph.",
    ]
    for item in conditions:
        print("  " + item)
    print()


def main() -> None:
    print("LRC14 DERIVED AP/GW TOURNAMENT NECESSARY-CONDITION ATLAS")
    print("=" * 78)
    print("External analogies:")
    print("  Tao-derived multiplicative functions: first coefficients of rigid")
    print("    multiplicative identities over dual-number extensions.")
    print("  Collatz/Farey lower-bound style: a cycle must live in a rational")
    print("    approximant corridor; wrong Farey child means escape from the apex.")
    print()
    print_row_table()
    print_even_slot_cycle()
    print_kuratowski_mediants()
    condition_tournament()
    print_necessary_conditions()
    print("[Summit readout]")
    print("  AP is the zero tangent.  GW is the only single dual-number tangent that")
    print("  preserves q=14, the unit binders, the antipodal derivative, and the")
    print("  Jacobsthal no-unit window.  Other petal edges are achieved tournament")
    print("  classes, but they are loose or arithmetically mislabelled.  A proof by")
    print("  tournament isomorphism must therefore use augmented classes:")
    print("    (tournament iso class, q-threshold, Farey index, residue-defect tangent).")


if __name__ == "__main__":
    main()
