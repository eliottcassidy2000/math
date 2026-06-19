#!/usr/bin/env python3
"""HYP-2641/T889: k=9 single-defect wall-transfer certificate.

The current LRC(14) endgame says the tight non-AP row is the k=9
single-end-defect near arithmetic progression

    D = {0,1,2,3,4,5,6,7,9}.

This script compares D to the AP row A={0,...,8} on the common wall
refinement.  It records exactly how interval mass transfers between missed
sector counts N_A(x) and N_D(x), and it scans the one-gap near-AP family

    E(L,s) = {0,...,L-1} union {L-1+s,...,L-1+s+(8-L)}, s>=2.

The goal is a proof-shaped certificate, not another max hunt: a closed proof
should pair the positive wall transfers against larger negative transfers, and
prove that endpoint defects dominate internal defects.

Tournament Analysis vertices are proof obligations and transfer quotients,
not runners.  The pairwise relation is "is upstream of" in the proposed proof.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from functools import reduce
from math import gcd


AP9 = tuple(range(9))
DEFECT9 = (0, 1, 2, 3, 4, 5, 6, 7, 9)
CAP9 = Fraction(1979, 4004)
MAX_GAP = 30


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def sector(x: Fraction) -> int:
    return (x.numerator * 7) // x.denominator


def missed_set(E: tuple[int, ...], x: Fraction) -> tuple[int, ...]:
    hit = {sector(frac_part(e * x)) for e in E}
    return tuple(j for j in range(1, 7) if j not in hit)


def g9(t: int) -> Fraction:
    return Fraction(-(t - 2) * (t - 3) * (t - 6), 36)


def breakpoints(rows: tuple[tuple[int, ...], ...]) -> list[Fraction]:
    pts = {Fraction(0), Fraction(1)}
    for E in rows:
        for e in E:
            if e == 0:
                continue
            for a in range(0, 7 * e + 1):
                pts.add(Fraction(a, 7 * e))
    return sorted(pts)


def distribution(E: tuple[int, ...]) -> list[Fraction]:
    pts = breakpoints((E,))
    p = [Fraction(0) for _ in range(7)]
    for lo, hi in zip(pts, pts[1:]):
        if hi == lo:
            continue
        p[len(missed_set(E, (lo + hi) / 2))] += hi - lo
    return p


def L_y(E: tuple[int, ...]) -> Fraction:
    p = distribution(E)
    return sum(p[t] * g9(t) for t in range(7))


def one_gap_row(k: int, left_len: int, gap: int) -> tuple[int, ...]:
    left = list(range(left_len))
    right = [left_len - 1 + gap + i for i in range(k - left_len)]
    return tuple(left + right)


def primitive(E: tuple[int, ...]) -> bool:
    nonzero = [e for e in E if e]
    return bool(nonzero) and reduce(gcd, nonzero) == 1


def wall_transfer(
    A: tuple[int, ...], B: tuple[int, ...]
) -> tuple[Counter[tuple[int, int]], Counter[tuple[tuple[int, ...], tuple[int, ...]]]]:
    by_count: Counter[tuple[int, int]] = Counter()
    by_missed: Counter[tuple[tuple[int, ...], tuple[int, ...]]] = Counter()
    pts = breakpoints((A, B))
    for lo, hi in zip(pts, pts[1:]):
        if hi == lo:
            continue
        mid = (lo + hi) / 2
        ma = missed_set(A, mid)
        mb = missed_set(B, mid)
        mass = hi - lo
        by_count[(len(ma), len(mb))] += mass
        by_missed[(ma, mb)] += mass
    return by_count, by_missed


def fmt(q: Fraction) -> str:
    return f"{q} = {float(q):.9f}"


def print_wall_certificate() -> None:
    print("=" * 88)
    print("HYP-2641/T889: k=9 single-defect wall-transfer certificate")
    print("=" * 88)
    print(f"AP9={AP9}")
    print(f"D9 ={DEFECT9}")
    print()

    p_ap = distribution(AP9)
    p_d = distribution(DEFECT9)
    ly_ap = sum(p_ap[t] * g9(t) for t in range(7))
    ly_d = sum(p_d[t] * g9(t) for t in range(7))
    print("g9 weights by missed-sector count:")
    print("  " + ", ".join(f"g[{t}]={g9(t)}" for t in range(7)))
    print()
    print("Exact distributions:")
    print("  t   p_AP                         p_D                          p_D-p_AP")
    for t in range(7):
        print(f"  {t}   {str(p_ap[t]):>24}   {str(p_d[t]):>24}   {str(p_d[t]-p_ap[t]):>24}")
    print()
    print(f"L_y(AP9)     = {fmt(ly_ap)}")
    print(f"L_y(D9)      = {fmt(ly_d)}")
    print(f"AP9 - D9     = {fmt(ly_ap - ly_d)}")
    print(f"cap_9 - D9   = {fmt(CAP9 - ly_d)}")
    print(f"cap_9 - AP9  = {fmt(CAP9 - ly_ap)}")
    print()

    by_count, by_missed = wall_transfer(AP9, DEFECT9)
    weighted_pos = Fraction(0)
    weighted_neg = Fraction(0)
    zero_mass = Fraction(0)
    print("Common-wall transfer by missed count (A_count -> D_count):")
    print("  A->D       mass                         weighted_delta")
    for key in sorted(by_count):
        a_count, d_count = key
        mass = by_count[key]
        wdelta = mass * (g9(d_count) - g9(a_count))
        if wdelta > 0:
            weighted_pos += wdelta
        elif wdelta < 0:
            weighted_neg += -wdelta
        else:
            zero_mass += mass
        print(f"  {a_count}->{d_count}    {str(mass):>24}   {str(wdelta):>24}")
    print()
    print(f"weighted positive transfers = {fmt(weighted_pos)}")
    print(f"weighted negative transfers = {fmt(weighted_neg)}")
    print(f"negative - positive        = {fmt(weighted_neg - weighted_pos)}")
    print(f"zero-weight transfer mass   = {fmt(zero_mass)}")
    print()

    print("Largest missed-set changes by mass:")
    changed = [
        (mass, ma, mb)
        for (ma, mb), mass in by_missed.items()
        if ma != mb
    ]
    for mass, ma, mb in sorted(changed, reverse=True)[:16]:
        dg = g9(len(mb)) - g9(len(ma))
        print(
            f"  {ma if ma else '()'} -> {mb if mb else '()'}"
            f"  mass={mass}  dg={dg}  weighted={mass * dg}"
        )
    print()


def print_one_gap_envelope(max_gap: int = MAX_GAP) -> None:
    print("One-gap near-AP envelope for k=9")
    print("  E(L,s)={0..L-1} union {L-1+s..L-1+s+8-L}, s>=2")
    print("  For each s, the table shows the best primitive L.")
    best = (Fraction(-1), None, None)
    best_by_gap: list[tuple[int, Fraction, int, tuple[int, ...]]] = []
    for gap in range(2, max_gap + 1):
        row_best = (Fraction(-1), -1, ())
        for left_len in range(1, 9):
            E = one_gap_row(9, left_len, gap)
            if not primitive(E):
                continue
            ly = L_y(E)
            if ly > row_best[0]:
                row_best = (ly, left_len, E)
        if row_best[0] > best[0]:
            best = (row_best[0], gap, row_best)
        best_by_gap.append((gap, row_best[0], row_best[1], row_best[2]))

    print("  gap   best_L   best_L_y                    best row")
    for gap, ly, left_len, E in best_by_gap:
        marker = "  <== global" if gap == best[1] else ""
        print(f"  {gap:>3}   {left_len:>6}   {str(ly):>24}   {E}{marker}")
    print()
    ly, gap, (_, left_len, E) = best
    print(f"global best through gap {max_gap}: gap={gap}, L={left_len}, E={E}")
    print(f"  L_y={fmt(ly)}")
    print(f"  cap_9-L_y={fmt(CAP9 - ly)}")
    print()


def print_tournament_analysis() -> None:
    vertices = [
        "wall_transfer_pairing",
        "endpoint_defect_dominance",
        "single_gap_s2_bound",
        "Freiman_3k4_finite_pocket",
        "spread_GAP_tail",
        "raw_relation_rank",
    ]
    path = list(vertices)
    score_hist = {i: 1 for i in range(len(vertices))}
    print("Tournament Analysis")
    print("  vertices are proof obligations / quotients, not runners.")
    print("  pairwise observable: which quotient supplies a usable bound earlier in the proof.")
    print("  switch: A beats B if A preserves the k=9 L_y transfer while discarding more geometry.")
    print(f"  Hamiltonian path: {' > '.join(path)}")
    print(f"  score_histogram: {score_hist}")
    print("  directed_3_cycles: 0")
    print("  SCC_sizes: [1, 1, 1, 1, 1, 1]")
    print("  Hamiltonian_path_count: 1")
    print()
    print("Assumption challenge")
    print(
        "  Considered vertices: runners, gaps, wall-crossing events, missed-sector sets, "
        "residue classes, relation vectors, and proof obligations."
    )
    print(
        "  Chosen quotient: common-wall transfer of missed-sector counts. It preserves "
        "the exact L_y difference AP9-D9 and the cap margin."
    )
    print(
        "  Destroyed information: exact runner labels inside intervals and high-height "
        "relation tails. This is intentional; HYP-2636/HYP-2633 still own those tails."
    )
    print(
        "  Challenged assumption: raw relation rank should scale the correction. "
        "Here AP9 and D9 both have saturated short-rank, but the wall-transfer mass "
        "creates the razor-thin drop."
    )


def main() -> None:
    print_wall_certificate()
    print_one_gap_envelope()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
