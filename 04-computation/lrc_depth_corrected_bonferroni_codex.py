#!/usr/bin/env python3
"""Exact audit for the activation-depth Bonferroni obligation.

The corrected Legendre/Venn picture suggests truncating the far-runner Newton
expansion after order 3.  This script records the guardrail: the order-3
truncation is not universal.  It becomes the right target only after the first
live far-packet layer is at depth at most 2; triple-active rows need an
order-4 target instead.
"""

from __future__ import annotations

from collections import Counter
from functools import cache
from itertools import combinations
from fractions import Fraction


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


@cache
def p0_exact(E_tuple: tuple[int, ...]) -> Fraction:
    """Exact measure that the nonzero speeds hit all six nonzero sectors."""
    E = tuple(sorted(set(e for e in E_tuple if e != 0)))
    if not E:
        return Fraction(0)

    breaks = {Fraction(0), Fraction(1)}
    for e in E:
        for j in range(8):
            m = 0
            while True:
                x = (Fraction(j, 7) + m) / e
                if x >= 1:
                    break
                breaks.add(x)
                m += 1

    total = Fraction(0)
    for lo, hi in zip(sorted(breaks), sorted(breaks)[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        sectors = set()
        for e in E:
            r = frac_part(e * mid)
            sectors.add((7 * r.numerator) // r.denominator)
        if set(range(1, 7)).issubset(sectors):
            total += hi - lo
    return total


def packet_totals(base: tuple[int, ...], far: tuple[int, ...]) -> tuple[Fraction, Fraction, list[Fraction]]:
    """Return p0(base), p0(base+far), and total Newton packets T_1..T_m."""
    n = len(far)
    values: dict[tuple[int, ...], Fraction] = {}
    for r in range(n + 1):
        for S in combinations(range(n), r):
            values[S] = p0_exact(tuple(sorted(base + tuple(far[i] for i in S))))

    totals: list[Fraction] = []
    for r in range(1, n + 1):
        packet_total = Fraction(0)
        for S in combinations(range(n), r):
            mixed = Fraction(0)
            for j in range(r + 1):
                for T in combinations(S, j):
                    mixed += ((-1) ** (r - j)) * values[T]
            packet_total += mixed
        totals.append(packet_total)
    return values[()], values[tuple(range(n))], totals


def first_live_depth(totals: list[Fraction]) -> int | None:
    for index, value in enumerate(totals, 1):
        if value:
            return index
    return None


def fmt(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


ROWS = [
    (
        "triple_active_counterexample",
        (0, 1, 2, 3),
        (16, 19, 22, 25, 28),
        "T1=T2=0, so order 3 is too early.",
    ),
    (
        "four_far_triple_active_edge",
        (0, 1, 2, 3),
        (12, 13, 14, 15),
        "Same activation defect with four far runners.",
    ),
    (
        "spread_far_high_depth_slack",
        (0, 3, 9),
        (27, 28, 42, 43, 47),
        "KPS S31u high-depth failure: Bonferroni-3 is zero, but p0 is far below cap.",
    ),
    (
        "edge_active_consecutive_base",
        (0, 1, 2, 3, 4),
        (20, 21, 22, 23, 24),
        "S31t-style edge-active row; order 3 is an upper bound.",
    ),
    (
        "edge_active_even_base",
        (0, 2, 4, 6, 8),
        (25, 27, 29, 31, 33),
        "Dilated/even base; order 3 remains an upper bound.",
    ),
    (
        "corner_edge_active_consecutive6",
        (0, 1, 2, 3, 4, 5),
        (30, 31, 32, 33, 34),
        "T1 is live, but the tail after order 3 is still nonpositive.",
    ),
]


def row_audit() -> list[dict[str, object]]:
    reports = []
    print("Activation-depth Bonferroni audit for LRC far packets")
    print("=" * 78)
    for name, base, far, note in ROWS:
        p_base, p_full, totals = packet_totals(base, far)
        r0 = first_live_depth(totals)
        tail_ge4 = sum(totals[3:], Fraction(0))
        tail_ge5 = sum(totals[4:], Fraction(0))
        order3_gap = p_base + sum(totals[:3], Fraction(0)) - p_full
        order4_gap = p_base + sum(totals[:4], Fraction(0)) - p_full
        report = {
            "name": name,
            "r0": r0,
            "tail_ge4": tail_ge4,
            "tail_ge5": tail_ge5,
            "order3_gap": order3_gap,
            "order4_gap": order4_gap,
        }
        reports.append(report)

        print(f"\n{name}")
        print(f"  base={base}")
        print(f"  far ={far}")
        print(f"  note={note}")
        print(f"  p0(base)={fmt(p_base)}")
        print(f"  p0(full)={fmt(p_full)}")
        print(f"  first live packet r0={r0}")
        for index, value in enumerate(totals, 1):
            print(f"  T_{index}={fmt(value)}")
        print(f"  tail T_>=4={fmt(tail_ge4)}")
        print(f"  tail T_>=5={fmt(tail_ge5)}")
        print(f"  order-3 gap (trunc-full)={fmt(order3_gap)}")
        print(f"  order-4 gap (trunc-full)={fmt(order4_gap)}")
    return reports


def print_synthesis(reports: list[dict[str, object]]) -> None:
    print("\nCorrected obligation")
    print("=" * 78)
    print("The universal Bonferroni-3 target is false.")
    bad = next(r for r in reports if r["name"] == "triple_active_counterexample")
    print(
        "  counterexample: r0=3 and order-3 gap="
        f"{fmt(bad['order3_gap'])}<0, while order-4 gap={fmt(bad['order4_gap'])}>0."
    )
    print()
    print("Sharp target:")
    print("  if r0<=2, prove T_>=4<=0 and close with order 3;")
    print("  if r0=3, prove T_>=5<=0 and close with order 4;")
    print("  if r0>=4, prove a slack bound p0<<cap or route to finite residual atlases.")
    print()
    print("This keeps the Legendre/Venn layer labels but adds the missing")
    print("activation-depth address before scalarizing the far-packet expansion.")


def print_tournament_analysis(reports: list[dict[str, object]]) -> None:
    print("\nTournament Analysis over proof obligations")
    print("=" * 78)
    vertices = {
        "activation_depth": {"preserves_predicate", "r0", "tail_sign", "mode_labels"},
        "order3_edge_active": {"r0<=2", "tail_ge4", "cap_target", "mode_labels"},
        "order4_triple_active": {"r0=3", "tail_ge5", "cap_target", "mode_labels"},
        "higher_depth_residual": {"r0>=4", "residual_atlas", "decorrelation"},
        "raw_bonferroni3": {"tail_ge4"},
        "raw_runner_vertices": set(),
    }

    names = list(vertices)
    scores = Counter({name: 0 for name in names})
    adj = {name: set() for name in names}
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            key_a = (len(vertices[a]), "preserves_predicate" in vertices[a], -i)
            key_b = (len(vertices[b]), "preserves_predicate" in vertices[b], -j)
            if key_a >= key_b:
                adj[a].add(b)
                scores[a] += 1
            else:
                adj[b].add(a)
                scores[b] += 1

    cycles3 = 0
    for a, b, c in combinations(names, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            cycles3 += 1
        if c in adj[a] and b in adj[c] and a in adj[b]:
            cycles3 += 1
    path = sorted(names, key=lambda name: (scores[name], len(vertices[name]), name), reverse=True)

    print("  vertices are proof obligations, not runners")
    print("  observable=(LRC predicate retained, first live layer, tail sign)")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={cycles3}")
    print("  Hamiltonian path=" + " > ".join(path))
    print()
    failures = [r["name"] for r in reports if r["order3_gap"] < 0]
    print(f"  raw_bonferroni3 failures in audited rows={failures}")
    print("  challenged assumption: all wide rows share the same truncation level.")


def main() -> None:
    reports = row_audit()
    print_synthesis(reports)
    print_tournament_analysis(reports)


if __name__ == "__main__":
    main()
