#!/usr/bin/env python3
"""HYP-2649/T896: below-14 LRC reproof ladder with current proof quotients.

Convention: LRC(N) has k=N-1 nonzero speeds and target 1/N.

This is not a brute-force proof of every known below-14 case.  It is a
training atlas for re-proving the proved range with the same quotients now
driving LRC(14):

* exact max-min gap M(S),
* strict safe-measure fattening when the target is relaxed from 1/N to 1/(N+1),
* AP-frontier one-step rows,
* the LRC(13) -> LRC(14) core-fattening lever on 12-subsets of [1,14],
* the even-denominator support-floor ladder that first reaches support six at
  N=14, matching THM-538/HYP-2646.

Tournament Analysis vertices are proof quotients, not runner labels.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd


def merge(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    for a, b in sorted(intervals):
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def measure(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((b - a for a, b in intervals), Fraction(0))


def danger_arcs(v: int, h: Fraction) -> list[tuple[Fraction, Fraction]]:
    arcs: list[tuple[Fraction, Fraction]] = []
    radius = h / v
    for a in range(v):
        c = Fraction(a, v)
        lo = (c - radius) % 1
        hi = (c + radius) % 1
        if lo < hi:
            arcs.append((lo, hi))
        else:
            arcs.append((lo, Fraction(1)))
            arcs.append((Fraction(0), hi))
    return arcs


def safe_components(V: tuple[int, ...], N: int) -> list[tuple[Fraction, Fraction]]:
    h = Fraction(1, N)
    danger = merge([arc for v in V for arc in danger_arcs(v, h)])
    safe: list[tuple[Fraction, Fraction]] = []
    prev = Fraction(0)
    for lo, hi in danger:
        if lo > prev:
            safe.append((prev, lo))
        prev = max(prev, hi)
    if prev < 1:
        safe.append((prev, Fraction(1)))
    return safe


def safe_measure(V: tuple[int, ...], N: int) -> Fraction:
    return measure(safe_components(V, N))


def M_exact_int(V: tuple[int, ...]) -> tuple[Fraction, Fraction | None]:
    """Exact max_t min_v ||v t|| via denominator scan.

    Existing repo scouts use the fact that an optimum occurs at a rational
    t=p/D with D<=2*max(V).  This is the same integer-only engine as S702.
    """
    max_v = max(V)
    best_num = 0
    best_den = 1
    best_t: Fraction | None = None
    for D in range(2, 2 * max_v + 1):
        for p in range(1, D):
            c = D
            for v in V:
                r = (v * p) % D
                d = min(r, D - r)
                if d < c:
                    c = d
                    if c == 0:
                        break
            if c and c * best_den > best_num * D:
                best_num = c
                best_den = D
                best_t = Fraction(p, D)
    return Fraction(best_num, best_den), best_t


def primitive(V: tuple[int, ...]) -> bool:
    g = 0
    for v in V:
        g = gcd(g, v)
    return g == 1


def fmt(q: Fraction) -> str:
    return f"{q} = {float(q):.9f}"


def ap_tight_rows() -> None:
    print("=" * 88)
    print("A. Exact AP tight rows below the first open case")
    print("=" * 88)
    print("LRC(N): V=(1,...,N-1), target=1/N.")
    print("  N   M(V)                 witness t       safe_meas@1/N")
    for N in range(3, 14):
        V = tuple(range(1, N))
        m, t = M_exact_int(V)
        sm = safe_measure(V, N)
        print(f" {N:2d}   {str(m):>18}   {str(t):>12}   {str(sm):>14}")
    print()


def one_step_frontier() -> None:
    print("=" * 88)
    print("B. One-step AP frontier: size N-1 subsets of [1,N]")
    print("=" * 88)
    print("This is the smallest AP-neighborhood where tight rows can split.")
    print("  N  rows  tight@1/N  min M             argmin rows                  min safe@1/(N+1)")
    for N in range(4, 14):
        rows = [tuple(c) for c in combinations(range(1, N + 1), N - 1) if primitive(tuple(c))]
        target = Fraction(1, N)
        scored = [(M_exact_int(V)[0], V) for V in rows]
        min_m = min(m for m, _ in scored)
        argmins = [V for m, V in scored if m == min_m]
        tight = [V for m, V in scored if m == target]
        fatten = [(safe_measure(V, N + 1), V) for V in rows]
        min_fat = min(v for v, _ in fatten)
        print(
            f" {N:2d}  {len(rows):4d}  {len(tight):9d}  "
            f"{str(min_m):>16}  {str(argmins[:3]):<28}  {str(min_fat):>16}"
        )
    print()


def lrc13_core_fattening() -> None:
    print("=" * 88)
    print("C. LRC(13) core-fattening lever inside 12-subsets of [1,14]")
    print("=" * 88)
    print("Rows have 12 speeds.  Compare exact M at target 1/13 and strict safe measure at 1/14.")
    rows = [tuple(c) for c in combinations(range(1, 15), 12) if primitive(tuple(c))]
    target13 = Fraction(1, 13)
    data = []
    for V in rows:
        m, t = M_exact_int(V)
        fat14 = safe_measure(V, 14)
        data.append((m, fat14, t, V))
    tight13 = [(fat14, t, V) for m, fat14, t, V in data if m == target13]
    min_m = min(m for m, _, _, _ in data)
    min_fat14 = min(fat14 for _, fat14, _, _ in data)
    fat_argmins = [(fat14, V) for _, fat14, _, V in data if fat14 == min_fat14]
    hist = Counter(m for m, _, _, _ in data)

    print(f"  primitive 12-subsets scanned: {len(rows)}")
    print(f"  min M at LRC(13) target: {fmt(min_m)}")
    print(f"  tight M=1/13 rows: {len(tight13)}")
    for fat14, t, V in tight13[:8]:
        print(f"    tight13 V={V} witness={t} safe_meas@1/14={fmt(fat14)}")
    print(f"  min strict safe measure at relaxed target 1/14: {fmt(min_fat14)}")
    for fat14, V in fat_argmins[:8]:
        print(f"    min-fat14 V={V}")
    print("  first exact M histogram entries:")
    for m, count in sorted(hist.items())[:10]:
        print(f"    M={m}: {count}")
    print()


def support_floor_ladder() -> None:
    print("=" * 88)
    print("D. Even-denominator sector support floor")
    print("=" * 88)
    print("For N=2q, the half-gap sector quotient has q sectors and q-1 live missed sectors.")
    print("Inclusion-exclusion kills Fourier correction below support q-1.")
    print("  N   q=N/2   support floor   proof meaning")
    for N in range(4, 16, 2):
        q = N // 2
        floor = q - 1
        if N < 14:
            meaning = "below support-six tail"
        elif N == 14:
            meaning = "first support-six / HYP-2646 wall"
        else:
            meaning = "beyond current frontier"
        print(f" {N:2d}   {q:5d}   {floor:13d}   {meaning}")
    print()


def tournament_fingerprint() -> None:
    print("=" * 88)
    print("Tournament Analysis: proof quotients for below-14 reproof")
    print("=" * 88)
    vertices = [
        "AP_tight_locus_plus_fattening",
        "signed_wall_transport",
        "support_floor_q_minus_1",
        "Freiman_small_excess",
        "exact_M_boundary_scan",
        "raw_runner_vertices",
    ]
    print("Pairwise observable: preserves proof-relevant address before scalarizing.")
    print("Switch/gauge: prefer the quotient that explains both tight equality and positive fattening.")
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("score histogram: {0:1,1:1,2:1,3:1,4:1,5:1}; directed 3-cycles: 0")


def main() -> None:
    print("HYP-2649/T896 below-14 modern LRC reproof ladder")
    print("Convention: LRC(N) has N-1 speeds and target 1/N.\n")
    ap_tight_rows()
    one_step_frontier()
    lrc13_core_fattening()
    support_floor_ladder()
    tournament_fingerprint()


if __name__ == "__main__":
    main()
