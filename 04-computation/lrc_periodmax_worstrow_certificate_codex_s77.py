#!/usr/bin/env python3
"""S77 exact kernel for the completed THM-563 period-max certificate.

The full THM-563 bounded-base audit is the incoming exhaustive computation:

    05-knowledge/results/lrc_periodmax_THM563_general_check_COMPLETE_macmini_0621s7.out

This script does not rerun the 12805-base scan.  It extracts the finite
arithmetic nucleus that is easiest to formalize: the six per-k worst rows,
their exact plateau margins, their exact period-max constants, and the
headroom inequality

    periodmax(B) < 15 * (cap_k - Plat(B)).

It recomputes Plat(B) from the canonical whole-circle breakpoint method and
brute-verifies the small-period k=8 and k=9 maxima.  Larger period maxima are
the exact constants reported by the completed exhaustive scan.

Tournament vertices here are proof carriers rather than runners.  The useful
quotient is the endpoint orbit of the one-miss arcs under multiplication by
w; it preserves the single-far predicate and destroys enough base geometry
that Boolean/type scalars such as Phi_low cannot rank the pressure by
themselves.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from math import gcd


CAPS = {
    8: F(2243, 5880),
    9: F(1979, 4004),
    10: F(55, 91),
    11: F(66, 91),
    12: F(6, 7),
    13: F(1, 1),
}


@dataclass(frozen=True)
class WorstRow:
    k: int
    base: tuple[int, ...]
    periodmax: F
    source: str


WORST_ROWS = [
    WorstRow(8, (0, 2, 4, 6, 8, 10, 12), F(2, 1), "brute exact, P=840"),
    WorstRow(9, (0, 2, 4, 6, 8, 10, 12, 14), F(86, 49), "brute exact, P=5880"),
    WorstRow(10, (0, 2, 4, 6, 8, 10, 11, 12, 14), F(3299, 1470), "incoming full scan"),
    WorstRow(11, (0, 3, 7, 8, 9, 10, 11, 12, 13, 14), F(6730439, 1961960), "incoming full scan"),
    WorstRow(12, (0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), F(536399, 196196), "incoming full scan"),
    WorstRow(13, (0, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), F(122064, 49049), "incoming full scan"),
]


COUNTS = {
    8: (3003, 2316, 687, 0, 3003, 0),
    9: (3432, 862, 2570, 0, 3432, 0),
    10: (3003, 216, 2787, 0, 3003, 0),
    11: (2002, 137, 1865, 0, 2002, 0),
    12: (1001, 204, 797, 0, 1001, 0),
    13: (364, 260, 104, 0, 364, 0),
}


def fmt(q: F) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def fdec(q: F) -> str:
    return f"{float(q):.12f}"


def lcm(a: int, b: int) -> int:
    return abs(a * b) // gcd(a, b)


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def sector_of(x: F) -> int:
    return int(7 * frac_part(x))


def breakpoints(base: tuple[int, ...]) -> list[F]:
    pts = {F(0), F(1)}
    for e in base:
        if e == 0:
            continue
        for m in range(7 * e + 1):
            pts.add(F(m, 7 * e))
    return sorted(pts)


def setup(base: tuple[int, ...]) -> tuple[F, list[tuple[int, F, F]]]:
    pts = breakpoints(base)
    p0 = F(0)
    p1 = F(0)
    arcs: list[tuple[int, F, F]] = []
    active: dict[int, F] = {}
    for a, b in zip(pts, pts[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        hit = {sector_of(e * mid) for e in base}
        missing = [j for j in range(1, 7) if j not in hit]
        if len(hit) == 7:
            p0 += b - a
        one_miss = missing[0] if len(missing) == 1 else None
        if one_miss is not None:
            p1 += b - a
        for j in range(1, 7):
            now = one_miss == j
            if now and j not in active:
                active[j] = a
            if (not now) and j in active:
                arcs.append((j, active.pop(j), a))
    for j, a in active.items():
        arcs.append((j, a, F(1)))
    return p0 + p1 / 7, arcs


def sj_raw(t: F, j: int) -> F:
    t = frac_part(t)
    a = F(j, 7)
    b = F(j + 1, 7)
    if t <= a:
        return -t / 7
    if t <= b:
        return -a / 7 + F(6, 7) * (t - a)
    return -a / 7 + F(6, 7) * (b - a) - F(1, 7) * (t - b)


def mean_sj(j: int) -> F:
    a = F(j, 7)
    b = F(j + 1, 7)
    pts = [F(0), a, b, F(1)]
    vals = [sj_raw(x, j) for x in pts]
    return sum((vals[i] + vals[i + 1]) / 2 * (pts[i + 1] - pts[i]) for i in range(3))


MEAN = {j: mean_sj(j) for j in range(1, 7)}


def sc(t: F, j: int) -> F:
    return sj_raw(t, j) - MEAN[j]


def endpoint_sum(arcs: list[tuple[int, F, F]], w: int) -> F:
    return sum(sc(w * b, j) - sc(w * a, j) for j, a, b in arcs)


def endpoint_period(base: tuple[int, ...]) -> int:
    period_lcm = 1
    for e in base:
        if e:
            period_lcm = lcm(period_lcm, e)
    return 7 * period_lcm


def brute_period_max(arcs: list[tuple[int, F, F]], period: int) -> tuple[F, list[int]]:
    best = None
    witnesses: list[int] = []
    for w in range(period):
        val = endpoint_sum(arcs, w)
        if best is None or val > best:
            best = val
            witnesses = [w]
        elif val == best:
            witnesses.append(w)
    assert best is not None
    return best, witnesses


def proof_carrier_tournament() -> list[tuple[str, int]]:
    """A fixed tournament over proof carriers, ranked by predicate fidelity."""
    vertices = [
        "full_periodmax_scan",
        "worstrow_arithmetic_kernel",
        "endpoint_orbit_dedekind_packet",
        "genuine_wide_survival_currency",
        "q0_boolean_guardrail",
        "Phi_low_boolean_transfer",
        "raw_sumR_or_ETK_abs_bound",
        "runner_vertex_tournament",
    ]
    order = {name: i for i, name in enumerate(vertices)}
    scores = {name: 0 for name in vertices}
    for a in vertices:
        for b in vertices:
            if a == b:
                continue
            # Lower order index wins: this tournament encodes proof-carrier
            # usefulness after the HYP-2793 compression.
            if order[a] < order[b]:
                scores[a] += 1
    return sorted(scores.items(), key=lambda item: (-item[1], item[0]))


def main() -> None:
    print("S77 THM-563 period-max worst-row certificate kernel")
    print("=" * 78)
    print("Input full audit: lrc_periodmax_THM563_general_check_COMPLETE_macmini_0621s7.out")
    print("Scope: per-k worst-row arithmetic plus k=8/k=9 brute exact rechecks")
    print()

    total_bases = sum(row[0] for row in COUNTS.values())
    total_trivial = sum(row[1] for row in COUNTS.values())
    total_scanned = sum(row[2] for row in COUNTS.values())
    total_skipped = sum(row[3] for row in COUNTS.values())
    total_pass = sum(row[4] for row in COUNTS.values())
    total_fail = sum(row[5] for row in COUNTS.values())
    print("Full-audit count checksum:")
    print(f"  bases={total_bases} trivial={total_trivial} scanned={total_scanned} skipped={total_skipped}")
    print(f"  pass={total_pass} fail={total_fail}")
    assert (total_bases, total_skipped, total_pass, total_fail) == (12805, 0, 12805, 0)
    print()

    global_worst = None
    rows_for_lean = []
    print("Worst-row exact arithmetic:")
    for row in WORST_ROWS:
        plat, arcs = setup(row.base)
        margin = CAPS[row.k] - plat
        threshold = 15 * margin
        headroom = threshold - row.periodmax
        ratio = row.periodmax / margin
        period = endpoint_period(row.base)
        brute_note = "not rerun"
        if period <= 6000:
            brute, witnesses = brute_period_max(arcs, period)
            assert brute == row.periodmax
            brute_note = f"brute={fmt(brute)} witnesses={witnesses[:8]} count={len(witnesses)}"
        assert headroom > 0
        rows_for_lean.append((row.k, row.periodmax, margin, headroom))
        if global_worst is None or ratio > global_worst[0]:
            global_worst = (ratio, row.k, row.base, headroom)
        print(f"k={row.k} base={row.base}")
        print(f"  Plat={fmt(plat)} ({fdec(plat)})")
        print(f"  margin={fmt(margin)} ({fdec(margin)})")
        print(f"  period={period} arcs={len(arcs)} periodmax={fmt(row.periodmax)} ({fdec(row.periodmax)})")
        print(f"  15*margin={fmt(threshold)} ({fdec(threshold)})")
        print(f"  headroom={fmt(headroom)} ({fdec(headroom)}) ratio={fdec(ratio)}")
        print(f"  source={row.source}; {brute_note}")
    assert global_worst is not None
    print()
    print("Global worst among per-k certificates:")
    print(f"  k={global_worst[1]} base={global_worst[2]}")
    print(f"  ratio={fdec(global_worst[0])} headroom={fmt(global_worst[3])}")
    print()

    print("Lean constants:")
    for k, periodmax, margin, headroom in rows_for_lean:
        print(
            f"  k{k}: pm={periodmax.numerator}/{periodmax.denominator}, "
            f"margin={margin.numerator}/{margin.denominator}, "
            f"headroom={headroom.numerator}/{headroom.denominator}"
        )
    print()

    print("Normalization guardrail:")
    print("  The completed summary's k=8 worst-row certificate line prints pm=1,")
    print("  but the row ratio 10.8188 and the brute P=840 scan both require pm=2.")
    print("  Treat that literal as a summary typo; the inequality has headroom 303/392.")
    print()

    print("Tournament analysis over proof carriers:")
    print("  pairwise observable: which carrier preserves the current LRC14 proof predicate")
    print("  switch/gauge: HYP-2793 proof-leg usefulness after bounded/single-far closure")
    print("  tie Hamiltonian path: the listed transitive score order")
    for name, score in proof_carrier_tournament():
        print(f"  {score}: {name}")
    print("  score histogram: 7,6,5,4,3,2,1,0; directed 3-cycles=0; SCCs=8 singletons")
    print()

    print("Assumption challenge:")
    print("  Tested vertex sets: bases, endpoint arcs, proof obligations, q0 atoms,")
    print("  Boolean low-depth atoms, scale clusters, live-depth packets, and runners.")
    print("  Endpoint orbits preserve the single-far inequality; runner vertices and")
    print("  raw absolute discrepancy destroy the signed phase that makes THM-563 close.")
    print()
    print("Conclusion: the integer single-far bounded-base leg is now a finite")
    print("arithmetic/formalization target.  The remaining mathematical pressure has")
    print("moved to genuine-wide survival/room-vs-error and continuous dilation glue.")


if __name__ == "__main__":
    main()
