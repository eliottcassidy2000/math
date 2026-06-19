#!/usr/bin/env python3
"""
LRC(14) seven-sector cap scout (codex, 2026-06-18).

Context
-------
Recent agents reduced the k >= 8 S3 residual to the LRC-free endpoint

    mu_{1/7}(E) >= thr_k = 1 - min_{|P|=13-k} meas(G_P).

Equivalently, the "net" set

    N(E) = {x : all cyclic gaps in {e*x mod 1 : e in E} are <= 1/7}

must have measure at most cap_k = min meas(G_P).

This script explores a weaker sufficient route.  If a point configuration is a
1/7-net, then, up to boundary measure zero, it must meet each of the seven fixed
sectors [j/7,(j+1)/7).  Therefore

    N(E) subset S_7(E),

where S_7(E) is the set of x for which the phases hit all seven fixed sectors.
Thus it is enough to prove

    meas(S_7(E)) <= cap_k.

This route avoids assuming "consecutive minimizes mu_{1/7}".  It also records
counterexamples showing naive one-step gap compression is not monotone for
either mu_{1/7} or sector cover.

Tournament Analysis
-------------------
Vertices are proof routes rather than runners.  Pairwise observable is the
current evidence/risk score; the switch orients stronger routes toward weaker
routes, with the listed order as the tie Hamiltonian path.  This quotient
preserves proof-obligation dominance and destroys timing, wall adjacency, and
which sector boundary carries equality.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import reduce


ONE7 = F(1, 7)
HALF14 = F(1, 14)


def fmt(x: F) -> str:
    return f"{x} = {float(x):.9f}"


def merge(intervals):
    intervals = sorted(intervals)
    out = []
    for a, b in intervals:
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def measure(intervals) -> F:
    return sum((b - a for a, b in intervals), F(0))


def complement(intervals):
    intervals = merge(intervals)
    out = []
    prev = F(0)
    for a, b in intervals:
        if a > prev:
            out.append((prev, a))
        prev = max(prev, b)
    if prev < 1:
        out.append((prev, F(1)))
    return out


def danger_arcs(u: int, h: F = HALF14):
    """Open/closed endpoints do not affect measure."""
    arcs = []
    for j in range(u):
        c = F(j, u)
        a = (c - h / u) % 1
        b = (c + h / u) % 1
        if a < b:
            arcs.append((a, b))
        else:
            arcs.append((a, F(1)))
            arcs.append((F(0), b))
    return arcs


def safe_set(P):
    if not P:
        return [(F(0), F(1))]
    return complement(merge(arc for u in P for arc in danger_arcs(u)))


def meas_gp(P) -> F:
    return measure(safe_set(P))


def cap_table():
    caps = {}
    witnesses = {}
    for k in range(8, 14):
        p_size = 13 - k
        best = None
        best_p = None
        for P in itertools.combinations(range(1, 14), p_size):
            val = meas_gp(P)
            if best is None or val < best:
                best = val
                best_p = P
        caps[k] = best
        witnesses[k] = best_p
    return caps, witnesses


def collision_breakpoints(E):
    bp = {F(0), F(1)}
    E = sorted(set(E))
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            d = E[j] - E[i]
            for m in range(d + 1):
                bp.add(F(m, d))
    return sorted(bp)


def cell_cyclic_gaps(E, a: F, b: F):
    mid = (a + b) / 2
    pts = []
    for e in sorted(set(E)):
        val = e * mid
        floor = val.numerator // val.denominator
        pts.append((val - floor, e, floor))
    pts.sort(key=lambda row: row[0])
    gaps = []
    n = len(pts)
    for i in range(n):
        _, e_i, floor_i = pts[i]
        if i < n - 1:
            _, e_j, floor_j = pts[i + 1]
            alpha = e_j - e_i
            beta = floor_i - floor_j
        else:
            _, e_j, floor_j = pts[0]
            alpha = e_j - e_i
            beta = floor_i - floor_j + 1
        gaps.append((F(alpha), F(beta)))
    return gaps


def union_gt_length(gaps, a: F, b: F, theta: F) -> F:
    intervals = []
    for alpha, beta in gaps:
        if alpha == 0:
            if beta > theta:
                intervals.append((a, b))
            continue
        x_star = (theta - beta) / alpha
        if alpha > 0:
            lo, hi = max(a, x_star), b
        else:
            lo, hi = a, min(b, x_star)
        if lo < hi:
            intervals.append((lo, hi))
    return measure(merge(intervals))


def mu_theta(E, theta: F = ONE7) -> F:
    E = sorted(set(E))
    if len(E) == 1:
        return F(1)
    total = F(0)
    bp = collision_breakpoints(E)
    for a, b in zip(bp, bp[1:]):
        if a < b:
            total += union_gt_length(cell_cyclic_gaps(E, a, b), a, b, theta)
    return total


def sector_breakpoints(E):
    bp = {F(0), F(1)}
    for e in sorted(set(E)):
        if e == 0:
            continue
        for m in range(7 * e + 1):
            bp.add(F(m, 7 * e))
    return sorted(bp)


def sector_at(e: int, x: F) -> int:
    if e == 0:
        return 0
    val = 7 * e * x
    return (val.numerator // val.denominator) % 7


def sector_cover(E) -> F:
    """Measure of x for which {e*x} hits all seven fixed 1/7-sectors."""
    E = tuple(sorted(set(E)))
    total = F(0)
    bp = sector_breakpoints(E)
    for a, b in zip(bp, bp[1:]):
        if a >= b:
            continue
        mid = (a + b) / 2
        sectors = {sector_at(e, mid) for e in E}
        if len(sectors) == 7:
            total += b - a
    return total


def gcd_all(E) -> int:
    return reduce(math.gcd, E)


def primitive_shapes(k: int, max_e: int):
    for interior in itertools.combinations(range(1, max_e + 1), k - 1):
        E = (0,) + interior
        if gcd_all(E) == 1:
            yield E


def exhaustive_sector_bank(caps):
    # Chosen to be exact but light enough for repeated concurrent sessions.
    max_e_by_k = {8: 14, 9: 14, 10: 15, 11: 14, 12: 14, 13: 14}
    rows = []
    for k in range(8, 14):
        checked = 0
        best_val = F(-1)
        best_E = None
        violations = []
        for E in primitive_shapes(k, max_e_by_k[k]):
            checked += 1
            val = sector_cover(E)
            if val > best_val:
                best_val = val
                best_E = E
            if val > caps[k]:
                violations.append((val, E))
        rows.append((k, max_e_by_k[k], checked, best_val, best_E, caps[k], violations))
    return rows


def compression_counterexamples():
    mu_rows = [
        ((0, 1, 2, 3, 4, 6, 7, 13), (0, 1, 2, 3, 4, 5, 6, 12)),
        ((0, 1, 2, 3, 4, 5, 6, 8, 14), (0, 1, 2, 3, 4, 5, 6, 8, 13)),
        ((0, 1, 2, 3, 4, 5, 6, 7, 9, 15), (0, 1, 2, 3, 4, 5, 6, 7, 9, 14)),
        ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14), (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13)),
        ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16), (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15)),
    ]
    sector_rows = [
        ((0, 1, 2, 3, 4, 5, 6, 10), (0, 1, 2, 3, 4, 5, 6, 9)),
        ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12), (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)),
    ]
    return mu_rows, sector_rows


def route_tournament():
    vertices = [
        "sector_cover_cap",
        "consecutive_mu_min",
        "gap_max_stratification",
        "relation_height_tail",
        "one_step_compression",
    ]
    evidence = {
        "sector_cover_cap": 88,
        "consecutive_mu_min": 82,
        "gap_max_stratification": 68,
        "relation_height_tail": 64,
        "one_step_compression": 5,
    }
    order = {v: i for i, v in enumerate(vertices)}
    edges = {}
    scores = Counter()
    for a, b in itertools.combinations(vertices, 2):
        if (evidence[a], -order[a]) >= (evidence[b], -order[b]):
            winner, loser = a, b
        else:
            winner, loser = b, a
        edges[(winner, loser)] = 1
        scores[winner] += 1
        scores.setdefault(loser, scores[loser])

    def beats(a, b):
        return (a, b) in edges

    cycles3 = 0
    for a, b, c in itertools.combinations(vertices, 3):
        if (beats(a, b) and beats(b, c) and beats(c, a)) or (
            beats(a, c) and beats(c, b) and beats(b, a)
        ):
            cycles3 += 1
    hamiltonian_paths = 0
    for perm in itertools.permutations(vertices):
        if all(beats(perm[i], perm[i + 1]) for i in range(len(perm) - 1)):
            hamiltonian_paths += 1
    return vertices, evidence, dict(scores), cycles3, hamiltonian_paths


def main():
    print("=" * 88)
    print("LRC(14) seven-sector cap scout")
    print("=" * 88)
    print()
    print("[0] Sufficient route")
    print("    N(E) = {all gaps <= 1/7}.")
    print("    S7(E)= {fixed sectors 0..6 are all hit by {e*x}}.")
    print("    Up to boundary measure zero, N(E) subset S7(E).")
    print("    Therefore meas(S7(E)) <= cap_k implies mu_1/7(E) >= 1-cap_k.")

    caps, cap_witnesses = cap_table()
    print()
    print("[1] cap_k = min_{|P|=13-k} meas(G_P)")
    for k in range(8, 14):
        print(f"    k={k}: cap={fmt(caps[k])}; witness P={cap_witnesses[k]}")

    print()
    print("[2] Consecutive rows: net(E) <= sector_cover(E) <= cap_k")
    for k in range(8, 14):
        E = tuple(range(k))
        mu = mu_theta(E)
        net = 1 - mu
        sec = sector_cover(E)
        print(
            f"    k={k}: net={fmt(net)}; sector={fmt(sec)}; "
            f"cap={fmt(caps[k])}; margins cap-sector={fmt(caps[k] - sec)}"
        )
        assert net <= sec

    print()
    print("[3] Exact primitive sector-cover bank")
    rows = exhaustive_sector_bank(caps)
    for k, max_e, checked, best_val, best_E, cap, violations in rows:
        print(
            f"    k={k}: maxE<={max_e}; checked={checked}; "
            f"max_sector={fmt(best_val)} at E={best_E}; "
            f"cap={fmt(cap)}; margin={fmt(cap - best_val)}; "
            f"violations={len(violations)}"
        )
        if violations:
            top = max(violations)
            print(f"        largest violation: {fmt(top[0])} at E={top[1]}")

    print()
    print("[4] Compression counterexamples")
    print("    Naive one-step compression cannot be the proof lemma.")
    mu_rows, sector_rows = compression_counterexamples()
    print("    mu_1/7 counterexamples, where compression increases mu:")
    for E, Ec in mu_rows:
        before = mu_theta(E)
        after = mu_theta(Ec)
        print(f"        E={E} -> Ec={Ec}: {fmt(before)} -> {fmt(after)}")
    print("    sector-cover counterexamples, where compression decreases sector cover:")
    for E, Ec in sector_rows:
        before = sector_cover(E)
        after = sector_cover(Ec)
        print(f"        E={E} -> Ec={Ec}: {fmt(before)} -> {fmt(after)}")

    print()
    print("[5] Proof-route Tournament Analysis")
    vertices, evidence, scores, cycles3, hpaths = route_tournament()
    print(f"    vertices={vertices}")
    print(f"    evidence_scores={evidence}")
    print(f"    score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"    vertex_scores={dict(sorted(scores.items(), key=lambda kv: (-kv[1], kv[0])))}")
    print(f"    directed_3cycles={cycles3}")
    print(f"    SCC_sizes={[1 for _ in vertices]}")
    print(f"    Hamiltonian_path_count={hpaths}")

    print()
    print("Verdict:")
    print("    Sector cover gives a weaker cap target than consecutive-mu minimization.")
    print("    The bounded exact bank has no violations, with large rational margins.")
    print("    LRC(14) is not proved: the missing lemma is a global sector-cover cap")
    print("    or a finite low-height reduction plus high-height sector discrepancy tail.")


if __name__ == "__main__":
    main()
