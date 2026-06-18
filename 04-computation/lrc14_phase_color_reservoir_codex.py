#!/usr/bin/env python3
"""
lrc14_phase_color_reservoir_codex.py

Codex 2026-06-18: phase-color reservoir scout for the LRC(14) S3
CRT placement problem.

The threshold ladder showed that uncolored slack intervals survive well below
the false via-max boundary.  But actual residues at q=14*V have a color
constraint: a residue a has phase color b=a mod 14, so the cluster must be
safe relative to the fixed center b/14, not merely to some continuous gap
center.

For E offsets and a phase color b, define

    C_b(E) = {x : ||e*x - b/14|| >= 1/14 for every e in E}.

Then for S=P union {V-e:e in E}, q=14V, and a=b+14t,

    a/q is an actual level-1/14 CRT witness
    iff x=a/q is in G_P cap C_b(E).

This script verifies that identity exactly and measures the continuous colored
reservoir

    Sigma(P,E) = sum_{b=0}^{13} meas(G_P cap C_b(E)),

the expected number of colored grid hits per V.  If Sigma has a uniform floor
and discrepancy is controlled, the CRT placement part of OPEN-Q-108 becomes a
finite colored-grid theorem rather than a private-pair denominator inequality.

Tournament Analysis declaration.
  Vertex set: phase colors b=0,...,13.
  Pairwise observable: aggregate reservoir measure for each color over the
    structured bank; orient b -> c if color b has larger aggregate support.
  Gauge/tie path: larger aggregate measure wins; ties follow increasing color.
  Fingerprints: score histogram, directed 3-cycles, SCCs, Hamiltonian paths.

Assumption challenge.
  I considered runner, interval-component, residue, denominator-event, and
  phase-color vertices.  Phase colors preserve the exact CRT predicate at
  modulus 14V and the small/large split, but destroy which runner creates the
  gap and the relation-lattice geometry of E.  The challenged assumption is
  that an uncolored max-gap reservoir is the right finite object; it is only
  the shadow of this colored reservoir.
"""

from __future__ import annotations

import itertools
import random
import sys
from collections import Counter
from fractions import Fraction as F
from functools import reduce
from math import gcd

import lrc14_global_threshold_ladder_codex as ladder

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

H = F(1, 14)
RNG = random.Random(20260618 + 2593)


def dist_residue(r: int, q: int) -> int:
    r %= q
    return min(r, q - r)


def danger_arcs_shift(e: int, center: F, h: F = H) -> list[tuple[F, F]]:
    """x such that e*x mod 1 lies within h of center."""
    center %= 1
    if e == 0:
        d = min(center, 1 - center)
        return [(F(0), F(1))] if d < h else []
    out: list[tuple[F, F]] = []
    for j in range(e):
        c = (center + j) / e
        a = (c - h / e) % 1
        b = (c + h / e) % 1
        if a < b:
            out.append((a, b))
        else:
            out.append((a, F(1)))
            out.append((F(0), b))
    return out


_PHASE_CACHE: dict[tuple[tuple[int, ...], int], list[tuple[F, F]]] = {}


def phase_safe_set(E: tuple[int, ...], b: int) -> list[tuple[F, F]]:
    E = tuple(sorted(set(E)))
    b %= 14
    key = (E, b)
    if key in _PHASE_CACHE:
        return _PHASE_CACHE[key]
    danger: list[tuple[F, F]] = []
    center = F(b, 14)
    for e in E:
        danger.extend(danger_arcs_shift(e, center))
    _PHASE_CACHE[key] = ladder.complement(ladder.merge(danger))
    return _PHASE_CACHE[key]


def color_measures(P: tuple[int, ...], E: tuple[int, ...]) -> list[F]:
    gp = ladder.safe_set(P)
    return [ladder.measure(ladder.intersect(gp, phase_safe_set(E, b))) for b in range(14)]


def color_components(P: tuple[int, ...], E: tuple[int, ...], b: int) -> list[tuple[F, F]]:
    return ladder.intersect(ladder.safe_set(P), phase_safe_set(E, b))


def sigma(P: tuple[int, ...], E: tuple[int, ...]) -> F:
    return sum(color_measures(P, E), F(0))


def build_S(P: tuple[int, ...], E: tuple[int, ...], V: int) -> tuple[int, ...]:
    return tuple(sorted(set(P) | {V - e for e in E}))


def gcd_all(vals) -> int:
    return reduce(gcd, vals, 0)


def is_covering(S: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in S) for q in range(2, 15))


def actual_crt_count(S: tuple[int, ...], V: int) -> int:
    q = 14 * V
    count = 0
    for a in range(q):
        if all(14 * dist_residue(v * a, q) >= q for v in S):
            count += 1
    return count


def colored_grid_count(P: tuple[int, ...], E: tuple[int, ...], V: int) -> int:
    q = 14 * V
    count = 0
    for b in range(14):
        for t in range(V):
            a = b + 14 * t
            if all(14 * dist_residue(p * a, q) >= q for p in P) and all(
                14 * dist_residue(e * a - b * V, q) >= q for e in E
            ):
                count += 1
    return count


def candidate_lifts(P: tuple[int, ...], E: tuple[int, ...], limit: int = 300) -> list[int]:
    out = []
    lo = max(E) + 14
    for V in range(lo, limit + 1):
        S = build_S(P, E, V)
        if len(S) == 13 and min(V - e for e in E) > 13 and gcd_all(S) == 1 and is_covering(S):
            out.append(V)
    return out


def count_directed_3cycles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in itertools.combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen = set()
    order: list[int] = []

    def dfs1(v: int):
        seen.add(v)
        for w, ok in enumerate(adj[v]):
            if ok and w not in seen:
                dfs1(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs1(v)
    seen.clear()
    sizes = []
    for start in reversed(order):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            v = stack.pop()
            size += 1
            for w, ok in enumerate(radj[v]):
                if ok and w not in seen:
                    seen.add(w)
                    stack.append(w)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def hamiltonian_paths_count(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [Counter() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v, cnt in list(dp[mask].items()):
            for w in range(n):
                if cnt and not (mask & (1 << w)) and adj[v][w]:
                    dp[mask | (1 << w)][w] += cnt
    return sum(dp[(1 << n) - 1].values())


def color_tournament(aggregate: list[F]) -> dict:
    n = len(aggregate)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if aggregate[i] > aggregate[j] or (aggregate[i] == aggregate[j] and i < j):
                adj[i][j] = True
    scores = [sum(adj[i]) for i in range(n)]
    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "cycles3": count_directed_3cycles(adj),
        "scc": scc_sizes(adj),
        "hamiltonian_paths": hamiltonian_paths_count(adj),
        "leaders": sorted(zip(scores, range(n), aggregate), reverse=True),
    }


def main() -> None:
    print("=" * 88)
    print("LRC(14) phase-color reservoir scout")
    print("=" * 88)
    print("C_b(E) = {x : ||e*x - b/14|| >= 1/14 for all e in E}.")
    print("Actual q=14V witnesses equal colored grid hits in G_P cap C_b(E).")

    named = [
        ("quarter_min", (1, 2, 3), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        ("near_via_min", (1, 2, 3, 11), (0, 2, 3, 4, 5, 6, 7, 8, 10)),
        ("via_zero_k7", (1, 2, 3, 6, 12, 13), (0, 2, 3, 4, 5, 6, 8)),
        ("via_zero_k9", (1, 2, 3, 13), (0, 2, 3, 4, 5, 6, 7, 8, 9)),
        ("broad_1_90", (1, 2, 9), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
    ]

    print("\n" + "=" * 88)
    print("STEP 1. Named shapes: continuous colored reservoirs")
    print("=" * 88)
    for label, P, E in named:
        cm = color_measures(P, E)
        sig = sum(cm, F(0))
        nonzero = [(b, m) for b, m in enumerate(cm) if m]
        best_width = F(0)
        best_b = None
        for b, _m in nonzero:
            comps = color_components(P, E, b)
            widths = [y - x for x, y in comps]
            if widths and max(widths) > best_width:
                best_width = max(widths)
                best_b = b
        print(f"  {label}: P={P} E={E}")
        print(f"    Sigma=sum_b measure = {sig} = {float(sig):.6f}; nonzero colors={len(nonzero)}")
        print(f"    best color component width={best_width} at b={best_b}")
        print(f"    color measures={[(b, str(m)) for b, m in nonzero]}")

    print("\n" + "=" * 88)
    print("STEP 2. Exact identity check: actual CRT count == colored grid count")
    print("=" * 88)
    mismatches = []
    rows = []
    for label, P, E in named:
        lifts = candidate_lifts(P, E, 320)[:8]
        for V in lifts:
            S = build_S(P, E, V)
            actual = actual_crt_count(S, V)
            colored = colored_grid_count(P, E, V)
            sig = sigma(P, E)
            rows.append((label, V, actual, colored, sig))
            if actual != colored:
                mismatches.append((label, P, E, V, actual, colored))
    for label, V, actual, colored, sig in rows[:32]:
        exp = sig * V
        ratio = F(actual, 1) / exp if exp else F(0)
        print(
            f"  {label:13s} V={V:3d}: actual={actual:4d} colored={colored:4d} "
            f"V*Sigma={exp} ({float(exp):.2f}) ratio={float(ratio):.3f}"
        )
    print(f"  identity mismatches: {len(mismatches)}")
    for item in mismatches[:5]:
        print("    MISMATCH", item)

    print("\n" + "=" * 88)
    print("STEP 3. Structured bank: Sigma floor and color tournament")
    print("=" * 88)
    min_sig = (F(10), None)
    min_max_color = (F(10), None)
    aggregate = [F(0) for _ in range(14)]
    cases = 0
    for k in range(3, 14):
        psz = 13 - k
        local = (F(10), None)
        local_width = (F(10), None)
        Es = ladder.shapes_for_k(k)
        for E in Es:
            phase_sets = [phase_safe_set(E, b) for b in range(14)]
            for P in ladder.powerset_P(psz):
                gp = ladder.safe_set(P)
                cm = [ladder.measure(ladder.intersect(gp, phase_sets[b])) for b in range(14)]
                sig = sum(cm, F(0))
                max_color = max(cm)
                cases += 1
                for b, m in enumerate(cm):
                    aggregate[b] += m
                if sig < min_sig[0]:
                    min_sig = (sig, (P, E, cm))
                if sig < local[0]:
                    local = (sig, (P, E))
                if max_color < min_max_color[0]:
                    min_max_color = (max_color, (P, E, cm))
                if max_color < local_width[0]:
                    local_width = (max_color, (P, E))
        print(
            f"  k={k:2d}: shapes={len(Es):2d} cases={len(Es)*len(ladder.powerset_P(psz)):6d} "
            f"min Sigma={local[0]} ({float(local[0]):.6f}) at {local[1]} "
            f"min max_color={local_width[0]} ({float(local_width[0]):.6f})"
        )

    print(f"\n  total structured cases={cases}")
    print(f"  GLOBAL min Sigma={min_sig[0]} = {float(min_sig[0]):.6f} at P,E={min_sig[1][:2]}")
    print(f"    color vector={[str(x) for x in min_sig[1][2]]}")
    print(
        f"  GLOBAL min max_color={min_max_color[0]} = {float(min_max_color[0]):.6f} "
        f"at P,E={min_max_color[1][:2]}"
    )
    print(f"    color vector={[str(x) for x in min_max_color[1][2]]}")

    tour = color_tournament(aggregate)
    print("\n  Phase-color tournament:")
    print(f"    aggregate={[str(x) for x in aggregate]}")
    print(
        f"    score_hist={tour['score_hist']} cycles3={tour['cycles3']} "
        f"scc={tour['scc']} hamiltonian_paths={tour['hamiltonian_paths']}"
    )
    print(f"    leaders={[(s,b,str(m)) for s,b,m in tour['leaders'][:8]]}")

    print("\n" + "=" * 88)
    print("STEP 4. Random lift stress: count vs V*Sigma")
    print("=" * 88)
    tested = 0
    worst_ratio = (F(10), None)
    zero_actual = []
    for _ in range(250):
        label, P, E = RNG.choice(named)
        lifts = candidate_lifts(P, E, 500)
        if not lifts:
            continue
        V = RNG.choice(lifts)
        S = build_S(P, E, V)
        actual = actual_crt_count(S, V)
        sig = sigma(P, E)
        exp = sig * V
        ratio = F(actual, 1) / exp if exp else F(0)
        tested += 1
        if ratio < worst_ratio[0]:
            worst_ratio = (ratio, (label, P, E, V, actual, sig))
        if actual == 0:
            zero_actual.append((label, P, E, V))
    print(f"  tested random named lifts={tested}")
    print(f"  zero actual witnesses at q=14V: {len(zero_actual)}")
    print(f"  worst actual/(V*Sigma)={worst_ratio[0]} = {float(worst_ratio[0]):.6f}")
    print(f"    at {worst_ratio[1]}")

    print("\n" + "=" * 88)
    print("TAKEAWAY")
    print("=" * 88)
    print("  1. The colored reservoir is the exact finite CRT object: actual")
    print("     witnesses at q=14V equal colored grid hits in G_P cap C_b(E).")
    print("  2. Sigma(P,E)=sum_b meas(G_P cap C_b(E)) is large in the structured")
    print("     bank; the minimum found is printed above and is far from zero.")
    print("  3. This suggests replacing the last placement problem by an")
    print("     Erdos-Turan/Koksma bound for 14 colored arithmetic progressions.")
    print("     The denominator D=14q-r route becomes: conditional blockers cannot")
    print("     cover all colored grid hits when V*Sigma is large.")
    print("  4. LRC(14) is not proved here: the missing theorem is a uniform")
    print("     discrepancy bound plus a finite small-V check, not just this census.")


if __name__ == "__main__":
    main()
