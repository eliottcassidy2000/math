#!/usr/bin/env python3
"""
LRC(14) empty-window moment lower certificates (codex-2026-06-18).

This is the complement-side analogue of the THM-534 seven-sector moment LP.
For the seven anchored open arcs

    W_j(E) = {x : (j/14, j/14 + 1/7) contains no frac(e*x), e in E},
    j = 0,...,6,

write R(x) = # {j : x in W_j(E)} and

    EWLB_A(E) = meas( union_j W_j(E) ) = P(R > 0).

If h(r) = sum_s z_s C(r,s) <= 1[r>0] for r=0,...,7, then

    EWLB_A(E) >= Phi_z(E) := sum_s z_s T_s(E),
    T_s(E) = E[C(R,s)].

The script records exact rational certificates that already clear the LRC(14)
threshold rows on the consecutive/AP shape and tests the resulting fixed
functional over bounded primitive banks.  The remaining proof target is the
scalar extremal statement Phi_z(E) >= thr_k for all E.

Tournament Analysis:
  vertices: the seven anchored empty-window regions W_0,...,W_6, not runners.
  observable: certificate load assigned by equal Shapley splitting of h(R(x))
    among the active empty windows at x.
  switch/gauge: orient i -> j when load_i >= load_j; exact ties are oriented
    by the cyclic Hamiltonian tie path 0 -> 1 -> ... -> 6.
  quotient preserves: which fixed regions carry the lower-bound certificate.
  quotient destroys: runner labels, cyclic adjacency of actual phase gaps,
    and boundary timing.  This challenges the assumption that LRC tournament
    vertices should be runners rather than loop regions/proof obligations.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd
import random
import sys


sys.stdout.reconfigure(line_buffering=True)


A = [F(j, 14) for j in range(7)]
WIDTH = F(1, 7)

# HYP-2602 / EWLB thresholds 1 - cap_k for k=8..12.
THR = {
    8: F(3637, 5880),
    9: F(2025, 4004),
    10: F(36, 91),
    11: F(25, 91),
    12: F(1, 7),
}

# Universal polynomial minorants h_R(r) <= 1[r>0], in the binomial basis C(r,s).
CERTS = {
    1: [F(0), F(1, 7)],
    2: [F(0), F(2, 5), F(-1, 10)],
    3: [F(0), F(5, 7), F(-3, 7), F(1, 7)],
    4: [F(0), F(16, 21), F(-11, 21), F(2, 7), F(-2, 21)],
}

# First degree that clears the consecutive/AP row.
ROW_DEGREE = {8: 4, 9: 3, 10: 3, 11: 2, 12: 1}


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def merge_intervals(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    intervals = sorted(intervals)
    merged: list[list[F]] = []
    for lo, hi in intervals:
        if hi <= lo:
            continue
        if not merged or lo > merged[-1][1]:
            merged.append([lo, hi])
        elif hi > merged[-1][1]:
            merged[-1][1] = hi
    return [(lo, hi) for lo, hi in merged]


def empty_intervals(E: tuple[int, ...], a0: F) -> list[tuple[F, F]]:
    """Return exact intervals where the anchored open arc is empty."""
    hits: list[tuple[F, F]] = []
    for e in E:
        if e == 0:
            continue
        ee = F(e)
        # x in ((a0+m)/e, (a0+WIDTH+m)/e) intersect [0,1).
        for m in range(-2, e + 3):
            lo = max(F(0), (a0 + m) / ee)
            hi = min(F(1), (a0 + WIDTH + m) / ee)
            if lo < hi:
                hits.append((lo, hi))

    hits = merge_intervals(hits)
    empty: list[tuple[F, F]] = []
    prev = F(0)
    for lo, hi in hits:
        if lo > prev:
            empty.append((prev, lo))
        if hi > prev:
            prev = hi
    if prev < 1:
        empty.append((prev, F(1)))
    return empty


def event_profile(E: tuple[int, ...], z: list[F] | None = None):
    """Exact distribution of R, its binomial moments, union, and z-loads."""
    events = [empty_intervals(E, a0) for a0 in A]
    bps = {F(0), F(1)}
    for event in events:
        for lo, hi in event:
            bps.add(lo)
            bps.add(hi)
    cuts = sorted(bps)

    T = [F(0)] * 8
    dist = [F(0)] * 8
    loads = [F(0)] * 7
    for x0, x1 in zip(cuts, cuts[1:]):
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        length = x1 - x0
        active = []
        for j, event in enumerate(events):
            if any(lo < xm < hi for lo, hi in event):
                active.append(j)
        r = len(active)
        dist[r] += length
        for s in range(8):
            T[s] += length * comb(r, s)
        if z is not None and r:
            h = sum(z[s] * comb(r, s) for s in range(len(z)))
            for j in active:
                loads[j] += length * h / r

    union = sum(dist[1:])
    ie_union = sum(((-1) ** (s + 1)) * T[s] for s in range(1, 8))
    return T, dist, union, ie_union, loads


def phi(T: list[F], z: list[F]) -> F:
    return sum(z[s] * T[s] for s in range(len(z)))


def h_values(z: list[F]) -> list[F]:
    return [sum(z[s] * comb(r, s) for s in range(len(z))) for r in range(8)]


def primitive(E: tuple[int, ...]) -> bool:
    g = 0
    for e in E:
        g = gcd(g, e)
    return g == 1


def directed_3cycles(adj: list[list[bool]]) -> int:
    total = 0
    n = len(adj)
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)

    def reach(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w, ok in enumerate(graph[v]):
                if ok and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    unused = set(range(n))
    sizes = []
    while unused:
        v = min(unused)
        comp = reach(v, adj) & reach(v, rev)
        sizes.append(len(comp))
        unused -= comp
    return sorted(sizes, reverse=True)


def ham_path_count(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def load_tournament(loads: list[F]):
    n = len(loads)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if loads[i] > loads[j] or (loads[i] == loads[j] and i < j):
                adj[i][j] = True
            else:
                adj[j][i] = True
    scores = [sum(row) for row in adj]
    hist = {s: scores.count(s) for s in sorted(set(scores))}
    return {
        "scores": scores,
        "score_hist": hist,
        "directed_3cycles": directed_3cycles(adj),
        "scc_sizes": scc_sizes(adj),
        "hamiltonian_paths": ham_path_count(adj),
    }


def print_minorants() -> None:
    print("=== universal lower minorants h_R(r) <= 1[r>0] ===")
    for R, z in CERTS.items():
        vals = h_values(z)
        ok = vals[0] <= 0 and all(v <= 1 for v in vals[1:])
        print(f"R={R}: z={[str(v) for v in z]}")
        print(f"     h(0..7)={[str(v) for v in vals]} valid={ok}")


def print_ap_rows() -> None:
    print("\n=== AP/consecutive exact rows ===")
    for k in range(8, 13):
        E = tuple(range(k))
        R = ROW_DEGREE[k]
        z = CERTS[R]
        T, dist, union, ie_union, loads = event_profile(E, z)
        bound = phi(T, z)
        print(
            f"k={k} R={R}: Phi={fmt(bound)} thr={fmt(THR[k])} "
            f"margin={fmt(bound - THR[k])} EWLB={fmt(union)} IEok={union == ie_union}"
        )
        print(f"     R-dist={[str(d) for d in dist]}")
        print(f"     load-tournament={load_tournament(loads)}")


def bounded_sweeps() -> None:
    print("\n=== bounded primitive sweeps for fixed certificates ===")
    max_by_k = {8: 14, 9: 13, 10: 13, 11: 13, 12: 13}
    for k in range(8, 13):
        R = ROW_DEGREE[k]
        z = CERTS[R]
        best = None
        below = 0
        count = 0
        for rest in combinations(range(1, max_by_k[k] + 1), k - 1):
            E = (0,) + rest
            if not primitive(E):
                continue
            T, _dist, union, _ie, _loads = event_profile(E)
            value = phi(T, z)
            count += 1
            if value < THR[k]:
                below += 1
            if best is None or value < best[0]:
                best = (value, E, union)
        print(
            f"k={k} R={R} maxE<={max_by_k[k]} count={count}: "
            f"minPhi={fmt(best[0])} at {best[1]} EWLB={fmt(best[2])} below_thr={below}"
        )


def random_stress() -> None:
    print("\n=== small random exact stress ===")
    rng = random.Random(2608)
    for k in range(8, 13):
        R = ROW_DEGREE[k]
        z = CERTS[R]
        best = None
        for _ in range(8):
            pts = {0}
            while len(pts) < k:
                pts.add(rng.randint(1, 35))
            E = tuple(sorted(pts))
            if not primitive(E):
                continue
            T, _dist, union, _ie, _loads = event_profile(E)
            value = phi(T, z)
            if best is None or value < best[0]:
                best = (value, E, union)
        print(
            f"k={k} R={R}: min random Phi={fmt(best[0])} at {best[1]} "
            f"EWLB={fmt(best[2])} margin={fmt(best[0] - THR[k])}"
        )


def main() -> None:
    print("LRC(14) empty-window moment lower certificates")
    print("A={j/14: j=0..6}, window width=1/7")
    print_minorants()
    print_ap_rows()
    bounded_sweeps()
    random_stress()
    print("\nDONE: LRC(14) not proved; reduced route is Phi_z(E) >= thr_k.")


if __name__ == "__main__":
    main()
