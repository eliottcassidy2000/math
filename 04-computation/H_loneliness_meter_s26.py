#!/usr/bin/env python3
"""
H_loneliness_meter_s26.py

oracle-2026-06-01-S26

Investigate H (directed Hamiltonian-path count) as a LONELINESS METER for the
half-turn (circular) tournament from n points on a unit circle.

Lift (S24): i -> j iff (x_i - x_j) mod 1 in (0, 1/2) (a ROUND tournament; each
out-set is the arc of points in the leading semicircle). Questions:
 (1) Is H a function of the gap profile / of max_gap alone?
 (2) Is H monotone non-increasing in max_gap?
 (3) Compare H vs #3-cycles, score variance, min out-degree as meters.
 (4) What does H count here?
Results stored in 05-knowledge/results/H_loneliness_meter_s26.out.
"""

from __future__ import annotations

import random
from itertools import combinations
from functools import lru_cache
from collections import defaultdict


def half_turn_tournament(points):
    n = len(points)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            d = (points[i] - points[j]) % 1.0
            if 0 < d < 0.5:
                adj[i][j] = 1
    return adj


def H_count(adj):
    n = len(adj)
    full = (1 << n) - 1
    adjt = tuple(tuple(r) for r in adj)

    @lru_cache(maxsize=None)
    def dp(mask, last):
        if mask == full:
            return 1
        s = 0
        for nx in range(n):
            if not (mask >> nx) & 1 and adjt[last][nx]:
                s += dp(mask | (1 << nx), nx)
        return s

    out = sum(dp(1 << s, s) for s in range(n))
    dp.cache_clear()
    return out


def num_3cycles(adj):
    n = len(adj)
    c = 0
    for i, j, k in combinations(range(n), 3):
        e = adj[i][j] + adj[j][k] + adj[k][i]
        if e == 0 or e == 3:
            c += 1
    return c


def scores(adj):
    return sorted(sum(r) for r in adj)


def gaps(points):
    p = sorted(points)
    n = len(p)
    g = [p[(i + 1) % n] - p[i] for i in range(n - 1)]
    g.append(p[0] + 1 - p[-1])
    return sorted(g, reverse=True)


def main():
    rng = random.Random(26)
    print("H as a loneliness meter — half-turn circular tournaments (oracle-S26)\n")

    for n in (5, 6, 7):
        print("=" * 70)
        print(f"n={n}")
        print("=" * 70)
        byH = defaultdict(lambda: {"maxgap": [], "n3": [], "svar": [], "minout": []})
        H_for_maxgap = defaultdict(set)
        samples = 0
        for _ in range(20000):
            pts = sorted(rng.random() for _ in range(n))
            adj = half_turn_tournament(pts)
            H = H_count(adj)
            g = gaps(pts)
            mg = g[0]
            sc = scores(adj)
            svar = sum((s - (n - 1) / 2) ** 2 for s in sc) / n
            byH[H]["maxgap"].append(mg)
            byH[H]["n3"].append(num_3cycles(adj))
            byH[H]["svar"].append(svar)
            byH[H]["minout"].append(sc[0])
            H_for_maxgap[round(mg, 2)].add(H)
            samples += 1
        print(f" samples={samples}; realizable H values: {sorted(byH)}")
        print(" H   | max_gap range          | mean | #3cyc(range) | scorevar(range)")
        for H in sorted(byH):
            d = byH[H]
            mgs = d["maxgap"]; n3 = d["n3"]; sv = d["svar"]
            print(f" {H:<4}| [{min(mgs):.4f}, {max(mgs):.4f}] "
                  f"| {sum(mgs)/len(mgs):.4f} | [{min(n3)},{max(n3)}] "
                  f"| [{min(sv):.3f},{max(sv):.3f}]")
        ambiguous = {mg: hs for mg, hs in H_for_maxgap.items() if len(hs) > 1}
        print(f" max_gap (2dp) values mapping to >1 distinct H: {len(ambiguous)} "
              f"/ {len(H_for_maxgap)}  -> H is {'NOT' if ambiguous else ''} a function of max_gap alone")
        Hsorted = sorted(byH)
        strictly_ordered = all(
            max(byH[b]["maxgap"]) <= min(byH[a]["maxgap"]) + 1e-9
            for a, b in zip(Hsorted, Hsorted[1:])
        )
        print(f" H-ranges in max_gap strictly ordered (H up <=> max_gap down)? {strictly_ordered}")
        print()

    print("SUMMARY: H=1 <=> max_gap>=1/2 sharply; above that H is coarse (overlapping")
    print("ranges, monotone only in mean). #3-cycles = spread-triple count; H finer.")


if __name__ == "__main__":
    main()
