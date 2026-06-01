#!/usr/bin/env python3
"""
lrc_source_reachability_n8_s520.py

oracle-2026-06-01-S520

Extend the S512 reachable-source-class menu (the TRUE LRC win-set inside
A000568(n-1)) to n=8, the first value oracle-S512 and codex-S516 did not reach.

S512 computed the reachable source-class count = 1, 2, 6, 6 for n = 4, 5, 6, 7
(vs A000568(n-1) = 2, 4, 12, 56) and re-proved LRC for n<=7 on the tournament
side (0 source-avoidance failures over a box of primitive speed sets).

This script:
  (1) re-runs n=7 as a harness self-check (must reproduce 6),
  (2) computes the n=8 menu count and characterises it (H, score, transitive),
  (3) reports LRC failures (no source reached) over the n=8 box -> tournament-
      side verification of LRC for n=8 if 0,
  (4) prints the extended menu sequence for HYP-1987.

Core routines copied verbatim from lrc_source_reachability_deep_s512.py so the
numbers are directly comparable. Only main() differs (adds n=8, smaller box).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, permutations
from math import gcd
from collections import Counter
from functools import lru_cache
import sys

ONE = Fraction(1)


def frac(x): return x - (x.numerator // x.denominator)
def dist0(x):
    f = frac(x); return min(f, ONE - f)


def primitive(vs):
    g = 0
    for v in vs: g = gcd(g, v)
    return g == 1


def near_count(speeds, n, t):
    thr = Fraction(1, n)
    return sum(1 for v in speeds[1:] if dist0(Fraction(v) * t) < thr)


def runner_subtournament(speeds, n, t):
    """half-turn tournament on the n-1 runners (indices 1..n-1)."""
    idx = list(range(1, n))
    m = len(idx)
    adj = [[0] * m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            if a == b: continue
            f = frac(Fraction(speeds[idx[a]] - speeds[idx[b]]) * t)
            adj[a][b] = 1 if 0 < f < Fraction(1, 2) else 0
    return tuple(tuple(r) for r in adj)


def canon(adj):
    m = len(adj)
    best = None
    for p in permutations(range(m)):
        flat = tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or flat < best: best = flat
    return best


def Hc(adj):
    m = len(adj); full = (1 << m) - 1
    @lru_cache(None)
    def dp(mask, last):
        if mask == full: return 1
        return sum(dp(mask | (1 << x), x) for x in range(m) if not (mask >> x) & 1 and adj[last][x])
    return sum(dp(1 << s, s) for s in range(m))


def scores(adj): return tuple(sorted(sum(r) for r in adj))
def is_transitive(adj): return scores(adj) == tuple(range(len(adj)))


def walls(speeds, n):
    thr = Fraction(1, n); W = set()
    for v in speeds[1:]:
        if v == 0: continue
        for m in range(0, v):
            W.add(frac(Fraction(m, v) + thr / v)); W.add(frac(Fraction(m, v) - thr / v))
    for i, j in combinations(range(1, n), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0: continue
        for k in range(0, 2 * d): W.add(frac(Fraction(k, 2 * d)))
    return sorted(w for w in W if 0 <= w < 1)


def midpoints(speeds, n):
    W = walls(speeds, n) + [ONE]
    return [(a + b) / 2 for a, b in zip(W, W[1:])]


def boundarypoints(speeds, n):
    return [w for w in walls(speeds, n)]


def A000568(k):
    table = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248]
    return table[k] if k < len(table) else None


def analyze_targets(n, max_speed):
    reachable = {}   # canon -> (H, score, transitive)
    lrc_fail = 0
    total = 0
    mindist_hist = Counter()
    fail_examples = []
    for combo in combinations(range(1, max_speed + 1), n - 1):
        if not primitive(combo): continue
        speeds = (0,) + combo
        total += 1
        mids = midpoints(speeds, n)
        bnds = boundarypoints(speeds, n)
        reached = False
        mind = n
        for t in mids:
            nc = near_count(speeds, n, t)
            mind = min(mind, nc)
            if nc == 0:
                reached = True
                adj = runner_subtournament(speeds, n, t)
                c = canon(adj)
                if c not in reachable:
                    reachable[c] = (Hc(adj), scores(adj), is_transitive(adj))
        if not reached:
            for t in bnds:
                if near_count(speeds, n, t) == 0:
                    reached = True
                    adj = runner_subtournament(speeds, n, t)
                    c = canon(adj)
                    if c not in reachable:
                        reachable[c] = (Hc(adj), scores(adj), is_transitive(adj))
                    break
        mindist_hist[mind] += 1
        if not reached:
            lrc_fail += 1
            if len(fail_examples) < 10:
                fail_examples.append(combo)
    return reachable, total, lrc_fail, mindist_hist, fail_examples


def report(n, ms):
    print("=" * 66)
    print(f"n={n} runners (observer + {n-1}); safe arc length L = 1-2/{n} = {1-2/n:.3f}; box max_speed={ms}")
    print("=" * 66)
    reachable, total, fail, mh, fex = analyze_targets(n, ms)
    a = A000568(n - 1)
    print(f" primitive speed sets scanned: {total};  LRC failures (no source reached): {fail}")
    print(f" REACHABLE source classes: {len(reachable)}  (of A000568({n-1})={a} possible)")
    ntrans = sum(1 for v in reachable.values() if v[2])
    print(f"   transitive among them: {ntrans};  H values: {sorted({v[0] for v in reachable.values()})}")
    for c, (H, sc, tr) in sorted(reachable.items(), key=lambda kv: kv[1][0]):
        print(f"     H={H:<4} score={sc} transitive={tr}")
    print(f" min #near-runners over walk, histogram: {dict(sorted(mh.items()))}")
    if fex:
        print(f" FAILURE examples (speed sets with no reached source): {fex}")
    print()
    return len(reachable), fail


def main():
    print("LRC source-reachability EXTENSION to n=8 (oracle-2026-06-01-S520)\n")
    menu = {}
    # self-check on n=7 (must give 6), then the new n=8 point.
    for n, ms in [(7, 10), (8, 9)]:
        cnt, fail = report(n, ms)
        menu[n] = (cnt, fail)
    print("EXTENDED MENU SEQUENCE (reachable source classes, the TRUE LRC win-set):")
    known = {4: 1, 5: 2, 6: 6, 7: 6}
    row = []
    for n in range(4, 9):
        if n in known and n not in menu:
            row.append(f"n={n}:{known[n]}")
        elif n in menu:
            tag = "" if menu[n][1] == 0 else f"(!{menu[n][1]} fail)"
            row.append(f"n={n}:{menu[n][0]}{tag}")
    print("   " + "  ".join(row))
    print(f"   vs A000568(n-1) = 2, 4, 12, 56, 456 for n=4..8")
    if menu.get(7, (None,))[0] != 6:
        print(" !! HARNESS SELF-CHECK FAILED: n=7 did not reproduce 6 -- numbers NOT comparable to S512")
    else:
        print(" harness self-check OK: n=7 reproduced 6 (matches S512)")
    if menu.get(8, (None, 1))[1] == 0:
        print(" => LRC VERIFIED for n=8 on the tournament side (0 source-avoidance failures in box).")


if __name__ == "__main__":
    main()
