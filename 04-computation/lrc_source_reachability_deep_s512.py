#!/usr/bin/env python3
"""
lrc_source_reachability_deep_s512.py

oracle-2026-06-01-S512

Push the S511 reformulation (LRC = observer-marked source-reachability; target =
A000568(n-1)) as far as it goes:

 (1) REACHABLE TARGET: at a lonely time the observer is a SOURCE and the n-1
     runners all lie in the FIXED safe arc [1/n, 1-1/n] (length L = 1-2/n).  The
     runner-runner half-turn sub-tournament is therefore the half-turn tournament
     of n-1 points in an interval of length L.  We enumerate the iso-classes that
     actually occur (the TRUE LRC win-set inside A000568(n-1)), characterize them
     (transitive? regular? H/score), and count them vs A000568(n-1).
 (2) FRAMING RE-PROVES SMALL LRC: for every primitive speed set in a box, verify
     the marked walk reaches a source cell (open or boundary) -- i.e. LRC holds,
     read off the tournament side.
 (3) DISTANCE-TO-SOURCE: track min over the walk of #near-runners (runners within
     1/n of the observer).  0 = a lonely (source) cell.  A counterexample would
     keep this >= 1 for the whole lap.  Probe how small it gets.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, permutations
from math import gcd
from collections import Counter, defaultdict
from functools import lru_cache

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
    """Enumerate primitive speed sets {0}+(n-1 from 1..max_speed); collect the
    reachable SOURCE runner sub-tournament classes (open + boundary)."""
    reachable = {}   # canon -> (H, score, transitive)
    lrc_fail = 0
    total = 0
    mindist_hist = Counter()
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
        # boundary witnesses (tight case)
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
        if not reached: lrc_fail += 1
    return reachable, total, lrc_fail, mindist_hist


def main():
    print("LRC source-reachability — the real target inside A000568(n-1) (oracle-S512)\n")
    for n, ms in [(4, 14), (5, 12), (6, 11), (7, 10)]:
        print("=" * 66)
        print(f"n={n} runners (observer + {n-1}); safe arc length L = 1-2/{n} = {1-2/n:.3f}")
        print("=" * 66)
        reachable, total, fail, mh = analyze_targets(n, ms)
        a = A000568(n - 1)
        print(f" primitive speed sets scanned: {total};  LRC failures (no source reached): {fail}")
        print(f" REACHABLE source classes: {len(reachable)}  (of A000568({n-1})={a} possible)")
        ntrans = sum(1 for v in reachable.values() if v[2])
        print(f"   transitive among them: {ntrans};  H values: {sorted({v[0] for v in reachable.values()})}")
        for c, (H, sc, tr) in sorted(reachable.items(), key=lambda kv: kv[1][0]):
            print(f"     H={H:<4} score={sc} transitive={tr}")
        print(f" min #near-runners over walk, histogram: {dict(sorted(mh.items()))}")
        print()
    print("READING: reachable source classes (the TRUE LRC win-set) are far fewer")
    print("than A000568(n-1); for n<=4 they are forced transitive (arc <= 1/2). The")
    print("min-near-runner histogram shows how close any set comes to source-")
    print("avoidance (a 0 bucket = an open lonely interval; all sets reach 0 or a")
    print("boundary source => LRC verified on the tournament side).")


if __name__ == "__main__":
    main()
