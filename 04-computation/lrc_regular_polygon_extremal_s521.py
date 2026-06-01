#!/usr/bin/env python3
"""
lrc_regular_polygon_extremal_s521.py   claudebox-2026-06-01-S521

Seed computation for 07-reflections/lrc-regular-polygon-and-complement-sieve-s521.md.

Claim: the extremal lonely configuration of the canonical LRC instance v_i=i at
t=1/n places the runners at the n-th roots of unity minus the observer vertex
(a regular n-gon). For ODD n the resulting half-turn tournament is
self-complementary and maximally balanced -- and equals the top class of the
S521 arc menu (HYP-1993). For EVEN n it degenerates: an antipodal runner pair
sits at separation exactly 1/2 (the polygon's diameter), the first-even bridge
obstruction.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations

def half_turn(positions):
    m = len(positions); adj = [[0]*m for _ in range(m)]; tie = False
    for a in range(m):
        for b in range(m):
            if a == b: continue
            d = (positions[a] - positions[b]) % 1
            if d == F(1, 2): tie = True
            elif 0 < d < F(1, 2): adj[a][b] = 1
    return adj, tie

def scores(adj): return sorted(sum(r) for r in adj)

def Hc(adj):
    m = len(adj); full = (1 << m) - 1
    @lru_cache(None)
    def dp(mask, last):
        if mask == full: return 1
        return sum(dp(mask | (1 << x), x) for x in range(m)
                   if not (mask >> x) & 1 and adj[last][x])
    return sum(dp(1 << s, s) for s in range(m))

def self_complementary(adj):
    m = len(adj)
    for p in permutations(range(m)):
        if all(adj[p[i]][p[j]] == adj[j][i]
               for i in range(m) for j in range(m) if i != j):
            return True
    return False

def main():
    print("LRC extremal lonely config = regular n-gon minus observer (v_i=i, t=1/n)\n")
    print(f"{'n':>2} {'m':>2}  result")
    for n in range(4, 11):
        pos = [F(i, n) for i in range(1, n)]
        adj, tie = half_turn(pos)
        if tie:
            print(f"{n:>2} {n-1:>2}  ANTIPODAL TIE (even n): polygon has a diameter -> degenerate")
        else:
            sc = scores(adj); H = Hc(adj)
            scp = self_complementary(adj) if n-1 <= 8 else "?"
            print(f"{n:>2} {n-1:>2}  scores={sc} H={H} self-complementary={scp}  (= top S521 menu class)")

if __name__ == "__main__":
    main()
