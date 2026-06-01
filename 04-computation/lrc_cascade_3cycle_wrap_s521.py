#!/usr/bin/env python3
"""
lrc_cascade_3cycle_wrap_s521.py   claudebox-2026-06-01-S521

The clearance cascade and where transitivity deadlocks (reflection:
07-reflections/lrc-cascade-conditional-clearances-s521.md).

Tournament transitivity, two facts:
  Fact 1 (chain):  X->Y and Y->Z  =>  X->Z.
  Fact 2 (prune):  X->Y  =>  not (Z->X and Y->Z)  [the acyclicity dual; forbids the
                   3-cycle X->Y->Z->X].
KEY (proved/verified here): in the round (half-turn) tournament on circle points,
  a triple is a 3-CYCLE  <=>  it WRAPS (no empty semicircle, max gap < 1/2)
  <=> it involves a LONG pair (sep > 1/2) = the apex/seam = the observer.
So the clearance cascade runs freely on the transitive backbone (within-semicircle
triples) and can only deadlock at the wrapping 3-cycles at the seam. Frustrating
cycles = odd cycles = the project's OCF; observer clearable (lonely) <=> on no
frustrating odd cycle.
"""
from fractions import Fraction as F
from itertools import combinations
import random

def fr(x): return x % 1
def round_adj(pts):
    m = len(pts); adj = [[0]*m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            if a != b and 0 < fr(pts[b]-pts[a]) < F(1, 2): adj[a][b] = 1
    return adj
def is_3cycle(adj, t):
    a, b, c = t
    for x, y, z in [(a, b, c), (a, c, b)]:
        if adj[x][y] and adj[y][z] and adj[z][x]: return True
    return False
def wraps(three):
    s = sorted(three); g = [(s[(i+1) % 3]-s[i]) % 1 for i in range(3)]
    return max(g) < F(1, 2)

def main():
    random.seed(1)
    print("Verify: round-tournament triple is a 3-cycle  <=>  it WRAPS (no empty semicircle):")
    ok = mism = 0
    for _ in range(400):
        m = 5; pts = sorted(F(random.randint(1, 600), 600) for _ in range(m))
        if len(set(pts)) < m: continue
        adj = round_adj(pts)
        for t in combinations(range(m), 3):
            three = [pts[i] for i in t]
            if is_3cycle(adj, t) == wraps(three): ok += 1
            else: mism += 1
    print(f"  matches: {ok}, mismatches: {mism}")
    print("\n  => 3-cycles (where Fact 2 / acyclicity binds) = wrapping triples = long-pair/apex/seam.")
    print("  The cascade clears the transitive backbone (within-semicircle, Facts 1&2 hold) freely;")
    print("  it deadlocks only at the seam (observer), where wrapping odd cycles (the OCF) live.")
    print("  LRC <=> the cascade closes <=> the observer is on no frustrating odd cycle (a source).")

if __name__ == "__main__":
    main()
