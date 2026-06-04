#!/usr/bin/env python3
"""
S628 — practical speedups: (1) tight-instance ENUMERATION via the gap_shells filter + structural
pruning (the task we kept timing out on, S621); (2) UNIT-DISTANCE counting/search via the lattice
neighbor rule and CM-rotation restriction (HYP-2230).
"""
from fractions import Fraction as Fr
from math import gcd, comb
from functools import reduce
import itertools, time, math, cmath, sys
sys.path.insert(0, '04-computation')
from lrc_fast_s628 import gap_brute, gap_shells

# ========================= (1) tight-instance enumeration: brute vs fast-filter
# convention here: n 'movers' = the speeds; tight at gap 1/(n+1) (S620/S621 convention).
def is_tight_brute(S):
    n = len(S); return gap_brute(S) == Fr(1, n+1)

def is_tight_fast(S):
    """reject fast if a shell witness already beats 1/(n+1) (=> loose); confirm survivors by brute.
       gap_shells with Mmax=2n+1 (=2(#speeds)+1) covers the witness orbit for this convention."""
    n = len(S); target = Fr(1, n+1)
    gs = gap_shells(S, n+1, Mmax=2*n+1)        # n+1 'runners' -> shells up to 2(n+1)-1=2n+1
    if gs > target:                            # provably loose -> not tight, no brute needed
        return False
    return gap_brute(S) == target              # confirm (handles high-shell residual)

def enum_tight(n, R, fast):
    delta = Fr(1, n+1); out = []
    test = is_tight_fast if fast else is_tight_brute
    for s in itertools.combinations(range(1, R+1), n):
        if reduce(gcd, s) != 1: continue
        if fast and max(s) > 2*n-1:            # THM-411 pruning: v_max <= 2n-1
            # still must include sets within bound; skip only if ALL > bound impossible -> just prune this s
            continue
        if test(list(s)): out.append(s)
    return out

print("(1) TIGHT-INSTANCE ENUMERATION  brute vs fast (gap_shells filter + v_max<=2n-1 prune)")
for n, R in [(5, 14), (6, 16)]:
    t0 = time.perf_counter(); B = enum_tight(n, R, False); tb = time.perf_counter()-t0
    t0 = time.perf_counter(); F = enum_tight(n, R, True);  tf = time.perf_counter()-t0
    same = set(B) == set(F)
    print(f"   n={n} R={R}: brute {tb:.2f}s -> {len(B)} tight ;  fast {tf:.2f}s -> {len(F)} tight ;"
          f"  same={same}  speedup x{tb/tf:.1f}")

# ========================= (2) unit distance: lattice O(n) count + CM-rotation search
w = cmath.exp(1j*math.pi/3)
def tri_count_brute(pts):                       # O(k^2)
    e = 0
    for a, b in itertools.combinations(pts, 2):
        if abs(abs(a-b)-1) < 1e-7: e += 1
    return e
def tri_count_fast(ij):                          # O(k): triangular-lattice neighbor rule
    S = set(ij); e = 0
    for (i, j) in ij:
        for di, dj in ((1,0),(0,1),(1,-1)):      # 3 of 6 (avoid double count)
            if (i+di, j+dj) in S: e += 1
    return e

print("\n(2) UNIT-DISTANCE count on triangular lattice  brute O(k^2) vs neighbor O(k)")
rng = __import__('random').Random(1)
def disk(R): return [(i,j) for i in range(-R,R+1) for j in range(-R,R+1) if i*i+i*j+j*j<=R*R]
L = disk(7)
tb = tf = 0.0; mism = 0
for _ in range(300):
    sub = rng.sample(L, 22)
    pts = [i + j*w for (i,j) in sub]
    t0=time.perf_counter(); cb = tri_count_brute(pts); tb += time.perf_counter()-t0
    t0=time.perf_counter(); cf = tri_count_fast(sub);  tf += time.perf_counter()-t0
    if cb != cf: mism += 1
print(f"   mismatches: {mism}/300 ;  brute {tb*1000:.1f} ms  neighbor {tf*1000:.1f} ms  speedup x{tb/tf:.1f}")

# CM-rotation restriction: candidate unit-distance-creating rotations = alpha/conj(alpha), bounded height
we = cmath.exp(2j*math.pi/3)
def cm_rotations(H):
    rots = {}
    for a in range(-H, H+1):
        for b in range(-H, H+1):
            z = a + b*we
            if abs(z) < 1e-9: continue
            ang = round(math.degrees(cmath.phase(z/z.conjugate())) % 360, 4)
            rots[ang] = (a, b)
    return sorted(rots)
print("\n   CM-rotation search space: angles alpha/abar of Z[w], |height|<=4:",
      f"{len(cm_rotations(4))} distinct rotations (vs a continuum) — search these only (HYP-2230).")
