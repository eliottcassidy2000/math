#!/usr/bin/env python3
"""
The simultaneous-coverage criterion (THM-1016)  — mac-mini-2026-07-18-S112
==========================================================================
Owner directive: attack the metric wall with tooth positions.

For a tight n-set with sheet number s, on-sheet part E=sU and off-sheet F:
  T_{w,j} = { tau : ||w(tau+j)/s|| <= L },  L = 1/(n+1)
  Cov(F,s,L) = intersect_{j in Z/s} union_{w in F} T_{w,j}
THEOREM (proved): G = {tau : phi_U(tau) > L} satisfies G subset Cov.  G is nonempty
(M(U) >= 1/(|U|+1) > L), so Cov = EMPTY excludes the configuration outright.

This CROSSES the counting wall of THM-1006 sec.H (capacity+primitivity satisfiable
for every val): here we measure how many capacity+primitivity-passing configs the
metric kills.
"""
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations

G = 300000
tau = np.arange(G) / G
L = 1.0 / 13

def cov_empty(F, s, L=L):
    inter = np.ones(G, dtype=bool)
    for j in range(s):
        un = np.zeros(G, dtype=bool)
        for w in F:
            x = (w * (tau + j) / s) % 1.0
            un |= (np.minimum(x, 1 - x) <= L + 1e-12)
        inter &= un
        if not inter.any(): return True
    return False

def maxint(F, s, L=L):
    inter = np.ones(G, dtype=bool)
    for j in range(s):
        un = np.zeros(G, dtype=bool)
        for w in F:
            x = (w * (tau + j) / s) % 1.0
            un |= (np.minimum(x, 1 - x) <= L + 1e-12)
        inter &= un
    if not inter.any(): return 0.0
    if inter.all(): return 1.0
    d = np.diff(np.concatenate(([0], inter.view(np.int8), [0])))
    st = np.where(d == 1)[0]; en = np.where(d == -1)[0]; runs = en - st
    best = runs.max()
    if inter[0] and inter[-1]: best = max(best, runs[0] + runs[-1])
    return best / G

def capacity_ok(F, s):
    tot = 0.0
    for w in F:
        g = gcd(w, s); D = s // g
        if D < 2: return False
        tot += (int(2 * D / 13) + 1) / D
    return tot >= 1.0

def prim_ok(F, s):
    return reduce(gcd, [gcd(w, s) for w in F]) == 1

print("=" * 72)
print("(1) smallest strict gain: s=2, F={1,3} — capacity OK, primitivity OK, Cov EMPTY")
print("=" * 72)
print(f"  capacity_ok={capacity_ok([1,3],2)}  prim_ok={prim_ok([1,3],2)}  cov_empty={cov_empty([1,3],2)}")
print("  closed form: (T_10 u T_30) = [0,2/13] u [24/39,28/39];")
print("               (T_11 u T_31) = [11/39,15/39] u [11/13,1)   -> intersection EMPTY")

print()
print("=" * 72)
print("(2) THE STRICT GAIN over capacity+primitivity (THM-1006 sec.H's irreducible class)")
print("=" * 72)
passed = 0; gain = []
for s in range(2, 9):
    pool = [w for w in range(1, 25) if w % s != 0]
    for r in (2, 3):
        for F in combinations(pool, r):
            F = list(F)
            if not capacity_ok(F, s) or not prim_ok(F, s): continue
            passed += 1
            if cov_empty(F, s): gain.append((s, F))
from collections import Counter
print(f"  pass capacity AND primitivity: {passed}")
print(f"  of those EXCLUDED by tooth positions (Cov empty): {len(gain)} ({100*len(gain)/passed:.1f}%)")
print(f"  by sheet number s: {dict(sorted(Counter(s for s,_ in gain).items()))}")

print()
print("=" * 72)
print("(3) what survives: maxint(Cov) decays like 1/w  => surviving configs need LARGE")
print("    off-sheet speeds and a correspondingly tiny on-sheet safe interval")
print("=" * 72)
for F in ([1,7],[1,11],[1,19],[1,31],[3,5]):
    print(f"  s=2 F={F}: maxint(Cov)={maxint(F,2):.5f}")
print()
print("VERDICT: metric kills SMALL off-sheet configs (55.8% of the counting-irreducible")
print("  class); LARGE off-sheet configs survive and need a height bound. Complementary.")
