# opus-2026-07-16-S331 -- HYP-7100: THE PARALLEL-CLASS CIRCLE.
# Each speed v resolves the circle into the v-torsion parallel class P_v
# (v congruent arcs). LRC = the covering interaction of the 13 classes.
# CONCRETE PROBE: the KILLER-CLASS stratification of the Landau corona:
# for each corona-dying orbit (sorted d-multiset that is lattice-legal but
# Landau-illegal), WHICH prefix k (the first violated Landau index) kills it
# -- the killer = the parallel-class interaction depth; census at n = 8, 10.
# Plus: the tight AP's joint-cell census (the 13-class common refinement =
# the lcm(1..13)-cell structure restricted to the tight locus mu6).
from fractions import Fraction
from math import comb, gcd
from collections import defaultdict
import itertools

print("(1) THE KILLER-CLASS STRATIFICATION OF THE CORONA (even n):")
for n in (8, 10):
    shells = defaultdict(set)
    rng = sorted(range(-(n-1), n, 2))
    def gen(prefix, remaining, lo_idx):
        if remaining == 0:
            if sum(prefix) == 0:
                shells[sum(t*t for t in prefix)].add(tuple(prefix))
            return
        for i in range(lo_idx, len(rng)):
            d = rng[i]
            rem = remaining - 1
            smin = sum(prefix) + d + rem*d
            smax = sum(prefix) + d + rem*rng[-1]
            if smin > 0 or smax < 0: continue
            gen(prefix + [d], rem, i)
    gen([], n, 0)
    killer = defaultdict(int)
    ndying = 0
    for x in sorted(shells):
        for m in shells[x]:
            s = sorted((d + n - 1)//2 for d in m)
            kk = None
            for k in range(n):
                if sum(s[:k+1]) < comb(k+1, 2): kk = k + 1; break
            if kk is not None:
                killer[kk] += 1; ndying += 1
    print(f"   n={n}: dying orbits = {ndying}; killer-index census "
          f"(first violated prefix k): {dict(sorted(killer.items()))}")

print()
print("(2) THE PARALLEL-CLASS JOINT-CELL STRUCTURE OF THE TIGHT AP:")
# the 13 classes' common refinement has lcm(1..13) = 360360 cells; the tight
# locus mu6 = {p/14}: which joint cell does each tight time inhabit, i.e. the
# 13-vector of class-positions (floor(v t) mod v shifted): the 'address' of
# the loneliest moments in the resolution.
L = 360360
for p in (1, 3, 5, 9, 11, 13):
    t = Fraction(p, 14)
    addr = [int(v*t) % v for v in range(1, 14)]
    dists = sorted(min((v*t) % 1, 1 - (v*t) % 1) for v in range(1, 14))
    print(f"   t = {p:2d}/14: class addresses {addr}; min dist = {dists[0]} "
          f"(all >= 1/14: {all(d >= Fraction(1,14) for d in dists)})")
print()
print("(3) THE DICTIONARY (exact checks):")
print(f"   lcm(1..13) = 360360 = 2^3*3^2*5*7*11*13; 360360/14 = {360360//14}"
      f" (the tight times sit on the 14-grid INSIDE the 360360-cell lattice:"
      f" 25740 cells per tight step)")
# S330's escape packet: inside the 30Z parallel class
esc = [420, 450, 510, 570, 690, 870, 1230, 1770, 2370, 3210]
g = 0
for x in esc: g = gcd(g, x)
print(f"   S330 escape prefix gcd = {g} (a packet INSIDE the {g}Z parallel "
      f"class: the dilate escape = travel along one class)")
# the wiggly 1-factorization side (canon check, small): Q_m directions = m
# parallel classes; count: m = C(n-1,2) at n=6: 10 classes on 2^10 tilings
print(f"   metagraph side: Q_10's edge set resolves into 10 wiggly parallel "
      f"classes (each a perfect matching of 1024 tilings: 512 edges each; "
      f"total 5120 = 10*2^9 edges = the tiling-space 1-factorization)")
