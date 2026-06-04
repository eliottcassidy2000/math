#!/usr/bin/env python3
"""S621 part 3 — crisp boolean verification of two claims about the tight family."""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, sys
sys.path.insert(0, '04-computation')
from lrc_tight_enum_s621 import enum_tight, norm, gap_exact   # lightweight only

def witnesses(speeds):
    V = [abs(v) for v in speeds]; g = gap_exact(speeds); cands = set()
    for i in range(len(V)):
        vi = V[i]
        for k in range(0, 2*vi+1):
            t = Fr(2*k+1, 2*vi)
            if 0 < t <= Fr(1, 2): cands.add(t)
        for j in range(i):
            vj = V[j]
            for d in (vi+vj, abs(vi-vj)):
                if d == 0: continue
                kk = 1
                while Fr(kk, d) <= Fr(1, 2):
                    cands.add(Fr(kk, d)); kk += 1
    return g, sorted(t for t in cands if min(norm(v*t) for v in V) == g)

def phi(m): return sum(1 for k in range(1, m) if gcd(k, m) == 1)

print("CLAIM D (Kravitz ladder: every primitive set with gap<1/n has g=s/(ns+1), s in Z+):")
for n in [3, 4, 5, 6]:
    R = 3*n+1; bad = []; tot = 0
    for s in itertools.combinations(range(1, R+1), n):
        if reduce(gcd, s) != 1: continue
        g = gap_exact(list(s))
        if g < Fr(1, n):
            tot += 1
            if (g/(1-n*g)).denominator != 1: bad.append((s, str(g)))
    print(f"  n={n}: {tot} sets with gap<1/n; OFF-ladder: {len(bad)} {bad[:3]}")

print("\nCLAIM B (witnesses of every tight set == (Z/m)^* orbit in (0,1/2], count phi(m)/2):")
for n in range(3, 8):
    m = n+1; orbit = sorted(Fr(k, m) for k in range(1, m) if gcd(k, m) == 1 and Fr(k, m) <= Fr(1, 2))
    allok = True
    for s in enum_tight(n, max(2*n+4, 14), K=120):
        _, W = witnesses(list(s))
        if W != orbit: allok = False; print("   MISMATCH", s, [str(x) for x in W])
    print(f"  n={n}: orbit size {len(orbit)} (phi(m)/2={phi(m)//2}); all tights match={allok}")
