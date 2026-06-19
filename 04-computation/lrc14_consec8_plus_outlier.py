#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) k=9: the family  E_W = {0,1,2,3,4,5,6,7} U {W}  (consec_8 + one outlier).
The brute force found these are the worst non-AP 9-sets. Scan W over a large range
to find the SUPREMUM L_y of this family, and confirm it stays < cap_9 = 0.49426.

Also scan consec_8 + outlier where the outlier is NEGATIVE (i.e. consec block not at
the bottom): E = {0,1,...,7} translated, with outlier below -> same as {0} U {b,b+1,...,b+7}.
By translation-invariance of the ORBIT mod 1?  NO: scaling x is the only invariance.
Translating E by a constant c multiplies nothing -- frac((e+c)x) shifts all points by
frac(cx), a RIGID ROTATION of the orbit on the circle.  N counts MISSED sectors; a rigid
rotation of all orbit points does NOT preserve which sectors are missed (sectors fixed).
So translation matters.  We therefore also test the outlier on the LEFT.

kind-pasteur-2026-06-19.
"""
import sys
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def N_at(E, x):
    hit = set()
    for e in E:
        v = e * x
        v = v - (v.numerator // v.denominator)
        hit.add((v.numerator * 7) // v.denominator)
    return sum(1 for j in range(1, 7) if j not in hit)

def dist_p(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        ae = abs(e)
        for a in range(0, 7*ae + 1):
            bps.add(Fraction(a, 7*ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        p[N_at(E, (lo+hi)/2)] += (hi-lo)
    return p

g = [Fraction(-(t-2)*(t-3)*(t-6), 36) for t in range(7)]
def L_y(E):
    p = dist_p(E)
    return sum(p[t]*g[t] for t in range(7))

cap = Fraction(49426, 100000)
print(f"cap_9 = {float(cap):.5f}\n")
print("Family A:  E = {0,1,2,3,4,5,6,7, W}, outlier on the RIGHT")
bestA = (Fraction(-1), None)
for W in range(9, 120):
    E = list(range(8)) + [W]
    Lv = L_y(E)
    if Lv > bestA[0]:
        bestA = (Lv, W)
print(f"  sup over W in [9,119]: L_y={float(bestA[0]):.6f} at W={bestA[1]}  margin={float(cap-bestA[0]):.6f}")

# print a few small W explicitly
for W in range(9, 20):
    print(f"    W={W}: L_y={float(L_y(list(range(8))+[W])):.6f}")

print("\nFamily B:  E = {0} U {b, b+1, ..., b+7}, consec block shifted up (outlier=0 on left)")
bestB = (Fraction(-1), None)
for b in range(2, 120):
    E = [0] + list(range(b, b+8))
    Lv = L_y(E)
    if Lv > bestB[0]:
        bestB = (Lv, b)
print(f"  sup over b in [2,119]: L_y={float(bestB[0]):.6f} at b={bestB[1]}  margin={float(cap-bestB[0]):.6f}")
for b in range(2, 12):
    print(f"    b={b}: L_y={float(L_y([0]+list(range(b,b+8)))):.6f}")
