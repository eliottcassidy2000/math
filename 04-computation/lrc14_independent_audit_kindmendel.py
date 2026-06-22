#!/usr/bin/env python3
"""Independent exact-rational re-derivation of the LRC(14) sector-route quantities,
from definitions, to verify the published values and the positivity linchpin.
kind-mendel-2026-06-22-S1."""
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
from itertools import combinations

def lcm(a,b): return a*b//gcd(a,b)
def norm(y):
    "||y|| = distance from rational y to nearest integer"
    f = y - (y.numerator//y.denominator)   # frac in [0,1)
    return min(f, 1-f)

def meas_GP(P):
    "meas{x in [0,1): ||p x|| >= 1/14 for all p in P}, exact."
    if not P: return F(1)
    bps=set([F(0),F(1)])
    for p in P:
        for n in range(p):
            for s in (F(-1,14*p),F(1,14*p)):
                x=(F(n,p)+s)%1
                bps.add(x)
    bps=sorted(bps)
    total=F(0)
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2
        if all(norm(p*mid)>=F(1,14) for p in P):
            total+=b-a
    return total

def p0(E):
    "meas{x: all 7 sectors [j/7,(j+1)/7) hit by {frac(e x): e in E}}, exact."
    pos=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in pos:
        for a in range(7*e+1):
            bps.add(F(a,7*e)%1 if F(a,7*e)<1 else F(0))
            bps.add(F(a,7*e))
    bps=sorted(x for x in bps if 0<=x<=1)
    total=F(0)
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid=(a+b)/2
        sectors=set()
        for e in E:
            f=(e*mid)%1
            sectors.add(int(7*f))   # floor since f<1 exact
        if sectors>=set(range(7)):
            total+=b-a
    return total

# --- CHECK 1: reproduce published per-size minima of meas(G_P), P subset {1..13} ---
published={1:F(6,7),2:F(66,91),3:F(55,91),4:F(1979,4004),5:F(2243,5880),
           6:F(3029,10780)}  # psz=1..6 (k=12..7); higher psz expensive, check 1..6
print("=== CHECK 1: min meas(G_P) over |P|=psz, P subset {1..13} ===")
allok=True
for psz in range(1,7):
    best=None
    for P in combinations(range(1,14),psz):
        m=meas_GP(P)
        if best is None or m<best[0]: best=(m,P)
    ok = best[0]==published[psz]
    allok&=ok
    print(f"psz={psz}: min={best[0]}={float(best[0]):.6f} at P={best[1]}  published={published[psz]}  {'OK' if ok else 'MISMATCH'}")
print("CHECK 1", "PASS" if allok else "FAIL")

# --- CHECK 2: p0(consec_8) should be 481/1470 ---
print("\n=== CHECK 2: p0(consec_k) ===")
v=p0(list(range(8)))
print(f"p0(consec_8)={v}={float(v):.6f}  published=481/1470={float(F(481,1470)):.6f}  {'OK' if v==F(481,1470) else 'MISMATCH'}")
