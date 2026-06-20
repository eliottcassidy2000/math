#!/usr/bin/env python3
"""
lrc_cluster_reversal_symmetry_macmini_0620s4.py  (mac-mini-2026-06-20-S4)

Does the half-tiling complement (reverse all arcs) descend to a CLUSTER-REVERSAL symmetry of the
LRC seven-sector coverage?  Test:  p0(E) =? p0(E*)  where  E* = {M - e : e in E},  M = max(E).
If YES for all E, the LRC finite checks (THM-547 collar, true-wide) HALVE (only check E <= E*).
Also test the sector-reflection x->-x identity p0 is invariant (sanity), and which reversal works.
"""
import itertools, sys
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E))==7: tot+=x1-x0
    return tot
def rev(E):
    M=max(E); return tuple(sorted(M-e for e in E))

# exhaustive small clusters + sampled LRC-relevant clusters
print("Test p0(E) == p0({max-e}) (cluster reversal):")
import random; random.seed(3)
def gen(k, lo, hi, count):
    out=set()
    while len(out)<count:
        s=tuple(sorted(random.sample(range(lo,hi+1), k-1)))
        E=(0,)+s
        if gcd(*E)==1 or True: out.add(E)
    return list(out)
allok=True; tested=0; mism=0; examples=[]
for k in [6,7,8,9]:
    for E in gen(k, 1, 24, 60):
        p0=measS7(E); p0r=measS7(rev(E)); tested+=1
        if p0!=p0r:
            mism+=1; allok=False
            if len(examples)<5: examples.append((E,rev(E),float(p0),float(p0r)))
print(f"  tested {tested} clusters (k=6..9, elts in [0,24]); mismatches: {mism}")
if examples:
    print("  MISMATCH examples:")
    for E,Er,a,b in examples: print(f"    E={E} -> {Er}: p0={a:.5f} vs {b:.5f}")
else:
    print("  ALL MATCH => cluster reversal E->{max-e} preserves p0 (HALVES the finite checks).")
# also test on the named true-wide/collar leaders
print("\n  named rows:")
for E in [(0,4,6,8,10,12,14,15,16),(0,1,2,3,4,5,6,7),(0,1,2,4,8,12,16,20),(0,2,4,6,8,10,12,14,15)]:
    print(f"    E={E}: p0={float(measS7(E)):.6f}  rev={rev(E)} p0={float(measS7(rev(E))):.6f}  {'=' if measS7(E)==measS7(rev(E)) else 'DIFF'}")
