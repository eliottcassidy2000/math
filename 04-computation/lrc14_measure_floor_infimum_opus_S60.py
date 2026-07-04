#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
IS THERE AN ACTUAL (UNIFORM POSITIVE) MEASURE FLOOR?  inf meas(lonely S) over primitive covering
families: bounded below, or -> 0 as magnitude grows?  (opus-2026-07-03-S60)

meas(lonely S) > 0 for all primitive covering S IS LRC(14). The question of a UNIFORM floor is:
inf meas > 0? Two heuristics disagree:
 * SLOPE bound (proved here): meas(lonely S) >= 2(M(S) - 1/14)/v_max, since F(t)=min_v||vt|| is
   v_max-Lipschitz and peaks at M(S). This -> 0 as v_max -> infinity (if M stays ~ tight).
 * DECORRELATION (THM-611) + RIGIDITY (kps HYP-4062, tight locus = AP, covering repelled by the
   covering-min gap 14/183): a single large runner decorrelates (meas -> (6/7)meas(rest), bounded),
   and covering can't stay near the tight AP -> suggests meas bounded below.
This settles it by DIRECT SEARCH: min meas over primitive covering families at growing magnitude caps.
If min meas -> 0, there is NO uniform floor (LRC holds per-family but inf meas = 0). If it stabilizes,
a uniform floor is plausible. Also verifies the slope bound.
"""
import sys, random, itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

BAND = Fr(1, 14)
def sa(v): return [((Fr(k)+BAND)/v, (Fr(k+1)-BAND)/v) for k in range(v)]
def inter(A, B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]);hi=min(A[i][1],B[j][1])
        if lo<hi:r.append((lo,hi))
        if A[i][1]<B[j][1]:i+=1
        else:j+=1
    return r
def region(S):
    S=sorted(set(S)); a=sa(S[0])
    for v in S[1:]:
        a=inter(a,sa(v))
        if not a: return []
    return a
def meas(S): return sum(h-l for l,h in region(S))
def cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def prim(S): return reduce(gcd,S)==1
def Mval(S):
    """exact M(S)=max_t min_v ||vt||: max over the safe-region arc structure at increasing band."""
    # M = largest b such that {min||vt||>=b} nonempty; = max over candidate t* of min_v ||v t*||.
    # candidate t* are at intersections k/v +- ...; use the fact M is attained at a rational a/q,
    # q = lcm-ish; here approximate via fine search is unsafe. Use: M = max over v,w,k,l of the
    # crossing points. Cheap robust proxy: sample the region endpoints of the 1/14-safe set is empty
    # for tight; instead do a moderate rational scan of t=a/Q for Q=lcm capped. For SEARCH we only
    # need meas; M is reported only for the few finalists via fine float scan.
    return None
def Mfloat(S, G=200003):
    import numpy as np
    t=(np.arange(G)+0.5)/G
    F=np.full(G, 1.0)
    for v in S:
        fr=(v*t)%1.0
        d=np.minimum(fr,1.0-fr)
        F=np.minimum(F,d)
    i=int(F.argmax()); return F[i], t[i]

def scan(cap, ntrials, seed=0):
    random.seed(seed)
    poolsmall=[x for x in range(1,cap+1)]
    best=(9.9,None)
    # structured: near-imprimitive (c*subset + primitivizer) + random
    trials=0
    while trials<ntrials:
        trials+=1
        S=random.sample(poolsmall,13)
        if not any(s%14==0 for s in S): continue
        if not cov(S) or not prim(S): continue
        m=float(meas(S))
        if 0<m<best[0]: best=(m,tuple(sorted(S)))
    # also throw in scaled near-AP structured families
    for c in range(2,cap//13+1):
        for drop in range(1,14):
            base=[c*j for j in range(1,14) if j!=drop]
            for add in [drop, c*drop+1, c*drop-1, 2*c-1, c+1]:
                S=sorted(set(base+[add]))
                if len(S)==13 and cov(S) and prim(S) and max(S)<=cap:
                    m=float(meas(S))
                    if 0<m<best[0]: best=(m,tuple(sorted(S)))
    return best

print("="*100)
print(" inf meas(lonely S) over PRIMITIVE COVERING families, by magnitude cap")
print("="*100)
print(f"  {'cap':>5} {'min meas':>12} {'v_max':>7} {'M(S)':>9} {'2(M-1/14)/vmax':>16} {'family'}")
import numpy as np
for cap in [20, 30, 45, 70, 110]:
    b=scan(cap, 40000, seed=1)
    if b[1] is None: print(f"  {cap:>5}   (none)"); continue
    S=list(b[1]); Mf,tf=Mfloat(S); vmax=max(S)
    slope=2*(Mf-1.0/14)/vmax
    print(f"  {cap:>5} {b[0]:>12.6f} {vmax:>7} {Mf:>9.5f} {slope:>16.6f}  {S}")

print("\n"+"="*100); print(" SLOPE-BOUND CHECK: meas >= 2(M-1/14)/v_max  (should hold for all)"); print("="*100)
for S in [[1,2,3,4,5,6,8,9,10,11,12,13,14],
          [2,4,5,6,8,10,14,16,18,20,22,24,26],
          sorted(set([j for j in range(1,13)])|{182})]:
    Mf,_=Mfloat(S); m=float(meas(S)); vmax=max(S); sl=2*(Mf-1.0/14)/vmax
    print(f"  S={S}\n    meas={m:.6f}  M={Mf:.5f}  v_max={vmax}  2(M-1/14)/vmax={sl:.6f}  bound holds: {m+1e-9>=sl}")

print("\n"+"="*100); print(" READING"); print("="*100)
print("  If min meas DECREASES with cap (toward 0) and the min-family has LARGE v_max => NO uniform")
print("  floor (inf meas=0; LRC holds per-family). If it STABILIZES => a uniform floor is plausible.")
print("DONE.")
