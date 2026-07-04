#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
IS inf_U gap(U) > 0 ?  (m=2,f=2 confinement as a uniform gap). opus-S64.
gap(U)=min_{odd w1,w2}(M(2U u {w1,w2})-1/14). The most likely place for gap->0 is NEAR the tight even
block 2*{1..13} (M=1/14): U = {1..13}\{i,j} (remove 2 even runners), add 2 odd tighteners. Compute the
MIN gap over all pairs {i,j} (and best odd tighteners); exact-verify the minimizer. If min gap is bounded
away from 0 => inf>0 plausible. If some {i,j} -> 0 => inf=0 (statement false).
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import numpy as np
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def Mf(S,G=60013):
    t=(np.arange(G)+0.5)/G;F=np.full(G,1.0)
    for v in S:
        fr=(v*t)%1.0;F=np.minimum(F,np.minimum(fr,1.0-fr))
    return float(F.max())
def norm(x):
    x=x-int(x);
    if x<0:x+=1
    return min(x,1-x)
def exact_M(S):
    S=sorted(set(S));cands=set()
    for v in S:
        for k in range(v):cands.add(Fr(2*k+1,2*v))
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for den in (S[i]+S[j],abs(S[i]-S[j])):
                if den:
                    for s in range(den):cands.add(Fr(s,den))
    b=Fr(0)
    for t in cands:
        val=min(norm(v*t) for v in S)
        if val>b:b=val
    return b
def prim(S):return reduce(gcd,S)==1
T=1.0/14
def best_gap(U, wmax=45):
    E=[2*u for u in U];odds=[w for w in range(1,wmax+1,2) if w not in E]
    best=(1.0,None)
    for w1,w2 in itertools.combinations(odds,2):
        S=E+[w1,w2]
        if not prim(S):continue
        M=Mf(S)
        if M>=T-1e-5 and M-T<best[0]:best=(M-T,(w1,w2))
    return best
print("MIN gap over near-block U={1..13}\{i,j} (E=2U), best odd tighteners (float; minimizer exact-verified)")
print("  %-16s %10s %14s" % ("removed {i,j}","gap(float)","best (w1,w2)"))
rows=[]
for i,j in itertools.combinations(range(1,14),2):
    U=[x for x in range(1,14) if x not in (i,j)]
    g,ws=best_gap(U)
    rows.append((g,(i,j),ws,U))
rows.sort()
for g,ij,ws,U in rows[:10]:
    print("  %-16s %10.6f %14s" % (str(ij),g,str(ws)))
print("  ... (showing 10 smallest of %d pairs)"%len(rows))
gmin,ij,ws,U=rows[0]
Se=sorted([2*u for u in U]+list(ws))
Me=exact_M(Se)
print(f"\n  SMALLEST gap: removed {ij}, U={U}")
print(f"    S={Se}  M(S) exact={Me}  gap exact={float(Me-Fr(1,14)):.6f}  =1/14? {Me==Fr(1,14)}  primitive={prim(Se)}")
print(f"    1/84={float(Fr(1,84)):.6f}")
# dilation check: does the smallest-gap family keep its gap under dilation (=> unbounded extremizer)?
print("\n  dilation of the smallest-gap shape (c*U), best tighteners re-searched:")
for c in [1,2,3]:
    Uc=[c*u for u in U]; g,ws=best_gap(Uc, wmax=min(60,12*max(Uc)))
    print(f"    c={c}: u_max={max(Uc):>3} gap(float)={g:.6f} (w1,w2)={ws}")
print("\nREADING: if the min near-block gap is bounded away from 0 (~>0.005) and stays so under dilation,")
print("  inf_U gap(U) > 0 is plausible; if it -> 0 for some {i,j}, inf=0 (uniform gap FALSE).")
print("DONE.")
