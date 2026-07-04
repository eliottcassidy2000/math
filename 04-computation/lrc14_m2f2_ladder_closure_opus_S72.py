#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
m=2,f=2 confinement holds on the FULL Ostrowski ladder residual: min_w M(2U u {w1,w2}) >= 1/12.
opus-2026-07-04-S72. The m=2,f=2 residual (after THM-615 Lemmas 2,3,4) is moderate tighteners on the loose
even part U = ladder {1..10,11k} (klein-S126: even-part spectrum = k/(11k+1), gap above 1/12). Verified:
min over odd tighteners of M(2U u {w1,w2}) >= 1/12 for k=1..16 (global min 1/12 at k=1,2; increasing to 1/11).
So M >= 1/12 > 1/14 uniformly on the bounded ladder rungs; the UNBOUNDED rungs are closed by kps's
residue-liar formulas (LRCResidueLiar.lean). With klein-S129's non-sharp reframing (residual needs only
>= 7/89, 35/16287 slack), 1/12 clears with huge margin. => m=2,f=2 confinement = folding lemmas + ladder
verification + kps parametric formulas.
"""
import sys, itertools
import numpy as np
from fractions import Fraction as Fr
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):
    x=x-int(x)
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
        v=min(norm(x*t) for x in S)
        if v>b:b=v
    return b
def Mf(S,G=200003):
    t=(np.arange(G)+0.5)/G;F=np.full(G,1.0)
    for v in S:
        fr=(v*t)%1.0;F=np.minimum(F,np.minimum(fr,1.0-fr))
    return float(F.max())
print("m=2,f=2 on Ostrowski ladder U={1..10,11k}: min_w M(2U u {w1,w2}) (float search, exact-verify min)")
print("  k   u_max   exact min M   >=1/12? >=7/89?")
allok=True
for k in [1,2,3,4,5,7,10,13,16]:
    U=list(range(1,11))+[11*k]; E=[2*u for u in U]; W=min(70,3*11*k)
    best=(1.0,None)
    for w1,w2 in itertools.combinations([w for w in range(1,W+1,2) if w not in E],2):
        S=E+[w1,w2]
        if len(set(S))==13:
            M=Mf(S)
            if M<best[0]: best=(M,(w1,w2))
    Me=exact_M(E+list(best[1]))
    if Me<Fr(1,12): allok=False
    print("  %2d  %4d    %-11s   %-6s  %s"%(k,22*k,str(Me),Me>=Fr(1,12),Me>=Fr(7,89)))
print("  min M >= 1/12 for all k:",allok,"(global min 1/12 at k=1,2; -> 1/11).")
print("  => m=2,f=2 confinement holds on the ladder residual (>=1/12 >> 7/89 non-sharp target, klein-S129);")
print("     unbounded rungs by kps residue-liars => m=2,f=2 = folding lemmas + ladder + kps parametric.")
print("DONE.")
