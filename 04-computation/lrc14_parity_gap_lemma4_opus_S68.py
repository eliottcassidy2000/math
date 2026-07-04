#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THM-615 Lemma 4 (PARITY GAP): the ODD-ness of the two tighteners forces a gap on the confinement core.
opus-2026-07-04-S68. For 11-runner U with M(U)>1/12, I0 = global-max component of {g_E>=1/12}, width L:
  same-mode extremity on I0 needs I0 inside a w1-half-integer arc AND a w2-integer arc; their centers
  c1=(2l+1)/(2w1), c2=k/w2 satisfy c1-c2 = [(2l+1)w2 - 2k w1]/(2 w1 w2), numerator ODD (both w odd) => >=1
  => |c1-c2|>=1/(2 w1 w2). But both arcs contain I0 => |c1-c2| <= (w1+w2)/(12 w1 w2) - L. So L <= (w1+w2-6)/(12 w1 w2).
  => if L > (w1+w2-6)/(12 w1 w2) then M >= 1/12; in particular w1+w2<=6 => M>=1/12 unconditionally.
Disposes the SMALL end (all w1+w2<=6). With Lemma 3 (large end), residual = moderate tighteners on near-AP.
"""
import sys
from fractions import Fraction as Fr
import random
import numpy as np
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
random.seed(7)
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
def L_width(U,band=1/12,G=400000):
    t=(np.arange(G)+0.5)/G; g=np.full(G,1.0)
    for u in U:
        fr=(2*u*t)%1.0; g=np.minimum(g,np.minimum(fr,1.0-fr))
    hi=g>=band; i0=int(g.argmax()); lo=i0; hh=i0
    while lo>0 and hi[lo-1]: lo-=1
    while hh<G-1 and hi[hh+1]: hh+=1
    return (hh-lo+1)/G

print("(A) w1+w2 <= 6 ({1,3},{1,5}) => M(2U u {w1,w2}) >= 1/12, all 11-runner U (M(U)>1/12):")
bad=tot=0
for _ in range(80):
    U=sorted(random.sample(range(1,20),11))
    if exact_M([2*u for u in U])<=Fr(1,12): continue
    for w1,w2 in [(1,3),(1,5)]:
        S=[2*u for u in U]+[w1,w2]
        if len(set(S))!=13: continue
        tot+=1
        if exact_M(S)<Fr(1,12): bad+=1
print("   tested %d, M<1/12 violations: %d"%(tot,bad))
print("(B) general L-condition L>(w1+w2-6)/(12 w1 w2) => M>=1/12 (spot check where it fires):")
bad2=tot2=0
for _ in range(60):
    U=sorted(random.sample(range(1,15),11)); MU=float(exact_M([2*u for u in U]))
    if MU<=1/12+1e-4: continue
    L=L_width(U)
    for w1,w2 in [(3,5),(5,7),(3,7),(1,9),(5,9),(7,9),(3,11)]:
        thr=(w1+w2-6)/(12*w1*w2)
        if L>thr:
            tot2+=1
            if exact_M([2*u for u in U]+[w1,w2])<Fr(1,12)-Fr(1,10**6): bad2+=1
print("   tested %d (L>thr cases), M<1/12 violations: %d"%(tot2,bad2))
print("\nRESIDUAL after Lemmas 2,3,4 + THM-616: MODERATE tighteners (6<w1+w2, both<=u_max/(6(M(U)-1/12)))")
print("  on the near-AP corner; for PRIMITIVE families near-AP => bounded (rigidity) => finite. The parity")
print("  gap (odd tighteners => odd=even coincidence >=1) is the new mechanism closing the small end.")
print("DONE.")
