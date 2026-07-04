#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
ORBIT-COVERING: f tighteners are USELESS unless they can COVER the m-orbit at the g-argmax.
opus-2026-07-04-S69. Generalizes THM-616 (f=1) to f < 7m/(m+7) (~ f<=6 large m).
At the g-argmax t0, M(mU u F) = max_t min(g(t),Phi_F(t)), Phi_F(t0)=max_j min_{w in F}||w(t0+j/m)||.
Each w is danger (<1/14) at <= m/7 + gcd(w,m) of the m orbit points; if sum < m the orbit is uncovered
=> some point all-safe => Phi_F(t0)>=1/14 => M >= min(M(U),1/14) > 1/14 (M(U)>=1/(e+1)). f coprime => need
f(m+7)<7m. m=2 => f<=1 (=THM-616); the HARD m=2,f=2 is exactly f=m (orbit just coverable = the parity case).
Also: m=2,f=2 confinement HOLDS on the Ostrowski ladder U={1..10,11k} (min_w M(2U u {w1,w2}) >= 1/12).
"""
import sys, itertools
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
print("(1) f<=6 tighteners are USELESS: M(mU u F)=M(U)=1/(e+1) (orbit uncovered):")
bad=0
for m in [8,10,12,20]:
    for f in [2,3,4,5,6]:
        e=13-f; U=list(range(1,e+1)); F=[]; c=1
        while len(F)<f:
            c+=1
            if c%m and c not in F and c not in [m*u for u in U]: F.append(c)
        S=[m*u for u in U]+F; M=exact_M(S)
        if M!=Fr(1,e+1): bad+=1; print("   dev m=%d f=%d M=%s vs 1/%d"%(m,f,M,e+1))
print("   M != 1/(e+1) deviations:",bad,"(0 => tighteners useless for f<=6, these U)")
print("(2) m=2,f=2 confinement holds on Ostrowski ladder U={1..10,11k}: min_w M(2U u {w1,w2}) >= 1/12:")
for k in [1,2,3,4,5,6]:
    U=list(range(1,11))+[11*k]; best=Fr(1)
    for w1,w2 in itertools.combinations(range(1,42,2),2):
        S=[2*u for u in U]+[w1,w2]
        if len(set(S))==13:
            M=exact_M(S)
            if M<best: best=M
    print("   k=%d M(U)=%-7s min_w M=%-7s >=1/12? %s"%(k,str(exact_M([2*u for u in U])),str(best),best>=Fr(1,12)))
print("DONE.")
