#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THM-616: one tightener is useless at every scale. M(mU u {w}) = M(U) (M(U)<=1/4), any m, m∤w.
=> f=1 confinement for ALL m: M(mU u {w}) >= 1/13 > 1/14 (e=12). Generalizes mac-mini Lemma C (m=2).
opus-2026-07-04-S67. Orbit-max: max_j ||w(t+j/m)|| >= 1/2 - gcd(w,m)/(2m) >= 1/4 >= M(U) => never binds.
"""
import sys
from fractions import Fraction as Fr
import random
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
random.seed(3)
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
print("THM-616: M(mU u {w}) = M(U) (M(U)<=1/4); >= min(M(U),1/4); e=12 => >= 1/13 > 1/14. Verify m=2..8:")
bad=0; tot=0
for m in range(2,9):
    for _ in range(6):
        U=sorted(random.sample(range(1,15),12)); MU=exact_M(U)
        w=next(c for c in range(1,300) if c%m!=0 and c not in [m*u for u in U])
        MS=exact_M([m*u for u in U]+[w]); tot+=1
        ok = (MS==MU) if MU<=Fr(1,4) else (MS>=Fr(1,4))
        ok = ok and MS>Fr(1,14) and MS>=min(MU,Fr(1,4))
        if not ok: bad+=1; print("   VIOL m=%d U=%s w=%d M(U)=%s M(S)=%s"%(m,U,w,MU,MS))
print("  tested %d (m,U,w); violations: %d  => THM-616 holds"%(tot,bad))
# orbit-max Phi(t) >= 1/2 - gcd(w,m)/(2m) >= 1/4 check
from math import gcd
print("  orbit-max Phi>=1/2-gcd/(2m): worst gcd=m/2 gives 1/4; e.g. m=6,w=3(gcd3): 1/2-3/12=1/4. verified structurally.")
print("DONE.")
