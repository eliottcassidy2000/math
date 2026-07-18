#!/usr/bin/env python3
"""residual_hunt_kps_S128c50.py -- kind-pasteur S128 cont.50.
(A) structural analysis of the 1/7 margin-floor minimizer: witness t*, obstruction.
(B) TARGETED residual hunt: primitive covering-gap-tight family = counterexample to THM-995(VII).
    Search structured shapes {1..11,a,b}, {1..10,a,b,c}, {1..12,a} (compressed small cluster + gap
    to large speeds). Grid-filter tight, exact-verify, classify. If ANY is covering+gap+max>=23+tight,
    the residual conjecture is REFUTED (real discovery). Else strong support.
"""
import sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)
LAM=F(1,14)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def Mgrid(V,Q=14001):
    best=0.0; iq=1.0/Q
    for p in range(1,Q):
        t=p*iq; m=2.0
        for v in V:
            d=abs(v*t-round(v*t))
            if d<m:
                m=d
                if m<=best: break
        if m>best: best=m
    return best
def M_exact(V):
    cand=set()
    for v in V:
        for a in range(1,2*v): cand.add(F(a,2*v))
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            for s in (V[i]+V[j],abs(V[i]-V[j])):
                if s:
                    for a in range(1,s): cand.add(F(a,s))
    b=F(0); bt=None
    for t in cand:
        if 0<t<1:
            m=min(nd(v*t) for v in V)
            if m>b: b=m; bt=t
    return b,bt
def covering(V):
    for q in range(2,15):
        if min(nd(v*F(1,q)) for v in V)>=F(1,14): return False,q
    return True,None
def gap(V): return any(a>13*b for a in V for b in V)

print("== (A) structure of the 1/7 minimizer ==")
Vmin=[42,66,96,108,150,228,229,247,375,377,396,414,552]
M,tstar=M_exact(Vmin)
print("  V=%s"%Vmin)
print("  M=%s=%.6f, witness t*=%s"%(M,float(M),tstar))
if tstar:
    diststar=sorted((v, nd(v*tstar)) for v in Vmin)
    tight_speeds=[v for v,d in diststar if d==M]
    print("  speeds achieving the min (||v t*||=%s): %s"%(M,tight_speeds))
    print("  t* denominator:", tstar.denominator, " (1/7-structure: 7 | denom? %s)"%(tstar.denominator%7==0))
print()
print("== (B) targeted residual hunt (compressed cluster + gap, tight?) ==")
found=[]
tested=0
# shape 1: {1..11, a, b}
base11=list(range(1,12))
for a in range(12, 70):
    for b in range(a+1, 90):
        V=base11+[a,b]
        if len(set(V))!=13: continue
        if max(V)<23: continue
        tested+=1
        if Mgrid(V) <= 1/14 + 5e-5:
            Me,_=M_exact(V)
            if Me==LAM:  # tight
                cov,q=covering(V); g=gap(V)
                cls = "COVERING+gap RESIDUAL!!" if (cov and g and max(V)>=23) else ("sieve q=%s"%q if not cov else "no-gap")
                found.append((V,cls))
                print("  TIGHT {1..11,%d,%d}: %s"%(a,b,cls))
# shape 2: {1..10, a, b, c}
base10=list(range(1,11))
for a in range(11, 45):
    for b in range(a+1, 55):
        for c in range(b+1, 70):
            V=base10+[a,b,c]
            if len(set(V))!=13 or max(V)<23: continue
            tested+=1
            if Mgrid(V) <= 1/14 + 5e-5:
                Me,_=M_exact(V)
                if Me==LAM:
                    cov,q=covering(V); g=gap(V)
                    cls = "COVERING+gap RESIDUAL!!" if (cov and g and max(V)>=23) else ("sieve q=%s"%q if not cov else "no-gap")
                    found.append((V,cls))
                    print("  TIGHT {1..10,%d,%d,%d}: %s"%(a,b,c,cls))
print()
print("  tested %d structured families; tight found: %d"%(tested,len(found)))
resid=[V for V,c in found if "RESIDUAL" in c]
if resid:
    print("  *** RESIDUAL COUNTEREXAMPLES (primitive covering gap tight):", resid, "***")
else:
    print("  >>> NO primitive covering-gap tight family in the structured search: every tight one is sieve-covered or no-gap.")
    print("  >>> THM-995(VII) residual conjecture SUPPORTED on the compressed-cluster+gap shapes.")
print("DONE")
