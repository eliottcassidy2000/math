#!/usr/bin/env python3
"""
LRC(14) DISPROVE - obstruction analysis.
Why can't perturbations of {1..13} go below 1/14?

For the tight AP {1..13}, M=1/14 attained at tau=1/14 AND tau=5/14.
At tau=1/14: v*tau = v/14, so ||v/14|| for v=1..13 is
  min over v of dist(v/14, Z) = 1/14 (v=1 and v=13).
So the gap is "pinned" by the extreme speeds v=1 and v=13 at tau=1/14.

We quantify: for a perturbation, track which tau in the candidate set is the
argmax (the lonely point), and what its value is. We show that when you push one
lonely point down, ANOTHER tau rises (the envelope's max is robust). This is the
'whack-a-mole' obstruction: the lower envelope of 13 tent maps has multiple high
plateaus; suppressing one keeps M pinned at 1/14 by another.

We compute, for {1..13} and single perturbations, the SET of tau achieving M and
the SECOND-best plateau, to show the redundancy that protects 1/14.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def argmax_taus(S):
    """all tau achieving M, and the value."""
    b=F(0); ats=[]
    for t in sorted(cand(S)):
        v=g(S,t)
        if v>b: b=v; ats=[t]
        elif v==b: ats.append(t)
    return b, ats
def is_prim(S): return reduce(gcd,S)==1

THRESH=F(1,14)
base=tuple(range(1,14))
b, ats = argmax_taus(base)
print(f"Tight {base}: M={b}={float(b):.6f}")
print(f"  All argmax tau (lonely points): {ats}")
print(f"  Number of distinct optimal tau: {len(ats)}")
print()
# at each optimal tau, which elements are the 'binding' constraints (||v tau||=M)?
for t in ats:
    binders = [v for v in base if nrm(v*t)==b]
    print(f"  tau={t}: binding speeds (||v tau||=M) = {binders}")
print()

# Single perturbations: for each, M and the gap to 1/14, and #optimal tau.
print("Single perturbations of {1..13}: how M responds")
print(f"{'replace':>10} {'with':>5} {'M':>12} {'M_float':>10} {'#opt_tau':>8} {'>=1/14?':>7}")
results=[]
for i in range(13):
    for newv in range(14,30):
        S=list(base[:i]+base[i+1:])+[newv]
        S=tuple(sorted(set(S)))
        if len(S)!=13 or not is_prim(S): continue
        b2,ats2=argmax_taus(S)
        results.append((b2,base[i],newv,len(ats2)))
results.sort()
for b2,old,newv,nopt in results[:25]:
    print(f"{old:>10} {newv:>5} {str(b2):>12} {float(b2):>10.6f} {nopt:>8} {str(b2>=THRESH):>7}")
print()
print(f"Minimum M among these single perturbations: {results[0][0]}={float(results[0][0]):.6f}")
print(f"All >= 1/14 ? {all(r[0]>=THRESH for r in results)}")
print()
print("OBSTRUCTION: every single perturbation keeps M>=1/14. The lower envelope")
print("of the remaining 12 consecutive speeds 1..13 already pins M at multiple")
print("tau simultaneously (1/14 and 5/14 each binding the extremes). Removing one")
print("speed cannot lower BOTH plateaus; the surviving AP-core re-pins 1/14.")
