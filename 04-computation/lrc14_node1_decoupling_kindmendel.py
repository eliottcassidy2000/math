#!/usr/bin/env python3
"""Node 1 swing: explicit finite-V decoupling for the boundary core {b,2b,..,12b,V}.
L >= (6/7)mu_12 - b*r_12/(7V), positive for V > kappa*b, kappa=r_12/(6 mu_12).
kind-mendel-2026-06-22-S3."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
def gl(xs): return reduce(lambda a,b:a*b//gcd(a,b),[x for x in xs if x],1)
def nrm(y): f=y-floor(y); return min(f,1-f)

def lonely_exact(S, thr=F(1,14)):
    S=[s for s in S if s]
    bps=set([F(0),F(1)])
    for s in S:
        for n in range(s+1):
            for sign in (F(-1),F(1)):
                x=F(n,s)+sign*thr/s
                if 0<=x<=1: bps.add(x)
    bps=sorted(bps); tot=F(0); flags=[]
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2; ok=all(nrm(s*mid)>=thr for s in S)
        flags.append(ok)
        if ok: tot+=b-a
    arcs=0; prev=False
    for ok in flags:
        if ok and not prev: arcs+=1
        prev=ok
    return tot,arcs

mu12,r12=lonely_exact(list(range(1,13)))
kappa=F(r12,6*mu12)
print(f"mu_12 = L({{1..12}}@1/14) = {mu12} = {float(mu12):.6f}")
print(f"r_12  (arc count)        = {r12}")
print(f"kappa = r_12/(6 mu_12)   = {kappa} = {float(kappa):.2f}   => need V > kappa*b")
print(f"limit floor (6/7)mu_12   = {float(F(6,7)*mu12):.6f}\n")

# confirm scale-invariance + arc multiplication: L({b..12b}) = mu_12, arcs = 12b
for b in [1,2,3]:
    m,a=lonely_exact([j*b for j in range(1,13)])
    print(f"  L({{{b}..{12*b}}})={float(m):.6f} (==mu_12:{m==mu12}) arcs={a} (==12b={12*b}:{a==12*b})")

print("\n=== Verify boundary core: L({b,..,12b,V}) >= (6/7)mu_12 - b*r12/(7V),  V>kappa*b => L>0 ===")
print(f"{'b':>2} {'V':>5} {'V/b':>6} {'V>kb':>5} {'L_exact':>9} {'lowerbnd':>9} {'bnd>0':>5} {'L>0':>4} {'valid':>5}")
for b in [1,2,3]:
    kb=float(kappa*b)
    for V in sorted({14*m for m in range(1,40)} | {14*m for m in [40,60,100]}):
        if V <= 12*b or V in {j*b for j in range(1,13)}: continue
        if V < kb-30 or V > kb+200: 
            if not (V in (70,84,140,420,840) and b==1): continue
        S=[j*b for j in range(1,13)]+[V]
        if len(set(S))<13: continue
        L,_=lonely_exact(S)
        bnd=F(6,7)*mu12-F(b*r12,7*V)
        print(f"{b:>2} {V:>5} {V/b:>6.1f} {str(V>kappa*b):>5} {float(L):>9.6f} {float(bnd):>9.6f} "
              f"{str(bnd>0):>5} {str(L>0):>4} {str(L>=bnd-F(1,10**12)):>5}")

print("\n=== The residual GAP: comparable V (V <= kappa*b) makes the elementary bound vacuous ===")
for b,V in [(1,14),(1,28),(2,28),(3,14),(1,42)]:
    S=[j*b for j in range(1,13)]+[V]
    if len(set(S))<13: continue
    L,_=lonely_exact(S); bnd=F(6,7)*mu12-F(b*r12,7*V)
    print(f"b={b} V={V:>3} V/b={V/b:>5.1f}: L_exact={float(L):.6f} (>0:{L>0})  bnd={float(bnd):+.4f} (vacuous:{bnd<=0})")
