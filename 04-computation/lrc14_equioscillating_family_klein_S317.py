#!/usr/bin/env python3
"""
lrc14_equioscillating_family_klein_S317.py
==========================================
klein-2026-07-18-S317 (owner: DFS the "169" thread for inspiration).

THE EQUIOSCILLATING LADDER.  Body {1..b}, killer k == -1 (mod q), witness t = m/q: then v=1 and v=k both
sit at distance m/q (equioscillation), and the body is safe iff q >= m(b+1). Smallest such q with b=12:
        V_m = {1..12} u {13m},     M(V_m) = m/(13m+1)     [exact; verified m=1..16,28,42]
M is INCREASING in m -> 1/13 from below; V_m is primitive for all m; and
        V_m covering  <=>  14 | 13m  <=>  14 | m.
CONSEQUENCE -- why the covering-minimum is 14/183:
  m=1  killer 13   M=1/14      the tight AP {1..13} (LRC(14) extremal), non-covering
  m=13 killer 169  M=13/170    BELOW the covering-min, but non-covering (169=13^2 has no factor 14)
  m=14 killer 182  M=14/183    COVERING -> THE COVERING-MINIMUM (THM-724)
So the covering-min is the FIRST RUNG OF A MONOTONE LADDER AT WHICH COVERING SWITCHES ON; everything below
is excluded solely by 14 nmid m. 169 marks the LAST NON-COVERING RUNG before it (and 169 = 13^2 is an
accident: 14*12+1 = 169 too, which is why the n=13 deep well {1..11,168} looked structural).
Re-derives THM-1014's tower from the parametric side: exceptions are V_m with 14|m = {1..12,182m'},
rho = 13m/12 forced past 13 by m >= 14.
NB the OTHER 169 is real: 169/2352 = (1/12)(13/14)^2 is the balance ladder's per-step loss at j=2, the
reason clustering broke telescoping -- fixed by the additive tail certificate (THM-1015).
"""
import numpy as np
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
def g(xs): return reduce(gcd,xs)
def iscov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def Mexact(A):
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A); bn,bd=0,1
    for dd in ds:
        m=np.arange(1,dd); r=(np.outer(m,Av))%dd
        mn=np.minimum(r,dd-r).min(axis=1); val=int(mn.max())
        if val*bd>bn*dd: bn,bd=val,dd
    return Fr(bn,bd)
if __name__=="__main__":
    print("  m | killer | predicted | exact     | covering")
    for m in list(range(1,17))+[28,42]:
        V=sorted(set(list(range(1,13))+[13*m]))
        if len(V)!=13: continue
        print(" %3d | %6d | %-9s | %-9s | %s"%(m,13*m,Fr(m,13*m+1),Mexact(V),iscov(V)))
    print("\ncovering-min = first rung with 14|m  =>  m=14, killer 182, M = 14/183")
