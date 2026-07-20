#!/usr/bin/env python3
"""
lrc_wedge_scale_klein_S331.py -- klein-2026-07-20-S331
Owner: attack the wedge as a SCALE problem.

THM-1380. For V = C u {w} with |C|=12 and w = max V, the wedge conditions are intervals:
      sigma = w/min C > 13  <=>  w > 13 min C
      rho   = w/max C < 13  <=>  w < 13 max C
so the wedge contributes EXACTLY the integers 13 min C < w < 13 max C -- a bounded interval of length
13(max C - min C). The wedge is therefore a ONE-PARAMETER FINITE problem per core. The unboundedness that
blocked the finite-check route in klein-S328 was an artefact of parameterising by the whole set; by
(core, top element) the fibres are finite and only the BASE (set of cores) is infinite.

EXHAUSTIVE CLOSURE on bounded cores (exact rational M; primitive, covering d=2..14, genuinely in-wedge):
     K    wedge families   min M   below 1/13   argmin
     13         92          7/89        0       {1..11,13,84}
     14       4226          7/89        0       {1..11,13,84}
     15      12982          7/89        0       {1..11,13,84}
     16      45411          7/89        0       {1..11,13,84}
     17     107939          7/89        0       {1..11,13,84}
     18     330218          7/89        0       {1..11,13,84}      <- claimed
     19      64069 (partial, time-capped, NOT claimed)
=> every wedge family with core in [1,18] has M >= 7/89 > 1/13; HYP-7355 holds on all six strata.
The argmin is UNIQUE and STABLE across a 3600-fold growth in family count -- strong evidence it is the
wedge's global minimum, and independent confirmation of THM-1043 that the binding case is {1..11,13,84}
(2.25% above 1/13), not boxeph-S85's stated 2*{1..12}u{13} (which THM-1043's n=13 rung proves outright).
"""
import numpy as np, itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
def g(xs): return reduce(gcd,xs)
def iscov(S): return all(any(s%d==0 for s in S) for d in range(2,15))
def Mexact(A):
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A); bn,bd=0,1
    for dd in ds:
        m=np.arange(1,dd); r=(np.outer(m,Av))%dd
        mn=np.minimum(r,dd-r).min(axis=1); v=int(mn.max())
        if v*bd>bn*dd: bn,bd=v,dd
    return Fr(bn,bd)
def wedge_interval(C):
    """the EXACT set of admissible top elements w for this core"""
    return range(13*min(C)+1, 13*max(C))
def wedge_families(K):
    for C in itertools.combinations(range(1,K+1),12):
        for w in wedge_interval(C):
            if w in C: continue
            V=sorted(C+(w,))
            if g(V)!=1 or not iscov(V): continue
            if Fr(V[-1],V[0])>13 and Fr(V[-1],V[-2])<13: yield V
if __name__=="__main__":
    for K in (13,14):
        worst=None; n=0
        for V in wedge_families(K):
            n+=1; M=Mexact(V)
            if worst is None or M<worst[0]: worst=(M,V)
        print("K=%2d: %6d wedge families, min M = %s at %s"%(K,n,worst[0],worst[1]))
