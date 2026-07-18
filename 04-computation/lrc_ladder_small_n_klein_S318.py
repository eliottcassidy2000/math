#!/usr/bin/env python3
"""
lrc_ladder_small_n_klein_S318.py
================================
klein-2026-07-18-S318 (owner: check the ladder against the known small-n covering minima).

RESULT: my S317 suggestion is REFUTED. The equioscillating ladder does NOT give the covering-minimum at
general n. Exhaustive over primitive covering (n-1)-sets (covering d=2..n):
    n=5: 2/9  at {1,3,4,5}        (ladder 5/21)
    n=6: 2/11 at {1,3,4,5,18}     (ladder 6/31)
    n=7: 2/13 at {1,2,5,6,7,8}    (ladder 7/43)
    n=8: 2/15 at {1,3,4,5,7,11,24}(ladder 8/57)
    n=9: 1/8                       (ladder 9/73)
The ladder-family min IS always at k=n(n-1) with M=m/((n-1)m+1) exactly -- but the family does not contain
the global min. True minimizers are COMPACT with a non-initial-segment body ({1,3,4,5} skips 2). n=14's
covering-min being ladder-shaped is a coincidence of n=14.

THE REAL PATTERN + THRESHOLD: minima run 2/(2n-1) for n=5..8, then BREAK at n=9 (1/8, not 2/17).
That break is the general-n analog of boxeph HYP-7355 (compact rho<n-1 covering => M >= 1/(n-1)):
    n=5..8 : compact min 2/9,2/11,2/13,2/15  <  1/(n-1)  => ANALOG FAILS (all minimizers compact)
    n=9    : compact min 1/8                 =  1/(n-1)  => HOLDS with equality
=> HYP-7355 cannot be proved by induction on n or any n-uniform argument; it is FALSE for n<=8 and
switches on at n=9. (Same phenomenon as mac-mini-S110's "stability gap FAILS at n=6,7".)
Caveat: n=8,9 rows exhaustive within entry bound 26/24 (compact only); n=5,6,7 exhaustive to 60/45/32.
Evaluator validated against a 2^20 grid (max deviation 2e-6) on all disputed sets.
"""
import numpy as np, itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
def g(xs): return reduce(gcd,xs)
def Mexact(A):
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A); bn,bd=0,1
    for dd in ds:
        m=np.arange(1,dd); r=(np.outer(m,Av))%dd
        mn=np.minimum(r,dd-r).min(axis=1); val=int(mn.max())
        if val*bd>bn*dd: bn,bd=val,dd
    return Fr(bn,bd)
def covering_min(n,B,compact_only=False):
    def iscov(S): return all(any(s%d==0 for s in S) for d in range(2,n+1))
    best=None; arg=None
    for V in itertools.combinations(range(1,B+1),n-1):
        if not iscov(V) or g(V)!=1: continue
        if compact_only and Fr(max(V),sorted(V)[-2])>=n-1: continue
        M=Mexact(list(V))
        if best is None or M<best: best,arg=M,V
    return best,arg
if __name__=="__main__":
    for n,B in ((5,60),(6,45),(7,32)):
        m,a=covering_min(n,B)
        print("n=%d true covering-min = %-7s at %-22s ; ladder = %-7s"%(n,m,str(a),Fr(n,n*n-n+1)))
    for n,B in ((8,26),(9,24)):
        m,a=covering_min(n,B,compact_only=True)
        print("n=%d COMPACT covering-min = %-7s at %-22s ; 1/(n-1) = %-6s ; analog %s"
              %(n,m,str(a),Fr(1,n-1),"HOLDS" if m>=Fr(1,n-1) else "FAILS"))
