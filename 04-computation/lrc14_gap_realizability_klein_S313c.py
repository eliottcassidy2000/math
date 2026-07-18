#!/usr/bin/env python3
"""
lrc14_gap_realizability_klein_S313c.py
======================================
klein-2026-07-17-S313c (owner: rule out (3,38) by integer realizability; keep working the full bound).

OUTCOME: (3,38) is NOT ruled out. Three routes tried, all fail -- documented so they aren't repeated.
  1. RESIDUE COVERING rules out nothing: B_c={k:|ck|_q<=val} admits greedy covers of size 7-10 <= 12
     at EVERY admissible (val,q), q<=149.
  2. MONOTONICITY is the wrong direction: A subset Allowed(a) => M(A) >= M(Allowed(a)) = 3/38, a LOWER
     bound; detecting a gap violation needs an UPPER bound M < 2/25.
  3. SEARCH finds no violation: random and greedy-descent 12-subsets of Allowed(a) bottom out at
     min M = 1/11 = 0.0909 >> 2/25 = 0.08.
NEW STRUCTURE: writing q = 13 val - s, the gap is exactly {val/(13val-s) : 1 <= s < val/2}. Since
val -> infinity is admissible (s fixed), the in-gap denominator list is finite IFF max(A) is bounded:
CRUX (C) is NOT a finite denominator check in general. The 3/38 -> >=1/11 jump on descending to 12
elements indicates the right route is an INVERSE/RIGIDITY theorem, not denominator elimination.
"""
import numpy as np
from fractions import Fraction as Fr
from math import gcd
lo,hi=Fr(1,13),Fr(2,25)
def Mexact(A):
    """exact M = max_t min_v ||v t||; maximizer sits at t=K/(v_i+v_j) (THM-999 Lemma A, general-M form)."""
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A); bn,bd=0,1
    for d in ds:
        m=np.arange(1,d); r=(np.outer(m,Av))%d
        mn=np.minimum(r,d-r).min(axis=1); val=int(mn.max())
        if val*bd > bn*d: bn,bd=val,d
    return Fr(bn,bd)
def gap_pairs(qmax): return [(v,q) for v in range(1,14) for q in range(2,qmax+1) if lo<Fr(v,q)<hi]
print("gap parametrization q = 13*val - s, 1 <= s < val/2:")
for val in range(3,9):
    ss=[s for s in range(1,(val+1)//2+1) if s<Fr(val,2)]
    print("  val=%d -> q =%s"%(val,[13*val-s for s in ss]))
print("\ngreedy descent from Allowed(a) down to 12 elements:")
for (val,q,B) in [(3,38,25),(3,38,30),(4,51,35)]:
    best=None
    for a in range(1,q//2+1):
        if gcd(a,q)!=1: continue
        A=[v for v in range(1,B+1) if min((v*a)%q,(q-(v*a)%q))>=val]
        if len(A)<12: continue
        while len(A)>12:
            cand=None
            for x in A:
                B2=[y for y in A if y!=x]; M2=Mexact(B2)
                if cand is None or M2<cand[0]: cand=(M2,B2)
            A=cand[1]
        M=Mexact(A)
        if best is None or M<best: best=M
    print("  (val,q,B)=(%d,%d,%d): greedy-min 12-subset M = %s = %.5f  %s"%(val,q,B,best,float(best),
          "IN GAP" if lo<best<hi else ">= 2/25 (no violation)"))
