#!/usr/bin/env python3
"""
lrc14_gap_reduction_klein_S313b.py
==================================
klein-2026-07-17-S313 (owner: prove the analytic stability bound non-AP => M >= 2/25).

That bound is CRUX (C) and OPEN. THM-1002 proves the BOUNDED case and pins the obstruction:
  (1) PAIR-SUM DENOMINATOR BOUND: the maximizer t*=a/q (lowest terms) satisfies q | (v_i+v_j),
      so q <= 2*max(A); and M(A) = val/q with val = min_v |v a|_q.
  (2) M=val/q in (1/13,2/25) <=> 12.5 val < q < 13 val; val=1 -> q in (12.5,13) EMPTY,
      val=2 -> q in (25,26) EMPTY. So a violation needs val>=3, q>=38.
  (3) THEOREM: max(A) <= 18 => q <= 36 < 38 => M = 1/13 or M >= 2/25.
  (4) NO EXTENSION via residues: the discrete covering condition (B_c={k:|ck|_q<=val} must cover Z/q
      with <=12 sets) is SATISFIABLE at every admissible (val,q) (greedy 7-10 <= 12) -> the obstruction
      is integer realizability, i.e. the hard LRC content.
  ERRATUM: the claim "crossing forces v == +-1 mod q" is FALSE (only when val=1); refuted 332/400.
Evidence: 0 gap violations in ~174k primitive 12-sets (119,998 with max<=36; 53,740 with max<=60).
"""
import numpy as np, random
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
def g(xs): return reduce(gcd,xs)
def Mexact(A):
    """exact M via the pair-sum ruler: maximizer sits at t=K/(v_i+v_j)."""
    ds=sorted({A[i]+A[j] for i in range(len(A)) for j in range(i,len(A))})
    Av=np.array(A); best=Fr(0); arg=None
    for d in ds:
        m=np.arange(1,d); r=(np.outer(m,Av))%d
        mn=np.minimum(r,d-r).min(axis=1)
        k=int(mn.argmax()); val=int(mn[k])
        if val>0 and Fr(val,d)>best: best=Fr(val,d); arg=(k+1,d)   # (witness numerator, denominator)
    return best,arg
lo,hi=Fr(1,13),Fr(2,25)
pairs=[(v,q) for v in range(1,12) for q in range(2,150) if lo<Fr(v,q)<hi]
print("in-gap (val,q), q<=149:",pairs[:12]); print("smallest q =",min(q for _,q in pairs),"=> max(A)>=19 needed")
print("THEOREM: max(A)<=18 => M=1/13 or M>=2/25   [q<=2max(A)<=36 < 38]")
def Bc(c,q,val): return frozenset(k for k in range(q) if min((c*k)%q,(q-(c*k)%q))<=val)
print("\ndiscrete covering test (can <=12 sets B_c cover Z/q?):")
for v,q in pairs[:8]:
    cands=sorted({Bc(c,q,v) for c in range(v,q-v+1)},key=lambda s:-len(s)); U=frozenset(range(q)); cov=set(); n=0
    for _ in range(15):
        b=max(cands,key=lambda s:len(s-cov))
        if not (b-cov): break
        cov|=b; n+=1
        if cov==U: break
    print("  (val,q)=(%2d,%3d) greedy cover=%2d %s"%(v,q,n,"<=12 => NOT ruled out" if cov==U and n<=12 else "-> impossible"))
print("\nsanity:", [(nm,str(Mexact(A)[0])) for nm,A in
      [("AP{1..12}",list(range(1,13))),("{1..11,24}",list(range(1,12))+[24])]])
