#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S29 (HYP-4592) -- the TWO-SIDED WITNESS COMPETITION.

A gap-witness is a NARROW two-sided target.  The canonical S_N intends clearance 3 at
q=3N+2 (the mediant 3/(3N+2), always strictly inside the open gap).  It is a genuine
gap-witness iff that intended witness (i) HOLDS clearance exactly 3 AND (ii) DOMINATES
all competitors at q' in {v_i+-v_j} u {2v_i}.  Two failure walls bracket the target:
  OVERSHOOT (competitor beats from above; N=12: q'=35 clears 3 => 3/35 > 2/25 top),
  DEGRADE  (intended clearance < 3 => M collapses to floor; N=31: q=95 clears 2 => 1/32).
Generic 6k+1 threads the slab; N=12 hits the top wall, N=31 the bottom.  This is the
(G) proof-target: show the OVERSHOOT wall is forced for EVERY covering 12-family.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def clear_at(S,q):
    """best clearance c and its a, at fixed denominator q (max over a of min-gap)."""
    best=0; ba=0
    for a in range(1,q):
        mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
        if mn>best: best=mn; ba=a
    return best,ba

def Mfull(S):
    S=sorted(set(S)); Q=set([2*v for v in S])
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0); bq=None
    for q in Q:
        c,_=clear_at(S,q)
        if F(c,q)>best: best=F(c,q); bq=q
    return best,bq

if __name__=="__main__":
    print("WITNESS COMPETITION: S_N is a gap-witness iff intended q=3N+2 holds clearance 3 AND no competitor overshoots.")
    print("Two failure modes: OVERSHOOT (M>gap-top, e.g. N=12) vs DEGRADE (intended clearance<3 => M<=floor, e.g. N=31).\n")
    print(f"  {'N':>3} {'gap window':>20} {'M(S_N)':>8} {'q=3N+2 clr':>11} {'M vs gap':>18}")
    for N in [7,12,13,19,25,31,37,43]:
        S=sorted(set(list(range(1,N-1))+[N,3*(N-1)])); S=[x//reduce(gcd,S) for x in S]
        lo,hi=F(1,N+1),F(2,2*N+1); M,bq=Mfull(S)
        cint,_=clear_at(S,3*N+2)   # clearance the INTENDED mediant witness actually holds
        if M>hi: verd="OVERSHOOT (too big)"
        elif M<=lo: verd="DEGRADE (<=floor)"
        elif lo<M<hi: verd="** IN GAP **"
        else: verd="below"
        print(f"  {N:>3} ({float(lo):.4f},{float(hi):.4f}) {str(M):>8} {str(cint)+'/'+str(3*N+2):>11} {verd:>18}")
    print("\n  N=12: intended q=38 holds clr 3 BUT competitor q=35 (=2+33) also clr 3 => 3/35 OVERSHOOTS gap-top 2/25.")
    print("  N=31: intended q=95 clearance DROPS to 2 (AP elt collides in width-2 hole) => M falls to floor 1/32.")
    print("  => gap-witness = narrow two-sided target; (G)-target: OVERSHOOT wall (factorable q'<38, clr>=3) forced for ALL covering 12-families.")
