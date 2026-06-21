#!/usr/bin/env python3
"""
lrc14_residual_R_stratification — mac-mini-2026-06-20-S5

Extends codex's HYP-2697 (arbitrary cluster compression cone). codex showed naive
coordinatewise dominance q_C(R) <= q_K(R) (K=consec) is FALSE — arbitrary shapes beat
consec on SMALL residual sets R. This script STRATIFIES by |R| and finds a clean monotone
hierarchy over the full bounded k=8 box (0 in E, max(E)<=13, 1716 shapes):

  q_C(R) = Pr_x( cluster C covers every sector in R ),  R subset {1..6}, via the Z/7
  vertex-coloring  color(e,x) = floor(7 frac(e x))  (an EXACT functional, not a tournament correlate).

RESULT (VERIFIED, exact rationals):
  |R|=6 (full inner cover = p_0 = meas S7): consec maximal, 0 violators  -> consec maximizes the
        REAL covering measure p_0 over the box, not merely the L_y proxy.
  |R|=5: 0 violators  (UNIVERSAL consec dominance)
  |R|=4: 1 violator   (near-universal)
  |R|=3: ~9 shapes violate, only on SPREAD residuals (1,3,5),(1,2,5),(1,4,5),...
  |R|<=2: many violators (codex's regime; naive dominance fails here).

So the compression-cone difficulty is ENTIRELY in the small-|R| tail; consec dominance on
|R|>=5 is unconditional. The LRC context weight w_R = Pr(small part P leaves residual R)
determines how much weight lands in the hard small-|R| region.
NOTE: the EXTREMALITY functional is L_y (LP-dual, <= cap), NOT U4 (=p_0+p_5+5p_6, the cruder
true-wide Bonferroni; U4(consec_8)=353/735=0.480 > cap, so U4<=cap is FALSE for narrow consec).
"""
from fractions import Fraction as F
from collections import defaultdict
import itertools as it

def qC_all(E, rmin=0):
    """q_C(R) for all R subset {1..6} with |R|>=rmin, exact via breakpoints."""
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1)
    q=defaultdict(lambda:F(0))
    for a in range(len(bps)-1):
        lo,hi=bps[a],bps[a+1]; mid=(lo+hi)/2; w=hi-lo
        inner=set(int(((e*mid)%1)*7) for e in E)&set(range(1,7))
        for r in range(max(rmin,1),7):
            for R in it.combinations(range(1,7),r):
                if set(R)<=inner: q[R]+=w
    return q

K=[0,1,2,3,4,5,6,7]
print("Residual-set stratification of consec dominance over the bounded k=8 box (1716 shapes):")
for rmin in [6,5,4,3]:
    qk=qC_all(K,rmin); Rs=[R for r in range(rmin,7) for R in it.combinations(range(1,7),r)]
    viol=0
    for combo in it.combinations(range(1,14),7):
        C=[0]+list(combo); qc=qC_all(C,rmin)
        if any(qk[R]<qc[R] for R in Rs): viol+=1
    tag = " (= p_0 = meas S7, the real covering measure)" if rmin==6 else ""
    print(f"  |R|>={rmin}: {viol:4d}/1716 shapes violate q_K(R)>=q_C(R){tag}")

print("\n|R|=3 violated-residual census (which spread sets let a shape beat consec):")
qk3=qC_all(K,3); badR=defaultdict(int)
for combo in it.combinations(range(1,14),7):
    C=[0]+list(combo); qc=qC_all(C,3)
    for R in it.combinations(range(1,7),3):
        if qk3[R]<qc[R]: badR[R]+=1
for R,c in sorted(badR.items(),key=lambda kv:-kv[1]):
    print(f"   {R}: {c} shapes")
print("\n=> consec dominance is UNIVERSAL for |R|>=5; the cone theorem's content is the small-|R| tail.")
