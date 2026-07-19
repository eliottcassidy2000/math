#!/usr/bin/env python3
"""hit_vs_sojourn_kps_S128c80.py -- kind-pasteur S128 cont.80 (part 2).
MY 4|e CRITERION IS REFUTED.  It correctly predicts which geodesics HIT the centre, but
hitting is NOT sufficient for maximal sojourn.
(1,2,7) has e = 4, and indeed at u = 3/4:  frac(-3/4)=1/4, frac(-3/2)=1/2, frac(-21/4)=3/4
-- it passes EXACTLY through c = (1/4,1/2,3/4).  Yet its sojourn is only ~7.3x |B|, not 28x.
WHY: sojourn ~ (number of hits per unit u) * (chord time), and the chord time scales as
2 rho / d_max.  So among directions that hit, the SLOWEST wins -- and (1,2,3) is the minimal
increasing triple, with (m,2m,3m) tied because m hits exactly compensate d_max = 3m."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def bad(g):
    cuts=[]
    for x in g:
        h=(7.0/6.0)*x
        a=max(h-1/6,0.0); b=min(h,1.0)
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=0.0; L=0.0
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        if b>cur: cur=b
    if 1.0-cur>L: L=1.0-cur
    return L<=1/6+1e-12
def sojourn(d,M=200000):
    return sum(1 for t in range(M) if bad(tuple((-di*t/M)%1.0 for di in d)))/M
print("### (1) does (1,2,7) really hit the centre? exact check at u = 3/4 ###")
u=F(3,4)
for d in [1,2,7]:
    print("   frac(-%d * 3/4) = %s"%(d,(-d*u)%1))
print("   target c = (1/4, 1/2, 3/4)  -> HIT: %s"%(tuple((-d*u)%1 for d in [1,2,7])==(F(1,4),F(1,2),F(3,4))))
print()
print("### (2) hitting does NOT give 28x -- sojourn vs d_max among CENTRE-HITTING directions ###")
print("  direction     hits c   d_max   sojourn     x|B|      sojourn*d_max")
B=0.003367
for d in [(1,2,3),(2,4,6),(3,6,9),(1,2,7),(1,6,7),(1,2,11),(1,6,11),(1,10,15),(2,4,14)]:
    u=F(3,4)
    hit = tuple((-di*u)%1 for di in d)==(F(1,4),F(1,2),F(3,4))
    if not hit:
        for t in range(1,400):
            uu=F(t,400)
            if tuple((-di*uu)%1 for di in d)==(F(1,4),F(1,2),F(3,4)): hit=True; break
    s=sojourn(d)
    print("  %-13s %-8s %-7d %-11.6f %-9.2f %.4f"%(str(d),"YES" if hit else "no",max(d),s,s/B,s*max(d)))
print()
print("### (3) the product sojourn * d_max is the invariant ###")
print("  for the proportional family (m,2m,3m) it should be constant:")
for m in [1,2,3,5,10]:
    d=(m,2*m,3*m); s=sojourn(d)
    print("    m=%-3d d_max=%-4d sojourn=%.6f   sojourn*d_max = %.4f"%(m,3*m,s,s*3*m))
print()
print("  CONCLUSION: 4|e predicts HITTING, not maximal sojourn.  Sojourn ~ hits * 2rho/d_max,")
print("  so among hitting directions the slowest wins; (1,2,3) is the minimal increasing")
print("  triple that hits, and (m,2m,3m) ties it because m hits compensate d_max = 3m exactly.")
print("DONE")
