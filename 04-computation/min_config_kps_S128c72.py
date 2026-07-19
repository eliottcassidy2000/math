#!/usr/bin/env python3
"""min_config_kps_S128c72.py -- kind-pasteur S128 cont.72.
FINDING THE MINIMISING CONFIGURATION.

CONSERVATIVE REDUCTION.  The four-comb theorem needs SOME surviving piece inside the
core-safe component to exceed 1/(7 k4).  Rather than track which gaps the core makes
available, minimise over ALL k1-gaps in [0,1]:

    W(k1,k2,k3,k4) := min over k1-gaps g of ( longest piece of g minus the k2,k3,k4 teeth ).

If  7*k4*W > 1  then the theorem holds for that quadruple no matter which gap the core
leaves, provided the component contains at least one full k1-gap -- which it does, since the
495-core atlas gives component length >= 1/70 and every legal k1 exceeds 13*max(P) >= 104.
So W is a genuinely conservative certificate, and the question is its MINIMUM over quadruples.
PRINT DATA ONLY."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def worst_gap(k1,k2,k3,k4):
    """min over k1-gaps of the longest surviving piece; exact rationals"""
    best=None
    for j in range(k1):
        A=F(j,k1)+F(1,14*k1); B=F(j+1,k1)-F(1,14*k1)
        cuts=[]
        for k in (k2,k3,k4):
            i=int(A*k)-1
            while F(i,k)-F(1,14*k) < B:
                x=F(i,k)-F(1,14*k); y=F(i,k)+F(1,14*k)
                if y>A and x<B: cuts.append((max(x,A),min(y,B)))
                i+=1
        cuts.sort(); cur=A; L=F(0)
        for x,y in cuts:
            if x>cur and x-cur>L: L=x-cur
            cur=max(cur,y)
        if B>cur and B-cur>L: L=B-cur
        if best is None or L<best[0]: best=(L,j)
    return best
print("### W = min over k1-gaps of longest surviving piece ; need 7*k4*W > 1 ###")
print("  killers                     7*k4*W      min-gap index")
rows=[]
CASES=[(157,158,159,160),(157,159,161,163),(371,374,377,379),(550,553,554,558),
       (157,158,160,163),(200,201,203,206),(300,301,302,303),(157,170,183,196),
       (163,167,173,179),(157,158,159,161),(211,212,214,215),(157,163,169,175)]
for ks in CASES:
    W,j=worst_gap(*ks)
    v=float(7*ks[3]*W)
    rows.append((v,ks,j))
    print("  %-27s %-11.5f %d"%(str(ks),v,j))
rows.sort()
print()
print("  minimum over these: %.5f at %s"%(rows[0][0],str(rows[0][1])))
print()
print("### systematic search: consecutive-type quadruples over a k1 range ###")
print("  looking for the smallest 7*k4*W")
glob=None
for k1 in range(157,340):
    for d2 in range(1,7):
        for d3 in range(d2+1,d2+7):
            for d4 in range(d3+1,d3+7):
                ks=(k1,k1+d2,k1+d3,k1+d4)
                W,j=worst_gap(*ks)
                v=float(7*ks[3]*W)
                if glob is None or v<glob[0]: glob=(v,ks,j)
    if k1%40==0: print("  ... k1=%d  running min %.5f at %s"%(k1,glob[0],str(glob[1])))
print()
print("  GLOBAL MIN over the searched family: 7*k4*W = %.5f"%glob[0])
print("    at killers %s (worst gap index %d)"%(str(glob[1]),glob[2]))
print("  theorem needs > 1 ; the three-tooth ratio target was 1.295 x equal-split")
W,_=worst_gap(*glob[1])
k1,k2,k3,k4=glob[1]
eq=(F(6,7*k1)-F(1,7*k2)-F(1,7*k3)-F(1,7*k4))/4
print("  at that configuration: W = %s ; equal-split = %s ; ratio W/eq = %.5f"%(W,eq,float(W/eq)))
print("DONE")
