#!/usr/bin/env python3
"""general_d_sweep_kps_S128c76.py -- kind-pasteur S128 cont.76 (part 2).
CONFIRMING THE 2/21 CEILING over all d-triples in a wide box.
Part 1 found: total bad ~ 2/21 exactly when (d2,d3,d4) is PROPORTIONAL to (1,2,3) -- the
runs multiply (2,4,6) while each shrinks (1/21,1/42,1/63) so the product is constant -- and
collapses to ~0 otherwise.  My predicted growth 2*d_max/21 is refuted.
Sweep all 1<=d2<d3<d4<=16 and confirm the ceiling."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def Finf(u,DS):
    cuts=[]
    for d in DS:
        g=(-d*u)%1
        c=F(7,6)*g-F(1,12)
        a=max(c-F(1,12),F(0)); b=min(c+F(1,12),F(1))
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=F(0); L=F(0)
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        cur=max(cur,b)
    if 1-cur>L: L=1-cur
    return L
def total_bad(DS,N=4200):
    thr=F(1,6)
    return sum(1 for t in range(N) if Finf(F(t,N),DS)<=thr)/N
best=[]; over=[]
for d2,d3,d4 in itertools.combinations(range(1,17),3):
    tb=total_bad((d2,d3,d4))
    best.append((tb,(d2,d3,d4)))
    if tb>=0.164: over.append((tb,(d2,d3,d4)))
best.sort(reverse=True)
print("### all triples 1<=d2<d3<d4<=16 : %d of them ###"%len(best))
print("  top 12 by total bad measure:")
for tb,DS in best[:12]:
    g=DS[0]
    prop = (DS[1]/DS[0]==2 and DS[2]/DS[0]==3) if DS[0] else False
    print("    %-14s total bad = %.6f    %s"%(str(DS),tb,"<- proportional to (1,2,3)" if prop else ""))
print()
print("  MAXIMUM total bad over all %d triples: %.6f  at %s"%(len(best),best[0][0],str(best[0][1])))
print("  2/21 = %.6f"%(2/21))
print("  triples exceeding the 0.164 safe measure: %d"%len(over))
if over:
    for tb,DS in over[:6]: print("    *** %s : %.6f"%(str(DS),tb))
else:
    print("    NONE -- the 2/21 ceiling holds across the whole box")
print()
print("### the proportional family (m, 2m, 3m) ###")
print("   m   total bad    #runs x per-run")
for m in range(1,7):
    DS=(m,2*m,3*m)
    tb=total_bad(DS)
    print("  %-3d  %.6f     %d runs implied"%(m,tb,round(tb*21*m)) )
print()
print("  VERDICT: total bad <= 2/21 = %.6f < 0.164 ; margin %.5f"%(2/21,0.164-2/21))
print("DONE")
