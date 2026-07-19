#!/usr/bin/env python3
"""standoff_test_kps_S128c79.py -- kind-pasteur S128 cont.79.
TESTING THE STANDOFF BOUND BEFORE TRYING TO PROVE IT.

THM-1150 conjectured: every non-proportional integer direction keeps standoff > rho ~ 0.041
from all six centres (permutations of (1/4,1/2,3/4)), so its geodesic MISSES B.
BUT a geodesic of direction d has length ~|d| in T^3, so for LARGE non-proportional d it
should equidistribute and come arbitrarily close to any point -- making standoff -> 0 and
the conjecture FALSE.  If so the maximiser claim must survive by a different route:
sojourn -> |B| = 0.0034 for generic d, against 2/21 = 0.0952 for d ~ (1,2,3), a factor 28.
Test standoff and sojourn as |d| grows.  PRINT DATA ONLY."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
CENTRES=[tuple(x/4 for x in p) for p in itertools.permutations([1,2,3])]
def bad_from_g(g):
    cuts=[]
    for x in g:
        h=(7.0/6.0)*x
        a=max(h-1.0/6.0,0.0); b=min(h,1.0)
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=0.0; L=0.0
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        if b>cur: cur=b
    if 1.0-cur>L: L=1.0-cur
    return L<=1.0/6.0+1e-12
def standoff(DS,M=200000):
    best=1.0
    for t in range(M):
        u=t/M
        g=tuple((-d*u)%1.0 for d in DS)
        for c in CENTRES:
            dd=max(min(abs(gi-ci),1.0-abs(gi-ci)) for gi,ci in zip(g,c))
            if dd<best: best=dd
    return best
def sojourn(DS,M=200000):
    c=0
    for t in range(M):
        u=t/M
        if bad_from_g(tuple((-d*u)%1.0 for d in DS)): c+=1
    return c/M
print("### standoff and sojourn as |d| grows, NON-proportional directions ###")
print("  direction              d_max   standoff    sojourn     vs |B|=0.003367")
random.seed(79)
rows=[]
for DS in [(1,2,4),(2,3,4),(3,5,7),(2,4,7),(5,9,14),(7,13,23),(11,19,37),
           (23,41,67),(41,73,131),(97,173,281),(211,367,593)]:
    so=standoff(DS); sj=sojourn(DS)
    rows.append((max(DS),so,sj,DS))
    print("  %-22s %-7d %-11.6f %-11.6f %.2f x"%(str(DS),max(DS),so,sj,sj/0.003367))
print()
print("### the proportional family for comparison ###")
for DS in [(1,2,3),(5,10,15),(20,40,60),(97,194,291)]:
    so=standoff(DS); sj=sojourn(DS)
    print("  %-22s %-7d %-11.6f %-11.6f %.2f x"%(str(DS),max(DS),so,sj,sj/0.003367))
print()
print("### verdict on the standoff conjecture ###")
mn=min(r[1] for r in rows)
print("  smallest standoff among the non-proportional directions tested: %.6f"%mn)
print("  THM-1150's rho = 0.0412")
print("  conjecture 'standoff > rho for all non-proportional d': %s"%(
    "HOLDS on this sample" if mn>0.0412 else "*** REFUTED -- standoff falls below rho"))
if mn<=0.0412:
    bad=[r for r in rows if r[1]<=0.0412]
    print("  witnesses:")
    for dmax,so,sj,DS in bad[:6]:
        print("    %-22s standoff %.6f  sojourn %.6f"%(str(DS),so,sj))
print()
print("  large-|d| sojourn should approach |B| = 0.003367 if equidistribution is the mechanism")
print("DONE")
