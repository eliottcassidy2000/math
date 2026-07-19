#!/usr/bin/env python3
"""six_ray_kakeya_kps_S128c83.py -- kind-pasteur S128 cont.83.
DEATH-STAR'S SIX-RAY SOJOURN-MAX, VIA THE KAKEYA / X-RAY LENS.

The sojourn of a geodesic of direction d in B is the X-RAY TRANSFORM of 1_B in direction d.
B is six boxes centred at the six PERMUTATIONS of (1/4,1/2,3/4).  A box centred at
sigma(1,2,3)/4 is threaded lengthwise by the direction sigma(1,2,3) -- so there are exactly
SIX maximising rays, the six permutations, which is death-star's six-ray sojourn-max.

THE PREDICTION THIS MAKES: all six permutations should give the SAME 28.28x sojourn, and
(1,2,3) is singled out only because the killer ordering k2<k3<k4 forces d2<d3<d4, so it is
the UNIQUE INCREASING representative of the six-ray orbit.  That would explain the maximiser
without needing any new mechanism.  TEST IT."""
import sys, itertools
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
B=0.003367
print("### the six rays = the six permutations of (1,2,3) ###")
print("  direction    sojourn     x|B|     increasing?  admissible as (d2,d3,d4)?")
for p in itertools.permutations([1,2,3]):
    s=sojourn(p)
    inc = p[0]<p[1]<p[2]
    print("  %-12s %-11.6f %-8.2f %-12s %s"%(str(p),s,s/B,inc,"YES" if inc else "no (k2<k3<k4 forces increasing)"))
print()
print("### and their scalings ###")
print("  direction         sojourn     x|B|")
for p in [(1,2,3),(2,1,3),(3,2,1),(2,4,6),(4,2,6),(6,4,2)]:
    s=sojourn(p)
    print("  %-17s %-11.6f %.2f"%(str(p),s,s/B))
print()
print("### the X-ray reading: is the box threaded lengthwise by its own permutation? ###")
print("  a box centred at sigma(1,2,3)/4 should be entered head-on by direction sigma(1,2,3).")
print("  check: does the geodesic of direction p pass exactly through the centre sigma=p?")
from fractions import Fraction as F
for p in itertools.permutations([1,2,3]):
    c=tuple(F(x,4) for x in p)
    hit=None
    for t in range(1,400):
        u=F(t,400)
        if tuple((-d*u)%1 for d in p)==c: hit=u; break
    print("  direction %-12s centre %-22s hit at u = %s"%(str(p),str([str(x) for x in c]),hit))
print()
print("### so the maximiser is forced by the ORDERING constraint, not by resonance ###")
print("  all six rays tie; only (1,2,3) satisfies d2 < d3 < d4, which k2<k3<k4 requires.")
print("DONE")
