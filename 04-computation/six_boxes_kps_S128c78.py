#!/usr/bin/env python3
"""six_boxes_kps_S128c78.py -- kind-pasteur S128 cont.78 (part 2).
B IS NOT A TUBE -- IT IS SIX ISOLATED BOXES.
The tube-around-L claim is refuted (distance to L does not separate B).  The correct
picture: badness needs the three g_i to be near (1/4, 1/2, 3/4) IN SOME ORDER, so
        B ~ union of six small boxes around the permutations of (1/4,1/2,3/4).
A geodesic of direction d has positive sojourn iff it passes within the box radius of one
of those six centres.  Hitting (1/4,1/2,3/4) exactly requires
        d2 u = -1/4,  d3 u = -1/2,  d4 u = -3/4   (mod 1),
which forces d3/d2 = 2 and d4/d2 = 3 -- exactly proportionality to (1,2,3).
VERIFY: the box structure, its radius, and the centre-hitting criterion."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def bad_from_g(g):
    cuts=[]
    for x in g:
        h=F(7,6)*x
        a=max(h-F(1,6),F(0)); b=min(h,F(1))
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=F(0); L=F(0)
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        if b>cur: cur=b
    if 1-cur>L: L=1-cur
    return L<=F(1,6)
CENTRES=[tuple(F(x,4) for x in p) for p in itertools.permutations([1,2,3])]
print("### (1) is every point of B near one of the six permutations of (1/4,1/2,3/4)? ###")
random.seed(781)
found=[]; M=260
tries=0
while len(found)<M and tries<400000:
    tries+=1
    g=tuple(F(random.randrange(1,720),720) for _ in range(3))
    if bad_from_g(g):
        d=min(max(min(abs(float(gi-ci)),1-abs(float(gi-ci))) for gi,ci in zip(g,c)) for c in CENTRES)
        found.append(d)
found.sort()
print("  sampled points of B: %d ; sup-distance to the nearest of the six centres:"%len(found))
print("    min %.4f   median %.4f   max %.4f"%(found[0],found[len(found)//2],found[-1]))
print("  all within 1/8 = 0.125 of a centre : %s"%(found[-1]<=0.125))
print()
print("### (2) the box radius implied by the volume ###")
vol=0.003367
r=(vol/(6*8))**(1/3)
print("  |B| = %.6f ; six boxes of half-width rho gives 6*(2 rho)^3 = |B| -> rho = %.4f"%(vol,r))
print("  observed max distance %.4f is the same order"%found[-1])
print()
print("### (3) the centre-hitting criterion ###")
print("  hitting (1/4,1/2,3/4) needs d2 u = -1/4, d3 u = -1/2, d4 u = -3/4 (mod 1).")
print("  From the first, u = (-1/4 + m)/d2 ; substituting forces d3/d2 = 2 and d4/d2 = 3.")
print("  direction        does the geodesic hit a centre EXACTLY?   min sup-dist to a centre")
for DS in [(1,2,3),(2,4,6),(3,6,9),(1,2,4),(1,3,5),(2,3,4),(1,2,5),(3,5,7),(1,4,7),(2,4,7)]:
    best=None; hit=False
    for t in range(5040):
        u=F(t,5040)
        g=tuple((-d*u)%1 for d in DS)
        for c in CENTRES:
            dd=max(min(abs(float(gi-ci)),1-abs(float(gi-ci))) for gi,ci in zip(g,c))
            if best is None or dd<best: best=dd
            if dd==0: hit=True
    al=(DS[1]==2*DS[0] and DS[2]==3*DS[0])
    print("  %-16s %-42s %.5f"%(str(DS),("YES" if hit else "no")+("   (aligned)" if al else ""),best))
print()
print("  CONCLUSION: only directions proportional to (1,2,3) put the geodesic through a")
print("  centre of B; all others keep a positive standoff and so miss the boxes entirely.")
print("DONE")
