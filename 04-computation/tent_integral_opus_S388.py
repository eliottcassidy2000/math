# opus-2026-07-17-S388 -- HYP-7740: THE TENT-INTEGRAL CRITERION.
#
# g(t) = min_v ||v t|| vanishes at every k/v and rises to a local maximum in each
# gap between consecutive zeros.  On a cell of length L the graph lies under the
# tent of height max, so area <= max * L / 2.  Therefore if EVERY local maximum
# were < 1/14,
#     int_0^1 g  <  (1/14) * (sum of L_i)/2  =  (1/14)*(1/2)  =  1/28.
# CONTRAPOSITIVE:   int g >= 1/28   =>   some local max >= 1/14   =>  LRC(14).
#
# A sufficient condition, pointwise in origin but computable as a single integral.
# THE HEURISTIC THAT MAKES IT INTERESTING: if the 13 values ||v t|| behaved like
# independent uniforms on [0,1/2], then E[min] = (1/2)/(13+1) = 1/28 EXACTLY --
# the criterion sits precisely on the boundary of the independence model, which
# is another appearance of the knife-edge this problem keeps producing.
from fractions import Fraction as F
import random
LAM=F(1,14)

def g_integral_exact(V):
    """exact int_0^1 min_v ||v t|| dt, by integrating over the arrangement."""
    # breakpoints: zeros k/v and peaks (2k+1)/(2v) of each wave
    bps=set([F(0),F(1)])
    for v in V:
        for k in range(v+1): bps.add(F(k,v))
        for k in range(v):   bps.add(F(2*k+1,2*v))
    # crossings of two waves: k/(v_i -+ v_j)
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            for q in (V[i]+V[j], abs(V[i]-V[j])):
                if q==0: continue
                for k in range(q+1): bps.add(F(k,q))
    pts=sorted(p for p in bps if 0<=p<=1)
    tot=F(0)
    def gval(t):
        best=None
        for v in V:
            r=(v*t)%1
            d=min(r,1-r)
            if best is None or d<best: best=d
        return best
    for a,b in zip(pts, pts[1:]):
        if b<=a: continue
        # g is linear on each cell: integrate via endpoints (trapezoid is EXACT)
        tot += (gval(a)+gval(b))*(b-a)/2
    return tot

THRESH=F(1,28)
print(f"    threshold 1/28 = {float(THRESH):.8f}")
print()
print("(1) THE TIGHT FAMILIES -- criterion must FAIL there (max is exactly 1/14)")
for nm,V in [("{1,...,13}",list(range(1,14))),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
    I=g_integral_exact(V)
    print(f"    {nm:16s} int g = {float(I):.8f}   {'CERTIFIES' if I>=THRESH else 'below threshold'}"
          f"   ratio to 1/28 = {float(I/THRESH):.4f}")

print()
print("(2) RANDOM FAMILIES -- how often does the criterion fire?")
random.seed(388)
fires=0; n=0; ratios=[]
for _ in range(12):
    V=sorted(random.sample(range(1,40),13))
    I=g_integral_exact(V); n+=1
    ratios.append(float(I/THRESH))
    if I>=THRESH: fires+=1
ratios.sort()
print(f"    {n} families: criterion fires {fires}/{n}")
print(f"    int g / (1/28): min {ratios[0]:.4f}, median {ratios[len(ratios)//2]:.4f}, max {ratios[-1]:.4f}")
print()
print("(3) THE INDEPENDENCE HEURISTIC: E[min of 13 uniforms on [0,1/2]] = (1/2)/14 = 1/28")
print("    so the criterion sits EXACTLY on the boundary of the independence model --")
print("    the same knife-edge that makes S1 = 13/7 and the k=7 ceiling appear.")
