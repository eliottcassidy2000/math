#!/usr/bin/env python3
"""tube_geometry_kps_S128c78.py -- kind-pasteur S128 cont.78.
THE MAXIMISER PROOF, GEOMETRIC FORM.

Write g_i = frac(-d_i u) in [0,1); the tooth right-edge is h_i = (7/6) g_i and the
configuration is BAD iff the longest surviving piece is <= 1/6.  Badness depends ONLY on
(g2,g3,g4), so there is a FIXED region B in the 3-torus, independent of d:
        B = { g in T^3 : configuration is bad } .
The map u -> (frac(-d2 u), frac(-d3 u), frac(-d4 u)) is a CLOSED GEODESIC of direction
(d2,d3,d4), so
        bad measure(d) = sojourn time of that geodesic in B.
CLAIM: B is a TUBE around the balanced line L = {(s,2s,3s)}.  A geodesic of direction d
therefore runs ALONG the tube axis when d ~ (1,2,3) -- maximal sojourn -- and crosses it
transversally otherwise -- sojourn suppressed by the crossing angle.
TESTS: (a) exact volume of B; (b) is B concentrated near L; (c) the balanced point.
PRINT DATA ONLY."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def bad_from_g(g):
    """g = (g2,g3,g4) in [0,1)^3 ; returns True if longest surviving piece <= 1/6"""
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
    return L<=F(1,6), L
print("### (a) the volume of B ###")
N=140
cnt=0; tot=0
for i in range(N):
    for j in range(N):
        for k in range(N):
            g=(F(2*i+1,2*N),F(2*j+1,2*N),F(2*k+1,2*N))
            tot+=1
            if bad_from_g(g)[0]: cnt+=1
print("  grid %d^3 : |B| = %d/%d = %.6f"%(N,cnt,tot,cnt/tot))
print("  compare: measured bad measure of the (1,2,3) geodesic = 0.0952 = 2/21")
print("  ratio (sojourn)/(volume) = %.1f  -> the geodesic is FAR from equidistributed"%(0.095238/(cnt/tot)))
print()
print("### (b) is B a tube around the line L = {(s,2s,3s)}? ###")
def dist_to_L(g):
    """min over s of the torus distance from g to (s,2s,3s)"""
    best=None
    M=2000
    for t in range(M):
        s=t/M
        d=0.0
        for c,x in zip((1,2,3),g):
            e=abs((c*s)%1.0-float(x)); e=min(e,1.0-e)
            d=max(d,e)
        if best is None or d<best: best=d
    return best
random.seed(78)
inB=[]; outB=[]
while len(inB)<120 or len(outB)<120:
    g=(F(random.randrange(1,2000),2000),F(random.randrange(1,2000),2000),F(random.randrange(1,2000),2000))
    b,_=bad_from_g(g)
    if b and len(inB)<120: inB.append(dist_to_L(g))
    elif (not b) and len(outB)<120: outB.append(dist_to_L(g))
inB.sort(); outB.sort()
print("  sup-distance to L, points INSIDE B : min %.4f  median %.4f  max %.4f"%(inB[0],inB[60],inB[-1]))
print("  sup-distance to L, points OUTSIDE B: min %.4f  median %.4f  max %.4f"%(outB[0],outB[60],outB[-1]))
print("  every point of B within %.4f of L : %s"%(max(inB),max(inB)<0.2))
print()
print("### (c) the balanced point sits on L ###")
gb=(F(1,4),F(1,2),F(3,4))
b,L0=bad_from_g(gb)
print("  g = (1/4, 1/2, 3/4) = (s,2s,3s) at s=1/4 ; bad: %s ; longest piece = %s"%(b,L0))
print("  its tooth right edges h = %s"%[str(F(7,6)*x) for x in gb])
print()
print("### (d) sojourn vs direction: does alignment with (1,2,3) maximise? ###")
def sojourn(DS,M=4200):
    c=0
    for t in range(M):
        u=F(t,M)
        g=tuple((-d*u)%1 for d in DS)
        if bad_from_g(g)[0]: c+=1
    return c/M
print("  direction        sojourn      aligned with (1,2,3)?")
for DS in [(1,2,3),(2,4,6),(3,6,9),(1,2,4),(1,3,5),(2,3,4),(1,2,5),(3,5,7),(1,4,7)]:
    al = (DS[1]==2*DS[0] and DS[2]==3*DS[0])
    print("  %-16s %-12.6f %s"%(str(DS),sojourn(DS),"YES" if al else "no"))
print("DONE")
