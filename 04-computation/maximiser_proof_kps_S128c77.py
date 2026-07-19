#!/usr/bin/env python3
"""maximiser_proof_kps_S128c77.py -- kind-pasteur S128 cont.77.
WHY d PROPORTIONAL TO (1,2,3) IS THE MAXIMISER -- the structural reason.

In the continuum model, tooth i occupies [h_i - 1/6, h_i] cap [0,1] with
        h_i = (7/6) * frac(-d_i u) .
BAD means all four surviving pieces are <= 1/6.  The teeth have total width 1/2, so the
four pieces total 1/2, and 'all <= 1/6' with total 1/2 forces them close to 1/8 each.

THE EXACTLY-BALANCED CONFIGURATION.  Pieces 1/8, tooth, 1/8, tooth, 1/8, tooth, 1/8 puts the
tooth RIGHT EDGES at
        1/8 + 1/6 = 7/24 ,   7/24 + 1/8 + 1/6 = 7/12 ,   7/12 + 1/8 + 1/6 = 7/8 ,
i.e. h = (7/24, 7/12, 7/8) -- in EXACT RATIO 1 : 2 : 3.
Since h_i = (7/6) frac(-d_i u), the balanced configuration needs frac(-d_i u) in ratio
1:2:3.  Sustaining that over an INTERVAL of u (which positive bad measure requires) forces
d proportional to (1,2,3); for other d the ratio can only hold at isolated u.
VERIFY the balanced-configuration arithmetic and the ratio structure of real bad sets."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def teeth(u,DS):
    out=[]
    for d in DS:
        h=F(7,6)*((-d*u)%1)
        a=max(h-F(1,6),F(0)); b=min(h,F(1))
        if b>a: out.append((a,b))
    return sorted(out)
def Finf(u,DS):
    cur=F(0); L=F(0)
    for a,b in teeth(u,DS):
        if a>cur and a-cur>L: L=a-cur
        if b>cur: cur=b
    if 1-cur>L: L=1-cur
    return L
print("### (1) the exactly-balanced configuration ###")
h=[F(7,24),F(7,12),F(7,8)]
print("  tooth right edges h = %s"%[str(x) for x in h])
print("  ratios h1:h2:h3 = 1 : %s : %s"%(h[1]/h[0],h[2]/h[0]))
pieces=[h[0]-F(1,6), h[1]-F(1,6)-h[0], h[2]-F(1,6)-h[1], 1-h[2]]
print("  resulting pieces = %s ; all equal 1/8: %s ; sum = %s"%(
    [str(p) for p in pieces],all(p==F(1,8) for p in pieces),sum(pieces)))
print("  longest piece 1/8 <= threshold 1/6 : %s  (so this configuration is BAD)"%(F(1,8)<=F(1,6)))
print()
print("### (2) does d ~ (1,2,3) realise it exactly, and others not? ###")
print("  d-triple      u realising h ratio 1:2:3        F(u)      bad?")
for DS in [(1,2,3),(2,4,6),(3,6,9),(1,2,4),(1,3,5),(2,3,4),(1,2,5)]:
    # want frac(-d_i u) = (1/4, 1/2, 3/4)*(6/7)*... solve for the (1,2,3) case: frac(-u)=1/4
    u=F(3,4)   # frac(-3/4) = 1/4
    hs=[F(7,6)*((-d*u)%1) for d in DS]
    r=[float(x/hs[0]) if hs[0] else 0 for x in hs]
    print("  %-13s h=%-28s ratio %-18s F=%-9s %s"%(
        str(DS),str([str(x) for x in hs]),str([round(x,3) for x in r]),str(Finf(u,DS)),Finf(u,DS)<=F(1,6)))
print()
print("### (3) the bad set of the balanced family, exactly ###")
for m in [1,2,3]:
    DS=(m,2*m,3*m)
    # bad interval from THM-1147: u in [5/(21m) + j/m, 2/(7m) + j/m] for j=0..m-1, plus mirrors
    lo=F(5,21*m); hi=F(2,7*m)
    print("  m=%d : predicted first run [%s, %s], width %s"%(m,lo,hi,hi-lo))
    print("        F at the two ends: %s and %s (threshold 1/6)"%(Finf(lo,DS),Finf(hi,DS)))
    print("        F just inside: %s"%Finf((lo+hi)/2,DS))
print()
print("### (4) for non-proportional d, can the 1:2:3 ratio persist on an interval? ###")
print("  d-triple      #u in a 2520-grid with h-ratio within 1%% of (1,2,3)")
N=2520
for DS in [(1,2,3),(2,4,6),(1,2,4),(1,3,5),(2,3,4),(1,2,5),(3,5,7)]:
    c=0
    for t in range(N):
        u=F(t,N)
        hs=[F(7,6)*((-d*u)%1) for d in DS]
        if hs[0]==0: continue
        r2=float(hs[1]/hs[0]); r3=float(hs[2]/hs[0])
        if abs(r2-2)<0.02 and abs(r3-3)<0.03: c+=1
    print("  %-13s %d / %d = %.4f"%(str(DS),c,N,c/N))
print("DONE")
