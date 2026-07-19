#!/usr/bin/env python3
"""whole_safe_set_kps_S128c73.py -- kind-pasteur S128 cont.73 (part 2).
THE FIX: use the WHOLE safe set, not one component.
The bad-j run is CONTIGUOUS, of width ~0.046*k1, i.e. it occupies an interval of length
~0.046 in t.  But S(P) has total measure 0.16-0.36 spread across [0,1].  So S(P) cannot be
swallowed by the bad zone, and the theorem only needs ONE good k1-gap anywhere in S(P).
TEST: for every 8-speed core and consecutive-type quadruple, is there a k1-gap that lies
inside S(P) AND outside the bad zone (i.e. a gap whose surviving piece beats 1/(7k4))?
PRINT DATA ONLY."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def safe_set(P):
    bps={F(0),F(1)}
    for p in P:
        for j in range(p+1):
            for s in (F(1,14*p),-F(1,14*p)):
                v=F(j,p)+s
                if 0<=v<=1: bps.add(v)
    B=sorted(bps); out=[]
    for i in range(len(B)-1):
        a,b=B[i],B[i+1]
        if b<=a: continue
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P): out.append((a,b))
    mg=[]
    for a,b in out:
        if mg and mg[-1][1]==a: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    return mg
def best_over_safe(P,ks):
    """max over k1-gaps lying inside S(P) of the longest surviving piece"""
    k1,k2,k3,k4=ks
    iv=safe_set(P); best=F(0); n=0
    for (A0,B0) in iv:
        j=int(A0*k1)-1
        while F(j,k1)-F(1,14*k1) < B0:
            A=F(j,k1)+F(1,14*k1); B=F(j+1,k1)-F(1,14*k1)
            if A>=A0 and B<=B0 and B>A:
                n+=1
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
                if L>best: best=L
            j+=1
    return best,n
C8=[sorted(c) for c in itertools.combinations(range(1,13),8)]
random.seed(73)
print("### using the WHOLE safe set: max over k1-gaps inside S(P) ###")
print("  killers                   cores   min over cores of 7*k4*best   #gaps(min core)   worst core")
allmin=None
for ks in [(157,158,159,160),(197,198,199,200),(317,318,319,320),(394,395,396,397),
           (300,301,302,303),(371,374,377,379)]:
    mn=None
    for P in C8:
        M=max(P)
        if ks[0]<=13*M: continue
        b,n=best_over_safe(P,ks)
        v=float(7*ks[3]*b)
        if mn is None or v<mn[0]: mn=(v,tuple(P),n)
    if mn is None: continue
    print("  %-25s %-7d %-29.5f %-17d %s"%(str(ks),len(C8),mn[0],mn[2],list(mn[1])))
    if allmin is None or mn[0]<allmin[0]: allmin=(mn[0],ks,mn[1])
print()
print("  MIN over all tested (core, quadruple) pairs: 7*k4*best = %.5f"%allmin[0])
print("    at killers %s, core %s"%(str(allmin[1]),list(allmin[2])))
print("  theorem needs > 1 : %s"%("HOLDS" if allmin[0]>1 else "FAILS"))
print()
print("### broader sweep over consecutive quadruples ###")
gm=None
for k1 in range(157,420,17):
    ks=(k1,k1+1,k1+2,k1+3)
    for P in random.sample(C8,40):
        M=max(P)
        if k1<=13*M: continue
        b,n=best_over_safe(P,ks)
        v=float(7*ks[3]*b)
        if gm is None or v<gm[0]: gm=(v,ks,tuple(P))
print("  min 7*k4*best over the sweep: %.5f at killers %s core %s"%(gm[0],str(gm[1]),list(gm[2])))
print("  VERDICT: %s"%("whole-safe-set conditioning HOLDS on this sweep" if gm[0]>1 else "FAILS"))
print("DONE")
