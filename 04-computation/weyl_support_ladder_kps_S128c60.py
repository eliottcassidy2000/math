#!/usr/bin/env python3
"""weyl_support_ladder_kps_S128c60.py -- kind-pasteur S128 cont.60.
THE SIGN LAW and the SUPPORT LADDER.
c_m = -sin(pi m/7)/(pi m):  NEGATIVE for m mod 14 in {1..6}, ZERO for 7|m, POSITIVE for
m mod 14 in {8..13}.  So a relation whose coefficients all have |m_j| <= 6 contributes with
sign (-1)^s, s = support.  That makes the Fourier expansion an ALTERNATING LADDER IN
SUPPORT -- the same alternation as the Bonferroni ladder of THM-930/935, but indexed by
RELATIONS instead of events.  Test the sign law, compute the ladder terms w_s, and see how
much of the exact correction lives at low support.  PRINT DATA ONLY."""
import sys, math, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
LAM=1/14.0
def c(m):
    if m==0: return 1-2*LAM
    return -math.sin(2*math.pi*m*LAM)/(math.pi*m)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def mu_exact(V):
    bps={F(0),F(1)}
    for p in V:
        for j in range(p+1):
            for s in (F(1,14*p),-F(1,14*p)):
                x=F(j,p)+s
                if 0<=x<=1: bps.add(x)
    B=sorted(bps); tot=F(0)
    for i in range(len(B)-1):
        a,b=B[i],B[i+1]
        if b<=a: continue
        if all(nd(p*((a+b)/2))>=F(1,14) for p in V): tot+=b-a
    return tot
print("### (1) the sign law ###")
ok=True
for m in range(1,60):
    if m%7==0:
        if abs(c(m))>1e-15: ok=False
        continue
    r=m%14
    want = -1 if 1<=r<=6 else 1
    if (c(m)<0)!=(want<0): ok=False; print("   violation at m=%d: c=%g want sign %d"%(m,c(m),want))
print("  c_m sign = -1 for m mod 14 in {1..6}, +1 for {8..13}, 0 for 7|m, checked m<60: %s"%ok)
print("  => a relation with all |m_j| <= 6 contributes with sign (-1)^support")
print()
def ladder(V,smax=4,Mmax=6,Mmax5=3):
    """signed ladder terms w_s (normalised by the main term), relations with m.V=0"""
    n=len(V); main=(1-2*LAM)**n
    W={}
    for s in range(2,smax+1):
        M = Mmax if s<=4 else Mmax5
        acc=0.0
        for idx in itertools.combinations(range(n),s):
            vs=[V[i] for i in idx]
            for mm in itertools.product([x for x in range(-M,M+1) if x!=0],repeat=s):
                if mm[0]<=0: continue
                if sum(mm[t]*vs[t] for t in range(s))!=0: continue
                pr=1.0
                for t in range(s): pr*=c(mm[t])
                acc+=2*pr/((1-2*LAM)**s)
        W[s]=acc
    return W,main
FAMS=[("tight {1..13}",list(range(1,14))),
      ("{1..12} u {169}",list(range(1,13))+[169]),
      ("{2..12} u {169,182}",list(range(2,13))+[169,182]),
      ("{2..12} u {1000,2000}",list(range(2,13))+[1000,2000]),
      ("{1..11} u {312,364}",list(range(1,12))+[312,364]),
      ]
print("### (2) the support ladder:  mu/main = 1 + w2 + w3 + w4 + ...  (signed) ###")
print("  family                          mu/main    1+w2      +w3       +w4       low-supp captures")
for name,V in FAMS:
    if len(set(V))!=13: continue
    m=float(mu_exact(V)); W,main=ladder(V)
    r=m/main
    p2=1+W[2]; p3=p2+W[3]; p4=p3+W[4]
    print("  %-31s %-10.5f %-9.5f %-9.5f %-9.5f %.1f%%"%(
        name,r,p2,p3,p4,100*(p4-1)/(r-1) if abs(r-1)>1e-12 else float('nan')))
print()
print("### (3) the individual ladder terms (signs should alternate: w2>0, w3<0, w4>0) ###")
print("  family                          w2         w3         w4")
for name,V in FAMS:
    if len(set(V))!=13: continue
    W,main=ladder(V)
    print("  %-31s %-10.5f %-10.5f %-10.5f"%(name,W[2],W[3],W[4]))
print()
print("### (4) does the alternating truncation BOUND mu?  (1+w2+w3 <= mu/main <= 1+w2 ?) ###")
for name,V in FAMS:
    if len(set(V))!=13: continue
    m=float(mu_exact(V)); W,main=ladder(V); r=m/main
    lo=1+W[2]+W[3]; hi=1+W[2]
    print("  %-31s lower %-9.5f  mu/main %-9.5f  upper %-9.5f   brackets: %s"%(
        name,lo,r,hi,lo<=r<=hi))
print("DONE")
