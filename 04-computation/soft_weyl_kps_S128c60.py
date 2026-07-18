#!/usr/bin/env python3
"""soft_weyl_kps_S128c60.py -- kind-pasteur S128 cont.60.
SOFT WEYL / FOURIER VIEW of mu(V) = meas{t : ||v t|| >= 1/14 for all v}.
Expand 1_{||x||>=lam} = sum_m c_m e(mx),  c_0 = 1-2lam,  c_m = -sin(2 pi m lam)/(pi m).
At lam = 1/14:  c_m = -sin(pi m / 7)/(pi m)  -- VANISHES exactly when 7 | m.
Orthogonality gives   mu = sum over m with m.v = 0 of prod_j c_{m_j},
so   mu = (1-2lam)^13 + (relation corrections).
MAIN TERM (6/7)^13 = 0.13537 > 0: the whole difficulty is the correction.
Measure the correction on real families -- including the TIGHT one, where mu = 0 exactly and
the correction must cancel the main term exactly.  PRINT DATA ONLY."""
import sys, math, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
LAM=F(1,14)
def c(m):
    if m==0: return 1-2*float(LAM)
    return -math.sin(2*math.pi*m*float(LAM))/(math.pi*m)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def mu_exact(V):
    """exact measure of the safe set, as a Fraction"""
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
        if all(nd(p*((a+b)/2))>=LAM for p in V): tot+=b-a
    return tot
print("### the coefficients at lam = 1/14 ###")
print("  m    c_m          |c_m|*(7/6)   note")
for m in range(0,16):
    note=""
    if m and m%7==0: note="<-- VANISHES (7 | m)"
    print("  %-4d %-12.8f %-13.8f %s"%(m,c(m),abs(c(m))*7/6,note))
print()
print("  MAIN TERM (6/7)^13 = %.8f"%((6/7)**13))
print("  the safe-indicator Fourier support at lam=1/14 avoids ALL multiples of 7.")
print()
print("### correction = mu - (6/7)^13, measured on real families ###")
print("  family                                   mu exact        mu float     correction     corr/main")
main=(6/7)**13
FAMS=[("tight {1..13}",list(range(1,14))),
      ("{1..12} u {169}",list(range(1,13))+[169]),
      ("{2..12} u {169,182}",list(range(2,13))+[169,182]),
      ("{2..12} u {169,338}",list(range(2,13))+[169,338]),
      ("{1..11} u {312,364}",list(range(1,12))+[312,364]),
      ("{2..12} u {1000,2000}",list(range(2,13))+[1000,2000]),
      ("{1,3,4,5,7,11,18,32}+pad",[1,3,4,5,7,11,18,32,101,103,107,109,113]),
      ]
for name,V in FAMS:
    if len(set(V))!=13: 
        print("  %-40s (not 13 distinct, skipped)"%name); continue
    m=mu_exact(V); mf=float(m)
    print("  %-40s %-15s %-12.8f %-14.8f %.4f"%(name,("%d/%d"%(m.numerator,m.denominator))[:15],mf,mf-main,(mf-main)/main))
print()
print("### the relation lattice: which relations carry the correction? ###")
print("  support-2 (ratio) and support-3 (additive) relations, |m_j| <= 6, weight = prod |c|*(7/6)")
def relweight(V,smax=3,Mmax=6):
    n=len(V); tot=0.0; rows=[]
    for s in range(2,smax+1):
        for idx in itertools.combinations(range(n),s):
            for mm in itertools.product([x for x in range(-Mmax,Mmax+1) if x!=0],repeat=s):
                if mm[0]<=0: continue          # take one of each +-m pair
                if sum(mm[t]*V[idx[t]] for t in range(s))!=0: continue
                w=1.0
                for t in range(s): w*=abs(c(mm[t]))*7/6
                tot+=2*w
                rows.append((2*w,[ (V[idx[t]],mm[t]) for t in range(s)]))
    rows.sort(reverse=True)
    return tot,rows
for name,V in FAMS[:5]:
    if len(set(V))!=13: continue
    tot,rows=relweight(V)
    print("  %-40s sum of |weights| (supp<=3) = %.5f   %s"%(
        name,tot,"< 1 OK" if tot<1 else ">= 1 (absolute bound fails)"))
    print("        heaviest:",[(round(w,5),r) for w,r in rows[:3]])
print("DONE")
