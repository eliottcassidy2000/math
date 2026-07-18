#!/usr/bin/env python3
"""weyl_dividing_line_kps_S128c61.py -- kind-pasteur S128 cont.61.
WHERE SOFT WEYL WORKS AND WHERE IT DIES -- reconciling four agents' concurrent results.
death-star THM-1037 PROVES the position lemma by soft Weyl; boxeph HYP-7505 closes density
off-resonance with it; boxeph S95 says it is the WRONG tool for one-line rigidity; my
THM-1061(V) shows the relation-support ladder DIVERGES.  Claim: these are consistent, and
the dividing line is whether ONE oscillating factor is tested against a FIXED set, or a
PRODUCT over all 13 speeds is expanded.
ONE SPEED: for G = union of C intervals, integration by parts gives
    |int_G e(k t) dt| = |sum_i [e(k b_i) - e(k a_i)]/(2 pi i k)| <= C/(pi |k|)   -- decays in k.
That IS my measure horn's boundary term and klein's 2dN/w, in Fourier clothing.
PRODUCT: expanding 13 factors gives the relation-support ladder, which diverges.
PRINT DATA ONLY."""
import sys, math, itertools, cmath
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
print("### (1) the one-speed soft Weyl bound  |int_G e(kt)| <= C/(pi k)  ###")
print("  core           C    mu        k       |int| exact    C/(pi k)     ratio")
for drop in [1,6,12]:
    P=[x for x in range(1,13) if x!=drop]
    iv=safe_set(P); C=len(iv); mu=float(sum(b-a for a,b in iv))
    for k in [157,500,2000,20000]:
        s=0j
        for a,b in iv:
            a=float(a); b=float(b)
            s+=(cmath.exp(2j*math.pi*k*b)-cmath.exp(2j*math.pi*k*a))/(2j*math.pi*k)
        bnd=C/(math.pi*k)
        print("  drop%-2d        %-4d %-9.5f %-7d %-13.8f %-12.8f %.4f"%(
            drop,C,mu,k,abs(s),bnd,abs(s)/bnd))
print()
print("### (2) that bound IS my measure horn's boundary term ###")
print("  geometric:  bad(k) within G  <=  mu/7 + 2C/(7k)")
print("  Fourier  :  |1_G hat(k)|     <=  C/(pi k)      -- same C/k shape, same origin (endpoints)")
print("  klein    :  loss             <=  2 d mu + 2 d N / w,  d = 1/14   -- same again")
print("  core           C    mu        k       true bad     mu/7+2C/(7k)   slack")
for drop in [1,12]:
    P=[x for x in range(1,13) if x!=drop]
    iv=safe_set(P); C=len(iv); mu=float(sum(b-a for a,b in iv))
    for k in [157,900,20000]:
        bad=0.0
        for a,b in iv:
            a=float(a); b=float(b)
            jlo=int(a*k)-1; jhi=int(b*k)+1; w=1.0/(14.0*k)
            for j in range(jlo,jhi+1):
                x=max(j/k-w,a); y=min(j/k+w,b)
                if y>x: bad+=y-x
        bnd=mu/7+2*C/(7*k)
        print("  drop%-2d        %-4d %-9.5f %-7d %-12.8f %-14.8f %.4f"%(drop,C,mu,k,bad,bnd,bnd/bad))
print()
print("### (3) the dividing line: one factor converges, the 13-fold product does not ###")
print("  one speed  : error term 2C/(7k) -> 0 as k grows.  C is fixed by the core.")
for drop in [1]:
    P=[x for x in range(1,13) if x!=drop]
    iv=safe_set(P); C=len(iv)
    for k in [157,1000,10000,100000,1000000]:
        print("     k=%-9d  error 2C/(7k) = %.3e"%(k,2*C/(7*k)))
print("  product    : relation-support ladder terms (THM-1061 V) GROW:")
print("     w2=+1.12  w3=-5.23  w4=+12.06   (tight family; partial sums 2.12, -3.11, +8.95 vs true 0)")
print()
print("  => soft Weyl is a ONE-FREQUENCY tool.  It proves death-star's position lemma, closes")
print("     boxeph's density off-resonance, and underwrites my measure horn and klein's")
print("     recursion -- all single-speed statements.  It does NOT bound the 13-fold product.")
print("DONE")
