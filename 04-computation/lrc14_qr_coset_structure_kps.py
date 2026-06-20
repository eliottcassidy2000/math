#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QR/NQR structure of the coset kernel D7 (NUMERIC scan, all 7776 cosets) +
the chi-twisted Gauss-sum collapse test.  Goal: does grouping the relation lattice
by the QR/NQR pattern of (n mod 7) make the signed sum cancel / telescope?

We compute D7(c) for every c in (F_7^*)^6 and ask:

  (i)  How does |D7(c)| depend on the QR/NQR pattern (#QR coords) of c?
  (ii) Sum over each QR/NQR *class* (cosets with a fixed QR/NQR pattern) -- does it
       vanish or collapse to a multiple of the Gauss sum g = i sqrt7?
  (iii) The KEY: the archimedean weights are real; the sum K(n) = arch(n) D7(n mod 7).
       Group n by residue class c.  Within a residue class, arch(n) is a REAL signed
       sum over the lattice.  D7(c) is a FIXED cyclotomic number.  So
          sum_{n: n=c mod7} K(n) = D7(c) * [ sum_{n=c} arch(n) ].
       The cancellation is then driven by (a) the arch lattice sum per residue, and
       (b) the cyclotomic D7(c).  We test whether sum_c D7(c) (unweighted) is a
       Gauss-sum multiple -- the cleanest possible collapse.
"""
import sys, itertools, cmath, math
from collections import defaultdict
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

Z=cmath.exp(2j*math.pi/7)
QR={1,2,4}; NQR={3,5,6}
def chi(m):
    m%=7
    return 0 if m==0 else (1 if m in QR else -1)

def D7_num(c):
    pref=1.0+0.0j
    for cj in c: pref*=(1 - Z**((-cj)%7))
    acc=0.0+0.0j
    for r in range(7):
        for T in itertools.combinations(range(1,7),r):
            p=1.0+0.0j
            for cj in c:
                # sigma_T(cj) = sum_{t in T} zeta^{-cj t}
                p*= sum(Z**((-cj*t)%7) for t in T)
            acc += ((-1)**r)*p
    return pref*acc

# faster: precompute sigma_T(m) for all T,m
def precompute():
    sig={}
    for r in range(7):
        for T in itertools.combinations(range(1,7),r):
            for m in range(1,7):
                sig[(T,m)]=sum(Z**((-m*t)%7) for t in T)
    return sig
SIG=precompute()
Tlist=[T for r in range(7) for T in itertools.combinations(range(1,7),r)]
SGN={T:(-1)**len(T) for T in Tlist}
PREF={m:(1 - Z**((-m)%7)) for m in range(1,7)}
def D7_fast(c):
    pref=1.0+0.0j
    for cj in c: pref*=PREF[cj]
    acc=0.0+0.0j
    for T in Tlist:
        p=1.0+0.0j
        for cj in c: p*=SIG[(T,cj)]
        acc+=SGN[T]*p
    return pref*acc

if __name__=="__main__":
    g=sum(chi(r)*Z**r for r in range(1,7))
    print(f"Gauss sum g={g:.6f}  (i sqrt7={1j*math.sqrt(7):.6f})  g^2={g*g:.4f}")

    # full scan
    S=0.0+0.0j
    by_pattern=defaultdict(complex)      # key = #QR coords
    by_chiprod=defaultdict(complex)      # key = chi(prod c_j) in {+1,-1}
    Schi_prod=0.0+0.0j                   # sum chi(prod) D7
    Schi_coord=[0.0+0.0j]*6              # sum chi(c_i) D7 for each coord i
    absmax=0.0; cnt=0
    # store magnitude stats by #QR
    magsum=defaultdict(float); magcnt=defaultdict(int)
    for c in itertools.product(range(1,7),repeat=6):
        d=D7_fast(c)
        S+=d; cnt+=1
        nqr=sum(1 for x in c if x in QR)
        by_pattern[nqr]+=d
        magsum[nqr]+=abs(d); magcnt[nqr]+=1
        cp=1
        for x in c: cp*=chi(x)
        by_chiprod[cp]+=d
        Schi_prod += cp*d
        for i in range(6): Schi_coord[i]+= chi(c[i])*d
        absmax=max(absmax,abs(d))

    print(f"\n#cosets={cnt}   max|D7|={absmax:.4f}")
    print(f"RAW   sum_c D7(c) = {S:.6f}     S/g = {S/g:.6f}")
    print(f"chi-prod twisted  sum chi(prod c) D7 = {Schi_prod:.6f}   /g = {Schi_prod/g:.6f}")
    for i in range(6):
        print(f"  chi(c_{i+1}) twisted sum = {Schi_coord[i]:.6f}  /g={Schi_coord[i]/g:.6f}")

    print(f"\nSum of D7 by #QR-coords pattern (0..6):")
    tot=0
    for k in range(7):
        print(f"   #QR={k}: sum={by_pattern[k]:.6f}   mean|D7|={magsum[k]/max(magcnt[k],1):.5f}  count={magcnt[k]}")
    print(f"\nSum by chi(prod c_j):  +1 -> {by_chiprod[1]:.6f}   -1 -> {by_chiprod[-1]:.6f}")
    print(f"   (+1)+(-1) = {by_chiprod[1]+by_chiprod[-1]:.6f}  = raw S")
    print(f"   (+1)-(-1) = {by_chiprod[1]-by_chiprod[-1]:.6f}  = chi-prod twisted")
