#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A — INNER (relation-code) leg of B_2, exact Weyl expansion.

GOAL: express the shape-dependent part of B_2 = 1 - S_1 + S_2 as a sum over the
relation code Lambda(E), with EXACT Fourier kernels, and show the consec win is
a relation weight-enumerator statement.

Single-sector miss probability via Fourier (EXACT, finite sum):
  Let chi_s(u) = 1_{u in [s/7,(s+1)/7) mod 1}, period-1 indicator of sector s.
  Its Fourier series: chi_s(u) = sum_{m in Z} a_{s,m} e(m u),
     a_{s,0}=1/7;  a_{s,m} = (1/(2 pi i m)) (e(-m s/7) - e(-m(s+1)/7)) for m!=0.
  "sector s NOT hit by e_i" = 1 - chi_s(e_i x).  "s missed by all" =
     prod_i (1 - chi_s(e_i x)).  P(s missed) = integral_0^1 prod_i(1-chi_s(e_i x)) dx.
  Expand the product: sum over subsets J of {1..k} of (-1)^{|J|}
     integral prod_{i in J} chi_s(e_i x) dx.  Each factor expands in m_i; the
     integral is nonzero iff sum_{i in J} m_i e_i = 0 (a RELATION supported on J).
  So P(s missed) = sum_{relations n supported anywhere} (-1)^{supp(n)} *
     prod_{i: n_i!=0} a_{s,n_i}   (where n_i is the freq on coord i; n=0 gives
     the iid term (-1)^0 * ... no, n=0 term: J=empty -> +1; plus each single
     coord J={i} with m_i=0 excluded since a_{s,0} included only via... ).
  Cleanest: P(s missed) = prod_i (1 - chi_s(e_i x)) integrated; FULLY expand and
  collect by the integer relation vector n (n_i = total freq on coord i):
     P(s missed) = sum_{n: sum n_i e_i =0} (-1)^{supp(n)} prod_{i:n_i!=0} a_{s,n_i}
                   + (iid contributions where n_i=0 contribute a_{s,0}=1/7 with the
                     binomial expansion of (1-1/7)).
  We compute P(s missed) DIRECTLY (exact, via occupancy law) to avoid the bookkeeping,
  AND we compute the RELATION-CODE truncation to show it converges and is graded by
  support.  The point: S_1 = sum_s P(s missed); S_2 = sum_{s<t} P(s,t missed); the
  relation-code corrections to S_1 are support>=2 (single-coord relations need
  n_i e_i=0 -> n_i=0 since e_i!=0), and to S_2 likewise -> B_2's shape-dependence
  is a SUPPORT>=2 relation weight-enumerator with OUTER Krawtchouk g_2 signs.
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict

def occupancy_law(E):
    E=sorted(set(E))
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1)
    occ=defaultdict(lambda:F(0)); pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        occ[frozenset(hit)]+=hi-lo; pi[len(hit)]+=hi-lo
    return dict(occ),pi

def S_r_list(E):
    occ,pi=occupancy_law(E)
    return [sum(pi[h]*comb(7-h,r) for h in range(8)) for r in range(8)],pi

# iid baselines for S_1, S_2 (k iid uniform sectors):
def iid_S1(k): # E[N_iid] = sum_{s=0}^{6} P_iid(s missed) = 7*(6/7)^k
    return F(7)*F(6,7)**k
def iid_S2(k): # sum_{s<t} P_iid(both missed) = C(7,2)*(5/7)^k
    return F(comb(7,2))*F(5,7)**k
def iid_measS7(k):
    Snk=sum((-1)**(7-i)*comb(7,i)*i**k for i in range(8))
    return F(Snk,7**k)

def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

# relation shells by support, box B
def shells(E,B):
    E=list(E); k=len(E); nz=[i for i,e in enumerate(E) if e!=0]
    cnt=defaultdict(int)
    for n in itertools.product(*([range(-B,B+1)]*k)):
        if all(v==0 for v in n): continue
        if sum(n[i]*E[i] for i in range(k))!=0: continue
        supp=sum(1 for i in nz if n[i]!=0)
        if supp==0: continue
        cnt[supp]+=1
    return dict(cnt)

if __name__=="__main__":
    print("="*78)
    print("INNER LEG: B_2 shape-dependence vs iid, and the SUPPORT>=2 relation code")
    print("="*78)
    print(" B_2 = 1 - S_1 + S_2.  iid part B_2^iid = 1 - iid_S1 + iid_S2 (shape-free).")
    print(" corr_B2 = B_2 - B_2^iid = -(S_1-iid_S1) + (S_2-iid_S2) is the relation-code part.")
    print(" Single-coord relations are impossible (n_i e_i=0, e_i!=0 -> n_i=0), so the")
    print(" leading relation shell is SUPPORT-2 (e.g. consec: 2*e_1=1*e_2, i.e. 2*1=1*2).\n")
    for k,W,B in [(8,12,2),(9,12,2)]:
        C=consec(k)
        b2iid=1-iid_S1(k)+iid_S2(k)
        Src,pic=S_r_list(C); B2c=1-Src[1]+Src[2]; corrB2c=B2c-b2iid
        shc=shells(C,B)
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        bank=[E for E in bank if primitive(E)]
        print(f" k={k} span<= {W} box B={B}: {len(bank)} shapes;  B_2^iid={float(b2iid):.6f}")
        print(f"   consec: B_2={float(B2c):.6f}  corr_B2=B_2-iid={float(corrB2c):+.6f}  shells={dict(sorted(shc.items()))}")
        # correlate corr_B2 with support-2 and support-3 shell counts
        xs2=[];xs3=[];ys=[]
        for E in bank:
            Sr,pi=S_r_list(list(E)); b2=1-Sr[1]+Sr[2]; corrb2=float(b2-b2iid)
            sh=shells(list(E),B)
            xs2.append(sh.get(2,0)); xs3.append(sh.get(3,0)); ys.append(corrb2)
        def pear(a,b):
            n=len(a);ma=sum(a)/n;mb=sum(b)/n
            va=sum((x-ma)**2 for x in a);vb=sum((x-mb)**2 for x in b)
            return float('nan') if va==0 or vb==0 else sum((x-ma)*(y-mb) for x,y in zip(a,b))/((va*vb)**.5)
        print(f"   Pearson(corr_B2, #support-2 relations) = {pear(xs2,ys):+.4f}")
        print(f"   Pearson(corr_B2, #support-3 relations) = {pear(xs3,ys):+.4f}")
        # is corr_B2 maximized by consec? (== B_2 maximized, already verified)
        beat=sum(1 for y in ys if y>float(corrB2c)+1e-12)
        print(f"   #shapes beating consec corr_B2 = {beat}  (consec-max: {beat==0})\n")
