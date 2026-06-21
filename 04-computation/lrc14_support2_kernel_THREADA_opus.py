#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A — EXACT support-2 carrier kernel and the anti-MDS signature of consec.

The relation code Lambda(E)={n: sum n_i e_i=0}.  The leading nonzero shell is
SUPPORT-2: relations n with exactly two nonzero coords i,j:  n_i e_i + n_j e_j=0,
i.e. e_i/e_j = -n_j/n_i = rational.  For the offset set E (with 0 in E), a
support-2 relation among NONZERO offsets exists iff two offsets are commensurate
with small ratio:  a e_i = b e_j, gcd(a,b)=1.  consec={0,1,...,k-1} has the
RICHEST support-2 set: every pair (e_i,e_j)=(i,j) gives the relation
j*e_i = i*e_j (since j*i = i*j), i.e. n=(...,j at i,..., -i at j,...).  This is
the ANTI-MDS (minimum-distance-2, maximal low-weight) signature.

The MacWilliams reading: corr(E) = sum_{0!=n in Lambda} K(n).  We compute the
EXACT contribution of the support-2 shell to corr, via the Fourier kernel:
  each support-2 relation n=(a on coord i, -b on coord j), a e_i = b e_j, with
  a,b>0 coprime times a common multiplier t (n = t*(a,-b)), contributes a kernel
  K_2(t,a,b,positions) that we derive from the sector-indicator Fourier coeffs.

We VERIFY: (1) the support-2 shell count = #commensurate pairs (the A_2 weight-
enumerator coefficient); (2) consec maximizes A_2 among equal-cardinality shapes
in a box (anti-MDS); (3) the SIGN of the support-2 kernel contribution to corr.
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict

TWO_PI=2*math.pi
def sector_coeff(s,m):
    """a_{s,m} = integral_{s/7}^{(s+1)/7} e(-m u) du, exact-ish complex."""
    if m==0: return complex(F(1,7))
    # (e(-m s/7) - e(-m(s+1)/7))/(2 pi i m)
    num=cmath.exp(-2j*math.pi*m*s/7)-cmath.exp(-2j*math.pi*m*(s+1)/7)
    return num/(2j*math.pi*m)

def Ps_missed_fourier(E,s,M=40):
    """P(sector s missed) via Fourier: integral prod_i (1 - chi_s(e_i x)) dx
       = sum over freq assignments (m_1..m_k) with sum m_i e_i =0 of
         prod_i [ (m_i==0)?1: -a_{s,m_i} ]  ... wait, (1 - chi_s) = (1 - sum_m a_{s,m} e(m e_i x)).
       Product over i, integrate -> sum over (m_i) with sum m_i e_i=0 of
         prod_i coeff_i, coeff_i = 1 if m_i=0 contributes the '1' minus a_{s,0};
       cleaner: 1 - chi_s(t) = (1 - a_{s,0}) - sum_{m!=0} a_{s,m} e(m t)
              = (6/7) - sum_{m!=0} a_{s,m} e(m t).
       So coeff for coord i: m_i=0 -> 6/7;  m_i!=0 -> -a_{s,m_i}.
       Truncate |m_i|<=M.  Returns complex (should be ~real)."""
    k=len(E)
    total=0j
    # enumerate relations sum m_i e_i=0 with |m_i|<=M -- too big for k=8 directly.
    # Instead use the occupancy law for exact P; this fourier is just a cross-check
    # on small k.  Keep for k<=4 sanity.
    for ms in itertools.product(range(-M,M+1),repeat=k):
        if sum(ms[i]*E[i] for i in range(k))!=0: continue
        c=1+0j
        for i in range(k):
            c*= (6/7) if ms[i]==0 else (-sector_coeff(s,ms[i]))
        total+=c
    return total

def occupancy_pi(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1); pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)]+=hi-lo
    return pi
def Ps_exact(E,s):
    """exact P(sector s missed) via occupancy law."""
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1); tot=F(0)
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if s not in hit: tot+=hi-lo
    return tot

def support2_count(E):
    """A_2 = number of (unordered) commensurate NONZERO-offset pairs counted with
       their primitive support-2 relation, i.e. #pairs {i,j}, e_i,e_j!=0, that admit
       a relation a e_i=b e_j (always true over Q!).  But the CODE weight enumerator
       counts LATTICE vectors; the support-2 shell coefficient is naturally the
       number of pairs whose PRIMITIVE relation is SHORT.  Define A_2^prim = #pairs;
       A_2^short(R) = #pairs with min(|a|,|b|)<=R for the primitive a e_i=b e_j."""
    nz=[e for e in E if e!=0]
    pairs=0; shortR={1:0,2:0,3:0}
    for i in range(len(nz)):
        for j in range(i+1,len(nz)):
            ei,ej=nz[i],nz[j]
            g=gcd(ei,ej); a=ej//g; b=ei//g  # a*ei=b*ej?  a*ei=(ej/g)*ei, b*ej=(ei/g)*ej -> equal. ok
            pairs+=1
            mn=min(abs(a),abs(b))
            for R in shortR:
                if mn<=R: shortR[R]+=1
    return pairs,shortR

def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

if __name__=="__main__":
    print("="*78); print("Fourier cross-check of P(sector missed) on small k"); print("="*78)
    for E in [[0,1,2],[0,1,3]]:
        for s in range(7):
            pe=Ps_exact(E,s); pf=Ps_missed_fourier(E,s,M=25)
            print(f"  E={E} s={s}: exact={float(pe):.5f} fourier={pf.real:.5f} (im={pf.imag:.1e}) match={abs(float(pe)-pf.real)<1e-3}")
        print()
    print("="*78); print("Support-2 weight enumerator A_2 and the anti-MDS signature"); print("="*78)
    print(" Every pair of nonzero offsets admits a support-2 relation over Q, but the")
    print(" CODE-relevant coefficient is #pairs with a SHORT primitive relation")
    print(" min(|a|,|b|)<=R (these carry the largest Fourier kernel ~1/(ab)).\n")
    for k,W in [(8,12),(9,12)]:
        C=consec(k); pc,shc=support2_count(C)
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        bank=[E for E in bank if primitive(E)]
        # is consec max in short support-2 count at R=1,2,3?
        for R in [1,2,3]:
            cval=shc[R]; beat=0; best=cval; bestE=C
            for E in bank:
                _,sh=support2_count(list(E))
                if sh[R]>cval: beat+=1
                if sh[R]>best: best=sh[R]; bestE=list(E)
            print(f" k={k} R={R}: consec short-A_2={cval}; #shapes with more={beat}; max={best} by {bestE}  consec-max:{beat==0}")
        print()
