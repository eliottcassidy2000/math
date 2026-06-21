#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A part 2 (opus): EXACT carrier kernel K(n) and the MacWilliams link.

measS7(E) = meas{ x in [0,1) : every sector j/7..(j+1)/7, j=0..6, is hit by
some frac(e_i x) }.  Write it as an integral of a product of sector-coverage
indicators and EXPAND in Fourier characters e(m x).  Then integrate over x:
only frequency-zero terms survive, i.e. integer combinations sum_i n_i e_i = 0.

Sector-coverage indicator for sector s (s=0..6):
   1_{s hit} = 1 - prod_{i}(1 - 1_{floor(7 frac(e_i x)) = s}).
The full coverage indicator = prod_{s=0}^{6} 1_{s hit}.  Expanding the product
over s and the inner products over i gives, after Fourier expansion of each
1_{floor(7 frac(e_i x))=s} = sum_m c_m e(m e_i x) with
   c_m = integral over the sector = (1/7) sinc-type Fourier coeff:
   1_{u in [s/7,(s+1)/7)} has Fourier coeff a_{s,m} = integral_{s/7}^{(s+1)/7} e(-m u) du
        = e(-m s/7) * (1 - e(-m/7)) / (2 pi i m)   for m!=0,  = 1/7 for m=0.

So 1_{floor(7 frac(t))=s} as function of t (period 1) = sum_m a_{s,m} e(m t).

After full expansion measS7 = sum over frequency assignments that net to a
RELATION sum_i n_i e_i = 0 of a product of a_{s,m} coefficients.  Grouping by the
integer relation vector n in Lambda(E):
     measS7(E) = sum_{n in Lambda(E)} K(n),
where K(0)=iid_k and K(n) for n!=0 is the EXACT carrier kernel.

Rather than derive K(n) symbolically (messy), we VERIFY the structure and the
MacWilliams link numerically-exactly via a DIFFERENT, fully rigorous route:
the depth law itself, sliced by which sectors are hit.  Then we connect the
two via the MacWilliams transform on the *sector occupancy code*.

KEY NEW OBJECT (the genuine MacWilliams pairing):
   For a fixed x, define the SECTOR OCCUPANCY VECTOR v(x) in {0,1}^7,
   v(x)_s = 1 if sector s is hit.  This is a codeword-like object in F_2^7.
   The DEPTH h(x) = wt(v(x)) (Hamming weight).  measS7 = P(v = all-ones).
   pi_E(h) = the WEIGHT DISTRIBUTION of the random vector v(x), x ~ Unif.
   The Bonferroni/coverage weights g_J are EXACTLY partial sums of the
   MacWilliams/Krawtchouk transform on F_2^7:
        [v = all-ones] = (1/2^7) sum_{u in F_2^7} (-1)^{<u, v+1>}
   and grouping u by weight gives the Krawtchouk kernel K_j^{(7)}.

So measS7 = P(N=0) is LITERALLY the MacWilliams dual read of the occupancy
weight distribution pi_E.  We verify the exact Krawtchouk/MacWilliams identity:
     P(N=0) = (1/2^7) sum_{w=0}^{7} (-1)^w * (number of u of weight w with
              <u,v>=parity...) ...  -> = sum_h pi_E(h) * (Krawtchouk read).
The CLEAN statement: P(all hit) = sum_{S subset sectors} (-1)^{|S|}
P(all sectors in S are MISSED) = inclusion-exclusion = sum_J (-1)^J S_J, and
each S_J = <g_J, pi_E> with g_J the partial Krawtchouk.  This is the
MacWilliams transform of the occupancy distribution restricted to the "missed"
poset.  VERIFY exactly.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb

def occupancy_law(E):
    """Exact: for x ~ Unif[0,1), the distribution of the occupancy SET (which
       sectors hit).  Returns dict frozenset(hit sectors) -> measure (Fraction).
       Also returns depth law pi[h]."""
    E=sorted(set(E))
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1)
    from collections import defaultdict
    occ=defaultdict(lambda:F(0)); pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2
        hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        occ[frozenset(hit)]+=hi-lo
        pi[len(hit)]+=hi-lo
    return dict(occ), pi

def S_J_via_IE(E,J):
    """S_J = sum_{|S|<=J, S subset of 7 sectors, S!=empty handled} ... the
       coverage IE truncation computed DIRECTLY from the occupancy law:
       S_r = sum_{|S|=r} P(all sectors in S are MISSED).
       S_0=1; measS7 = sum_{r=0}^{7} (-1)^r S_r  (full IE)."""
    occ,pi=occupancy_law(E)
    # P(all sectors in S missed) = measure of x whose hit-set avoids S
    # = sum over occ sets T with T ∩ S = empty of occ[T].
    sectors=range(7)
    Sr=[F(0)]*8
    for r in range(8):
        tot=F(0)
        for Sset in itertools.combinations(sectors,r):
            Ss=set(Sset)
            m=F(0)
            for T,meas in occ.items():
                if Ss.isdisjoint(T): m+=meas
            tot+=m
        Sr[r]=tot
    # truncation
    return sum((-1)**r*Sr[r] for r in range(J+1)), Sr, pi

def measS7(E):
    occ,pi=occupancy_law(E); return pi[7]

if __name__=="__main__":
    print("="*78)
    print("EXACT MacWilliams / inclusion-exclusion identity for measS7")
    print("="*78)
    print(" measS7 = sum_{r=0}^{7} (-1)^r S_r,  S_r = sum_{|S|=r} P(sectors S all missed).")
    print(" And S_r = <C(7-h,r), pi_E(h)>  (the r-th elementary symmetric / Krawtchouk read).")
    print(" Verify S_r = sum_h pi_E(h) C(7-h, r) and full IE = pi_E(7) = measS7.\n")
    for E in [list(range(8)),[0,2,3,4,5,6,7,8],[0,1,2,3,5,7,9,11]]:
        full,Sr,pi=S_J_via_IE(E,7)
        # check S_r = sum_h pi[h] C(7-h,r)
        ok=all(Sr[r]==sum(pi[h]*comb(7-h,r) for h in range(8)) for r in range(8))
        print(f" E={E}")
        print(f"   S_r (r=0..7) = {[str(x) for x in Sr]}")
        print(f"   S_r == <C(7-h,r),pi> for all r:  {ok}")
        print(f"   full IE sum(-1)^r S_r = {full} = {float(full):.6f}   measS7=pi[7]={float(pi[7]):.6f}   match={full==pi[7]}")
        # the EVEN Bonferroni truncations B_2,B_4 and measS7=B_6
        for J in [2,4,6]:
            BJ=sum((-1)**r*Sr[r] for r in range(J+1))
            print(f"     B_{J} = {float(BJ):.6f}")
        print()
