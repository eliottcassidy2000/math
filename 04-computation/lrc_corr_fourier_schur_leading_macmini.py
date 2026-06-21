#!/usr/bin/env python3
"""
lrc_corr_fourier_schur_leading_macmini.py
THREAD A, part (iii)/deliverable: VERIFY the offset-relation-lattice Fourier identity for corr(E),
and isolate the SCHUR-TRIPLE (support-3 additive-triangle) LEADING TERM.

corr(E) = measS7(E) - iid_k
        = sum_{ 0 != n in Lambda(E) }  W(n)        (offset-relation-lattice sum)
where for a frequency vector n=(n_e)_{e in E} with sum_e n_e e = 0,
   W(n) = sum_{T subset {0..6}} (-1)^|T| prod_e a_{n_e}(T),
   a_m(T) = m-th Fourier coeff of g_T(t)=1[floor(7 frac t) not in T]
          = (1/7) sum_{j not in T} sinc-type:  a_m(T)= integral_0^1 g_T(t) e(-m t) dt.
   For m=0: a_0(T)=(7-|T|)/7. For m!=0: a_m(T)= sum_{j not in T} c_{j,m}, with
   c_{j,m}= integral over sector j of e(-mt) = e(-m j/7)*(1-e(-m/7))/(2 pi i m)  (a fixed 7th-root weight).

This script:
 (A) DIRECTLY verifies corr(E) == sum over Lambda(E) Fourier products, on small E (numerically,
     truncating |n_e|<=N), confirming the lattice picture.
 (B) Isolates the SUPPORT-3 (Schur-triple) part: relations n with exactly 3 nonzero entries
     {+1,+1,-1} on a Schur triple a+b=c. This is the LEADING cross-block term.
 (C) Tests whether corr ~ const + (leading-Schur-term) so that |corr - iid-floor part| is
     CONTROLLED by the Schur-triple Fourier sum -> a bound for the multi-block case.
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

TWOPI=2*math.pi
def e_(t): return cmath.exp(2j*math.pi*t)

# exact measS7 for ground truth
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for mm in range(0,7*e+1): bps.add(F(mm,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            if e==0: continue
            secs.add(int(((e*xm)%1)*7))
        if 0 in E: secs.add(0)
        if len(secs)==7: total+=x1-x0
    return total

from math import comb
def iid_k(k):
    # iid floor of COVER (all 7 sectors hit) = 7! S(k,7)/7^k (surjective occupancy)
    return sum((-1)**j*comb(7,j)*F(7-j,7)**k for j in range(8))

# Fourier coeff a_m(T): integral_0^1 1[sector(t) not in T] e(-m t) dt
def a_m(m, T):
    Tset=set(T)
    if m==0:
        return complex((7-len(Tset))/7.0, 0.0)
    s=0j
    for j in range(7):
        if j in Tset: continue
        # sector j = [j/7,(j+1)/7); integral e(-m t) dt = [e(-m t)/(-2pi i m)]_{j/7}^{(j+1)/7}
        hi=e_(-m*(j+1)/7); lo=e_(-m*j/7)
        s += (hi-lo)/(-2j*math.pi*m)
    return s

def corr_lattice(E, N):
    """Sum over Lambda(E)\{0} truncated to |n_e|<=N of Fourier product, summed over T-subsets.
       Returns complex (should be ~ real == corr)."""
    E=sorted(set(E)); k=len(E)
    # precompute a_m(T) for needed m and all T
    allT=[]
    for r in range(0,8):
        for T in itertools.combinations(range(7), r):
            allT.append((T,(-1)**r))
    total=0j
    ranges=[range(-N,N+1)]*k
    for n in itertools.product(*ranges):
        if all(v==0 for v in n): continue
        if sum(n[i]*E[i] for i in range(k))!=0: continue
        # this n is a relation; add its weighted Fourier product
        w=0j
        for T,sgn in allT:
            prod=1.0+0j
            for i in range(k):
                prod*=a_m(n[i],T)
                if prod==0: break
            w+=sgn*prod
        total+=w
    return total

def schur_relations(E):
    """All support-3 relations n with entries {+1,+1,-1} on a Schur triple a+b=c (a,b,c in E)."""
    E=sorted(set(E)); idx={e:i for i,e in enumerate(E)}; k=len(E); rels=[]
    Eset=set(E)
    for a,b in itertools.combinations_with_replacement(E,2):
        c=a+b
        if c in Eset and c!=a and c!=b and a!=b:  # distinct triple
            n=[0]*k; n[idx[a]]+=1; n[idx[b]]+=1; n[idx[c]]-=1
            rels.append(tuple(n))
    # dedup
    return list(set(rels))

def corr_schur_term(E):
    """The support-3 (Schur-triple) part of the lattice sum: sum over Schur relations and their
       negatives of the Fourier-product weight. The LEADING cross-block correction."""
    E=sorted(set(E)); k=len(E)
    allT=[]
    for r in range(0,8):
        for T in itertools.combinations(range(7), r):
            allT.append((T,(-1)**r))
    rels=schur_relations(E)
    # include negatives (n and -n both in lattice)
    allrels=set(rels)
    for n in rels:
        allrels.add(tuple(-v for v in n))
    total=0j
    for n in allrels:
        w=0j
        for T,sgn in allT:
            prod=1.0+0j
            for i in range(k):
                prod*=a_m(n[i],T)
                if prod==0: break
            w+=sgn*prod
        total+=w
    return total

if __name__=="__main__":
    print("="*92)
    print("(A) VERIFY corr(E) == lattice Fourier sum (truncated |n_e|<=N)")
    print("="*92)
    tests=[("tri[0,1,2,3]",[0,1,2,3]),("[0,1,3]",[0,1,3]),("[0,1,2,3,4]",[0,1,2,3,4]),
           ("twoblock[0,1,2,3,10,11]"[:24],[0,1,2,3,10,11])]
    for name,E in tests:
        exact=float(measS7(E)-iid_k(len(set(E))))
        for N in (3,6,10):
            lat=corr_lattice(E,N)
            print(f"  {name:<24} N={N:>2}  lattice={lat.real:+.5f} (im {lat.imag:+.1e})  exact_corr={exact:+.5f}  diff={lat.real-exact:+.5f}")
        print()

    print("="*92)
    print("(B/C) SCHUR-TRIPLE LEADING TERM vs full corr (does support-3 dominate the correction?)")
    print("="*92)
    pop=[("consec[0..3]",[0,1,2,3]),("consec[0..4]",[0,1,2,3,4]),("consec[0..5]",list(range(6))),
         ("Sidon[0,1,3,7]",[0,1,3,7]),("dissoc[0,1,3,7,15]",[0,1,3,7,15]),
         ("twoblock[0,1,2,50,51,52]",[0,1,2,50,51,52]),
         ("twoblock_far[0,1,2,500,501,502]",[0,1,2,500,501,502])]
    print(f"  {'E':<32}{'#Schur':>7}{'schurTerm':>11}{'fullCorr':>11}{'iidFloorPart':>13}")
    for name,E in pop:
        sch=corr_schur_term(E).real
        full=float(measS7(E)-iid_k(len(set(E))))
        nrel=len(schur_relations(E))
        # "iid floor part" here = full corr minus schur term = the residual (support>=4 + support-2)
        print(f"  {name:<32}{nrel:>7}{sch:>+11.5f}{full:>+11.5f}{full-sch:>+13.5f}")
    print("\n  If schurTerm tracks fullCorr (residual small), the support-3 additive-triangle")
    print("  Fourier sum is the LEADING controllable term -> handle on the multi-block carrier error.")
