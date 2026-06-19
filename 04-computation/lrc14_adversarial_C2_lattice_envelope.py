#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL TEST of CANDIDATE 2's quantitative claim:
  "support>=6 forces a 6-fold lattice sum; Cauchy-Schwarz over 6 nonzero coords gives
   prod (sum 1/n_j^2)^{1/2} = (pi^2/6)^3 < 17, FINITE."

The claim conflates two index sets:
  (A) the FREE product Z^6 \ {0-coords} : sum over THAT of prod 1/|n_j| ... but with L^2
      (1/n^2) it would be (pi^2/3)^6 or so -- a product, finite, but IRRELEVANT.
  (B) the actual contributing set = the RELATION LATTICE Lambda(E) intersect {support>=6}.
      The kernel K(n) is summed over n in Lambda(E), a rank-(k-2) sublattice of Z^{k-1}.

The divergence MISTAKE-078 found (partial sum 7.42 at 10^5, growing like log) is the
support->=6 ABSOLUTE envelope sum_{n in Lambda(E)} prod_j |s_T(n_j)| with |s_T(n)|<=c1/|n|.
Candidate 2 claims this is actually a CONVERGENT 6-fold object via 1/n^2 L^2 decay.

KEY POINT WE TEST EXACTLY: on the relation lattice Lambda(E), the constraint sum n_j e_j = 0
TIES the coordinates together. A 1-parameter sublattice (e.g. scaling a single primitive
support-6 relation by t=1,2,3,...) gives n = t*v, and prod_j |s_T(n_j)| ~ prod c1/|t v_j|
= C/t^6, whose sum over t CONVERGES. BUT there are MANY primitive support-6 relations, and
their COUNT within |n|<=L grows polynomially. The real question: does sum over ALL
support-6 lattice points of prod c1/|n_j| converge?

We compute, for an EXACT wide set E, the partial envelope
  ENV(L) = sum_{0!=n in Lambda(E), |n|_inf<=L, support>=6} prod_{j in supp} c1/|n_j|
and watch ENV(L) vs L. If it keeps growing (esp. ~ powers of log or polynomial), Candidate 2's
"(pi^2/6)^3<17" is REFUTED as stated (Cauchy-Schwarz across coords does not apply to a
lattice-tied sum, and the per-coordinate 1/n^2 cannot be factored out).

We ALSO compute the EXACT signed correction (measS7 - M7) for the same E to show what the
TRUE quantity is (it is small and bounded -- that's the point, but the ENVELOPE is not the tool).
"""
import sys, itertools, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None

c1 = 0.6973026  # sup |sin(pi n /7)|/(pi |n|) scale (the per-coord envelope constant)

def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            y=(e*xm)%1; secs.add((y.numerator*7)//y.denominator)
        if len(secs)==7: total+=x1-x0
    return total

def M7(k):
    s=F(0)
    for t in range(0,7):
        s += F((-1)**t * math.comb(6,t)) * F(7-t,7)**(k-1)
    return s

def env_partial(E, L):
    """sum over n in Lambda(E)\{0}, |n|_inf<=L, with >=6 nonzero-and-not-7-mult coords,
       of prod over those coords of c1/|n_j|.  E has 0 as first element; the relation is
       sum_{j>=1} n_j e_j = 0 over the NONZERO offsets (e_0=0 contributes a free factor 1)."""
    Enz=[e for e in E if e!=0]
    kk=len(Enz)
    # enumerate n in [-L,L]^kk with sum n_j e_j = 0; count envelope of support>=6 (non-7-mult) coords
    total=0.0
    cnt=0
    for ns in itertools.product(range(-L,L+1), repeat=kk):
        if sum(n*e for n,e in zip(ns,Enz))!=0: continue
        nz=[n for n in ns if n!=0 and n%7!=0]
        if len(nz)<6: continue
        if len(ns)-ns.count(0) != len([n for n in ns if n!=0]):  # paranoia
            pass
        # only the genuinely-contributing coords (nonzero, not mult of 7) carry s_T factor;
        # a coord that is nonzero but mult-of-7 has s_T=0 -> whole term 0, skip it.
        if any(n!=0 and n%7==0 for n in ns):
            continue
        p=1.0
        for n in nz: p*=c1/abs(n)
        total+=p; cnt+=1
    return total, cnt

if __name__=="__main__":
    print("="*80)
    print("CANDIDATE 2 TEST: does the support>=6 lattice ENVELOPE sum converge?")
    print("="*80)
    print(f"  Candidate 2 claims the bound is (pi^2/6)^3 = {(math.pi**2/6)**3:.4f} < 17 (FINITE).")
    print(f"  (that is sum 1/n^2 = pi^2/6 = {math.pi**2/6:.5f}, cubed for 6 coords paired)")
    print()
    # Use sets with at least 7 nonzero offsets (k>=8) but SMALL so enumeration is feasible.
    # A genuinely WIDE primitive set and the AP, both k=8 (7 nonzero offsets -> rank-6 lattice).
    tests = [
        ("AP k=8 {0..7}", list(range(8))),
        ("wide k=8 {0,1,2,3,4,5,6,30}", [0,1,2,3,4,5,6,30]),
        ("wide k=8 {0,1,3,7,12,20,30,44} Sidon", [0,1,3,7,12,20,30,44]),
    ]
    for name,E in tests:
        g=measS7(E); m7=M7(len(E)); corr=g-m7
        print(f"  {name}: TRUE signed correction (measS7-M7) = {float(corr):+.5f}")
        for L in [7,10,13,16]:
            ev,cnt = env_partial(E,L)
            print(f"     |n|_inf<={L:>2}: support>=6 ENVELOPE sum = {ev:8.4f}  ({cnt} lattice pts)")
        print()
    print("DIAGNOSIS:")
    print(" - If ENVELOPE keeps growing with L (does not stabilize well below 17), then the")
    print("   '(pi^2/6)^3 < 17' Cauchy-Schwarz argument is a CATEGORY ERROR: the sum is over a")
    print("   RANK-6 LATTICE (coords tied by sum n_j e_j=0), not a free product Z^6, so you")
    print("   cannot factor sum 1/n_j^2 across coordinates.")
