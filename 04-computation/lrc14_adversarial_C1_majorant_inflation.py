#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL TEST of CANDIDATE 1 (Selberg-Beurling bandlimited per-factor majorant).

Candidate 1's plan: replace each factor 1_{B_T}(e_i x) by a degree-N Vaaler/Beurling
MAJORANT V^+ (finite Fourier support |n|<=N).  Then
  meas(S7(E)) <= integral prod V^+  =  sum over Lambda(E) cap [-N,N]^{k-1} of prod hatV^+(n_j),
a FINITE lattice sum, uniform in E (box size N depends only on geometry, not E).

THE FATAL FLAW (we make it exact). The Vaaler majorant inflates the 0-mode of each factor:
  hatV^+(0) = |B_T| + 1/(N+1).
The product of k such factors has 0-mode (|B_T|+1/(N+1))^? -> the main term M7 is replaced by
  M7_infl(k,N) = sum_T (-1)^|T| (1 - |T|/7 + 1/(N+1))^{k-1}   [the k-1 nonzero offsets carry it].
Even with a PERFECT high-mode majorant (all n!=0 modes EXACT), this main-term inflation alone
must stay below cap_k.  We check EXACTLY (Fractions) whether M7_infl(k,N) <= cap_k, and at what N.

ADDITIONAL FATAL FLAW (signed product). meas(S7)=sum_T (-1)^|T| (prod of nonneg factors).
For T with (-1)^|T|=+1 you need an UPPER bound on the product -> use V^+ on each factor.
For T with (-1)^|T|=-1 you need a LOWER bound on the (nonneg) product -> use V^- (minorant).
But V^- can be NEGATIVE on parts of the torus, so prod V^- is NOT a valid lower bound on a
product of indicators in general (a product of lower bounds bounds the product below only if all
factors are nonneg, which V^- need not be). We FLAG this; the inflation test alone already decides.
"""
import sys, math, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None

caps={8:F(2243,5880)}  # exact cap_8 from canon (THM-535)
cap_float={8:0.38153,9:0.49426,10:0.6044}

def M7_infl_exact(k,N):
    """exact rational main-term inflation with 0-mode bumped by 1/(N+1)."""
    bump=F(1,N+1)
    s=F(0)
    for t in range(0,7):
        base = F(7-t,7) + bump
        s += F((-1)**t * math.comb(6,t)) * base**(k-1)
    return s

def M7_exact(k):
    s=F(0)
    for t in range(0,7):
        s += F((-1)**t * math.comb(6,t)) * F(7-t,7)**(k-1)
    return s

if __name__=="__main__":
    print("="*80)
    print("CANDIDATE 1 TEST: best-case main-term inflation of the Beurling-Selberg majorant")
    print("="*80)
    print("Even a PERFECT high-mode majorant inflates the main term to M7_infl(k,N).")
    print("For the certificate to work we need M7_infl(k,N) <= cap_k (necessary condition).\n")
    for k in [8,9,10]:
        capk = caps.get(k)
        capk_f = float(capk) if capk is not None else cap_float[k]
        m7=float(M7_exact(k))
        print(f"  k={k}: M7(exact)={m7:.5f}  cap_{k}={capk_f:.5f}  (headroom for inflation = {capk_f-m7:.5f})")
        # find smallest N making M7_infl <= cap (exact comparison when cap exact)
        firstok=None
        for N in [3,6,12,24,48,96,200,500,1000,5000]:
            mi=M7_infl_exact(k,N)
            ok = (mi<=capk) if capk is not None else (float(mi)<=capk_f)
            print(f"      N={N:>5}: M7_infl={float(mi):.5f}  <= cap? {ok}")
            if ok and firstok is None: firstok=N
        print(f"      -> inflation drops below cap at N>= {firstok}\n")

    print("="*80)
    print("BUT: the inflation test is NECESSARY, not sufficient. The REAL killer is that the")
    print("FINITE-box lattice sum at that N must reproduce the TRUE correction delta_k=0.303 (k=8)")
    print("to within (cap - meas) = 0.054.  The earlier Q1 diagnostic (BS_bandlimited_signed)")
    print("showed the SIGNED truncated lattice sum does NOT converge to delta_k even at N=11")
    print("(it sat near -0.03 vs the true +0.303). A majorant only makes the bound LARGER than")
    print("the already-non-converging signed sum, so the box must be enormous AND the inflation")
    print("re-enters. The two requirements (small inflation -> large N; finite box -> tractable N)")
    print("are in tension, and crucially NO ONE HAS EXHIBITED a single (E,N) where the EXACT")
    print("majorant integral <= cap with a verified finite box. The mechanism is GROUNDED but")
    print("UNEXECUTED -- same status as the Minkowski count (MISTAKE-078).")
