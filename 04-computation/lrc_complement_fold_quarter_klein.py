#!/usr/bin/env python3
"""
lrc_complement_fold_quarter_klein.py  --  klein-2026-07-01-S79

Apply the tournament FOLD thinking to LRC-14, abstractly.
Tournament (S78): cube folds by <sigma=complement, flip> = Klein-four -> FULL->HALF->QUARTER; an obstruction
(SC covering excess) DISSOLVES in the fold-adapted (half-address) coordinates (a coordinate artifact).

LRC analogue:
 - COMPLEMENT-FOLD: the antipode iota: t->1-t (S55) is the LRC complement; M(S), the danger pattern, and the
   lonely measure are ALL iota-invariant, so LRC folds to the half-circle [0,1/2].
 - QUARTER-FOLD: the phase space is Z/Phi6, Phi6=3*61 (Eisenstein primes). The complement -1 mod Phi6 = the
   antipode = 182 CRT-FACTORS: 182=(-1 mod3,-1 mod61)=62*121, product of PARTIAL complements. So
   <62,121>=Klein-four folds Z/Phi6 -> a 'quarter' of the phase space, exactly as <sigma,flip> did the cube.
 - THE TEST (does the S78 dissolution happen?): does the covering-min obstruction DISSOLVE under the CRT
   quarter-fold, i.e., does the binding factor over the primes 3,61? Compute M(construction) restricted to
   each modulus q: if the construction binds ONLY at the FULL Phi6 (not at 3 or 61), the obstruction is
   metric-irreducible (does NOT dissolve) -- unlike the tournament (coordinate artifact). That would EXPLAIN
   why LRC is harder: the covering-min lives at the composite Phi6, invisible to the CRT quarter.
"""
from fractions import Fraction as F
import numpy as np

def Mrestrict(S, q):
    """max over t=a/q of min_v ||v t|| -- the loneliness available at modulus q (exact)."""
    best=F(0)
    for a in range(q):
        m=min(min((v*a)%q, q-((v*a)%q)) for v in S)   # ||v a/q|| * q  (integer numerator)
        val=F(m,q)
        if val>best: best=val
    return best

if __name__=="__main__":
    n=14; Phi6=n*n-n+1; C=list(range(1,n-1))+[n*(n-1)]; Mc=F(n,Phi6)
    print(f"n={n} Phi6={Phi6}=3*61; construction C={C}; covering-min M_C=n/Phi6={Mc}={float(Mc):.5f}")

    print("\n(1) COMPLEMENT-FOLD (iota: t->1-t). M(S)=max_t min_v||vt|| is iota-invariant (||v(1-t)||=||vt||),")
    print("    so the whole problem folds to the half-circle [0,1/2]; the 2 binding atoms {t*,1-t*} fold to 1.")
    N=200000; t=np.arange(N)/N
    G=np.full(N,1.0)
    for v in C:
        x=(v*t)%1.0; G=np.minimum(G,np.minimum(x,1-x))
    # iota-symmetry check: G(t)=G(1-t)
    iota_sym=np.allclose(G, G[::-1] if N%2==0 else G, atol=1e-9) or np.max(np.abs(G-np.concatenate([[G[0]],G[1:][::-1]])))<1e-6
    print(f"    iota-symmetry of the loneliness profile G(t)=G(1-t): {np.max(np.abs(G[1:]-G[1:][::-1]))<1e-9} (fold to [0,1/2] valid)")

    print("\n(2) QUARTER-FOLD group: complement -1 mod Phi6 = 182 = 62*121 (partial complements mod 3, mod 61);")
    print(f"    <62,121> = {{1,62,121,182}} = Klein-four in (Z/{Phi6})*; iota=182=(-1 mod3,-1 mod61). This folds")
    print(f"    the phase space Z/Phi6 into a QUARTER, exactly as <sigma,flip> folded the tournament cube.")

    print("\n(3) THE TEST -- does the covering-min obstruction DISSOLVE under the CRT quarter (binding factor over 3,61)?")
    print(f"    M(construction) restricted to each modulus q (loneliness available at q):")
    for q in [2,3,6,7,14,61,183, Phi6]:
        mr=Mrestrict(C,q)
        flag = " <-- = M_C (the DEEP binding)" if mr==Mc else (" (covered/shallow)" if mr<Mc else "")
        print(f"      q={q:>4}: max_t min_v||vt|| = {mr} = {float(mr):.5f}{flag}")
    # is the binding ONLY at Phi6 (composite), not at the prime factors 3, 61?
    m3=Mrestrict(C,3); m61=Mrestrict(C,61); mP=Mrestrict(C,Phi6)
    print(f"    => modulus 3: {float(m3):.4f}; modulus 61: {float(m61):.4f}; modulus Phi6=3*61: {float(mP):.4f}=M_C")
    print(f"    The construction COVERS the prime moduli 3 and 61 (small M there) and binds ONLY at the COMPOSITE Phi6.")
    print(f"    => the covering-min is METRIC-IRREDUCIBLE at Phi6: it does NOT factor over the CRT quarter (primes 3,61).")
    print(f"    So the obstruction does NOT dissolve under the quarter-fold (unlike the tournament SC-cover, a coord artifact).")
    print(f"    This is the analytic non-factorization of THM-503 (L is NOT an Euler product), seen as: the covering-min")
    print(f"    lives at the DEEP composite modulus Phi6 (the 2nd CF convergent), invisible to the prime CRT factors.")
    print("\n    WHAT THE COMPLEMENT-FOLD DOES BUY: halves the problem to [0,1/2]; lonely measure iota-even => Verblunsky")
    print("    REAL (S66); the parity lemma (S55, odd D => #lonely even) is the iota-fold fixed-point count. The HALF")
    print("    fold is useful; the QUARTER (CRT) fold is a real group symmetry but the metric obstruction is Phi6-deep.")
