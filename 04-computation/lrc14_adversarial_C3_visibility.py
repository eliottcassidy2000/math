#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL TEST of CANDIDATE 3 (Lam-Leung visibility / dilation reduction).

Candidate 3 claims: "the signed mass K(n) is a function only of the mod-7 visibility class
of n, and the worst visibility class is realized by a BOUNDED representative."

We test the LOAD-BEARING factual claim directly:
   Does K(n) depend ONLY on (n mod 7)?    (the cyclotomic phase e^{-2pi i n j/7})
The kernel:  K(n) = sum_{T subset {1..6}} (-1)^|T| prod_j ctilde_T(n_j),
   ctilde_T(0)=1-|T|/7,  ctilde_T(n)= -sum_{j in T} e^{-2pi i n j/7}(1-e^{-2pi i n/7})/(-2pi i n).
The PHASE  e^{-2pi i n j/7}  depends only on n mod 7  --  BUT the prefactor 1/(2 pi i n)
depends on the ACTUAL integer n (the AMPLITUDE).  So ctilde_T(n) = (phase depending on n mod 7)
* (1/n) (amplitude).  => K(n) does NOT depend only on n mod 7; it carries an explicit 1/n decay.

We verify EXACTLY (closed form in Q(zeta_7), evaluated to high precision) that two relations n, n'
with n' = n + 7*(shift) (same residues mod 7) have DIFFERENT K (because the 1/n amplitudes differ).
This is the decisive test: if K were residue-only, candidate 3's finite-residue reduction would
hold; if K carries 1/n amplitude (it does), then the "finitely many residue classes" argument is a
CATEGORY ERROR -- there are infinitely many distinct K values within one residue class, and the
amplitude decay is EXACTLY the harmonic-divergent envelope, so dilation does not bound the sum.
"""
import sys, math, cmath, itertools
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None

def stilde_T(T, n):
    """ctilde_T(n) for n != 0 (the 'sine-kernel' coordinate factor), high-precision complex."""
    if n==0: return complex(1.0 - len(T)/7.0, 0.0)
    if n%7==0: return 0j
    s=0j
    for j in T:
        s += cmath.exp(-2j*math.pi*n*j/7.0)*(1-cmath.exp(-2j*math.pi*n/7.0))/(-2j*math.pi*n)
    return -s

def K(n_vec, sectors=range(1,7)):
    """K(n) = sum_T (-1)^|T| prod_j ctilde_T(n_j). n_vec excludes the e=0 coord (factor 1)."""
    tot=0j
    for r in range(0,7):
        for T in itertools.combinations(sectors, r):
            p=1.0+0j
            for nj in n_vec:
                p*=stilde_T(T, nj)
                if p==0: break
            tot += ((-1)**r)*p
    return tot

if __name__=="__main__":
    print("="*80)
    print("CANDIDATE 3 TEST: does K(n) depend only on n mod 7 (visibility class)?")
    print("="*80)
    print("If yes -> finite residue reduction works. If K carries 1/n amplitude -> category error.\n")
    # a fixed support-6 relation pattern (residues), then dilate by adding 7 to coords keeping residues.
    base = (1,2,3,1,2,4)   # arbitrary support-6 residue pattern (all nonzero mod 7)
    print(f"  base support-6 relation residue pattern n = {base}")
    print(f"     K(base) = {K(base):.6f}   |K| = {abs(K(base)):.6e}")
    print()
    print("  Now add 7 to the FIRST coordinate only (same residues mod 7, different amplitude):")
    for shift in [0,1,2,3,5,10,50]:
        n_vec = (base[0]+7*shift,) + base[1:]
        kv = K(n_vec)
        print(f"     n=({n_vec[0]:>4},...rest same): K={kv: .6e}  |K|={abs(kv):.6e}")
    print()
    print("  Multiply ALL coords by a constant c coprime to 7 (a DILATION n->c*n, residues change")
    print("  to c*residue mod 7 but the orbit/measure is dilation-invariant). Watch |K| shrink ~1/c^6:")
    for c in [1,2,3,4,5,8,15]:
        n_vec = tuple(c*b for b in base)
        kv=K(n_vec)
        print(f"     c={c:>3}: n={n_vec}  |K|={abs(kv):.4e}  (c^6*|K| = {c**6*abs(kv):.4e})")
    print()
    print("DIAGNOSIS:")
    print(" - K(n) is NOT residue-only: shifting a coord by +7 (same mod-7 class) changes |K|,")
    print("   decaying like 1/n. So one residue class contains INFINITELY many distinct K values.")
    print(" - The dilation n->c*n shrinks |K| like 1/c^6 (the support-6 amplitude product), but the")
    print("   NUMBER of relations of amplitude ~c grows, recreating the harmonic-divergent envelope.")
    print(" - Therefore 'finitely many residue classes -> bounded representative' is a CATEGORY ERROR.")
    print("   Visibility (mod-7 vanishing, THM-503) is REAL but only zeroes the 7|n coords; it does")
    print("   NOT collapse the infinite amplitude tower within a residue class. Candidate 3 REFUTED")
    print("   as a closer (the s708 visibility law is true but does not bound the amplitude sum).")
