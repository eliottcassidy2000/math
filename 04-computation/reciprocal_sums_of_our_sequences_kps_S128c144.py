#!/usr/bin/env python3
"""
RECIPROCAL SUMS of the project's integer sequences — each sum(1/a_n) is a SUB-SERIES of the
harmonic series. kind-pasteur-2026-07-21-S128c144.  Owner: reciprocal of an integer sequence = a
subset of the harmonic numbers; 1+1/2+1/3+1/4+1/5 > 2 already, while sum 1/T_n = 2 exactly; study
our sequences' reciprocal sums extensively and extend.

Core provable result: the FIGURATE LADDER.  For k>=2,  sum_{n>=k} 1/C(n,k) = k/(k-1)
  (telescoping 1/C(n,k) = k/(k-1)[1/C(n-1,k-1) - 1/C(n,k-1)]).
 k=2 (arcs=triangular) -> 2 (the owner's value); k=3 (tetrahedral) -> 3/2; ... -> 1 as k->inf.
 The HARMONIC series is the k=1 rung C(n,1)=n -> DIVERGES: the ladder's divergent origin.
So the project's DIMENSIONAL ladder vertices(n) -> arcs(C(n,2)) -> C(n,3) -> ... is EXACTLY the
figurate reciprocal ladder crossing from divergence (harmonic, dim 1) to convergence (2, dim 2).
"""
import mpmath as mp
mp.mp.dps = 40
from fractions import Fraction

# ---------------- sequence generators (as callables giving many terms) ----------------
def binom(n,k):
    from math import comb
    return comb(n,k)

SEQS = {}   # name -> (list of a_n, index note, growth tag)
# figurate / simplicial
SEQS["vertices n (dim1)"]      = ([n for n in range(1,4000)], "n>=1", "linear")
SEQS["arcs C(n,2)=T (dim2)"]   = ([binom(n,2) for n in range(2,2000)], "n>=2", "poly-2")
SEQS["tetrahedral C(n,3)"]     = ([binom(n,3) for n in range(3,1500)], "n>=3", "poly-3")
SEQS["C(n,4) (dim4)"]          = ([binom(n,4) for n in range(4,1200)], "n>=4", "poly-4")
SEQS["C(n,5) (dim5)"]          = ([binom(n,5) for n in range(5,1000)], "n>=5", "poly-5")
SEQS["squares n^2"]            = ([n*n for n in range(1,4000)], "n>=1", "poly-2")
SEQS["var(lam^2)=2C(n,3)"]     = ([2*binom(n,3) for n in range(3,1500)], "n>=3", "poly-3")
SEQS["central binom C(2n,n)"]  = ([binom(2*n,n) for n in range(0,60)], "n>=0", "exp-4^n")
SEQS["Catalan"]                = ([binom(2*n,n)//(n+1) for n in range(0,60)], "n>=0", "exp")
# 2-adic / lacunary (theta-type)
SEQS["2^C(n,2) (lab.tourn)"]   = ([2**binom(n,2) for n in range(1,40)], "n>=1", "lacunary")
SEQS["2^C(n-1,2) (switch cls)"]= ([2**binom(n-1,2) for n in range(1,40)], "n>=1", "lacunary")
SEQS["2^n (Cayley-Dickson)"]   = ([2**n for n in range(1,80)], "n>=1", "geometric")
SEQS["2^n-1 (Mersenne)"]       = ([2**n-1 for n in range(1,80)], "n>=1", "geometric")
# census / combinatorial (super-exponential)
SEQS["A000568 tournaments"]    = ([1,1,1,2,4,12,56,456,6880,191536,9733056,903753248], "n>=1", "super-exp")
SEQS["A002854 even graphs V(E)"]=([1,1,2,3,7,16,54,243,2038,33120,1182004], "n>=1", "super-exp")
SEQS["A000571 score seqs"]     = ([1,1,1,2,4,9,22,59,167,490,1486,4639,14805,48107], "n>=1", "exp~4^n")
SEQS["A000182 tangent"]        = ([1,2,16,272,7936,353792,22368256,1903757312], "n>=1", "super-exp")
SEQS["A000364 secant/Euler"]   = ([1,1,5,61,1385,50521,2702765,199360981], "n>=0", "super-exp")
SEQS["factorial n!"]           = ([1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800], "n>=0","super-exp")
SEQS["Fibonacci"]              = ([1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765], "n>=1","exp-phi")
# arithmetic / "edge" sequences from the H-work
SEQS["odd numbers (H parity)"] = ([2*n-1 for n in range(1,4000)], "n>=1", "linear")
SEQS["H-spectrum odds\\{7,21}"] = ([x for x in (2*n-1 for n in range(1,4000)) if x not in (7,21)], "n>=1", "linear")
SEQS["pentagonal n(3n-1)/2"]   = ([n*(3*n-1)//2 for n in range(1,3000)], "n>=1", "poly-2")

def rec_sum(vals):
    s = mp.mpf(0)
    for a in vals:
        if a>0: s += mp.mpf(1)/mp.mpf(a)
    return s

print("="*84)
print("RECIPROCAL SUMS  sum(1/a_n)  — each is a sub-series of the harmonic series")
print("="*84)
print(f"{'sequence':30s} {'index':7s} {'growth':11s} {'sum 1/a_n':>22s}  note")
print("-"*84)
KNOWN = {
 "arcs C(n,2)=T (dim2)": ("= 2 exactly (figurate k=2)"),
 "tetrahedral C(n,3)": ("= 3/2 exactly (figurate k=3)"),
 "C(n,4) (dim4)": ("= 4/3 exactly (figurate k=4)"),
 "C(n,5) (dim5)": ("= 5/4 exactly (figurate k=5)"),
 "var(lam^2)=2C(n,3)": ("= 3/4 (= (1/2)*3/2)"),
 "squares n^2": ("= pi^2/6 (Basel)"),
 "factorial n!": ("= e (n>=0)"),
 "2^n (Cayley-Dickson)": ("= 1 (n>=1)"),
 "vertices n (dim1)": ("DIVERGES (harmonic) — the edge"),
 "odd numbers (H parity)": ("DIVERGES (~ (1/2)ln)"),
 "H-spectrum odds\\{7,21}": ("DIVERGES (removing 2 terms keeps divergence)"),
}
for name,(vals,idx,growth) in SEQS.items():
    if growth=="linear":
        note = KNOWN.get(name,"DIVERGES")
        # show partial sum to N to display slow growth
        partial = rec_sum(vals[:len(vals)])
        print(f"{name:30s} {idx:7s} {growth:11s} {mp.nstr(partial,10):>22s}  {note} (partial, N={len(vals)})")
    else:
        s = rec_sum(vals)
        note = KNOWN.get(name,"")
        print(f"{name:30s} {idx:7s} {growth:11s} {mp.nstr(s,18):>22s}  {note}")

print()
print("="*84)
print("(1) THE FIGURATE LADDER  sum_{n>=k} 1/C(n,k) = k/(k-1)  — VERIFIED + telescoping proof")
print("="*84)
for k in range(2,9):
    s = rec_sum([binom(n,k) for n in range(k, 4000)])
    exact = Fraction(k,k-1)
    print(f"  k={k}: sum 1/C(n,k) = {mp.nstr(s,15)}   vs  k/(k-1) = {exact} = {float(exact):.10f}  match={abs(float(s)-float(exact))<1e-9}")
print("  telescoping: 1/C(n,k) = k/(k-1) [1/C(n-1,k-1) - 1/C(n,k-1)]  (identity, proven in THM)")
print("  k=1 (harmonic, vertices): k/(k-1) -> infinity  = the DIVERGENT origin of the ladder.")
print("  k->inf: k/(k-1) -> 1.  Ladder of sums: inf, 2, 3/2, 4/3, 5/4, ... -> 1.")

print()
print("="*84)
print("(2) THE 2-ADIC / LACUNARY (theta) sums — partial theta values")
print("="*84)
for name in ["2^C(n,2) (lab.tourn)","2^C(n-1,2) (switch cls)"]:
    vals=SEQS[name][0]; s=rec_sum(vals)
    print(f"  {name}: sum = {mp.nstr(s,25)}")
# Jacobi theta reference: sum_{n>=1} q^{n(n-1)/2} at q=1/2  (arcs 2^{-C(n,2)})
q=mp.mpf(1)/2
theta_arc = mp.nsum(lambda n: q**(n*(n-1)//2 if False else mp.mpf(n)*(n-1)/2), [1, mp.inf])
print(f"  cross-check sum_ n>=1  2^(-n(n-1)/2) = {mp.nstr(theta_arc,25)}  (partial theta at q=1/2)")
# Euler pentagonal product at q=1/2:  prod (1-2^-n) = signed pentagonal theta
prodE = mp.mpf(1)
for n in range(1,200): prodE *= (1-mp.mpf(2)**(-n))
print(f"  Euler prod_(n>=1)(1-2^-n) = {mp.nstr(prodE,25)}  (= signed pentagonal theta, pentagonal-number thm)")

print()
print("="*84)
print("(3) NAMED CONSTANTS the sums hit (fingerprints)")
print("="*84)
checks = [
 ("Fibonacci", "reciprocal Fibonacci const psi ~ 3.359885666"),
 ("Catalan", "sum 1/Catalan = 2 + 4pi/(9 sqrt3)?"),
 ("central binom C(2n,n)", "= 4/3 + 2 pi sqrt3 / 27"),
 ("2^n-1 (Mersenne)", "Erdos-Borwein constant ~ 1.606695"),
 ("A000568 tournaments", "converges fast (super-exp)"),
]
for name,desc in checks:
    s=rec_sum(SEQS[name][0]); print(f"  {name:26s} = {mp.nstr(s,18)}   {desc}")
# explicit closed-form checks
print(f"    check central-binom vs 4/3+2pi*sqrt3/27 = {mp.nstr(mp.mpf(4)/3 + 2*mp.pi*mp.sqrt(3)/27,18)}")
print(f"    check squares vs pi^2/6 = {mp.nstr(mp.pi**2/6,18)}")
print(f"    check factorial vs e = {mp.nstr(mp.e,18)}")
print("DONE.")
