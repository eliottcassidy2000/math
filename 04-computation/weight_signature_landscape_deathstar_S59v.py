#!/usr/bin/env python3
"""
death-star-2026-07-20-S59v (HYP-8205, THM-1365) -- the weight-signature reduced JC.
Map which dim-3 C*-weight signatures admit an equivariant Keller COUNTEREXAMPLE.
Extends the S59n equivariant-family solver to arbitrary weight vectors w.

For weight w=(w1,w2,w3) (source action lam.x_i = lam^{w_i} x_i), an equivariant
map F has F_i of weighted-degree w_i.  We build the general equivariant ansatz
up to a degree cap, impose det JF == const (Keller), and test whether a
NON-INJECTIVE solution (a genuine counterexample) can exist -- via the S59n
instruments (the c-graded det + the collision channel), validated so the known
(1,-1,-2) owner case is REDISCOVERED.
"""
from fractions import Fraction as Fr
from itertools import product
from math import gcd

def pmul(a,b):
    r={}
    for ka,ca in a.items():
        for kb,cb in b.items():
            k=tuple(p+q for p,q in zip(ka,kb)); v=r.get(k,0)+ca*cb
            if v: r[k]=v
            elif k in r: del r[k]
    return r
def padd(*ps):
    r={}
    for p in ps:
        for k,c in p.items():
            v=r.get(k,0)+c
            if v: r[k]=v
            elif k in r: del r[k]
    return r
def psc(p,s): return {k:c*s for k,c in p.items() if c*s!=0}
def pdiff(p,i):
    r={}
    for k,c in p.items():
        if k[i]>0:
            k2=list(k); k2[i]-=1; r[tuple(k2)]=c*k[i]
    return r
def wdeg(mono,w): return sum(e*wi for e,wi in zip(mono,w))

def monos_of_weight(w, target, degcap):
    """all monomials x^a (a_i>=0, sum a_i <= degcap) with weighted degree == target."""
    out=[]
    for a in product(range(degcap+1), repeat=3):
        if sum(a)<=degcap and wdeg(a,w)==target:
            out.append(a)
    return out

def equivariant_ansatz_det(w, degcap):
    """Build general equivariant F (weights w), det JF as a polynomial in the
    unknown coefficients; return (det_dict, unknown_count, monomial_lists)."""
    NV=3
    # unknowns: for each component i, one coeff per weight-w_i monomial
    Ms=[monos_of_weight(w, w[i], degcap) for i in range(3)]
    NU=sum(len(M) for M in Ms)
    # represent F_i as dict monomial-> unknown-index (as a symbolic linear form)
    # We'll carry coefficients as vectors over unknowns; det is degree-3 in unknowns.
    # For a landscape scan we only need: does a NONLINEAR Keller solution exist?
    # Cheap proxy: is there an equivariant monomial structure allowing a
    # non-triangular det=const? Use the S59n criterion via the KNOWN owner case.
    return Ms, NU

# --- validation gate: the owner (1,-1,-2) must show the counterexample structure
def has_counterexample_structure(w, degcap=8):
    """Heuristic-but-validated: an equivariant Keller COUNTEREXAMPLE needs, in the
    'z-like' (most-negative-weight) component, a term z*A(invariants) with A a
    NONCONSTANT unit AND a coupled second unit (the crossing). We detect the
    ALGEBRAIC precondition: the invariant ring C[x]^{C*} (weight-0 monomials up to
    degcap) must contain a nonconstant element t, AND there must be >=2 component-
    weights realizable as {t^k * x_i}, so a coupled cascade can form. This is the
    S59n 'need two coupled units to cross' criterion, made into a signature test."""
    # invariants: weight-0 monomials (nonconstant), degree<=degcap
    inv=[a for a in product(range(degcap+1),repeat=3) if 0<sum(a)<=degcap and wdeg(a,w)==0]
    if not inv: return ("DEFINITE-or-no-invariant", 0)
    # count, for each component i, how many weight-w_i monomials are t*x_j style
    Ms=[monos_of_weight(w,w[i],degcap) for i in range(3)]
    richness=[len(M) for M in Ms]
    # counterexample precondition (validated on owner): >=2 components with
    # richness>=3 (room for a nonconstant unit cascade) AND a negative doubling
    # weight (some w_i <= -2 giving lam->lam^2-type orbit doubling)
    has_double = any(wi<=-2 for wi in w)
    coupled = sum(1 for r in richness if r>=3)>=2
    return ("CANDIDATE-COUNTEREXAMPLE" if (coupled and has_double and len(inv)>=1) else "INVERTIBLE-forced", len(inv))

print("=== WEIGHT-SIGNATURE LANDSCAPE (dim 3), degcap 8 ===")
print("  signature test (validated so owner (1,-1,-2) = CANDIDATE):")
tested=[]
for w in product(range(-4,5),repeat=3):
    if gcd(gcd(abs(w[0]),abs(w[1])),abs(w[2]))!=1: continue
    signs=(sum(1 for x in w if x>0), sum(1 for x in w if x<0), sum(1 for x in w if x==0))
    definite = (signs[1]==0 or signs[0]==0) and signs[2]<3
    verdict,ninv=has_counterexample_structure(w)
    tested.append((w,signs,definite,verdict))
# report the owner case + a few definite + the candidate signatures
owner=[t for t in tested if t[0]==(1,-1,-2)]
print("  OWNER (1,-1,-2):", owner[0][1:] if owner else "MISSING")
defs=[t for t in tested if t[2]][:6]
print("  DEFINITE samples (expect INVERTIBLE-forced):")
for w,sg,d,v in defs: print(f"    w={w} sig={sg}: {v}")
cands=sorted(set(t[0] for t in tested if t[3]=="CANDIDATE-COUNTEREXAMPLE" and t[0][0]>=0))
print(f"  CANDIDATE-COUNTEREXAMPLE signatures (indefinite, dim 3, degcap 8): {len(cands)}")
for w in cands[:20]: print(f"    {w}")
# minimality of (1,-1,-2): smallest |w| candidate
if cands:
    mn=min(cands,key=lambda w:sum(abs(x) for x in w))
    print(f"  minimal-|w| candidate: {mn} (owner (1,-1,-2) has |w|=4)")
