#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_threadA_sawtooth_kpswf5.py  (kind-pasteur 2026-06-21, THREAD A part 3)

We PROVED computationally:
  edge marginal  f_P(a,1) = 2 a (P-a)   (exact, all primes 3..19).
Now find the full f_P(a,b).  The cell-discrepancy is an L1 sum of a sawtooth, so the
natural closed form is a DEDEKIND-SUM / double-sawtooth.  We rederive f_P(a,b) from
first principles and test against the data.

Derivation recap (residue-only):
  rowdef(a,b) = sum_{j=0}^{P-1} | P*N_j - a*q |,
  N_j = #{ z in [q j, q(j+1)) : z mod P < a },  where q is any rep with q mod P = b.
But by the residue-only theorem this depends only on (a mod P, b mod P).  Take q = b
(the smallest positive rep, when b>=1) and p with p mod P = a, p//P large enough.
Actually rowdef only needs a=p%P and q (the ACTUAL q, not just q%P), but the VALUE is
constant across the residue class -- so we can take the *smallest* in-window q with
q%P=b and any compatible p.  Cleanest: take the limit/representative directly.

Equivalent clean object (substitute z = q j + s, s in [0,q)):
  Over a full period the count N_j is governed by the fractional parts {a*?}.  Empirically
  f_P(a,1)=2a(P-a) is the variance-like term sum_x B_2-ish.  We TEST the double form:

  CANDIDATE D (double sawtooth / 'twisted' parabola):
     f_P(a,b) = 2 * sum_{k=1}^{P-1} ((k a / P)) * ((k b / P)) * (-P)   ... (Dedekind-type)
  where ((x)) = x-floor(x)-1/2 is the sawtooth.  We'll just compute several natural
  Dedekind-type sums and match.

  Also CANDIDATE M (min/lattice):
     f_P(a,b) relates to # lattice points / the area 2ab(P-a)(P-b)/P plus a defect.
"""
from math import gcd, floor
from fractions import Fraction as Fr


def saw(x):
    """((x)) sawtooth: x - floor(x) - 1/2, with ((integer))=0."""
    fx = x - floor(x)
    if fx == 0:
        return Fr(0)
    return Fr(fx) - Fr(1, 2)


def rowdef_residue(a, b, P):
    """exact rowdef value for residue class (a,b), a,b in 1..P-1 (we use 1..(P-1)/2)."""
    for Qm in range(1, 60):
        q = b + P * Qm
        for Pm in range(1, 60):
            p = a + P * Pm
            if p <= q or p >= P * q:
                continue
            if gcd(p, q) != 1:
                continue
            base = p // P
            rr = p % P
            rd = 0
            for j in range(P):
                Nj = sum(1 for z in range(q * j, q * (j + 1)) if z % P < rr)
                rd += abs(P * Nj - rr * q)
            return rd
    return None


def dedekind_like(a, b, P):
    """A Dedekind-sum-style double sawtooth:  -2P * sum_{k=1}^{P-1} ((k a/P))((k b/P))."""
    s = Fr(0)
    for k in range(1, P):
        s += saw(Fr(k * a, P)) * saw(Fr(k * b, P))
    return -2 * P * s


def dedekind_like2(a, b, P):
    """variant: 4 * sum ((ka/P))((kb/P)) * ... try +2P."""
    s = Fr(0)
    for k in range(1, P):
        s += saw(Fr(k * a, P)) * saw(Fr(k * b, P))
    return s  # raw, we'll inspect the ratio


def main():
    primes = [3, 5, 7, 11, 13, 17, 19]
    print("Build f_P and test Dedekind double-sawtooth closed forms")
    print("=" * 70)
    F = {}
    for P in primes:
        h = (P - 1) // 2
        f = {(a, b): rowdef_residue(a, b, P) for a in range(1, h + 1) for b in range(1, h + 1)}
        F[P] = (f, h)

    # Test the raw Dedekind double sum s(a,b)=sum_k ((ka/P))((kb/P)) and find scalar.
    print("Ratio  f_P(a,b) / [ sum_k ((ka/P))((kb/P)) ]  (look for constant = -2P?):")
    for P in primes:
        f, h = F[P]
        ratios = set()
        bad = False
        for a in range(1, h + 1):
            for b in range(1, h + 1):
                s = dedekind_like2(a, b, P)
                if s == 0:
                    bad = True
                    continue
                ratios.add(Fr(f[(a, b)]) / s)
        print(f"   P={P:2d}: ratios={sorted(ratios)}  (want single value){'  [has zero denom]' if bad else ''}")

    # Direct test: f_P(a,b) == -2P * sum_k ((ka/P))((kb/P)) ?
    print("\nDirect test  f_P(a,b) == -2P * sum_k ((ka/P))((kb/P)):")
    allok = True
    for P in primes:
        f, h = F[P]
        ok = True
        for a in range(1, h + 1):
            for b in range(1, h + 1):
                if f[(a, b)] != dedekind_like(a, b, P):
                    ok = False
        print(f"   P={P:2d}: {'MATCH' if ok else 'NO'}")
        allok = allok and ok
    print(f"   ==> {'UNIFORM DEDEKIND MATCH' if allok else 'Dedekind form fails'}")

    # If it fails, fit f_P over the FULL range a,b in 1..P-1 (not just half) to see the
    # natural symmetric object, then reduce.
    print("\nFull-range (a,b in 1..P-1) Dedekind check for P=7 to see exact relation:")
    P = 7
    for a in range(1, P):
        row = []
        for b in range(1, P):
            ded = dedekind_like(a, b, P)
            row.append(str(ded))
        print(f"   a={a}: " + " ".join(f"{x:>5s}" for x in row))
    print("   (compare to reduced f via ||.|| symmetry)")

    # Alternative closed form: the parabola-product MINUS a correction.
    # f(a,1)=2a(P-a) exact.  Maybe f(a,b)=2*min-based.  Test:
    #   f(a,b) =? 2*( a*(P-b) + b*(P-a) ) ... for a<=b?  check edges: b=1 -> 2(a(P-1)+(P-a))
    #   no.  Test f(a,b) = 2*a*b*(P-a)*(P-b)/P + corr.  Print area term & defect.
    print("\nArea term A=2ab(P-a)(P-b)/P and defect f-A (P=7,11):")
    for P in (7, 11, 13):
        f, h = F[P]
        print(f"  P={P}:")
        for a in range(1, h + 1):
            row = []
            for b in range(1, h + 1):
                A = Fr(2 * a * b * (P - a) * (P - b), P)
                row.append(f[(a, b)] - A)
            print(f"    a={a}: defects = {[str(x) for x in row]}")


if __name__ == "__main__":
    main()
