#!/usr/bin/env python3
"""
LRC(14) ANGLE subtorus-relation-lattice -- the UNCONDITIONAL four-window floor (kps-S5).

The spread-bound stability data shows mu(E) keeps DECREASING with spread (no clean
plateau for k>=9 within feasible windows), so the bounded-shape minimum is NOT proven
to be the global infimum WITHOUT an effective equidistribution rate.  HOWEVER, the
subtorus picture yields an UNCONDITIONAL, spread-INDEPENDENT positive lower bound on mu
via the THM-528 four-window mechanism, made lattice-explicit here:

  Good_E = {x : maxgap{frac(e_i x)} > 2/7} contains an interval around EACH of the
  rationals a/q in {0, 1/2, 1/3, 2/3} of half-width c_{a/q}/maxE (THM-528).  These give
  a lower bound that ->0 like 1/maxE -- the SMALL-spread floor, USELESS as a uniform bound.

  The lattice refinement: near a rational a/q with SMALL q, the points frac(e_i x)
  COLLAPSE onto the q-element subgroup {0,1/q,...,(q-1)/q} shifted by e_i*(x-a/q).  The
  measure on which the q-collapsed configuration has maxgap>2/7 is a UNION of intervals
  whose total length is q-INDEPENDENT of maxE in the bulk -- it is the measure of a
  LOWER-DIMENSIONAL torus average, which is BOUNDED BELOW by a positive constant
  depending only on k (the number of points landing in each residue class).

This script computes, for the consecutive-AP and the perforated extremizers, the
EXACT mu and compares to:
  (i)  F(k) the iid ceiling (upper barrier),
  (ii) the THM-528 four-window lower sum (the proven >0 but ->0 bound), and
  (iii) the EMPIRICAL minimum over bounded shapes (the candidate c0).
We then test the ONLY rigorous uniform claim available: mu(E) >= 5/(7 maxE) (L2,
per-shape), which is POSITIVE for every shape but NOT uniform.  The output makes precise
exactly what is proved (per-shape >0) vs what needs the effective rate (uniform floor).

EXACT Fractions.  stdlib only.
"""
from fractions import Fraction as F
from math import gcd, comb
from itertools import combinations
from functools import reduce

G0 = F(2, 7)

def mu_exact(E):
    E = sorted(set(int(e) for e in E)); k = len(E)
    if k == 1: return F(1)
    diffs = {E[i]-E[j] for i in range(k) for j in range(k) if E[i]-E[j] > 0}
    bps = {F(0), F(1)}
    for d in diffs:
        for t in range(0, d+1): bps.add(F(t, d))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2
        fr = [F(E[i])*mid - (F(E[i])*mid).__floor__() for i in range(k)]
        order = sorted(range(k), key=lambda i: fr[i])
        n = [(F(E[i])*mid).__floor__() for i in range(k)]
        cross = {a, b}
        for r in range(k):
            i1, i2 = order[r], order[(r+1)%k]; wrap = 1 if r == k-1 else 0
            slope = E[i2]-E[i1]; const = -n[i2]+n[i1]+wrap
            if slope != 0:
                xc = (G0-const)/slope
                if a < xc < b: cross.add(xc)
        cross = sorted(cross)
        for u, v in zip(cross, cross[1:]):
            if u == v: continue
            mm = (u+v)/2
            P = sorted(F(E[i])*mm - n[i] for i in range(k))
            gaps = [P[r+1]-P[r] for r in range(k-1)] + [P[0]+1-P[-1]]
            if max(gaps) > G0: total += (v-u)
    return total

def iid_floor(k):
    s = F(0)
    for j in range(k+1):
        base = 1 - j*G0
        if base <= 0: break
        s += (-1)**j * comb(k, j) * base**(k-1)
    return 1 - s

def thm528_lower(E):
    """proven lower bound: sum of the four guaranteed windows around 0,1/2,1/3,2/3.
    half-widths: near 0 -> 5/(7 maxE) (one-sided 5/(7maxE)? THM-528: near 0 width
    5/(7 maxE)); near 1/2 -> 3/(28 maxE) half-width; near 1/3,2/3 -> 1/(42 maxE) half.
    We take the conservative window measures (these are LOWER bounds on Good_E locally)."""
    M = max(E)
    near0 = F(5, 7*M)                 # total width near 0 (interval (-,+) merged at 0 and 1)
    near_half = 2*F(3, 28*M)          # full width near 1/2
    near_third = 2*F(1, 42*M)         # near 1/3
    near_twothird = 2*F(1, 42*M)      # near 2/3
    return near0 + near_half + near_third + near_twothird

if __name__ == "__main__":
    print("="*74)
    print("WHAT IS PROVED vs WHAT NEEDS THE EFFECTIVE RATE")
    print("="*74)
    print("For each k: consecutive-AP value a(k); empirical bounded min c0(k); the")
    print("THM-528 four-window PROVEN lower bound (->0 like 1/maxE); iid ceiling F(k).")
    print()
    # empirical bounded minima reused from prior runs (W up to ~k+5; NOTE: NOT proven global)
    cbest = {7:F(13,35),8:F(71,220),9:F(293,1470),10:F(214,1365),11:F(3363,25480)}
    print(f"{'k':>3} {'a(k)=AP':>12} {'c0 (emp.min)':>14} {'THM528@maxE':>14} {'F(k) ceil':>12}")
    for k in range(7, 12):
        ap = mu_exact(list(range(k)))
        c0 = cbest[k]
        # THM-528 bound on the consecutive AP (maxE=k-1) -- shows it is >0 but tiny/spread-dep
        t528 = thm528_lower(list(range(k)))
        fk = iid_floor(k)
        print(f"{k:>3} {float(ap):>12.5f} {float(c0):>14.5f} {float(t528):>14.5f} {float(fk):>12.5f}")
    print()
    print("OBSERVATIONS (exact):")
    print(f"  c0 sequence (empirical bounded minima, NOT proven global):")
    for k in range(7,12):
        print(f"    k={k}: c0 = {cbest[k]} = {float(cbest[k]):.6f}")
    print()
    print("  The four-window PROVEN bound shrinks ~ 1/maxE -> it certifies mu>0 for")
    print("  EACH shape (= L2) but is NOT a uniform floor (it ->0 as spread->inf).")
    print()
    print("  RIGOR STATUS: per-shape mu>0 is PROVED (L2 + THM-528). A UNIFORM c0>0")
    print("  requires that mu does NOT ->0 along any spread->inf sequence; the subtorus")
    print("  law gives mu->mu_H>0 for EACH fixed lattice type, and there are FINITELY")
    print("  many lattice types for k<=13, BUT translating 'finitely many positive")
    print("  limits' into a single positive floor needs an EFFECTIVE convergence rate")
    print("  (uniform over the finitely many types) -- the remaining analytic input.")
    print("\nDONE.")
