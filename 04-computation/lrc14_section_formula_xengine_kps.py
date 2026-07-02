#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYP-3956: the SECTION FORMULA + the x-side gap-surplus engine. kind-pasteur-2026-07-01-S33.

RESTATEMENT (Lean-facing). THM-599's object A(U) = E_c[L^c(U)] satisfies, by doing the c-integral
in CLOSED FORM per x-section (the section {c : ||c - u_i x|| >= h for all i} has measure exactly
the gap surplus of the phase multiset):

    (SF)   A(U) = Int_0^1 F_U(x) dx,     F_U(x) = Sum over circular gaps g of {frac(u_i x)} of (g - w)^+.

All d-fold band overlaps likewise become ONE-CIRCLE integrals of explicit piecewise-linear
integrands (pair: Int (w - ||m x||)^+ dx = w^2; AP-triple: Int (w - 2||m x||)^+ dx = w^2/2 = 2h^2).
No subtorus Haar, no Pontryagin duality: the Lean inputs are Fubini on T^2, the AddCircle ball
volume, and the measure-preservation of x -> m*x — all in mathlib.

ENUMERATION. F_U is piecewise-AFFINE in x. Its breakpoints occur only where the phase ORDER changes
(u_i x = u_j x mod 1: x = j/m, m = u_i - u_j) or a gap crosses w (x = (j +- w)/m = (7j +- 1)/(7m)).
Hence ALL breakpoints lie in (1/M) Z with

    M = 7 * lcm{ |u_i - u_j| }        (w = 1/7; h = 1/14 enters only through w = 2h).

Between consecutive multiples of 1/M the integrand is affine, so the EXACT UNIFORM-GRID TRAPEZOID

    (GT)   A(U) = (1/M) * Sum_{r=0}^{M-1} [ F_U(r/M + 0) + F_U((r+1)/M - 0) ] / 2
                = (1/M) * Sum_{r=0}^{M-1} F_U( (2r+1) / (2M) )        (midpoint form)

is EXACT — a finite sum of rationals. In Lean this is a Finset.sum of an explicit function; the
only analytic lemma needed (once, generic) is "integral of a piecewise-affine function = its
midpoint/trapezoid sum on a refining grid".

TESTS:
 X1: closed-form section identities — pair w^2, AP-triple 2h^2 — by the x-engine (exact rationals).
 X2: the x-engine REPRODUCES the entire S32 c-engine ledger k = 2..13 EXACTLY (rational equality),
     with timing comparison (breakpoints over differences only vs Sum 4uu').
 X3: SCALE demo — the exact A(U) census over ALL primitive 7-subsets of {1..13} (1716 patterns
     before primitivization) in one run: min, argmin, full-histogram head. (Infeasible for the
     c-engine in this budget; the x-engine's M = 7*lcm(diffs <= 12) stays small.)
"""
import sys, itertools, time
from fractions import Fraction as Fr
from math import gcd, floor
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

w = Fr(1, 7)   # = 2h, h = 1/14

def F_U(U, x):
    """gap surplus of the phase multiset {u x mod 1} at rational x — exact."""
    ph = sorted((u*x) % 1 for u in U)
    k = len(ph)
    tot = Fr(0)
    for i in range(k):
        g = (ph[(i+1) % k] - ph[i]) % 1
        if g > w: tot += g - w
    return tot

def lcm(a, b): return a*b // gcd(a, b)

def A_xengine(U):
    """exact A(U) by the x-BREAKPOINT sum: breakpoints of F_U lie in
       { j/m } ∪ { (7j±1)/(7m) } over the pairwise differences m — F_U is affine between
       consecutive breakpoints, so Sum length × F(midpoint) is EXACT.
       (The uniform 1/M grid (GT) is the Lean-facing STATEMENT; this is the same identity
       evaluated on the coarsest complete breakpoint set.)"""
    U = sorted(set(U))
    bps = {Fr(0), Fr(1)}
    for i in range(len(U)):
        for j in range(i+1, len(U)):
            m = U[j] - U[i]
            for jj in range(m):
                bps.add(Fr(jj, m))
                bps.add(Fr(7*jj + 1, 7*m))
                bps.add(Fr((7*jj - 1) % (7*m), 7*m))
    bps = sorted(bps)
    tot = Fr(0)
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i+1]
        if hi > lo:
            tot += (hi - lo) * F_U(U, (lo + hi)/2)
    return tot, len(bps)

# ============================================================================================
print("="*100)
print(" X1: CLOSED-FORM SECTION IDENTITIES (exact)")
print("="*100)
# pair: A(0, m aside) — directly integrate the pair-overlap integrand via the engine on {0, m}? The
# pair OVERLAP itself = w^2 is (ii); here we confirm the section formula gives A(pair) = 1 - 2w + w^2.
for U in [(1,2), (3,10), (5,37)]:
    A, M = A_xengine(U)
    print(f"  A{U} = {A}  (= (1-w)^2 = 36/49? {A == Fr(36,49)})   grid M = {M}")
for U in [(1,2,3), (10,24,38)]:
    A, M = A_xengine(U)
    print(f"  A{U} = {A}  (= 61/98? {A == Fr(61,98)})   grid M = {M}")

print("\n" + "="*100)
print(" X2: x-ENGINE REPRODUCES THE S32 LEDGER EXACTLY (rational equality) + timing")
print("="*100)
ledger = { 2: ((5,37), Fr(36,49)), 3: ((10,24,38), Fr(61,98)), 4: ((5,6,7,8), Fr(11,21)),
           5: ((3,5,6,7,9), Fr(127,294)), 6: ((2,4,6,7,8,10), Fr(141,392)),
           7: ((3,5,6,7,8,9,11), Fr(173,588)), 8: ((1,3,5,6,7,8,9,11), Fr(169,686)),
           9: ((1,3,4,5,6,7,8,9,11), Fr(4267,20580)), 10: ((1,3,4,5,6,7,8,9,11,13), Fr(4279,24696)),
           11: ((2,3,4,5,6,7,8,9,10,12,14), Fr(9451,61740)),
           12: ((4,8,9,10,12,13,14,15,16,17,22,24), Fr(4615489,35315280)),
           13: ((3,7,8,9,11,13,15,16,17,19,23,24,27), Fr(144699147,1267426160)) }
allok = True
for k in range(2, 14):
    U, want = ledger[k]
    t0 = time.time()
    A, M = A_xengine(U)
    dt = time.time() - t0
    ok = (A == want)
    allok &= ok
    print(f"  k={k:2d}: x-engine A = {str(A):>24}   == c-engine? {ok}   (M = {M}, {dt:.2f}s)")
print(f"  ALL k=2..13 REPRODUCED EXACTLY: {allok}")

print("\n" + "="*100)
print(" X3: SCALE DEMO — exact census over ALL 7-subsets of {1..13} (primitivized)")
print("="*100)
t0 = time.time()
seen = {}
for C in itertools.combinations(range(1, 14), 7):
    g = 0
    for c in C: g = gcd(g, c)
    U = tuple(sorted(c//g for c in C))
    key = tuple(u - U[0] for u in U)     # translation-invariant canonical pattern
    if key in seen: continue
    A, nb = A_xengine(U)
    seen[key] = A
mins = sorted(set(seen.values()))
print(f"  patterns evaluated (after primitivize + translation dedup): {len(seen)}   ({time.time()-t0:.1f}s)")
print(f"  MIN A = {mins[0]} = {float(mins[0]):.6f}")
argm = [k for k, v in seen.items() if v == mins[0]][:4]
print(f"  argmin difference-patterns (first 4): {argm}")
print(f"  next values: " + ", ".join(f"{str(v)} ({float(v):.4f})" for v in mins[1:4]))
MP = Fr(14249, 252252)
print(f"  min vs witnessMP: x{float(mins[0]/MP):.2f}   ALL >= MP: {mins[0] >= MP}")
print("DONE.")
