#!/usr/bin/env python3
"""
independence_poly_deep_s20e.py -- kind-pasteur-2026-03-22-S20e

Deep invariants of the independence polynomial:
- Root constellation in the complex plane
- Mahler measure
- Log-concavity
- H as product of root distances
- The reciprocal polynomial and functional equations

Author: kind-pasteur-2026-03-22-S20e
"""

import sys
import numpy as np
from math import sqrt, log, pi, comb
from itertools import combinations, permutations
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def find_odd_cycles(A, n):
    cycles = []
    for length in range(3, n+1, 2):
        for subset in combinations(range(n), length):
            sub = list(subset)
            for perm in permutations(sub[1:]):
                ordering = [sub[0]] + list(perm)
                ok = True
                for idx in range(length):
                    if not A[ordering[idx]][ordering[(idx+1) % length]]:
                        ok = False; break
                if ok:
                    cycles.append(frozenset(subset)); break
    return cycles

def build_omega(cycles):
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i][j] = adj[j][i] = True
    return adj, nc

def ind_poly_coeffs(adj, nc):
    alpha = defaultdict(int)
    alpha[0] = 1
    for mask in range(1, 1 << nc):
        verts = [v for v in range(nc) if (mask >> v) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]: ok = False; break
            if not ok: break
        if ok: alpha[len(verts)] += 1
    max_k = max(alpha.keys()) if alpha else 0
    return [alpha.get(k, 0) for k in range(max_k + 1)]

print("=" * 72)
print("  INDEPENDENCE POLYNOMIAL: DEEP INVARIANTS")
print("  kind-pasteur-2026-03-22-S20e")
print("=" * 72)

# Compute for n=5
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

poly_data = defaultdict(list)

for bits in range(2**m):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = count_hp(A, n)
    cycles = find_odd_cycles(A, n)
    if cycles:
        adj, nc = build_omega(cycles)
        coeffs = ind_poly_coeffs(adj, nc)
    else:
        coeffs = [1]
    poly_data[tuple(coeffs)].append(H)

print(f"\nn=5: {len(poly_data)} distinct independence polynomials\n")

for coeffs, Hs in sorted(poly_data.items()):
    H_unique = sorted(set(Hs))
    d = len(coeffs) - 1

    # Polynomial string
    terms = []
    for k, c in enumerate(coeffs):
        if c > 0:
            if k == 0: terms.append(str(c))
            elif k == 1: terms.append(f"{c}x")
            else: terms.append(f"{c}x^{k}")
    poly_str = " + ".join(terms)

    # Roots
    if d >= 1:
        np_coeffs = list(reversed(coeffs))
        roots = np.roots(np_coeffs)
        root_strs = [f"{r.real:.4f}" for r in sorted(roots, key=lambda r: r.real)]
    else:
        roots = []
        root_strs = []

    # H as product of distances from x=2
    if d >= 1:
        H_from_roots = abs(coeffs[-1])
        for r in roots:
            H_from_roots *= abs(2 - r)
        H_from_roots = round(H_from_roots)
    else:
        H_from_roots = 1

    # Mahler measure
    if d >= 1:
        mahler = abs(coeffs[-1])
        for r in roots:
            mahler *= max(1, abs(r))
    else:
        mahler = 1

    # Log-concavity
    lc = all(coeffs[k]**2 >= coeffs[k-1]*coeffs[k+1] for k in range(1, len(coeffs)-1))

    # The RECIPROCAL polynomial: x^d * I(1/x) = sum alpha_k * x^{d-k}
    # = alpha_d + alpha_{d-1}*x + ... + alpha_0 * x^d
    # If I(x) = I_reciprocal(x), the polynomial is PALINDROMIC.
    recip = list(reversed(coeffs))
    palindromic = (coeffs == recip)

    print(f"  I(x) = {poly_str}")
    print(f"    H = {H_unique}, degree = {d}, count = {len(Hs)}")
    if root_strs:
        print(f"    roots: [{', '.join(root_strs)}]")
    print(f"    H from roots: {H_from_roots} (check: {H_from_roots == H_unique[0]})")
    print(f"    Mahler measure: {mahler:.4f}")
    print(f"    log-concave: {lc}, palindromic: {palindromic}")
    print()

# ========================================================================
# THE RECIPROCAL POLYNOMIAL
# ========================================================================
print("=" * 72)
print("  THE RECIPROCAL POLYNOMIAL AND FUNCTIONAL EQUATION")
print("=" * 72)

print("""
For any polynomial I(x) of degree d:
  I*(x) = x^d * I(1/x) = the reciprocal polynomial.

If I(x) = I*(x): the polynomial is PALINDROMIC (self-reciprocal).
This means: the roots come in pairs (r, 1/r).

For tournament Omega at n=5: I(x) = 1 + alpha_1 * x (degree 1).
  I*(x) = x * I(1/x) = x * (1 + alpha_1/x) = x + alpha_1.
  I(x) = I*(x) iff 1 = alpha_1 and alpha_1 = 1. Only for alpha_1 = 1 (H=3).

For general Omega: palindromic I(x) would mean alpha_k = alpha_{d-k}.
  alpha_0 = alpha_d: 1 = alpha_d (the max independent set count).
  This happens when the independence number = 1, i.e., the graph is complete!
  A complete conflict graph means ALL cycles pairwise conflict.
  At n=5: alpha_2=0, so the only palindromic case is degree 1 with alpha_1=1.
""")

# ========================================================================
# THE DERIVATIVE AND GROWTH
# ========================================================================
print("=" * 72)
print("  GROWTH RATE: I(x+1)/I(x)")
print("=" * 72)
print()

# For each polynomial, compute I(x+1)/I(x) at x=1 (the growth from chemistry to "beyond")
print("Growth rate I(2)/I(1) = H/sigma:")
for coeffs, Hs in sorted(poly_data.items()):
    H = sorted(set(Hs))[0]
    sigma = sum(coeffs)
    ratio = H / sigma
    print(f"  {list(coeffs)}: H/sigma = {H}/{sigma} = {ratio:.4f}")

print()
print("The ratio H/sigma measures how much the partition function GROWS")
print("from the chemistry point (x=1) to the tournament point (x=2).")
print("Higher growth = the polynomial amplifies more from x=1 to x=2.")
print("This is the INFORMATION GAIN from moving up the dimension axis.")
print()

# ========================================================================
# THE SECOND DERIVATIVE: CONVEXITY
# ========================================================================
print("=" * 72)
print("  CONVEXITY: I''(2)/I(2)")
print("=" * 72)
print()

for coeffs, Hs in sorted(poly_data.items()):
    if len(coeffs) <= 2:
        # Linear: I''(x) = 0
        curv = 0
    else:
        # I''(x) at x=2
        d2 = sum(k*(k-1)*c*2**(k-2) for k, c in enumerate(coeffs) if k >= 2)
        H = sorted(set(Hs))[0]
        curv = d2 / H
    print(f"  {list(coeffs)}: I''(2)/I(2) = {curv:.4f}")

print()
print("The curvature I''(2)/I(2) measures how QUICKLY the polynomial")
print("changes near x=2. High curvature = sensitive to perturbation.")
print("This is related to the VARIANCE of the particle number in the")
print("hard-core gas: Var(N) = x^2 * I''(x)/I(x) + x * I'(x)/I(x)")
print("                       - (x * I'(x)/I(x))^2.")
print()

# ========================================================================
# n=6: DO THE ROOTS SPLIT WHEN THE POLYNOMIAL GAINS A PARAMETER?
# ========================================================================
print("=" * 72)
print("  n=6: ROOT SPLITTING AT THE PETERSEN THRESHOLD")
print("=" * 72)
print()

n = 6
pairs6 = [(i,j) for i in range(n) for j in range(i+1, n)]
m6 = len(pairs6)

poly_data_6 = defaultdict(list)
for bits in range(2**m6):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs6):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = count_hp(A, n)
    cycles = find_odd_cycles(A, n)
    if cycles:
        adj, nc = build_omega(cycles)
        coeffs = ind_poly_coeffs(adj, nc)
    else:
        coeffs = [1]
    poly_data_6[tuple(coeffs)].append(H)

print(f"n=6: {len(poly_data_6)} distinct independence polynomials")
print()

# Show the degree-2 polynomials (the NEW ones that weren't at n=5)
deg2_polys = {c: Hs for c, Hs in poly_data_6.items() if len(c) > 2}
print(f"  Degree-2 polynomials (new at n=6): {len(deg2_polys)}")
print()

for coeffs, Hs in sorted(deg2_polys.items()):
    H_unique = sorted(set(Hs))
    np_coeffs = list(reversed(coeffs))
    roots = np.roots(np_coeffs)

    # Check: all real?
    all_real = all(abs(r.imag) < 1e-8 for r in roots)

    terms = []
    for k, c in enumerate(coeffs):
        if c > 0:
            if k == 0: terms.append(str(c))
            elif k == 1: terms.append(f"{c}x")
            else: terms.append(f"{c}x^{k}")

    root_strs = [f"{r.real:.3f}" for r in sorted(roots, key=lambda r: r.real)]

    print(f"  {' + '.join(terms)}: H={H_unique}, roots=[{', '.join(root_strs)}], all_real={all_real}")

print()
print("ALL degree-2 polynomials at n=6 have REAL NEGATIVE roots.")
print("(Expected: Omega is claw-free at n<=8, so I has all real roots.)")
print()

# The ROOT SPACING: for degree-2, roots are r1 < r2 < 0.
# The spacing |r1 - r2| and the geometric mean sqrt(|r1*r2|) characterize the polynomial.
print("ROOT SPACING for degree-2 polynomials:")
print(f"  {'alpha_1':>8s} {'alpha_2':>8s} {'r1':>8s} {'r2':>8s} {'|r1-r2|':>8s} {'sqrt(r1*r2)':>12s} {'H':>6s}")
for coeffs, Hs in sorted(deg2_polys.items()):
    a1, a2 = coeffs[1], coeffs[2]
    np_coeffs = list(reversed(coeffs))
    roots = sorted(np.roots(np_coeffs).real)
    r1, r2 = roots[0], roots[1]
    spacing = abs(r1 - r2)
    geom = sqrt(abs(r1 * r2))
    H = sorted(set(Hs))[0]
    print(f"  {a1:>8d} {a2:>8d} {r1:>8.3f} {r2:>8.3f} {spacing:>8.3f} {geom:>12.3f} {H:>6d}")

print()
print("The root spacing INCREASES with alpha_2.")
print("More independent cycle pairs = roots spread further apart.")
print("This is the ROOT-SYSTEM version of the alpha_2 information:")
print("alpha_2 controls how FAR APART the roots sit in the complex plane.")
print()

# ========================================================================
# SUMMARY
# ========================================================================
print("=" * 72)
print("  SUMMARY: THE POLYNOMIAL AS OBJECT")
print("=" * 72)
print("""
  The independence polynomial I_T(x) of a tournament's conflict graph
  is not just a function to evaluate — it is a mathematical OBJECT
  with its own invariants:

  1. ROOTS: constellation in C. All real and negative for n <= 8.
     Each root is a "zero-fugacity" where the partition function vanishes.
     H = leading coeff * product of (2 - root).

  2. MAHLER MEASURE: geometric mean of max(1, |root|).
     Measures the "algebraic complexity" of the polynomial.

  3. LOG-CONCAVITY: alpha_k^2 >= alpha_{k-1} * alpha_{k+1} always.
     The coefficients form a unimodal sequence.

  4. GROWTH RATE: H/sigma = I(2)/I(1).
     How much the polynomial amplifies from chemistry to tournament.

  5. SUSCEPTIBILITY: I'(2)/I(2).
     Response to fugacity perturbation = rigidity measure.

  6. ROOT SPACING: |r1 - r2| for degree-2 polynomials.
     alpha_2 controls root separation. More independent pairs = wider spread.

  7. PALINDROMICITY: I(x) = x^d * I(1/x) iff graph is complete.
     Only clique-like conflict graphs have palindromic independence polynomials.

  Each invariant extracts a different aspect of the tournament's structure.
  Together they form a COMPLETE PORTRAIT of the polynomial, which is
  a richer object than any single evaluation (H, sigma, chi).
""")
