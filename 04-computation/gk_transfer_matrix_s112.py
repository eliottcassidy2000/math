#!/usr/bin/env python3
"""
gk_transfer_matrix_s112.py — Transfer matrix proof that g_k is cubic for k>=3
kind-pasteur-2026-03-15-S112

PROVED: E[prod_{j in S} Z_j] = 2^c / (n)_L where c=#components, L=|S|.

So g_k(m) = (1/2) * sum over k-matchings of P_{m+2k-1} of 2^{c(M)}.

A k-matching of path P_N (N vertices, N-1 edges) is a set of k non-adjacent
edges. Two matched edges at positions i and i+2 have connected position sets.

Transfer matrix for sum over matchings weighted by 2^{#components} * x^{#edges}:

States:
  A: not in cluster (last match >= 3 edges ago or never)
  B: just matched (current edge is in matching)
  C: one edge after match (can extend cluster at next position)

Transfer matrix M(x) for processing each edge:
  From A: don't match -> A, match -> B (new cluster, factor 2x)
  From B: can't match (adjacent) -> C
  From C: don't match -> A (cluster ends), match -> B (extend, factor x)

  M = [[1, 2x, 0],
       [0,  0, 1],
       [1,  x, 0]]

The total weight for N-1 edges starting from state A is:
  f(N,x) = sum over final states of [M^{N-1}]_{A,*} applied to start vector [1,0,0]

Coefficient of x^k in f(N,x) is sum over k-matchings of 2^c.
Then g_k(m) = (1/2) * [x^k] f(m+2k, x) where N = m+2k-1+1 = m+2k.
Wait: the path P_{n-1} has n-1 vertices and n-2 edges. N = n-1, edges = n-2.
So we apply M^{n-2} times, and n = m+2k, so n-2 = m+2k-2.

Let me just compute [x^k] of (sum of entries of M^{n-2} applied to [1,0,0])
for various n, extract g_k(m), and verify the polynomial structure.
"""

import numpy as np
from fractions import Fraction
from math import comb

def transfer_matrix_gk(n_max, k_max):
    """Compute g_k(m) = (1/2) * [x^k] f(n,x) for various n.

    f(n,x) = sum of entries of M(x)^{n-2} * [1,0,0]^T
    """
    results = {}  # results[(k, m)] = g_k(m)

    for n in range(3, n_max + 1):
        num_edges = n - 2

        # Track polynomials in x as coefficient vectors
        # State vector: 3 polynomials (one per state), each as list of coefficients
        # poly[i] = coefficient of x^i

        max_deg = min(k_max, num_edges // 2)

        # Initialize: start in state A
        # state_A = [1] (constant 1), state_B = [0], state_C = [0]
        state = [[Fraction(1)], [Fraction(0)], [Fraction(0)]]

        def poly_add(a, b):
            result = [Fraction(0)] * max(len(a), len(b))
            for i, v in enumerate(a):
                result[i] += v
            for i, v in enumerate(b):
                result[i] += v
            return result

        def poly_scale(a, c):
            return [v * c for v in a]

        def poly_shift(a):
            """Multiply polynomial by x (shift coefficients right)."""
            return [Fraction(0)] + list(a)

        for step in range(num_edges):
            A, B, C = state

            # New states after processing this edge:
            # new_A = A (don't match from A) + C (don't match from C, cluster ends)
            new_A = poly_add(A, C)

            # new_B = 2x*A (match from A, new cluster) + x*C (match from C, extend cluster)
            new_B = poly_add(poly_scale(poly_shift(A), Fraction(2)),
                            poly_shift(C))

            # new_C = B (must not match from B)
            new_C = list(B)

            # Truncate to max_deg
            for poly in [new_A, new_B, new_C]:
                while len(poly) > max_deg + 1:
                    poly.pop()

            state = [new_A, new_B, new_C]

        # Total: sum over all states
        total = poly_add(poly_add(state[0], state[1]), state[2])

        for k in range(1, min(max_deg + 1, len(total))):
            coeff = total[k]
            m = n - 2*k
            if m >= 0:
                gk = coeff / 2  # g_k(m) = (1/2) * [x^k]
                results[(k, m)] = gk

    return results

print("="*70)
print("g_k(m) FROM TRANSFER MATRIX")
print("="*70)

results = transfer_matrix_gk(30, 10)

for k in range(1, 11):
    print(f"\nk={k}:")
    vals = [(m, results[(k,m)]) for m in range(0, 25) if (k,m) in results]
    for m, gk in vals:
        print(f"  g_{k}({m}) = {gk}")

# Now verify against known g_k polynomials
print("\n" + "="*70)
print("VERIFY AGAINST KNOWN POLYNOMIALS")
print("="*70)

def g_known(k, m):
    if k == 1: return Fraction(m)
    if k == 2: return Fraction(m * m)
    if k == 3: return Fraction(2*m**3 + m, 3)
    if k == 4: return Fraction(10*m**3 - 33*m**2 + 50*m - 24, 3)
    if k == 5: return Fraction(388*m**3 - 2040*m**2 + 3431*m - 1776, 3)
    return None

for k in range(1, 6):
    print(f"\nk={k}:")
    for m in range(0, 15):
        if (k, m) not in results:
            continue
        tm = results[(k, m)]
        kn = g_known(k, m)
        if kn is not None:
            match = "OK" if tm == kn else f"MISMATCH (tm={tm}, known={kn})"
            print(f"  m={m}: transfer={tm}, known={kn} {match}")

# NOW: fit polynomials to g_k for k=1..10 and determine degree
print("\n" + "="*70)
print("POLYNOMIAL DEGREE ANALYSIS")
print("="*70)

def fit_degree(values):
    """Given [(m, g_k(m))], determine the degree of the polynomial.
    Use forward differences."""
    if not values:
        return -1
    # Need consecutive m values starting from some m0
    ms = [m for m, _ in values]
    vs = [v for _, v in values]

    # Forward differences
    diffs = list(vs)
    degree = 0
    for d in range(1, len(diffs)):
        new_diffs = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
        if all(v == 0 for v in new_diffs):
            break
        degree = d
        diffs = new_diffs

    return degree

for k in range(1, 11):
    vals = [(m, results[(k,m)]) for m in range(0, 25) if (k,m) in results]
    if len(vals) < 3:
        continue

    deg = fit_degree(vals)

    # Also compute forward differences to show structure
    ms = [m for m, _ in vals]
    vs = [v for _, v in vals]

    # Check if degree is exactly 3 for k >= 3
    print(f"\nk={k}: degree = {deg}, {len(vals)} data points")
    print(f"  g_{k}(0..5) = {[str(results.get((k,m), '?')) for m in range(6)]}")

    # Forward differences at order deg
    diffs = list(vs)
    for d in range(deg):
        diffs = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
    if diffs:
        print(f"  {deg}-th differences: {diffs[:6]} (should be constant)")

# Extract the cubic coefficients for k >= 3
print("\n" + "="*70)
print("CUBIC COEFFICIENTS: 3*g_k(m) = a*m^3 + b*m^2 + c*m + d")
print("="*70)

for k in range(3, 11):
    vals = [(m, results[(k,m)]) for m in range(0, 25) if (k,m) in results]
    if len(vals) < 4:
        continue

    # g_k(m) = (a*m^3 + b*m^2 + c*m + d) / 3
    # Use m=0,1,2,3 to solve for a,b,c,d
    g0 = results.get((k, 0), Fraction(0))
    g1 = results.get((k, 1), Fraction(0))
    g2 = results.get((k, 2), Fraction(0))
    g3 = results.get((k, 3), Fraction(0))

    d = 3 * g0
    c = 3 * g1 - d
    # g(2) = (8a + 4b + 2c + d)/3 => 3g(2) = 8a + 4b + 2c + d
    # g(3) = (27a + 9b + 3c + d)/3 => 3g(3) = 27a + 9b + 3c + d
    # From g(2): 8a + 4b = 3g2 - 2c - d
    # From g(3): 27a + 9b = 3g3 - 3c - d

    rhs2 = 3*g2 - 2*c - d
    rhs3 = 3*g3 - 3*c - d

    # 8a + 4b = rhs2
    # 27a + 9b = rhs3
    # 9*rhs2 - 4*rhs3 = 72a + 36b - 108a - 36b = -36a
    a = (4*rhs3 - 9*rhs2) / 36
    b = (rhs2 - 8*a) / 4

    print(f"\nk={k}: 3*g_{k}(m) = {a}*m^3 + {b}*m^2 + {c}*m + {d}")

    # Verify at m=4,5
    for m in [4, 5, 6, 10]:
        pred = (a*m**3 + b*m**2 + c*m + d) / 3
        actual = results.get((k, m), None)
        if actual is not None:
            match = "OK" if pred == actual else f"FAIL (pred={pred})"
            print(f"  m={m}: {match}")

print("\nDone!")
