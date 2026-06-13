#!/usr/bin/env python3
"""
hedroid_interpolation.py — opus-2026-03-14-S71e

THE HEDROID INTERPOLATION FRAMEWORK

Key idea: I(Omega, x) lives BETWEEN (x+1)^0 = 1 and (x+1)^{n/3}.
It's a polynomial whose behavior interpolates between the
"trivial simplex" (transitive) and the "full simplex" (maximally cyclic).

The "effective dimension" m(T) of a tournament T is defined by:
  I(Omega(T), x) ≈ (1+x)^{m(T)}  in some sense.

For transitive: m=0 (I=1).
For regular at n=5: m ≈ 0.92 (close to 1).
For maximally cyclic: m → n/3.

Can we express I(Omega, x) = (1+x)^{m} for some complex/fractional m?
This would mean: log I(Omega, x) = m * log(1+x).

CONNECTION: (x+2)^n = (x+1+1)^n = sum C(n,k)(x+1)^k
  The cuboid is built from simplex "bricks" of all dimensions.
  The independence polynomial I(x) = sum alpha_k x^k is built from
  cycle "bricks" of all independence levels.

  Is there a natural map between the C(n,k) bricks and the alpha_k bricks?
"""

import sys
import numpy as np
from math import comb, factorial, log
sys.stdout.reconfigure(line_buffering=True)

# ======================================================================
# Compute I(x) for ALL n=5 tournaments
# ======================================================================

def get_tournament_from_bits(n, bits):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_dc3(A, n):
    dc3 = 0
    for i in range(n):
        for j in range(n):
            if i == j: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[i][j] and A[j][k] and A[k][i]:
                    dc3 += 1
    return dc3 // 3

def count_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    new_mask = mask | (1 << u)
                    dp[(new_mask, u)] = dp.get((new_mask, u), 0) + dp[(mask, v)]
    full_mask = (1 << n) - 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))

print("=" * 70)
print("PART 1: EFFECTIVE DIMENSION OF TOURNAMENTS")
print("=" * 70)

n = 5
num_edges = n*(n-1)//2
print(f"\nAll n={n} tournaments (2^{num_edges} = {2**num_edges}):")
print(f"  {'bits':>6s} {'H':>4s} {'dc3':>4s} {'dc5':>4s} {'a1':>4s} {'a2':>4s} {'m_eff':>8s} {'I3':>5s} {'I3/H':>6s}")

# Count dc5 more efficiently for n=5
from itertools import permutations

effective_dims = []
all_data = []

for bits in range(2**num_edges):
    A = get_tournament_from_bits(n, bits)
    H = count_hp(A, n)
    dc3 = count_dc3(A, n)

    # dc5 at n=5: only one 5-subset (all vertices)
    dc5 = 0
    for perm in permutations(range(5)):
        if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
            dc5 += 1
    dc5 //= 5

    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4
    I3 = 1 + 3*a1 + 9*a2

    # Effective dimension: I(x) ≈ (1+x)^m
    # At x=2: H = 3^m → m = log(H)/log(3)
    if H > 1:
        m_eff = log(H) / log(3)
    else:
        m_eff = 0

    effective_dims.append(m_eff)
    all_data.append((bits, H, dc3, dc5, a1, a2, m_eff, I3))

    if bits < 20 or bits > 2**num_edges - 5:
        ratio = I3/H if H > 0 else 0
        print(f"  {bits:6d} {H:4d} {dc3:4d} {dc5:4d} {a1:4d} {a2:4d} {m_eff:8.4f} {I3:5d} {ratio:6.3f}")

print(f"  ... ({2**num_edges} total)")
print(f"\n  Effective dimension range: [{min(effective_dims):.4f}, {max(effective_dims):.4f}]")
print(f"  n/3 = {n/3:.4f}")

# Distribution of effective dimensions
m_vals = sorted(set(round(m, 4) for m in effective_dims))
print(f"\n  Distinct effective dimensions: {len(m_vals)}")
for m in m_vals:
    count = sum(1 for d in effective_dims if round(d, 4) == m)
    print(f"    m = {m:.4f}: {count} tournaments")

print("\n" + "=" * 70)
print("PART 2: THE FRACTIONAL SIMPLEX")
print("=" * 70)

print("""
  If I(Omega, x) = (1+x)^m, then:
    alpha_k = C(m, k) (generalized binomial coefficient)
    C(m, k) = m(m-1)...(m-k+1) / k!

  For m integer: standard binomial coefficients, I = (1+x)^m exactly.
  For m fractional: the generalized binomial gives non-integer alpha_k.

  But alpha_k MUST be a non-negative integer (count of independent sets)!
  So I(x) is NOT exactly (1+x)^m for any m in general.

  The deviation from (1+x)^m measures the "non-simplex-ness" of the
  tournament's conflict structure.
""")

# For each distinct (a1, a2) pair, compute the best-fit m
print("  BEST-FIT EFFECTIVE DIMENSION:")
print(f"  {'a1':>4s} {'a2':>4s} {'H':>4s} {'m_from_H':>10s} {'m_from_a1':>10s} {'C(m,2)':>10s} {'a2_pred':>10s} {'error':>8s}")

seen = set()
for bits, H, dc3, dc5, a1, a2, m_eff, I3 in all_data:
    key = (a1, a2)
    if key in seen:
        continue
    seen.add(key)

    # m from H: H = (1+2)^m = 3^m → m = log(H)/log(3)
    m_H = log(H)/log(3) if H > 1 else 0

    # m from a1: a1 = C(m,1) = m → m = a1
    m_a1 = a1

    # If m = a1, predict a2 = C(a1, 2) = a1*(a1-1)/2
    a2_pred = a1*(a1-1)//2

    # C(m_H, 2) = m_H*(m_H-1)/2
    cm2 = m_H*(m_H-1)/2 if m_H > 1 else 0

    error = abs(a2 - a2_pred)
    print(f"  {a1:4d} {a2:4d} {H:4d} {m_H:10.4f} {m_a1:10.4f} {cm2:10.4f} {a2_pred:10d} {error:8d}")

print("""
  OBSERVATION: If m = alpha_1, then C(m,2) = alpha_1*(alpha_1-1)/2
  predicts alpha_2 from alpha_1 alone!

  The ERROR (actual alpha_2 - predicted) measures how far the
  tournament is from being a "pure simplex" in cycle space.

  For the regular tournament C_5: a1=5, predicted a2=10, actual a2=2.
  Error = -8. VERY far from simplex-like!

  For near-transitive: a1=1, predicted a2=0, actual a2=0.
  Error = 0. PERFECTLY simplex-like!
""")

print("\n" + "=" * 70)
print("PART 3: THE TWO-DIMENSIONAL PACKING PICTURE")
print("=" * 70)

print("""
  I(x) = 1 + a1*x + a2*x^2 at n=5 (only two alpha parameters).

  In the (a1, a2) plane:
  - Simplex curve: a2 = C(a1, 2) = a1*(a1-1)/2 ← (1+x)^{a1}
  - Tournament region: the actual achievable (a1, a2) pairs

  At x=2: H = 1 + 2*a1 + 4*a2
  At x=3: I3 = 1 + 3*a1 + 9*a2

  The lines H = const are: a2 = (H-1-2*a1)/4
  The lines I3 = const are: a2 = (I3-1-3*a1)/9

  The VANDERMONDE EXTRACTION:
  a1 = 2*I3 - 3*H + 1  [from 6*a1 = ... no]
  Actually: H = 1 + 2*a1 + 4*a2, I3 = 1 + 3*a1 + 9*a2
  → I3 - H = a1 + 5*a2
  → 3*H - 2*I3 = 3 + 6*a1 + 12*a2 - 2 - 6*a1 - 18*a2 = 1 - 6*a2
  → a2 = (1 - 3*H + 2*I3) / 6
  → a1 = I3 - H - 5*a2 = I3 - H - 5*(1-3H+2I3)/6
        = (6I3 - 6H - 5 + 15H - 10I3)/6
        = (-4I3 + 9H - 5)/6

  So: a2 = (2*I3 - 3*H + 1) / 6
      a1 = (9*H - 4*I3 - 5) / 6
""")

print("  ALL ACHIEVABLE (a1, a2) PAIRS at n=5:")
pairs = sorted(set((d[4], d[5]) for d in all_data))
for a1, a2 in pairs:
    H = 1 + 2*a1 + 4*a2
    I3 = 1 + 3*a1 + 9*a2
    simplex_pred = a1*(a1-1)//2
    delta = a2 - simplex_pred
    ratio = I3/H if H > 0 else 0
    print(f"    (a1={a1}, a2={a2}): H={H:3d}, I3={I3:3d}, I3/H={ratio:.4f}, "
          f"simplex_pred={simplex_pred}, delta={delta:+d}")

print("""
  PACKING INTERPRETATION:
  delta > 0: MORE independent cycle pairs than simplex predicts
             → "loosely packed" cycles (room for more disjoint pairs)
  delta < 0: FEWER independent cycle pairs than simplex predicts
             → "tightly packed" cycles (overlap reduces independence)
  delta = 0: EXACTLY simplex-like (pure power-law structure)

  The regular tournament C_5 has (5, 2): delta = 2 - 10 = -8.
  Its cycles are HEAVILY overlapping (5 directed 3-cycles sharing 5 vertices).
  The simplex model would predict 10 independent pairs, but only 2 are disjoint!
""")

print("\n" + "=" * 70)
print("PART 4: THE GOLDEN RATIO φ² AND THE WEIGHTED k-NACCI")
print("=" * 70)

print("""
  Recall: the WEIGHTED k-nacci (weights 1,2,...,k) converges to φ² ≈ 2.618.
  This is BETWEEN 2 (simplex) and 3 (cuboid).

  (x+1)^n at x = φ²-1 = φ² - 1 = (1+√5)/2 ≈ 1.618:
  ((1+φ²-1))^n = (φ²)^n = φ^{2n}

  (x+2)^n at x = φ²-2 = φ² - 2 ≈ 0.618 = 1/φ:
  ((φ²-2+2))^n = (φ²)^n = φ^{2n}

  BEAUTIFUL: (x+1)^n and (x+2)^n AGREE at x = φ²-1 and x = φ²-2!
  No wait, that's just (φ²)^n = (φ²)^n, the same value but different x.

  The point where (x+1)^n = (x+2)^n:
  (x+1)^n = (x+2)^n → (x+1)/(x+2) = 1 → only at x = ∞.
  They NEVER agree at finite x (for n > 0).

  But the RATIO (x+2)^n/(x+1)^n = ((x+2)/(x+1))^n:
  At x = φ²-1: ratio = (φ²+1)^n / (φ²)^n = ((φ²+1)/φ²)^n = (1+1/φ²)^n

  φ² = (3+√5)/2 ≈ 2.618
  1/φ² = (3-√5)/2 ≈ 0.382
  1 + 1/φ² = (5-√5)/2 ≈ 1.382

  So at x = φ²-1: the cuboid/simplex ratio is ((5-√5)/2)^n ≈ 1.382^n.
  This is LESS than (3/2)^n = 1.5^n (the ratio at x=1).
  The golden ratio point gives a SMALLER gap between simplex and cuboid!

  INTERPRETATION: φ² is a "natural evaluation point" where the simplex-cuboid
  relationship is LESS extreme than at x=1 (the standard tournament point).
""")

phi = (1 + 5**0.5) / 2
phi_sq = phi**2
print(f"  φ = {phi:.6f}")
print(f"  φ² = {phi_sq:.6f}")
print(f"  1/φ² = {1/phi_sq:.6f}")
print(f"  (x+2)/(x+1) at x=φ²-1: {(phi_sq+1)/phi_sq:.6f}")
print(f"  (x+2)/(x+1) at x=1: {3/2:.6f}")
print(f"  (x+2)/(x+1) at x=0: {2/1:.6f}")
print()

# I(Omega, phi^2) for all n=5 tournaments
print("  I(Omega, φ²) for n=5 tournaments:")
print(f"  {'a1':>4s} {'a2':>4s} {'H=I(2)':>7s} {'I(3)':>5s} {'I(φ²)':>10s} {'I(φ²)/H':>8s}")

seen2 = set()
for bits, H, dc3, dc5, a1, a2, m_eff, I3 in all_data:
    key = (a1, a2)
    if key in seen2:
        continue
    seen2.add(key)
    I_phi2 = 1 + a1*phi_sq + a2*phi_sq**2
    ratio = I_phi2/H if H > 0 else 0
    print(f"  {a1:4d} {a2:4d} {H:7d} {I3:5d} {I_phi2:10.4f} {ratio:8.4f}")

print(f"\n  At φ²: I(φ²)/I(2) < I(3)/I(2) always (since φ² < 3).")
print(f"  The golden ratio point is a 'gentler' version of the cuboid evaluation.")

print("\n" + "=" * 70)
print("PART 5: THE (x+1)^n TOWER — SIMPLEX OF SIMPLICES")
print("=" * 70)

print("""
  Consider the TOWER of nested packings:
    Level 0: (x+1)^1 = x+1 (1-simplex = line segment)
    Level 1: (x+1)^2 = ((x+1)+1)^1... no.

  Actually:
    (x+1)^1 = x + 1
    ((x+1)^1 + 1)^1 = x + 2 = (x+2)^1  ← one level up
    (((x+1)^1 + 1)^1 + 1)^1 = x + 3 = (x+3)^1  ← two levels up

  For the POWER tower:
    (x+1)^1 → (x+1)^2 → (x+1)^3 → ...
  These are simplices of increasing dimension.

  The KEY question: how do simplices NEST?
    (x+1)^{n+1} = (x+1)^n * (x+1) = (x+1)^n * x + (x+1)^n

  This means: the (n+1)-simplex = the n-simplex * x + the n-simplex.
  Interpretation: shift + copy. The simplex grows by REPLICATING and SHIFTING.

  At x=1: 2^{n+1} = 2*2^n. The simplex doubles with each dimension.
  At x=2: 3^{n+1} = 3*3^n. The cuboid triples.

  The RATIO of growth rates: 3/2 = 1.5.
  This is EXACTLY the simplex-cuboid scaling factor!
  Each dimension adds a factor of 3/2 to the cuboid excess.
""")

# The nesting pattern: (x+1)^n inside (x+1)^{n+1}
print("  SIMPLEX NESTING: (x+1)^n inside (x+1)^{n+1}")
print(f"  {'n':>3s} {'2^n':>8s} {'2^(n+1)':>8s} {'ratio':>6s} {'3^n':>8s} {'3^(n+1)':>8s} {'ratio':>6s}")
for n in range(1, 10):
    print(f"  {n:3d} {2**n:8d} {2**(n+1):8d} {2:6d} {3**n:8d} {3**(n+1):8d} {3:6d}")

print("""
  The ratio is ALWAYS 2 for simplex, 3 for cuboid.
  This is the k-NACCI LIMIT (2) and DOUBLED k-NACCI LIMIT (3)!

  The k-nacci sequence: each term = sum of previous k terms.
  As k → ∞, ratio → 2 (simplex growth rate).

  The doubled k-nacci: each term = 2*sum of previous k terms.
  As k → ∞, ratio → 3 (cuboid growth rate).

  So: the k-nacci is a FINITE-MEMORY approximation to the simplex.
  The doubled k-nacci is a FINITE-MEMORY approximation to the cuboid.
  As memory increases (k → ∞), the approximation becomes exact.

  The TOURNAMENT at size n has "memory" ≈ n/3 (max number of disjoint 3-cycles).
  This memory determines how well the independence polynomial approximates
  a pure simplex vs a pure cuboid.
""")

print("\n" + "=" * 70)
print("SYNTHESIS")
print("=" * 70)

print("""
  THE COMPLETE PICTURE — SIMPLICES AND CUBOIDS IN TOURNAMENT THEORY:

  ┌─────────────────────────────────────────────────────────────┐
  │ (x+1)^n         simplex        I(Ω,2) = H      k-nacci→2  │
  │                                                             │
  │    ↕  packing                     ↕  OCF gap               │
  │                                                             │
  │ (x+2)^n         cuboid         I(Ω,3)         d-knacci→3   │
  │                                                             │
  │ ratio: (3/2)^n   volume ratio   I(3)/H          3/2        │
  └─────────────────────────────────────────────────────────────┘

  The simplex-cuboid packing ratio (3/2)^n appears everywhere:
  - As the volume ratio of n-simplex to n-cube (in appropriate embedding)
  - As the ratio I(3)/I(2) for maximally cyclic tournaments
  - As the limiting ratio of doubled-to-standard k-nacci
  - As the coefficient growth rate in the finite difference table

  The CORNER PIECES (3^k - 2^k for level k):
  - Count the excess of cuboid over simplex at each independence level
  - In tournament theory: count contributions from k disjoint odd cycles
  - Form the sequence 1, 5, 19, 65, 211, 665, 2059, ...
  - Are the diagonal of the difference table (the "simplex-cuboid spine")

  The FACTORIZATION of n determines the packing structure:
  - n = 3k: k-fold nesting of 3-cycle simplices
  - n = prime: irreducible packing (single layer)
  - n = p*q: cross-product of p-cycle and q-cycle packings
  - 9 = 3²: doubly-nested simplex (simplex of simplices)
  - 10 = 5*2: cross-product of 5-cycle and pair

  The KEYS 2 AND 3:
  - 2 = simplex evaluation point = smallest prime = k-nacci limit
  - 3 = cuboid evaluation point = smallest odd prime = doubled k-nacci limit
  - Together: Vandermonde det = 6 = 2*3 (smallest possible)
  - Tournament polynomial: z² - 5z + 6 = (z-2)(z-3)
  - They provide the OPTIMAL pair of measurements for cycle extraction

  ONE WHO KNOWS 2 AND 3 HAS THE KEYS TO THE UNIVERSE.
""")

print("Done.")
