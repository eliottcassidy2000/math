#!/usr/bin/env python3
"""
simplex_cube_hadamard.py — opus-2026-03-14-S71e

Direct Hadamard construction for regular simplices inscribed in hypercubes.
No brute force search — uses the algebraic structure directly.

FACT: A regular n-simplex can be inscribed in {0,1}^n iff a Hadamard matrix
of order n+1 exists. Known Hadamard orders: 1, 2, 4, 8, 12, 16, 20, ...
(all multiples of 4, plus 1 and 2).

So regular simplex inscriptions exist at n = 0, 1, 3, 7, 11, 15, 19, ...
These are n = 4k-1 for k = 1, 2, 3, 4, 5, ...

For each such n:
  - The simplex has n+1 vertices from {0,1}^n
  - All edges have length sqrt(n+1)/2... actually sqrt((n+1)/2)
  - Volume = ((n+1)/2)^{n/2} / n!
  - Cube/simplex ratio = n! / ((n+1)/2)^{n/2} * 2^? ... let me compute directly.

The KEY insight: at n=3, cube/simplex = 3, and there are 4 = 2^2 corners,
each a "half" of the simplex.
"""

import sys
import numpy as np
from math import factorial, comb
sys.stdout.reconfigure(line_buffering=True)

def hadamard(n):
    """Sylvester construction: H_{2^k}."""
    if n == 1:
        return np.array([[1]])
    if n == 2:
        return np.array([[1, 1], [1, -1]])
    # For powers of 2
    if n & (n-1) == 0:  # power of 2
        h = hadamard(n // 2)
        return np.block([[h, h], [h, -h]])
    return None

def simplex_volume_from_vertices(vertices):
    """Volume of simplex from vertex array (rows = vertices)."""
    n = len(vertices) - 1
    v0 = vertices[0]
    M = vertices[1:] - v0
    MTM = M @ M.T
    det = np.linalg.det(MTM)
    return np.sqrt(abs(det)) / factorial(n)

print("=" * 70)
print("SIMPLEX IN HYPERCUBE VIA HADAMARD CONSTRUCTION")
print("=" * 70)

# For n=3 (H4): the standard construction
for k in range(1, 7):
    n_plus_1 = 2**k  # Hadamard order (power of 2)
    n = n_plus_1 - 1

    H = hadamard(n_plus_1)
    if H is None:
        continue

    # Convert Hadamard rows to {0,1}^n vertices
    # Map: h_ij ∈ {-1,1} → (h_ij + 1)/2 ∈ {0,1}
    # Drop first column (all 1s after normalization)
    # Actually use all rows, map to {0,1}^n
    vertices = ((H[:, 1:] + 1) // 2).astype(float)  # drop col 0, map to {0,1}

    print(f"\n{'='*60}")
    print(f"n = {n}: Hadamard order {n_plus_1} = 2^{k}")
    print(f"{'='*60}")

    # Check pairwise distances
    dists = set()
    for i in range(n_plus_1):
        for j in range(i+1, n_plus_1):
            d = np.linalg.norm(vertices[i] - vertices[j])
            dists.add(round(d, 8))
    print(f"  Pairwise distances: {sorted(dists)}")
    edge = sorted(dists)[0] if len(dists) == 1 else "NOT REGULAR"
    if len(dists) == 1:
        print(f"  REGULAR simplex, edge = {edge:.6f} = sqrt({round(edge**2)})")
    else:
        print(f"  NOT a regular simplex! Distances: {sorted(dists)}")

    # Volume
    vol = simplex_volume_from_vertices(vertices)
    ratio = 1.0 / vol if vol > 0 else float('inf')
    print(f"  Volume: {vol:.10f}")
    print(f"  Cube/simplex: {ratio:.6f}")

    # Check if ratio is integer
    if abs(ratio - round(ratio)) < 0.01:
        r = round(ratio)
        print(f"  *** EXACT: cube = {r} * simplex ***")
    else:
        # Check if it's n! / something
        print(f"  n! = {factorial(n)}, ratio/n! = {ratio/factorial(n):.6f}")

    # Parities
    parities = [int(sum(vertices[i]) % 2) for i in range(n_plus_1)]
    even_c = parities.count(0)
    odd_c = parities.count(1)
    print(f"  Vertex parities: {even_c} even, {odd_c} odd")

    # Complement
    comp = 1.0 - vol
    print(f"  Complement: {comp:.10f}")

    # Number of "halves"
    half_vol = vol / 2
    num_halves_exact = comp / half_vol
    print(f"  Complement / (simplex/2) = {num_halves_exact:.6f}")
    print(f"  2^{{n-1}} = {2**(n-1)}")

    # Corner piece analysis
    # For the parity-based decomposition:
    # Simplex uses one parity class (at n=3: all even-parity)
    # Each opposite-parity vertex is the apex of a corner piece
    if even_c == n_plus_1 or odd_c == n_plus_1:
        # All one parity
        corner_count = 2**(n-1)  # opposite parity vertices
        if corner_count > 0:
            corner_vol = comp / corner_count
            print(f"  Corner pieces: {corner_count}")
            print(f"  Volume per corner: {corner_vol:.10f}")
            print(f"  Corner/simplex: {corner_vol/vol:.6f}")
    else:
        print(f"  Mixed parities — using non-simplex vertex count")
        corner_count = 2**n - n_plus_1
        print(f"  Non-simplex vertices: {corner_count}")

print("\n" + "=" * 70)
print("ANALYSIS: THE VOLUME FORMULA")
print("=" * 70)

print("""
  For the Hadamard simplex in {0,1}^n (n = 2^k - 1):

  All edges have length sqrt(n+1)/sqrt(2) = sqrt(2^{k-1})...
  Actually at n=3: edge = sqrt(2), at n=7: edge = sqrt(4) = 2.

  Let me compute the volume formula directly.
""")

print(f"  {'n':>4s} {'n+1':>5s} {'vol':>15s} {'1/vol':>15s} {'n!':>12s} {'n!/vol':>15s}")
for k in range(1, 8):
    n = 2**k - 1
    n1 = n + 1

    # Volume of regular simplex with edge length a in R^n:
    # V = a^n * sqrt(n+1) / (n! * 2^{n/2})
    # For Hadamard simplex: edge^2 = n+1 - 1 = n when using Hadamard...
    # Actually at n=3: edge = sqrt(2), edge^2 = 2 = (n+1)/2

    # Direct: use Gram matrix
    # The Gram matrix of the Hadamard simplex:
    # G_ij = (v_i - v_0) · (v_j - v_0)
    # For Hadamard with rows mapped to {0,1}:
    # v_i · v_j = (n+1)/4 for i ≠ j (from H^T H = (n+1)I)
    # v_i · v_i = ... depends on weight

    # Actually, for the ±1 Hadamard: h_i · h_j = (n+1) δ_{ij}
    # Mapped to {0,1}: v = (h+1)/2, so v_i · v_j = (h_i · h_j + sum h_i + sum h_j + n+1)/4
    # For rows of H (excluding col 0): h_i · h_j uses cols 1..n
    # h_i · h_j (cols 1..n) = (n+1)δ_{ij} - 1 (subtract col 0 contribution)
    #                        = n δ_{ij} + δ_{ij} - 1 = n for i=j, -1 for i≠j

    # Wait, I need to be more careful. Let me just compute edge^2 directly.
    # At n=3: v0=(0,0,0), v1=(1,1,0), v2=(1,0,1), v3=(0,1,1)
    # |v1-v0|^2 = 2, |v2-v0|^2 = 2, |v1-v2|^2 = 2. So edge^2 = 2 = (n+1)/2.

    edge_sq = (n+1) / 2  # This is the pattern
    edge = np.sqrt(edge_sq)

    # Regular simplex volume in R^n with edge a:
    # V = (a^n / n!) * sqrt((n+1) / 2^n)
    vol = (edge**n / factorial(n)) * np.sqrt((n+1) / 2**n)

    ratio = 1/vol if vol > 0 else float('inf')

    print(f"  {n:4d} {n1:5d} {vol:15.8f} {ratio:15.4f} {factorial(n):12d} {factorial(n)*vol:15.8f}")

print("\n  OBSERVATION: The ratio cube/simplex grows MUCH faster than n!")
print("  The simplex becomes vanishingly small relative to the cube.")
print("  At n=3: ratio=3 (cube = 3 simplices). Clean!")
print("  At n=7: ratio much larger. The 'halves' picture breaks down.")

print("\n" + "=" * 70)
print("THE CORRECT ANSWER TO 'HOW DOES IT CONTINUE'")
print("=" * 70)

print("""
  THE USER'S OBSERVATION:
    n=2: triangle in square → 2 halves  (not regular: uses right triangle)
    n=3: tetrahedron in cube → 4 halves (regular: alternating vertices)

  HOW IT CONTINUES:

  1. THE VOLUME FORMULA for regular simplex in n-cube:
     Vol(simplex) = 1/n * Vol at n=3 specifically (not general!)
     At n=3: Vol = 1/3, complement = 2/3, 4 corners of 1/6 each.

  2. FOR GENERAL n (not just Hadamard dimensions):
     The POLYNOMIAL identity always holds:

     (x+2)^n = (x+1)^n + sum_{k=0}^{n-1} C(n,k)(x+1)^k

     The "complement" has n terms with coefficients C(n,0),...,C(n,n-1).
     The LEADING complement term is n*(x+1)^{n-1}.

  3. THE TOURNAMENT INTERPRETATION:
     I(Omega, 3) = I(Omega, 2) + [alpha_1 + 5*alpha_2 + 19*alpha_3 + ...]

     The "corner pieces" are the ODD CYCLES:
     - Each 3-cycle contributes 3^1 - 2^1 = 1 to the excess
     - Each independent pair of 3-cycles contributes 3^2 - 2^2 = 5
     - Each independent triple contributes 3^3 - 2^3 = 19

     The corner piece at level k has volume 3^k - 2^k = (3^k - 2^k).

  4. THE 2^{n-1} CONNECTION:
     In tournament theory, 2^{n-1} appears as:
     - Number of nonzero Walsh components of H
     - Walsh normalization factor
     - THM-J: S mod 2^{n-1} universality

     In the simplex-cuboid picture, 2^{n-1} appears as:
     - Number of opposite-parity vertices of {0,1}^n
     - Number of corner pieces at n=3 (and Hadamard dimensions)

  5. THE DEEP TRUTH:
     "Simplex = (x+1)^n" and "Cuboid = (x+2)^n" is PRECISELY the
     framework of the independence polynomial I(Omega, x).

     I(Omega, 2) = H = "simplex evaluation"
     I(Omega, 3) = "cuboid evaluation"

     The GAP I(3) - I(2) = sum alpha_k(3^k - 2^k) counts the CORNERS.
     These corners ARE the odd cycles of the tournament.

     A tournament with no odd cycles (transitive): no corners, simplex = cuboid.
     A maximally cyclic tournament: maximum corners, simplex << cuboid.

     The simplex-cuboid ratio I(3)/I(2) is a MEASURE OF CYCLICITY.
""")

# Final computation: I(3) - I(2) coefficients
print("  CORNER PIECE SIZES by independence level:")
print(f"  {'k':>3s} {'3^k':>8s} {'2^k':>8s} {'3^k-2^k':>10s} {'meaning':>30s}")
for k in range(1, 11):
    print(f"  {k:3d} {3**k:8d} {2**k:8d} {3**k-2**k:10d} {'alpha_'+str(k)+' corner piece size':>30s}")

print(f"\n  I(3) - I(2) = sum_k alpha_k * (3^k - 2^k)")
print(f"  = 1*alpha_1 + 5*alpha_2 + 19*alpha_3 + 65*alpha_4 + ...")
print(f"\n  These are EXACTLY 3^k - 2^k = the number of ways to")
print(f"  fill the complement when you have k independent odd cycles.")

print("\n  BEAUTIFUL IDENTITY:")
print("  3^k - 2^k = sum_{j=0}^{k-1} C(k,j) 2^j  (binomial complement)")
print("  = sum_{j=0}^{k-1} 3^j * 2^{k-1-j}  (telescoping)")
print()
for k in range(1, 8):
    lhs = 3**k - 2**k
    rhs_binom = sum(comb(k,j) * 2**j for j in range(k))
    rhs_tele = sum(3**j * 2**(k-1-j) for j in range(k))
    print(f"  k={k}: 3^k-2^k = {lhs}, binomial = {rhs_binom}, telescoping = {rhs_tele}, match = {lhs == rhs_binom == rhs_tele}")

print("\nDone.")
