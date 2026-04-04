#!/usr/bin/env python3
"""
Searching for a closed-form expression for H(x₁,...,xₘ).
opus-2026-04-03-S28

What we know:
  H(0,...,0) = 1  (transitive)
  Linear coeff of tile j = 2^(sⱼ-1)
  Cross-end pairs: c₂ = 2^(s₁+s₂-2)
  Same-end adjacent pairs: c₂ = -2

The linear+cross-end pattern suggests:
  H ≈ prod over some structure of (1 + 2^(s-1) * x_j)

If H were a PRODUCT: H = prod_j (1 + a_j x_j),
then the order-k coefficient of x_{j₁}...x_{jₖ} = prod a_{jᵢ}.
The order-2 coefficient of x_j x_k = a_j * a_k.

For the cross-end pair CF at n=5: a_C * a_F = 2^1 * 2^1 = 4 = c₂. ✓
For the same-end pair AB at n=5: a_A * a_B = 2^3 * 2^2 = 32. But c₂ = -2. ✗

So H is NOT a simple product. But what about a product along "chains"?

The cross-end pairs form CHAINS through the staircase:
  (3,1) → through vertex 3 → (5,3) → through vertex 5 → ...

Each chain is a path from bottom-left to top-right of the staircase.
Could H be a SUM over chains, where each chain contributes a product?

This is reminiscent of a PATH INTEGRAL through the staircase.

Actually, there's a classical interpretation: the Hamiltonian paths of
the tournament correspond to paths in the staircase where you "choose"
which arcs to use. The transfer matrix formulation should give H as a
matrix product.

Let me check: does the transfer matrix factorize along skip layers?
"""

import numpy as np
from collections import defaultdict

def tournament_from_bits(n, bits_list):
    T = [[0]*n for _ in range(n)]
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    for i in range(n, 1, -1):
        T[i-1][i-2] = 1
    for idx, (x, y) in enumerate(tiles):
        if idx < len(bits_list):
            if bits_list[idx] == 0:
                T[x-1][y-1] = 1
            else:
                T[y-1][x-1] = 1
    return T

def count_hp(T, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

# The transfer matrix approach:
# H(T) = sum over permutations sigma of prod_{i=1}^{n-1} T[sigma(i)][sigma(i+1)]
# = 1^T M^{n-1} 1 where M is the adjacency matrix

# For our tiling model, M depends on the tile bits.
# M[i][j] = 1 if arc i→j exists.
# Each tile bit controls one entry of M.

# But M is an n×n matrix, and H = sum_v1 sum_v2≠v1 ... sum_vn≠v1..vn-1 prod T[vi][vi+1]
# This is NOT simply a matrix power.

# The CORRECT transfer matrix uses the DP state:
# dp[S][v] where S is a bitmask of visited vertices and v is the current vertex.
# This is exponential in n.

# However, we can think of H as a PERMANENT-like quantity:
# H = permanent of a certain matrix related to T.

# Actually, H(T) for a tournament on n vertices equals:
#   sum_{sigma in S_n} prod_{i=1}^{n-1} T[sigma(i), sigma(i+1)]
# This is NOT quite a permanent (permanent is prod T[i, sigma(i)]).
# It's the "path permanent" or "rook polynomial on a path."

# Let me investigate the structure more directly.
# For n=5, I have the complete multilinear polynomial (35 terms).
# Can I write it in a compact form?

n = 5
m = 6
tiles = get_tiles(n)

# Compute all H values
H_vals = {}
for bits in range(2**m):
    b = [(bits >> i) & 1 for i in range(m)]
    T = tournament_from_bits(n, b)
    H_vals[bits] = count_hp(T, n)

# Multilinear coefficients
poly = {}
for S in range(2**m):
    order = bin(S).count('1')
    coeff = 0
    for T_bits in range(2**m):
        if T_bits & S != T_bits:
            continue
        diff = bin(S & ~T_bits).count('1')
        coeff += ((-1) ** diff) * H_vals[T_bits]
    if abs(coeff) > 1e-10:
        poly[S] = int(coeff)

print("="*70)
print(f"COMPLETE MULTILINEAR POLYNOMIAL: n={n}, m={m}")
print("="*70)

print(f"\nTiles: A=(5,1) B=(4,1) C=(3,1) D=(5,2) E=(4,2) F=(5,3)")
print(f"Skips: A=4    B=3    C=2    D=3    E=2    F=2")

print(f"\nH = ", end="")
terms = []
for S in sorted(poly, key=lambda s: (bin(s).count('1'), s)):
    coeff = poly[S]
    support = [chr(65+i) for i in range(m) if (S >> i) & 1]
    var = ''.join(support) if support else '1'
    if coeff > 0:
        terms.append(f"{coeff}{var}")
    else:
        terms.append(f"({coeff}){var}")

print(' + '.join(terms[:20]))
if len(terms) > 20:
    print(f"  + ... ({len(terms) - 20} more terms)")

# Total number of terms
print(f"\nTotal terms: {len(poly)}")

# Can we factor this?
# Group terms by which "chain" they belong to.
# Chain 1: A-B-C (bottom row, y=1)
# Chain 2: A-D-F (through 5→2→5→3)... no, chains are harder to define.

# Let me try a different approach: VERTEX-BASED decomposition.
# Each vertex v appears in some tiles. Can we write H in terms of
# "vertex contributions"?

# Vertex 1: tiles A=(5,1), B=(4,1), C=(3,1) — all lower vertex = 1
# Vertex 2: tiles D=(5,2), E=(4,2) — lower vertex = 2
# Vertex 3: tiles C=(3,1), F=(5,3) — vertex 3 in both
# Vertex 4: tiles B=(4,1), E=(4,2) — upper vertex = 4
# Vertex 5: tiles A=(5,1), D=(5,2), F=(5,3) — upper vertex = 5

# The "skip chain" from vertex 1 to vertex 5:
# 1 → (via C) → 3 → (via F) → 5
# This chain uses tiles C and F (cross-end through vertex 3 and 5).
# c₂(C,F) = 2^(2+2-2) = 4. And the product (1+2C)(1+2F) = 1+2C+2F+4CF.
# Compare with actual: H has constant 1, linear C=2, F=2, CF=4. MATCH!

# Check: does (1+2C)(1+2F) give the right contribution for ALL terms involving C and F?
# At n=5, the full polynomial has:
# Terms with C but not F: 2C, -2AC, -2BC, 2CD, 0CE, -2ACE, ...
# So (1+2C)(1+2F) doesn't capture everything — there are cross-terms with other tiles.

# But the PURE C-F interaction is perfectly captured by the product (1+2C)(1+2F).
# This suggests a partial product structure along chains.

# Let me try: the "anti-diagonal chains" through the staircase.
# Each Hamiltonian path through the vertex set can be described by
# a sequence of arcs. The arcs that are "tile arcs" (not base-path arcs)
# form the non-trivial part. These tile arcs form PATHS through the staircase.

# THE KEY INSIGHT: H(T) counts labeled paths in a graph that depends on T.
# The transfer matrix IS the adjacency matrix of T.
# H = sum_sigma prod T[sigma(i), sigma(i+1)] = permanent of a path matrix.

# For the tiling model: T[x-1][y-1] = x_tile for tile (x,y), and
# T[i][i-1] = 1 for base-path arcs. All other arcs are zero (we're
# using the TILING model, not the full tournament).
# Wait no — the full tournament has ALL arcs. T[i][j] = 1 for exactly
# one of {(i,j), (j,i)} for each pair.

# The point is: the tiling bits enter the adjacency matrix LINEARLY.
# H is a permanent-like function of a matrix with linear entries.
# The permanent of a matrix with linear entries is a polynomial of
# degree at most n-1 (number of rows minus 1 for path permanent).

# And the degree we observe is ≤ n-2 or so. This is the "column count" bound.

print(f"\n{'='*70}")
print(f"TESTING PARTIAL PRODUCT ALONG ANTI-DIAGONAL CHAINS")
print(f"{'='*70}")

# For n=5, the anti-diagonal chains are:
# Path 1→3→5 using tiles C=(3,1) and F=(5,3): product (1+2C)(1+2F) = 1+2C+2F+4CF
# Path 1→4→... : not a clean chain
# Path 2→4→5: not tiles...

# Actually the relevant chains are formed by cross-end tile pairs.
# Cross-end (x₁=v=y₂): tile (v,y₁) and tile (x₂,v).
# The "through" path goes: y₁ → v → x₂.

# At n=5:
# Through v=3: C=(3,1) → F=(5,3) → [end at 5]
# Through v=4: B=(4,1) → E=(4,2) → [wait, (4,2) has upper=4]
#   B=(4,1): lower=1, upper=4. E=(4,2): lower=2, upper=4.
#   These share upper vertex 4 — same-end, not cross-end!

# Cross-end chains through middle vertices:
# v=3: need (v, y₁) and (x₂, v) = (3, y₁) and (x₂, 3)
#   (3,1) and (5,3): cross-end through 3 → chain C-F
# v=5: need (5, y) and (x, 5)
#   No tiles with upper=5... wait (5,1), (5,2), (5,3) all have upper=5.
#   And no tiles have lower=5 (since y < n-1 = 4). So no cross-end through 5.

# Interesting: at n=5, the only cross-end chain is C-F through vertex 3.
# The product (1+2C)(1+2F) contributes 4 terms: 1, 2C, 2F, 4CF.
# The full polynomial has 35 terms.

# The other terms come from:
# - Individual tiles: 8A, 4B, 4D (not in any cross-end chain)
# - Same-end interference: -2AB, -2AC, -2AD, -2AF, -2BC, -2BE, -2DE, -2DF
# - Higher-order terms

# So the product structure is PARTIAL — it captures the cross-end cooperation
# but not the full polynomial.

# OBSERVATION: The polynomial can be written as:
# H = (1 + 2C)(1 + 2F)(1 + ...) + correction terms

# What if we write H = (1 + Σ 2^(sⱼ-1) xⱼ) + cross-terms?
# The "naive" model H_naive = 1 + Σ 2^(sⱼ-1) xⱼ is the linearization.
# This gives H_naive(1,...,1) = 1 + 2+4+2+4+2+2 = 17 but actual H(1,...,1) = 9.
# The corrections are LARGE.

# Let me compute the "chain product" for n=5:
# The three vertices 1, 3, 5 form a "chain" through the staircase.
# Tiles: C=(3,1), F=(5,3) form the cross-end chain.
# Tiles A=(5,1) is the apex — a direct jump from 1 to 5.

# Is H related to counting PATHS through the staircase?
# A path from bottom-left (vertex 1) to top-right (vertex n) through
# the staircase, using tiles as steps...

# This is getting into deep territory. Let me just verify a specific
# structural formula.

# FORMULA ATTEMPT for n=5:
# H = 1 + 8A + 4B + 2C + 4D + 2E + 2F
#     - 2(AB + AC + AD + AF + BC + BE + DE + DF)
#     + 4(AE + BD + CF) + 2(BF + CD)
#     - 4(BCD + BDF) - 2(ABD + ABE + ACE + ADE + BDE + AEF)
#     + 2(ACDE + BCDE)

# This doesn't simplify to a product.

# BUT: look at the pattern by which vertices are "active":
# The AE pair: tiles (5,1) and (4,2). These are DISJOINT in vertices.
# c₂(AE) = +4 = 2^(4+2-2) = 2^4 = 16? No, we computed c₂=4 earlier.
# Wait, +4 was the Möbius coefficient in the multilinear polynomial.
# The WALSH coefficient was different.

print(f"\nLinear model: H_lin = 1 + Σ 2^(sⱼ-1) xⱼ")
for bits in [0, 63, 42, 21]:
    b = [(bits >> i) & 1 for i in range(m)]
    H_actual = H_vals[bits]
    H_linear = 1 + sum(2**(tiles[j][0]-tiles[j][1]-1) * b[j] for j in range(m))
    print(f"  bits={bits:06b}: H_actual={H_actual}, H_linear={H_linear}, "
          f"error={H_actual-H_linear}")

# The errors are systematic: the linear model OVERCOUNTS.
# The corrections (order 2+) are NEGATIVE for overlapping tiles.

print(f"\nKey structural observation:")
print(f"  The multilinear polynomial has a 'tree of corrections' structure.")
print(f"  Order 0: 1 (transitive)")
print(f"  Order 1: +2^(s-1) per tile (new HPs from single flips)")
print(f"  Order 2: corrections for overlap/cooperation")
print(f"  Order 3: corrections to corrections")
print(f"  Order 4: highest-order corrections")
print(f"  This is an INCLUSION-EXCLUSION over tile subsets!")
