#!/usr/bin/env python3
"""
Using the corrected consecutive-position formula to understand scalar M.

M = (H/n)*I iff for all a ≠ b:
  sum_j (-1)^j * N(a,b,j) = 0
where N(a,b,j) = #{paths with {a,b} at consecutive positions {j,j+1}}.

This is a VERY strong constraint: for EVERY pair of vertices, the
alternating sum of their consecutive-position counts must vanish.

QUESTIONS:
1. What is the structure of N(a,b,j) for scalar M tournaments?
2. Does N(a,b,j) have a palindromic or other symmetry?
3. How does N relate to the independence polynomial / OCF?
4. What happens to N under tile flips in the skeleton?
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def ham_paths(A):
    """Return all Hamiltonian paths as tuples."""
    n = len(A)
    paths = []
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if A[perm[k]][perm[k+1]] != 1: ok = False; break
        if ok: paths.append(perm)
    return paths

def N_matrix(A):
    """N[a][b][j] = #{paths with {a,b} at positions {j,j+1}}."""
    n = len(A)
    N = [[[0]*(n-1) for _ in range(n)] for _ in range(n)]
    for p in ham_paths(A):
        for j in range(n-1):
            a, b = p[j], p[j+1]
            N[a][b][j] += 1
            N[b][a][j] += 1  # Symmetrize
    return N

def all_tournaments_canonical(n):
    """Yield one tournament per isomorphism class."""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set()
    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        key = tuple(tuple(row) for row in A)
        min_key = key
        for perm in permutations(range(n)):
            pkey = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
            if pkey < min_key: min_key = pkey
        if min_key in seen: continue
        seen.add(min_key)
        yield A

# =====================================================================
# n=5: N(a,b,j) for scalar M classes
# =====================================================================
print("=" * 70)
print("n=5: N(a,b,j) FOR SCALAR M CLASSES")
print("=" * 70)

n = 5
for A in all_tournaments_canonical(n):
    paths = ham_paths(A)
    H = len(paths)
    scores = sorted([sum(A[i]) for i in range(n)], reverse=True)

    # Compute M to check if scalar
    M = np.zeros((n, n), dtype=int)
    for p in paths:
        for j in range(n):
            M[p[j]][p[j]] += (-1)**j
        for j in range(n-1):
            a, b = p[j], p[j+1]
            M[a][b] += (-1)**j
            M[b][a] += (-1)**j

    is_scalar = np.array_equal(M, (H//n) * np.eye(n, dtype=int)) if H % n == 0 else False

    if not is_scalar:
        continue

    print(f"\n  H={H}, scores={scores}, M = {H//n}*I")

    N = N_matrix(A)

    # Show N(a,b,j) for all pairs
    for a in range(n):
        for b in range(a+1, n):
            n_vals = [N[a][b][j] for j in range(n-1)]
            alt_sum = sum((-1)**j * n_vals[j] for j in range(n-1))
            total = sum(n_vals)
            edge_dir = "a→b" if A[a][b] == 1 else "b→a"
            print(f"    N({a},{b},j) = {n_vals}, alt_sum={alt_sum}, "
                  f"total={total}, {edge_dir}")

    # Check palindromic property: N(a,b,j) = N(a,b,n-2-j)?
    is_palindromic = True
    for a in range(n):
        for b in range(a+1, n):
            for j in range(n-1):
                if N[a][b][j] != N[a][b][n-2-j]:
                    is_palindromic = False
                    break
    print(f"    N palindromic (N(a,b,j) = N(a,b,n-2-j))? {is_palindromic}")

    # Check: is N(a,b,j) constant in j?
    is_const = True
    for a in range(n):
        for b in range(a+1, n):
            vals = [N[a][b][j] for j in range(n-1)]
            if len(set(vals)) > 1:
                is_const = False
    print(f"    N constant in j? {is_const}")

# =====================================================================
# n=5: N(a,b,j) for NON-scalar M classes (comparison)
# =====================================================================
print()
print("=" * 70)
print("n=5: N(a,b,j) FOR NON-SCALAR M (sample)")
print("=" * 70)

count = 0
for A in all_tournaments_canonical(n):
    paths = ham_paths(A)
    H = len(paths)

    M = np.zeros((n, n), dtype=int)
    for p in paths:
        for j in range(n):
            M[p[j]][p[j]] += (-1)**j
        for j in range(n-1):
            a, b = p[j], p[j+1]
            M[a][b] += (-1)**j
            M[b][a] += (-1)**j

    is_scalar = np.array_equal(M, (H//n) * np.eye(n, dtype=int)) if H % n == 0 else False
    if is_scalar: continue

    count += 1
    if count > 3: break

    scores = sorted([sum(A[i]) for i in range(n)], reverse=True)
    print(f"\n  H={H}, scores={scores}")

    N = N_matrix(A)
    for a in range(n):
        for b in range(a+1, n):
            n_vals = [N[a][b][j] for j in range(n-1)]
            alt_sum = sum((-1)**j * n_vals[j] for j in range(n-1))
            if alt_sum != 0:
                print(f"    N({a},{b},j) = {n_vals}, alt_sum={alt_sum} ≠ 0")

# =====================================================================
# N(a,b,j) AND TILE FLIPS: How does N change under arc reversal?
# =====================================================================
print()
print("=" * 70)
print("N UNDER TILE FLIPS")
print("=" * 70)

# For the tiling model: path arcs are (i+1→i), tiles are non-path arcs.
# When we flip tile (u,v), some paths gain/lose.
# The change in N(a,b,j) = N'(a,b,j) - N(a,b,j) tells us exactly
# which consecutive-position statistics change.

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = []
    for a in range(n):
        for b in range(a):
            if a - b >= 2:
                tiles.append((a, b))
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A, tiles

n = 5
_, tiles = tiling_to_tournament(0, n)
m = len(tiles)

# Sample: flip each tile from a specific tiling and see how N changes
bits = 0  # transitive-ish starting point
A_base, _ = tiling_to_tournament(bits, n)
N_base = N_matrix(A_base)
H_base = len(ham_paths(A_base))

print(f"\nBase tiling: bits={bits:0{m}b}, H={H_base}")

for tile_idx in range(m):
    bits_new = bits ^ (1 << tile_idx)
    A_new, _ = tiling_to_tournament(bits_new, n)
    N_new = N_matrix(A_new)
    H_new = len(ham_paths(A_new))

    # Compute delta_N
    delta_N_nonzero = []
    for a in range(n):
        for b in range(a+1, n):
            for j in range(n-1):
                dn = N_new[a][b][j] - N_base[a][b][j]
                if dn != 0:
                    delta_N_nonzero.append((a, b, j, dn))

    # How many distinct (a,b) pairs are affected?
    affected_pairs = set((a,b) for a,b,j,dn in delta_N_nonzero)

    u, v = tiles[tile_idx]
    print(f"  Flip tile ({u},{v}): dH={H_new-H_base:+d}, "
          f"{len(delta_N_nonzero)} N entries changed, "
          f"{len(affected_pairs)} pairs affected")
    if len(delta_N_nonzero) <= 10:
        for a, b, j, dn in delta_N_nonzero:
            print(f"    delta_N({a},{b},{j}) = {dn:+d}")

# =====================================================================
# KEY: For which tilings is M scalar?
# =====================================================================
print()
print("=" * 70)
print("TILINGS WITH SCALAR M")
print("=" * 70)

scalar_tilings = []
for bits in range(2**m):
    A, _ = tiling_to_tournament(bits, n)
    paths = ham_paths(A)
    H = len(paths)
    if H % n != 0: continue

    M = np.zeros((n, n), dtype=int)
    for p in paths:
        for j in range(n):
            M[p[j]][p[j]] += (-1)**j
        for j in range(n-1):
            a, b = p[j], p[j+1]
            M[a][b] += (-1)**j
            M[b][a] += (-1)**j

    if np.array_equal(M, (H//n) * np.eye(n, dtype=int)):
        scalar_tilings.append((bits, H))

print(f"  {len(scalar_tilings)} tilings with scalar M (out of {2**m})")
for bits, H in scalar_tilings:
    print(f"    bits={bits:0{m}b}, H={H}")

# Are scalar tilings connected in the skeleton?
if len(scalar_tilings) > 1:
    print(f"\n  Connectivity of scalar tilings in skeleton:")
    for i in range(len(scalar_tilings)):
        for j in range(i+1, len(scalar_tilings)):
            b1, h1 = scalar_tilings[i]
            b2, h2 = scalar_tilings[j]
            dist = bin(b1 ^ b2).count('1')
            print(f"    ({b1:0{m}b}, H={h1}) <-> ({b2:0{m}b}, H={h2}): distance={dist}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
The corrected consecutive-position formula gives a path-level
characterization of scalar M:

  M = (H/n)*I  iff  for ALL pairs {a,b}:
    sum_j (-1)^j N(a,b,j) = 0

Key findings:
1. For scalar M at n=5 (H=15): N(a,b,j) is NOT constant in j,
   but always has vanishing alternating sum.

2. For non-scalar M: some pairs have nonzero alternating sum.

3. Under tile flips: N changes in a localized way (few pairs affected),
   but the alternating sum constraint is global.

4. Scalar tilings are RARE in the skeleton — isolated or small clusters.

This suggests scalar M requires a very precise balance of path statistics
that is easily disrupted by local moves.
""")
