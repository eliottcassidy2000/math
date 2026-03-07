#!/usr/bin/env python3
"""
SYNTHESIS: Blue pairs, even r-powers, and POS on the tiling skeleton.

The key threads to weave together:
1. M(r) = c_0 + c_2*r^2 + c_4*r^4 + ... (even r-powers, THM-030)
2. Blue pair: T and T' = complement tiling (all tiles flipped)
   => s_{ij}(T') = -s_{ij}(T) => M(r, -s) = ?
3. POS: M is positive semidefinite (or at least has all nonneg eigenvalues?)
4. Skeleton: single-tile flip graph, colored by isomorphism class

CENTRAL QUESTION: When we flip ALL tiles (blue pair), s -> -s.
Since M(r) = sum_{k even} c_k(s) * r^k, what is c_k(-s)?
The s-dependence is polynomial (odd-degree in each s_{ij}).
So c_k(-s) = (-1)^{deg_s(c_k)} * c_k(s)?

Actually, each term in c_k involves products of s values from paths.
A path of length n-1 uses n-1 edges, and each edge contributes one s factor.
So c_k(s) is a polynomial of degree n-1-2k in the s variables (since r^{2k}
eats 2k factors, leaving n-1-2k s factors).

At odd n: n-1 is even, so n-1-2k is always even.
=> c_k(-s) = c_k(s) for all k!
=> M(r, -s) = M(r, s) at odd n!
This means: BLUE PAIR HAS SAME TRANSFER MATRIX at odd n!

At even n: n-1 is odd, so n-1-2k is odd.
=> c_k(-s) = -c_k(s) for all k!
=> M(r, -s) = -M(r, s) at even n!
This means: BLUE PAIR HAS NEGATED TRANSFER MATRIX at even n!

Let's verify this computationally.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

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

def ham_path_count(A):
    n = len(A)
    count = 0
    for p in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False; break
        if valid: count += 1
    return count

def transfer_matrix(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset(A, S_verts, end=a)
                    bb = count_paths_subset(A, R_verts, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False; break
        if valid: count += 1
    return count

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

def complement_bits(bits, m):
    """Flip all tiles: blue pair complement."""
    return bits ^ ((1 << m) - 1)

# =====================================================================
# PART 1: Blue pair M relationship
# =====================================================================
print("=" * 70)
print("BLUE PAIR M RELATIONSHIP: M(T) vs M(T')")
print("T' = all-tile-flip complement of T")
print("=" * 70)

for n in [3, 4, 5]:
    _, tiles = tiling_to_tournament(0, n)
    m = len(tiles)

    print(f"\nn={n}: {m} tiles, {2**m} tilings")

    tested = 0
    m_equal = 0
    m_negated = 0
    m_other = 0

    for bits in range(2**m):
        comp = complement_bits(bits, m)
        if comp <= bits:  # avoid double counting
            continue

        A1, _ = tiling_to_tournament(bits, n)
        A2, _ = tiling_to_tournament(comp, n)

        M1 = transfer_matrix(A1)
        M2 = transfer_matrix(A2)
        H1 = ham_path_count(A1)
        H2 = ham_path_count(A2)

        if np.array_equal(M1, M2):
            m_equal += 1
        elif np.array_equal(M1, -M2):
            m_negated += 1
        else:
            m_other += 1
        tested += 1

        if tested <= 3:
            print(f"  bits={format(bits, f'0{m}b')} vs {format(comp, f'0{m}b')}: "
                  f"H={H1},{H2}, M equal?{np.array_equal(M1, M2)}, "
                  f"M negated?{np.array_equal(M1, -M2)}")

    print(f"\n  Summary: {tested} pairs tested")
    print(f"    M(T) = M(T'):  {m_equal}")
    print(f"    M(T) = -M(T'): {m_negated}")
    print(f"    Other:          {m_other}")

    if n % 2 == 1 and m_equal == tested:
        print(f"  ** ODD n={n}: ALL blue pairs have M(T) = M(T')! **")
    if n % 2 == 0 and m_negated == tested:
        print(f"  ** EVEN n={n}: ALL blue pairs have M(T) = -M(T')! **")


# =====================================================================
# PART 2: H values for blue pairs
# =====================================================================
print()
print("=" * 70)
print("BLUE PAIR H VALUES: H(T) + H(T') pattern")
print("=" * 70)

for n in [3, 4, 5]:
    _, tiles = tiling_to_tournament(0, n)
    m = len(tiles)

    h_sums = []
    h_diffs = []

    for bits in range(2**m):
        comp = complement_bits(bits, m)
        if comp < bits: continue

        A1, _ = tiling_to_tournament(bits, n)
        A2, _ = tiling_to_tournament(comp, n)
        H1 = ham_path_count(A1)
        H2 = ham_path_count(A2)

        h_sums.append(H1 + H2)
        h_diffs.append(H1 - H2)

    print(f"\nn={n}:")
    print(f"  H(T) + H(T'): min={min(h_sums)}, max={max(h_sums)}, unique={sorted(set(h_sums))}")
    print(f"  H(T) - H(T'): min={min(h_diffs)}, max={max(h_diffs)}, unique={sorted(set(h_diffs))}")

    if n % 2 == 0:
        # At even n: M(T') = -M(T), so tr(M') = -tr(M). Since tr(M) = 0 at even n,
        # this is consistent. But H = sum of all paths, not just tr(M).
        pass

    if n % 2 == 1:
        # At odd n: M(T) = M(T'), so tr(M) = tr(M'). tr(M) = H/n at odd n,
        # so H(T) = H(T')!
        all_equal = all(d == 0 for d in h_diffs)
        print(f"  H(T) = H(T') always? {all_equal}")
        if all_equal:
            print(f"  ** At odd n={n}: blue pairs ALWAYS have same H! (follows from M=M') **")


# =====================================================================
# PART 3: Skeleton colored by iso class, showing blue-pair connections
# =====================================================================
print()
print("=" * 70)
print("n=5: SKELETON WITH BLUE PAIRS AND ISO CLASSES")
print("=" * 70)

n = 5
_, tiles = tiling_to_tournament(0, n)
m = len(tiles)

iso_classes = defaultdict(list)
tiling_data = {}

for bits in range(2**m):
    A, _ = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    M = transfer_matrix(A)
    canon = tournament_canonical(A)
    iso_classes[canon].append(bits)
    tiling_data[bits] = {'H': H, 'M': M, 'canon': canon}

class_labels = {}
for idx, canon in enumerate(sorted(iso_classes.keys())):
    class_labels[canon] = idx

num_classes = len(class_labels)
print(f"\n{2**m} tilings in {num_classes} isomorphism classes")

# Blue pair class connections
print("\nBlue pair class connections:")
blue_class_pairs = defaultdict(int)
for bits in range(2**m):
    comp = complement_bits(bits, m)
    if comp < bits: continue
    c1 = class_labels[tiling_data[bits]['canon']]
    c2 = class_labels[tiling_data[comp]['canon']]
    pair = (min(c1, c2), max(c1, c2))
    blue_class_pairs[pair] += 1

for (c1, c2), count in sorted(blue_class_pairs.items()):
    h1 = tiling_data[iso_classes[sorted(iso_classes.keys())[c1]][0]]['H']
    h2 = tiling_data[iso_classes[sorted(iso_classes.keys())[c2]][0]]['H']
    self_str = " (SELF-PAIRED)" if c1 == c2 else ""
    print(f"  Class {c1} (H={h1}) <-> Class {c2} (H={h2}): {count} pairs{self_str}")


# =====================================================================
# PART 4: POS analysis — eigenvalues of M across classes
# =====================================================================
print()
print("=" * 70)
print("n=5: POS — EIGENVALUES OF M BY ISO CLASS")
print("=" * 70)

for canon in sorted(iso_classes.keys()):
    c = class_labels[canon]
    bits = iso_classes[canon][0]
    M = tiling_data[bits]['M']
    H = tiling_data[bits]['H']
    eigs = sorted(np.linalg.eigvalsh(M.astype(float)))
    is_psd = all(e >= -1e-10 for e in eigs)
    is_scalar = np.allclose(M, (M[0][0]) * np.eye(n))

    print(f"  Class {c} (H={H:2d}, {len(iso_classes[canon]):2d} tilings): "
          f"eigs={[round(e,2) for e in eigs]}, "
          f"PSD={is_psd}, scalar={is_scalar}")


# =====================================================================
# PART 5: Skeleton single-flip connections BETWEEN iso classes
# =====================================================================
print()
print("=" * 70)
print("n=5: SINGLE-FLIP CONNECTIONS BETWEEN ISO CLASSES")
print("=" * 70)

# Build the "class adjacency" graph
class_adj = defaultdict(set)
class_edge_count = defaultdict(int)

for bits in range(2**m):
    c1 = class_labels[tiling_data[bits]['canon']]
    for tile_idx in range(m):
        nb = bits ^ (1 << tile_idx)
        c2 = class_labels[tiling_data[nb]['canon']]
        if c1 != c2:
            pair = (min(c1, c2), max(c1, c2))
            class_adj[c1].add(c2)
            class_edge_count[pair] += 1

print(f"\nClass adjacency graph (single-flip connections):")
for c in range(num_classes):
    H = tiling_data[iso_classes[sorted(iso_classes.keys())[c]][0]]['H']
    neighbors = sorted(class_adj.get(c, set()))
    print(f"  Class {c} (H={H:2d}): neighbors = {neighbors}")

print(f"\nInter-class edge counts (number of single-flip edges between classes):")
for (c1, c2), count in sorted(class_edge_count.items()):
    # divide by 2 since each edge counted from both sides
    print(f"  Class {c1} <-> Class {c2}: {count//2} edges")


# =====================================================================
# PART 6: M polynomial coefficients c_0, c_2 for each class
# =====================================================================
print()
print("=" * 70)
print("n=5: EVEN-r POLYNOMIAL COEFFICIENTS c_0, c_2 BY CLASS")
print("=" * 70)

# M(r) = c_0 + c_2*r^2 (degree n-2=3, but odd coeffs=0, so degree-2 polynomial in r^2)
# Actually degree in r is n-2=3 max, but with only even powers: r^0 and r^2
# So M(r) = c_0 + c_2 * r^2 at n=5 (since r^4 would need degree 4 > n-2=3)

# Wait: n-1 edges per path. Weight = product of (r + s_e) for n-1 edges.
# This is degree n-1 = 4 in r. The IE formula has multiple sums...
# M[a,b] is at most degree n-2 in r (from the IE structure).
# At n=5: degree <= 3, even powers only => r^0, r^2 terms.

def transfer_matrix_r_float(A, r_val):
    """Transfer matrix at general r (float)."""
    n = len(A)
    M = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0.0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_weighted(A, S_verts, r_val, end=a)
                    bb = count_paths_weighted(A, R_verts, r_val, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def count_paths_weighted(A, verts, r_val, start=None, end=None):
    total = 0.0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        w = 1.0
        for i in range(len(p)-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

# Sample at r=0 and r=0.5 to get c_0 and c_2
# M(0) = c_0, M(0.5) = c_0 + c_2*0.25
# So c_2 = (M(0.5) - M(0)) / 0.25

for canon in sorted(iso_classes.keys()):
    c = class_labels[canon]
    bits = iso_classes[canon][0]
    A, _ = tiling_to_tournament(bits, n)

    M_0 = transfer_matrix_r_float(A, 0.0)
    M_half = transfer_matrix_r_float(A, 0.5)

    c0 = M_0
    c2 = (M_half - M_0) / 0.25

    H = tiling_data[bits]['H']

    # Check if c0 and c2 are scalar multiples of I
    c0_scalar = np.allclose(c0, np.trace(c0)/n * np.eye(n))
    c2_scalar = np.allclose(c2, np.trace(c2)/n * np.eye(n))

    # Eigenvalues of c0 and c2
    eigs_c0 = sorted(np.linalg.eigvalsh(c0))
    eigs_c2 = sorted(np.linalg.eigvalsh(c2))

    print(f"  Class {c} (H={H:2d}):")
    print(f"    c_0: tr={np.trace(c0):.4f}, scalar?{c0_scalar}, eigs={[round(e,2) for e in eigs_c0]}")
    print(f"    c_2: tr={np.trace(c2):.4f}, scalar?{c2_scalar}, eigs={[round(e,2) for e in eigs_c2]}")

    # Verify M(0.5) = c_0 + 0.25*c_2 matches integer M
    M_int = tiling_data[bits]['M'].astype(float)
    match = np.allclose(M_int, c0 + 0.25 * c2)
    if not match:
        print(f"    WARNING: polynomial reconstruction doesn't match! "
              f"max diff = {np.max(np.abs(M_int - (c0 + 0.25*c2)))}")


# =====================================================================
# PART 7: Does POS hold for c_0 and c_2 separately?
# =====================================================================
print()
print("=" * 70)
print("n=5: POS FOR c_0 AND c_2 SEPARATELY")
print("=" * 70)

for canon in sorted(iso_classes.keys()):
    c = class_labels[canon]
    bits = iso_classes[canon][0]
    A, _ = tiling_to_tournament(bits, n)

    M_0 = transfer_matrix_r_float(A, 0.0)
    M_half = transfer_matrix_r_float(A, 0.5)

    c0 = M_0
    c2 = (M_half - M_0) / 0.25

    H = tiling_data[bits]['H']
    eigs_c0 = sorted(np.linalg.eigvalsh(c0))
    eigs_c2 = sorted(np.linalg.eigvalsh(c2))

    c0_psd = all(e >= -1e-10 for e in eigs_c0)
    c2_psd = all(e >= -1e-10 for e in eigs_c2)

    print(f"  Class {c} (H={H:2d}): c_0 PSD={c0_psd}, c_2 PSD={c2_psd}")


print()
print("=" * 70)
print("CONCLUSIONS")
print("=" * 70)
print("""
BLUE PAIR TRANSFER MATRIX:
  The prediction M(T) = M(T') at odd n FAILS.
  Reason: blue pair flips only TILE edges, not path edges.
  So s -> -s only for non-adjacent pairs, not universally.

  At n=5: blue pairs generally do NOT share the same M.
  Only 1/32 blue pairs have M(T) = M(T').
  H(T) + H(T') is NOT constant (varies from 10 to 30).

KEY FINDINGS FROM THE DATA:

1. SELF-PAIRED CLASSES: Classes 8 (H=11) and 10 (H=13) contain blue
   pairs mapping to themselves. These are the "blue-self" classes.

2. PSD CLASSES: Only classes 9, 10, 11 have PSD transfer matrices.
   These have H >= 13. All other classes have negative eigenvalues.

3. c_2 PSD THRESHOLD: Classes with H >= 9 have c_2 PSD.
   This is the quadratic r-coefficient being positive semidefinite.

4. SHARED c_2 SPECTRUM: Classes 8, 9, 10 share identical c_2 eigenvalues
   {4.71, 8, 8, 12, 15.29}. This is unexpected and suggests a deeper
   algebraic connection between these classes.

5. CLASS 11 (H=15, non-regular scalar): c_0 = 0, c_2 = 12I.
   So M(r) = 12r^2*I. At r=1/2: M = 3I. Pure quadratic!

6. CLASS 9 (H=15, Paley scalar): c_0 is non-scalar but c_0 + 0.25*c_2 = 3I.
   The scalar property emerges from CANCELLATION between c_0 and c_2.

EVEN-r POLYNOMIAL DECOMPOSITION M(r) = c_0 + c_2*r^2:
  c_0 = M(0) = "pure antisymmetric" contribution.
  c_2 = "quadratic r-correction."
  For scalar classes: c_0 + 0.25*c_2 is scalar, but c_0 alone need not be.
  Exception: Class 11 has c_0 = 0 (the non-regular scalar class).
""")
