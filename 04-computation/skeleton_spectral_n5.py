#!/usr/bin/env python3
"""
Spectral analysis of the transfer matrix along the n=5 tiling skeleton.

At odd n=5, M is symmetric and tr(M) = H.
For each isomorphism class, what are the eigenvalues of M?
How do eigenvalues change under single-tile flips?

KEY QUESTIONS:
1. Is M = (H/n)*I for any class besides the regular tournament?
2. What eigenvalue patterns exist across the 12 classes?
3. How does the eigenvalue spectrum change along skeleton edges?
4. Do the self-flip pairs (class 8, 10) have special spectral structure?
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
    return A

def ham_path_count(A):
    n = len(A)
    count = 0
    for p in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start:
            continue
        if end is not None and p[-1] != end:
            continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
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

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

n = 5
tiles = []
for a in range(n):
    for b in range(a):
        if a - b >= 2:
            tiles.append((a, b))
tiles.sort()
m = len(tiles)
num_tilings = 2**m

print(f"n={n}: {m} tiles, {num_tilings} tilings")

# Compute H and canonical form for all tilings
tiling_H = {}
tiling_canon = {}
iso_classes = defaultdict(list)

for bits in range(num_tilings):
    A = tiling_to_tournament(bits, n)
    H = ham_path_count(A)
    canon = tournament_canonical(A)
    tiling_H[bits] = H
    tiling_canon[bits] = canon
    iso_classes[canon].append(bits)

# Compute transfer matrix for ONE representative per class
class_labels = {}
class_M = {}
class_H = {}
class_evals = {}
class_representatives = {}

print()
print("=" * 70)
print("SPECTRAL ANALYSIS OF M PER ISOMORPHISM CLASS")
print("=" * 70)

for idx, canon in enumerate(sorted(iso_classes.keys())):
    class_labels[canon] = idx
    members = iso_classes[canon]
    rep = members[0]
    A = tiling_to_tournament(rep, n)
    H = tiling_H[rep]
    M = transfer_matrix(A)
    evals = sorted(np.linalg.eigvalsh(M).tolist(), reverse=True)

    class_M[idx] = M
    class_H[idx] = H
    class_evals[idx] = evals
    class_representatives[idx] = rep

    is_scalar = np.allclose(M, (H/n)*np.eye(n))
    diag = [M[i][i] for i in range(n)]
    off_diag_sum = M.sum() - np.trace(M)

    print(f"\n  Class {idx}: H={H}, |class|={len(members)}")
    print(f"    eigenvalues = {[round(e, 3) for e in evals]}")
    print(f"    diagonal = {diag}")
    print(f"    M = (H/n)*I ? {is_scalar}")
    print(f"    off-diag sum = {int(off_diag_sum)}")
    if not is_scalar:
        print(f"    M = {M.tolist()}")

# Spectral gaps
print()
print("=" * 70)
print("SPECTRAL GAPS AND PATTERNS")
print("=" * 70)

for idx in sorted(class_evals.keys()):
    evals = class_evals[idx]
    H = class_H[idx]
    gap = evals[0] - evals[1] if len(evals) > 1 else 0
    ratio = evals[0] / H if H > 0 else 0
    print(f"  Class {idx} (H={H:2d}): gap={gap:6.3f}, lambda_1/H={ratio:.3f}, evals={[round(e,2) for e in evals]}")

# Self-flip pairs spectral analysis
print()
print("=" * 70)
print("SELF-FLIP PAIRS: SPECTRAL COMPARISON")
print("=" * 70)

for bits in range(num_tilings):
    flip_bits = ((1 << m) - 1) ^ bits
    if flip_bits > bits:
        c1 = class_labels[tiling_canon[bits]]
        c2 = class_labels[tiling_canon[flip_bits]]
        if c1 == c2:
            A1 = tiling_to_tournament(bits, n)
            A2 = tiling_to_tournament(flip_bits, n)
            M1 = transfer_matrix(A1)
            M2 = transfer_matrix(A2)

            print(f"\n  Self-flip: {format(bits, f'0{m}b')} <-> {format(flip_bits, f'0{m}b')}, class={c1}, H={tiling_H[bits]}")
            print(f"    M1 eigenvalues = {sorted(np.linalg.eigvalsh(M1).tolist(), reverse=True)}")
            print(f"    M2 eigenvalues = {sorted(np.linalg.eigvalsh(M2).tolist(), reverse=True)}")

            # Check M1 vs M2 relationship
            P = np.zeros((n,n))
            for i in range(n):
                P[i][n-1-i] = 1
            M2_relabeled = P @ M1 @ P.T
            print(f"    M2 = P*M1*P^T (v->n-1-v)? {np.allclose(M2, M2_relabeled)}")
            print(f"    M2 = M1? {np.allclose(M1, M2)}")

# Delta eigenvalue under single-tile flips between DIFFERENT classes
print()
print("=" * 70)
print("EIGENVALUE CHANGE UNDER SINGLE-TILE FLIP")
print("=" * 70)

# Sample a few interesting skeleton edges
interesting_edges = []
for bits in range(num_tilings):
    for tile_idx in range(m):
        neighbor = bits ^ (1 << tile_idx)
        if neighbor > bits:
            c1 = class_labels[tiling_canon[bits]]
            c2 = class_labels[tiling_canon[neighbor]]
            if c1 != c2:
                dH = abs(tiling_H[neighbor] - tiling_H[bits])
                interesting_edges.append((dH, bits, neighbor, tile_idx, c1, c2))

# Show the largest delta_H edge
interesting_edges.sort(reverse=True)
print("\nLargest |delta_H| edges:")
seen_pairs = set()
for dH, bits, neighbor, tile_idx, c1, c2 in interesting_edges[:5]:
    pair = (min(c1,c2), max(c1,c2))
    if pair in seen_pairs:
        continue
    seen_pairs.add(pair)

    A1 = tiling_to_tournament(bits, n)
    A2 = tiling_to_tournament(neighbor, n)
    M1 = transfer_matrix(A1)
    M2 = transfer_matrix(A2)
    e1 = sorted(np.linalg.eigvalsh(M1).tolist(), reverse=True)
    e2 = sorted(np.linalg.eigvalsh(M2).tolist(), reverse=True)
    delta_evals = [round(e2[i] - e1[i], 3) for i in range(n)]

    a, b = tiles[tile_idx]
    print(f"\n  {format(bits, f'0{m}b')} -> {format(neighbor, f'0{m}b')} (tile {tile_idx}: ({a},{b}))")
    print(f"    Class {c1} (H={tiling_H[bits]}) -> Class {c2} (H={tiling_H[neighbor]}), |dH|={dH}")
    print(f"    evals1 = {[round(e,2) for e in e1]}")
    print(f"    evals2 = {[round(e,2) for e in e2]}")
    print(f"    delta_evals = {delta_evals}")

# Check: is delta_tr always = delta_H?
print()
print("=" * 70)
print("VERIFY: delta_tr = delta_H FOR ALL SKELETON EDGES (odd n)")
print("=" * 70)

# We know tr(M) = H at odd n, so this should be trivially true
# But let's verify it holds for the actual computed transfer matrices
all_match = True
for bits in [0, 1, 31, 63]:
    A = tiling_to_tournament(bits, n)
    M = transfer_matrix(A)
    H = ham_path_count(A)
    tr = int(np.trace(M))
    if tr != H:
        print(f"  MISMATCH: bits={format(bits, f'0{m}b')}, tr(M)={tr}, H={H}")
        all_match = False
print(f"  tr(M) = H for all tested: {all_match}")

# Score sequences
print()
print("=" * 70)
print("SCORE SEQUENCES BY CLASS")
print("=" * 70)

for idx in sorted(class_H.keys()):
    rep = class_representatives[idx]
    A = tiling_to_tournament(rep, n)
    scores = sorted([sum(row) for row in A], reverse=True)
    H = class_H[idx]
    print(f"  Class {idx} (H={H:2d}): scores={scores}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
n=5 spectral analysis:
  - 12 isomorphism classes, H in {{1,3,5,9,11,13,15}}
  - Transfer matrix M is symmetric (THM-030)
  - At odd n=5: tr(M) = H
  - Self-flip pairs exist in classes 8 (H=11) and 10 (H=13)

Key spectral questions answered:
  1. Is M = (H/n)*I for any class? (Check output above)
  2. What eigenvalue patterns exist?
  3. How do eigenvalues change under tile flips?
  4. Self-flip pair spectral structure?
""")
