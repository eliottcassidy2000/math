#!/usr/bin/env python3
"""
Transfer matrix on the GS subspace.

The transfer matrix M for a tournament T with backbone 0->1->...->n-1
has entries M[S,b][S',b'] where S is a set of visited vertices, b is the
current endpoint, and we extend by one step.

In the tiling model, the tournament is determined by the tiling bits.
The W-polynomial W(r) = sum_P prod(r + s_e) where s_e = A[p_i,p_{i+1}] - 1/2.

For GS tilings, the tournament has additional symmetry.
Question: does the transfer matrix restricted to GS tilings give the
skeleton adjacency matrix?

A simpler approach: the "H-vector" h_i = H(T_i) for each SC class.
Is there a matrix M such that M * h = lambda * h for some eigenvalue?
Or does the skeleton adjacency A relate to some transformation of H-values?

Let's explore: the GS flip graph maps GS tiling T to flip(T).
H(T) = #Hamiltonian paths of the tournament defined by tiling T.
For a GS tiling T in SC class i, flip(T) is in SC class j (j != i at odd n).
So H(T_i) tilings of class i get sent to class j, contributing to H(T_j).

But WAIT: flip is NOT class-preserving (different GS tilings of same class
can flip to different classes). So the "transfer" is many-to-many.

Let's compute: for each pair (i,j) of SC classes, count how many
GS tilings of class i flip to class j. Call this K[i][j].
This IS the weighted adjacency matrix of the skeleton!

Then: K[i][j] = #{GS tilings in class i that flip to class j}.
And degree(i) = sum_j K[i][j] = #{GS tilings in class i} (since all
flip to some class j at odd n).

The row sums of K equal the GS class sizes, which are always odd.

Question: does K have any nice algebraic structure?

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations, combinations
from collections import defaultdict

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def converse(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]

def is_SC(A, n):
    return canonical(A, n) == canonical(converse(A, n), n)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_ham_paths(A, n):
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def dp(mask, last):
        if mask == (1 << n) - 1:
            return 1
        total = 0
        for v in range(n):
            if not (mask & (1 << v)) and A[last][v]:
                total += dp(mask | (1 << v), v)
        return total
    return sum(dp(1 << v, v) for v in range(n))

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs = []
    fixed = []
    seen = set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen:
            continue
        ti, tj = n-1-j, n-1-i
        if ti > tj:
            ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx)
            seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx))
            seen.add(idx)
            seen.add(tidx)
    return pairs, fixed

def gen_gs_tilings(n, pairs, fixed):
    gs_dof = len(pairs) + len(fixed)
    result = []
    for free_val in range(2**gs_dof):
        bits = 0
        for k, (idx1, idx2) in enumerate(pairs):
            if (free_val >> k) & 1:
                bits |= (1 << idx1) | (1 << idx2)
        for k, fidx in enumerate(fixed):
            if (free_val >> (len(pairs) + k)) & 1:
                bits |= (1 << fidx)
        result.append(bits)
    return result

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

n = 5
m = num_tiling_bits(n)
print(f"TRANSFER MATRIX ANALYSIS at n={n}")
print(f"{'='*60}")

pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

# Build class database
canon_db = {}
class_list = []
bits_to_class = {}

for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    c = canonical(A, n)
    if c not in canon_db:
        canon_db[c] = len(class_list)
        class_list.append({
            'rep': A, 'sc': is_SC(A, n),
            'scores': score_seq(A, n),
            't3': count_3cycles(A, n),
            'H': count_ham_paths(A, n),
            'gs_tilings': set(), 'tilings': set()
        })
    idx = canon_db[c]
    class_list[idx]['tilings'].add(bits)
    bits_to_class[bits] = idx

for bits in gs_tilings:
    idx = bits_to_class[bits]
    class_list[idx]['gs_tilings'].add(bits)

sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
sc_to_local = {v: i for i, v in enumerate(sc_indices)}
nsc = len(sc_indices)

# Build K matrix
K = [[0]*nsc for _ in range(nsc)]
for bits in gs_tilings:
    c_from = bits_to_class[bits]
    flipped = flip_tiling(bits, m)
    c_to = bits_to_class[flipped]
    if c_from in sc_to_local and c_to in sc_to_local:
        K[sc_to_local[c_from]][sc_to_local[c_to]] += 1

# K is the skeleton adjacency matrix
print(f"\n  K matrix (rows/cols = SC classes in order {sc_indices}):")
header = "   " + " ".join(f"{sc_indices[j]:3d}" for j in range(nsc))
print(header)
for i in range(nsc):
    row = f"{sc_indices[i]:3d}" + " ".join(f"{K[i][j]:3d}" for j in range(nsc))
    sc = sc_indices[i]
    print(f"{row}  | GS={len(class_list[sc]['gs_tilings'])}, H={class_list[sc]['H']}, t3={class_list[sc]['t3']}")

# Check: K = K^T (symmetric)?
is_sym = all(K[i][j] == K[j][i] for i in range(nsc) for j in range(nsc))
print(f"\n  Symmetric: {is_sym}")

# Check: row sums = GS class sizes
for i, sc in enumerate(sc_indices):
    row_sum = sum(K[i])
    gs_size = len(class_list[sc]['gs_tilings'])
    assert row_sum == gs_size, f"Row sum {row_sum} != GS size {gs_size} for class {sc}"
print(f"  Row sums = GS class sizes: TRUE")

# K^2: what does K*K look like?
try:
    import numpy as np
    K_np = np.array(K, dtype=float)
    K2 = K_np @ K_np
    print(f"\n  K^2 matrix:")
    header2 = "   " + " ".join(f"{sc_indices[j]:4d}" for j in range(nsc))
    print(header2)
    for i in range(nsc):
        row2 = f"{sc_indices[i]:3d}" + " ".join(f"{int(K2[i][j]):4d}" for j in range(nsc))
        print(row2)

    # K^2 diagonal = sum of squares of K entries in each row
    # = number of "round trips" from class i back to class i via GS flip
    print(f"\n  K^2 diagonal (# GS round-trip paths):")
    for i, sc in enumerate(sc_indices):
        c = class_list[sc]
        print(f"    Class {sc} (t3={c['t3']}, H={c['H']}): K^2[i,i] = {int(K2[i][i])}")

    # Eigenvalues of K
    eigenvalues = np.linalg.eigvalsh(K_np)
    print(f"\n  K eigenvalues: {[f'{e:.4f}' for e in sorted(eigenvalues)[::-1]]}")

    # Check: sum of H*H values = sum of K^2?
    total_K2 = np.sum(K2)
    total_HH = sum(class_list[sc]['H']**2 for sc in sc_indices)
    print(f"\n  sum(K^2) = {total_K2}")
    print(f"  sum(H^2) for SC classes = {total_HH}")

    # Check: is K related to the GS projection of the Hamiltonian path matrix?
    # The total number of GS tilings is 2^gs_dof
    # For each GS tiling, it belongs to exactly one SC class
    # The "GS H-vector" g_i = #GS tilings in class i
    g = np.array([len(class_list[sc]['gs_tilings']) for sc in sc_indices])
    h = np.array([class_list[sc]['H'] for sc in sc_indices])
    print(f"\n  g (GS counts): {g}")
    print(f"  h (total H): {h}")
    print(f"  g/h (GS fraction): {np.round(g/h, 4)}")

    # K is symmetric with row sums g. Is K a doubly stochastic matrix after normalization?
    # K_normalized = diag(1/g) * K has row sums 1
    K_norm = K_np / g[:, np.newaxis]
    print(f"\n  K_norm row sums: {K_norm.sum(axis=1)}")
    print(f"  K_norm col sums: {K_norm.sum(axis=0)}")

    # Check: does K_norm have eigenvector proportional to g?
    Kg = K_np @ g
    print(f"\n  K * g = {Kg}")
    print(f"  Is K*g proportional to g? {np.allclose(Kg / g, Kg[0] / g[0])}")

except ImportError:
    print("\n  numpy not available")

print("\nDONE")
