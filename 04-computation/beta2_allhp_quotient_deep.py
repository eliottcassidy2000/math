"""
Deep analysis: R = B₁ / <Euler + ALL HP telescopes>.

At n=5: dim(R) = 0 for 600 tournaments, dim(R) = 2 for 120.
Question: Are the 120 exactly the #bad=3 tournaments?
And what IS the 2D residual?

opus-2026-03-09-S51j
"""
import numpy as np
from collections import Counter
from itertools import permutations

def tournament_from_bits(n, bits):
    A = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def get_tts(A, n):
    tts = []
    for a in range(n):
        for b in range(n):
            if b == a or A[a][b] == 0: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] == 1 and A[a][c] == 1:
                    tts.append((a,b,c))
    return tts

def boundary2_matrix(A, n, edges, tts):
    edge_idx = {e: i for i, e in enumerate(edges)}
    mat = np.zeros((len(edges), len(tts)), dtype=float)
    for j, (a,b,c) in enumerate(tts):
        if (b,c) in edge_idx: mat[edge_idx[(b,c)], j] += 1
        if (a,c) in edge_idx: mat[edge_idx[(a,c)], j] -= 1
        if (a,b) in edge_idx: mat[edge_idx[(a,b)], j] += 1
    return mat

def beta1(A, n):
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    if not tts: return len(edges) - (n-1)
    mat = boundary2_matrix(A, n, edges, tts)
    r = np.linalg.matrix_rank(mat, tol=1e-8)
    return len(edges) - (n-1) - r

def find_all_hps(A, n, max_hps=500):
    hps = []
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            hps.append(list(perm))
            if len(hps) >= max_hps:
                break
    return hps


# ============================================================
# CORRELATION: dim(R) vs #bad
# ============================================================
print("=" * 60)
print("CORRELATION: dim(B₁/<Euler + ALL HPs>) vs #bad")
print("=" * 60)

n = 5
ne = n*(n-1)//2
results = []

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)
    rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)

    # Count bad vertices
    bad_count = 0
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            bad_count += 1

    # Compute quotient with ALL HPs
    generators = [np.sum(mat, axis=1)]  # Euler
    hps = find_all_hps(A, n, max_hps=200)
    for hp in hps:
        hp_tel = np.zeros(len(edges))
        for i in range(len(hp)-2):
            a, b, c = hp[i], hp[i+1], hp[i+2]
            if A[a][b] == 1 and A[b][c] == 1 and A[a][c] == 1:
                if (a,b,c) in tts:
                    j = tts.index((a,b,c))
                    hp_tel += (-1)**i * mat[:, j]
        generators.append(hp_tel)

    gen_mat = np.array(generators).T
    U_B, S_B, _ = np.linalg.svd(mat, full_matrices=True)
    B1_proj = U_B[:, :rank_B1]
    gen_in_B1 = B1_proj.T @ gen_mat
    rank_gen = np.linalg.matrix_rank(gen_in_B1, tol=1e-8)
    dim_R = rank_B1 - rank_gen

    results.append((bits, bad_count, dim_R, len(hps)))

# Cross-tabulate
print(f"\nn=5 exhaustive ({len(results)} tournaments with β₁=0):")
cross = Counter()
for _, bad, dimR, _ in results:
    cross[(bad, dimR)] += 1
print(f"\n(#bad, dim_R) → count:")
for (bad, dimR), cnt in sorted(cross.items()):
    print(f"  #bad={bad}, dim(R)={dimR}: {cnt}")

# Are the dim(R)=2 cases exactly #bad=3?
dimR2_bad3 = sum(1 for _, bad, dimR, _ in results if dimR == 2 and bad == 3)
dimR2_total = sum(1 for _, _, dimR, _ in results if dimR == 2)
bad3_total = sum(1 for _, bad, _, _ in results if bad == 3)
print(f"\ndim(R)=2 AND #bad=3: {dimR2_bad3}")
print(f"dim(R)=2 total: {dimR2_total}")
print(f"#bad=3 total: {bad3_total}")
if dimR2_total == bad3_total and dimR2_bad3 == dimR2_total:
    print("*** PERFECT MATCH: dim(R)=2 ⟺ #bad=3 ***")
else:
    print(f"NOT perfect match")


# ============================================================
# WHAT IF: R = B₁ / <all HP-consecutive individual TTs>
# (not alternating sums, but the individual TT boundaries)
# ============================================================
print("\n" + "=" * 60)
print("ALT: B₁ / <all individual HP-consecutive TT boundaries>")
print("=" * 60)

dim_R2_hist = Counter()
cross2 = Counter()

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)
    rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)

    bad_count = 0
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            bad_count += 1

    # Collect ALL individual HP-consecutive TT boundaries
    hps = find_all_hps(A, n, max_hps=200)
    tt_indices_used = set()
    for hp in hps:
        for i in range(len(hp)-2):
            a, b, c = hp[i], hp[i+1], hp[i+2]
            if A[a][b] == 1 and A[b][c] == 1 and A[a][c] == 1:
                if (a,b,c) in tts:
                    tt_indices_used.add(tts.index((a,b,c)))

    # How many TTs are HP-consecutive vs total?
    generators = [mat[:, j] for j in tt_indices_used]
    if not generators:
        dim_R2_hist[rank_B1] += 1
        cross2[(bad_count, rank_B1)] += 1
        continue

    gen_mat = np.array(generators).T
    U_B, S_B, _ = np.linalg.svd(mat, full_matrices=True)
    B1_proj = U_B[:, :rank_B1]
    gen_in_B1 = B1_proj.T @ gen_mat
    rank_gen = np.linalg.matrix_rank(gen_in_B1, tol=1e-8)
    dim_R = rank_B1 - rank_gen

    dim_R2_hist[dim_R] += 1
    cross2[(bad_count, dim_R)] += 1

print(f"n=5: dim(R) = {dict(sorted(dim_R2_hist.items()))}")
print(f"Cross-tab:")
for (bad, dimR), cnt in sorted(cross2.items()):
    print(f"  #bad={bad}, dim(R)={dimR}: {cnt}")

# What fraction of TTs are HP-consecutive?
print("\n--- HP-consecutive TT coverage ---")
coverage = []
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue
    tts = get_tts(A, n)
    hps = find_all_hps(A, n, max_hps=200)
    tt_set = set()
    for hp in hps:
        for i in range(len(hp)-2):
            a, b, c = hp[i], hp[i+1], hp[i+2]
            if A[a][b] == 1 and A[b][c] == 1 and A[a][c] == 1:
                if (a,b,c) in tts:
                    tt_set.add(tts.index((a,b,c)))
    coverage.append(len(tt_set) / max(1, len(tts)))

print(f"HP-consecutive TT coverage: min={min(coverage):.2%}, max={max(coverage):.2%}, mean={np.mean(coverage):.2%}")
if min(coverage) == 1.0:
    print("*** ALL TTs are HP-consecutive for SOME HP! ***")
    print("So individual HP-consecutive TTs = ALL TTs = full B₁")
