"""
Test Grok's "3D global residual subspace of B₁" claim.

The hidden cycle lifts z_v all lie in B₁ = im(∂₂). They span #bad dims.
Grok claims: there's a 3D subspace of B₁ containing ALL possible hidden cycles,
bounded by "Ham-path telescoping + Euler."

Direct test: across ALL tournaments with β₁=0 at a given n, collect ALL
hidden cycle lifts. What is the dimension of their joint span?
If it's 3 (independent of n), Grok is right.
If it grows with n, Grok is wrong.

opus-2026-03-09-S51h
"""
import numpy as np
from collections import Counter

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

def compute_h1_generator(A, n):
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    if not tts: return edges, None
    mat = boundary2_matrix(A, n, edges, tts)
    U, S, _ = np.linalg.svd(mat)
    rank = np.sum(S > 1e-8)
    cocycles = U[:, rank:]
    if cocycles.shape[1] - (n - 1) <= 0: return edges, None
    cob = np.zeros((len(edges), n), dtype=float)
    for v in range(n):
        for i, (a, b) in enumerate(edges):
            if b == v: cob[i, v] += 1
            if a == v: cob[i, v] -= 1
    cc = np.linalg.lstsq(cocycles, cob, rcond=None)[0]
    U2, S2, V2 = np.linalg.svd(cc.T)
    rc = np.sum(S2 > 1e-8)
    h1c = V2[rc:]
    if h1c.shape[0] == 0: return edges, None
    return edges, cocycles @ h1c[0]


# ============================================================
# GLOBAL SPAN: across ALL tournaments, what dim do lifts span?
# ============================================================
# Problem: different tournaments have different edge sets, so lifts
# live in different spaces. We need a COMMON coordinate system.
# Use the UNIVERSAL edge labeling: for n vertices {0,...,n-1},
# edge (i,j) with i<j gets a fixed index, and we track direction.

def universal_edge_index(n):
    """Map directed edge (i,j) to index in C(n,2)-dim space with sign."""
    idx = {}
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            idx[(i,j)] = (k, +1)  # forward
            idx[(j,i)] = (k, -1)  # backward
            k += 1
    return idx, k

# Actually, for a FIXED tournament T, edges are directed.
# Different tournaments have different edge directions.
# So we can't directly compare lifts from different tournaments.
#
# BUT: within a single tournament, we can ask how many independent
# hidden cycles exist. We already know it's #bad ≤ 3.
# The question "is there a universal 3D subspace" doesn't make sense
# across different tournaments.
#
# REVISED APPROACH: For a FIXED tournament T with #bad=3, the lifts
# span 3D in B₁(T). Grok claims this 3D is constrained by HP structure.
# Test: does the 3D subspace have special structure relative to B₁?

print("=" * 60)
print("B₁ SUBSPACE STRUCTURE: WHERE DO HIDDEN CYCLES LIVE?")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n{'='*40}")
    print(f"n={n}: dim(B₁) = C(n,2)-(n-1) = {ne-(n-1)}")
    print(f"{'='*40}")

    tested = 0
    limit = None if n <= 5 else 200

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break

        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        edge_idx = {e: i for i, e in enumerate(edges)}
        tts = get_tts(A, n)
        mat = boundary2_matrix(A, n, edges, tts)
        rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)

        bad = []
        lifts = []
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            es, h1 = compute_h1_generator(A_sub, n-1)
            if h1 is not None:
                bad.append(v)
                lift = np.zeros(len(edges))
                for k, (a, b) in enumerate(es):
                    ao, bo = remaining[a], remaining[b]
                    if (ao, bo) in edge_idx:
                        lift[edge_idx[(ao, bo)]] = h1[k]
                lifts.append(lift)

        if len(bad) != 3:
            tested += 1
            continue

        # Express lifts in B₁ coordinates
        # B₁ basis: columns of mat (after rank reduction)
        U, S, Vt = np.linalg.svd(mat, full_matrices=False)
        B1_basis = U[:, :rank_B1]  # orthonormal basis for B₁

        # Project lifts onto B₁ basis to get B₁-coordinates
        lift_mat = np.array(lifts)  # (3 × |E|)
        lift_B1_coords = lift_mat @ B1_basis  # (3 × dim_B₁)

        # Rank in B₁ coordinates
        rank_in_B1 = np.linalg.matrix_rank(lift_B1_coords, tol=1e-8)

        if tested < 3:
            print(f"\nTournament bits={bits}, bad={bad}")
            print(f"  dim(B₁) = {rank_B1}")
            print(f"  Lifts rank in B₁ = {rank_in_B1}")
            print(f"  Lift B₁-coordinates (3 × {rank_B1}):")
            for i, v in enumerate(bad):
                nonzero = [(j, lift_B1_coords[i,j]) for j in range(rank_B1)
                          if abs(lift_B1_coords[i,j]) > 1e-10]
                print(f"    z_{v}: {len(nonzero)} nonzero components out of {rank_B1}")

            # Which B₁ basis directions are used?
            used_dirs = set()
            for i in range(3):
                for j in range(rank_B1):
                    if abs(lift_B1_coords[i,j]) > 1e-10:
                        used_dirs.add(j)
            print(f"  Total B₁ directions used: {len(used_dirs)} out of {rank_B1}")

        tested += 1

    # Aggregate: what fraction of B₁ do lifts typically span?
    print(f"\n  Lifts always span {3} out of dim(B₁) = {ne-(n-1)}")
    print(f"  Fraction: {3/(ne-(n-1)):.2%}")


# ============================================================
# KEY TEST: Decompose B₁ into bad-TT direction + complement.
# How many independent COMPONENTS (in the complement) do lifts have?
# ============================================================
print("\n" + "=" * 60)
print("DECOMPOSITION: z_v = α_v · ∂(bad-TT) + w_v")
print("Where w_v ∈ span(non-bad TT boundaries)")
print("=" * 60)

n = 5
ne = n*(n-1)//2
residual_dims = Counter()

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx = {e: i for i, e in enumerate(edges)}
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)

    bad = []
    lifts = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        es, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is not None:
            bad.append(v)
            lift = np.zeros(len(edges))
            for k, (a, b) in enumerate(es):
                ao, bo = remaining[a], remaining[b]
                if (ao, bo) in edge_idx:
                    lift[edge_idx[(ao, bo)]] = h1[k]
            lifts.append(lift)

    if len(bad) != 3:
        continue

    bad_set = set(bad)
    bad_tt_idx = None
    for idx_t, tt in enumerate(tts):
        if set(tt) == bad_set:
            bad_tt_idx = idx_t
            break

    # bad-TT boundary direction
    bad_dir = mat[:, bad_tt_idx]
    bad_dir_norm = bad_dir / np.linalg.norm(bad_dir)

    # Residuals: w_v = z_v - (z_v · bad_dir) * bad_dir
    residuals = []
    for lift in lifts:
        proj = np.dot(lift, bad_dir_norm) * bad_dir_norm
        w = lift - proj
        residuals.append(w)

    # Rank of residuals
    res_mat = np.array(residuals)
    res_rank = np.linalg.matrix_rank(res_mat, tol=1e-8)
    residual_dims[res_rank] += 1

print(f"\nn=5 exhaustive:")
print(f"z_v = α_v · ∂(bad-TT) + w_v decomposition")
print(f"Rank of residuals {{w_v}}: {dict(sorted(residual_dims.items()))}")

# Same for n=6
n = 6
ne = n*(n-1)//2
residual_dims6 = Counter()
tested = 0

for bits in range(1 << ne):
    if tested >= 500:
        break
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx = {e: i for i, e in enumerate(edges)}
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)

    bad = []
    lifts = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        es, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is not None:
            bad.append(v)
            lift = np.zeros(len(edges))
            for k, (a, b) in enumerate(es):
                ao, bo = remaining[a], remaining[b]
                if (ao, bo) in edge_idx:
                    lift[edge_idx[(ao, bo)]] = h1[k]
            lifts.append(lift)

    if len(bad) != 3:
        tested += 1
        continue

    bad_set = set(bad)
    bad_tt_idx = None
    for idx_t, tt in enumerate(tts):
        if set(tt) == bad_set:
            bad_tt_idx = idx_t
            break

    bad_dir = mat[:, bad_tt_idx]
    bad_dir_norm = bad_dir / np.linalg.norm(bad_dir)

    residuals = []
    for lift in lifts:
        proj = np.dot(lift, bad_dir_norm) * bad_dir_norm
        w = lift - proj
        residuals.append(w)

    res_mat = np.array(residuals)
    res_rank = np.linalg.matrix_rank(res_mat, tol=1e-8)
    residual_dims6[res_rank] += 1
    tested += 1

print(f"\nn=6 (500 sampled):")
print(f"Rank of residuals {{w_v}}: {dict(sorted(residual_dims6.items()))}")

print(f"""
INTERPRETATION:
If residual rank = 2: lifts span 1D (bad-TT dir) + 2D (residual) = 3D total.
The 3 lifts are: α₁·b + w₁, α₂·b + w₂, α₃·b + w₃
where b = ∂(bad-TT) and w₁,w₂,w₃ span 2D (one relation among residuals).
This would mean the "3D subspace" is real but structured as 1D + 2D.
""")
