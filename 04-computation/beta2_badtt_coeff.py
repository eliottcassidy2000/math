"""
Test Grok's B₁-proof claim: every z_v = Σ α_{v,TT} ∂(TT) has α≠0 on
bad-vertex TT, and the bad-TT boundaries span only 3D in B₁.

Grok says: {∂(bad-TT)} span ≤3D due to "global Euler + Rédei HP 2-cycles".

Key tests:
1. Is α_{v,bad-TT} ≠ 0 for ALL bad vertices v?
2. How many "bad TTs" are there? (3 bad vertices form exactly 1 TT)
3. Does ∂(bad-TT) span a small subspace of B₁?
4. What fraction of z_v's B₁-component comes from bad-TT vs good TTs?

NOTE: There is only ONE bad-TT (the TT formed by the 3 bad vertices).
So "bad-TT boundaries span ≤ 3D" is trivially true — it's 1D (one vector).
The question is whether z_v has nonzero coefficient on this one vector.

opus-2026-03-09-S51g
"""
import numpy as np
from itertools import permutations
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
    U, S, Vt = np.linalg.svd(mat)
    rank = np.sum(S > 1e-8)
    cocycles = U[:, rank:]
    if cocycles.shape[1] - (n - 1) <= 0: return edges, None
    coboundary_mat = np.zeros((len(edges), n), dtype=float)
    for v in range(n):
        for i, (a, b) in enumerate(edges):
            if b == v: coboundary_mat[i, v] += 1
            if a == v: coboundary_mat[i, v] -= 1
    cob_coords = np.linalg.lstsq(cocycles, coboundary_mat, rcond=None)[0]
    U_c, S_c, Vt_c = np.linalg.svd(cob_coords.T)
    rank_cob = np.sum(S_c > 1e-8)
    h1_coords = Vt_c[rank_cob:]
    if h1_coords.shape[0] == 0: return edges, None
    return edges, cocycles @ h1_coords[0]


# ============================================================
# TEST 1: Bad-TT coefficient in z_v decomposition
# ============================================================
print("=" * 60)
print("BAD-TT COEFFICIENT IN z_v = Σ α_j ∂₂(TT_j)")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")

    nonzero_count = 0
    zero_count = 0
    total_decomps = 0
    coeff_ratios = []  # |α_bad| / max|α_j|
    tested = 0
    limit = None if n <= 5 else 500

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break

        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        edge_idx_T = {e: i for i, e in enumerate(edges_T)}
        tts = get_tts(A, n)
        mat = boundary2_matrix(A, n, edges_T, tts)

        # Find bad vertices
        bad = []
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            if beta1(A_sub, n-1) == 1:
                bad.append(v)

        if len(bad) != 3:
            tested += 1
            continue

        # Find bad-TT index
        bad_set = set(bad)
        bad_tt_idx = None
        for idx_t, tt in enumerate(tts):
            if set(tt) == bad_set:
                bad_tt_idx = idx_t
                break
        if bad_tt_idx is None:
            tested += 1
            continue

        # For each bad vertex, decompose z_v as TT-boundary combination
        for v in bad:
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            edges_sub, h1 = compute_h1_generator(A_sub, n-1)
            if h1 is None:
                continue

            lift = np.zeros(len(edges_T))
            for k, (a_sub, b_sub) in enumerate(edges_sub):
                a_orig = remaining[a_sub]
                b_orig = remaining[b_sub]
                if (a_orig, b_orig) in edge_idx_T:
                    lift[edge_idx_T[(a_orig, b_orig)]] = h1[k]

            # Solve: lift = mat @ coeffs (least squares, should be exact)
            coeffs, _, _, _ = np.linalg.lstsq(mat, lift, rcond=None)
            residual = np.linalg.norm(mat @ coeffs - lift)
            assert residual < 1e-8, f"Decomposition failed: residual={residual}"

            alpha_bad = coeffs[bad_tt_idx]
            max_alpha = np.max(np.abs(coeffs))

            total_decomps += 1
            if abs(alpha_bad) > 1e-10:
                nonzero_count += 1
                coeff_ratios.append(abs(alpha_bad) / max_alpha)
            else:
                zero_count += 1

        tested += 1

    print(f"Total decompositions: {total_decomps}")
    print(f"α_bad ≠ 0: {nonzero_count} ({100*nonzero_count/max(1,total_decomps):.1f}%)")
    print(f"α_bad = 0: {zero_count} ({100*zero_count/max(1,total_decomps):.1f}%)")
    if coeff_ratios:
        print(f"|α_bad|/max|α|: mean={np.mean(coeff_ratios):.4f}, min={np.min(coeff_ratios):.4f}")


# ============================================================
# TEST 2: Decomposition uniqueness — is it unique?
# ============================================================
print("\n" + "=" * 60)
print("IS THE TT-BOUNDARY DECOMPOSITION UNIQUE?")
print("=" * 60)

print("""
When redundancy > 0, the decomposition z_v = Σ α_j ∂(TT_j) is NOT unique.
The null space of ∂₂ (right null space) gives freedom in choosing coefficients.
Grok's α_{v,bad-TT}≠0 claim might depend on the choice of decomposition.

Question: Is α_bad ALWAYS nonzero, or can we choose coefficients to make it 0?
This is equivalent to: can z_v be written ONLY using non-bad TT boundaries?
""")

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")

    can_avoid = 0
    must_use = 0
    tested = 0
    limit = None if n <= 5 else 500

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break

        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        edge_idx_T = {e: i for i, e in enumerate(edges_T)}
        tts = get_tts(A, n)
        mat = boundary2_matrix(A, n, edges_T, tts)

        bad = []
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            if beta1(A_sub, n-1) == 1:
                bad.append(v)

        if len(bad) != 3:
            tested += 1
            continue

        bad_set = set(bad)
        bad_tt_idx = None
        for idx_t, tt in enumerate(tts):
            if set(tt) == bad_set:
                bad_tt_idx = idx_t
                break
        if bad_tt_idx is None:
            tested += 1
            continue

        # Remove bad-TT column: can we still express each z_v?
        non_bad_cols = [j for j in range(len(tts)) if j != bad_tt_idx]
        mat_no_bad = mat[:, non_bad_cols]

        for v in bad:
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            edges_sub, h1 = compute_h1_generator(A_sub, n-1)
            if h1 is None:
                continue

            lift = np.zeros(len(edges_T))
            for k, (a_sub, b_sub) in enumerate(edges_sub):
                a_orig = remaining[a_sub]
                b_orig = remaining[b_sub]
                if (a_orig, b_orig) in edge_idx_T:
                    lift[edge_idx_T[(a_orig, b_orig)]] = h1[k]

            # Try to express lift using ONLY non-bad TT boundaries
            coeffs, _, _, _ = np.linalg.lstsq(mat_no_bad, lift, rcond=None)
            residual = np.linalg.norm(mat_no_bad @ coeffs - lift)

            if residual < 1e-8:
                can_avoid += 1
            else:
                must_use += 1

        tested += 1

    print(f"z_v CAN avoid bad-TT: {can_avoid}")
    print(f"z_v MUST use bad-TT: {must_use}")
    print(f"MUST use rate: {100*must_use/max(1,must_use+can_avoid):.1f}%")

    if must_use > 0:
        print(f"*** z_v REQUIRES bad-TT boundary — consistent with RC! ***")


# ============================================================
# TEST 3: Classify TTs by how many bad vertices they contain
# ============================================================
print("\n" + "=" * 60)
print("TT-BOUNDARY SUBSPACE STRUCTURE")
print("=" * 60)

n = 5
ne = n*(n-1)//2
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue
    bad = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            bad.append(v)
    if len(bad) != 3:
        continue

    bad_set = set(bad)
    edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges_T, tts)

    # Classify TTs
    tt_classes = {}
    for j, tt in enumerate(tts):
        bc = sum(1 for v in tt if v in bad_set)
        tp = f"{bc}-bad"
        if tp not in tt_classes:
            tt_classes[tp] = []
        tt_classes[tp].append(j)

    print(f"\nExample: bits={bits}, bad={bad}, #TTs={len(tts)}")
    for tp, indices in sorted(tt_classes.items()):
        cols = mat[:, indices]
        rank = np.linalg.matrix_rank(cols, tol=1e-8)
        print(f"  {tp} TTs ({len(indices)}): ∂₂-rank = {rank}")
        for j in indices:
            print(f"    TT{tts[j]}")

    # Full rank and its decomposition
    full_rank = np.linalg.matrix_rank(mat, tol=1e-8)
    print(f"  Full ∂₂ rank: {full_rank}")
    print(f"  dim(B₁) = {full_rank}")

    # Can non-bad TTs (0-bad + 1-bad + 2-bad) span all of B₁?
    non3_indices = [j for j, tt in enumerate(tts) if sum(1 for v in tt if v in bad_set) < 3]
    mat_non3 = mat[:, non3_indices]
    rank_non3 = np.linalg.matrix_rank(mat_non3, tol=1e-8)
    print(f"  Rank without 3-bad TT: {rank_non3}")
    print(f"  Rank drop when removing bad-TT: {full_rank - rank_non3}")

    break
