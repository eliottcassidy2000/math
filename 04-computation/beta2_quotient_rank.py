"""
Test Grok's C₁/B₁ quotient argument for |bad|≤3.

Grok claims: 4 lifts in C₁/B₁ (dim n-1) satisfy alt-sum relation,
so rank ≤ 3. But we showed 3 lifts ALSO have alt-sum ∈ im(∂₂)!
So the same argument would give rank ≤ 2 for 3 lifts.

Key test: what is the ACTUAL rank of the lifts' images in C₁/B₁?
If it's < #bad, then C₁/B₁ loses information, and the argument
can't distinguish #bad from #bad-1.

opus-2026-03-09-S51f
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
    U, S, Vt = np.linalg.svd(mat)
    rank = np.sum(S > 1e-8)
    cocycles = U[:, rank:]
    b1 = cocycles.shape[1] - (n - 1)
    if b1 <= 0: return edges, None
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
    h1_gen = cocycles @ h1_coords[0]
    return edges, h1_gen


# ============================================================
# RANK OF LIFTS IN C₁ vs C₁/B₁
# ============================================================
print("=" * 60)
print("RANK OF HIDDEN CYCLE LIFTS: C₁ vs C₁/B₁")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    print(f"dim(C₁) = {n*(n-1)//2}")
    # When β₁=0: dim(B₁) = C(n,2) - (n-1), so dim(C₁/B₁) = n-1
    print(f"dim(C₁/B₁) = {n-1} when β₁=0")

    rank_pairs = Counter()  # (rank_C1, rank_quotient) -> count
    tested = 0
    limit = 1000 if n == 6 else None

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
        B1 = mat  # columns of ∂₂ span B₁

        # Get bad vertices and their lifts
        bad = []
        lifts = []
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            edges_sub, h1 = compute_h1_generator(A_sub, n-1)
            if h1 is not None:
                bad.append(v)
                lift = np.zeros(len(edges_T))
                for k, (a_sub, b_sub) in enumerate(edges_sub):
                    a_orig = remaining[a_sub]
                    b_orig = remaining[b_sub]
                    if (a_orig, b_orig) in edge_idx_T:
                        lift[edge_idx_T[(a_orig, b_orig)]] = h1[k]
                lifts.append(lift)

        if not lifts:
            tested += 1
            continue

        nbad = len(bad)

        # Rank in C₁
        lift_mat = np.array(lifts)  # (#bad × |E|)
        rank_C1 = np.linalg.matrix_rank(lift_mat, tol=1e-8)

        # Rank in C₁/B₁
        # Project lifts onto complement of B₁
        # B₁ = column space of mat (∂₂)
        rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)
        U_B, S_B, _ = np.linalg.svd(mat, full_matrices=True)
        # Projector onto B₁
        P_B = U_B[:, :rank_B1] @ U_B[:, :rank_B1].T
        # Project lifts onto complement of B₁
        lifts_mod_B1 = lift_mat - (lift_mat @ P_B.T)  # project out B₁ component
        # Wait, need to think about this correctly
        # Each lift z_v ∈ C₁. Its image [z_v] in C₁/B₁.
        # Rank of {[z_v]} = rank of {z_v mod B₁}
        # = rank of {z_v} projected onto B₁^⊥
        P_B_perp = np.eye(len(edges_T)) - P_B
        lifts_proj = (P_B_perp @ lift_mat.T).T  # (#bad × |E|)
        rank_quotient = np.linalg.matrix_rank(lifts_proj, tol=1e-8)

        rank_pairs[(nbad, rank_C1, rank_quotient)] += 1
        tested += 1

    print(f"Tested: {tested}")
    print(f"(#bad, rank_C₁, rank_{'{C₁/B₁}'}) → count:")
    for key in sorted(rank_pairs.keys()):
        nbad, rc, rq = key
        print(f"  #bad={nbad}, rank_C₁={rc}, rank_C₁/B₁={rq}: {rank_pairs[key]}")


# ============================================================
# CRITICAL: What ARE the lifts mod B₁?
# ============================================================
print("\n" + "=" * 60)
print("LIFTS MOD B₁: DETAILED STRUCTURE (n=5)")
print("=" * 60)

n = 5
ne = n*(n-1)//2
example_count = 0

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    bad = []
    lifts = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        edges_sub, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is not None:
            bad.append(v)
            lift = np.zeros(n*(n-1)//2)
            edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
            edge_idx_T = {e: i for i, e in enumerate(edges_T)}
            for k, (a_sub, b_sub) in enumerate(edges_sub):
                a_orig = remaining[a_sub]
                b_orig = remaining[b_sub]
                if (a_orig, b_orig) in edge_idx_T:
                    lift[edge_idx_T[(a_orig, b_orig)]] = h1[k]
            lifts.append(lift)

    if len(bad) != 3:
        continue

    example_count += 1
    if example_count > 1:
        break

    edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx_T = {e: i for i, e in enumerate(edges_T)}
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges_T, tts)
    rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)
    U_B, S_B, _ = np.linalg.svd(mat, full_matrices=True)
    P_B = U_B[:, :rank_B1] @ U_B[:, :rank_B1].T
    P_B_perp = np.eye(len(edges_T)) - P_B

    print(f"\nTournament bits={bits}, bad={bad}")
    print(f"dim(B₁) = {rank_B1}, dim(C₁/B₁) = {len(edges_T) - rank_B1}")

    lift_mat = np.array(lifts)
    lifts_proj = (P_B_perp @ lift_mat.T).T

    for i, v in enumerate(bad):
        print(f"\n  z_{v} mod B₁ (nonzero components):")
        for j, e in enumerate(edges_T):
            if abs(lifts_proj[i, j]) > 1e-10:
                print(f"    [{e}] = {lifts_proj[i, j]:.6f}")
        print(f"  norm in B₁⊥: {np.linalg.norm(lifts_proj[i]):.6f}")
        print(f"  norm in B₁: {np.linalg.norm(lift_mat[i] - lifts_proj[i]):.6f}")

    # Check: pairwise angles in C₁/B₁
    print(f"\n  Pairwise angles in C₁/B₁:")
    for i in range(3):
        for j in range(i+1, 3):
            a = lifts_proj[i]
            b = lifts_proj[j]
            cos = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-15)
            print(f"    z_{bad[i]} · z_{bad[j]} / norms = {cos:.6f}")

    # What is the rank?
    rank_q = np.linalg.matrix_rank(lifts_proj, tol=1e-8)
    print(f"\n  Rank of 3 lifts mod B₁ = {rank_q}")
    print(f"  (If rank=2, then ONE relation exists among the 3 lifts mod B₁)")

    # Check: is there an explicit relation?
    if rank_q == 0:
        print(f"  *** ALL lifts are in B₁! Quotient C₁/B₁ sees nothing. ***")
        # Verify: each lift is in im(∂₂)
        for i, v in enumerate(bad):
            coeff, _, _, _ = np.linalg.lstsq(mat, lift_mat[i], rcond=None)
            residual = np.linalg.norm(mat @ coeff - lift_mat[i])
            print(f"    z_{v} ∈ im(∂₂): residual = {residual:.2e}")


# ============================================================
# AGGREGATE: rank in C₁/B₁ for all #bad values
# ============================================================
print("\n" + "=" * 60)
print("AGGREGATE: rank mod B₁ vs #bad (n=5 exhaustive, n=6 sampled)")
print("=" * 60)

# Already printed above, but let's also check:
# If rank_quotient < #bad, then the quotient argument is WEAKER than C₁ rank.
# Grok's argument would need rank_quotient = #bad to work.

# The data already shows this. Let me also check: what is dim(C₁/B₁)?
for n in [5]:
    ne = n*(n-1)//2
    print(f"\n--- n={n}: dim(C₁/B₁) = {n-1} ---")
    print(f"If #bad lifts have rank r in C₁/B₁, then in C₁/B₁ we need r ≤ {n-1}")
    print(f"But if rank_quotient < #bad, the quotient loses info.")
    print(f"For Grok's argument: need rank_quotient = #bad, and #bad ≤ dim(C₁/B₁) = {n-1}")
    print(f"So the quotient bound gives #bad ≤ {n-1}, not #bad ≤ 3!")
    print(f"At n=5: #bad ≤ 4, at n=6: #bad ≤ 5, etc.")
    print(f"This is WEAKER than #bad ≤ 3!")
