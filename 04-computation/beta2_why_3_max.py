"""
WHY is the maximum hidden cycle dimension 3?

The lift dimension = #bad (each bad vertex contributes 1 independent direction).
So the question reduces to: why is #bad ≤ 3?

Approach: check if 4 bad vertices would force a LINEAR DEPENDENCY among lifts.
At n=7+, check tournaments with many bad vertices.

Key insight: each hidden cycle z_v is supported on edges NOT touching v.
For 4 bad vertices v₁,v₂,v₃,v₄, the supports overlap heavily.
Maybe the overlap forces a dependency?

Also: check Σ_v β₁(T\v) at n=7 exhaustively (impossible, 2^21).
Use sampling instead.

opus-2026-03-09-S51d
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
# STRUCTURAL ANALYSIS: Support overlaps of hidden cycles
# ============================================================
print("=" * 60)
print("SUPPORT STRUCTURE OF HIDDEN CYCLES")
print("=" * 60)

n = 5
ne = n*(n-1)//2
example_count = 0

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx_T = {e: i for i, e in enumerate(edges_T)}

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

    example_count += 1
    if example_count > 2:
        break

    print(f"\nTournament bits={bits}, bad={bad}")

    for v in bad:
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]

        edges_sub, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is None:
            continue

        # Lift to T's edge space
        lift = np.zeros(len(edges_T))
        for k, (a_sub, b_sub) in enumerate(edges_sub):
            a_orig = remaining[a_sub]
            b_orig = remaining[b_sub]
            if (a_orig, b_orig) in edge_idx_T:
                lift[edge_idx_T[(a_orig, b_orig)]] = h1[k]

        support = [edges_T[i] for i in range(len(edges_T)) if abs(lift[i]) > 1e-10]
        zero_edges = [edges_T[i] for i in range(len(edges_T)) if abs(lift[i]) < 1e-10]
        print(f"  Hidden cycle z_{v}:")
        print(f"    Support: {len(support)} edges (0 on {len(zero_edges)} edges touching v={v})")
        print(f"    Zero edges: {zero_edges}")
        for i, e in enumerate(edges_T):
            if abs(lift[i]) > 1e-10:
                print(f"      z_{v}({e}) = {lift[i]:.4f}")


# ============================================================
# KEY TEST: At n=7, sample and check Σ_v β₁(T\v)
# ============================================================
print("\n" + "=" * 60)
print("n=7,8: SAMPLED Σ_v β₁(T\\v) DISTRIBUTION")
print("=" * 60)

import random
random.seed(42)

for n in [7, 8]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    sum_hist = Counter()
    samples = 5000 if n <= 7 else 2000

    for _ in range(samples):
        bits = random.randint(0, (1 << ne) - 1)
        A = tournament_from_bits(n, bits)
        b1 = beta1(A, n)

        s = 0
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            s += beta1(A_sub, n-1)

        sum_hist[(b1, s)] += 1

    print(f"(β₁(T), Σ_v β₁(T\\v)) → count:")
    for key in sorted(sum_hist.keys()):
        b1, s = key
        print(f"  β₁={b1}, Σ={s}: {sum_hist[key]}")

    # Check maximum sum when β₁=0
    max_sum_b0 = max((s for (b1, s) in sum_hist.keys() if b1 == 0), default=0)
    print(f"Max Σ when β₁=0: {max_sum_b0}")
    if max_sum_b0 <= 3:
        print(f"*** HYP-282 HOLDS at n={n} (sampled) ***")
    else:
        print(f"*** HYP-282 VIOLATED at n={n}! ***")


# ============================================================
# ALGEBRAIC CONSTRAINT: Why can't 4 hidden cycles be independent?
# ============================================================
print("\n" + "=" * 60)
print("WHY ≤ 3: ALGEBRAIC STRUCTURE")
print("=" * 60)

print("""
Key observations:
1. Each hidden cycle z_v has ZERO on all edges touching v
2. z_v is supported on C(n-1,2) edges (the edges of T\\v)
3. For 4 bad vertices v₁,...,v₄, each z_vᵢ misses n-1 edges

But the total edge space has C(n,2) = n(n-1)/2 dimensions.
The 4 lifts together span ≤ C(n,2) dimensions.

CRITICAL: z_v is NOT arbitrary — it's an H₁ cocycle of T\\v.
It satisfies ALL TT constraints of T\\v.
When lifted, it satisfies TT constraints for ALL TTs not involving v.
The ONLY TT constraint it might violate is those involving v.

For β₁(T)=0, ALL cocycles of T are coboundaries.
So z_v (lifted) is NOT a cocycle of T.
The "defect" of z_v = its failure to satisfy TTs involving v.

For 4 bad vertices: each z_vᵢ fails at TTs involving vᵢ.
BUT: the TTs involving v₁,...,v₄ overlap heavily.
This might force a dependency among the defects.

Let's check: for each lift z_v, compute ⟨z_v, ∂₂(τ)⟩ for each TT τ.
""")

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

    edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx_T = {e: i for i, e in enumerate(edges_T)}
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges_T, tts)

    print(f"\nExample: bits={bits}, bad={bad}")
    print(f"TTs: {tts}")

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

        # Compute defect: ⟨z_v, ∂₂(τ)⟩ for each TT
        defects = mat.T @ lift  # (#TTs,)
        nonzero_defects = [(tts[j], defects[j]) for j in range(len(tts)) if abs(defects[j]) > 1e-10]
        print(f"  z_{v} defects (nonzero only):")
        for tt, d in nonzero_defects:
            involves_v = v in tt
            print(f"    TT{tt}: {d:.4f} {'(involves v)' if involves_v else '(NOT involving v!)'}")

    break  # Just one example
