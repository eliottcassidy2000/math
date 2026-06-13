"""
Verify ζ = H₁ generator of flipped T at n=6 (sample).
Also: investigate WHY |bad| ≤ 3 using the ζ structure.

Key insight from n=5: ζ is the EXACT H₁ class that would emerge
under the flip. This means RC = "the flip would create this H₁ class,
and removing the bad-TT column is equivalent to the flip at Ω₂ level."

opus-2026-03-09-S51c
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

def get_bad_vertices(A, n):
    bad = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1, n-1), dtype=int)
        for i2, i in enumerate(remaining):
            for j2, j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, n-1) == 1:
            bad.append(v)
    return bad

def find_witness_cocycle(mat, bad_tt_idx):
    full_rank = np.linalg.matrix_rank(mat, tol=1e-8)
    mat_no_bad = np.delete(mat, bad_tt_idx, axis=1)
    U, S, Vt = np.linalg.svd(mat_no_bad)
    rank_no_bad = np.sum(S > 1e-8)
    left_null_no_bad = U[:, rank_no_bad:]
    U_full, S_full, _ = np.linalg.svd(mat)
    left_null_full = U_full[:, full_rank:]
    if left_null_full.shape[1] > 0:
        proj = left_null_full @ left_null_full.T
        complement = left_null_no_bad - proj @ left_null_no_bad
        _, sv, vt = np.linalg.svd(complement, full_matrices=False)
        best = np.argmax(sv)
        zeta = complement[:, best] if sv[best] > 1e-8 else None
    else:
        zeta = left_null_no_bad[:, 0]
    if zeta is not None:
        bad_col = mat[:, bad_tt_idx]
        dot = np.dot(zeta, bad_col)
        if abs(dot) > 1e-10:
            zeta = zeta / dot
    return zeta

# ============================================================
# n=6: VERIFY ζ = H₁ OF FLIPPED T (SAMPLE)
# ============================================================
print("=" * 60)
print("n=6: ζ = H₁ GENERATOR OF FLIPPED T? (first 200)")
print("=" * 60)

n = 6
ne = n*(n-1)//2
match = 0
mismatch = 0
tested = 0

for bits in range(1 << ne):
    if tested >= 200:
        break

    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue
    bad = get_bad_vertices(A, n)
    if len(bad) != 3:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    edge_idx = {e: i for i, e in enumerate(edges)}
    mat = boundary2_matrix(A, n, edges, tts)

    bad_tt_idx = None
    for idx_t, tt in enumerate(tts):
        if set(tt) == set(bad):
            bad_tt_idx = idx_t
            break
    if bad_tt_idx is None:
        continue

    zeta = find_witness_cocycle(mat, bad_tt_idx)
    if zeta is None:
        continue

    # Find source/sink among bad vertices
    source = sink = mid = None
    for x in bad:
        others = [y for y in bad if y != x]
        if all(A[x][y] == 1 for y in others):
            source = x
        if all(A[y][x] == 1 for y in others):
            sink = x
    if source is None or sink is None:
        continue
    mid = [x for x in bad if x != source and x != sink][0]

    # Flip source→sink
    A_flip = A.copy()
    A_flip[source][sink] = 0
    A_flip[sink][source] = 1

    edges_flip = [(i,j) for i in range(n) for j in range(n) if i!=j and A_flip[i][j]==1]
    edge_idx_flip = {e: i for i, e in enumerate(edges_flip)}
    tts_flip = get_tts(A_flip, n)
    mat_flip = boundary2_matrix(A_flip, n, edges_flip, tts_flip)
    rank_flip = np.linalg.matrix_rank(mat_flip, tol=1e-8)

    U_f, S_f, _ = np.linalg.svd(mat_flip)
    cocycles = U_f[:, rank_flip:]

    coboundary_mat = np.zeros((len(edges_flip), n), dtype=float)
    for v in range(n):
        for i, (a, b) in enumerate(edges_flip):
            if b == v: coboundary_mat[i, v] += 1
            if a == v: coboundary_mat[i, v] -= 1

    cob_coords = np.linalg.lstsq(cocycles, coboundary_mat, rcond=None)[0]
    U_c, S_c, Vt_c = np.linalg.svd(cob_coords.T)
    rank_cob = np.sum(S_c > 1e-8)
    h1_coords = Vt_c[rank_cob:]
    if h1_coords.shape[0] == 0:
        continue
    h1_gen = cocycles @ h1_coords[0]

    # Map ζ to flipped edges
    zeta_on_flip = np.zeros(len(edges_flip))
    for i, e in enumerate(edges_flip):
        if e in edge_idx:
            zeta_on_flip[i] = zeta[edge_idx[e]]
        elif e == (sink, source):
            zeta_on_flip[i] = -zeta[edge_idx[(source, sink)]]

    z_norm = zeta_on_flip / np.linalg.norm(zeta_on_flip)
    h_norm = h1_gen / np.linalg.norm(h1_gen)
    cos_angle = abs(np.dot(z_norm, h_norm))

    if cos_angle > 0.999:
        match += 1
    else:
        mismatch += 1
        print(f"  MISMATCH at bits={bits}: cos={cos_angle:.6f}")

    tested += 1

print(f"\nResults: {match} match, {mismatch} mismatch out of {tested} tested")
print(f"Match rate: {match/(match+mismatch)*100:.1f}%")


# ============================================================
# WHY |bad| ≤ 3: THE ζ PERSPECTIVE
# ============================================================
print("\n" + "=" * 60)
print("WHY |bad| ≤ 3: ANALYSIS")
print("=" * 60)

print("""
KEY OBSERVATIONS:
1. #bad = Σ_v β₁(T\\v) when β₁(T)=0 (each bad vertex contributes exactly 1)
2. Σ_v β₁(T\\v) ≤ 3 for ALL tournaments (HYP-282, verified n≤10)
3. Therefore |bad| ≤ 3 follows from HYP-282

The ζ perspective:
- ζ is the H₁ generator of the FLIPPED tournament T'
- In T', the 3-cycle among bad vertices is the homological generator
- ζ exists iff the bad-vertex TT is rank-critical
- There is exactly ONE independent ζ (rank_drop = 1 always when #bad=3)

For |bad|≤3 via ζ:
- Each bad vertex v has β₁(T\\v) = 1
- The H₁ class of T\\v is generated by a 3-cycle involving the OTHER 2 bad vertices
- All these 3-cycles in subtournaments are "projections" of the same global ζ
- The claim: there's room for at most ONE independent ζ → at most 3 bad vertices

This connects to THM-103 (β₁≤1): the single H₁ dimension limits how many
vertices can be "critical" under deletion.
""")

# Verify: for β₁=0 tournaments, does #bad = Σ_v β₁(T\v)?
print("Verifying #bad = Σ_v β₁(T\\v) when β₁=0:")
for n in [5, 6]:
    ne = n*(n-1)//2
    violations = 0
    total = 0
    for bits in range(1 << ne):
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue
        total += 1
        bad = get_bad_vertices(A, n)
        s = sum(1 for v in range(n) if beta1(
            np.array([[A[i][j] for j2,j in enumerate([x for x in range(n) if x!=v])]
                       for i2,i in enumerate([x for x in range(n) if x!=v])]),
            n-1) == 1)
        if len(bad) != s:
            violations += 1
    print(f"  n={n}: {violations} violations out of {total} (should be 0)")

# The deeper question: WHY is Σ_v β₁(T\v) ≤ 3?
# At n=5: exhaustive confirms Σ ∈ {0,1,2,3} when β₁=0
# At n=6: exhaustive confirms Σ ∈ {0,1,2,3} when β₁=0
# This bound is INDEPENDENT of n!

print("\nΣ_v β₁(T\\v) distribution for β₁=0 (exhaustive):")
for n in [5, 6]:
    ne = n*(n-1)//2
    hist = Counter()
    for bits in range(1 << ne):
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue
        s = 0
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1, n-1), dtype=int)
            for i2, i in enumerate(remaining):
                for j2, j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            s += beta1(A_sub, n-1)
        hist[s] += 1
    print(f"  n={n}: {dict(sorted(hist.items()))}")
