"""
Test: Is the bad-vertex transitive triple rank-critical at n=6?
At n=6, #bad=3 tournaments have redundancy 2-3 (not zero like n=5).
Question: Is the bad-vertex TT still rank-critical despite this?

Bad vertex = v where beta1(T\v) = 1
Bad-vertex TT = the unique TT formed by the 3 bad vertices (if transitive)

opus-2026-03-09-S51
"""
import numpy as np
from itertools import combinations, permutations

def tournament_from_bits(n, bits):
    A = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def get_tts(A, n):
    """Get all transitive triples"""
    tts = []
    for a in range(n):
        for b in range(n):
            if b == a: continue
            if A[a][b] == 0: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] == 1 and A[a][c] == 1:
                    tts.append((a,b,c))
    return tts

def boundary2_matrix(A, n, edges, tts):
    """Build ∂₂ restricted to Ω₂ (transitive triples)"""
    edge_idx = {e: i for i, e in enumerate(edges)}
    mat = np.zeros((len(edges), len(tts)), dtype=float)
    for j, (a,b,c) in enumerate(tts):
        # ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)
        if (b,c) in edge_idx: mat[edge_idx[(b,c)], j] += 1
        if (a,c) in edge_idx: mat[edge_idx[(a,c)], j] -= 1
        if (a,b) in edge_idx: mat[edge_idx[(a,b)], j] += 1
    return mat

def beta1(A, n):
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    if len(tts) == 0:
        return len(edges) - (n-1)  # all of Z1
    mat = boundary2_matrix(A, n, edges, tts)
    r = np.linalg.matrix_rank(mat, tol=1e-8)
    # beta1 = dim(Z1) - dim(B1) = (C(n,2) - (n-1)) - r
    return len(edges) - (n-1) - r

def get_bad_vertices(A, n):
    """Bad vertex = v where beta1(T\v) = 1"""
    bad = []
    for v in range(n):
        # Delete vertex v
        remaining = [i for i in range(n) if i != v]
        m = n - 1
        A_sub = np.zeros((m, m), dtype=int)
        for i2, i in enumerate(remaining):
            for j2, j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        if beta1(A_sub, m) == 1:
            bad.append(v)
    return bad

n = 6
ne = n*(n-1)//2
total = 1 << ne

count_tested = 0
bad3_count = 0
bad_tt_rc = 0
bad_tt_not_rc = 0
bad_not_tt = 0  # bad vertices don't form a TT

for bits in range(total):
    A = tournament_from_bits(n, bits)
    b1 = beta1(A, n)
    if b1 != 0:
        continue

    bad = get_bad_vertices(A, n)
    if len(bad) != 3:
        continue

    bad3_count += 1

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    edge_idx = {e: i for i, e in enumerate(edges)}

    # Check if bad vertices form a TT
    a, b, c = bad
    bad_tt_idx = None
    for perm in permutations(bad):
        x, y, z = perm
        if A[x][y] == 1 and A[y][z] == 1 and A[x][z] == 1:
            # (x,y,z) is a TT
            try:
                bad_tt_idx = tts.index((x,y,z))
            except ValueError:
                pass
            break

    if bad_tt_idx is None:
        bad_not_tt += 1
        continue

    # Build full ∂₂ matrix and check rank
    mat = boundary2_matrix(A, n, edges, tts)
    full_rank = np.linalg.matrix_rank(mat, tol=1e-8)

    # Remove bad-vertex TT column and check rank
    mat_reduced = np.delete(mat, bad_tt_idx, axis=1)
    reduced_rank = np.linalg.matrix_rank(mat_reduced, tol=1e-8)

    redundancy = len(tts) - full_rank

    if reduced_rank < full_rank:
        bad_tt_rc += 1
    else:
        bad_tt_not_rc += 1

    if bad3_count <= 5 or bad_tt_not_rc > 0:
        print(f"T{bad3_count}: bad={bad}, #TTs={len(tts)}, rank={full_rank}, "
              f"redundancy={redundancy}, bad-TT RC={reduced_rank < full_rank}")

print(f"\n=== SUMMARY n={n} ===")
print(f"Total with beta1=0 and #bad=3: {bad3_count}")
print(f"Bad vertices form TT: {bad_tt_rc + bad_tt_not_rc}")
print(f"Bad vertices NOT TT: {bad_not_tt}")
print(f"Bad-vertex TT is rank-critical: {bad_tt_rc}")
print(f"Bad-vertex TT is NOT rank-critical: {bad_tt_not_rc}")
