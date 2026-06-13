"""
Since Grok's projection-to-V argument fails (mixed span full 3D),
investigate WHY the bad-vertex TT is still rank-critical.

The issue: projecting to V forgets structure. The bad-TT column is
dependent on mixed TTs in V-projection, but NOT in the FULL edge space.
So the linear independence must come from the non-bad components.

Key question: What makes the bad-TT column unique in the full ∂₂ matrix?

opus-2026-03-09-S51
"""
import numpy as np
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


# Pick a specific n=5 tournament with #bad=3
n = 5
for bits in range(1 << (n*(n-1)//2)):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue
    bad = get_bad_vertices(A, n)
    if len(bad) != 3:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    edge_idx = {e: i for i, e in enumerate(edges)}

    # Find bad TT
    bad_tt_idx = None
    for idx_t, tt in enumerate(tts):
        if set(tt) == set(bad):
            bad_tt_idx = idx_t
            break

    mat = boundary2_matrix(A, n, edges, tts)
    full_rank = np.linalg.matrix_rank(mat, tol=1e-8)

    print(f"Tournament bits={bits}")
    print(f"  Bad vertices: {bad}")
    print(f"  Edges: {edges}")
    print(f"  TTs: {tts}")
    print(f"  Bad-TT index: {bad_tt_idx}, TT={tts[bad_tt_idx]}")
    print(f"  #TTs={len(tts)}, rank={full_rank}, redundancy={len(tts)-full_rank}")
    print()

    # Show the full ∂₂ matrix
    print("  Full ∂₂ matrix (rows=edges, cols=TTs):")
    for i, e in enumerate(edges):
        row = mat[i, :]
        nonzero = [(j, int(row[j])) for j in range(len(tts)) if abs(row[j]) > 1e-10]
        if nonzero:
            print(f"    {e}: " + ", ".join(f"TT{j}({tts[j]})={v}" for j, v in nonzero))

    # Classify TTs
    print(f"\n  TT classification:")
    for j, tt in enumerate(tts):
        bad_count = sum(1 for v in tt if v in bad)
        tp = 'BAD' if bad_count == 3 else ('MIXED' if bad_count == 2 else 'GOOD')
        col = mat[:, j]
        nonzero_edges = [(edges[i], int(col[i])) for i in range(len(edges)) if abs(col[i]) > 1e-10]
        rc_marker = " ← BAD-TT" if j == bad_tt_idx else ""
        print(f"    TT{j} {tt} [{tp}]: " +
              ", ".join(f"{e}={v}" for e, v in nonzero_edges) + rc_marker)

    # Now: can the bad-TT column be written as a linear combination of others?
    # If RC, it can't. Find the null space of mat to see the dependency.
    # The left null space tells us which row constraints separate bad-TT from others.

    # Find the cocycle that "sees" the bad-TT uniquely
    # A 1-cocycle w ∈ ker(∂₂^T) satisfies w · (each column of ∂₂) = 0
    # dim(ker(∂₂^T)) = n-1+β₁ = n-1 (since β₁=0)

    # Remove bad-TT column, find left null space
    mat_no_bad = np.delete(mat, bad_tt_idx, axis=1)
    # Left null space of mat_no_bad = ker(mat_no_bad^T)
    U, S, Vt = np.linalg.svd(mat_no_bad)
    tol = 1e-8
    rank_no_bad = np.sum(S > tol)
    # Left null space = columns of U corresponding to zero singular values
    left_null = U[:, rank_no_bad:]  # shape: (|edges|, dim_null)

    print(f"\n  rank(∂₂) = {full_rank}")
    print(f"  rank(∂₂ without bad-TT) = {rank_no_bad}")
    print(f"  dim(left null of ∂₂\\bad-TT) = {left_null.shape[1]}")

    # The extra null vector (compared to full ∂₂) is the one that
    # vanishes on all non-bad-TT columns but NOT on bad-TT column
    # Find it: it's in left_null but not in left_null of full mat
    U_full, S_full, _ = np.linalg.svd(mat)
    left_null_full = U_full[:, full_rank:]

    # The "witness" vector w: in left null of ∂₂\bad-TT but not of ∂₂
    # Project left_null onto complement of left_null_full
    # This gives us the direction that sees the bad-TT
    if left_null.shape[1] > left_null_full.shape[1]:
        # Find the extra direction
        # Project left_null_full out of left_null
        proj = left_null_full @ left_null_full.T  # projector onto left null of full
        complement = left_null - proj @ left_null
        # SVD to find nonzero directions
        _, sv, vt = np.linalg.svd(complement)
        witness_idx = np.argmax(sv)
        witness = complement[:, 0]  # take first column
        # Normalize
        witness = witness / np.linalg.norm(witness)

        # Check: w · (bad-TT column) should be nonzero
        bad_col = mat[:, bad_tt_idx]
        w_dot_bad = np.dot(witness, bad_col)

        print(f"\n  Witness cocycle w (sees bad-TT but not others):")
        for i, e in enumerate(edges):
            if abs(witness[i]) > 1e-10:
                print(f"    w{e} = {witness[i]:.4f}")
        print(f"  w · bad-TT-col = {w_dot_bad:.4f}")

        # Verify: w · (every other column) should be ~0
        max_other = max(abs(np.dot(witness, mat[:, j])) for j in range(len(tts)) if j != bad_tt_idx)
        print(f"  max |w · other-col| = {max_other:.2e}")

    print()
    break  # Just first example

# Now check: what's special about bad-TT column in the FULL space?
# At n=5 with redundancy=0, ALL TTs are RC. The projection argument
# fails because it only looks at 3 dimensions out of 10.
# The real reason bad-TT is RC is simply that ALL TTs are RC when
# redundancy=0 (square matrix, all columns needed).

print("\n" + "="*60)
print("KEY INSIGHT: WHY PROJECTION ARGUMENT FAILS")
print("="*60)
print("""
At n=5: redundancy=0 means ALL columns are needed (trivially RC).
Grok's projection is unnecessary — the full matrix is already tight.

At n=6: redundancy=2-3, so NOT all TTs are RC.
The question is WHY the bad-TT is STILL RC.

The V-projection argument fails because mixed TTs span all of V.
Mixed TTs touching each of the 3 bad edges provide e_xy, e_yz, e_xz
individually, spanning all 3D. No star-cocycle constraint prevents this.

The real mechanism must involve the FULL boundary structure, not just
the bad-edge projection.
""")

# Check: at n=6, which TTs are NOT rank-critical?
print("="*60)
print("n=6: WHICH TTs ARE rank-critical? (detailed example)")
print("="*60)

n = 6
count = 0
for bits in range(1 << (n*(n-1)//2)):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue
    bad = get_bad_vertices(A, n)
    if len(bad) != 3:
        continue

    count += 1
    if count > 1:  # Just first example
        break

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)
    full_rank = np.linalg.matrix_rank(mat, tol=1e-8)

    print(f"  Bad: {bad}, #TTs={len(tts)}, rank={full_rank}, redundancy={len(tts)-full_rank}")

    # Check each TT for rank-criticality
    rc_list = []
    nrc_list = []
    for j, tt in enumerate(tts):
        mat_j = np.delete(mat, j, axis=1)
        rj = np.linalg.matrix_rank(mat_j, tol=1e-8)
        bad_count = sum(1 for v in tt if v in bad)
        tp = 'BAD' if bad_count == 3 else ('MIXED' if bad_count == 2 else 'GOOD')
        is_rc = rj < full_rank
        if is_rc:
            rc_list.append((j, tt, tp))
        else:
            nrc_list.append((j, tt, tp))

    print(f"  Rank-critical ({len(rc_list)}):")
    for j, tt, tp in rc_list:
        print(f"    TT{j} {tt} [{tp}]")
    print(f"  NOT rank-critical ({len(nrc_list)}):")
    for j, tt, tp in nrc_list:
        print(f"    TT{j} {tt} [{tp}]")

    # Key: Are the non-RC TTs always GOOD? Always MIXED?
    nrc_types = [tp for _, _, tp in nrc_list]
    rc_types = [tp for _, _, tp in rc_list]
    print(f"\n  RC types: {dict(zip(*np.unique(rc_types, return_counts=True)))}")
    print(f"  Non-RC types: {dict(zip(*np.unique(nrc_types, return_counts=True)))}")


# Larger statistics: at n=6, what types are non-RC?
print("\n" + "="*60)
print("n=6: AGGREGATE — which TT types are RC vs non-RC?")
print("="*60)

from collections import Counter
rc_type_counts = Counter()
nrc_type_counts = Counter()
total = 0

for bits in range(1 << (n*(n-1)//2)):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue
    bad = get_bad_vertices(A, n)
    if len(bad) != 3:
        continue
    total += 1

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)
    full_rank = np.linalg.matrix_rank(mat, tol=1e-8)

    for j, tt in enumerate(tts):
        bad_count = sum(1 for v in tt if v in bad)
        tp = 'BAD' if bad_count == 3 else ('MIXED' if bad_count == 2 else 'GOOD')
        mat_j = np.delete(mat, j, axis=1)
        rj = np.linalg.matrix_rank(mat_j, tol=1e-8)
        if rj < full_rank:
            rc_type_counts[tp] += 1
        else:
            nrc_type_counts[tp] += 1

print(f"Total tournaments: {total}")
print(f"RC counts by type: {dict(rc_type_counts)}")
print(f"Non-RC counts by type: {dict(nrc_type_counts)}")
print(f"\nPer-tournament averages:")
for tp in ['BAD', 'MIXED', 'GOOD']:
    rc = rc_type_counts.get(tp, 0)
    nrc = nrc_type_counts.get(tp, 0)
    print(f"  {tp}: {rc/total:.2f} RC, {nrc/total:.2f} non-RC (total {(rc+nrc)/total:.2f})")
