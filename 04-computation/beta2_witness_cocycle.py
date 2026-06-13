"""
Deep analysis of the witness cocycle ζ that certifies bad-TT rank-criticality.

Questions:
1. Is ζ unique (up to scaling)? — should be dim=1 extra null space
2. Is ζ the H₁ generator of the FLIPPED tournament?
3. Does ζ have a clean algebraic form?
4. How does ζ-multiplicity relate to |bad|?
5. Can we prove |bad|≤3 from ζ structure?

opus-2026-03-09-S51c
"""
import numpy as np
from itertools import permutations, combinations
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
    """Find the unique (up to scaling) 1-form ζ that:
    - vanishes on all non-bad-TT columns of ∂₂
    - is nonzero on the bad-TT column
    """
    full_rank = np.linalg.matrix_rank(mat, tol=1e-8)
    mat_no_bad = np.delete(mat, bad_tt_idx, axis=1)

    # Left null space of mat_no_bad
    U, S, Vt = np.linalg.svd(mat_no_bad)
    rank_no_bad = np.sum(S > 1e-8)
    left_null_no_bad = U[:, rank_no_bad:]

    # Left null space of full mat
    U_full, S_full, _ = np.linalg.svd(mat)
    left_null_full = U_full[:, full_rank:]

    # ζ is the direction in left_null_no_bad that's NOT in left_null_full
    # Project out the full null space
    if left_null_full.shape[1] > 0:
        proj = left_null_full @ left_null_full.T
        complement = left_null_no_bad - proj @ left_null_no_bad
        # SVD to find the nonzero direction
        _, sv, vt = np.linalg.svd(complement, full_matrices=False)
        # Take the direction with largest singular value
        best = np.argmax(sv)
        zeta = complement[:, best] if sv[best] > 1e-8 else None
    else:
        zeta = left_null_no_bad[:, 0]

    if zeta is not None:
        # Normalize so ⟨ζ, ∂_bad⟩ = 1
        bad_col = mat[:, bad_tt_idx]
        dot = np.dot(zeta, bad_col)
        if abs(dot) > 1e-10:
            zeta = zeta / dot

    return zeta


# ============================================================
# DETAILED ANALYSIS AT n=5
# ============================================================
print("=" * 60)
print("n=5: WITNESS COCYCLE ζ ANALYSIS")
print("=" * 60)

n = 5
ne = n*(n-1)//2
examples_shown = 0

for bits in range(1 << ne):
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

    # Find bad TT
    bad_tt_idx = None
    for idx_t, tt in enumerate(tts):
        if set(tt) == set(bad):
            bad_tt_idx = idx_t
            break

    zeta = find_witness_cocycle(mat, bad_tt_idx)

    if examples_shown < 3:
        print(f"\nTournament bits={bits}, bad={bad}")
        print(f"  Bad-TT: {tts[bad_tt_idx]}")
        print(f"  ζ (normalized so ⟨ζ,∂_bad⟩=1):")
        for i, e in enumerate(edges):
            if abs(zeta[i]) > 1e-10:
                print(f"    ζ({e}) = {zeta[i]:.6f}")

        # Check: is ζ a 1-cocycle of the ORIGINAL tournament?
        # A 1-cocycle satisfies ζ·(∂₂ column) = 0 for all TTs
        residuals = [np.dot(zeta, mat[:, j]) for j in range(len(tts))]
        print(f"  ⟨ζ, ∂₂ cols⟩: {[f'{r:.4f}' for r in residuals]}")
        print(f"  ζ is a 1-cocycle of T: {all(abs(r) < 1e-8 for r in residuals)}")
        print(f"  ζ is NOT a 1-cocycle (only of T\\bad-TT): correct")

        # Now compute the FLIPPED tournament and its H₁ generator
        # Flip: source→sink edge of bad-vertex TT
        bad_set = set(bad)
        # Find which bad vertex is source (beats both others) and sink
        for x in bad:
            others = [y for y in bad if y != x]
            if all(A[x][y] == 1 for y in others):
                source = x
            if all(A[y][x] == 1 for y in others):
                sink = x

        A_flip = A.copy()
        A_flip[source][sink] = 0
        A_flip[sink][source] = 1

        edges_flip = [(i,j) for i in range(n) for j in range(n) if i!=j and A_flip[i][j]==1]
        tts_flip = get_tts(A_flip, n)
        mat_flip = boundary2_matrix(A_flip, n, edges_flip, tts_flip)

        b1_flip = beta1(A_flip, n)
        print(f"\n  Flipped tournament: β₁={b1_flip}")

        if b1_flip == 1:
            # Find H₁ generator
            rank_flip = np.linalg.matrix_rank(mat_flip, tol=1e-8)
            U_f, S_f, _ = np.linalg.svd(mat_flip)
            left_null_flip = U_f[:, rank_flip:]

            # The 1-cocycles are left null space of ∂₂
            # H¹ = cocycles / coboundaries
            # Coboundary: δf(a,b) = f(b) - f(a)
            # dim(cocycles) = |E| - rank(∂₂) = C(n,2) - rank
            # dim(coboundaries) = n-1
            # β₁ = dim(cocycles) - (n-1)

            # Left null space of mat_flip gives ALL cocycles
            cocycles = left_null_flip  # columns are cocycle basis

            # Coboundary basis: δ(e_v) for each v
            edge_idx_flip = {e: i for i, e in enumerate(edges_flip)}
            coboundary_mat = np.zeros((len(edges_flip), n), dtype=float)
            for v in range(n):
                for i, (a, b) in enumerate(edges_flip):
                    if b == v:
                        coboundary_mat[i, v] += 1
                    if a == v:
                        coboundary_mat[i, v] -= 1

            # Project cocycles onto complement of coboundaries
            # H¹ rep = cocycle not in coboundary span
            cob_in_null = left_null_flip.T @ coboundary_mat  # project coboundaries into cocycle space
            # Find the direction in cocycle space orthogonal to coboundary projections
            if cocycles.shape[1] > 0:
                # cocycle space coords of coboundaries
                cob_coords = np.linalg.lstsq(cocycles, coboundary_mat, rcond=None)[0]
                # H¹ generator: cocycle direction not spanned by coboundary coords
                U_c, S_c, Vt_c = np.linalg.svd(cob_coords.T)
                rank_cob = np.sum(S_c > 1e-8)
                # kernel of cob_coords.T gives H¹ direction in cocycle-coordinate space
                h1_coords = Vt_c[rank_cob:]  # rows are kernel basis
                if h1_coords.shape[0] > 0:
                    h1_gen = cocycles @ h1_coords[0]
                    # Normalize
                    h1_gen = h1_gen / np.max(np.abs(h1_gen))

                    print(f"  H₁ generator of flipped T:")
                    for i, e in enumerate(edges_flip):
                        if abs(h1_gen[i]) > 1e-10:
                            print(f"    h₁({e}) = {h1_gen[i]:.6f}")

                    # Compare ζ (on original edges) with h₁ (on flipped edges)
                    # The edges differ only in (source,sink) ↔ (sink,source)
                    print(f"\n  Edge comparison (orig vs flipped):")
                    print(f"    Original has ({source},{sink}), flipped has ({sink},{source})")

                    # Map ζ to flipped edges: same for shared edges,
                    # ζ(source,sink) in original → -ζ for (sink,source) in flipped?
                    zeta_on_flip = np.zeros(len(edges_flip))
                    for i, e in enumerate(edges_flip):
                        if e in edge_idx:
                            zeta_on_flip[i] = zeta[edge_idx[e]]
                        elif e == (sink, source):
                            # This edge was (source, sink) in original
                            zeta_on_flip[i] = -zeta[edge_idx[(source, sink)]]

                    # Check if ζ (mapped) is proportional to h₁
                    # Normalize both
                    z_norm = zeta_on_flip / np.linalg.norm(zeta_on_flip)
                    h_norm = h1_gen / np.linalg.norm(h1_gen)
                    cos_angle = abs(np.dot(z_norm, h_norm))
                    print(f"    cos(angle between ζ_mapped and h₁) = {cos_angle:.6f}")
                    print(f"    ζ IS h₁ of flipped T: {cos_angle > 0.999}")

        examples_shown += 1

    if examples_shown >= 3:
        break


# ============================================================
# AGGREGATE: Is ζ always the H₁ generator of flipped T?
# ============================================================
print("\n" + "=" * 60)
print("n=5: IS ζ ALWAYS THE H₁ GENERATOR OF FLIPPED T?")
print("=" * 60)

n = 5
ne = n*(n-1)//2
match_count = 0
mismatch_count = 0

for bits in range(1 << ne):
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

    zeta = find_witness_cocycle(mat, bad_tt_idx)
    if zeta is None:
        continue

    # Find source/sink
    for x in bad:
        others = [y for y in bad if y != x]
        if all(A[x][y] == 1 for y in others):
            source = x
        if all(A[y][x] == 1 for y in others):
            sink = x

    # Flip
    A_flip = A.copy()
    A_flip[source][sink] = 0
    A_flip[sink][source] = 1

    edges_flip = [(i,j) for i in range(n) for j in range(n) if i!=j and A_flip[i][j]==1]
    tts_flip = get_tts(A_flip, n)
    edge_idx_flip = {e: i for i, e in enumerate(edges_flip)}
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
        match_count += 1
    else:
        mismatch_count += 1

print(f"ζ = H₁ generator of flipped T: {match_count}")
print(f"ζ ≠ H₁ generator: {mismatch_count}")
print(f"Match rate: {match_count/(match_count+mismatch_count)*100:.1f}%")


# ============================================================
# ζ MULTIPLICITY AND |bad| BOUND
# ============================================================
print("\n" + "=" * 60)
print("n=5,6: ζ-MULTIPLICITY AND |bad| RELATIONSHIP")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")

    # For each tournament with β₁=0, compute:
    # - |bad|
    # - dim of "extended left null space" (cocycles that kill all non-bad-vertex TTs)
    # This "ζ-multiplicity" should relate to #bad

    bad_count_hist = Counter()
    total = 0

    for bits in range(1 << ne):
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        bad = get_bad_vertices(A, n)
        nbad = len(bad)

        if nbad == 0:
            bad_count_hist[(nbad, 0)] += 1
            total += 1
            continue

        edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        tts = get_tts(A, n)
        mat = boundary2_matrix(A, n, edges, tts)
        full_rank = np.linalg.matrix_rank(mat, tol=1e-8)

        # Find ALL TTs involving only bad vertices
        bad_set = set(bad)
        bad_tt_indices = [j for j, tt in enumerate(tts) if all(v in bad_set for v in tt)]

        if not bad_tt_indices:
            bad_count_hist[(nbad, 0)] += 1
            total += 1
            continue

        # Remove ALL bad-vertex TTs and check rank drop
        non_bad_cols = [j for j in range(len(tts)) if j not in bad_tt_indices]
        mat_non_bad = mat[:, non_bad_cols] if non_bad_cols else np.zeros((len(edges), 0))
        rank_non_bad = np.linalg.matrix_rank(mat_non_bad, tol=1e-8) if mat_non_bad.shape[1] > 0 else 0

        rank_drop = full_rank - rank_non_bad

        # "ζ-multiplicity" = how many independent witness cocycles?
        # = rank_drop = number of dimensions the bad-vertex TTs contribute independently

        bad_count_hist[(nbad, rank_drop)] += 1
        total += 1

    print(f"Total β₁=0 tournaments: {total}")
    print(f"(#bad, rank_drop of bad-TTs) → count:")
    for key in sorted(bad_count_hist.keys()):
        nbad, rd = key
        print(f"  #bad={nbad}, rank_drop={rd}: {bad_count_hist[key]}")


# ============================================================
# n=5,6: CHECK Sum_v beta1(T\v) VALUES
# ============================================================
print("\n" + "=" * 60)
print("n=5,6: Sum_v β₁(T\\v) DISTRIBUTION")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    sum_hist = Counter()

    for bits in range(1 << ne):
        A = tournament_from_bits(n, bits)
        b1 = beta1(A, n)

        s = 0
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1, n-1), dtype=int)
            for i2, i in enumerate(remaining):
                for j2, j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            s += beta1(A_sub, n-1)

        sum_hist[(b1, s)] += 1

    print(f"(β₁(T), Σ_v β₁(T\\v)) → count:")
    for key in sorted(sum_hist.keys()):
        b1, s = key
        print(f"  β₁={b1}, Σ={s}: {sum_hist[key]}")
