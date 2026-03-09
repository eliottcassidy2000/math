"""
Test Grok's algebraic proof of HYP-282: Σ_v β₁(T\v) ≤ 3.

Grok's claim: local H₁ classes from each T\v embed into a quotient space
of dimension ≤ 3. Test this by:

1. For each bad vertex v (β₁(T\v)=1), compute the H₁ generator z_v of T\v
2. Lift z_v to a 1-form on T (extend by assigning values to edges touching v)
3. Check: do these lifts span a space of dimension ≤ 3?

Also test: what IS the obstruction space? When β₁(T)=0, Z¹(T) = coboundaries.
The restriction res_v: Z¹(T) → Z¹(T\v) is NOT surjective for bad v.
The "missing" cocycles of T\v form the cokernel of res_v.

opus-2026-03-09-S51d
"""
import numpy as np
from itertools import combinations
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

def compute_cocycles(A, n):
    """Return basis for Z¹(T) = left null space of ∂₂."""
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    if not tts:
        return edges, np.eye(len(edges))
    mat = boundary2_matrix(A, n, edges, tts)
    U, S, Vt = np.linalg.svd(mat)
    rank = np.sum(S > 1e-8)
    cocycles = U[:, rank:]  # left null space
    return edges, cocycles

def compute_h1_generator(A, n):
    """Return H¹ generator (cocycle not coboundary), or None if β₁=0."""
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    if not tts:
        return edges, None
    mat = boundary2_matrix(A, n, edges, tts)
    U, S, Vt = np.linalg.svd(mat)
    rank = np.sum(S > 1e-8)
    cocycles = U[:, rank:]
    b1 = cocycles.shape[1] - (n - 1)
    if b1 <= 0:
        return edges, None

    # Coboundary matrix
    coboundary_mat = np.zeros((len(edges), n), dtype=float)
    for v in range(n):
        for i, (a, b) in enumerate(edges):
            if b == v: coboundary_mat[i, v] += 1
            if a == v: coboundary_mat[i, v] -= 1

    # Project coboundaries into cocycle space
    cob_coords = np.linalg.lstsq(cocycles, coboundary_mat, rcond=None)[0]
    U_c, S_c, Vt_c = np.linalg.svd(cob_coords.T)
    rank_cob = np.sum(S_c > 1e-8)
    h1_coords = Vt_c[rank_cob:]
    if h1_coords.shape[0] == 0:
        return edges, None
    h1_gen = cocycles @ h1_coords[0]
    return edges, h1_gen

def beta1(A, n):
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    if not tts: return len(edges) - (n-1)
    mat = boundary2_matrix(A, n, edges, tts)
    r = np.linalg.matrix_rank(mat, tol=1e-8)
    return len(edges) - (n-1) - r

def subtournament(A, n, v):
    """Return A with vertex v removed."""
    remaining = [i for i in range(n) if i != v]
    m = n - 1
    A_sub = np.zeros((m, m), dtype=int)
    for i2, i in enumerate(remaining):
        for j2, j in enumerate(remaining):
            A_sub[i2][j2] = A[i][j]
    return A_sub, remaining


# ============================================================
# APPROACH 1: Lift H₁ generators to full edge space
# ============================================================
print("=" * 60)
print("APPROACH 1: LIFTING H₁ GENERATORS")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    lift_dims = Counter()
    tested = 0

    for bits in range(1 << ne):
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        # Find bad vertices
        bad = []
        h1_gens = {}  # v -> h1 generator of T\v (in T\v's edge space)
        for v in range(n):
            A_sub, remaining = subtournament(A, n, v)
            edges_sub, h1 = compute_h1_generator(A_sub, n-1)
            if h1 is not None:
                bad.append(v)
                h1_gens[v] = (edges_sub, h1, remaining)

        if not bad:
            continue

        # Lift each h1_gen to the full edge space of T
        edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        edge_idx_T = {e: i for i, e in enumerate(edges_T)}

        lifts = []
        for v in bad:
            edges_sub, h1, remaining = h1_gens[v]
            lift = np.zeros(len(edges_T))
            for k, (a_sub, b_sub) in enumerate(edges_sub):
                # Map back to original vertex labels
                a_orig = remaining[a_sub]
                b_orig = remaining[b_sub]
                if (a_orig, b_orig) in edge_idx_T:
                    lift[edge_idx_T[(a_orig, b_orig)]] = h1[k]
            # Edges touching v get value 0 (natural extension)
            lifts.append(lift)

        # Check span dimension of lifts
        if lifts:
            lift_mat = np.array(lifts)  # (#bad × |E|)
            dim = np.linalg.matrix_rank(lift_mat, tol=1e-8)
            lift_dims[dim] += 1

        tested += 1

    print(f"Tested: {tested}")
    print(f"Dimension of lift span: {dict(sorted(lift_dims.items()))}")


# ============================================================
# APPROACH 2: Restriction map cokernel
# ============================================================
print("\n" + "=" * 60)
print("APPROACH 2: RESTRICTION MAP ANALYSIS")
print("=" * 60)

print("""
When β₁(T)=0: Z¹(T) = coboundaries (dim n-1).
Restriction res_v: Z¹(T) → Z¹(T\\v) maps cocycles to cocycles.
For bad v: β₁(T\\v)=1, so Z¹(T\\v) has dim (n-2)+1 = n-1.
Coboundaries of T\\v have dim n-2.
res_v maps (n-1)-dim space into (n-1)-dim space.
H¹(T\\v) = 1-dim, so res_v CANNOT cover the H¹ direction.

Key: dim(coker(res_v)) = ?
""")

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    coker_dims = Counter()
    tested = 0

    for bits in range(1 << ne):
        if n == 6 and tested >= 500:
            break

        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        edge_idx_T = {e: i for i, e in enumerate(edges_T)}

        # Z¹(T) basis
        _, cocycles_T = compute_cocycles(A, n)  # columns = cocycle basis

        for v in range(n):
            A_sub, remaining = subtournament(A, n, v)
            b1_sub = beta1(A_sub, n-1)
            if b1_sub == 0:
                continue

            edges_sub = [(i,j) for i in range(n-1) for j in range(n-1)
                        if i!=j and A_sub[i][j]==1]
            _, cocycles_sub = compute_cocycles(A_sub, n-1)

            # Build restriction matrix: maps cocycles of T to edges of T\v
            # res_v(z)(a,b) = z(a,b) for a,b ≠ v
            # In coordinates: project cocycles_T onto edges not touching v
            edge_map = {}  # edge index in T\v → edge index in T
            for k, (a_sub, b_sub) in enumerate(edges_sub):
                a_orig = remaining[a_sub]
                b_orig = remaining[b_sub]
                if (a_orig, b_orig) in edge_idx_T:
                    edge_map[k] = edge_idx_T[(a_orig, b_orig)]

            # Restriction: select rows of cocycles_T corresponding to T\v edges
            res_mat = np.zeros((len(edges_sub), cocycles_T.shape[1]))
            for k, idx_t in edge_map.items():
                res_mat[k, :] = cocycles_T[idx_t, :]

            # Image of restriction in edge space of T\v
            im_res = res_mat  # columns = restricted cocycles
            rank_im = np.linalg.matrix_rank(im_res, tol=1e-8)

            # Z¹(T\v) dimension
            dim_Z1_sub = cocycles_sub.shape[1]

            # Check: does im(res_v) land inside Z¹(T\v)?
            # Project im_res columns onto cocycles_sub
            # residual should be zero if res maps cocycles to cocycles
            proj = cocycles_sub @ cocycles_sub.T  # projector onto Z¹(T\v)
            residual = im_res - proj @ im_res
            max_residual = np.max(np.abs(residual))

            coker_dim = dim_Z1_sub - rank_im
            coker_dims[(n, coker_dim)] += 1

        tested += 1

    print(f"Tested: {tested}")
    for key in sorted(coker_dims.keys()):
        nn, cd = key
        if nn == n:
            print(f"  coker(res_v) dim = {cd}: {coker_dims[key]} instances")


# ============================================================
# APPROACH 3: The "3D obstruction space" — does it exist?
# ============================================================
print("\n" + "=" * 60)
print("APPROACH 3: OBSTRUCTION SPACE DIMENSION")
print("=" * 60)

print("""
For β₁(T)=0, each bad vertex v contributes a "hidden cycle" z_v ∈ Z¹(T\\v)
that is NOT in im(res_v). The question: do these hidden cycles live in a
space of dimension ≤ 3?

More precisely: consider the UNION of all Z¹(T\\v) for bad v.
Map each to a common ambient space (e.g., R^{all edges of T}).
How many independent dimensions do the hidden cycles span?
""")

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    hidden_dims = Counter()
    tested = 0

    for bits in range(1 << ne):
        if n == 6 and tested >= 1000:
            break

        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        edge_idx_T = {e: i for i, e in enumerate(edges_T)}

        # Z¹(T) in full edge space
        _, cocycles_T = compute_cocycles(A, n)

        hidden_cycles = []  # lifted to T's edge space

        for v in range(n):
            A_sub, remaining = subtournament(A, n, v)
            b1_sub = beta1(A_sub, n-1)
            if b1_sub == 0:
                continue

            # Get H¹ generator of T\v
            edges_sub, h1_gen = compute_h1_generator(A_sub, n-1)
            if h1_gen is None:
                continue

            # Lift to T's edge space (zero on edges touching v)
            lift = np.zeros(len(edges_T))
            for k, (a_sub, b_sub) in enumerate(edges_sub):
                a_orig = remaining[a_sub]
                b_orig = remaining[b_sub]
                if (a_orig, b_orig) in edge_idx_T:
                    lift[edge_idx_T[(a_orig, b_orig)]] = h1_gen[k]

            hidden_cycles.append(lift)

        if hidden_cycles:
            hc_mat = np.array(hidden_cycles)
            dim = np.linalg.matrix_rank(hc_mat, tol=1e-8)
            hidden_dims[dim] += 1
        else:
            hidden_dims[0] += 1

        tested += 1

    print(f"Tested: {tested}")
    print(f"Dimension of hidden cycle span: {dict(sorted(hidden_dims.items()))}")
    max_dim = max(d for d in hidden_dims.keys())
    print(f"Maximum dimension: {max_dim}")
    if max_dim <= 3:
        print(f"*** CONFIRMS |bad| ≤ 3 mechanism: hidden cycles span ≤ 3D ***")
    else:
        print(f"*** EXCEEDS 3D — Grok's claim needs refinement ***")
