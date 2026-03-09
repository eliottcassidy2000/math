"""
Test Grok's "global LES" claim for HYP-282.

Grok claims: a global LES
  0 → ker(δ₃) → H₃(Δ) → ⊕_v H₃(Δ_v) → coker → 0
has rank(ker δ₃) = 3 exactly, independent of bad vertices.

This would be a LES involving the full path complex Δ(T) and all vertex
deletions Δ(T\v). Let's test:

1. Does such a sequence exist? What is δ₃?
2. Compute H_k(T, T\v) for all v via standard LES
3. Check if any natural global map has fixed-rank kernel = 3
4. Compute the "global restriction map" res: H_k(T) → ⊕_v H_k(T\v)
   and check dim(ker), dim(coker)

For β₂=0 context: the relevant level is k=1 (hidden cycles).
For β₃ context (kind-pasteur's work): k=3.

opus-2026-03-09-S51k
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

def boundary1_matrix(n, edges):
    """∂₁: C₁ → C₀. ∂₁(a,b) = b - a."""
    mat = np.zeros((n, len(edges)), dtype=float)
    for j, (a, b) in enumerate(edges):
        mat[b, j] += 1
        mat[a, j] -= 1
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
# TEST 1: Global restriction map res: H_1(T) → ⊕_v H_1(T\v)
# ============================================================
print("=" * 60)
print("TEST 1: Global restriction map at H₁ level")
print("When β₁(T)=0: H₁(T)=0, so res is trivially 0→⊕_v H₁(T\\v)")
print("This gives NO constraint on Σ β₁(T\\v).")
print("=" * 60)

# When β₁(T) = 0, H₁(T) = 0. The map res: 0 → ⊕_v H₁(T\v) has
# trivial kernel. There's no "ker(δ₃)" to speak of at this level.
print("\nβ₁(T)=0 ⟹ H₁(T)=0 ⟹ res: 0→⊕H₁(T\\v) has ker=0, coker=⊕H₁(T\\v)")
print("No 3D constraint possible from this direction.\n")


# ============================================================
# TEST 2: Relative homology H_2(T, T\v) at each vertex
# ============================================================
print("=" * 60)
print("TEST 2: H₂(T, T\\v) dimensions for each vertex v")
print("LES: H₂(T) → H₂(T,T\\v) → H₁(T\\v) → H₁(T)")
print("When β₂=β₁=0: H₂(T,T\\v) ≅ H₁(T\\v)")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    tested = 0
    limit = None if n <= 5 else 300
    rel_dim_sums = Counter()

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        bad_count = 0
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            if beta1(A_sub, n-1) == 1:
                bad_count += 1

        # H₂(T,T\v) ≅ H₁(T\v), so Σ dim H₂(T,T\v) = Σ β₁(T\v) = #bad
        rel_dim_sums[(bad_count, bad_count)] += 1  # trivially equal
        tested += 1

    print(f"  Σ_v dim H₂(T,T\\v) = Σ_v β₁(T\\v) = #bad (trivially)")
    print(f"  Distribution of #bad: {dict(Counter(k[0] for k,v in rel_dim_sums.items() for _ in range(v)))}")


# ============================================================
# TEST 3: The ACTUAL sequence Grok might mean
# At the CHAIN level: restriction C_k(T) → ⊕_v C_k(T\v)
# ============================================================
print("\n" + "=" * 60)
print("TEST 3: Chain-level restriction map")
print("res: Z₁(T) → ⊕_v Z₁(T\\v) / ⊕_v B₁(T\\v)")
print("= Z₁(T) → ⊕_v H₁(T\\v)")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    tested = 0
    limit = None if n <= 5 else 200
    ker_dims = Counter()
    coker_dims = Counter()

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        edge_idx_T = {e: i for i, e in enumerate(edges_T)}
        tts_T = get_tts(A, n)
        mat_T = boundary2_matrix(A, n, edges_T, tts_T)
        d1_T = boundary1_matrix(n, edges_T)

        # Z₁(T) = ker(∂₁) in C₁(T)
        # dim Z₁ = |E| - rank(∂₁) = C(n,2) - (n-1)
        rank_d1 = np.linalg.matrix_rank(d1_T, tol=1e-8)
        dim_Z1 = len(edges_T) - rank_d1

        # B₁(T) = im(∂₂) in Z₁(T)
        rank_d2 = np.linalg.matrix_rank(mat_T, tol=1e-8)

        # For each v, compute H₁(T\v) and the restriction map
        # The restriction takes a 1-cycle z ∈ Z₁(T) and restricts
        # it to edges not involving v, then projects to H₁(T\v)
        bad = []
        h1_gens = []
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            es_sub, h1 = compute_h1_generator(A_sub, n-1)
            if h1 is not None:
                bad.append(v)
                h1_gens.append((v, remaining, es_sub, h1))

        if len(bad) == 0:
            tested += 1
            continue

        # Build the global restriction map
        # res: Z₁(T) → ⊕_v H₁(T\v)
        # For each bad vertex v, H₁(T\v) is 1-dimensional
        # The map sends z ∈ Z₁(T) to (⟨z|_{T\v}, h1_gen_v⟩)_v
        # where the inner product detects the H₁ component

        # First, get Z₁(T) basis
        U_d1, S_d1, Vt_d1 = np.linalg.svd(d1_T, full_matrices=True)
        Z1_basis = Vt_d1[rank_d1:].T  # |E| × dim_Z1 (right null space)

        # For each bad v: the h1 generator as a vector in C₁(T) restricted to T\v edges
        # We need: for z ∈ Z₁(T), what's its H₁(T\v) component?
        # Restrict z to T\v edges, then project out B₁(T\v)

        # Build restriction matrices
        res_rows = []
        for v, remaining, es_sub, h1 in h1_gens:
            # h1 lives in C₁(T\v). Embed into C₁(T).
            edge_idx_sub = {}
            for k, (a,b) in enumerate(es_sub):
                ao, bo = remaining[a], remaining[b]
                if (ao, bo) in edge_idx_T:
                    edge_idx_sub[k] = edge_idx_T[(ao, bo)]

            h1_embed = np.zeros(len(edges_T))
            for k in range(len(h1)):
                if k in edge_idx_sub:
                    h1_embed[edge_idx_sub[k]] = h1[k]

            # The H₁ component of z restricted to T\v is ⟨z, h1_embed⟩ / ⟨h1, h1⟩
            # (since H₁ is 1D, this is a scalar)
            # In Z₁ coordinates:
            row = Z1_basis.T @ h1_embed  # dim_Z1 vector
            res_rows.append(row)

        if not res_rows:
            tested += 1
            continue

        # res_matrix: dim_Z1 → #bad
        res_mat = np.array(res_rows)  # (#bad × dim_Z1)
        rank_res = np.linalg.matrix_rank(res_mat, tol=1e-8)
        ker_dim = dim_Z1 - rank_res  # dim of {z ∈ Z₁ : invisible to all H₁(T\v)}
        coker_dim = len(bad) - rank_res

        ker_dims[(len(bad), ker_dim)] += 1
        coker_dims[(len(bad), coker_dim)] += 1

        if tested < 3 and len(bad) == 3:
            print(f"\n  bits={bits}, bad={bad}")
            print(f"  dim(Z₁) = {dim_Z1}, #bad = {len(bad)}")
            print(f"  rank(res) = {rank_res}")
            print(f"  ker(res) = {ker_dim}")
            print(f"  coker(res) = {coker_dim}")

        tested += 1

    print(f"\n  ker(res) by (#bad, dim): {dict(sorted(ker_dims.items()))}")
    print(f"  coker(res) by (#bad, dim): {dict(sorted(coker_dims.items()))}")


# ============================================================
# TEST 4: Same but restricted to B₁(T) instead of Z₁(T)
# res: B₁(T) → ⊕_v H₁(T\v)
# ============================================================
print("\n" + "=" * 60)
print("TEST 4: Restriction from B₁(T) → ⊕_v H₁(T\\v)")
print("When β₁=0: B₁=Z₁, so this should match TEST 3")
print("=" * 60)

# When β₁(T) = 0: Z₁ = B₁ (homology vanishes), so B₁ = Z₁.
# The maps are the same! This is consistent.
print("When β₁(T)=0: B₁(T) = Z₁(T), so res: B₁→⊕H₁ = res: Z₁→⊕H₁")
print("No new information.\n")


# ============================================================
# TEST 5: Is there a FIXED dimension among all tournaments?
# Check: dim(Z₁), dim(B₁), rank(res), ker(res), coker(res)
# ============================================================
print("=" * 60)
print("TEST 5: Which quantities are FIXED (independent of tournament)?")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    tested = 0
    limit = None if n <= 5 else 300

    quantities = {'dim_Z1': set(), 'dim_B1': set(), 'rank_res': set(),
                  'ker_res': set(), 'coker_res': set(), '#bad': set(),
                  'rank_res_from_B1': set()}

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        edge_idx_T = {e: i for i, e in enumerate(edges_T)}
        tts_T = get_tts(A, n)
        mat_T = boundary2_matrix(A, n, edges_T, tts_T)
        d1_T = boundary1_matrix(n, edges_T)

        rank_d1 = np.linalg.matrix_rank(d1_T, tol=1e-8)
        dim_Z1 = len(edges_T) - rank_d1
        rank_d2 = np.linalg.matrix_rank(mat_T, tol=1e-8)

        quantities['dim_Z1'].add(dim_Z1)
        quantities['dim_B1'].add(rank_d2)

        # Get bad vertices and their H₁ generators
        bad = []
        h1_embeds = []
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            es_sub, h1 = compute_h1_generator(A_sub, n-1)
            if h1 is not None:
                bad.append(v)
                h1_embed = np.zeros(len(edges_T))
                for k, (a, b) in enumerate(es_sub):
                    ao, bo = remaining[a], remaining[b]
                    if (ao, bo) in edge_idx_T:
                        h1_embed[edge_idx_T[(ao, bo)]] = h1[k]
                h1_embeds.append(h1_embed)

        quantities['#bad'].add(len(bad))

        if len(bad) == 0:
            tested += 1
            continue

        # Restriction map using Z₁ basis
        U_d1, S_d1, Vt_d1 = np.linalg.svd(d1_T, full_matrices=True)
        Z1_basis = Vt_d1[rank_d1:].T  # right null space

        res_mat = np.array([Z1_basis.T @ h for h in h1_embeds])
        rank_res = np.linalg.matrix_rank(res_mat, tol=1e-8)
        ker_res = dim_Z1 - rank_res
        coker_res = len(bad) - rank_res

        quantities['rank_res'].add(rank_res)
        quantities['ker_res'].add(ker_res)
        quantities['coker_res'].add(coker_res)

        tested += 1

    for name, vals in sorted(quantities.items()):
        fixed = "FIXED" if len(vals) == 1 else "VARIES"
        print(f"  {name}: {sorted(vals)} [{fixed}]")

    print(f"  (dim_Z₁ = C(n,2)-(n-1) = {n*(n-1)//2-(n-1)} always)")
    print(f"  (dim_B₁ = dim_Z₁ - β₁ = {n*(n-1)//2-(n-1)} when β₁=0)")


# ============================================================
# TEST 6: Mayer-Vietoris / Čech approach
# Consider the "global" boundary map
# ∂: C₂(T) → C₁(T)
# and restrict to the "local" part involving each vertex
# ============================================================
print("\n" + "=" * 60)
print("TEST 6: Vertex-star decomposition of C₂(T)")
print("Each TT (a,b,c) belongs to the star of vertex a.")
print("Group TT boundaries by their star-vertex.")
print("=" * 60)

n = 5
ne = n*(n-1)//2
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)
    rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)

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

    # Group TTs by star-vertex (first vertex)
    star_tts = {v: [] for v in range(n)}
    for j, (a,b,c) in enumerate(tts):
        star_tts[a].append(j)

    print(f"\nExample: bits={bits}, bad={bad}")
    for v in range(n):
        cols = star_tts[v]
        if cols:
            sub_mat = mat[:, cols]
            rank_v = np.linalg.matrix_rank(sub_mat, tol=1e-8)
        else:
            rank_v = 0
        bad_label = " [BAD]" if v in bad else ""
        print(f"  Star(v={v}): {len(cols)} TTs, rank={rank_v}{bad_label}")

    # What's the rank of non-star-bad TTs?
    non_bad_star = []
    for v in range(n):
        if v not in bad:
            non_bad_star.extend(star_tts[v])
    rank_non_bad_star = np.linalg.matrix_rank(mat[:, non_bad_star], tol=1e-8)

    bad_star = []
    for v in bad:
        bad_star.extend(star_tts[v])
    rank_bad_star = np.linalg.matrix_rank(mat[:, bad_star], tol=1e-8)

    print(f"  rank(non-bad stars) = {rank_non_bad_star}")
    print(f"  rank(bad stars) = {rank_bad_star}")
    print(f"  rank(all) = {rank_B1}")
    print(f"  codim(non-bad stars in B₁) = {rank_B1 - rank_non_bad_star}")

    break


# ============================================================
# TEST 7: The number 3 = C(n,2) - dim(B₁) - 1?
# At n=5: C(5,2) - 6 - 1 = 3. Coincidence?
# At n=6: C(6,2) - 10 - 1 = 4. NOT 3!
# ============================================================
print("\n" + "=" * 60)
print("TEST 7: Where does '3' come from?")
print("=" * 60)
for n in range(4, 9):
    dim_Z1 = n*(n-1)//2 - (n-1)
    print(f"  n={n}: dim(Z₁)=dim(B₁)={dim_Z1}, n-2={n-2}, C(n-1,2)={((n-1)*(n-2))//2}")

print("""
Observation:
  - dim(B₁) = C(n,2) - (n-1) = (n-1)(n-2)/2
  - This is C(n-1, 2).
  - #bad ≤ 3 holds for ALL n tested (5-10).
  - 3 is NOT C(n-1,2) - something for general n.
  - 3 = C(3,2) = the number of edges in the bad triple.
  - Or: 3 = the number of vertices in a triangle.
  - The constraint is NOT from a fixed-dim linear algebra argument
    unless the target space is truly 3D for all n.
""")
