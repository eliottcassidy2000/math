"""
Test Grok's CORRECTED quotient: R = B₁ / <Euler, HP 2-cycles> (NO stars).

Grok says stars "annihilate" (are relations), not generators. The quotient
by ONLY Euler + HP 2-cycles should be exactly 3D.

"HP 2-cycles" = for each Rédei HP, the alternating telescoping sum gives
relations among TT boundaries. We test multiple interpretations:
1. HP-consecutive TTs (individual TT boundaries along the HP)
2. HP alternating sum (Σ (-1)^i ∂(TT_{p_i,p_{i+1},p_{i+2}}))
3. ALL Rédei HPs' telescoping relations (not just one HP)
4. Euler = sum of all ∂₂ columns

opus-2026-03-09-S51j
"""
import numpy as np
from collections import Counter
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

def find_hamiltonian_path(A, n):
    """Greedy HP construction for tournaments."""
    path = [0]
    for v in range(1, n):
        inserted = False
        for i in range(len(path)):
            if A[v][path[i]] == 1:
                path.insert(i, v)
                inserted = True
                break
        if not inserted:
            path.append(v)
    return path

def find_all_hps_small(A, n, max_hps=50):
    """Find multiple HPs by permutation search (only feasible for small n)."""
    hps = []
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            hps.append(list(perm))
            if len(hps) >= max_hps:
                break
    return hps


# ============================================================
# TEST 1: R = B₁ / <Euler only>
# ============================================================
print("=" * 60)
print("TEST 1: R = B₁ / <Euler> (just sum of all ∂₂ columns)")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    dim_R_hist = Counter()
    tested = 0
    limit = None if n <= 5 else 500

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        tts = get_tts(A, n)
        mat = boundary2_matrix(A, n, edges, tts)
        rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)

        euler = np.sum(mat, axis=1).reshape(-1, 1)
        U_B, S_B, _ = np.linalg.svd(mat, full_matrices=True)
        B1_proj = U_B[:, :rank_B1]
        euler_in_B1 = B1_proj.T @ euler
        rank_euler = np.linalg.matrix_rank(euler_in_B1, tol=1e-8)
        dim_R = rank_B1 - rank_euler
        dim_R_hist[dim_R] += 1
        tested += 1

    print(f"  n={n}: tested={tested}, dim(R)={dict(sorted(dim_R_hist.items()))}")


# ============================================================
# TEST 2: R = B₁ / <Euler + single HP telescope>
# ============================================================
print("\n" + "=" * 60)
print("TEST 2: R = B₁ / <Euler + HP telescope (single HP)>")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    dim_R_hist = Counter()
    tested = 0
    limit = None if n <= 5 else 500

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        tts = get_tts(A, n)
        mat = boundary2_matrix(A, n, edges, tts)
        rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)

        hp = find_hamiltonian_path(A, n)
        generators = []

        # Euler
        generators.append(np.sum(mat, axis=1))

        # HP telescope: alternating sum of consecutive TT boundaries
        hp_tel = np.zeros(len(edges))
        for i in range(len(hp)-2):
            a, b, c = hp[i], hp[i+1], hp[i+2]
            if A[a][b] == 1 and A[b][c] == 1 and A[a][c] == 1:
                tt = (a, b, c)
                if tt in tts:
                    j = tts.index(tt)
                    hp_tel += (-1)**i * mat[:, j]
        generators.append(hp_tel)

        # Also add individual HP-consecutive TT boundaries
        for i in range(len(hp)-2):
            a, b, c = hp[i], hp[i+1], hp[i+2]
            if A[a][b] == 1 and A[b][c] == 1 and A[a][c] == 1:
                if (a,b,c) in tts:
                    j = tts.index((a,b,c))
                    generators.append(mat[:, j])

        gen_mat = np.array(generators).T
        U_B, S_B, _ = np.linalg.svd(mat, full_matrices=True)
        B1_proj = U_B[:, :rank_B1]
        gen_in_B1 = B1_proj.T @ gen_mat
        rank_gen = np.linalg.matrix_rank(gen_in_B1, tol=1e-8)
        dim_R = rank_B1 - rank_gen
        dim_R_hist[dim_R] += 1
        tested += 1

    print(f"  n={n}: tested={tested}, dim(R)={dict(sorted(dim_R_hist.items()))}")


# ============================================================
# TEST 3: R = B₁ / <Euler + ALL HP telescopes>
# ============================================================
print("\n" + "=" * 60)
print("TEST 3: R = B₁ / <Euler + ALL HP telescopes> (n=5 only)")
print("=" * 60)

n = 5
ne = n*(n-1)//2
dim_R_hist = Counter()
tested = 0

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)
    rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)

    generators = []
    generators.append(np.sum(mat, axis=1))  # Euler

    # ALL HPs
    hps = find_all_hps_small(A, n, max_hps=200)
    for hp in hps:
        # Alternating telescope
        hp_tel = np.zeros(len(edges))
        for i in range(len(hp)-2):
            a, b, c = hp[i], hp[i+1], hp[i+2]
            if A[a][b] == 1 and A[b][c] == 1 and A[a][c] == 1:
                if (a,b,c) in tts:
                    j = tts.index((a,b,c))
                    hp_tel += (-1)**i * mat[:, j]
        generators.append(hp_tel)

        # Individual consecutive TTs
        for i in range(len(hp)-2):
            a, b, c = hp[i], hp[i+1], hp[i+2]
            if A[a][b] == 1 and A[b][c] == 1 and A[a][c] == 1:
                if (a,b,c) in tts:
                    j = tts.index((a,b,c))
                    generators.append(mat[:, j])

    if not generators:
        tested += 1
        continue

    gen_mat = np.array(generators).T
    U_B, S_B, _ = np.linalg.svd(mat, full_matrices=True)
    B1_proj = U_B[:, :rank_B1]
    gen_in_B1 = B1_proj.T @ gen_mat
    rank_gen = np.linalg.matrix_rank(gen_in_B1, tol=1e-8)
    dim_R = rank_B1 - rank_gen
    dim_R_hist[dim_R] += 1
    tested += 1

print(f"  n={n}: tested={tested}, dim(R)={dict(sorted(dim_R_hist.items()))}")


# ============================================================
# TEST 4: What DOES give dim(R) = 3 exactly?
# ============================================================
print("\n" + "=" * 60)
print("TEST 4: WHAT subspace of B₁ gives codimension 3?")
print("=" * 60)

print("""
If dim(R) = 3 means dim(B₁) - dim(generators ∩ B₁) = 3, then at n=5
we need generators spanning dim(B₁)-3 = 6-3 = 3 dimensional subspace of B₁.
What natural 3D subspace of B₁ works?
""")

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

    # Bad vertices
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

    print(f"\nExample: bits={bits}, bad={bad}")
    print(f"  dim(B₁) = {rank_B1}, need codim-3 subspace = {rank_B1-3}D")

    # What if the "quotient" generators are the HIDDEN CYCLE LIFTS themselves?
    # Then R = B₁ / <z_v : v bad> has dim = rank_B1 - rank(lifts in B₁)
    # We know lifts span 3D => dim(R) = rank_B1 - 3 = 3 at n=5
    # This is trivially true but maybe what Grok means?
    print(f"  If R = B₁ / <hidden lifts>: dim = {rank_B1} - 3 = {rank_B1-3}")

    # The COMPLEMENT: what is the 3D subspace NOT spanned by non-bad TTs?
    bad_set = set(bad)
    bad_tt_indices = [j for j, tt in enumerate(tts) if set(tt) == bad_set]
    non_bad_indices = [j for j in range(len(tts)) if j not in bad_tt_indices]

    mat_non_bad = mat[:, non_bad_indices]
    rank_non_bad = np.linalg.matrix_rank(mat_non_bad, tol=1e-8)
    print(f"  rank(non-bad TTs) = {rank_non_bad}")
    print(f"  codim of non-bad in B₁ = {rank_B1 - rank_non_bad}")

    # What about TTs that DON'T involve any bad vertex?
    pure_good_indices = [j for j, tt in enumerate(tts) if not any(v in bad_set for v in tt)]
    if pure_good_indices:
        mat_pure = mat[:, pure_good_indices]
        rank_pure = np.linalg.matrix_rank(mat_pure, tol=1e-8)
        print(f"  rank(pure-good TTs, no bad vertex) = {rank_pure}")
        print(f"  codim = {rank_B1 - rank_pure}")

    break


# ============================================================
# TEST 5: ker(i_*) dimension — Grok's claim it's 3D
# ============================================================
print("\n" + "=" * 60)
print("TEST 5: dim(ker(i_*)) for inclusion T\\v → T")
print("Grok claims ker(i_*: H_1(T\\v) → H_1(T)) gives 3D global space")
print("=" * 60)

# This doesn't make sense directly — i_* goes from H_k(T\v) to H_k(T).
# For β₂: i_*: H_2(T\v) → H_2(T). But β₂=0 makes both sides trivial.
# For β₁: i_*: H_1(T\v) → H_1(T). When β₁(T)=0, target is 0.
# So ker(i_*) = H_1(T\v) which has dim = β₁(T\v) ∈ {0,1}.

# Maybe Grok means: over ALL vertices v, the total dimension?
# Σ_v dim(ker(i_*: H_1(T\v) → H_1(T))) = Σ_v β₁(T\v) = #bad ≤ 3

# Or maybe at the chain level: ker of the chain map C_k(T\v) → C_k(T)?

print("\nFor β₁(T)=0: ker(i_*: H_1(T\\v) → H_1(T)) = H_1(T\\v)")
print("So dim(ker) = β₁(T\\v) = 0 or 1")
print("Σ_v ker(i_*) = Σ_v β₁(T\\v) = #bad ≤ 3")
print("This is just HYP-282 restated. Not a new mechanism.")

# But let's check: is there a 3D COMBINED kernel space at the chain level?
for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    tested = 0
    limit = None if n <= 5 else 200

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break
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
            tested += 1
            continue

        if tested < 2:
            print(f"  bits={bits}, bad={bad}")
            print(f"  Σ_v β₁(T\\v) = {len(bad)} = #bad")
            print(f"  ker(i_*) per vertex: each 0 or 1, sum = {len(bad)}")
            print(f"  This IS the |bad|≤3 question repackaged")
        tested += 1

    print(f"  Total tested: {tested}")
