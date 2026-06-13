"""
Test Grok's i_*-injectivity claim for β₂=0 proof.

Grok claims:
1. ζ_bad (witness cocycles) map to ker(i_*) which is "3D global span"
2. α=±1/√2 full-support forces their images to be basis of ker
3. |bad|≤3 follows from pigeon-hole

The inclusion i: T\v → T induces i_*: H_k(T\v) → H_k(T).
For β₁: i_*: H_1(T\v) → H_1(T).
When β₁(T)=0, target H_1(T)=0, so ker(i_*) = H_1(T\v).
Each ker is 0 or 1-dimensional. No "3D global space" possible.

BUT: maybe at chain level there's a common target space?
Test: embed all H_1(T\v) into C_1(T) and check span.

opus-2026-03-09-S51j
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
# TEST: Chain-level embedding of hidden cycles
# ============================================================
print("=" * 60)
print("CHAIN-LEVEL i_* ANALYSIS: Hidden cycles in C₁(T)")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n{'='*40}")
    print(f"n={n}")
    print(f"{'='*40}")

    tested = 0
    limit = None if n <= 5 else 300
    span_hist = Counter()
    B1_span_hist = Counter()

    for bits in range(1 << ne):
        if limit and tested >= limit:
            break
        A = tournament_from_bits(n, bits)
        if beta1(A, n) != 0:
            continue

        edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
        edge_idx = {e: i for i, e in enumerate(edges)}
        tts = get_tts(A, n)
        mat = boundary2_matrix(A, n, edges, tts)
        rank_B1 = np.linalg.matrix_rank(mat, tol=1e-8)

        # Find bad vertices and their hidden cycle lifts
        bad = []
        lifts = []
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            A_sub = np.zeros((n-1,n-1), dtype=int)
            for i2,i in enumerate(remaining):
                for j2,j in enumerate(remaining):
                    A_sub[i2][j2] = A[i][j]
            es, h1 = compute_h1_generator(A_sub, n-1)
            if h1 is not None:
                bad.append(v)
                # Embed into C₁(T)
                lift = np.zeros(len(edges))
                for k, (a, b) in enumerate(es):
                    ao, bo = remaining[a], remaining[b]
                    if (ao, bo) in edge_idx:
                        lift[edge_idx[(ao, bo)]] = h1[k]
                lifts.append(lift)

        if len(bad) == 0:
            tested += 1
            continue

        # Span in C₁(T)
        lift_mat = np.array(lifts)
        span_C1 = np.linalg.matrix_rank(lift_mat, tol=1e-8)
        span_hist[(len(bad), span_C1)] += 1

        # Project onto B₁
        U_B, S_B, _ = np.linalg.svd(mat, full_matrices=True)
        B1_proj = U_B[:, :rank_B1]
        lift_in_B1 = lift_mat @ B1_proj
        span_B1 = np.linalg.matrix_rank(lift_in_B1, tol=1e-8)
        B1_span_hist[(len(bad), span_B1)] += 1

        if tested < 3 and len(bad) == 3:
            print(f"\n  bits={bits}, bad={bad}")
            print(f"  Lift span in C₁: {span_C1}")
            print(f"  Lift span in B₁: {span_B1}")
            print(f"  dim(B₁) = {rank_B1}")

            # Check: are lifts linearly independent?
            # If they span exactly #bad, they're independent
            print(f"  Independent: {'YES' if span_C1 == len(bad) else 'NO'}")

            # Check support overlap
            for i, v in enumerate(bad):
                support = [e for j, e in enumerate(edges) if abs(lifts[i][j]) > 1e-10]
                involves_bad = [e for e in support if any(b in e for b in bad if b != v)]
                print(f"    z_{v}: {len(support)} edges, {len(involves_bad)} involve other bad vertices")

        tested += 1

    print(f"\n  C₁ span by (#bad, span): {dict(sorted(span_hist.items()))}")
    print(f"  B₁ span by (#bad, span): {dict(sorted(B1_span_hist.items()))}")


# ============================================================
# KEY TEST: Is there a COMMON 3D subspace containing all lifts
# across multiple tournaments?
# ============================================================
print("\n" + "=" * 60)
print("COMMON SUBSPACE TEST: Do lifts from different tournaments")
print("share a common ≤3D subspace of C₁(Kₙ)?")
print("=" * 60)

n = 5
ne = n*(n-1)//2

# Use universal edge ordering: (i,j) for i<j, with sign
# Edge (i,j) for i→j in T gets weight +w, (j,i) gets -w at same position
univ_edges = [(i,j) for i in range(n) for j in range(i+1, n)]
univ_idx = {}
for k, (i,j) in enumerate(univ_edges):
    univ_idx[(i,j)] = (k, +1)
    univ_idx[(j,i)] = (k, -1)
dim_univ = len(univ_edges)

all_univ_lifts = []
all_bad_counts = []

for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx = {e: i for i, e in enumerate(edges)}

    bad = []
    lifts = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        es, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is not None:
            bad.append(v)
            # Embed into UNIVERSAL edge space
            univ_lift = np.zeros(dim_univ)
            for k, (a, b) in enumerate(es):
                ao, bo = remaining[a], remaining[b]
                if (ao, bo) in univ_idx:
                    idx, sgn = univ_idx[(ao, bo)]
                    univ_lift[idx] += sgn * h1[k]
            all_univ_lifts.append(univ_lift)
    all_bad_counts.append(len(bad))

if all_univ_lifts:
    lift_mat = np.array(all_univ_lifts)
    total_span = np.linalg.matrix_rank(lift_mat, tol=1e-8)
    print(f"\nn=5: {len(all_univ_lifts)} lifts from {sum(1 for b in all_bad_counts if b>0)} tournaments")
    print(f"Total span in universal C(5,2)=10 edge space: {total_span}")
    print(f"If ≤3: there IS a universal 3D subspace")
    print(f"If >3: there is NOT")

    # Incremental: how fast does span grow?
    seen = np.zeros((0, dim_univ))
    ranks = []
    for lift in all_univ_lifts:
        seen = np.vstack([seen, lift.reshape(1,-1)])
        ranks.append(np.linalg.matrix_rank(seen, tol=1e-8))
    print(f"Rank growth: {ranks[:20]}...")
    print(f"Final rank: {ranks[-1]}")


# ============================================================
# GROK'S SPECIFIC CLAIM: α=±1/√2 forces basis
# ============================================================
print("\n" + "=" * 60)
print("α=±1/√2 ANALYSIS: Does this constrain anything?")
print("=" * 60)

print("""
Grok claims: Each ζ_bad has full support and α=±1/√2 coefficient on
the bad-TT boundary. Since ||∂(bad-TT)|| = √3, we have:
  ζ · ∂(bad-TT) = ±1/√2 * √3 = ±√(3/2)

But this is about the INNER PRODUCT of a 1-form with a boundary.
It constrains the projection, not the full vector.

The 3 lifts z_v (one per bad vertex) are each:
  z_v = (±1/√2) * ∂(bad-TT) + residual_v
  where residual_v ⊥ ∂(bad-TT)

The residuals span 2D (we showed rank 3 total = 1D bad-TT + 2D residual).
The "pigeonhole" would need: residuals live in a fixed 2D space.
But they DON'T live in a fixed space across tournaments!
""")

n = 5
ne = n*(n-1)//2

# Check residual structure for each 3-bad tournament
residual_dims = Counter()
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx = {e: i for i, e in enumerate(edges)}
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)

    bad = []
    lifts = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        es, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is not None:
            bad.append(v)
            lift = np.zeros(len(edges))
            for k, (a, b) in enumerate(es):
                ao, bo = remaining[a], remaining[b]
                if (ao, bo) in edge_idx:
                    lift[edge_idx[(ao, bo)]] = h1[k]
            lifts.append(lift)

    if len(bad) != 3:
        continue

    # Bad-TT direction
    bad_set = set(bad)
    bad_tt_idx = None
    for idx_t, tt in enumerate(tts):
        if set(tt) == bad_set:
            bad_tt_idx = idx_t
            break

    bad_dir = mat[:, bad_tt_idx]
    bad_dir_norm = bad_dir / np.linalg.norm(bad_dir)

    # Residuals
    residuals = []
    for lift in lifts:
        proj = np.dot(lift, bad_dir_norm) * bad_dir_norm
        residuals.append(lift - proj)

    res_mat = np.array(residuals)
    res_rank = np.linalg.matrix_rank(res_mat, tol=1e-8)
    residual_dims[res_rank] += 1

print(f"n=5 exhaustive: residual rank distribution = {dict(sorted(residual_dims.items()))}")
print(f"(Residuals = z_v minus bad-TT-direction projection)")
print()

# Check if residuals have any special relationship
# Like: sum of residuals = 0?
sum_zero = 0
sum_nonzero = 0
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx = {e: i for i, e in enumerate(edges)}
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)

    bad = []
    lifts = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        es, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is not None:
            bad.append(v)
            lift = np.zeros(len(edges))
            for k, (a, b) in enumerate(es):
                ao, bo = remaining[a], remaining[b]
                if (ao, bo) in edge_idx:
                    lift[edge_idx[(ao, bo)]] = h1[k]
            lifts.append(lift)

    if len(bad) != 3:
        continue

    # Check sum of lifts
    s = sum(lifts)
    if np.linalg.norm(s) < 1e-8:
        sum_zero += 1
    else:
        sum_nonzero += 1

print(f"Sum of 3 lifts = 0: {sum_zero} tournaments")
print(f"Sum of 3 lifts ≠ 0: {sum_nonzero} tournaments")
print(f"(If sum=0 always, then lifts satisfy 1 linear relation, span ≤ 2)")

# Check signed sum: z_a - z_b + z_c
signed_sum_in_B1 = 0
signed_sum_outside = 0
for bits in range(1 << ne):
    A = tournament_from_bits(n, bits)
    if beta1(A, n) != 0:
        continue

    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx = {e: i for i, e in enumerate(edges)}
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges, tts)

    bad = []
    lifts = []
    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        es, h1 = compute_h1_generator(A_sub, n-1)
        if h1 is not None:
            bad.append(v)
            lift = np.zeros(len(edges))
            for k, (a, b) in enumerate(es):
                ao, bo = remaining[a], remaining[b]
                if (ao, bo) in edge_idx:
                    lift[edge_idx[(ao, bo)]] = h1[k]
            lifts.append(lift)

    if len(bad) != 3:
        continue

    # Try all 8 sign patterns
    for s0 in [-1, 1]:
        for s1 in [-1, 1]:
            for s2 in [-1, 1]:
                combo = s0 * lifts[0] + s1 * lifts[1] + s2 * lifts[2]
                if np.linalg.norm(combo) < 1e-8:
                    pass  # counted above

print("\nCONCLUSION:")
print("Grok's i_*-injectivity argument reduces to |bad|=Σ β₁(T\\v) ≤ 3.")
print("But WHY Σ β₁(T\\v) ≤ 3 is exactly HYP-282 — the open question.")
print("α=±1/√2 is a real structural finding but doesn't constrain the COUNT.")
print("The 'pigeonhole on ker(i_*)' has no well-defined 3D target space.")
