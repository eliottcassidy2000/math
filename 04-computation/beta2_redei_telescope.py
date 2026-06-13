"""
Test Grok's Rédei telescoping argument for HYP-282.

Claim: Given a HP P of tournament T, lift each local ζ_v to a global cochain
via P-extension. Then for any 4 bad vertices v₁,...,v₄, the alternating sum
Σ (-1)^i ζ_{v_i} telescopes to an exact boundary on P-triples.

Test:
1. Find HP of T
2. For each bad vertex v, compute the "P-extended lift" of ζ_v
3. Check if alternating sums of 4 lifts vanish

But first: we've established #bad ≤ 3 empirically. So we can't test
with 4 bad vertices directly. Instead, test whether:
- The 3 hidden cycles z_v₁, z_v₂, z_v₃ have any special relationship via the HP
- Any alternating sum of 3 lifts has a clean algebraic structure
- The HP ordering constrains which vertices can be bad

opus-2026-03-09-S51e
"""
import numpy as np
from itertools import permutations
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

def find_hamiltonian_path(A, n):
    """Find a HP using greedy insertion (works for tournaments)."""
    path = [0]
    for v in range(1, n):
        # Insert v into path
        inserted = False
        for i in range(len(path)):
            if A[v][path[i]] == 1:  # v beats path[i]
                path.insert(i, v)
                inserted = True
                break
        if not inserted:
            path.append(v)
    # Verify
    for i in range(len(path)-1):
        assert A[path[i]][path[i+1]] == 1, f"Not a HP: {path[i]}→{path[i+1]}"
    return path

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
# WHERE ARE BAD VERTICES IN THE HAMILTONIAN PATH?
# ============================================================
print("=" * 60)
print("BAD VERTEX POSITIONS IN HAMILTONIAN PATH")
print("=" * 60)

for n in [5, 6]:
    ne = n*(n-1)//2
    print(f"\n--- n={n} ---")
    pos_hist = Counter()  # positions of bad vertices in HP
    gap_hist = Counter()  # gaps between consecutive bad vertices
    consec_hist = Counter()  # are bad vertices consecutive in HP?
    total = 0

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

        if not bad:
            continue

        hp = find_hamiltonian_path(A, n)
        hp_pos = {v: i for i, v in enumerate(hp)}

        bad_positions = sorted([hp_pos[v] for v in bad])
        for p in bad_positions:
            pos_hist[(n, p)] += 1

        if len(bad) == 3:
            # Check if bad vertices are consecutive in HP
            bp = bad_positions
            consec = (bp[2] - bp[0] == 2)
            consec_hist[(n, consec)] += 1

            # Gaps
            gaps = tuple(bp[i+1] - bp[i] for i in range(2))
            gap_hist[(n, gaps)] += 1

        total += 1

    print(f"Total with bad verts: {total}")
    print(f"Bad vertex HP positions (position → count):")
    for key in sorted(pos_hist.keys()):
        nn, p = key
        if nn == n:
            print(f"  pos {p}: {pos_hist[key]}")

    if any(k[0] == n for k in consec_hist):
        print(f"Bad vertices consecutive in HP (#bad=3):")
        for key in sorted(consec_hist.keys()):
            nn, c = key
            if nn == n:
                print(f"  consecutive={c}: {consec_hist[key]}")

    if any(k[0] == n for k in gap_hist):
        print(f"Gap patterns between bad vertices (#bad=3):")
        for key in sorted(gap_hist.keys()):
            nn, g = key
            if nn == n:
                print(f"  gaps={g}: {gap_hist[key]}")


# ============================================================
# P-EXTENSION LIFT TEST
# ============================================================
print("\n" + "=" * 60)
print("P-EXTENSION LIFT: STRUCTURE OF HIDDEN CYCLES ALONG HP")
print("=" * 60)

n = 5
ne = n*(n-1)//2
example_count = 0

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

    hp = find_hamiltonian_path(A, n)
    edges_T = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    edge_idx_T = {e: i for i, e in enumerate(edges_T)}

    example_count += 1
    if example_count > 2:
        break

    print(f"\nTournament bits={bits}, bad={bad}")
    print(f"HP: {' → '.join(str(v) for v in hp)}")
    hp_pos = {v: i for i, v in enumerate(hp)}
    print(f"Bad positions in HP: {[hp_pos[v] for v in bad]}")

    # For each bad vertex, show the hidden cycle values on HP edges
    hp_edges = [(hp[i], hp[i+1]) for i in range(n-1)]
    print(f"HP edges: {hp_edges}")

    lifts = []
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
        lifts.append(lift)

        # Show values on HP edges
        hp_vals = [lift[edge_idx_T[e]] if e in edge_idx_T else 0 for e in hp_edges]
        print(f"  z_{v} on HP edges: {[f'{x:.4f}' for x in hp_vals]}")

    # Check: alternating sum of 3 lifts
    alt_sum = lifts[0] - lifts[1] + lifts[2]
    print(f"\n  Alternating sum z_{bad[0]} - z_{bad[1]} + z_{bad[2]}:")
    for i, e in enumerate(edges_T):
        if abs(alt_sum[i]) > 1e-10:
            print(f"    ({e}): {alt_sum[i]:.6f}")
    print(f"  Norm of alternating sum: {np.linalg.norm(alt_sum):.6f}")

    # Check: is alternating sum a coboundary of T?
    coboundary_mat = np.zeros((len(edges_T), n), dtype=float)
    for v in range(n):
        for i, (a, b) in enumerate(edges_T):
            if b == v: coboundary_mat[i, v] += 1
            if a == v: coboundary_mat[i, v] -= 1

    # Solve: alt_sum = coboundary_mat @ f?
    f, residuals, _, _ = np.linalg.lstsq(coboundary_mat, alt_sum, rcond=None)
    is_coboundary = np.linalg.norm(coboundary_mat @ f - alt_sum) < 1e-8
    print(f"  Alternating sum is coboundary: {is_coboundary}")

    # Check: is alt_sum in im(∂₂)?
    tts = get_tts(A, n)
    mat = boundary2_matrix(A, n, edges_T, tts)
    coeff, _, _, _ = np.linalg.lstsq(mat, alt_sum, rcond=None)
    is_boundary = np.linalg.norm(mat @ coeff - alt_sum) < 1e-8
    print(f"  Alternating sum is in im(∂₂): {is_boundary}")

    # Check pairwise sums
    print(f"\n  Pairwise sums:")
    for i in range(3):
        for j in range(i+1, 3):
            s = lifts[i] + lifts[j]
            f2, _, _, _ = np.linalg.lstsq(coboundary_mat, s, rcond=None)
            is_cob = np.linalg.norm(coboundary_mat @ f2 - s) < 1e-8
            print(f"    z_{bad[i]} + z_{bad[j]}: norm={np.linalg.norm(s):.4f}, coboundary={is_cob}")


# ============================================================
# KEY TEST: Can we construct a HYPOTHETICAL 4th bad vertex?
# ============================================================
print("\n" + "=" * 60)
print("HYPOTHETICAL 4th BAD VERTEX: DEPENDENCY TEST")
print("=" * 60)

print("""
Since #bad is always ≤ 3, we can't directly test with 4 bad vertices.
Instead: for each tournament with #bad=3, check whether adding ANY
non-bad vertex's hidden cycle would create a dependency.

For non-bad vertex w (β₁(T\\w)=0): the "hidden cycle" z_w = 0.
So trivially, {z_v1, z_v2, z_v3, z_w=0} has rank 3 ≤ 3.

The REAL question: is there a structural reason why β₁(T\\w) CAN'T be 1
for a 4th vertex w when 3 vertices are already bad?

Approach: for #bad=3 tournaments, check the RANK of ∂₂ restricted to
T\\w for each non-bad w. How close is it to having β₁=1?
""")

n = 5
ne = n*(n-1)//2
margin_hist = Counter()
total = 0

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
        b1v = beta1(A_sub, n-1)
        if b1v == 1:
            bad.append(v)

    if len(bad) != 3:
        continue
    total += 1

    # For each non-bad vertex, check "margin to being bad"
    for w in range(n):
        if w in bad:
            continue
        remaining = [i for i in range(n) if i != w]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]

        # β₁(T\w) = 0. How many TTs would need to be removed to make β₁=1?
        edges_sub = [(i,j) for i in range(n-1) for j in range(n-1) if i!=j and A_sub[i][j]==1]
        tts_sub = get_tts(A_sub, n-1)
        mat_sub = boundary2_matrix(A_sub, n-1, edges_sub, tts_sub)
        rank_sub = np.linalg.matrix_rank(mat_sub, tol=1e-8)
        redundancy = len(tts_sub) - rank_sub
        # β₁ = C(n-1,2) - (n-2) - rank = 6 - 3 - rank
        # For β₁=1 need rank = 6-3-1 = 2 (at n-1=4)
        # Current rank = 3 (since β₁=0 means rank = C(3,2)-(4-1)+1-0 = 3-3+1 = hm wait
        # dim Z₁ = C(4,2) - (4-1) = 6-3 = 3
        # β₁ = dim Z₁ - rank = 3 - rank
        # β₁ = 0 → rank = 3
        # For β₁ = 1 → rank = 2, so need to DROP rank by 1
        margin_hist[redundancy] += 1

print(f"\nTotal #bad=3 tournaments: {total}")
print(f"Redundancy of T\\w for non-bad w (margin to β₁=1):")
print(f"  (redundancy = #TTs - rank in T\\w)")
for key in sorted(margin_hist.keys()):
    print(f"  redundancy={key}: {margin_hist[key]} instances")

print(f"\n  redundancy=0 means ALL TTs of T\\w are rank-critical")
print(f"  Higher redundancy = harder to create β₁=1 by removing 1 TT")


# ============================================================
# n=6: SAME TEST — margin for non-bad vertices
# ============================================================
print("\n--- n=6 (first 500) ---")
n = 6
ne = n*(n-1)//2
margin_hist6 = Counter()
total6 = 0

for bits in range(1 << ne):
    if total6 >= 500:
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
        continue
    total6 += 1

    for w in range(n):
        if w in bad:
            continue
        remaining = [i for i in range(n) if i != w]
        A_sub = np.zeros((n-1,n-1), dtype=int)
        for i2,i in enumerate(remaining):
            for j2,j in enumerate(remaining):
                A_sub[i2][j2] = A[i][j]
        edges_sub = [(i,j) for i in range(n-1) for j in range(n-1) if i!=j and A_sub[i][j]==1]
        tts_sub = get_tts(A_sub, n-1)
        mat_sub = boundary2_matrix(A_sub, n-1, edges_sub, tts_sub)
        rank_sub = np.linalg.matrix_rank(mat_sub, tol=1e-8)
        redundancy = len(tts_sub) - rank_sub
        margin_hist6[redundancy] += 1

print(f"Total #bad=3 tournaments: {total6}")
print(f"Redundancy of T\\w for non-bad w:")
for key in sorted(margin_hist6.keys()):
    print(f"  redundancy={key}: {margin_hist6[key]} instances")
