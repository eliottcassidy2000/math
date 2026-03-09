"""
Test Grok's algebraic proof of rank-criticality of bad-vertex TT.

Idea: Project ∂₂|Ω₂ columns onto V = span{bad-edges} (3D subspace).
- Bad-TT → full support vector r = e_xy + e_yz - e_xz
- Mixed TTs (2 bad verts) → ±single e_ij
- Good TTs (0-1 bad verts) → 0
Claim: Mixed TTs span ≤ 2D hyperplane in V (star-cocycle constraint),
  and r lies outside → bad-TT column independent of others → RC always.

opus-2026-03-09-S51 (testing Grok's proposal)
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

def test_grok_proof(A, n, verbose=False):
    """Test Grok's projection argument for a single tournament."""
    bad = get_bad_vertices(A, n)
    if len(bad) != 3:
        return None

    b1 = beta1(A, n)
    if b1 != 0:
        return None

    x, y, z = bad
    edges = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]==1]
    tts = get_tts(A, n)
    edge_idx = {e: i for i, e in enumerate(edges)}

    # Find the bad-vertex TT (the one with all 3 bad verts, transitive)
    bad_tt_idx = None
    bad_tt = None
    for idx_t, tt in enumerate(tts):
        if set(tt) == set(bad):
            bad_tt_idx = idx_t
            bad_tt = tt
            break

    if bad_tt_idx is None:
        return None  # shouldn't happen

    # Identify the 3 bad-edges (directed edges among bad vertices)
    bad_edges = []
    for i in bad:
        for j in bad:
            if i != j and A[i][j] == 1:
                bad_edges.append((i,j))
    # Should be exactly 3 edges (tournament on 3 vertices)
    assert len(bad_edges) == 3, f"Expected 3 bad edges, got {len(bad_edges)}"

    bad_edge_idx = {e: i for i, e in enumerate(bad_edges)}

    # Project each TT column onto V = span{bad_edges}
    # ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)
    # Projection: keep only terms that are bad-edges
    projections = []
    tt_types = []  # 'bad', 'mixed', 'good'
    for tt in tts:
        a, b, c = tt
        bad_count = sum(1 for v in tt if v in bad)

        proj = np.zeros(3)
        # boundary terms
        terms = [(b,c,1), (a,c,-1), (a,b,1)]
        for (u,v,sign) in terms:
            if (u,v) in bad_edge_idx:
                proj[bad_edge_idx[(u,v)]] += sign

        projections.append(proj)
        if bad_count == 3:
            tt_types.append('bad')
        elif bad_count == 2:
            tt_types.append('mixed')
        else:
            tt_types.append('good')

    # Analyze projections
    bad_proj = None
    mixed_projs = []
    good_projs = []

    for i, (proj, tp) in enumerate(zip(projections, tt_types)):
        if tp == 'bad':
            bad_proj = proj
        elif tp == 'mixed':
            mixed_projs.append(proj)
        else:
            good_projs.append(proj)

    # Check Grok's claims:
    # 1. Bad-TT has full support (all 3 components nonzero)
    bad_support = np.sum(np.abs(bad_proj) > 1e-10)

    # 2. Good TTs project to 0
    good_zero = all(np.max(np.abs(p)) < 1e-10 for p in good_projs) if good_projs else True

    # 3. Mixed TTs project to single-edge vectors (support = 1)
    mixed_supports = [np.sum(np.abs(p) > 1e-10) for p in mixed_projs]
    mixed_all_single = all(s <= 1 for s in mixed_supports)

    # 4. Mixed TTs span ≤ 2D
    if mixed_projs:
        mixed_mat = np.array(mixed_projs).T  # 3 × #mixed
        mixed_rank = np.linalg.matrix_rank(mixed_mat, tol=1e-8)
    else:
        mixed_rank = 0

    # 5. Bad proj outside mixed span?
    if mixed_projs:
        augmented = np.column_stack([mixed_mat, bad_proj.reshape(3,1)])
        aug_rank = np.linalg.matrix_rank(augmented, tol=1e-8)
        bad_outside = aug_rank > mixed_rank
    else:
        bad_outside = True

    result = {
        'bad_support': bad_support,
        'good_zero': good_zero,
        'mixed_all_single': mixed_all_single,
        'mixed_rank': mixed_rank,
        'bad_outside': bad_outside,
        'n_mixed': len(mixed_projs),
        'n_good': len(good_projs),
        'bad_proj': bad_proj,
    }

    if verbose:
        print(f"  Bad vertices: {bad}, Bad-TT: {bad_tt}")
        print(f"  Bad edges: {bad_edges}")
        print(f"  Bad-TT proj: {bad_proj} (support={bad_support})")
        print(f"  Good TTs all zero: {good_zero} ({len(good_projs)} good TTs)")
        print(f"  Mixed TTs all single-edge: {mixed_all_single} ({len(mixed_projs)} mixed TTs)")
        print(f"  Mixed rank: {mixed_rank}")
        print(f"  Bad proj outside mixed span: {bad_outside}")
        for i, (proj, tp) in enumerate(zip(projections, tt_types)):
            if tp == 'mixed' and np.max(np.abs(proj)) > 1e-10:
                print(f"    Mixed TT {tts[i]}: proj={proj}")

    return result


# ============================================================
# EXHAUSTIVE TEST
# ============================================================
for n in [5, 6]:
    print(f"\n{'='*60}")
    print(f"n={n}: EXHAUSTIVE TEST OF GROK'S RC PROOF")
    print(f"{'='*60}")

    ne = n*(n-1)//2
    total_tested = 0
    counts = Counter()
    failures = []
    first_verbose = True

    for bits in range(1 << ne):
        A = tournament_from_bits(n, bits)
        result = test_grok_proof(A, n, verbose=(first_verbose and n <= 6))
        if result is None:
            continue

        total_tested += 1
        if first_verbose and n <= 6:
            first_verbose = False

        # Track all claims
        if not result['good_zero']:
            counts['good_nonzero'] += 1
        if not result['mixed_all_single']:
            counts['mixed_not_single'] += 1
        if result['mixed_rank'] > 2:
            counts['mixed_rank_3'] += 1
        if not result['bad_outside']:
            counts['bad_inside_span'] += 1
            failures.append(bits)
        if result['bad_support'] != 3:
            counts['bad_not_full'] += 1

        counts[f'mixed_rank_{result["mixed_rank"]}'] += 1
        counts[f'bad_support_{result["bad_support"]}'] += 1

    print(f"\nTotal tournaments with β₁=0, #bad=3: {total_tested}")
    print(f"\nGrok's claims:")
    print(f"  1. Bad-TT has full support (3): {counts.get('bad_not_full', 0)} violations")
    print(f"  2. Good TTs project to 0: {counts.get('good_nonzero', 0)} violations")
    print(f"  3. Mixed TTs are single-edge: {counts.get('mixed_not_single', 0)} violations")
    print(f"  4. Mixed rank ≤ 2: {counts.get('mixed_rank_3', 0)} violations")
    print(f"  5. Bad proj outside mixed span: {counts.get('bad_inside_span', 0)} violations")
    print(f"\nMixed rank distribution: " +
          ", ".join(f"rank={k.split('_')[-1]}: {v}" for k,v in sorted(counts.items()) if k.startswith('mixed_rank_')))
    print(f"Bad support distribution: " +
          ", ".join(f"support={k.split('_')[-1]}: {v}" for k,v in sorted(counts.items()) if k.startswith('bad_support_')))

    if failures:
        print(f"\n*** FAILURES at bits: {failures[:5]} ***")
        A = tournament_from_bits(n, failures[0])
        test_grok_proof(A, n, verbose=True)
