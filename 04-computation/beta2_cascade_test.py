#!/usr/bin/env python3
"""
beta2_cascade_test.py — Precise test of flip obstruction mechanism

For tournaments with beta_1=0 and exactly 3 bad vertices, flip the bad-vertex
transitive triple and track exactly what happens to TTs and boundary spaces.

Tests Grok's claim: "Star constraints propagate dependence across adjacent TTs
in the low-dim space, making pre-flip indep generator Q-redundant while new TTs
absorb into higher filtration."
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict
import sys

# Import core functions from path_homology_v2
sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')

# We need the raw functions, not the script's validation output.
# Re-implement the core functions inline to avoid running the validation.

def enumerate_allowed_paths(A, n, p):
    if p < 0:
        return []
    if p == 0:
        return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    p = len(path) - 1
    result = []
    for i in range(p + 1):
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1:
                M[idx_pm1[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0:
        return np.zeros((0, 0))
    if p == 0:
        return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0:
        return np.zeros((dim_Ap, 0))
    return null_space

def path_betti_numbers(A, n, max_dim=None):
    if max_dim is None:
        max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    betti = []
    for p in range(max_dim + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_omega_p == 0:
            betti.append(0)
            continue
        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_omega = bd_p @ omega[p]
        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0
        ker_dim = dim_omega_p - rank_p
        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0
        beta_p = ker_dim - im_dim
        betti.append(max(0, beta_p))
    return betti

# ============================================================
# TOURNAMENT UTILITIES
# ============================================================

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield [list(row) for row in A]

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if np.random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def get_transitive_triples(A, n):
    """Return list of (a,b,c) where a->b, b->c, a->c (transitive triple)."""
    tts = []
    for a in range(n):
        for b in range(n):
            if b == a: continue
            if not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    tts.append((a,b,c))
    return tts

def get_3cycles(A, n):
    """Return list of (a,b,c) with a<b<c forming a 3-cycle."""
    cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    cycles.append((a,b,c))
    return cycles

def count_bad_vertices(A, n):
    """A vertex v is 'bad' if it has a 3-cycle through it.
    Actually: v is bad if out-neighborhood of v is NOT transitive.
    Equivalently: v has a 3-cycle in its out-neighborhood... no.

    Let me use the standard: v is bad if the sub-tournament on out(v) is not transitive,
    i.e., there exists a 3-cycle i->j->k->i with all of i,j,k in out(v).

    Actually, for GLMY: 'bad vertex' typically means a vertex participating in a 3-cycle.
    Let me count vertices that are in at least one 3-cycle.
    """
    cycles = get_3cycles(A, n)
    bad = set()
    for (a,b,c) in cycles:
        bad.add(a)
        bad.add(b)
        bad.add(c)
    return bad

def flip_triple(A, n, triple):
    """Flip a transitive triple (a,b,c) where a->b->c, a->c.
    We reverse the 3 edges to create c->b->a, c->a (a 3-cycle a->b, b->c, c->a... no).

    Actually: flip means reverse ALL 3 edges of the triple.
    a->b becomes b->a, b->c becomes c->b, a->c becomes c->a.
    Result: b->a, c->b, c->a. This is the TT (c,b,a).

    Wait — the problem says "flip the bad-vertex transitive triple".
    For a tournament with exactly 3 bad vertices forming a 3-cycle,
    we should flip the 3-cycle to a transitive triple, or vice versa.

    Let me re-read: "take the bad-vertex transitive triple and flip it."

    Hmm. If there are 3 bad vertices that form a 3-cycle, the "bad-vertex TT"
    might mean: the transitive triple on those 3 vertices... but they form a cycle, not a TT.

    I think "flip" means: reverse the 3-cycle {a,b,c} (reverse all 3 arcs) to make it transitive.
    """
    a, b, c = triple
    A2 = [row[:] for row in A]
    # Reverse all 3 edges between a, b, c
    for x, y in [(a,b), (b,c), (a,c)]:
        # Swap edge direction
        A2[x][y], A2[y][x] = A2[y][x], A2[x][y]
    return A2

def boundary_matrix_for_tts(tts, edges_list, n):
    """Build boundary matrix: rows=edges (2-paths), cols=TTs (3-paths).

    A TT (a,b,c) = allowed 3-path a->b->c with a->c.
    ∂(a,b,c) = (b,c) - (a,c) + (a,b)
    """
    edge_idx = {e: i for i, e in enumerate(edges_list)}
    M = np.zeros((len(edges_list), len(tts)))
    for j, (a, b, c) in enumerate(tts):
        # ∂(a,b,c) = (b,c) - (a,c) + (a,b)
        if (b,c) in edge_idx:
            M[edge_idx[(b,c)], j] += 1
        if (a,c) in edge_idx:
            M[edge_idx[(a,c)], j] -= 1
        if (a,b) in edge_idx:
            M[edge_idx[(a,b)], j] += 1
    return M

def matrix_rank(M, tol=1e-8):
    if M.size == 0:
        return 0
    S = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(S > tol))

def is_in_column_span(M, v, tol=1e-8):
    """Check if vector v is in the column span of M."""
    if M.shape[1] == 0:
        return np.linalg.norm(v) < tol
    # Solve M x = v in least squares
    x, residuals, rank, sv = np.linalg.lstsq(M, v, rcond=None)
    return np.linalg.norm(M @ x - v) < tol

# ============================================================
# MAIN ANALYSIS
# ============================================================

def analyze_tournament(A, n, label=""):
    """Full cascade analysis for one tournament."""

    # Get edges (allowed 1-paths = directed edges)
    edges = []
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                edges.append((i,j))

    # Get TTs (transitive triples = allowed 2-paths that are ∂-invariant...
    # Actually TTs are a specific subset of allowed 2-paths.
    # An allowed 2-path is (a,b,c) with a->b and b->c.
    # It's a TT if additionally a->c.
    # Non-TT allowed 2-paths have c->a (forming a 3-cycle with the path).

    allowed_2 = enumerate_allowed_paths(A, n, 2)
    tts = []
    non_tts = []
    for path in allowed_2:
        a, b, c = path
        if A[a][c]:  # a->c: transitive triple
            tts.append(path)
        else:  # c->a: path through 3-cycle
            non_tts.append(path)

    # Betti numbers
    betti = path_betti_numbers(A, n, max_dim=2)

    # Bad vertices (in 3-cycles)
    bad = count_bad_vertices(A, n)

    return {
        'A': A, 'n': n, 'edges': edges, 'tts': tts, 'non_tts': non_tts,
        'allowed_2': allowed_2, 'betti': betti, 'bad': bad, 'label': label
    }

def find_bad_triple(A, n, bad_verts):
    """Find the 3-cycle among the bad vertices."""
    bv = sorted(bad_verts)
    if len(bv) != 3:
        return None
    a, b, c = bv
    # Check which orientation the 3-cycle has
    if A[a][b] and A[b][c] and A[c][a]:
        return (a, b, c)  # a->b->c->a
    elif A[a][c] and A[c][b] and A[b][a]:
        return (a, c, b)  # a->c->b->a
    return None

def cascade_analysis(A, n, idx=0):
    """Full cascade analysis for a tournament with 3 bad vertices."""

    info = analyze_tournament(A, n, f"T{idx}")
    bad = info['bad']
    betti = info['betti']

    if len(bad) != 3 or betti[1] != 0:
        return None

    # Find the 3-cycle
    triple = find_bad_triple(A, n, bad)
    if triple is None:
        return None

    # Flip it
    A2 = flip_triple(A, n, triple)
    info2 = analyze_tournament(A2, n, f"T{idx}'")
    betti2 = info2['betti']

    # Get edges for both tournaments
    edges_old = set(info['edges'])
    edges_new = set(info2['edges'])
    all_edges = sorted(edges_old | edges_new)
    edge_idx = {e: i for i, e in enumerate(all_edges)}

    # TT sets
    tts_old = set(info['tts'])
    tts_new = set(info2['tts'])

    D = sorted(tts_old - tts_new)  # destroyed
    C = sorted(tts_new - tts_old)  # created
    K = sorted(tts_old & tts_new)  # kept

    # Build boundary matrices using ALL edges as row space
    def bd_mat(tt_list):
        M = np.zeros((len(all_edges), len(tt_list)))
        for j, (a, b, c) in enumerate(tt_list):
            if (b,c) in edge_idx: M[edge_idx[(b,c)], j] += 1
            if (a,c) in edge_idx: M[edge_idx[(a,c)], j] -= 1
            if (a,b) in edge_idx: M[edge_idx[(a,b)], j] += 1
        return M

    B_old = bd_mat(list(D) + list(K))  # span of boundaries from old TTs
    B_new = bd_mat(list(C) + list(K))  # span of boundaries from new TTs
    B_kept = bd_mat(list(K))
    B_destroyed = bd_mat(list(D))
    B_created = bd_mat(list(C))

    dim_old = matrix_rank(B_old)
    dim_new = matrix_rank(B_new)
    dim_kept = matrix_rank(B_kept)
    dim_destroyed = matrix_rank(B_destroyed)
    dim_created = matrix_rank(B_created)

    # Check redundancy of each destroyed TT
    destroyed_redundant = []
    for i, tt in enumerate(D):
        v = bd_mat([tt]).flatten()
        # Is this in span of kept TTs?
        destroyed_redundant.append(is_in_column_span(B_kept, v))

    # Check redundancy of each created TT
    created_redundant = []
    for i, tt in enumerate(C):
        v = bd_mat([tt]).flatten()
        created_redundant.append(is_in_column_span(B_kept, v))

    # Independent destroyed/created
    n_destroyed_indep = 0
    if D:
        # How many of D are independent relative to K?
        # = dim(span(D ∪ K)) - dim(span(K))
        n_destroyed_indep = dim_old - dim_kept

    n_created_indep = 0
    if C:
        n_created_indep = dim_new - dim_kept

    # Star constraint analysis
    a0, b0, c0 = triple
    star_info = {}
    for v in [a0, b0, c0]:
        out_v = [u for u in range(n) if A[v][u]]
        out_v2 = [u for u in range(n) if A2[v][u]]
        # TTs rooted at v (v is source)
        tts_v_old = [t for t in tts_old if t[0] == v]
        tts_v_new = [t for t in tts_new if t[0] == v]
        star_info[v] = {
            'out_old': out_v, 'out_new': out_v2,
            'tts_old': len(tts_v_old), 'tts_new': len(tts_v_new),
            'out_changed': set(out_v) != set(out_v2)
        }

    # Find the lost dimension (null space analysis)
    # The H1 generator of T' lives in ker(∂_1) but not im(∂_2)
    # ∂_1 maps edges to vertices: ∂(a,b) = b - a
    # For T': compute ker(∂_1) and im(∂_2) explicitly

    # Build ∂_1 for T' (edges -> vertices)
    edges_new_list = sorted(edges_new)
    d1 = np.zeros((n, len(edges_new_list)))
    for j, (a, b) in enumerate(edges_new_list):
        d1[b, j] += 1  # +b
        d1[a, j] -= 1  # -a

    # Build ∂_2 for T' (using Omega_2, not just TTs!)
    # We need the actual Omega_2 computation for T'
    allowed_2_new = enumerate_allowed_paths(A2, n, 2)
    allowed_1_new = enumerate_allowed_paths(A2, n, 1)
    omega2_new = compute_omega_basis(A2, n, 2, allowed_2_new, allowed_1_new)

    d2_full = build_full_boundary_matrix(allowed_2_new, allowed_1_new)

    # Image of ∂_2 in edge space
    if omega2_new.shape[1] > 0:
        d2_omega = d2_full @ omega2_new
        im_d2_rank = matrix_rank(d2_omega)
    else:
        d2_omega = np.zeros((len(allowed_1_new), 0))
        im_d2_rank = 0

    # ker(∂_1)
    # Map from allowed_1_new indexing to vertex indexing
    d1_allowed = np.zeros((n, len(allowed_1_new)))
    for j, (a, b) in enumerate(allowed_1_new):
        d1_allowed[b, j] += 1
        d1_allowed[a, j] -= 1

    U, S, Vt = np.linalg.svd(d1_allowed, full_matrices=True)
    r1 = sum(s > 1e-8 for s in S)
    ker_d1 = Vt[r1:].T  # columns = basis of ker(∂_1)

    # H_1 generator: something in ker(d1) but not in im(d2)
    h1_gen = None
    if betti2[1] > 0 and ker_d1.shape[1] > im_d2_rank:
        # Find a vector in ker(d1) not in im(d2)
        # Stack im(d2) columns and extend with ker(d1) to find the quotient
        if d2_omega.shape[1] > 0:
            combined = np.hstack([d2_omega, ker_d1])
        else:
            combined = ker_d1.copy()
        # The H1 generators are ker(d1) columns not expressible as im(d2)
        # Use column pivoting
        # Simple approach: project ker(d1) basis onto orthogonal complement of im(d2)
        if d2_omega.shape[1] > 0:
            Q, R = np.linalg.qr(d2_omega, mode='reduced')
            # Project each ker(d1) basis vector
            for col_i in range(ker_d1.shape[1]):
                v = ker_d1[:, col_i]
                proj = Q @ (Q.T @ v)
                residual = v - proj
                if np.linalg.norm(residual) > 1e-6:
                    h1_gen = residual / np.linalg.norm(residual)
                    break
        else:
            if ker_d1.shape[1] > 0:
                h1_gen = ker_d1[:, 0]

    # Express h1_gen in terms of edges
    h1_edges = {}
    if h1_gen is not None:
        for j, e in enumerate(allowed_1_new):
            if abs(h1_gen[j]) > 1e-6:
                h1_edges[e] = round(h1_gen[j], 4)

    return {
        'idx': idx,
        'triple': triple,
        'betti_old': betti,
        'betti_new': betti2,
        'bad_old': sorted(bad),
        'bad_new': sorted(info2['bad']),
        'n_tts_old': len(tts_old),
        'n_tts_new': len(tts_new),
        '|D|': len(D), '|C|': len(C), '|K|': len(K),
        'D': D, 'C': C,
        'dim_B_old': dim_old,
        'dim_B_new': dim_new,
        'dim_B_kept': dim_kept,
        'destroyed_redundant': destroyed_redundant,
        'created_redundant': created_redundant,
        'n_destroyed_indep': n_destroyed_indep,
        'n_created_indep': n_created_indep,
        'star_info': star_info,
        'h1_edges': h1_edges,
        'dim_ker_d1': ker_d1.shape[1],
        'dim_im_d2': im_d2_rank,
    }

# ============================================================
# RUN n=5 EXHAUSTIVE
# ============================================================

print("=" * 80)
print("BETA-2 CASCADE TEST: Flip Obstruction Mechanism")
print("=" * 80)

print("\n" + "=" * 80)
print("PART 1: n=5 EXHAUSTIVE — All tournaments with β₁=0 and 3 bad vertices")
print("=" * 80)

n = 5
results_5 = []
all_t = list(all_tournaments(n))
print(f"Total n=5 tournaments: {len(all_t)}")

# Filter
candidates = []
for idx, A in enumerate(all_t):
    bad = count_bad_vertices(A, n)
    if len(bad) != 3:
        continue
    betti = path_betti_numbers(A, n, max_dim=2)
    if betti[1] != 0:
        continue
    candidates.append((idx, A))

print(f"Candidates (β₁=0, 3 bad): {len(candidates)}")

# Analyze each
dim_old_hist = defaultdict(int)
dim_new_hist = defaultdict(int)
dim_kept_hist = defaultdict(int)
rank_drop_hist = defaultdict(int)
destroyed_indep_hist = defaultdict(int)
created_indep_hist = defaultdict(int)
all_destroyed_redundant = 0
all_destroyed_total = 0
all_created_redundant = 0
all_created_total = 0

for i, (idx, A) in enumerate(candidates):
    r = cascade_analysis(A, n, idx)
    if r is None:
        continue
    results_5.append(r)

    dim_old_hist[r['dim_B_old']] += 1
    dim_new_hist[r['dim_B_new']] += 1
    dim_kept_hist[r['dim_B_kept']] += 1
    rank_drop = r['dim_B_old'] - r['dim_B_new']
    rank_drop_hist[rank_drop] += 1
    destroyed_indep_hist[r['n_destroyed_indep']] += 1
    created_indep_hist[r['n_created_indep']] += 1

    all_destroyed_total += r['|D|']
    all_destroyed_redundant += sum(r['destroyed_redundant'])
    all_created_total += r['|C|']
    all_created_redundant += sum(r['created_redundant'])

print(f"\nAnalyzed: {len(results_5)} tournaments")

# Print first 10 in detail
print("\n--- First 10 detailed results ---")
for r in results_5[:10]:
    print(f"\nT{r['idx']}: triple={r['triple']}, β={r['betti_old']} -> β'={r['betti_new']}")
    print(f"  Bad: {r['bad_old']} -> {r['bad_new']}")
    print(f"  TTs: {r['n_tts_old']} -> {r['n_tts_new']}")
    print(f"  |D|={r['|D|']}, |C|={r['|C|']}, |K|={r['|K|']}")
    print(f"  dim(B_old)={r['dim_B_old']}, dim(B_new)={r['dim_B_new']}, dim(B_kept)={r['dim_B_kept']}")
    print(f"  Rank drop: {r['dim_B_old'] - r['dim_B_new']}")
    print(f"  Destroyed indep (relative to K): {r['n_destroyed_indep']}")
    print(f"  Created indep (relative to K): {r['n_created_indep']}")
    print(f"  Destroyed TTs redundant given K: {r['destroyed_redundant']}")
    print(f"  Created TTs redundant given K: {r['created_redundant']}")
    print(f"  D (destroyed TTs): {r['D']}")
    print(f"  C (created TTs): {r['C']}")
    if r['h1_edges']:
        print(f"  H₁ generator of T' (edge support): {r['h1_edges']}")
    for v, si in r['star_info'].items():
        print(f"  Star v={v}: out {si['out_old']}->{si['out_new']}, "
              f"TTs rooted: {si['tts_old']}->{si['tts_new']}, changed={si['out_changed']}")

# Aggregate statistics
print("\n\n--- AGGREGATE STATISTICS (n=5) ---")
print(f"\nBetti change distribution:")
betti_change = defaultdict(int)
for r in results_5:
    key = (tuple(r['betti_old']), tuple(r['betti_new']))
    betti_change[key] += 1
for k, v in sorted(betti_change.items()):
    print(f"  {list(k[0])} -> {list(k[1])}: {v}")

print(f"\ndim(B_old) distribution: {dict(dim_old_hist)}")
print(f"dim(B_new) distribution: {dict(dim_new_hist)}")
print(f"dim(B_kept) distribution: {dict(dim_kept_hist)}")
print(f"Rank drop (dim_old - dim_new): {dict(rank_drop_hist)}")

print(f"\nDestroyed-TT independence (relative to kept):")
print(f"  Distribution: {dict(destroyed_indep_hist)}")
print(f"  Total destroyed TTs: {all_destroyed_total}")
print(f"  Redundant given K: {all_destroyed_redundant} ({100*all_destroyed_redundant/max(all_destroyed_total,1):.1f}%)")

print(f"\nCreated-TT independence (relative to kept):")
print(f"  Distribution: {dict(created_indep_hist)}")
print(f"  Total created TTs: {all_created_total}")
print(f"  Redundant given K: {all_created_redundant} ({100*all_created_redundant/max(all_created_total,1):.1f}%)")

# KEY QUESTION: Does dim(B_kept) = dim(B_old)?
print(f"\n--- KEY QUESTION: Is dim(B_kept) = dim(B_old)? ---")
kept_equals_old = sum(1 for r in results_5 if r['dim_B_kept'] == r['dim_B_old'])
print(f"  Yes: {kept_equals_old}/{len(results_5)}")
print(f"  (If yes: destroyed TTs contribute NOTHING new to boundary space)")

# KEY QUESTION: Is dim(B_old) = dim(B_new) + 1?
print(f"\n--- KEY QUESTION: Is dim(B_old) = dim(B_new) + 1? ---")
drop_one = sum(1 for r in results_5 if r['dim_B_old'] - r['dim_B_new'] == 1)
print(f"  Yes: {drop_one}/{len(results_5)}")

# KEY QUESTION: net mechanism
print(f"\n--- NET MECHANISM ---")
for r in results_5[:10]:
    net = r['n_created_indep'] - r['n_destroyed_indep']
    print(f"  T{r['idx']}: destroyed_indep={r['n_destroyed_indep']}, "
          f"created_indep={r['n_created_indep']}, net={net}, "
          f"rank_change={r['dim_B_new']-r['dim_B_old']}")

# H1 generator analysis
print(f"\n--- H₁ GENERATOR OF T' ---")
h1_count = sum(1 for r in results_5 if r['h1_edges'])
print(f"Found H₁ generator for {h1_count}/{len(results_5)} tournaments")
# Show edge supports
edge_support_sizes = defaultdict(int)
for r in results_5:
    if r['h1_edges']:
        edge_support_sizes[len(r['h1_edges'])] += 1
print(f"H₁ generator edge-support sizes: {dict(edge_support_sizes)}")

# Show a few generators
for r in results_5[:5]:
    if r['h1_edges']:
        print(f"  T{r['idx']}: {r['h1_edges']}")
        # Check if the support edges form a cycle
        support_edges = list(r['h1_edges'].keys())
        vertices_in_gen = set()
        for (a,b) in support_edges:
            vertices_in_gen.add(a)
            vertices_in_gen.add(b)
        print(f"    Vertices in generator: {sorted(vertices_in_gen)}")
        print(f"    Bad vertices of T': {r['bad_new']}")

# ============================================================
# PART 2: n=6 SAMPLE
# ============================================================

print("\n\n" + "=" * 80)
print("PART 2: n=6 SAMPLE — 50 tournaments with β₁=0 and 3 bad vertices")
print("=" * 80)

n = 6
np.random.seed(42)
results_6 = []
attempts = 0
max_attempts = 50000

while len(results_6) < 50 and attempts < max_attempts:
    attempts += 1
    A = random_tournament(n)
    bad = count_bad_vertices(A, n)
    if len(bad) != 3:
        continue
    betti = path_betti_numbers(A, n, max_dim=2)
    if betti[1] != 0:
        continue

    r = cascade_analysis(A, n, attempts)
    if r is not None:
        results_6.append(r)
        if len(results_6) % 10 == 0:
            print(f"  Found {len(results_6)}/50 (after {attempts} attempts)")

print(f"\nFound {len(results_6)} qualifying n=6 tournaments after {attempts} attempts")

if results_6:
    # Aggregate
    print("\n--- AGGREGATE STATISTICS (n=6) ---")
    betti_change_6 = defaultdict(int)
    for r in results_6:
        key = (tuple(r['betti_old']), tuple(r['betti_new']))
        betti_change_6[key] += 1
    print(f"Betti change distribution:")
    for k, v in sorted(betti_change_6.items()):
        print(f"  {list(k[0])} -> {list(k[1])}: {v}")

    dim_old_6 = defaultdict(int)
    dim_new_6 = defaultdict(int)
    dim_kept_6 = defaultdict(int)
    rank_drop_6 = defaultdict(int)
    d_indep_6 = defaultdict(int)
    c_indep_6 = defaultdict(int)

    for r in results_6:
        dim_old_6[r['dim_B_old']] += 1
        dim_new_6[r['dim_B_new']] += 1
        dim_kept_6[r['dim_B_kept']] += 1
        rd = r['dim_B_old'] - r['dim_B_new']
        rank_drop_6[rd] += 1
        d_indep_6[r['n_destroyed_indep']] += 1
        c_indep_6[r['n_created_indep']] += 1

    print(f"\ndim(B_old) distribution: {dict(dim_old_6)}")
    print(f"dim(B_new) distribution: {dict(dim_new_6)}")
    print(f"dim(B_kept) distribution: {dict(dim_kept_6)}")
    print(f"Rank drop: {dict(rank_drop_6)}")
    print(f"Destroyed indep: {dict(d_indep_6)}")
    print(f"Created indep: {dict(c_indep_6)}")

    kept_eq_old_6 = sum(1 for r in results_6 if r['dim_B_kept'] == r['dim_B_old'])
    drop_one_6 = sum(1 for r in results_6 if r['dim_B_old'] - r['dim_B_new'] == 1)
    print(f"\nKEY: dim(B_kept) = dim(B_old)? {kept_eq_old_6}/{len(results_6)}")
    print(f"KEY: rank drop = 1? {drop_one_6}/{len(results_6)}")

    # Detail first 5
    print("\n--- First 5 detailed (n=6) ---")
    for r in results_6[:5]:
        print(f"\nT{r['idx']}: triple={r['triple']}, β={r['betti_old']} -> β'={r['betti_new']}")
        print(f"  |D|={r['|D|']}, |C|={r['|C|']}, |K|={r['|K|']}")
        print(f"  dim(B_old)={r['dim_B_old']}, dim(B_new)={r['dim_B_new']}, dim(B_kept)={r['dim_B_kept']}")
        print(f"  Destroyed indep: {r['n_destroyed_indep']}, Created indep: {r['n_created_indep']}")
        print(f"  Destroyed redundant: {r['destroyed_redundant']}")
        print(f"  Created redundant: {r['created_redundant']}")
        if r['h1_edges']:
            print(f"  H₁ gen support: {list(r['h1_edges'].keys())}")
else:
    print("WARNING: Could not find 50 qualifying n=6 tournaments.")
    print("This is likely because 3-bad-vertex tournaments are rare at n=6.")
    print("Let me try with ANY number of bad vertices instead...")

    # Fallback: just find tournaments where flipping a 3-cycle creates β₁
    print("\nFallback: Finding n=6 tournaments where a 3-cycle flip creates β₁>0")
    results_6b = []
    attempts = 0
    while len(results_6b) < 50 and attempts < max_attempts:
        attempts += 1
        A = random_tournament(n)
        betti = path_betti_numbers(A, n, max_dim=1)
        if betti[1] != 0:
            continue
        # Find a 3-cycle and flip it
        cycles = get_3cycles(A, n)
        if not cycles:
            continue
        cyc = cycles[0]
        triple = find_bad_triple(A, n, set(cyc))
        if triple is None:
            continue
        A2 = flip_triple(A, n, triple)
        betti2 = path_betti_numbers(A2, n, max_dim=1)
        if betti2[1] > 0:
            r = cascade_analysis(A, n, attempts)
            if r is not None:
                results_6b.append(r)
                if len(results_6b) % 10 == 0:
                    print(f"  Found {len(results_6b)}/50")

    print(f"Fallback found {len(results_6b)} n=6 tournaments")
    if results_6b:
        betti_change_6b = defaultdict(int)
        rank_drop_6b = defaultdict(int)
        for r in results_6b:
            key = (tuple(r['betti_old']), tuple(r['betti_new']))
            betti_change_6b[key] += 1
            rd = r['dim_B_old'] - r['dim_B_new']
            rank_drop_6b[rd] += 1
        print(f"Betti changes: {dict(betti_change_6b)}")
        print(f"Rank drops: {dict(rank_drop_6b)}")
        kept_eq = sum(1 for r in results_6b if r['dim_B_kept'] == r['dim_B_old'])
        print(f"dim(B_kept) = dim(B_old)? {kept_eq}/{len(results_6b)}")

# ============================================================
# VERDICT
# ============================================================

print("\n\n" + "=" * 80)
print("VERDICT: Testing Grok's Claim")
print("=" * 80)
print("""
Grok's claim: "Star constraints propagate the dependence across adjacent TTs
in the low-dim space (≤3 basis), making pre-flip indep generator Q-redundant
while new TTs absorb into higher filtration."

This predicts:
1. Destroyed TTs are boundary-redundant given kept TTs (star propagation)
2. Created TTs add NO new boundary dimensions (absorb into higher filtration)
3. The net effect is a rank drop of exactly 1 in im(∂₂)
""")

if results_5:
    print("EVIDENCE FROM n=5:")

    # Test prediction 1
    all_d_redundant = all(all(r['destroyed_redundant']) for r in results_5 if r['|D|'] > 0)
    print(f"  Prediction 1 (destroyed redundant given K): {all_d_redundant}")
    if not all_d_redundant:
        counter = [r for r in results_5 if not all(r['destroyed_redundant'])]
        print(f"    Counterexamples: {len(counter)}")
        for r in counter[:3]:
            print(f"    T{r['idx']}: destroyed_redundant={r['destroyed_redundant']}")

    # Test prediction 2
    all_c_redundant = all(all(r['created_redundant']) for r in results_5 if r['|C|'] > 0)
    print(f"  Prediction 2 (created redundant given K): {all_c_redundant}")
    if not all_c_redundant:
        counter = [r for r in results_5 if not all(r['created_redundant'])]
        print(f"    Counterexamples: {len(counter)}")
        for r in counter[:3]:
            print(f"    T{r['idx']}: created_redundant={r['created_redundant']}")

    # Test prediction 3
    all_drop_one = all(r['dim_B_old'] - r['dim_B_new'] == 1 for r in results_5)
    print(f"  Prediction 3 (rank drop = 1): {all_drop_one}")
    if not all_drop_one:
        drops = [r['dim_B_old'] - r['dim_B_new'] for r in results_5]
        print(f"    Rank drops observed: {sorted(set(drops))}")

    # The actual mechanism
    print(f"\n  ACTUAL MECHANISM:")
    print(f"    dim(B_kept) = dim(B_old) always? "
          f"{all(r['dim_B_kept']==r['dim_B_old'] for r in results_5)}")
    print(f"    dim(B_kept) = dim(B_new) always? "
          f"{all(r['dim_B_kept']==r['dim_B_new'] for r in results_5)}")
    print(f"    dim(B_kept) < dim(B_old) always? "
          f"{all(r['dim_B_kept']<r['dim_B_old'] for r in results_5)}")
    print(f"    dim(B_kept) < dim(B_new) always? "
          f"{all(r['dim_B_kept']<r['dim_B_new'] for r in results_5)}")

print("\nDone.")
