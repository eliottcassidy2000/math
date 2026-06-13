#!/usr/bin/env python3
"""
Test Grok's claim about bad-vertex TTs and rank of ∂₂.

Claim: "Bad verts transitive => their TT-boundaries form basis for residual
im(d₂)^⊥ in the Q-resolution. With |bad|>=4, pigeon on star constraints +
global H₁=0 forces at least one isolated bad TT whose ∂ survives in ker(d₂^*),
injecting β₁>=1."

Tests:
1. For β₁=0 with k bad verts: TT structure among bad verts, rank-drop on removal
2. For β₁=1: bad vert count, SC sub-tournament, TT structure
3. Pigeonhole: Is the bad-vertex TT uniquely rank-critical?
4. Dual test: ker(∂₂^T) structure
5. Hypothetical 4th bad vertex analysis at n=6

opus-2026-03-09
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
import sys

# ---- Core path homology functions (from path_homology_v2.py) ----

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
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
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

# ---- Tournament utilities ----

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield [list(row) for row in A]

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = [[0]*(n-1) for _ in range(n-1)]
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B, keep

def find_bad_vertices(A, n):
    """Bad vertex v: β₁(T)=0 but β₁(T\v)=1."""
    bad = []
    for v in range(n):
        B, _ = delete_vertex(A, n, v)
        betti = path_betti_numbers(B, n-1, max_dim=1)
        if len(betti) > 1 and betti[1] >= 1:
            bad.append(v)
    return bad

def is_transitive_subtournament(A, verts):
    """Check if vertices form a transitive tournament."""
    vs = sorted(verts)
    k = len(vs)
    if k <= 1:
        return True
    # Transitive = acyclic = has a unique topological order
    # Check: no 3-cycle among the vertices
    for a, b, c in combinations(vs, 3):
        if A[a][b] and A[b][c] and A[c][a]:
            return False
        if A[a][c] and A[c][b] and A[b][a]:
            return False
    return True

def has_3cycle(A, verts):
    """Check if there's a 3-cycle among the given vertices."""
    vs = sorted(verts)
    for a, b, c in combinations(vs, 3):
        if A[a][b] and A[b][c] and A[c][a]:
            return True
        if A[a][c] and A[c][b] and A[b][a]:
            return True
    return False

def find_transitive_triples(A, verts):
    """Find all transitive triples (a,b,c) with a->b->c and a->c among verts."""
    tts = []
    for a, b, c in permutations(verts, 3):
        if A[a][b] and A[b][c] and A[a][c]:
            # Canonical: sorted source first
            if a < c:  # avoid counting (a,b,c) and (c,b,a) etc
                tts.append((a, b, c))
    # Deduplicate: each TT {a,b,c} with a->b->c, a->c appears once
    # Actually let's be more careful: a transitive triple on {x,y,z} has
    # exactly one vertex that beats both others. Let's find unique triples.
    seen = set()
    unique_tts = []
    for a, b, c in tts:
        key = frozenset([a, b, c])
        if key not in seen:
            seen.add(key)
            unique_tts.append((a, b, c))  # a->b->c, a->c
    return unique_tts

def get_omega2_data(A, n):
    """Get Ω₂ basis and boundary matrix ∂₂: Ω₂ → Ω₁."""
    allowed = {}
    for p in range(-1, 4):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)

    omega2 = compute_omega_basis(A, n, 2, allowed[2], allowed[1])
    omega1 = compute_omega_basis(A, n, 1, allowed[1], allowed[0])

    bd2_full = build_full_boundary_matrix(allowed[2], allowed[1])

    # ∂₂ restricted to Ω₂, mapping to A₁
    bd2_omega = bd2_full @ omega2

    return allowed, omega2, omega1, bd2_full, bd2_omega

def rank_of(M, tol=1e-8):
    if M.shape[0] == 0 or M.shape[1] == 0:
        return 0
    sv = np.linalg.svd(M, compute_uv=False)
    return int(sum(s > tol for s in sv))

# ========== MAIN ANALYSIS ==========

def analyze_tournament(A, n, verbose=False):
    """Full analysis of one tournament."""
    betti = path_betti_numbers(A, n, max_dim=2)
    b1 = betti[1] if len(betti) > 1 else 0

    bad = find_bad_vertices(A, n)
    k_bad = len(bad)

    is_trans_bad = is_transitive_subtournament(A, bad) if k_bad >= 2 else True
    has_3cyc_bad = has_3cycle(A, bad) if k_bad >= 3 else False

    tts_bad = find_transitive_triples(A, bad) if k_bad >= 3 else []

    return {
        'betti': betti,
        'b1': b1,
        'bad': bad,
        'k_bad': k_bad,
        'is_trans_bad': is_trans_bad,
        'has_3cyc_bad': has_3cyc_bad,
        'tts_bad': tts_bad,
    }

def rank_drop_test(A, n, info):
    """Test 1 & 3: Does removing a bad-TT from Ω₂ drop rank(∂₂)?"""
    allowed, omega2, omega1, bd2_full, bd2_omega = get_omega2_data(A, n)

    full_rank = rank_of(bd2_omega)
    dim_omega2 = omega2.shape[1] if omega2.ndim == 2 else 0

    allowed2 = allowed[2]
    allowed1 = allowed[1]
    idx2 = {path: i for i, path in enumerate(allowed2)}

    results = {}

    # For each TT among bad vertices, check if it's in Ω₂ and if removing it drops rank
    for tt in info['tts_bad']:
        a, b, c = tt
        # The TT (a,b,c) with a->b->c, a->c corresponds to the allowed 2-path (a,b,c)
        # if a->b and b->c are edges. Check:
        if A[a][b] and A[b][c]:
            path_abc = (a, b, c)
        else:
            # Try other orderings that form allowed paths
            path_abc = None
            for p in permutations([a, b, c]):
                if A[p[0]][p[1]] and A[p[1]][p[2]]:
                    path_abc = p
                    break

        if path_abc is None or path_abc not in idx2:
            results[tt] = {'in_allowed2': False}
            continue

        col_idx = idx2[path_abc]

        # Express this path in Ω₂ coordinates
        # path_abc as standard basis vector in A₂
        e = np.zeros(len(allowed2))
        e[col_idx] = 1.0

        # Check if it's in Ω₂: project onto omega2
        if dim_omega2 == 0:
            results[tt] = {'in_omega2': False}
            continue

        # omega2 columns span Ω₂. Check if e is in column space.
        # Actually, the path itself need not be in Ω₂; Ω₂ elements are
        # LINEAR COMBINATIONS of allowed paths. What we want is:
        # Does removing this specific TT (as a potential Ω₂ generator) change things?

        # Better approach: find the boundary of this path in A₁
        bd_path = bd2_full[:, col_idx]  # boundary of (a,b,c)

        # Check: is this path an element of Ω₂?
        # Path (a,b,c) ∈ Ω₂ iff all faces of (a,b,c) are in A₁ (i.e., allowed)
        # faces: (b,c), (a,c), (a,b)
        faces = [(b,c), (a,c), (a,b)] if path_abc == (a,b,c) else None
        if faces is None:
            p = path_abc
            faces = [(p[1],p[2]), (p[0],p[2]), (p[0],p[1])]

        allowed1_set = set(allowed1)
        all_faces_allowed = all(f in allowed1_set for f in faces)

        # For a TT (a,b,c) with a->b, b->c, a->c:
        # faces are (b,c), (a,c), (a,b) — all are directed edges, so all allowed
        # This means (a,b,c) ∈ Ω₂ (it's a "doubly-transitive" path)

        results[tt] = {
            'path': path_abc,
            'in_allowed2': True,
            'all_faces_allowed': all_faces_allowed,
            'in_omega2_direct': all_faces_allowed,
        }

        if not all_faces_allowed:
            continue

        # Now test: remove this column from Ω₂ and check rank drop
        # The path is in A₂ and also in Ω₂ (since all faces allowed).
        # Express it in Ω₂ coordinates: solve omega2 @ c = e
        # If omega2 is orthonormal (from SVD), c = omega2^T @ e
        c_coords = omega2.T @ e
        if np.linalg.norm(c_coords) < 1e-8:
            results[tt]['in_omega2'] = False
            continue

        results[tt]['in_omega2'] = True

        # Remove this direction from Ω₂: project out c_coords direction
        c_hat = c_coords / np.linalg.norm(c_coords)
        # New Ω₂ basis: omega2 minus the c_hat direction
        # omega2_new = omega2 @ (I - c_hat c_hat^T)
        proj = np.eye(dim_omega2) - np.outer(c_hat, c_hat)
        omega2_reduced = omega2 @ proj

        # Compute rank of ∂₂ on reduced space
        bd2_reduced = bd2_full @ omega2_reduced
        reduced_rank = rank_of(bd2_reduced)

        results[tt]['full_rank'] = full_rank
        results[tt]['reduced_rank'] = reduced_rank
        results[tt]['rank_dropped'] = (reduced_rank < full_rank)

    # Also test ALL TTs in the tournament for rank-criticality
    all_tts = find_transitive_triples(A, list(range(n)))
    n_rank_critical = 0
    n_total_tts = len(all_tts)
    rank_critical_tts = []

    for tt in all_tts:
        a, b, c = tt
        # Find allowed path
        path_abc = None
        for p in permutations([a, b, c]):
            if A[p[0]][p[1]] and A[p[1]][p[2]]:
                path_abc = p
                break
        if path_abc is None or path_abc not in idx2:
            continue

        col_idx = idx2[path_abc]
        e = np.zeros(len(allowed2))
        e[col_idx] = 1.0

        # Check faces allowed
        p = path_abc
        faces = [(p[1],p[2]), (p[0],p[2]), (p[0],p[1])]
        allowed1_set = set(allowed1)
        if not all(f in allowed1_set for f in faces):
            continue

        c_coords = omega2.T @ e
        if np.linalg.norm(c_coords) < 1e-8:
            continue

        c_hat = c_coords / np.linalg.norm(c_coords)
        proj = np.eye(dim_omega2) - np.outer(c_hat, c_hat)
        omega2_reduced = omega2 @ proj
        bd2_reduced = bd2_full @ omega2_reduced
        reduced_rank = rank_of(bd2_reduced)

        if reduced_rank < full_rank:
            n_rank_critical += 1
            rank_critical_tts.append(frozenset(tt))

    return results, full_rank, dim_omega2, n_rank_critical, n_total_tts, rank_critical_tts

def dual_test(A, n):
    """Test 4: Compute ker(∂₂^T) and related dual quantities."""
    allowed, omega2, omega1, bd2_full, bd2_omega = get_omega2_data(A, n)

    # ∂₂^T: A₁^* → A₂^*
    bd2_T = bd2_omega.T  # transpose of ∂₂|_{Ω₂}

    # ker(∂₂^T)
    if bd2_T.shape[0] > 0 and bd2_T.shape[1] > 0:
        U, S, Vt = np.linalg.svd(bd2_T, full_matrices=True)
        r = sum(s > 1e-8 for s in S)
        dim_ker_bd2T = bd2_T.shape[1] - r  # columns minus rank
    else:
        dim_ker_bd2T = 0

    dim_omega2_val = omega2.shape[1] if omega2.ndim == 2 else 0
    rank_bd2 = rank_of(bd2_omega)

    return {
        'dim_omega2': dim_omega2_val,
        'rank_bd2': rank_bd2,
        'dim_ker_bd2T': dim_ker_bd2T,
        'bd2_omega_shape': bd2_omega.shape,
    }

# ========== RUN ANALYSIS ==========

print("=" * 70)
print("GROK'S BAD-VERTEX TT BASIS CLAIM — EXHAUSTIVE TEST")
print("=" * 70)

for n in [5, 6]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    # Collect statistics
    stats_b1_0 = defaultdict(list)  # k_bad -> list of info dicts
    stats_b1_1 = []

    count = 0
    for A in all_tournaments(n):
        count += 1
        if count % 500 == 0:
            print(f"  ... processed {count} tournaments", file=sys.stderr)

        info = analyze_tournament(A, n)

        if info['b1'] == 0:
            stats_b1_0[info['k_bad']].append((A, info))
        else:
            stats_b1_1.append((A, info))

    print(f"\nTotal tournaments: {count}")
    print(f"β₁=0: {sum(len(v) for v in stats_b1_0.values())}")
    print(f"β₁=1: {len(stats_b1_1)}")

    # ======== TEST 1: β₁=0, by number of bad vertices ========
    print(f"\n--- TEST 1: β₁=0 tournaments, bad vertex structure ---")
    for k in sorted(stats_b1_0.keys()):
        entries = stats_b1_0[k]
        print(f"\n  k_bad = {k}: {len(entries)} tournaments")

        if k >= 2:
            n_trans = sum(1 for _, info in entries if info['is_trans_bad'])
            print(f"    Bad verts transitive: {n_trans}/{len(entries)}")

        if k >= 3:
            n_3cyc = sum(1 for _, info in entries if info['has_3cyc_bad'])
            tt_counts = [len(info['tts_bad']) for _, info in entries]
            print(f"    Bad verts have 3-cycle: {n_3cyc}/{len(entries)}")
            print(f"    TT count among bad verts: min={min(tt_counts)}, max={max(tt_counts)}, mean={np.mean(tt_counts):.2f}")

        # Detailed rank-drop analysis for k=3 (or k>=2 with TTs)
        if k >= 3 and len(entries) > 0:
            # Sample up to 50 for rank analysis (expensive)
            sample = entries[:50] if len(entries) > 50 else entries
            print(f"\n    Rank-drop analysis (sample of {len(sample)}):")

            n_bad_tt_rank_drops = 0
            n_bad_tt_tested = 0
            all_rank_critical_counts = []
            all_total_tt_counts = []
            bad_tt_is_unique_critical = 0

            for A_t, info in sample:
                rd_results, full_rank, dim_o2, n_rc, n_tt, rc_tts = rank_drop_test(A_t, n, info)

                for tt, rd in rd_results.items():
                    if rd.get('in_omega2', False):
                        n_bad_tt_tested += 1
                        if rd.get('rank_dropped', False):
                            n_bad_tt_rank_drops += 1

                all_rank_critical_counts.append(n_rc)
                all_total_tt_counts.append(n_tt)

                # Check if bad TTs are among the rank-critical ones
                bad_tt_sets = [frozenset(tt) for tt in info['tts_bad']]
                bad_in_critical = sum(1 for btt in bad_tt_sets if btt in rc_tts)
                if n_rc > 0 and bad_in_critical == n_rc:
                    bad_tt_is_unique_critical += 1

            print(f"    Bad TTs in Ω₂ that drop rank: {n_bad_tt_rank_drops}/{n_bad_tt_tested}")
            print(f"    Total rank-critical TTs: min={min(all_rank_critical_counts)}, max={max(all_rank_critical_counts)}, mean={np.mean(all_rank_critical_counts):.2f}")
            print(f"    Total TTs in tournament: min={min(all_total_tt_counts)}, max={max(all_total_tt_counts)}, mean={np.mean(all_total_tt_counts):.2f}")
            print(f"    Bad TT is the ONLY critical TT: {bad_tt_is_unique_critical}/{len(sample)}")

    # ======== TEST 2: β₁=1 tournaments ========
    print(f"\n--- TEST 2: β₁=1 tournaments ---")
    if stats_b1_1:
        bad_counts = [info['k_bad'] for _, info in stats_b1_1]
        print(f"  Count: {len(stats_b1_1)}")
        print(f"  Bad vertex counts: min={min(bad_counts)}, max={max(bad_counts)}, mean={np.mean(bad_counts):.2f}")

        bc_dist = defaultdict(int)
        for bc in bad_counts:
            bc_dist[bc] += 1
        print(f"  Distribution: {dict(sorted(bc_dist.items()))}")

        # Check sub-tournament structure of bad vertices
        n_trans_bad = 0
        n_sc_bad = 0
        tt_counts_b1 = []

        for A_t, info in stats_b1_1:
            if info['k_bad'] >= 2:
                if info['is_trans_bad']:
                    n_trans_bad += 1
                if info['has_3cyc_bad']:
                    n_sc_bad += 1
            tt_counts_b1.append(len(info['tts_bad']))

        n_with_multiple_bad = sum(1 for _, info in stats_b1_1 if info['k_bad'] >= 2)
        if n_with_multiple_bad > 0:
            print(f"  Among those with ≥2 bad verts ({n_with_multiple_bad}):")
            print(f"    Bad sub-tournament transitive: {n_trans_bad}")
            print(f"    Bad sub-tournament has 3-cycle (SC): {n_sc_bad}")

        # Rank-drop for β₁=1
        sample = stats_b1_1[:30] if len(stats_b1_1) > 30 else stats_b1_1
        print(f"\n  Rank-drop analysis for β₁=1 (sample of {len(sample)}):")
        for idx, (A_t, info) in enumerate(sample[:5]):
            if info['k_bad'] >= 3 and len(info['tts_bad']) > 0:
                rd_results, full_rank, dim_o2, n_rc, n_tt, rc_tts = rank_drop_test(A_t, n, info)
                print(f"    T#{idx}: k_bad={info['k_bad']}, tts_bad={len(info['tts_bad'])}, rank(∂₂)={full_rank}, dim(Ω₂)={dim_o2}, rank-critical TTs={n_rc}/{n_tt}")
                for tt, rd in rd_results.items():
                    print(f"      Bad TT {tt}: in_Ω₂={rd.get('in_omega2','?')}, rank_drop={rd.get('rank_dropped','?')}")
    else:
        print(f"  No β₁=1 tournaments at n={n}")

    # ======== TEST 3: Pigeonhole — uniqueness of bad TT ========
    print(f"\n--- TEST 3: Pigeonhole — bad-TT uniqueness (k_bad=3, β₁=0) ---")
    if 3 in stats_b1_0 and len(stats_b1_0[3]) > 0:
        sample = stats_b1_0[3][:50] if len(stats_b1_0[3]) > 50 else stats_b1_0[3]

        bad_tt_special = 0
        bad_tt_one_of_many = 0
        bad_tt_not_critical = 0

        for A_t, info in sample:
            rd_results, full_rank, dim_o2, n_rc, n_tt, rc_tts = rank_drop_test(A_t, n, info)

            bad_tt_sets = [frozenset(tt) for tt in info['tts_bad']]
            bad_critical = any(btt in rc_tts for btt in bad_tt_sets)

            if bad_critical and n_rc == 1:
                bad_tt_special += 1
            elif bad_critical and n_rc > 1:
                bad_tt_one_of_many += 1
            else:
                bad_tt_not_critical += 1

        print(f"  Bad TT is uniquely rank-critical: {bad_tt_special}/{len(sample)}")
        print(f"  Bad TT is critical but not unique: {bad_tt_one_of_many}/{len(sample)}")
        print(f"  Bad TT is NOT rank-critical: {bad_tt_not_critical}/{len(sample)}")
    else:
        print(f"  No k_bad=3 tournaments with β₁=0")

    # ======== TEST 4: Dual test ========
    print(f"\n--- TEST 4: Dual test — ker(∂₂^T) ---")
    # Sample from β₁=0
    for k in sorted(stats_b1_0.keys()):
        entries = stats_b1_0[k]
        sample = entries[:10] if len(entries) > 10 else entries
        dim_kers = []
        for A_t, info in sample:
            dt = dual_test(A_t, n)
            dim_kers.append(dt['dim_ker_bd2T'])
        if dim_kers:
            print(f"  β₁=0, k_bad={k}: dim(ker(∂₂^T)) = {set(dim_kers)} (from {len(sample)} samples)")

    # β₁=1
    if stats_b1_1:
        sample = stats_b1_1[:10]
        dim_kers = []
        for A_t, info in sample:
            dt = dual_test(A_t, n)
            dim_kers.append(dt['dim_ker_bd2T'])
        print(f"  β₁=1: dim(ker(∂₂^T)) = {set(dim_kers)} (from {len(sample)} samples)")

    # ======== TEST 5: Hypothetical 4th bad vertex (n=6 only) ========
    if n == 6:
        print(f"\n--- TEST 5: Hypothetical 4th bad vertex analysis ---")
        # Among β₁=0 with k_bad=3, find cases where a 4th vertex
        # "almost" looks bad (highest β₁(T\v) sensitivity)
        if 3 in stats_b1_0:
            sample = stats_b1_0[3][:20]
            print(f"  Checking {len(sample)} tournaments with k_bad=3...")

            for idx, (A_t, info) in enumerate(sample[:5]):
                bad = info['bad']
                non_bad = [v for v in range(n) if v not in bad]

                # For each non-bad vertex, check deletion β₁
                # (should be 0 since vertex is "good")
                # Check: if we form 4-vertex set bad ∪ {v}, is it transitive?
                print(f"\n    T#{idx}: bad={bad}")
                for v in non_bad:
                    ext_set = bad + [v]
                    trans = is_transitive_subtournament(A_t, ext_set)
                    tts_ext = find_transitive_triples(A_t, ext_set)

                    # How many TTs in the extended set?
                    # A transitive 4-tournament has C(4,3)=4 TTs
                    # A non-transitive one has fewer
                    print(f"      + v={v}: bad∪{{v}} transitive={trans}, #TTs={len(tts_ext)}")

                    if trans:
                        # This is the interesting case for Grok's argument
                        # Among the 4 TTs, which are "isolated"?
                        # (meaning: the TT doesn't share edges with non-bad-TTs)
                        print(f"        TTs in extended set: {tts_ext}")

                        # Check which TTs among the 4 bad verts would be rank-critical
                        # This is hypothetical: the 4th vertex isn't actually bad
                        # But we can ask: would a 4-vertex transitive bad set
                        # create more rank-critical TTs than can be accommodated?

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
