#!/usr/bin/env python3
"""
CORRECTED cascade test: flip ONLY the transitive edge source->sink
among the 3 bad vertices (b1-deletion bad) to create a 3-cycle.

Definition: vertex v is "bad" if b1(T\\v) > 0.
For tournaments with beta1(T)=0 and exactly 3 bad vertices:
  - The 3 bad vertices form a transitive triple (confirmed)
  - Flip ONLY source->sink to sink->source (creating 3-cycle among bad verts)
  - Track TT changes, boundary spaces, H1 generator

Bug in previous version: reversed ALL 3 edges instead of just one.
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict
import sys
import os
import random

sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')

# ============================================================
# CORE FUNCTIONS (inlined to avoid validation print)
# ============================================================

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

# ============================================================
# DETAILED HOMOLOGY
# ============================================================

def get_detailed_homology(A, n, p_target):
    """Get detailed homology info at dimension p_target."""
    max_dim = p_target + 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    p = p_target
    dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0

    bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
    bd_p_omega = bd_p @ omega[p]

    U, S, Vt = np.linalg.svd(bd_p_omega, full_matrices=True)
    rank_p = sum(s > 1e-8 for s in S)
    ker_basis_omega = Vt[rank_p:].T  # in Omega_p coords
    ker_basis_Ap = omega[p] @ ker_basis_omega  # in A_p coords

    dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
    if dim_omega_p1 > 0:
        bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
        bd_p1_omega = bd_p1 @ omega[p+1]
        U2, S2, Vt2 = np.linalg.svd(bd_p1_omega, full_matrices=True)
        im_rank = sum(s > 1e-8 for s in S2)
        im_basis_Ap = U2[:, :im_rank]
    else:
        im_basis_Ap = np.zeros((len(allowed[p]), 0))
        im_rank = 0

    beta_p = ker_basis_omega.shape[1] - im_rank

    # H_p generators
    if beta_p > 0 and im_rank > 0:
        Q_im, _ = np.linalg.qr(im_basis_Ap)
        Q_im = Q_im[:, :im_rank]
        ker_proj = ker_basis_Ap - Q_im @ (Q_im.T @ ker_basis_Ap)
        norms = np.linalg.norm(ker_proj, axis=0)
        gen_idx = np.where(norms > 1e-8)[0]
        generators_Ap = ker_proj[:, gen_idx]
    elif beta_p > 0:
        generators_Ap = ker_basis_Ap
    else:
        generators_Ap = np.zeros((len(allowed[p]), 0))

    return {
        'allowed_p': allowed[p],
        'allowed_pm1': allowed[p-1] if p > 0 else [],
        'omega_basis': omega[p],
        'ker_basis': ker_basis_Ap,
        'im_basis': im_basis_Ap,
        'generators': generators_Ap,
        'beta': beta_p,
        'dim_omega': dim_omega_p,
        'dim_ker': ker_basis_omega.shape[1],
        'dim_im': im_rank,
    }

# ============================================================
# TOURNAMENT UTILITIES
# ============================================================

def generate_all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(edges):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield [row[:] for row in A]

def subtournament(A, n, exclude_v):
    """Return (A', n', vertex_map) for T\\v."""
    verts = [i for i in range(n) if i != exclude_v]
    n2 = len(verts)
    A2 = [[0]*n2 for _ in range(n2)]
    for i in range(n2):
        for j in range(n2):
            A2[i][j] = A[verts[i]][verts[j]]
    return A2, n2, verts

def find_bad_vertices(A, n):
    """Bad vertex v: b1(T\\v) > 0."""
    bad = []
    for v in range(n):
        A2, n2, _ = subtournament(A, n, v)
        b = path_betti_numbers(A2, n2, max_dim=1)
        if b[1] > 0:
            bad.append(v)
    return bad

def compute_TT_set(A, n):
    """Set of transitive triples as frozensets."""
    tt_set = set()
    for a, b, c in combinations(range(n), 3):
        cyc1 = A[a][b] and A[b][c] and A[c][a]
        cyc2 = A[a][c] and A[c][b] and A[b][a]
        if not cyc1 and not cyc2:
            tt_set.add(frozenset([a, b, c]))
    return tt_set

def format_chain(vec, paths, threshold=1e-8):
    """Format a chain vector as human-readable expression."""
    terms = []
    for i, coeff in enumerate(vec):
        if abs(coeff) > threshold:
            sign = '+' if coeff > 0 else '-'
            c = abs(coeff)
            path_str = '->'.join(str(v) for v in paths[i])
            if abs(c - round(c)) < 1e-6:
                c = int(round(c))
                if c == 1:
                    terms.append(f"{sign}({path_str})")
                else:
                    terms.append(f"{sign}{c}({path_str})")
            else:
                terms.append(f"{sign}{c:.4f}({path_str})")
    if not terms:
        return "0"
    result = ' '.join(terms)
    if result.startswith('+'):
        result = result[1:]
    return result

# ============================================================
# MAIN ANALYSIS
# ============================================================

def analyze_one(A, n, verbose=True, label=""):
    """Analyze a single tournament: find bad verts, flip, track everything."""

    bad_verts = find_bad_vertices(A, n)
    if len(bad_verts) != 3:
        return None

    # Check β₁(T) = 0
    betti_orig = path_betti_numbers(A, n, max_dim=2)
    if betti_orig[1] != 0:
        return None

    # The 3 bad vertices should form a transitive triple
    a, b, c = bad_verts
    cyc1 = A[a][b] and A[b][c] and A[c][a]
    cyc2 = A[a][c] and A[c][b] and A[b][a]
    if cyc1 or cyc2:
        print(f"  WARNING: bad vertices {bad_verts} form a 3-CYCLE, not transitive triple!")
        return None

    # Find source, mid, sink
    scores = {v: sum(A[v][w] for w in bad_verts if w != v) for v in bad_verts}
    source = [v for v in bad_verts if scores[v] == 2][0]
    sink = [v for v in bad_verts if scores[v] == 0][0]
    mid = [v for v in bad_verts if scores[v] == 1][0]

    # Verify: source->mid, mid->sink, source->sink
    assert A[source][mid] == 1, f"Expected {source}->{mid}"
    assert A[mid][sink] == 1, f"Expected {mid}->{sink}"
    assert A[source][sink] == 1, f"Expected {source}->{sink}"

    # TTs before flip
    tt_before = compute_TT_set(A, n)

    # CORRECT FLIP: ONLY source->sink becomes sink->source
    A_new = [row[:] for row in A]
    A_new[source][sink] = 0
    A_new[sink][source] = 1

    # Verify only 2 entries changed
    changes = sum(1 for i in range(n) for j in range(n) if A[i][j] != A_new[i][j])
    assert changes == 2, f"Expected 2 entry changes, got {changes}"

    # Verify 3-cycle created: source->mid->sink->source
    assert A_new[source][mid] == 1, "source->mid preserved"
    assert A_new[mid][sink] == 1, "mid->sink preserved"
    assert A_new[sink][source] == 1, "sink->source (new)"

    # TTs after flip
    tt_after = compute_TT_set(A_new, n)

    destroyed = tt_before - tt_after
    created = tt_after - tt_before
    kept = tt_before & tt_after

    # Betti after flip
    betti_new = path_betti_numbers(A_new, n, max_dim=2)

    # Detailed H1 for both
    h1_orig = get_detailed_homology(A, n, 1)
    h1_new = get_detailed_homology(A_new, n, 1)

    result = {
        'source': source, 'mid': mid, 'sink': sink,
        'bad_verts': bad_verts,
        'beta1_orig': betti_orig[1],
        'beta1_new': betti_new[1],
        'betti_orig': betti_orig,
        'betti_new': betti_new,
        'tt_before': len(tt_before),
        'tt_after': len(tt_after),
        'destroyed': destroyed,
        'created': created,
        'kept': kept,
        'n_destroyed': len(destroyed),
        'n_created': len(created),
        'n_kept': len(kept),
        'orig_dim_omega1': h1_orig['dim_omega'],
        'orig_dim_ker1': h1_orig['dim_ker'],
        'orig_dim_im1': h1_orig['dim_im'],
        'new_dim_omega1': h1_new['dim_omega'],
        'new_dim_ker1': h1_new['dim_ker'],
        'new_dim_im1': h1_new['dim_im'],
    }

    if h1_new['generators'].shape[1] > 0:
        gen = h1_new['generators'][:, 0]
        max_abs = np.max(np.abs(gen))
        if max_abs > 1e-10:
            gen = gen / max_abs
        result['h1_generator_str'] = format_chain(gen, h1_new['allowed_p'])
        result['h1_generator_vec'] = gen
        result['h1_allowed_1paths'] = h1_new['allowed_p']

    # Also track what specific boundary direction is lost
    # Compare im(∂₂) old vs new
    if h1_orig['dim_im'] > h1_new['dim_im']:
        # There's a rank drop. Find which direction was lost.
        # Project old image basis vectors onto new image space
        im_old = h1_orig['im_basis']
        im_new = h1_new['im_basis']
        if im_new.shape[1] > 0:
            Q_new, _ = np.linalg.qr(im_new)
            Q_new = Q_new[:, :h1_new['dim_im']]
            # Find old image vectors NOT in new image
            residuals = im_old - Q_new @ (Q_new.T @ im_old)
            norms = np.linalg.norm(residuals, axis=0)
            lost_idx = np.where(norms > 1e-6)[0]
            if len(lost_idx) > 0:
                # The lost boundary direction
                lost_dir = residuals[:, lost_idx[0]]
                lost_dir = lost_dir / np.max(np.abs(lost_dir))
                # Need to use the ORIGINAL allowed 1-paths for formatting
                result['lost_boundary_str'] = format_chain(lost_dir, h1_orig['allowed_p'])
        else:
            # All boundary is lost (unlikely)
            result['lost_boundary_str'] = "(all boundary lost)"

    return result

def print_result(r, label=""):
    print(f"\n{'='*60}")
    print(f"Tournament {label}")
    print(f"  Bad vertices: {r['bad_verts']} (source={r['source']}, mid={r['mid']}, sink={r['sink']})")
    print(f"  Flip: {r['source']}->{r['sink']} becomes {r['sink']}->{r['source']}")
    print(f"  Creates 3-cycle: {r['source']}->{r['mid']}->{r['sink']}->{r['source']}")
    print(f"  beta1: {r['beta1_orig']} -> {r['beta1_new']}")
    print(f"  Full Betti before: {r['betti_orig']}")
    print(f"  Full Betti after:  {r['betti_new']}")
    print(f"  TT count: {r['tt_before']} -> {r['tt_after']}")
    print(f"  Destroyed: {r['n_destroyed']}, Created: {r['n_created']}, Kept: {r['n_kept']}")
    if r['n_destroyed'] <= 10:
        for tt in sorted(r['destroyed'], key=lambda s: tuple(sorted(s))):
            verts = sorted(tt)
            print(f"    D: {{{','.join(map(str,verts))}}}")
    if r['n_created'] <= 10:
        for tt in sorted(r['created'], key=lambda s: tuple(sorted(s))):
            verts = sorted(tt)
            print(f"    C: {{{','.join(map(str,verts))}}}")
    print(f"  Omega1 dim: {r['orig_dim_omega1']} -> {r['new_dim_omega1']}")
    print(f"  ker(d1) dim: {r['orig_dim_ker1']} -> {r['new_dim_ker1']}")
    print(f"  im(d2) dim:  {r['orig_dim_im1']} -> {r['new_dim_im1']}")
    if 'h1_generator_str' in r:
        print(f"  H1 generator: {r['h1_generator_str']}")
    if 'lost_boundary_str' in r:
        print(f"  Lost boundary direction: {r['lost_boundary_str']}")

def main():
    print("="*70)
    print("CORRECTED CASCADE TEST")
    print("Bad vertex = v where b1(T\\v) > 0")
    print("Flip ONLY source->sink (transitive edge) to create 3-cycle")
    print("="*70)

    # ==================== n=5 ====================
    print("\n" + "#"*70)
    print("# n=5: ALL tournaments with beta1=0 and 3 bad vertices")
    print("#"*70)

    n = 5
    count = 0
    results_n5 = []

    beta1_new_counts = defaultdict(int)
    destroyed_counts = defaultdict(int)
    created_counts = defaultdict(int)

    for A in generate_all_tournaments(n):
        betti = path_betti_numbers(A, n, max_dim=1)
        if betti[1] != 0:
            continue
        bad = find_bad_vertices(A, n)
        if len(bad) != 3:
            continue

        count += 1
        r = analyze_one(A, n, verbose=(count <= 5), label=f"n5-{count}")
        if r is None:
            continue
        results_n5.append(r)

        beta1_new_counts[r['beta1_new']] += 1
        destroyed_counts[r['n_destroyed']] += 1
        created_counts[r['n_created']] += 1

        if count <= 5:
            print_result(r, label=f"n5-{count}")

        if count % 20 == 0:
            print(f"  ... processed {count} tournaments", flush=True)

    print(f"\n{'='*70}")
    print(f"n=5 SUMMARY: {count} tournaments analyzed")
    print(f"{'='*70}")
    print(f"beta1 after flip: {dict(beta1_new_counts)}")
    print(f"# destroyed TTs: {dict(destroyed_counts)}")
    print(f"# created TTs: {dict(created_counts)}")

    # Dimension changes
    omega_changes = defaultdict(int)
    ker_changes = defaultdict(int)
    im_changes = defaultdict(int)
    for r in results_n5:
        omega_changes[(r['orig_dim_omega1'], r['new_dim_omega1'])] += 1
        ker_changes[(r['orig_dim_ker1'], r['new_dim_ker1'])] += 1
        im_changes[(r['orig_dim_im1'], r['new_dim_im1'])] += 1
    print(f"Omega1 dim (before->after): {dict(omega_changes)}")
    print(f"ker(d1) dim (before->after): {dict(ker_changes)}")
    print(f"im(d2) dim (before->after): {dict(im_changes)}")

    # Destroyed TT overlap with bad triple
    print(f"\nDESTROYED TT PATTERNS (overlap with bad triple):")
    destroy_overlap = defaultdict(int)
    create_overlap = defaultdict(int)
    for r in results_n5:
        bad_set = set(r['bad_verts'])
        for tt in r['destroyed']:
            overlap = len(tt & bad_set)
            destroy_overlap[overlap] += 1
        for tt in r['created']:
            overlap = len(tt & bad_set)
            create_overlap[overlap] += 1
    for k in sorted(destroy_overlap):
        print(f"  Destroyed, {k} verts overlap: {destroy_overlap[k]}")
    for k in sorted(create_overlap):
        print(f"  Created, {k} verts overlap: {create_overlap[k]}")

    # Which destroyed TTs contain source-sink?
    print(f"\nDESTROYED TTs containing source AND sink:")
    src_sink_count = 0
    src_only = 0
    sink_only = 0
    neither = 0
    for r in results_n5:
        for tt in r['destroyed']:
            has_src = r['source'] in tt
            has_snk = r['sink'] in tt
            if has_src and has_snk:
                src_sink_count += 1
            elif has_src:
                src_only += 1
            elif has_snk:
                sink_only += 1
            else:
                neither += 1
    print(f"  Both source+sink: {src_sink_count}")
    print(f"  Source only: {src_only}")
    print(f"  Sink only: {sink_only}")
    print(f"  Neither: {neither}")

    # Rank drop analysis
    print(f"\nRANK DROP ANALYSIS:")
    rank_drops = defaultdict(int)
    for r in results_n5:
        drop = r['orig_dim_im1'] - r['new_dim_im1']
        rank_drops[drop] += 1
    print(f"  im(d2) rank drop: {dict(rank_drops)}")

    omega_inc = defaultdict(int)
    for r in results_n5:
        inc = r['new_dim_omega1'] - r['orig_dim_omega1']
        omega_inc[inc] += 1
    print(f"  Omega1 dim increase: {dict(omega_inc)}")

    ker_inc = defaultdict(int)
    for r in results_n5:
        inc = r['new_dim_ker1'] - r['orig_dim_ker1']
        ker_inc[inc] += 1
    print(f"  ker(d1) dim increase: {dict(ker_inc)}")

    # H1 generators
    print(f"\nH1 GENERATORS (first 10 with beta1'=1):")
    gen_count = 0
    for r in results_n5:
        if 'h1_generator_str' in r and gen_count < 10:
            gen_count += 1
            print(f"  T{gen_count} (src={r['source']},mid={r['mid']},snk={r['sink']}): {r['h1_generator_str']}")
            if 'lost_boundary_str' in r:
                print(f"    Lost boundary: {r['lost_boundary_str']}")

    # ==================== n=6 ====================
    print("\n\n" + "#"*70)
    print("# n=6: 50 tournaments with beta1=0 and 3 bad vertices")
    print("#"*70)

    n = 6
    count6 = 0
    results_n6 = []

    beta1_new_counts6 = defaultdict(int)
    destroyed_counts6 = defaultdict(int)
    created_counts6 = defaultdict(int)

    for A in generate_all_tournaments(n):
        if count6 >= 50:
            break

        betti = path_betti_numbers(A, n, max_dim=1)
        if betti[1] != 0:
            continue
        bad = find_bad_vertices(A, n)
        if len(bad) != 3:
            continue

        count6 += 1
        r = analyze_one(A, n, verbose=(count6 <= 3), label=f"n6-{count6}")
        if r is None:
            continue
        results_n6.append(r)

        beta1_new_counts6[r['beta1_new']] += 1
        destroyed_counts6[r['n_destroyed']] += 1
        created_counts6[r['n_created']] += 1

        if count6 <= 3:
            print_result(r, label=f"n6-{count6}")

        if count6 % 10 == 0:
            print(f"  ... processed {count6}/50 tournaments", flush=True)

    print(f"\n{'='*70}")
    print(f"n=6 SUMMARY: {count6} tournaments analyzed")
    print(f"{'='*70}")
    print(f"beta1 after flip: {dict(beta1_new_counts6)}")
    print(f"# destroyed TTs: {dict(destroyed_counts6)}")
    print(f"# created TTs: {dict(created_counts6)}")

    omega_changes6 = defaultdict(int)
    ker_changes6 = defaultdict(int)
    im_changes6 = defaultdict(int)
    for r in results_n6:
        omega_changes6[(r['orig_dim_omega1'], r['new_dim_omega1'])] += 1
        ker_changes6[(r['orig_dim_ker1'], r['new_dim_ker1'])] += 1
        im_changes6[(r['orig_dim_im1'], r['new_dim_im1'])] += 1
    print(f"Omega1 dim (before->after): {dict(omega_changes6)}")
    print(f"ker(d1) dim (before->after): {dict(ker_changes6)}")
    print(f"im(d2) dim (before->after): {dict(im_changes6)}")

    # Destroyed TT patterns
    print(f"\nDESTROYED TT PATTERNS:")
    destroy_overlap6 = defaultdict(int)
    create_overlap6 = defaultdict(int)
    for r in results_n6:
        bad_set = set(r['bad_verts'])
        for tt in r['destroyed']:
            overlap = len(tt & bad_set)
            destroy_overlap6[overlap] += 1
        for tt in r['created']:
            overlap = len(tt & bad_set)
            create_overlap6[overlap] += 1
    for k in sorted(destroy_overlap6):
        print(f"  Destroyed, {k} verts overlap: {destroy_overlap6[k]}")
    for k in sorted(create_overlap6):
        print(f"  Created, {k} verts overlap: {create_overlap6[k]}")

    # Source-sink in destroyed
    print(f"\nDESTROYED TTs containing source AND sink:")
    src_sink6 = 0
    src_only6 = 0
    sink_only6 = 0
    neither6 = 0
    for r in results_n6:
        for tt in r['destroyed']:
            has_src = r['source'] in tt
            has_snk = r['sink'] in tt
            if has_src and has_snk:
                src_sink6 += 1
            elif has_src:
                src_only6 += 1
            elif has_snk:
                sink_only6 += 1
            else:
                neither6 += 1
    print(f"  Both source+sink: {src_sink6}")
    print(f"  Source only: {src_only6}")
    print(f"  Sink only: {sink_only6}")
    print(f"  Neither: {neither6}")

    rank_drops6 = defaultdict(int)
    for r in results_n6:
        drop = r['orig_dim_im1'] - r['new_dim_im1']
        rank_drops6[drop] += 1
    print(f"\nRANK DROP: im(d2) drop = {dict(rank_drops6)}")

    omega_inc6 = defaultdict(int)
    ker_inc6 = defaultdict(int)
    for r in results_n6:
        omega_inc6[r['new_dim_omega1'] - r['orig_dim_omega1']] += 1
        ker_inc6[r['new_dim_ker1'] - r['orig_dim_ker1']] += 1
    print(f"Omega1 increase: {dict(omega_inc6)}")
    print(f"ker(d1) increase: {dict(ker_inc6)}")

    # H1 generators for n=6
    print(f"\nH1 GENERATORS (first 5):")
    gen_count6 = 0
    for r in results_n6:
        if 'h1_generator_str' in r and gen_count6 < 5:
            gen_count6 += 1
            print(f"  T{gen_count6}: {r['h1_generator_str']}")
            if 'lost_boundary_str' in r:
                print(f"    Lost boundary: {r['lost_boundary_str']}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
