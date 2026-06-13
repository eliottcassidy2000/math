#!/usr/bin/env python3
"""
beta2_dim_drop_n6.py — Investigate Omega_2 dimension drop when flipping
a transitive-triple arc to create a 3-cycle.

Hypothesis (verified at n=5):
  Flipping a TT arc a->c to c->a in a beta_1=0 tournament causes:
    - dim(Omega_2) drops by exactly 1
    - dim(Z_1) stays the same
    - beta_1 goes from 0 to 1
"""
import sys
import numpy as np
from itertools import combinations
import time, random

def flush_print(*args, **kwargs):
    print(*args, **kwargs, flush=True)

# ===== Core path homology functions =====

def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
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
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

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
    if dim_Ap == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0: return np.zeros((dim_Ap, 0))
    return null_space

def compute_betti_and_dims(A, n, max_dim=2):
    allowed = {}
    for p in range(-1, max_dim + 2):
        allowed[p] = [] if p < 0 else enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    result = {}
    for p in range(max_dim + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_omega_p > 0:
            bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
            bd_p_omega = bd_p @ omega[p]
            if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
                S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
                rank_p = sum(s > 1e-8 for s in S_p)
            else: rank_p = 0
            ker_dim = dim_omega_p - rank_p
        else: ker_dim = 0; rank_p = 0

        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else: im_dim = 0

        result[p] = {
            'dim_A': len(allowed[p]),
            'dim_Omega': dim_omega_p,
            'dim_Z': ker_dim,
            'dim_B': im_dim,
            'betti': max(0, ker_dim - im_dim)
        }
    return result

def adj_matrix(n, edges_mask):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if edges_mask & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def find_one_tt(A, n):
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    return (a, b, c)
    return None

def flip_arc(A, a, c):
    Af = [row[:] for row in A]
    Af[a][c] = 0
    Af[c][a] = 1
    return Af


# ===== PART 0: Check Omega_2 vs A_2 relationship =====
def check_omega2_vs_A2():
    flush_print("=" * 70)
    flush_print("PART 0: Omega_2 vs A_2 for tournaments")
    flush_print("=" * 70)

    for n in [4, 5, 6]:
        num_edges = n*(n-1)//2
        total = 2**num_edges
        eq_count = 0
        neq_count = 0
        neq_examples = []
        for mask in range(total):
            A = adj_matrix(n, mask)
            info = compute_betti_and_dims(A, n, max_dim=2)
            if info[2]['dim_Omega'] == info[2]['dim_A']:
                eq_count += 1
            else:
                neq_count += 1
                if len(neq_examples) < 3:
                    neq_examples.append((mask, info[2]['dim_Omega'], info[2]['dim_A']))
        flush_print(f"  n={n}: Omega_2=A_2 in {eq_count}/{total}, Omega_2!=A_2 in {neq_count}/{total}")
        if neq_examples:
            for m, o2, a2 in neq_examples:
                flush_print(f"    Example: mask={m}, dim_Omega_2={o2}, |A_2|={a2}, gap={a2-o2}")
    flush_print()


# ===== PART 1: n=5 verification =====
def run_n5():
    n = 5
    num_edges = n*(n-1)//2
    total = 2**num_edges

    flush_print("=" * 70)
    flush_print(f"PART 1: n={n} EXHAUSTIVE — TT flip analysis ({total} tournaments)")
    flush_print("=" * 70)

    t0 = time.time()

    d_omega2_dist = {}
    d_z1_dist = {}
    d_b1_dist = {}
    beta1_after_dist = {}
    candidates = 0
    detailed = []

    for mask in range(total):
        A = adj_matrix(n, mask)
        info_orig = compute_betti_and_dims(A, n, max_dim=2)
        if info_orig[1]['betti'] != 0:
            continue

        tt = find_one_tt(A, n)
        if tt is None:
            continue

        a, b, c = tt
        Af = flip_arc(A, a, c)
        info_flip = compute_betti_and_dims(Af, n, max_dim=2)

        d_omega2 = info_flip[2]['dim_Omega'] - info_orig[2]['dim_Omega']
        d_z1 = info_flip[1]['dim_Z'] - info_orig[1]['dim_Z']
        d_b1 = info_flip[1]['dim_B'] - info_orig[1]['dim_B']
        new_beta1 = info_flip[1]['betti']

        d_omega2_dist[d_omega2] = d_omega2_dist.get(d_omega2, 0) + 1
        d_z1_dist[d_z1] = d_z1_dist.get(d_z1, 0) + 1
        d_b1_dist[d_b1] = d_b1_dist.get(d_b1, 0) + 1
        beta1_after_dist[new_beta1] = beta1_after_dist.get(new_beta1, 0) + 1
        candidates += 1

        if len(detailed) < 10:
            detailed.append({
                'mask': mask, 'tt': tt,
                'orig_Omega2': info_orig[2]['dim_Omega'], 'flip_Omega2': info_flip[2]['dim_Omega'],
                'orig_A2': info_orig[2]['dim_A'], 'flip_A2': info_flip[2]['dim_A'],
                'orig_Z1': info_orig[1]['dim_Z'], 'flip_Z1': info_flip[1]['dim_Z'],
                'orig_B1': info_orig[1]['dim_B'], 'flip_B1': info_flip[1]['dim_B'],
                'd_omega2': d_omega2, 'd_z1': d_z1, 'd_b1': d_b1,
                'beta1_after': new_beta1,
            })

    elapsed = time.time() - t0
    flush_print(f"  {candidates} TT flips from beta_1=0 tournaments in {elapsed:.1f}s")
    flush_print(f"  dim(Omega_2) change: {dict(sorted(d_omega2_dist.items()))}")
    flush_print(f"  dim(Z_1) change: {dict(sorted(d_z1_dist.items()))}")
    flush_print(f"  dim(B_1) change: {dict(sorted(d_b1_dist.items()))}")
    flush_print(f"  beta_1 after: {dict(sorted(beta1_after_dist.items()))}")

    flush_print(f"\n  Detailed (first 10):")
    for d in detailed:
        flush_print(f"    mask={d['mask']}, TT={d['tt']}: "
              f"Omega2: {d['orig_Omega2']}->{d['flip_Omega2']} (d={d['d_omega2']}), "
              f"|A2|: {d['orig_A2']}->{d['flip_A2']} (d={d['flip_A2']-d['orig_A2']}), "
              f"Z1: {d['orig_Z1']}->{d['flip_Z1']} (d={d['d_z1']}), "
              f"B1: {d['orig_B1']}->{d['flip_B1']} (d={d['d_b1']}), "
              f"beta1: 0->{d['beta1_after']}")
    flush_print()


# ===== PART 2: n=6 exhaustive =====
def run_n6():
    n = 6
    num_edges = n*(n-1)//2
    total = 2**num_edges

    flush_print("=" * 70)
    flush_print(f"PART 2: n={n} EXHAUSTIVE — TT flip analysis ({total} tournaments)")
    flush_print("=" * 70)

    t0 = time.time()

    # First: find all beta_1=0 tournaments (fast pass with max_dim=1)
    b1_zero_masks = []
    for mask in range(total):
        A = adj_matrix(n, mask)
        info = compute_betti_and_dims(A, n, max_dim=1)
        if info[1]['betti'] == 0:
            b1_zero_masks.append(mask)

    elapsed = time.time() - t0
    flush_print(f"  Found {len(b1_zero_masks)}/{total} with beta_1=0 in {elapsed:.1f}s")

    # Second: for each, flip ONE TT and compute full homology
    d_omega2_dist = {}
    d_z1_dist = {}
    d_b1_dist = {}
    beta1_after_dist = {}
    candidates = 0
    detailed = []

    t1 = time.time()
    for idx, mask in enumerate(b1_zero_masks):
        A = adj_matrix(n, mask)
        tt = find_one_tt(A, n)
        if tt is None:
            continue

        a, b, c = tt
        Af = flip_arc(A, a, c)

        info_orig = compute_betti_and_dims(A, n, max_dim=2)
        info_flip = compute_betti_and_dims(Af, n, max_dim=2)

        d_omega2 = info_flip[2]['dim_Omega'] - info_orig[2]['dim_Omega']
        d_z1 = info_flip[1]['dim_Z'] - info_orig[1]['dim_Z']
        d_b1 = info_flip[1]['dim_B'] - info_orig[1]['dim_B']
        new_beta1 = info_flip[1]['betti']

        d_omega2_dist[d_omega2] = d_omega2_dist.get(d_omega2, 0) + 1
        d_z1_dist[d_z1] = d_z1_dist.get(d_z1, 0) + 1
        d_b1_dist[d_b1] = d_b1_dist.get(d_b1, 0) + 1
        beta1_after_dist[new_beta1] = beta1_after_dist.get(new_beta1, 0) + 1
        candidates += 1

        if len(detailed) < 15:
            detailed.append({
                'mask': mask, 'tt': tt,
                'orig_Omega2': info_orig[2]['dim_Omega'], 'flip_Omega2': info_flip[2]['dim_Omega'],
                'orig_A2': info_orig[2]['dim_A'], 'flip_A2': info_flip[2]['dim_A'],
                'orig_Z1': info_orig[1]['dim_Z'], 'flip_Z1': info_flip[1]['dim_Z'],
                'orig_B1': info_orig[1]['dim_B'], 'flip_B1': info_flip[1]['dim_B'],
                'd_omega2': d_omega2, 'd_z1': d_z1, 'd_b1': d_b1,
                'beta1_after': new_beta1,
            })

        if idx % 4000 == 0:
            elapsed = time.time() - t1
            flush_print(f"  Processing: {idx}/{len(b1_zero_masks)}, {elapsed:.1f}s")

    elapsed = time.time() - t0
    flush_print(f"\n  Total: {candidates} TT flips in {elapsed:.1f}s")

    flush_print(f"\n--- RESULTS (n={n}) ---")
    flush_print(f"  dim(Omega_2) change: {dict(sorted(d_omega2_dist.items()))}")
    flush_print(f"  dim(Z_1) change: {dict(sorted(d_z1_dist.items()))}")
    flush_print(f"  dim(B_1) change: {dict(sorted(d_b1_dist.items()))}")
    flush_print(f"  beta_1 after: {dict(sorted(beta1_after_dist.items()))}")

    # Hypothesis checks
    if len(d_omega2_dist) == 1 and -1 in d_omega2_dist:
        flush_print(f"\n  *** CONFIRMED: dim(Omega_2) ALWAYS drops by exactly 1 ***")
    else:
        flush_print(f"\n  *** dim(Omega_2) drop varies — hypothesis needs refinement ***")

    if len(d_z1_dist) == 1 and 0 in d_z1_dist:
        flush_print(f"  *** CONFIRMED: dim(Z_1) always stays the same ***")
    else:
        flush_print(f"  *** dim(Z_1) changes — hypothesis needs refinement ***")

    if len(beta1_after_dist) == 1 and 1 in beta1_after_dist:
        flush_print(f"  *** CONFIRMED: beta_1 always becomes 1 ***")
    else:
        flush_print(f"  *** beta_1 varies — hypothesis needs refinement ***")

    flush_print(f"\n  Detailed (first 15):")
    for d in detailed:
        flush_print(f"    mask={d['mask']}, TT={d['tt']}: "
              f"Omega2: {d['orig_Omega2']}->{d['flip_Omega2']} (d={d['d_omega2']}), "
              f"|A2|: {d['orig_A2']}->{d['flip_A2']} (d={d['flip_A2']-d['orig_A2']}), "
              f"Z1: {d['orig_Z1']}->{d['flip_Z1']} (d={d['d_z1']}), "
              f"B1: {d['orig_B1']}->{d['flip_B1']} (d={d['d_b1']}), "
              f"beta1: 0->{d['beta1_after']}")


# ===== PART 3: What happens to |A_2| vs dim(Omega_2)? =====
def analyze_A2_vs_Omega2():
    """Since Omega_2 != A_2, the gap matters. Analyze how it changes on flip."""
    n = 6
    num_edges = n*(n-1)//2
    total = 2**num_edges

    flush_print(f"\n\n{'='*70}")
    flush_print(f"PART 3: |A_2| vs dim(Omega_2) gap analysis on flip (n={n})")
    flush_print(f"{'='*70}")

    t0 = time.time()

    gap_before_dist = {}
    gap_after_dist = {}
    dA2_dist = {}
    dOmega2_dist = {}
    dgap_dist = {}
    count = 0

    for mask in range(total):
        A = adj_matrix(n, mask)
        info = compute_betti_and_dims(A, n, max_dim=1)
        if info[1]['betti'] != 0:
            continue

        tt = find_one_tt(A, n)
        if tt is None:
            continue

        a, b, c = tt
        Af = flip_arc(A, a, c)

        info_orig = compute_betti_and_dims(A, n, max_dim=2)
        info_flip = compute_betti_and_dims(Af, n, max_dim=2)

        gap_before = info_orig[2]['dim_A'] - info_orig[2]['dim_Omega']
        gap_after = info_flip[2]['dim_A'] - info_flip[2]['dim_Omega']
        dA2 = info_flip[2]['dim_A'] - info_orig[2]['dim_A']
        dOmega2 = info_flip[2]['dim_Omega'] - info_orig[2]['dim_Omega']
        dgap = gap_after - gap_before

        gap_before_dist[gap_before] = gap_before_dist.get(gap_before, 0) + 1
        gap_after_dist[gap_after] = gap_after_dist.get(gap_after, 0) + 1
        dA2_dist[dA2] = dA2_dist.get(dA2, 0) + 1
        dOmega2_dist[dOmega2] = dOmega2_dist.get(dOmega2, 0) + 1
        dgap_dist[dgap] = dgap_dist.get(dgap, 0) + 1
        count += 1

    elapsed = time.time() - t0
    flush_print(f"  {count} flips in {elapsed:.1f}s")
    flush_print(f"  |A_2| - dim(Omega_2) BEFORE flip: {dict(sorted(gap_before_dist.items()))}")
    flush_print(f"  |A_2| - dim(Omega_2) AFTER flip:  {dict(sorted(gap_after_dist.items()))}")
    flush_print(f"  Delta |A_2|:         {dict(sorted(dA2_dist.items()))}")
    flush_print(f"  Delta dim(Omega_2):  {dict(sorted(dOmega2_dist.items()))}")
    flush_print(f"  Delta gap:           {dict(sorted(dgap_dist.items()))}")


# ===== PART 4: n=7 sampled =====
def run_n7(num_samples=100):
    n = 7
    num_edges = n*(n-1)//2

    flush_print(f"\n\n{'='*70}")
    flush_print(f"PART 4: n={n} SAMPLED ({num_samples} beta_1=0 tournaments)")
    flush_print(f"{'='*70}")

    random.seed(42)
    t0 = time.time()

    d_omega2_dist = {}
    d_z1_dist = {}
    d_b1_dist = {}
    beta1_after_dist = {}
    candidates = 0
    attempts = 0

    while candidates < num_samples and attempts < 50000:
        mask = random.randint(0, 2**num_edges - 1)
        A = adj_matrix(n, mask)
        attempts += 1

        info_check = compute_betti_and_dims(A, n, max_dim=1)
        if info_check[1]['betti'] != 0:
            continue

        tt = find_one_tt(A, n)
        if tt is None:
            continue

        a, b, c = tt
        Af = flip_arc(A, a, c)

        info_orig = compute_betti_and_dims(A, n, max_dim=2)
        info_flip = compute_betti_and_dims(Af, n, max_dim=2)

        d_omega2 = info_flip[2]['dim_Omega'] - info_orig[2]['dim_Omega']
        d_z1 = info_flip[1]['dim_Z'] - info_orig[1]['dim_Z']
        d_b1 = info_flip[1]['dim_B'] - info_orig[1]['dim_B']
        new_beta1 = info_flip[1]['betti']

        d_omega2_dist[d_omega2] = d_omega2_dist.get(d_omega2, 0) + 1
        d_z1_dist[d_z1] = d_z1_dist.get(d_z1, 0) + 1
        d_b1_dist[d_b1] = d_b1_dist.get(d_b1, 0) + 1
        beta1_after_dist[new_beta1] = beta1_after_dist.get(new_beta1, 0) + 1
        candidates += 1

        if candidates % 20 == 0:
            elapsed = time.time() - t0
            flush_print(f"  Progress: {candidates}/{num_samples}, {elapsed:.1f}s")

    elapsed = time.time() - t0
    flush_print(f"\n  {candidates} tournaments ({attempts} attempts) in {elapsed:.1f}s")
    flush_print(f"\n--- RESULTS (n={n}) ---")
    flush_print(f"  dim(Omega_2) change: {dict(sorted(d_omega2_dist.items()))}")
    flush_print(f"  dim(Z_1) change: {dict(sorted(d_z1_dist.items()))}")
    flush_print(f"  dim(B_1) change: {dict(sorted(d_b1_dist.items()))}")
    flush_print(f"  beta_1 after: {dict(sorted(beta1_after_dist.items()))}")

    if len(d_omega2_dist) == 1 and -1 in d_omega2_dist:
        flush_print(f"\n  *** n={n}: CONFIRMED — dim(Omega_2) always drops by exactly 1 ***")
    else:
        flush_print(f"\n  *** n={n}: dim(Omega_2) drop varies ***")

    if len(d_z1_dist) == 1 and 0 in d_z1_dist:
        flush_print(f"  *** n={n}: CONFIRMED — dim(Z_1) always stays the same ***")
    else:
        flush_print(f"  *** n={n}: dim(Z_1) changes ***")

    if len(beta1_after_dist) == 1 and 1 in beta1_after_dist:
        flush_print(f"  *** n={n}: CONFIRMED — beta_1 always becomes 1 ***")
    else:
        flush_print(f"  *** n={n}: beta_1 varies ***")


# ===== PART 5: Algebraic: why does dim(Omega_2) drop by 1? =====
def algebraic_why():
    flush_print(f"\n\n{'='*70}")
    flush_print("PART 5: ALGEBRAIC ANALYSIS — Why dim(Omega_2) drops by 1")
    flush_print(f"{'='*70}")

    flush_print("""
RECALL: For tournaments, A_1 = all directed edges, so ALL 1-paths (a,b)
are allowed (every pair has an arc). Thus every face of a 2-path is a
valid 1-path, meaning the Omega_2 constraint ("all faces in A_1") is
trivially satisfied.

WAIT: Omega_2 is about 2-CHAINS (elements of A_2). A 2-path is (x,y,z).
Its faces are (y,z), (x,z), (x,y).
  - (y,z): must be in A_1, i.e., y->z is an edge. YES (it's in the TT def).
  - (x,y): must be in A_1, i.e., x->y is an edge. YES.
  - (x,z): must be in A_1, i.e., x->z OR z->x is an edge.

For TOURNAMENTS, one of x->z or z->x exists. But A_1 = {(u,v): u->v in T}.
The face (x,z) is the 1-path from x to z, which requires x->z as an edge.
If z->x instead, then (x,z) is NOT in A_1!

So Omega_2 constraint: for TT (x,y,z) with x->y, y->z, x->z, the face
(x,z) requires x->z, which IS true by definition of TT. But...

Actually, let me recheck: the boundary of (x,y,z) is:
  d(x,y,z) = (y,z) - (x,z) + (x,y)
All three faces are directed paths. For a TT: x->y, y->z, x->z.
  (y,z): y->z ✓ in A_1
  (x,z): x->z ✓ in A_1 (since it's a TT)
  (x,y): x->y ✓ in A_1

So ALL faces of a TT are in A_1. Therefore Omega_2 = A_2 = span of TTs.

But PART 0 showed Omega_2 ≠ A_2!!! Let me debug this.
""")

    # Debug: find a specific case where Omega_2 != A_2
    n = 4
    flush_print(f"  Debugging at n={n}:")
    for mask in range(2**(n*(n-1)//2)):
        A = adj_matrix(n, mask)
        info = compute_betti_and_dims(A, n, max_dim=2)
        if info[2]['dim_Omega'] != info[2]['dim_A']:
            flush_print(f"    mask={mask}: dim_Omega_2={info[2]['dim_Omega']}, |A_2|={info[2]['dim_A']}")
            flush_print(f"    Adjacency: {A}")
            allowed_2 = enumerate_allowed_paths(A, n, 2)
            allowed_1 = enumerate_allowed_paths(A, n, 1)
            flush_print(f"    A_2 (TTs): {allowed_2}")
            flush_print(f"    A_1 (edges): {allowed_1}")
            allowed_1_set = set(allowed_1)
            for tt in allowed_2:
                faces_status = []
                for sign, face in boundary_coeffs(tt):
                    if len(set(face)) == len(face):
                        faces_status.append((face, face in allowed_1_set))
                flush_print(f"    TT {tt}: faces = {faces_status}")
            # Also check what Omega_2 basis looks like
            omega2 = compute_omega_basis(A, n, 2, allowed_2, allowed_1)
            flush_print(f"    Omega_2 basis shape: {omega2.shape}")
            if omega2.shape[1] > 0:
                for col in range(omega2.shape[1]):
                    vec = omega2[:, col]
                    nonzero = [(allowed_2[i], vec[i]) for i in range(len(vec)) if abs(vec[i]) > 1e-10]
                    flush_print(f"    Omega_2 basis[{col}]: {nonzero}")
            break

    flush_print()

    # Hmm, actually I think the issue is that A_2 != just TTs.
    # A_2 = allowed 2-paths = sequences (x,y,z) with x->y AND y->z (edges in T).
    # NOT requiring x->z! A TT additionally requires x->z.
    # But the original code uses "enumerate_allowed_paths" which just follows edges.
    flush_print("  IMPORTANT REALIZATION: A_2 = {(x,y,z) : x->y, y->z in T}")
    flush_print("  This includes BOTH transitive triples (x->z) AND 3-cycles (z->x)!")
    flush_print("  So |A_2| > |TT| in general.")
    flush_print()
    flush_print("  Omega_2 = {u in A_2 : d(u) in A_1}.")
    flush_print("  d(x,y,z) = (y,z) - (x,z) + (x,y).")
    flush_print("  Face (x,z) is NOT an allowed 1-path if z->x (not x->z)!")
    flush_print("  So 3-cycle paths (x,y,z) with z->x have a non-allowed face (x,z).")
    flush_print("  This means: Omega_2 is a PROPER subspace of A_2 when 3-cycles exist.")
    flush_print()

    # Now let's understand: if we flip a TT arc a->c to c->a:
    # - We DESTROY some paths in A_2 (those needing a->c)
    # - We CREATE some paths in A_2 (those using c->a)
    # - Some destroyed paths were in 3-cycles (face (x,z) with z->x not allowed)
    # - Some created paths form new 3-cycles
    # - The Omega_2 constraint changes because the "non-allowed faces" change

    flush_print("  For the flip a->c → c->a:")
    flush_print("  - 2-paths DESTROYED: any (x,y,z) using arc a->c as x->y or y->z")
    flush_print("  - 2-paths CREATED: any (x,y,z) using arc c->a as x->y or y->z")
    flush_print("  - The face (a,c) was in A_1 (since a->c), now (a,c) NOT in A_1 (c->a)")
    flush_print("  - The face (c,a) was NOT in A_1, now IS in A_1")
    flush_print("  - This changes the Omega_2 constraints!")


if __name__ == '__main__':
    check_omega2_vs_A2()
    run_n5()
    run_n6()
    analyze_A2_vs_Omega2()
    run_n7(num_samples=100)
    algebraic_why()
