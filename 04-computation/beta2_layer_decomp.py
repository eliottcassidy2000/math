#!/usr/bin/env python3
"""
beta2_layer_decomp.py — Layer decomposition of the β₂=0 proof

Key question: WHY does DT+cancel always fill Z₂?

Approach: Decompose into clear algebraic layers.
Layer 1: DT_f boundaries span im(∂₃^simp) inside Z₂
Layer 2: DT_c boundaries fill remaining simplicial cycles (flag H₂ holes)
Layer 3: Cancel pairs fill remaining NT directions

This script tests:
1. What is the rank contribution of each layer?
2. What fraction of Z₂ directions are pure TT vs NT vs mixed?
3. 4-vertex sub-tournament DT path structure

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations, permutations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


print("=" * 70)
print("β₂=0 LAYER DECOMPOSITION ANALYSIS")
print("=" * 70)

# =============================================================
# PART 1: Rank contribution of each DT/cancel layer
# =============================================================
print("\nPART 1: Layer-by-layer rank contributions")
print("-" * 50)

for n in [5, 6]:
    m = n*(n-1)//2
    total = 1 << m

    layer_stats = []
    fail_count = 0

    t0 = time.time()
    for bits in range(total):
        A = build_adj(n, bits)

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)

        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

        d2 = dim_om(om2)
        if d2 == 0 or not ap3:
            continue

        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        z2_dim = d2 - rk2
        if z2_dim == 0:
            continue

        bd3 = build_full_boundary_matrix(ap3, ap2)
        ap3_list = [tuple(p) for p in ap3]

        # Classify A₃ paths
        dtf_idx = []
        dtc_idx = []
        nondt_idx = []

        for i, p in enumerate(ap3_list):
            a, b, c, d = p
            if A[a][c] and A[b][d]:
                if A[a][d]:
                    dtf_idx.append(i)
                else:
                    dtc_idx.append(i)
            else:
                nondt_idx.append(i)

        # Build cancel vectors
        bad_groups = defaultdict(list)
        for i in nondt_idx:
            a, b, c, d = ap3_list[i]
            if not A[a][c]: bad_groups[('02', a, c)].append(i)
            if not A[b][d]: bad_groups[('13', b, d)].append(i)

        cancel_vecs_a3 = []
        for key, indices in bad_groups.items():
            for j in range(1, len(indices)):
                v = np.zeros(len(ap3_list))
                v[indices[0]] = 1; v[indices[j]] = -1
                cancel_vecs_a3.append(v)

        # Compute boundary images in Ω₂ coordinates
        def get_coords(col_indices=None, vecs=None):
            if col_indices is not None and len(col_indices) > 0:
                cols = np.column_stack([bd3[:, i] for i in col_indices])
            elif vecs is not None and len(vecs) > 0:
                V = np.column_stack(vecs)
                cols = bd3 @ V
            else:
                return np.zeros((d2, 0))
            return np.linalg.lstsq(om2, cols, rcond=None)[0]

        coords_f = get_coords(col_indices=dtf_idx)
        coords_c = get_coords(col_indices=dtc_idx)
        coords_cancel = get_coords(vecs=cancel_vecs_a3)

        rk_f = np.linalg.matrix_rank(coords_f, tol=1e-8) if coords_f.shape[1] > 0 else 0

        # Combined f+c
        parts_fc = [p for p in [coords_f, coords_c] if p.shape[1] > 0]
        rk_fc = np.linalg.matrix_rank(np.hstack(parts_fc), tol=1e-8) if parts_fc else 0

        # All
        parts_all = [p for p in [coords_f, coords_c, coords_cancel] if p.shape[1] > 0]
        rk_all = np.linalg.matrix_rank(np.hstack(parts_all), tol=1e-8) if parts_all else 0

        contrib_c = rk_fc - rk_f
        contrib_cancel = rk_all - rk_fc

        layer_stats.append((z2_dim, rk_f, contrib_c, contrib_cancel, rk_all))

        if rk_all < z2_dim:
            fail_count += 1

    elapsed = time.time() - t0
    print(f"\nn={n}: {len(layer_stats)} nontrivial tournaments ({elapsed:.0f}s)")

    if not layer_stats:
        continue

    print(f"  Z₂ dims: {Counter([s[0] for s in layer_stats])}")
    print(f"  DT_f rank: {Counter([s[1] for s in layer_stats])}")
    print(f"  DT_c incremental: {Counter([s[2] for s in layer_stats])}")
    print(f"  Cancel incremental: {Counter([s[3] for s in layer_stats])}")

    deficit_dist = Counter()
    for z, rf, cc, ccl, ra in layer_stats:
        deficit_dist[z - ra] += 1
    print(f"  Deficit (Z₂ - total rank): {dict(sorted(deficit_dist.items()))}")

    need_c = sum(1 for _, rf, cc, _, _ in layer_stats if cc > 0)
    need_cancel = sum(1 for _, _, _, ccl, _ in layer_stats if ccl > 0)
    print(f"  Need DT_c: {need_c}/{len(layer_stats)} ({100*need_c/len(layer_stats):.1f}%)")
    print(f"  Need cancel: {need_cancel}/{len(layer_stats)} ({100*need_cancel/len(layer_stats):.1f}%)")
    print(f"  Failures: {fail_count}")


# =============================================================
# PART 2: TT vs NT content of Z₂ directions
# =============================================================
print(f"\n{'='*70}")
print("PART 2: TT vs NT content of Z₂ directions at n=5")
print("-" * 50)

n = 5
m = n*(n-1)//2
directions = []

for bits in range(1 << m):
    A = build_adj(n, bits)

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

    d2 = dim_om(om2)
    if d2 == 0:
        continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0:
        continue

    ap2_list = [tuple(p) for p in ap2]
    tt_mask = np.array([1.0 if A[a][c] else 0.0 for a, b, c in ap2_list])

    U, S, Vt = np.linalg.svd(coords2, full_matrices=True)
    z2_basis = Vt[rk2:]

    for k in range(z2_dim):
        z2_dir_a2 = om2 @ z2_basis[k]
        tt_norm = np.linalg.norm(z2_dir_a2 * tt_mask)
        total_norm = np.linalg.norm(z2_dir_a2)
        if total_norm > 1e-10:
            tt_frac = tt_norm / total_norm
        else:
            tt_frac = 0
        directions.append(tt_frac)

print(f"Total Z₂ directions: {len(directions)}")
print(f"TT fraction: min={min(directions):.4f}, max={max(directions):.4f}, mean={np.mean(directions):.4f}")

pure_tt = sum(1 for f in directions if f > 0.99)
pure_nt = sum(1 for f in directions if f < 0.01)
mixed = len(directions) - pure_tt - pure_nt
print(f"Pure TT (>99%): {pure_tt}, Pure NT (<1%): {pure_nt}, Mixed: {mixed}")

# Histogram of TT fractions
bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01]
hist, _ = np.histogram(directions, bins=bins)
print(f"\nTT fraction histogram:")
for i in range(len(hist)):
    bar = '#' * (hist[i] // 5)
    print(f"  [{bins[i]:.1f}, {bins[i+1]:.1f}): {hist[i]:5d}  {bar}")


# =============================================================
# PART 3: 4-vertex sub-tournament DT path structure
# =============================================================
print(f"\n{'='*70}")
print("PART 3: 4-vertex subtournament DT path counts")
print("-" * 50)

# Enumerate all 4-vertex tournament types
for bits4 in range(1 << 6):  # C(4,2)=6 bits
    A4 = build_adj(4, bits4)
    n_dtf = 0
    n_dtc = 0
    n_nondt = 0

    for perm in permutations(range(4)):
        a, b, c, d = perm
        if A4[a][b] and A4[b][c] and A4[c][d]:
            if A4[a][c] and A4[b][d]:
                if A4[a][d]:
                    n_dtf += 1
                else:
                    n_dtc += 1
            else:
                n_nondt += 1

    c3 = sum(1 for trip in combinations(range(4), 3)
             if sum(A4[trip[i]][trip[j]] for i in range(3) for j in range(3) if i!=j) == 3)

    scores = tuple(sorted(sum(A4[i][j] for j in range(4) if j!=i) for i in range(4)))

    # Only print unique score sequences
    if bits4 == 0 or True:
        pass

# Actually, let me just do it by score sequence
score_types = defaultdict(list)
for bits4 in range(1 << 6):
    A4 = build_adj(4, bits4)
    scores = tuple(sorted(sum(A4[i][j] for j in range(4) if j!=i) for i in range(4)))

    n_dtf = n_dtc = n_nondt = 0
    for perm in permutations(range(4)):
        a, b, c, d = perm
        if A4[a][b] and A4[b][c] and A4[c][d]:
            if A4[a][c] and A4[b][d]:
                if A4[a][d]: n_dtf += 1
                else: n_dtc += 1
            else:
                n_nondt += 1

    c3 = sum(1 for trip in combinations(range(4), 3)
             if sum(A4[trip[i]][trip[j]] for i in range(3) for j in range(3) if i!=j) == 3)

    key = (scores, c3)
    if not score_types[key]:
        score_types[key] = (n_dtf, n_dtc, n_nondt)

print(f"{'scores':<20} {'c3':>3} {'DT_f':>5} {'DT_c':>5} {'nonDT':>6} {'|A3|':>5}")
for key in sorted(score_types.keys()):
    scores, c3 = key
    dtf, dtc, ndt = score_types[key]
    print(f"  {str(scores):<18} {c3:>3} {dtf:>5} {dtc:>5} {ndt:>6} {dtf+dtc+ndt:>5}")


# =============================================================
# PART 4: KEY TEST — Does ∂₃(DT_c) always have the right
# structure to fill simplicial holes?
# =============================================================
print(f"\n{'='*70}")
print("PART 4: DT_c boundary structure — does it always cover simp holes?")
print("-" * 50)

n = 5
m = n*(n-1)//2

dtc_always_fills_simp = True
total_with_simp_hole = 0

for bits in range(1 << m):
    A = build_adj(n, bits)

    # Quick flag H₂ check
    tt_triples = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c in (a,b) or not A[b][c]: continue
                if A[a][c]:
                    tt_triples.append((a, b, c))

    if not tt_triples:
        continue

    arcs = [(i, j) for i in range(n) for j in range(n) if i != j and A[i][j]]
    arc_idx = {a: i for i, a in enumerate(arcs)}

    d2_simp = np.zeros((len(arcs), len(tt_triples)))
    for j, (a, b, c) in enumerate(tt_triples):
        if (b, c) in arc_idx: d2_simp[arc_idx[(b, c)], j] += 1
        if (a, c) in arc_idx: d2_simp[arc_idx[(a, c)], j] -= 1
        if (a, b) in arc_idx: d2_simp[arc_idx[(a, b)], j] += 1

    rk_d2_simp = np.linalg.matrix_rank(d2_simp, tol=1e-8)
    z2_simp = len(tt_triples) - rk_d2_simp

    # Transitive 4-subsets
    trans4 = []
    for quad in combinations(range(n), 4):
        has_cycle = False
        for trip in combinations(quad, 3):
            total = sum(A[trip[i]][trip[j]] for i in range(3) for j in range(3) if i != j)
            if total == 3:
                has_cycle = True
                break
        if not has_cycle:
            sc = [(sum(A[v][w] for w in quad if w != v), v) for v in quad]
            sc.sort(reverse=True)
            trans4.append(tuple(v for _, v in sc))

    tt_idx = {t: i for i, t in enumerate(tt_triples)}
    if trans4:
        d3_simp = np.zeros((len(tt_triples), len(trans4)))
        for j, (a, b, c, d) in enumerate(trans4):
            faces = [(1, (b,c,d)), (-1, (a,c,d)), (1, (a,b,d)), (-1, (a,b,c))]
            for sign, face in faces:
                if face in tt_idx:
                    d3_simp[tt_idx[face], j] += sign
        rk_d3_simp = np.linalg.matrix_rank(d3_simp, tol=1e-8)
    else:
        d3_simp = np.zeros((len(tt_triples), 0))
        rk_d3_simp = 0

    h2_simp = z2_simp - rk_d3_simp

    if h2_simp == 0:
        continue

    total_with_simp_hole += 1

    # Now check: do DT_c boundaries cover the simplicial hole?
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap3:
        dtc_always_fills_simp = False
        continue

    bd3 = build_full_boundary_matrix(ap3, enumerate_allowed_paths(A, n, 2))
    ap3_list = [tuple(p) for p in ap3]

    dtc_idx = [i for i, p in enumerate(ap3_list)
               if A[p[0]][p[2]] and A[p[1]][p[3]] and not A[p[0]][p[3]]]

    if not dtc_idx:
        dtc_always_fills_simp = False
        continue

    # Get DT_c boundaries restricted to TT triples
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap2_list = [tuple(p) for p in ap2]
    tt_in_ap2 = {t: ap2_list.index(t) for t in tt_triples if t in ap2_list}

    # DT_c boundary projected onto TT components
    bd3_dtc = np.column_stack([bd3[:, i] for i in dtc_idx])
    tt_proj = np.zeros((len(tt_triples), len(dtc_idx)))
    for j_tt, triple in enumerate(tt_triples):
        if triple in tt_in_ap2:
            idx_a2 = tt_in_ap2[triple]
            tt_proj[j_tt, :] = bd3_dtc[idx_a2, :]

    # Find the simplicial hole directions (null of d2_simp not in im(d3_simp))
    U, S_vals, Vt = np.linalg.svd(d2_simp, full_matrices=True)
    z2_simp_basis = Vt[rk_d2_simp:].T

    if d3_simp.shape[1] > 0:
        proj = z2_simp_basis.T @ d3_simp
        U2, S2, _ = np.linalg.svd(proj, full_matrices=True)
        uncovered_start = sum(s > 1e-8 for s in S2)
        hole_dirs = z2_simp_basis @ U2[:, uncovered_start:]
    else:
        hole_dirs = z2_simp_basis

    # Check if DT_c TT-projections cover the hole directions
    dtc_tt_proj = tt_proj  # each column is a DT_c boundary in TT coords
    # Test: are hole_dirs in column span of dtc_tt_proj?
    for k in range(hole_dirs.shape[1]):
        hole = hole_dirs[:, k]
        x, res, _, _ = np.linalg.lstsq(dtc_tt_proj, hole, rcond=None)
        err = np.max(np.abs(dtc_tt_proj @ x - hole))
        if err > 1e-6:
            dtc_always_fills_simp = False
            if total_with_simp_hole <= 3:
                print(f"  FAIL: bits={bits}, hole not in DT_c TT-span, err={err:.2e}")

print(f"\nTotal with simplicial holes: {total_with_simp_hole}")
print(f"DT_c always fills simplicial holes: {dtc_always_fills_simp}")


print("\nDone.")
