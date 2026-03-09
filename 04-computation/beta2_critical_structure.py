#!/usr/bin/env python3
"""
beta2_critical_structure.py — Deep structure of critical filler vertices

For each tournament T with Sigma = sum h2_rel(v):
- Which vertices are critical (h2_rel = 1)?
- What is the filling 2-chain sigma_v with d2(sigma_v) = z_v?
- Do different critical vertices share 2-paths in their fillings?
- How do the filling 2-chains relate to each other?

Key hypothesis: the filling 2-chains for different critical vertices
must OVERLAP, limiting how many can coexist. If we can show they
share a "bottleneck", we get the bound Sigma <= 3.

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
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

def compute_h2rel(A, n, v, ap_T, om_T, bd_T):
    """Compute h2_rel(T, T\\v) using dimensions."""
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

    ap0s = enumerate_allowed_paths(A_sub, n1, 0)
    ap1s = enumerate_allowed_paths(A_sub, n1, 1)
    ap2s = enumerate_allowed_paths(A_sub, n1, 2)
    om1s = compute_omega_basis(A_sub, n1, 1, ap1s, ap0s)
    om2s = compute_omega_basis(A_sub, n1, 2, ap2s, ap1s) if ap2s else np.zeros((0,0))
    d1s = dim_om(om1s)
    d2s = dim_om(om2s)

    if d1s == 0:
        return 0, 0, None

    bd1s = build_full_boundary_matrix(ap1s, ap0s)
    rk1s = np.linalg.matrix_rank(bd1s @ om1s, tol=1e-8)
    z1s = d1s - rk1s

    if d2s > 0:
        bd2s = build_full_boundary_matrix(ap2s, ap1s)
        bd2oms = np.linalg.lstsq(om1s, bd2s @ om2s, rcond=None)[0]
        b1s = np.linalg.matrix_rank(bd2oms, tol=1e-8)
    else:
        b1s = 0

    beta1_sub = z1s - b1s

    if beta1_sub == 0:
        return 0, 0, None

    # Compute inclusion map i*: H1(T\v) -> H1(T)
    # H1(T\v) generator: find a 1-cycle not in B1(T\v)
    U, S, Vt = np.linalg.svd(bd1s @ om1s, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z1_basis = Vt[rk:, :]  # in Om1(T\v) coords

    if b1s > 0:
        # Project out B1
        b1_in_z1 = z1_basis @ bd2oms
        U_b, S_b, _ = np.linalg.svd(b1_in_z1, full_matrices=True)
        b_rk = sum(s > 1e-8 for s in S_b)
        h1_gen = U_b[:, b_rk:][:, 0]  # in Z1 coords
    else:
        h1_gen = z1_basis[0]  # first Z1 vector

    # H1 generator in A1(T\v) coords
    h1_A1_sub = om1s @ (z1_basis.T @ h1_gen)

    # Embed into A1(T)
    ap1_T = ap_T[1]
    ap1_T_list = [tuple(p) for p in ap1_T]
    h1_T = np.zeros(len(ap1_T_list))
    remap = {i: others[i] for i in range(n1)}
    for i, p in enumerate(ap1s):
        path_T = (remap[p[0]], remap[p[1]])
        if path_T in ap1_T_list:
            h1_T[ap1_T_list.index(path_T)] = h1_A1_sub[i]

    # Check if h1_T is in B1(T)
    om1_T = om_T[1]
    om2_T = om_T[2]
    d2_T = dim_om(om2_T)

    h1_om = np.linalg.lstsq(om1_T, h1_T, rcond=None)[0]

    if d2_T > 0:
        bd2_T = bd_T[2]
        bd2om_T = np.linalg.lstsq(om1_T, bd2_T @ om2_T, rcond=None)[0]
        preimage = np.linalg.lstsq(bd2om_T, h1_om, rcond=None)[0]
        resid = np.linalg.norm(bd2om_T @ preimage - h1_om)
        if resid < 1e-6:
            # h1 is a boundary in T — h2_rel = 1
            # The filling 2-chain is om2_T @ preimage (in A2(T) coords)
            filling = om2_T @ preimage
            return 1, beta1_sub, filling
        else:
            return 0, beta1_sub, None
    else:
        return 0, beta1_sub, None


print("=" * 70)
print("CRITICAL FILLER STRUCTURE")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# For Sigma=3 tournaments: analyze critical filler structure
sigma3_count = 0
filling_overlap_stats = []

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    # Pre-compute T data
    ap_T = {}
    om_T = {}
    bd_T = {}
    for p in range(5):
        ap_T[p] = enumerate_allowed_paths(A, n, p)
    om_T[1] = compute_omega_basis(A, n, 1, ap_T[1], ap_T[0])
    om_T[2] = compute_omega_basis(A, n, 2, ap_T[2], ap_T[1]) if ap_T[2] else np.zeros((0,0))
    bd_T[2] = build_full_boundary_matrix(ap_T[2], ap_T[1])

    # Compute h2_rel for each vertex
    critical_vertices = []
    fillings = {}
    for v in range(n):
        h2r, b1s, filling = compute_h2rel(A, n, v, ap_T, om_T, bd_T)
        if h2r > 0:
            critical_vertices.append(v)
            fillings[v] = filling

    sigma = len(critical_vertices)
    if sigma != 3:
        continue

    sigma3_count += 1
    if sigma3_count > 5:
        continue  # Only detail first 5

    print(f"\nT bits={bits}, scores={scores}, Sigma=3")
    print(f"  Critical vertices: {critical_vertices}")

    ap2_list = [tuple(p) for p in ap_T[2]]

    for v in critical_vertices:
        filling = fillings[v]
        fill_paths = []
        for i, p in enumerate(ap2_list):
            if abs(filling[i]) > 1e-8:
                has_v = v in p
                tt = "TT" if A[p[0]][p[2]] else "NT"
                fill_paths.append((p, round(filling[i], 4), has_v, tt))

        n_vpaths = sum(1 for _, _, hv, _ in fill_paths if hv)
        n_nonv = sum(1 for _, _, hv, _ in fill_paths if not hv)
        print(f"  v={v} (d+={scores[v]}): {len(fill_paths)} filling paths "
              f"({n_vpaths} v-paths, {n_nonv} non-v)")
        for p, c, hv, tt in fill_paths:
            print(f"    {c:+.4f} * {p} [{'V' if hv else ' '}, {tt}]")

    # Overlap analysis: do fillings share 2-paths?
    for v1, v2 in combinations(critical_vertices, 2):
        f1 = fillings[v1]
        f2 = fillings[v2]
        # Find common nonzero paths
        common = 0
        for i in range(len(ap2_list)):
            if abs(f1[i]) > 1e-8 and abs(f2[i]) > 1e-8:
                common += 1
        supp1 = sum(1 for x in f1 if abs(x) > 1e-8)
        supp2 = sum(1 for x in f2 if abs(x) > 1e-8)
        filling_overlap_stats.append((common, supp1, supp2))
        print(f"  Overlap({v1},{v2}): {common} shared paths (|f1|={supp1}, |f2|={supp2})")

print(f"\nTotal Sigma=3 tournaments: {sigma3_count}")
print(f"Overlap statistics (all Sigma=3 pairs):")
if filling_overlap_stats:
    overlaps = [o for o, _, _ in filling_overlap_stats]
    print(f"  Mean overlap: {np.mean(overlaps):.2f}")
    print(f"  Min/Max overlap: {min(overlaps)}/{max(overlaps)}")

# Key question: can we show that 4 critical fillers would require
# more 2-paths than exist?
print(f"\n{'='*70}")
print("COUNTING ARGUMENT")
print("=" * 70)

# For each n=5 tournament: count available v-triples per vertex
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    # Count TT triples through each vertex
    tt_through = {}
    for v in range(n):
        count_v = 0
        for p in ap_T[2]:
            if A[p[0]][p[2]] and v in p:  # TT path through v
                count_v += 1
        tt_through[v] = count_v

    # Count how many 2-paths are ONLY through specific vertices
    # (i.e., path (a,b,c) where the only critical vertex is one specific v)

    break  # Just show example

print(f"\nFirst tournament TT-triples through each vertex:")
print(f"  {tt_through}")

print("\nDone.")
