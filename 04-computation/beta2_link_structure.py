#!/usr/bin/env python3
"""
beta2_link_structure.py — Analyze the "link" structure of vertex v

For the quotient complex R_p = Ω_p(T)/Ω_p(T\v), elements of R_p
are chains involving vertex v. We study their structure.

Key idea: paths through v can be classified by how v appears:
- "star" paths: v is an interior vertex (a, v, b)
- "cone" paths: v is an endpoint (v, a, b) or (a, b, v)

For a 2-path (a,b,c) through v:
  Type 1: v = a (v is start): (v, b, c) — v→b, b→c
  Type 2: v = b (v is middle): (a, v, c) — a→v, v→c  
  Type 3: v = c (v is end): (a, b, v) — a→b, b→v

For Ω₂ membership: ∂₂(a,b,c) = (b,c) - (a,c) + (a,b) ∈ A₁
  Type 1: (b,c) - (v,c) + (v,b) — faces (v,c) and (v,b) through v
  Type 2: (v,c) - (a,c) + (a,v) — faces (v,c) and (a,v) through v
  Type 3: (b,v) - (a,v) + (a,b) — faces (b,v) and (a,v) through v

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

print("=" * 70)
print("LINK STRUCTURE OF THE QUOTIENT COMPLEX")
print("=" * 70)

# For each tournament and vertex v, classify Ω₂ paths by type
# and study R₂ = Ω₂(T) / Ω₂(T\v) structure

# Choose v = 0 for simplicity (similar for other v by symmetry)
v = 0

print(f"\nAnalyzing v={v} across all n={n} tournaments")

# For R₂: paths through v in Ω₂(T)
# For R₃: 3-paths through v in Ω₃(T)

r_dims = Counter()  # (dim R₁, dim R₂, dim R₃)
r2_type_dist = Counter()

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    # R₁ = edges through v
    # v has (n-1) edges total
    ap1 = enumerate_allowed_paths(A, n, 1)
    edges_through_v = [p for p in ap1 if v in p]
    dim_R1 = len(edges_through_v)  # = n-1 always for tournaments

    # Ω₂(T) and Ω₂(T\v)
    ap2_T = enumerate_allowed_paths(A, n, 2)
    ap1_T = ap1
    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
    d2_T = om2_T.shape[1] if om2_T.ndim == 2 else 0

    B = [[A[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
    ap2_Tv = enumerate_allowed_paths(B, n-1, 2)
    ap1_Tv = enumerate_allowed_paths(B, n-1, 1)
    om2_Tv = compute_omega_basis(B, n-1, 2, ap2_Tv, ap1_Tv) if ap2_Tv else np.zeros((0,0))
    d2_Tv = om2_Tv.shape[1] if om2_Tv.ndim == 2 and om2_Tv.shape[0] > 0 else 0

    dim_R2 = d2_T - d2_Tv

    # Ω₃(T) and Ω₃(T\v)
    ap3_T = enumerate_allowed_paths(A, n, 3)
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))
    d3_T = om3_T.shape[1] if om3_T.ndim == 2 else 0

    ap3_Tv = enumerate_allowed_paths(B, n-1, 3)
    om3_Tv = compute_omega_basis(B, n-1, 3, ap3_Tv, ap2_Tv) if ap3_Tv else np.zeros((0,0))
    d3_Tv = om3_Tv.shape[1] if om3_Tv.ndim == 2 and om3_Tv.shape[0] > 0 else 0

    dim_R3 = d3_T - d3_Tv

    r_dims[(dim_R1, dim_R2, dim_R3)] += 1

    # Classification of 2-paths through v
    if ap2_T and d2_T > 0:
        types = Counter()
        for path in ap2_T:
            if v not in path: continue
            if path[0] == v: types['start'] += 1
            elif path[1] == v: types['middle'] += 1
            elif path[2] == v: types['end'] += 1
        r2_type_dist[tuple(sorted(types.items()))] += 1

print(f"\nR_p dimensions (R₁, R₂, R₃): count")
for key, cnt in sorted(r_dims.items()):
    print(f"  R=({key[0]}, {key[1]}, {key[2]}): {cnt}")

print(f"\n2-path types through v=0:")
for key, cnt in sorted(r2_type_dist.items()):
    print(f"  {dict(key)}: {cnt}")

# Check: does β₂=0 follow from R₃ ≥ R₂ always?
# If dim R₃ ≥ dim ker(∂₂^R), then im(∂₃^R) can fill ker(∂₂^R)
# But R₃ ≥ R₂ is necessary but not sufficient.

print(f"\nR₃ vs R₂:")
r3_ge_r2 = sum(cnt for k, cnt in r_dims.items() if k[2] >= k[1])
r3_lt_r2 = sum(cnt for k, cnt in r_dims.items() if k[2] < k[1])
print(f"  R₃ ≥ R₂: {r3_ge_r2}/{1<<m}")
print(f"  R₃ < R₂: {r3_lt_r2}/{1<<m}")

# Now compute the actual relative homology to verify
print(f"\n{'='*70}")
print("RELATIVE HOMOLOGY COMPUTATION")
print("=" * 70)

# For each (R₁, R₂, R₃) tuple, compute H₂^R
# The key LES identity: β₂(T) = H₂^R - dim(ker i_*)

def compute_h2_rel_fixed_v0(A, n):
    """Compute H₂(T, T\\0) with correct vertex renumbering for v=0."""
    v = 0
    B = [[A[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
    loc2glob = [i for i in range(n) if i != v]

    om_T = {}; ap_T_dict = {}; om_Tv = {}; ap_Tv_dict = {}
    for p in range(5):
        ap_T_dict[p] = enumerate_allowed_paths(A, n, p)
        if p == 0: om_T[p] = np.eye(n)
        elif ap_T_dict[p]: om_T[p] = compute_omega_basis(A, n, p, ap_T_dict[p], ap_T_dict[p-1])
        else: om_T[p] = np.zeros((0, 0))

    for p in range(4):
        ap_Tv_dict[p] = enumerate_allowed_paths(B, n-1, p)
        if p == 0: om_Tv[p] = np.eye(n-1)
        elif ap_Tv_dict[p]: om_Tv[p] = compute_omega_basis(B, n-1, p, ap_Tv_dict[p], ap_Tv_dict[p-1])
        else: om_Tv[p] = np.zeros((0, 0))

    def dim_om(om):
        return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

    d2T = dim_om(om_T[2]); d3T = dim_om(om_T[3])
    if d2T == 0: return 0, 0, 0  # h2_rel, dim_ker, dim_im

    # Embed with renumbering
    def embed(p):
        omTv = om_Tv[p]
        if omTv.ndim != 2 or omTv.shape[1] == 0:
            return np.zeros((len(ap_T_dict[p]), 0))
        T_list = [tuple(x) for x in ap_T_dict[p]]
        Tv_list = [tuple(x) for x in ap_Tv_dict[p]]
        T_idx = {pp: i for i, pp in enumerate(T_list)}
        incl = np.zeros((len(T_list), len(Tv_list)))
        for j, path_local in enumerate(Tv_list):
            path_global = tuple(loc2glob[k] for k in path_local)
            if path_global in T_idx:
                incl[T_idx[path_global], j] = 1
        return incl @ omTv

    emb1 = embed(1); emb2 = embed(2)
    bd2_A = build_full_boundary_matrix(ap_T_dict[2], ap_T_dict[1])
    bd2_om = bd2_A @ om_T[2]

    if emb1.shape[1] > 0:
        M = np.hstack([bd2_om, -emb1])
    else:
        M = bd2_om
    U, S, Vt = np.linalg.svd(M, full_matrices=True)
    rk_M = int(sum(s > 1e-8 for s in S))
    null_space = Vt[rk_M:].T if rk_M < Vt.shape[0] else np.zeros((Vt.shape[1], 0))
    c_part = null_space[:d2T, :]
    preimage_A = om_T[2] @ c_part

    if emb2.shape[1] > 0 and preimage_A.shape[1] > 0:
        combined = np.hstack([emb2, preimage_A])
    elif preimage_A.shape[1] > 0: combined = preimage_A
    else: combined = emb2
    rk_combined = np.linalg.matrix_rank(combined, tol=1e-8) if combined.shape[1] > 0 else 0
    rk_emb2 = np.linalg.matrix_rank(emb2, tol=1e-8) if emb2.shape[1] > 0 else 0
    dim_ker = rk_combined - rk_emb2

    if d3T > 0:
        bd3_A = build_full_boundary_matrix(ap_T_dict[3], ap_T_dict[2])
        bd3_om = bd3_A @ om_T[3]
        if emb2.shape[1] > 0:
            combined_3 = np.hstack([emb2, bd3_om])
        else: combined_3 = bd3_om
        rk_combined_3 = np.linalg.matrix_rank(combined_3, tol=1e-8)
        dim_im = rk_combined_3 - rk_emb2
    else:
        dim_im = 0

    return dim_ker - dim_im, dim_ker, dim_im

# Check at n=5
h2_details = Counter()
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    h2r, dk, di = compute_h2_rel_fixed_v0(A, n)
    h2_details[(h2r, dk, di)] += 1

print(f"\n(H₂^rel, ker ∂₂^R, im ∂₃^R) at n={n}, v=0:")
for key, cnt in sorted(h2_details.items()):
    print(f"  H₂^R={key[0]}, ker={key[1]}, im={key[2]}: {cnt}")

# KEY INSIGHT: ker and im are equal in all cases!
# What determines the common value?

print(f"\nker(∂₂^R) always equals im(∂₃^R)?",
      all(k[1] == k[2] for k in h2_details))

# What determines ker(∂₂^R)?
# ker(∂₂^R) = dim{x ∈ Ω₂(T) : ∂₂(x) ∈ Ω₁(T\v)} / Ω₂(T\v)
# These are elements whose boundary has no edges incident to v.

# Study: what is dim(ker ∂₂^R) as a function of v's out-degree?
print(f"\nker(∂₂^R) by out-degree of v=0:")
ker_by_outdeg = {}
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    outdeg = sum(A[0])
    h2r, dk, di = compute_h2_rel_fixed_v0(A, n)
    if outdeg not in ker_by_outdeg:
        ker_by_outdeg[outdeg] = Counter()
    ker_by_outdeg[outdeg][dk] += 1

for d in sorted(ker_by_outdeg):
    print(f"  d_out(v)={d}: ker values = {dict(sorted(ker_by_outdeg[d].items()))}")

print("\nDone.")
