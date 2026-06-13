#!/usr/bin/env python3
"""
beta2_why_sigma_leq3.py — WHY is Σ β₁(T\v) ≤ 3 when β₁(T) = 0?

The filling space analysis:
When β₁(T) = 0 and β₁(T\v) = 1, vertex v "fills" the cycle in T\v.
The filling uses v-paths from Ω₂(T) (not in Ω₂(T\v)).

Key idea: the filling chains for different critical vertices must
lie in DIFFERENT parts of Ω₂(T). If they compete for the same 
Ω₂ resources, we get a dimension bound.

More precisely: the fillings σ_v ∈ Ω₂(T) with ∂₂(σ_v) = z_v 
are constrained by:
  1. σ_v must use v-paths (otherwise z_v ∈ B₁(T\v), contradicting β₁=1)
  2. The boundaries ∂₂(σ_v) are linearly independent in Z₁

This script analyzes:
1. The dimension of the "filling space" for each critical vertex
2. The intersection structure of filling spaces
3. Whether a rank argument bounds the number of fillers

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations
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

def compute_beta1(A, n):
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d1 = dim_om(om1)
    if d1 == 0: return 0, None, None, None, None
    bd1 = build_full_boundary_matrix(ap1, ap0)
    rk1 = np.linalg.matrix_rank(bd1 @ om1, tol=1e-8)
    z1 = d1 - rk1
    d2 = dim_om(om2)
    if d2 > 0:
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2om = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
        b1 = np.linalg.matrix_rank(bd2om, tol=1e-8)
    else:
        b1 = 0
        bd2 = None
        bd2om = None
    return z1 - b1, om1, om2, ap1, ap2


# ============================================================
# Analysis: Filling space structure at n=5
# ============================================================
print("=" * 70)
print("FILLING SPACE ANALYSIS AT n=5")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# For Σ=3 tournaments: analyze filling space dimensions
sigma3_examples = []

for bits in range(total):
    A = build_adj(n, bits)
    
    # Compute T data
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    
    d2 = dim_om(om2)
    if d2 == 0: continue
    
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2om = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
    
    beta1_T = 0  # Quick check
    bd1 = build_full_boundary_matrix(ap1, ap0)
    rk1 = np.linalg.matrix_rank(bd1 @ om1, tol=1e-8)
    d1 = dim_om(om1)
    b1 = np.linalg.matrix_rank(bd2om, tol=1e-8) if d2 > 0 else 0
    beta1_T = (d1 - rk1) - b1
    
    if beta1_T != 0: continue
    
    # Check each vertex deletion
    critical = []
    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
        b1v, om1s, om2s, ap1s, ap2s = compute_beta1(A_sub, n-1)
        if b1v == 1:
            critical.append(v)
    
    sigma = len(critical)
    if sigma != 3: continue
    
    if len(sigma3_examples) >= 3:
        continue
    
    # Deep analysis of filling space
    # For each critical v: find the H₁ generator z_v, then the filling σ_v
    
    # First: Z₁ and B₁ of T
    U1, S1, Vt1 = np.linalg.svd(bd1 @ om1, full_matrices=True)
    rk1_T = sum(s > 1e-8 for s in S1)
    z1_basis = Vt1[rk1_T:, :]  # In Ω₁ coords
    z1_dim = z1_basis.shape[0]
    
    # B₁ = im(∂₂|Ω₂) in Ω₁ coords
    b1_basis = bd2om  # Each column is ∂₂(ω₂_j) in Ω₁ coords
    b1_rk = np.linalg.matrix_rank(b1_basis, tol=1e-8)
    
    print(f"\nT bits={bits}, scores={[sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]}")
    print(f"  dim Ω₁={d1}, dim Z₁={z1_dim}, dim B₁={b1_rk}, dim Ω₂={d2}")
    print(f"  Critical vertices: {critical}")
    
    fillings = {}
    h1_gens = {}
    
    for v in critical:
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
        
        ap0s = enumerate_allowed_paths(A_sub, n-1, 0)
        ap1s = enumerate_allowed_paths(A_sub, n-1, 1)
        ap2s = enumerate_allowed_paths(A_sub, n-1, 2)
        om1s = compute_omega_basis(A_sub, n-1, 1, ap1s, ap0s)
        om2s = compute_omega_basis(A_sub, n-1, 2, ap2s, ap1s) if ap2s else np.zeros((0,0))
        
        bd1s = build_full_boundary_matrix(ap1s, ap0s)
        U_s, S_s, Vt_s = np.linalg.svd(bd1s @ om1s, full_matrices=True)
        rk_s = sum(s > 1e-8 for s in S_s)
        z1s_basis = Vt_s[rk_s:, :]  # Z₁(T\v) in Ω₁(T\v) coords
        
        d2s = dim_om(om2s)
        if d2s > 0:
            bd2s = build_full_boundary_matrix(ap2s, ap1s)
            bd2oms = np.linalg.lstsq(om1s, bd2s @ om2s, rcond=None)[0]
            b1s_rk = np.linalg.matrix_rank(bd2oms, tol=1e-8)
            # Project out B₁(T\v) from Z₁(T\v)
            if b1s_rk > 0:
                b1s_in_z1 = z1s_basis @ bd2oms
                U_b, S_b, _ = np.linalg.svd(b1s_in_z1, full_matrices=True)
                b_rk = sum(s > 1e-8 for s in S_b)
                h1_gen = U_b[:, b_rk:][:, 0]
            else:
                h1_gen = z1s_basis[0]
        else:
            h1_gen = z1s_basis[0]
        
        # H₁ generator in A₁(T\v) coords
        h1_A1_sub = om1s @ (z1s_basis.T @ h1_gen)
        
        # Embed into Ω₁(T)
        ap1_T_list = [tuple(p) for p in ap1]
        h1_T = np.zeros(len(ap1_T_list))
        remap = {i: others[i] for i in range(n-1)}
        for i, p in enumerate(ap1s):
            path_T = (remap[p[0]], remap[p[1]])
            if path_T in ap1_T_list:
                h1_T[ap1_T_list.index(path_T)] = h1_A1_sub[i]
        
        # Express in Ω₁ coords
        h1_om = np.linalg.lstsq(om1, h1_T, rcond=None)[0]
        
        # Find the filling: σ ∈ Ω₂(T) with ∂₂(σ) = z_v
        # bd2om @ σ_om = h1_om (in Ω₁ coords)
        sigma_om = np.linalg.lstsq(bd2om, h1_om, rcond=None)[0]
        resid = np.linalg.norm(bd2om @ sigma_om - h1_om)
        
        # σ in A₂(T) coords
        sigma_A2 = om2 @ sigma_om
        
        fillings[v] = sigma_A2
        h1_gens[v] = h1_T
        
        # Decompose σ: which paths are v-paths vs non-v-paths?
        ap2_list = [tuple(p) for p in ap2]
        v_support = sum(1 for i, p in enumerate(ap2_list) if abs(sigma_A2[i]) > 1e-8 and v in p)
        nonv_support = sum(1 for i, p in enumerate(ap2_list) if abs(sigma_A2[i]) > 1e-8 and v not in p)
        
        print(f"  v={v}: filling has {v_support} v-paths, {nonv_support} non-v-paths, resid={resid:.2e}")
    
    # Key test: are the H₁ generators linearly independent?
    h1_matrix = np.column_stack([h1_gens[v] for v in critical])
    h1_rk = np.linalg.matrix_rank(h1_matrix, tol=1e-8)
    print(f"  H₁ generators rank: {h1_rk} (of {len(critical)} critical vertices)")
    
    # Are the fillings linearly independent?
    fill_matrix = np.column_stack([fillings[v] for v in critical])
    fill_rk = np.linalg.matrix_rank(fill_matrix, tol=1e-8)
    print(f"  Fillings rank: {fill_rk} (of {len(critical)} fillings)")
    
    # What's the rank of the v-path part of each filling?
    for v in critical:
        v_part = np.zeros_like(fillings[v])
        nonv_part = np.zeros_like(fillings[v])
        for i, p in enumerate(ap2_list):
            if v in p:
                v_part[i] = fillings[v][i]
            else:
                nonv_part[i] = fillings[v][i]
        print(f"  v={v}: ||v-part||={np.linalg.norm(v_part):.4f}, ||nonv-part||={np.linalg.norm(nonv_part):.4f}")
    
    # Critical test: the ESSENTIAL v-paths
    # σ_v = σ_v^(v-paths) + σ_v^(non-v-paths)
    # ∂₂(σ_v) = z_v, where z_v is supported on non-v edges
    # But ∂₂(σ_v^(non-v-paths)) is also supported on non-v edges
    # So ∂₂(σ_v^(v-paths)) must compensate
    
    sigma3_examples.append(bits)


# ============================================================  
# Test: for n=6, is the H₁ generators rank always ≤ 3?
# ============================================================
print(f"\n{'='*70}")
print("n=6: H₁ GENERATOR RANK ANALYSIS")
print("=" * 70)

n = 6
m = n*(n-1)//2
total = 1 << m
t0 = time.time()

rank_dist = Counter()

for bits in range(total):
    A = build_adj(n, bits)
    
    # Quick β₁ check
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d1 = dim_om(om1)
    d2 = dim_om(om2)
    
    bd1 = build_full_boundary_matrix(ap1, ap0)
    rk1 = np.linalg.matrix_rank(bd1 @ om1, tol=1e-8)
    b1 = np.linalg.matrix_rank(np.linalg.lstsq(om1, build_full_boundary_matrix(ap2, ap1) @ om2, rcond=None)[0], tol=1e-8) if d2 > 0 else 0
    beta1_T = (d1 - rk1) - b1
    
    if beta1_T != 0:
        continue
    
    # Count critical vertices
    ncrit = 0
    for v in range(n):
        others = [i for i in range(n) if i != v]
        A_sub = [[A[others[i]][others[j]] for j in range(n-1)] for i in range(n-1)]
        b1v, _, _, _, _ = compute_beta1(A_sub, n-1)
        if b1v > 0:
            ncrit += 1
    
    rank_dist[ncrit] += 1
    
    if (bits+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=6: {total} tournaments in {elapsed:.0f}s (β₁=0 only)")
print(f"  #{'{critical v}'} distribution: {dict(sorted(rank_dist.items()))}")
print(f"  Max critical vertices: {max(rank_dist.keys()) if rank_dist else 0}")

print("\nDone.")
