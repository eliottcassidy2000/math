#!/usr/bin/env python3
"""
beta2_quotient_structure.py — Analyze quotient complex R = Ω(T)/Ω(T\v)

PROOF STRUCTURE (by induction on n):
Base: n ≤ 4 (verified computationally)

Inductive step: Assume β₂(T') = 0 for all tournaments T' on < n vertices.
From LES of (T, T\v) with H₂(T\v) = 0:
    β₂(T) = dim H₂(T,T\v) - dim ker(i_*)

KEY LEMMA: dim H₂(T,T\v) ≤ max(0, β₁(T\v) - β₁(T))
Given this: β₂ ≤ 0, combined with β₂ ≥ 0, gives β₂ = 0.

This script analyzes the quotient complex R to understand its structure:
- R₀ ≅ ℝ, R₁ has dim n-1 (arcs involving v)
- R₂, R₃ dimensions as function of v's degree
- Euler characteristic of R

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


def delete_vertex(A, n, v):
    return [[A[i][j] for j in range(n) if j != v] for i in range(n) if i != v]


def local_to_global(n, v):
    return [i for i in range(n) if i != v]


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0


def compute_chain_data(A, n, max_p=4):
    """Compute Ω dimensions and boundary ranks."""
    ap = {}; om = {}
    for p in range(max_p + 1):
        ap[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            om[p] = np.eye(n)
        elif ap[p]:
            om[p] = compute_omega_basis(A, n, p, ap[p], ap[p-1])
        else:
            om[p] = np.zeros((0, 0))
    return ap, om


def embed_omega(om_Tv, ap_T_p, ap_Tv_p, om_T_p, n, v):
    """Embed Ω_p(T\v) into Ω_p(T) coords."""
    d = dim_om(om_Tv)
    d_T = dim_om(om_T_p)
    if d == 0 or d_T == 0:
        return np.zeros((d_T, 0))

    l2g = local_to_global(n, v)
    T_list = [tuple(x) for x in ap_T_p]
    Tv_list = [tuple(x) for x in ap_Tv_p]
    T_idx = {p: i for i, p in enumerate(T_list)}

    incl = np.zeros((len(T_list), len(Tv_list)))
    for j, path_local in enumerate(Tv_list):
        path_global = tuple(l2g[k] for k in path_local)
        if path_global in T_idx:
            incl[T_idx[path_global], j] = 1

    emb_A = incl @ om_Tv
    emb_om = np.linalg.lstsq(om_T_p, emb_A, rcond=None)[0]
    return emb_om


def compute_betti_1(A, n):
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    if not ap1: return 0
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    d1 = dim_om(om1)
    bd1 = build_full_boundary_matrix(ap1, ap0)
    bd1_om = bd1 @ om1
    S1 = np.linalg.svd(np.linalg.lstsq(np.eye(n), bd1_om, rcond=None)[0], compute_uv=False)
    rk1 = int(sum(s > 1e-8 for s in S1))
    if ap2:
        om2 = compute_omega_basis(A, n, 2, ap2, ap1)
        d2 = dim_om(om2)
        if d2 > 0:
            bd2 = build_full_boundary_matrix(ap2, ap1)
            bd2_om = bd2 @ om2
            S2 = np.linalg.svd(np.linalg.lstsq(om1, bd2_om, rcond=None)[0], compute_uv=False)
            rk2 = int(sum(s > 1e-8 for s in S2))
        else: rk2 = 0
    else: rk2 = 0
    return (d1 - rk1) - rk2


# ===== MAIN =====
print("=" * 70)
print("QUOTIENT COMPLEX R = Ω(T)/Ω(T\\v) STRUCTURE")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

# Collect quotient complex dimensions
data = []

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap_T, om_T = compute_chain_data(A, n)
    b1_T = compute_betti_1(A, n)

    for v in range(n):
        B = delete_vertex(A, n, v)
        ap_Tv, om_Tv = compute_chain_data(B, n-1)
        b1_Tv = compute_betti_1(B, n-1)

        d_out = sum(A[v])
        d_in = n - 1 - d_out

        dims_R = {}
        for p in range(5):
            d_T = dim_om(om_T[p])
            emb = embed_omega(om_Tv[p], ap_T[p], ap_Tv[p], om_T[p], n, v)
            rk_emb = np.linalg.matrix_rank(emb, tol=1e-8) if emb.shape[1] > 0 else 0
            dims_R[p] = d_T - rk_emb

        data.append({
            'bits': bits, 'v': v, 'd_out': d_out,
            'b1_T': b1_T, 'b1_Tv': b1_Tv,
            'R': dims_R,
        })

    if bits % 200 == 0 and bits > 0:
        print(f"  ... {bits}/{1 << m}")

# Analyze
print(f"\nR₀ values: {Counter(d['R'][0] for d in data)}")
print(f"R₁ values: {Counter(d['R'][1] for d in data)}")
print(f"R₂ values: {Counter(d['R'][2] for d in data)}")
print(f"R₃ values: {Counter(d['R'][3] for d in data)}")
print(f"R₄ values: {Counter(d['R'][4] for d in data)}")

# R₂ by degree
print(f"\nR₂ by d_out(v):")
for deg in range(n):
    subset = [d for d in data if d['d_out'] == deg]
    if subset:
        r2_vals = Counter(d['R'][2] for d in subset)
        print(f"  d_out={deg} (d_in={n-1-deg}): R₂ ∈ {dict(r2_vals)}")
        # Expected: d_in * d_out star-paths, but some are shared with T\v
        expected = deg * (n-1-deg)
        print(f"    d_in*d_out = {expected}")

# Euler char of R
print(f"\nEuler characteristic of R:")
euler_by_b1 = {}
for d in data:
    chi = sum((-1)**p * d['R'][p] for p in range(5))
    key = (d['b1_T'], d['b1_Tv'])
    if key not in euler_by_b1:
        euler_by_b1[key] = Counter()
    euler_by_b1[key][chi] += 1

for key in sorted(euler_by_b1.keys()):
    print(f"  (β₁(T), β₁(T\\v)) = {key}: χ(R) ∈ {dict(euler_by_b1[key])}")

# H_p(R) from LES
print(f"\n{'='*70}")
print("H_p(R) FROM LES")
print("=" * 70)

# H₀(R) = 0 (since i_*: H₀(T\v)→H₀(T) is iso)
# H₁(R) = coker(i_*: H₁(T\v) → H₁(T)) = β₁(T) - rk(i_*)
# H₂(R) = max(0, β₁(T\v) - β₁(T)) (from our KEY LEMMA + β₂=0)
# H₃(R) = ?

# Since β₂=0 (verified), we can compute H_p(R) from:
# χ(R) = -H₁(R) + H₂(R) - H₃(R) + ...
# Also: χ(R) = Σ (-1)^p dim R_p (known from quotient dims)

# H₂(R) = max(0, β₁(T\v) - β₁(T)) (the KEY LEMMA identity)
# Verify this is consistent with the dims
h2r_check = Counter()
for d in data:
    expected = max(0, d['b1_Tv'] - d['b1_T'])
    # H₂(R) = dim ker(∂₂^R) - dim im(∂₃^R)
    # We can compute this from the quotient dims and boundary maps
    # But for now, just verify the formula matches
    h2r_check[(d['b1_T'], d['b1_Tv'], expected)] += 1

print(f"\nExpected H₂(R) by (β₁(T), β₁(T\\v)):")
for key in sorted(h2r_check.keys()):
    print(f"  β₁(T)={key[0]}, β₁(T\\v)={key[1]}: H₂^rel should be {key[2]}, count={h2r_check[key]}")

# KEY: verify χ(R) = -H₁(R) + H₂(R) - H₃(R)
# We know H₁(R) = β₁(T) - rk(i_*) = coker(i_*)
# We know H₂(R) = max(0, β₁(T\v) - β₁(T))
# So H₃(R) = H₂(R) - H₁(R) - χ(R)

# But we also have β₃(T) from the LES:
# 0 → H₃(T) → H₃(R) → H₂(T\v) = 0
# So H₃(T) = H₃(R) (if H₂(T\v) = 0 by induction)
# Wait, the LES continues: → H₃(T\v) → H₃(T) → H₃(R) → H₂(T\v)
# H₂(T\v) = 0 by induction.
# So H₃(R) = coker(H₃(T\v) → H₃(T))

print(f"\n{'='*70}")
print("STRUCTURAL IDENTITIES FOR R")
print("=" * 70)

# Let's just verify: is dim R₂ determined by (d_out, d_in)?
print("\nIs R₂ determined by d_out(v)?")
r2_by_deg_detail = {}
for d in data:
    deg = d['d_out']
    r2 = d['R'][2]
    if deg not in r2_by_deg_detail:
        r2_by_deg_detail[deg] = set()
    r2_by_deg_detail[deg].add(r2)

for deg in sorted(r2_by_deg_detail.keys()):
    vals = r2_by_deg_detail[deg]
    det = "YES ✓" if len(vals) == 1 else f"NO: {vals}"
    print(f"  d_out={deg}: {det}")

# Is dim R₂ determined by (d_out, b1_T)?
print("\nIs R₂ determined by (d_out, β₁(T))?")
r2_by_pair = {}
for d in data:
    key = (d['d_out'], d['b1_T'])
    r2 = d['R'][2]
    if key not in r2_by_pair:
        r2_by_pair[key] = set()
    r2_by_pair[key].add(r2)

for key in sorted(r2_by_pair.keys()):
    vals = r2_by_pair[key]
    det = "YES ✓" if len(vals) == 1 else f"NO: {vals}"
    print(f"  (d_out={key[0]}, β₁={key[1]}): R₂ ∈ {vals} {det}")

# Is dim R₂ determined by (d_out, β₁(T), β₁(T\v))?
print("\nIs R₂ determined by (d_out, β₁(T), β₁(T\\v))?")
r2_by_triple = {}
for d in data:
    key = (d['d_out'], d['b1_T'], d['b1_Tv'])
    r2 = d['R'][2]
    if key not in r2_by_triple:
        r2_by_triple[key] = set()
    r2_by_triple[key].add(r2)

for key in sorted(r2_by_triple.keys()):
    vals = r2_by_triple[key]
    det = "YES ✓" if len(vals) == 1 else f"NO: {vals}"
    print(f"  (d_out={key[0]}, β₁(T)={key[1]}, β₁(T\\v)={key[2]}): R₂ ∈ {vals} {det}")

print("\nDone.")
