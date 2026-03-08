#!/usr/bin/env python3
"""
beta2_connecting_map.py — Analyze Ω₂ structure to find tournament-specific β₂=0 argument

The LES reformulation β₂ = H₂^rel - dim(ker i_*) is circular.
We need a DIRECT argument. Key questions:
1. Is Ω₂ spanned entirely by TT (transitive triple) paths?
2. What do 2-cycles look like in tournaments?
3. Why does tournament completeness kill β₂?

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
print("Ω₂ STRUCTURAL ANALYSIS FOR TOURNAMENTS")
print("=" * 70)

# Q1: Is Ω₂ spanned entirely by TT paths?
print("\nQ1: Is Ω₂ = TT-span in tournaments?")
nt_in_omega = 0
total_checked = 0
tt_dim_vs_omega = Counter()

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap2 = enumerate_allowed_paths(A, n, 2)
    ap1 = enumerate_allowed_paths(A, n, 1)
    if not ap2: continue
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    d2 = om2.shape[1] if om2.ndim == 2 else 0
    if d2 == 0: continue
    total_checked += 1

    # Count TT paths
    tt_indices = [i for i, (a,b,c) in enumerate(ap2) if A[a][c]]
    nt_indices = [i for i, (a,b,c) in enumerate(ap2) if not A[a][c]]
    n_tt = len(tt_indices)

    # Check if Ω₂ basis has any NT components
    has_nt = False
    for col in range(d2):
        v = om2[:, col]
        for i in nt_indices:
            if abs(v[i]) > 1e-10:
                has_nt = True
                break
        if has_nt: break

    if has_nt:
        nt_in_omega += 1
    tt_dim_vs_omega[(n_tt, d2)] += 1

print(f"  Tournaments with NT in Ω₂: {nt_in_omega}/{total_checked}")

if nt_in_omega == 0:
    print("  *** Ω₂ = TT-span ALWAYS in n=5 tournaments ***")
    print("  This means dim(Ω₂) = |TT| always!")

print("\n  (|TT|, dim Ω₂): count")
for key, cnt in sorted(tt_dim_vs_omega.items()):
    eq = "=" if key[0] == key[1] else "≠"
    print(f"    TT={key[0]}, Ω₂={key[1]} ({eq}): {cnt}")

# Q2: Explicit 2-cycles
print(f"\n{'='*70}")
print("Q2: What do 2-cycles look like?")
print("=" * 70)

cycle_info = Counter()
example_count = 0

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap2 = enumerate_allowed_paths(A, n, 2)
    ap1 = enumerate_allowed_paths(A, n, 1)
    if not ap2: continue
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    d2 = om2.shape[1] if om2.ndim == 2 else 0
    if d2 == 0: continue

    # Boundary map
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    om1 = compute_omega_basis(A, n, 1, ap1, enumerate_allowed_paths(A, n, 0))
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    U, S, Vt = np.linalg.svd(coords2, full_matrices=True)
    rk = int(sum(s > 1e-8 for s in S))
    ker_dim = d2 - rk

    if ker_dim > 0:
        cycle_info[ker_dim] += 1
        if example_count < 5:
            example_count += 1
            ker_vecs = Vt[rk:].T
            scores = sorted(sum(A[i]) for i in range(n))
            print(f"\n  T#{bits} (scores={scores}), Ω₂={d2}, Z₂={ker_dim}:")

            for col in range(min(ker_dim, 2)):
                c_om = ker_vecs[:, col]
                c_A = om2 @ c_om
                terms = [(tuple(ap2[i]), round(c_A[i],4))
                         for i in range(len(c_A)) if abs(c_A[i]) > 1e-8]
                print(f"    z{col} = ", end="")
                for path, coeff in terms:
                    a, b, c = path
                    print(f"  {coeff:+.3f}*({a}{b}{c})", end="")
                print()

# Q3: Dimension analysis — Z₂ vs B₂
print(f"\n{'='*70}")
print("Q3: Z₂ = B₂ verification and dimension analysis")
print("=" * 70)

dim_data = []
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap2 = enumerate_allowed_paths(A, n, 2)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap2: continue

    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    d2 = om2.shape[1] if om2.ndim == 2 else 0
    if d2 == 0: continue

    om1 = compute_omega_basis(A, n, 1, ap1, enumerate_allowed_paths(A, n, 0))
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk_bd2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    ker_bd2 = d2 - rk_bd2

    rk_bd3 = 0
    d3 = 0
    if ap3:
        om3 = compute_omega_basis(A, n, 3, ap3, ap2)
        d3 = om3.shape[1] if om3.ndim == 2 else 0
        if d3 > 0:
            bd3 = build_full_boundary_matrix(ap3, ap2)
            bd3_om = bd3 @ om3
            coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
            rk_bd3 = np.linalg.matrix_rank(coords3, tol=1e-8)

    dim_data.append({
        'bits': bits,
        'd2': d2, 'd3': d3,
        'rk_bd2': rk_bd2, 'rk_bd3': rk_bd3,
        'Z2': ker_bd2, 'B2': rk_bd3,
        'b2': ker_bd2 - rk_bd3,
        'tt': sum(1 for a,b,c in ap2 if A[a][c]),
        'nt': sum(1 for a,b,c in ap2 if not A[a][c]),
    })

print(f"\nDimension statistics ({len(dim_data)} tournaments):")
z2_dist = Counter(d['Z2'] for d in dim_data)
b2_dist = Counter(d['B2'] for d in dim_data)
beta2_dist = Counter(d['b2'] for d in dim_data)

print(f"  Z₂ distribution: {dict(sorted(z2_dist.items()))}")
print(f"  B₂ distribution: {dict(sorted(b2_dist.items()))}")
print(f"  β₂ distribution: {dict(sorted(beta2_dist.items()))}")

# Check Z₂ = B₂ always
all_zero = all(d['b2'] == 0 for d in dim_data)
print(f"  β₂ = 0 always: {all_zero}")

# Relation between Z₂ and |TT|
print(f"\n  Z₂ as function of |TT|:")
z2_by_tt = {}
for d in dim_data:
    tt = d['tt']
    z2 = d['Z2']
    if tt not in z2_by_tt:
        z2_by_tt[tt] = Counter()
    z2_by_tt[tt][z2] += 1

for tt in sorted(z2_by_tt):
    print(f"    TT={tt}: Z₂ = {dict(sorted(z2_by_tt[tt].items()))}")

# Check if Z₂ = |TT| - rk(∂₂)
# rk(∂₂) should relate to Ω₁ structure
print(f"\n  rk(∂₂) as function of Ω₂ = TT:")
rk_by_d2 = {}
for d in dim_data:
    d2 = d['d2']
    rk = d['rk_bd2']
    if d2 not in rk_by_d2:
        rk_by_d2[d2] = Counter()
    rk_by_d2[d2][rk] += 1

for d2 in sorted(rk_by_d2):
    print(f"    Ω₂={d2}: rk(∂₂) = {dict(sorted(rk_by_d2[d2].items()))}")

# Check: is rk(∂₂) = Ω₂ - Z₂ always equal to some simple expression?
print(f"\n  Is Z₂ predictable from |TT|?")
# Z₂ = dim Ω₂ - rk(∂₂), and dim Ω₂ = |TT|
# So rk(∂₂) = |TT| - Z₂
# rk(∂₂) should equal something related to Ω₁ or edges

for d in dim_data[:10]:
    print(f"    TT={d['tt']}, Ω₂={d['d2']}, rk∂₂={d['rk_bd2']}, Z₂={d['Z2']}, "
          f"B₂={d['B2']}, β₂={d['b2']}")

print("\nDone.")
