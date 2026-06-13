#!/usr/bin/env python3
"""
beta2_z2_b2_structure.py — Direct analysis of why Z₂ = B₂

Study the explicit structure of 2-cycles (Z₂) and 2-boundaries (B₂)
in tournaments. For β₂=0, every 2-cycle must be a 2-boundary.

Key questions:
1. What do 2-cycles look like? Are they always "balanced flows"?
2. What produces 2-boundaries? Which 3-chains map to which 2-chains?
3. Is there an explicit pairing/involution that proves Z₂ = B₂?

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
print("Z₂ = B₂ STRUCTURE ANALYSIS")
print("=" * 70)

# For each tournament, compute:
# - Explicit Z₂ basis (2-cycles in Ω₂)
# - Explicit B₂ basis (images of ∂₃ from Ω₃)
# - The "matching" between them

# Focus on tournaments with interesting Z₂

# First: study Ω₃ structure
print("\n--- Ω₃ (DT 4-paths) structure ---")
print("A 4-path (a,b,c,d) is allowed if a→b, b→c, c→d.")
print("It's in Ω₃ if its boundary is in Ω₂.")
print("∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)")
print("Faces: three 2-paths, each must be in A₂ for ∂₃ ∈ Ω₂.")
print()

# At n=5, enumerate DT (doubly-transitive) 4-paths vs Ω₃
dt_vs_omega3 = Counter()

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    ap3 = enumerate_allowed_paths(A, n, 3)
    ap2 = enumerate_allowed_paths(A, n, 2)
    if not ap3: continue
    om3 = compute_omega_basis(A, n, 3, ap3, ap2)
    d3 = om3.shape[1] if om3.ndim == 2 else 0

    # Count DT paths: (a,b,c,d) with a→c AND b→d
    n_dt = sum(1 for a,b,c,d in ap3 if A[a][c] and A[b][d])
    dt_vs_omega3[(n_dt, d3)] += 1

print("(|DT|, dim Ω₃): count")
for key, cnt in sorted(dt_vs_omega3.items()):
    eq = "=" if key[0] == key[1] else "≠"
    print(f"  DT={key[0]}, Ω₃={key[1]} ({eq}): {cnt}")

# Study the surplus = dim(Ω₃) - dim(Z₂) = β₃ + rk(∂₃)
# Actually surplus = dim Ω₃ - dim Z₂ where Z₂ = ker(∂₂)
# No: rk(∂₃) = B₂ = Z₂ (since β₂=0), and dim Ω₃ - rk(∂₃) = ker(∂₃) - β₃... 
# Let me just compute things.

print(f"\n{'='*70}")
print("DETAILED CYCLE-BOUNDARY MATCHING")
print("=" * 70)

# Take specific tournaments and look at the matching
for bits in [0, 2, 4, 8, 15]:
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    scores = sorted(sum(A[i]) for i in range(n))

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
    U2, S2, Vt2 = np.linalg.svd(coords2, full_matrices=True)
    rk_bd2 = int(sum(s > 1e-8 for s in S2))
    z2_dim = d2 - rk_bd2

    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))
    d3 = om3.shape[1] if om3.ndim == 2 else 0

    if d3 > 0:
        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd3_om = bd3 @ om3
        coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
        rk_bd3 = np.linalg.matrix_rank(coords3, tol=1e-8)
    else:
        rk_bd3 = 0

    print(f"\nT#{bits} scores={scores}: Ω₂={d2}, Z₂={z2_dim}, Ω₃={d3}, B₂={rk_bd3}, β₂={z2_dim-rk_bd3}")

    # Explicit 2-cycles
    if z2_dim > 0:
        z2_basis = Vt2[rk_bd2:].T  # in Ω₂ coords
        z2_A = om2 @ z2_basis  # in A₂ coords

        print(f"  2-cycles (Z₂ basis):")
        for col in range(z2_dim):
            c = z2_A[:, col]
            terms = []
            for i in range(len(c)):
                if abs(c[i]) > 1e-8:
                    a,b,cc = ap2[i]
                    tt = "TT" if A[a][cc] else "NT"
                    terms.append(f"{c[i]:+.3f}*({a}{b}{cc})[{tt}]")
            print(f"    z{col}: {' '.join(terms)}")

    # Explicit 3-boundaries
    if d3 > 0 and rk_bd3 > 0:
        # Get B₂ basis (images of Ω₃ under ∂₃ that are independent)
        b2_om = coords3  # Ω₂ coords of ∂₃(Ω₃)
        b2_A = om2 @ b2_om  # A₂ coords

        print(f"  ∂₃(Ω₃) images (first {min(rk_bd3, 3)}):")
        for col in range(min(d3, 3)):
            c = b2_A[:, col]
            terms = []
            for i in range(len(c)):
                if abs(c[i]) > 1e-8:
                    a,b,cc = ap2[i]
                    tt = "TT" if A[a][cc] else "NT"
                    terms.append(f"{c[i]:+.3f}*({a}{b}{cc})[{tt}]")
            print(f"    ∂₃(e{col}): {' '.join(terms)}")

        # Show the Ω₃ elements
        print(f"  Ω₃ basis:")
        for col in range(min(d3, 3)):
            c = om3[:, col]
            terms = []
            for i in range(len(c)):
                if abs(c[i]) > 1e-8:
                    a,b,cc,d = ap3[i]
                    dt = "DT" if A[a][cc] and A[b][d] else "nDT"
                    terms.append(f"{c[i]:+.3f}*({a}{b}{cc}{d})[{dt}]")
            print(f"    e{col}: {' '.join(terms)}")

# Key question: is there a FORMULA for Z₂ in terms of the tournament?
print(f"\n{'='*70}")
print("Z₂ DIMENSION FORMULA")
print("=" * 70)

# Z₂ = ker(∂₂|_{Ω₂})
# ∂₂ maps Ω₂ → Ω₁ = A₁
# rk(∂₂|_{Ω₂}) = dim Ω₂ - dim Z₂
# 
# Since Ω₂ = TT paths (mostly), rk(∂₂|_{Ω₂}) relates to how
# independent the TT boundary images are.
#
# ∂₂(a,b,c) for TT: a→b, b→c, a→c
# = (b,c) - (a,c) + (a,b)
# These are "triangle flows" — adding flow on the triangle (a,b,c).
#
# The rank of the ∂₂ map = number of independent triangle flows.
# This is related to the CYCLE SPACE of the underlying undirected graph.
# For K_n (complete graph): cycle space has dim C(n,2) - (n-1) = C(n,2) - n + 1.
# At n=5: 10 - 4 = 6. And we see rk(∂₂) ∈ {5, 6}.
# Z₂ = Ω₂ - rk(∂₂).

# At n=5: Ω₂ ∈ {8,9,10}
# rk(∂₂) ∈ {5, 6}
# Z₂ = Ω₂ - rk(∂₂) ∈ {3, 4, 5}

# What determines rk(∂₂)?
# The images of TT paths are triangle flows.
# The rank depends on which triangles are "transitive" (TT).
# 
# Observation: if ALL C(n,3) triangles are TT (only for transitive tournament),
# then rk(∂₂) = C(n,2) - n + 1 (full cycle space rank).
# At n=5: rk = 6 for transitive.
# 
# When some triangles are 3-cycles (NT paths): those are NOT in Ω₂,
# so the rank can drop.

# Check: is rk(∂₂) = C(n,2) - n + 1 when Ω₂ = |TT|?
# At n=5: C(5,2) - 4 = 6
rk_analysis = Counter()
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

    om1 = compute_omega_basis(A, n, 1, ap1, enumerate_allowed_paths(A, n, 0))
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk = np.linalg.matrix_rank(coords, tol=1e-8)

    n_tt = sum(1 for a,b,c in ap2 if A[a][c])
    t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if not (A[i][j] and A[j][k] and A[i][k]) and
                not (A[j][i] and A[i][k] and A[j][k]) and
                not (A[k][i] and A[i][j] and A[k][j]) and
                not (A[k][j] and A[j][i] and A[k][i]))
    # Count 3-cycles
    c3 = 0
    for i in range(n):
        for j in range(n):
            if i == j: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[i][j] and A[j][k] and A[k][i]:
                    c3 += 1
    c3 //= 3  # each 3-cycle counted 3 times

    rk_analysis[(d2, rk, c3)] += 1

print("\n(dim Ω₂, rk ∂₂, c₃): count")
for key, cnt in sorted(rk_analysis.items()):
    z2 = key[0] - key[1]
    print(f"  Ω₂={key[0]}, rk∂₂={key[1]}, c₃={key[2]}, Z₂={z2}: {cnt}")

print("\nDone.")
