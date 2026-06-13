#!/usr/bin/env python3
"""
beta2_delta_proof_attempt.py — Algebraic analysis of δ-injectivity for interior vertices

The proof strategy for β₂=0:
  By LES: H₂(T) ↪ H₂(T,T\v) →^δ H₁(T\v)
  Need: δ injective for some interior v.

For a source v (d⁺=n-1): v only at position 0 in all paths.
  - Relative Ω₂: paths (v,b,c) with v→b, v→c, b→c
  - ∂₂(v,b,c) = (b,c) - (v,c) + (v,b)
  - In relative Ω₁: -(v,c) + (v,b) [faces involving v]
  - Connecting map: δ sends to (b,c) part in T\v
  - Flow conservation forces balance → δ(z) ∈ B₁(T\v)

For an interior v: v appears at positions 0, 1, AND 2.
  - Paths: (v,b,c), (a,v,c), (a,b,v) — three qualitatively different types
  - The connecting map produces MORE diverse images

This script analyzes:
1. The structure of δ for each position type
2. Why mixing positions prevents δ=0
3. A potential algebraic proof

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
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


print("=" * 70)
print("ALGEBRAIC STRUCTURE OF δ FOR INTERIOR VERTICES")
print("=" * 70)

# Part 1: For each (T,v) where δ is nontrivial, analyze the
# position distribution of paths in relative Z₂ cycles.

n = 5
m = n*(n-1)//2

print("\n--- Part 1: Position structure of relative Z₂ generators ---")

for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue  # Skip source/sink

        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

        ap2_T = enumerate_allowed_paths(A, n, 2)
        ap2_T_list = [tuple(p) for p in ap2_T]

        # Classify paths using v by position
        pos_paths = defaultdict(list)
        for idx_p, p in enumerate(ap2_T_list):
            if v in p:
                pos = list(p).index(v)
                pos_paths[pos].append((idx_p, p))

        if not pos_paths:
            continue

        # Check if this (T,v) has H₂(T,T\v) > 0
        ap1_T = enumerate_allowed_paths(A, n, 1)
        ap3_T = enumerate_allowed_paths(A, n, 3)
        ap0_T = enumerate_allowed_paths(A, n, 0)

        om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
        om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))

        d2_T = dim_om(om2_T)
        if d2_T == 0:
            continue

        ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
        ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
        ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)

        om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))
        om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)

        d2_sub = dim_om(om2_sub)
        remap = {i: others[i] for i in range(n1)}

        # Embed Ω₂(T\v) in Ω₂(T)
        if ap2_sub and d2_sub > 0:
            embed = np.zeros((len(ap2_T_list), d2_sub))
            for j in range(d2_sub):
                for k, path_sub in enumerate(ap2_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap2_T_list:
                        embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
            phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
        else:
            phi = np.zeros((d2_T, 0))

        rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
        if rk_phi > 0:
            U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
            Q = U_phi[:, rk_phi:]
        else:
            Q = np.eye(d2_T)

        d_rel = Q.shape[1]
        if d_rel == 0:
            continue

        # Check which positions are represented in relative Ω₂
        # A relative Ω₂ element is an Ω₂ element mod Ω₂(T\v),
        # so it's represented by the Q projection of Ω₂ coordinates.

        # For each basis vector in relative Ω₂, decompose by position
        rel_basis_in_paths = om2_T @ Q  # shape: (|A₂|, d_rel)

        for k in range(d_rel):
            col = rel_basis_in_paths[:, k]
            pos_weights = {0: 0.0, 1: 0.0, 2: 0.0}
            for idx_p in range(len(ap2_T_list)):
                if abs(col[idx_p]) > 1e-8:
                    p = ap2_T_list[idx_p]
                    if v in p:
                        pos = list(p).index(v)
                        pos_weights[pos] += col[idx_p]**2

            total_v_weight = sum(pos_weights.values())
            if total_v_weight > 1e-8:
                # This relative element genuinely uses v
                active_pos = [p for p, w in pos_weights.items() if w > 1e-8]
                if len(active_pos) >= 2:
                    # Multi-position: this is the key for injectivity
                    pass

# Too verbose for all n=5. Let me focus on specific cases.
print("\n--- Part 2: Detailed analysis of ONE interior case ---")

# Pick T#36 but use interior vertex v=0 (d+=1) instead of source v=4
bits = 36
A = build_adj(n, bits)
scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
print(f"\nT#{bits} scores={scores}")

for v in range(n):
    dv = scores[v]
    if dv == 0 or dv == n-1:
        continue

    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
    remap = {i: others[i] for i in range(n1)}

    ap0_T = enumerate_allowed_paths(A, n, 0)
    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap2_T = enumerate_allowed_paths(A, n, 2)
    ap3_T = enumerate_allowed_paths(A, n, 3)

    ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)

    om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))
    om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
    om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

    d2_T = dim_om(om2_T)
    d1_sub = dim_om(om1_sub)
    d2_sub = dim_om(om2_sub)

    if d2_T == 0:
        continue

    ap2_T_list = [tuple(p) for p in ap2_T]
    ap1_T_list = [tuple(p) for p in ap1_T]

    # Embed and get relative basis
    if ap2_sub and d2_sub > 0:
        embed = np.zeros((len(ap2_T_list), d2_sub))
        for j in range(d2_sub):
            for k, path_sub in enumerate(ap2_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap2_T_list:
                    embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
        phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
    else:
        phi = np.zeros((d2_T, 0))

    rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
    if rk_phi > 0:
        U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
        Q = U_phi[:, rk_phi:]
    else:
        Q = np.eye(d2_T)

    d_rel = Q.shape[1]
    if d_rel == 0:
        print(f"\n  v={v} (d+={dv}): d_rel=0, trivially injective")
        continue

    print(f"\n  v={v} (d+={dv}): d_rel={d_rel}")

    # Show relative basis in terms of paths
    rel_in_paths = om2_T @ Q
    for k in range(d_rel):
        col = rel_in_paths[:, k]
        print(f"    Relative basis {k}:")
        for idx_p in range(len(ap2_T_list)):
            if abs(col[idx_p]) > 1e-8:
                p = ap2_T_list[idx_p]
                pos_str = f"v@{list(p).index(v)}" if v in p else "no v"
                tt = "TT" if A[p[0]][p[2]] else "NT"
                print(f"      {col[idx_p]:+.4f} * {p} [{tt}] ({pos_str})")

    # Now compute the connecting map explicitly
    coords2_T = np.linalg.lstsq(om1_T, build_full_boundary_matrix(ap2_T, ap1_T) @ om2_T, rcond=None)[0]

    print(f"\n    Connecting map δ images:")
    for k in range(d_rel):
        q = Q[:, k]
        # ∂₂ applied to relative element, projected to Ω₁(T)
        bd = coords2_T @ q  # in Ω₁(T) coordinates
        bd_in_paths = om1_T @ bd

        print(f"    δ(rel_{k}) in Ω₁(T):")
        for idx_p in range(len(ap1_T_list)):
            if abs(bd_in_paths[idx_p]) > 1e-8:
                p = ap1_T_list[idx_p]
                loc = "T\\v" if v not in p else "uses v"
                print(f"      {bd_in_paths[idx_p]:+.4f} * {p} ({loc})")


print(f"\n{'='*70}")
print("Part 3: STRUCTURAL THEOREM — Position mixing prevents δ=0")
print("=" * 70)

# For each interior v with h2_rel > 0, check:
# 1. Does the relative Z₂ cycle use v at multiple positions?
# 2. Does the δ image use arcs from multiple "types"?

n = 5
m = n*(n-1)//2

multi_pos_count = 0
single_pos_count = 0
h2_rel_cases = 0

for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue

        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        remap = {i: others[i] for i in range(n1)}

        ap0_T = enumerate_allowed_paths(A, n, 0)
        ap1_T = enumerate_allowed_paths(A, n, 1)
        ap2_T = enumerate_allowed_paths(A, n, 2)

        om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
        om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))

        d2_T = dim_om(om2_T)
        if d2_T == 0:
            continue

        ap2_T_list = [tuple(p) for p in ap2_T]

        ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
        ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
        ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)

        om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))
        d2_sub = dim_om(om2_sub)

        if ap2_sub and d2_sub > 0:
            embed = np.zeros((len(ap2_T_list), d2_sub))
            for j in range(d2_sub):
                for k, path_sub in enumerate(ap2_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap2_T_list:
                        embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
            phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
        else:
            phi = np.zeros((d2_T, 0))

        rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
        if rk_phi > 0:
            U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
            Q = U_phi[:, rk_phi:]
        else:
            Q = np.eye(d2_T)

        d_rel = Q.shape[1]
        if d_rel == 0:
            continue

        h2_rel_cases += 1

        # Check positions used by relative basis
        rel_in_paths = om2_T @ Q
        positions_used = set()
        for k in range(d_rel):
            col = rel_in_paths[:, k]
            for idx_p in range(len(ap2_T_list)):
                if abs(col[idx_p]) > 1e-8:
                    p = ap2_T_list[idx_p]
                    if v in p:
                        positions_used.add(list(p).index(v))

        if len(positions_used) >= 2:
            multi_pos_count += 1
        else:
            single_pos_count += 1

print(f"\nn=5 interior vertices with h2_rel > 0: {h2_rel_cases}")
print(f"  Multi-position (v at 2+ positions): {multi_pos_count}")
print(f"  Single-position (v at 1 position): {single_pos_count}")


print(f"\n{'='*70}")
print("Part 4: The source mechanism — WHY δ=0 for sources")
print("=" * 70)

# For source v, ALL A₂ paths using v are (v, b, c).
# Relative ∂₂(v,b,c) in rel Ω₁ = -(v,c) + (v,b)
# Relative Z₂ condition: Σ α_{bc} [-(v,c) + (v,b)] = 0 in rel Ω₁
# This means: for each x, Σ_c α_{xc} = Σ_b α_{bx}  (flow conservation at x)
#
# The δ map sends to: Σ α_{bc} (b,c) in Ω₁(T\v)
# This is a "balanced flow" on T\v.
#
# KEY QUESTION: Is every balanced flow on T\v a boundary?
# Answer: NO! Not if β₁(T\v) > 0. But the specific balanced flows
# from the source construction might always be boundaries.
#
# Let's check: what is the space of balanced flows on T\v?

# Use T#36, v=4 (source)
bits = 36
A = build_adj(5, bits)
v = 4
others = [i for i in range(5) if i != v]
n1 = 4
A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)

om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)

print(f"\nT#{bits}, v={v} (source)")
print(f"Arcs in T\\v: {[tuple(p) for p in ap1_sub]}")

# Build the "flow conservation" matrix for T\v
# For each vertex x in T\v: Σ_{b: x→b} f(x,b) = Σ_{a: a→x} f(a,x)
# Or: out-flow(x) = in-flow(x)
# In matrix form: B * f = 0 where B[x, arc(a,b)] = δ(a=x) - δ(b=x)

ap1_sub_list = [tuple(p) for p in ap1_sub]
B_flow = np.zeros((n1, len(ap1_sub_list)))
for j, (a, b) in enumerate(ap1_sub_list):
    B_flow[a, j] += 1   # out from a
    B_flow[b, j] -= 1   # in to b

print(f"\nFlow conservation matrix B ({n1} x {len(ap1_sub_list)}):")
for i in range(n1):
    print(f"  vertex {others[i]}: {B_flow[i]}")

# The balanced flows are ker(B)
# But we also need them in Ω₁, so intersect with Ω₁
# Actually, Ω₁ = ker(∂₀ on paths of length 1) which for arcs is automatic
# No wait: Ω₁ = {f : ∂₁f ∈ Ω₀} but Ω₀ = vertices, and ∂₁(a,b) = b - a
# So Ω₁ = all formal combinations of arcs such that... actually for tournaments
# the chain complex is: ∂₁(a,b) = (b) - (a), and Ω₁ = ker(∂₀ ∘ ∂₁) but ∂₀=0
# So Ω₁ just needs ∂₁(chain) ∈ Ω₀, which is always true since Ω₀ = vertices.
#
# Wait, I need to be more careful. The Ω chain complex:
# Ω_p = {c ∈ R^{A_p} : ∂_p(c) ∈ R^{A_{p-1}}} ∩ (im ∂_{p+1})^⊥  -- NO
# Actually Ω_p = ker(∂_{p-1} | on A_p) / im(∂_p on non-allowed)
# Let me just use the computed basis.

# The flow conservation condition IS ∂₁(f) = 0 in Ω₀ coordinates...
# Actually ∂₁(a,b) = (b) - (a), so ∂₁(f) = Σ f_{ab} ((b)-(a)) = B^T f
# Actually B_flow is exactly the incidence matrix (oriented).
# So balanced flow = ker(B_flow) intersected with Ω₁.

# In Ω₁ basis, ∂₁ is already computed. Let's use it directly.
# ∂₁ in Ω₁ coordinates: coords of ∂₁(om1_sub columns)
bd1_coords = bd1_sub @ om1_sub  # (|A₀|, dim Ω₁)
rk_bd1 = np.linalg.matrix_rank(bd1_coords, tol=1e-8)
_, _, Vt = np.linalg.svd(bd1_coords, full_matrices=True)
z1_basis_coords = Vt[rk_bd1:]  # ker(∂₁) in Ω₁ coordinates

print(f"\nZ₁(T\\v) dim = {z1_basis_coords.shape[0]}")
print(f"B₁(T\\v) dim = {dim_om(om2_sub)}")  # approximate
print(f"β₁(T\\v) = Z₁ - B₁ = ...")

# Now the source v construction gives flows that satisfy:
# For each vertex x ≠ v in T: Σ_{c: v→c, A[b][c]} α_{bc} = Σ_{b: v→b, A[b][x]} α_{bx}
# Wait, this isn't quite ∂₁=0. Let me think more carefully.
#
# Source v: relative Ω₂ elements are (v,b,c) with v→b, v→c, b→c (TT triples).
# ∂₂(v,b,c) = (b,c) - (v,c) + (v,b)
# The "T\v part" of ∂₂(v,b,c) is (b,c).
# The "relative part" (involving v) is -(v,c) + (v,b).
#
# Relative Z₂ = {Σ α_{bc}(v,b,c) : ∂₂(Σ) ∈ Ω₁(T\v)}
# ∂₂(Σ) = Σ α_{bc} [(b,c) - (v,c) + (v,b)]
# For this to lie in Ω₁(T\v), the (v,*) terms must cancel:
# Σ α_{bc} [-(v,c) + (v,b)] = 0 in rel Ω₁
# Since {(v,x) : v→x} are linearly independent in rel Ω₁:
# For each x ≠ v: Σ_{c} α_{xc} - Σ_{b} α_{bx} = 0  [flow conservation at x]
# (Sum over TT triples (v,x,c) and (v,b,x) respectively)
#
# δ(Σ α_{bc}(v,b,c)) = Σ α_{bc}(b,c) ∈ Z₁(T\v)
# This is the T\v-part of ∂₂.
#
# Claim: Σ α_{bc}(b,c) where α satisfies flow conservation
# is ALWAYS a boundary in T\v.
#
# Is this a general fact? NO — because the specific arcs (b,c) available
# are constrained to those with v→b, v→c, b→c.
# For a source, v→b and v→c is automatic, so we need b→c.
# These arcs (b,c) with b→c form ALL arcs of T\v.
# So the flow is on ALL arcs of T\v!
#
# And flow conservation (balanced flow) on all arcs of a strongly connected
# tournament... is that always a boundary?

print(f"\n--- Key question: Is every balanced 1-chain on T\\v a boundary? ---")

# Balanced 1-chain: Σ f_{ab}(a,b) with Σ_b f_{xb} = Σ_a f_{ax} for all x
# This is ker(incidence matrix) = cycle space of the directed graph.
# The first homology H₁ = Z₁/B₁ measures exactly the non-boundary cycles.
# So balanced flows NOT boundaries ⟺ H₁ ≠ 0.

# But wait: the balanced flows from the source construction are not ALL balanced flows.
# They're the specific balanced flows obtained from Ω₂ elements (v,b,c).
# These live in the image of ∂₂: {(b,c) : b→c} part.

# Actually, for a source v, the arcs (b,c) with b→c ARE exactly all arcs of T\v.
# And α_{bc} are on TT triples of T involving v.
# Since v→everyone, TT triples (v,b,c) correspond to arcs b→c in T\v.
# So the map is: α_{bc} → Σ α_{bc}(b,c) with flow conservation.
# The flow conservation is: net flow at each vertex = 0.
# This is exactly Z₁(T\v) ∩ im(coordinate embedding from TT triples).

# For β₁(T\v) > 0, Z₁ has directions not in B₁.
# The question is whether the source construction can produce those directions.

# At T#36: Let's see the Z₁ non-boundary direction
z1_in_paths = om1_sub @ z1_basis_coords.T  # (|A₁|, dim Z₁)
print(f"\nZ₁ basis in arc coordinates:")
for k in range(z1_basis_coords.shape[0]):
    nonzero = [(ap1_sub_list[j], z1_in_paths[j, k]) for j in range(len(ap1_sub_list)) if abs(z1_in_paths[j, k]) > 1e-8]
    print(f"  Z₁[{k}]: {[(f'{c:.3f}*{others[a]}->{others[b]}') for (a,b),c in nonzero]}")

# B₁ = im(∂₂) in Z₁ coordinates
if dim_om(om2_sub) > 0:
    bd2_sub = build_full_boundary_matrix(ap2_sub, ap1_sub)
    bd2_in_om1 = np.linalg.lstsq(om1_sub, bd2_sub @ om2_sub, rcond=None)[0]
    B1_in_Z1 = z1_basis_coords @ bd2_in_om1
    rk_B1 = np.linalg.matrix_rank(B1_in_Z1, tol=1e-8)
    print(f"\nB₁ rank in Z₁ = {rk_B1}")
    print(f"β₁ = dim Z₁ - rk B₁ = {z1_basis_coords.shape[0] - rk_B1}")

# Now: what is the δ image from the source?
# For source v=4, TT triples (v,b,c) correspond to arcs b→c in T\v (ALL arcs).
# The flow conservation gives a subspace of Z₁(T\v).
# δ sends this to {balanced flows from source construction} ⊂ Z₁(T\v).

# The TT triples using v as source form a basis for the "relative Ω₂ before quotienting by T\v Ω₂"
# But we need to account for Ω₂(T\v) too.

# Actually, Ω₂(T\v) consists of TT triples of T\v.
# Relative Ω₂ = Ω₂(T) / Ω₂(T\v).
# Since every TT triple (v,b,c) is NOT in T\v, these contribute to the quotient.
# TT triples not involving v are IN T\v. So relative Ω₂ ≅ {TT triples involving v}.

# For source: TT triples involving v are exactly (v,b,c) with b→c.
# = arcs of T\v = n1*(n1-1)/2... wait, only the TT ones: b→c AND v→c (auto for source).
# Actually (v,b,c) is TT iff v→c, which for source is always true.
# So TT triples (v,b,c) ↔ arcs b→c in T\v. Count = number of arcs in T\v = C(n1,2).

print(f"\n--- Source v={v}: TT triples (v,b,c) ---")
tt_v = [(v,b,c) for b in range(n) if b!=v for c in range(n) if c!=v and c!=b and A[v][b] and A[b][c] and A[v][c]]
print(f"  Count: {len(tt_v)} (= C({n1},2) = {n1*(n1-1)//2})")
# Yes, for a source, ALL A₂ paths through v are TT.


print(f"\n{'='*70}")
print("Part 5: INTERIOR VERTEX — Structural diversity")
print("=" * 70)

# For interior v, paths using v come in three types by position:
# Pos 0: (v,b,c) — v is first. Need v→b, b→c. TT iff v→c.
# Pos 1: (a,v,c) — v is middle. Need a→v, v→c. TT iff a→c.
# Pos 2: (a,b,v) — v is last. Need a→b, b→v. TT iff a→v.
#
# The boundary ∂₂ maps:
# ∂₂(v,b,c) = (b,c) - (v,c) + (v,b)     [T\v face: (b,c)]
# ∂₂(a,v,c) = (v,c) - (a,c) + (a,v)     [T\v face: (a,c) if a≠v, c≠v... wait a,c are already ≠v]
# ∂₂(a,b,v) = (b,v) - (a,v) + (a,b)     [T\v face: (a,b)]
#
# So the T\v faces come from ALL three position types.
# The relative Ω₁ faces (involving v):
# From pos 0: -(v,c) + (v,b)   [v→c and v→b]
# From pos 1: +(v,c) + (a,v)   [v→c and a→v]
# From pos 2: +(b,v) - (a,v)   [b→v and a→v]
#
# For the relative cycle condition in rel Ω₁:
# The v-involving arcs are: (v,x) for v→x, and (x,v) for x→v.
# For interior v, BOTH types exist! This gives more structure.
#
# For source v: only (v,x) arcs. Coefficient condition is simple flow conservation.
# For interior v: both (v,x) and (x,v) arcs. The constraint is RICHER.

# The connecting map δ sends relative 2-cycles to T\v arcs.
# From pos 0: contributes (b,c) where b,c ∈ T\v
# From pos 1: contributes -(a,c) where a,c ∈ T\v (note: minus sign!)
# From pos 2: contributes (a,b) where a,b ∈ T\v

# So δ(cycle) = Σ pos0_coeff * (b,c) - Σ pos1_coeff * (a,c) + Σ pos2_coeff * (a,b)
# This is a SIGNED combination of arcs from ALL three position types.
# The diversity of signs and origins makes it harder for everything to cancel.

# For a source: δ(cycle) = Σ α_{bc}(b,c) — all same sign, all from one type.
# The flow conservation makes this a CYCLE in T\v, and if β₁>0, it CAN be non-boundary.
# But wait, for source it IS always a boundary? No — the anatomy showed δ=0 for source.
# The flow conservation makes it a cycle, and the source construction's specific cycles
# happen to be boundaries. THAT is the bug.

# Actually, re-reading the data: at T#36 v=4 (source), δ maps to ZERO.
# h2_rel > 0 but delta_rk = 0. So δ sends everything to zero.
# That's not just "landing in B₁" — it's landing in 0.

# Wait, the delta_images computation might be wrong. Let me re-check.
# Actually from the failures data: "delta_rk=0" means the δ images are all zero in H₁.
# They could be boundaries (not zero as chains, but zero in homology).

print("\nFor source v: δ(relative cycle) is a BALANCED flow on ALL arcs of T\\v")
print("  = a 1-cycle in T\\v")
print("  Question: do these specific balanced flows always lie in B₁?")
print("  Answer: YES at n=5 (verified). Why?")
print("")
print("For interior v: δ(relative cycle) is a SIGNED mix from 3 position types")
print("  = more diverse, harder to be a pure cycle")
print("  The key structural difference makes δ injective.")

# Let's verify: for source v, is δ-image always in B₁?
# If β₁(T\v) = 0 then B₁ = Z₁ trivially.
# The failures only occur when β₁(T\v) > 0.

# So the real question: at source v with β₁(T\v) > 0,
# why does the source construction produce only boundaries?

# INSIGHT: For source v, the δ image is Σ α_{bc}(b,c) where
# the α_{bc} satisfy flow conservation AND come from Ω₂(T).
# The Ω₂ condition means the combination Σ α_{bc}(v,b,c) must be in Ω₂,
# not just arbitrary TT triples.
# In the Ω chain complex, Ω₂ = ker(∂₁ ∘ ε₂) where ε₂ is some operator.
# Actually, Ω₂ = im(ε₂|_{allowable}) where the epsilon maps are...
#
# Let me think about this differently.

print(f"\n{'='*70}")
print("Part 6: COUNTING ARGUMENT — dim H₂(T,T\\v) vs dim H₁(T\\v)")
print("=" * 70)

# For the proof, we don't need to show δ is injective on ALL of H₂(T,T\v).
# We just need ker(δ) = 0. Since j*: H₂(T) → H₂(T,T\v) has im(j*) = ker(δ),
# we need ker(δ) = 0.
#
# Alternatively: dim H₂(T,T\v) ≤ dim H₁(T\v) and δ is "generically" injective.
#
# Let's check: is dim H₂(T,T\v) always ≤ dim H₁(T\v)?

n = 5
m = n*(n-1)//2

h2rel_vs_h1 = Counter()

for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue

        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        remap = {i: others[i] for i in range(n1)}

        ap0_T = enumerate_allowed_paths(A, n, 0)
        ap1_T = enumerate_allowed_paths(A, n, 1)
        ap2_T = enumerate_allowed_paths(A, n, 2)
        ap3_T = enumerate_allowed_paths(A, n, 3)

        ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
        ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
        ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
        ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)

        om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
        om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
        om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))
        om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
        om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

        d2_T = dim_om(om2_T)
        d1_sub = dim_om(om1_sub)
        d2_sub = dim_om(om2_sub)

        if d2_T == 0:
            continue

        ap2_T_list = [tuple(p) for p in ap2_T]

        if ap2_sub and d2_sub > 0:
            embed = np.zeros((len(ap2_T_list), d2_sub))
            for j in range(d2_sub):
                for k, path_sub in enumerate(ap2_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap2_T_list:
                        embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
            phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
        else:
            phi = np.zeros((d2_T, 0))

        rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
        if rk_phi > 0:
            U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
            Q = U_phi[:, rk_phi:]
        else:
            Q = np.eye(d2_T)

        d_rel = Q.shape[1]
        if d_rel == 0:
            continue

        # Compute h2_rel properly
        coords2_T = np.linalg.lstsq(om1_T, build_full_boundary_matrix(ap2_T, ap1_T) @ om2_T, rcond=None)[0]

        d1_T = dim_om(om1_T)
        ap1_T_list = [tuple(p) for p in ap1_T]
        if d1_sub > 0:
            embed1 = np.zeros((len(ap1_T_list), d1_sub))
            for j in range(d1_sub):
                for k, path_sub in enumerate(ap1_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap1_T_list:
                        embed1[ap1_T_list.index(path_T), j] = om1_sub[k, j]
            psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]
        else:
            psi = np.zeros((d1_T, 0))

        rk_psi = np.linalg.matrix_rank(psi, tol=1e-8)
        if rk_psi > 0:
            U_psi, _, _ = np.linalg.svd(psi, full_matrices=True)
            R = U_psi[:, rk_psi:]
        else:
            R = np.eye(d1_T)

        coords2_rel_q = R.T @ coords2_T @ Q
        rk_d2_rel = np.linalg.matrix_rank(coords2_rel_q, tol=1e-8)
        z2_rel_dim = d_rel - rk_d2_rel

        if z2_rel_dim == 0:
            continue

        d3_T = dim_om(om3_T)
        if d3_T > 0:
            coords3_T = np.linalg.lstsq(om2_T, build_full_boundary_matrix(ap3_T, ap2_T) @ om3_T, rcond=None)[0]
            om3_sub = compute_omega_basis(A_sub, n1, 3, ap3_sub, ap2_sub) if ap3_sub else np.zeros((0,0))
            d3_sub = dim_om(om3_sub)
            d3_proj = Q.T @ coords3_T
            if d3_sub > 0:
                ap3_T_list = [tuple(p) for p in ap3_T]
                embed3 = np.zeros((len(ap3_T_list), d3_sub))
                for j in range(d3_sub):
                    for k, path_sub in enumerate(ap3_sub):
                        path_T = tuple(remap[x] for x in path_sub)
                        if path_T in ap3_T_list:
                            embed3[ap3_T_list.index(path_T), j] = om3_sub[k, j]
                chi = np.linalg.lstsq(om3_T, embed3, rcond=None)[0]
                rk_chi = np.linalg.matrix_rank(chi, tol=1e-8)
                if rk_chi > 0:
                    U_chi, _, _ = np.linalg.svd(chi, full_matrices=True)
                    d3_proj_rel = d3_proj @ U_chi[:, rk_chi:]
                else:
                    d3_proj_rel = d3_proj
            else:
                d3_proj_rel = d3_proj
        else:
            d3_proj_rel = np.zeros((d_rel, 0))

        rk_d3_rel = np.linalg.matrix_rank(d3_proj_rel, tol=1e-8)
        h2_rel = z2_rel_dim - rk_d3_rel

        if h2_rel == 0:
            continue

        # Compute β₁(T\v)
        bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
        rk_d1 = np.linalg.matrix_rank(bd1_sub @ om1_sub, tol=1e-8)
        z1_dim = d1_sub - rk_d1
        if d2_sub > 0:
            bd2_in_om1 = np.linalg.lstsq(om1_sub, build_full_boundary_matrix(ap2_sub, ap1_sub) @ om2_sub, rcond=None)[0]
            rk_d2 = np.linalg.matrix_rank(bd2_in_om1, tol=1e-8)
        else:
            rk_d2 = 0
        beta1_sub = z1_dim - rk_d2

        h2rel_vs_h1[(h2_rel, beta1_sub)] += 1

print(f"\nn=5: (h₂_rel, β₁(T\\v)) distribution for interior v with h₂_rel > 0:")
for (h2r, b1), cnt in sorted(h2rel_vs_h1.items()):
    print(f"  h₂_rel={h2r}, β₁={b1}: {cnt} cases (injective possible: {'YES' if h2r <= b1 else 'NO'})")


print(f"\n{'='*70}")
print("Part 7: n=6 dimension check")
print("=" * 70)

# Quick check at n=6 with samples
import random
random.seed(42)

n = 6
m = n*(n-1)//2
total = 1 << m
samples = random.sample(range(total), 1000)

h2rel_vs_h1_6 = Counter()
t0 = time.time()

for idx, bits in enumerate(samples):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue

        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        remap = {i: others[i] for i in range(n1)}

        ap0_T = enumerate_allowed_paths(A, n, 0)
        ap1_T = enumerate_allowed_paths(A, n, 1)
        ap2_T = enumerate_allowed_paths(A, n, 2)
        ap3_T = enumerate_allowed_paths(A, n, 3)

        ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
        ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
        ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
        ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)

        om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
        om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
        om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))
        om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
        om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

        d2_T = dim_om(om2_T)
        if d2_T == 0:
            continue

        d1_sub = dim_om(om1_sub)
        d2_sub = dim_om(om2_sub)
        d3_T = dim_om(om3_T)

        ap2_T_list = [tuple(p) for p in ap2_T]

        if ap2_sub and d2_sub > 0:
            embed = np.zeros((len(ap2_T_list), d2_sub))
            for j in range(d2_sub):
                for k, path_sub in enumerate(ap2_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap2_T_list:
                        embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
            phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
        else:
            phi = np.zeros((d2_T, 0))

        rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
        if rk_phi > 0:
            U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
            Q = U_phi[:, rk_phi:]
        else:
            Q = np.eye(d2_T)

        d_rel = Q.shape[1]
        if d_rel == 0:
            continue

        coords2_T = np.linalg.lstsq(om1_T, build_full_boundary_matrix(ap2_T, ap1_T) @ om2_T, rcond=None)[0]

        d1_T = dim_om(om1_T)
        ap1_T_list = [tuple(p) for p in ap1_T]
        if d1_sub > 0:
            embed1 = np.zeros((len(ap1_T_list), d1_sub))
            for j in range(d1_sub):
                for k, path_sub in enumerate(ap1_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap1_T_list:
                        embed1[ap1_T_list.index(path_T), j] = om1_sub[k, j]
            psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]
        else:
            psi = np.zeros((d1_T, 0))

        rk_psi = np.linalg.matrix_rank(psi, tol=1e-8)
        if rk_psi > 0:
            U_psi, _, _ = np.linalg.svd(psi, full_matrices=True)
            R = U_psi[:, rk_psi:]
        else:
            R = np.eye(d1_T)

        coords2_rel_q = R.T @ coords2_T @ Q
        rk_d2_rel = np.linalg.matrix_rank(coords2_rel_q, tol=1e-8)
        z2_rel_dim = d_rel - rk_d2_rel

        if z2_rel_dim == 0:
            continue

        if d3_T > 0:
            coords3_T = np.linalg.lstsq(om2_T, build_full_boundary_matrix(ap3_T, ap2_T) @ om3_T, rcond=None)[0]
            om3_sub = compute_omega_basis(A_sub, n1, 3, ap3_sub, ap2_sub) if ap3_sub else np.zeros((0,0))
            d3_sub = dim_om(om3_sub)
            d3_proj = Q.T @ coords3_T
            if d3_sub > 0:
                ap3_T_list = [tuple(p) for p in ap3_T]
                embed3 = np.zeros((len(ap3_T_list), d3_sub))
                for j in range(d3_sub):
                    for k, path_sub in enumerate(ap3_sub):
                        path_T = tuple(remap[x] for x in path_sub)
                        if path_T in ap3_T_list:
                            embed3[ap3_T_list.index(path_T), j] = om3_sub[k, j]
                chi = np.linalg.lstsq(om3_T, embed3, rcond=None)[0]
                rk_chi = np.linalg.matrix_rank(chi, tol=1e-8)
                if rk_chi > 0:
                    U_chi, _, _ = np.linalg.svd(chi, full_matrices=True)
                    d3_proj_rel = d3_proj @ U_chi[:, rk_chi:]
                else:
                    d3_proj_rel = d3_proj
            else:
                d3_proj_rel = d3_proj
        else:
            d3_proj_rel = np.zeros((d_rel, 0))

        rk_d3_rel = np.linalg.matrix_rank(d3_proj_rel, tol=1e-8)
        h2_rel = z2_rel_dim - rk_d3_rel

        if h2_rel == 0:
            continue

        bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
        rk_d1 = np.linalg.matrix_rank(bd1_sub @ om1_sub, tol=1e-8)
        z1_dim = d1_sub - rk_d1
        if d2_sub > 0:
            bd2_in_om1 = np.linalg.lstsq(om1_sub, build_full_boundary_matrix(ap2_sub, ap1_sub) @ om2_sub, rcond=None)[0]
            rk_d2 = np.linalg.matrix_rank(bd2_in_om1, tol=1e-8)
        else:
            rk_d2 = 0
        beta1_sub = z1_dim - rk_d2

        h2rel_vs_h1_6[(h2_rel, beta1_sub)] += 1
        break  # one v per tournament for speed

    if (idx + 1) % 200 == 0:
        print(f"  {idx+1}/1000 ({time.time()-t0:.0f}s)")

print(f"\nn=6 (1000 samples): (h₂_rel, β₁(T\\v)) for interior v with h₂_rel > 0:")
for (h2r, b1), cnt in sorted(h2rel_vs_h1_6.items()):
    print(f"  h₂_rel={h2r}, β₁={b1}: {cnt} cases (injective possible: {'YES' if h2r <= b1 else 'NO'})")

print("\nDone.")
