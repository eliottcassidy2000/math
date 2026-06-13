#!/usr/bin/env python3
"""
beta2_cone_filling_formula.py - Test the "front+back cone" filling formula

Key observation from explicit filling at n=5:
The filling w for swap cycle z = sum M_{ab}[(a,b,v)-(v,a,b)] uses
3-paths with v ONLY at positions 0 and 3.

Hypothesis: w = (1/c) * [Cone_front(z) + Cone_back(z)]
where:
  Cone_front(a,b,v) = (?, a, b, v)  -- prepend some vertex to (a,b,v)
  Cone_back(v,a,b) = (v, a, b, ?)   -- append some vertex to (v,a,b)

But which vertex? For (a,b,v), we prepend g with g->a.
For (v,a,b), we append c with b->c.

Actually, let's define:
  h_front(a,b,v) = sum_{g: g->a, g!=b,v} (g,a,b,v)
  h_back(v,a,b) = sum_{c: b->c, c!=v,a} (v,a,b,c)

Then maybe w = (1/c) * [h_front(z_type2) + h_back(z_type0)]

Let me test this systematically.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


# ============================================================
# PART 1: Explicit cone construction at n=5
# ============================================================
print("=" * 70)
print("CONE FILLING FORMULA TEST")
print("=" * 70)

n = 5

# Use the specific example from the filling: T#12, v=0
bits = 12
A = build_adj(n, bits)
v = 0
P = [a for a in range(n) if a != v and A[v][a] == 1]
Q = [b for b in range(n) if b != v and A[b][v] == 1]

print(f"T#{bits}, v={v}, P={P}, Q={Q}")
print(f"Adjacency:")
for i in range(n):
    row = ''.join(str(A[i][j]) for j in range(n))
    print(f"  {i}: {row}")

# The swap cycle: z = (0,4,2)-(0,4,1)-(0,3,2)+(0,3,1) + (3,2,0)-(3,1,0)-(4,2,0)+(4,1,0)
# Type 0 (v at front): (0,4,2), (0,4,1), (0,3,2), (0,3,1) with coeffs +,-,-,+
# Type 2 (v at back): (3,2,0), (3,1,0), (4,2,0), (4,1,0) with coeffs +,-,-,+

# Cone_back: for each (v,a,b) in swap, create (v,a,b,c) for c with b->c, c!=v,a
print(f"\nCone_back construction:")
# Type 0 paths: (v,a,b) = (0,4,2), (0,4,1), (0,3,2), (0,3,1)
type0_paths = [(0,4,2,+0.5), (0,4,1,-0.5), (0,3,2,-0.5), (0,3,1,+0.5)]
for (v_, a, b, coeff) in type0_paths:
    successors = [c for c in range(n) if c != v_ and c != a and A[b][c] == 1]
    print(f"  ({v_},{a},{b}) [coeff={coeff:+.2f}]: b={b}, b->{{c}} = {successors}")
    for c in successors:
        print(f"    -> ({v_},{a},{b},{c}) coeff={coeff:+.4f}")

# Cone_front: for each (a,b,v) in swap, create (g,a,b,v) for g with g->a, g!=b,v
print(f"\nCone_front construction:")
type2_paths = [(3,2,0,+0.5), (3,1,0,-0.5), (4,2,0,-0.5), (4,1,0,+0.5)]
for (a, b, v_, coeff) in type2_paths:
    predecessors = [g for g in range(n) if g != v_ and g != b and A[g][a] == 1]
    print(f"  ({a},{b},{v_}) [coeff={coeff:+.2f}]: a={a}, {{g}}->a = {predecessors}")
    for g in predecessors:
        print(f"    -> ({g},{a},{b},{v_}) coeff={coeff:+.4f}")


# ============================================================
# PART 2: Build the cone 3-chain and compare with actual filling
# ============================================================
print(f"\n{'='*70}")
print("COMPARISON: CONE vs ACTUAL FILLING")
print("=" * 70)

paths2 = enumerate_allowed_paths(A, n, 2)
paths3 = enumerate_allowed_paths(A, n, 3)
paths1 = enumerate_allowed_paths(A, n, 1)
omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
path2_idx = {p: i for i, p in enumerate(paths2)}
path3_idx = {p: i for i, p in enumerate(paths3)}

# Build swap cycle z
M = {(3,1): -0.5, (3,2): 0.5, (4,1): 0.5, (4,2): -0.5}
z = np.zeros(len(paths2))
for (a, b), coeff in M.items():
    if (a,b,v) in path2_idx: z[path2_idx[(a,b,v)]] += coeff
    if (v,a,b) in path2_idx: z[path2_idx[(v,a,b)]] -= coeff

# Build cone_back: for (v,a,b) with coeff c, add c*(v,a,b,c') for each c' with b->c'
w_cone = np.zeros(len(paths3))
for (a, b), coeff in M.items():
    # Type 0: (v,a,b) has coefficient -coeff (negative in swap cycle)
    c_type0 = -coeff
    successors = [c for c in range(n) if c != v and c != a and A[b][c] == 1]
    for c in successors:
        path = (v, a, b, c)
        if path in path3_idx:
            w_cone[path3_idx[path]] += c_type0

    # Type 2: (a,b,v) has coefficient +coeff
    c_type2 = coeff
    predecessors = [g for g in range(n) if g != v and g != b and A[g][a] == 1]
    for g in predecessors:
        path = (g, a, b, v)
        if path in path3_idx:
            w_cone[path3_idx[path]] += c_type2

# Check d_3(w_cone) vs z
D3 = build_full_boundary_matrix(paths3, paths2)
d3_w_cone = D3 @ w_cone
err = np.max(np.abs(d3_w_cone - z))

print(f"Cone 3-chain d_3 error: {err:.2e}")

if err > 1e-6:
    print(f"  Cone does NOT exactly give z!")
    print(f"  Difference d_3(w_cone) - z:")
    for j, path in enumerate(paths2):
        diff = d3_w_cone[j] - z[j]
        if abs(diff) > 1e-10:
            print(f"    [{path}]: cone={d3_w_cone[j]:+.4f}, z={z[j]:+.4f}, diff={diff:+.4f}")

    # What IS d_3(w_cone)?
    print(f"\n  d_3(w_cone):")
    for j, path in enumerate(paths2):
        if abs(d3_w_cone[j]) > 1e-10:
            print(f"    [{path}]: {d3_w_cone[j]:+.4f}")

    # And the actual filling:
    D3_omega = D3 @ omega3
    w_actual, _, _, _ = np.linalg.lstsq(D3_omega, z, rcond=None)
    w_actual_A3 = omega3 @ w_actual
    print(f"\n  Actual filling:")
    for j, path in enumerate(paths3):
        if abs(w_actual_A3[j]) > 1e-10:
            print(f"    [{path}]: {w_actual_A3[j]:+.4f}")

    # What's the discrepancy? d_3(w_cone) - z = ?
    # This is a 2-chain in A_2. Is it a boundary? Is it in Omega_2?
    disc = d3_w_cone - z
    # Check if disc is in Omega_2
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    proj, _, _, _ = np.linalg.lstsq(omega2, disc, rcond=None)
    disc_proj = omega2 @ proj
    omega_err = np.max(np.abs(disc - disc_proj))
    print(f"\n  Discrepancy in Omega_2? err={omega_err:.2e}")

    # Is w_cone in Omega_3?
    proj3, _, _, _ = np.linalg.lstsq(omega3, w_cone, rcond=None)
    w_cone_proj = omega3 @ proj3
    omega3_err = np.max(np.abs(w_cone - w_cone_proj))
    print(f"  w_cone in Omega_3? err={omega3_err:.2e}")

else:
    print(f"  CONE EXACTLY FILLS THE SWAP CYCLE!")

    # But is w_cone in Omega_3?
    proj3, _, _, _ = np.linalg.lstsq(omega3, w_cone, rcond=None)
    w_cone_proj = omega3 @ proj3
    omega3_err = np.max(np.abs(w_cone - w_cone_proj))
    print(f"  w_cone in Omega_3? err={omega3_err:.2e}")


# ============================================================
# PART 3: Try different normalizations / averaging
# ============================================================
print(f"\n{'='*70}")
print("NORMALIZED CONE CONSTRUCTIONS")
print("=" * 70)

# Maybe: w = (1/k) * [cone_front + cone_back] for some k related to degree
# Or: average over ALL possible insertion vertices

# For a 2-path (a,b,c), the cone operator inserts v at ALL positions:
# K_0: (v, a, b, c) if v->a
# K_1: (a, v, b, c) if a->v and v->b
# K_2: (a, b, v, c) if b->v and v->c
# K_3: (a, b, c, v) if c->v

# For swap cycle z = sum M_{ab}[(a,b,v) - (v,a,b)]:
# Apply K at all positions and see what we get.

# For (a,b,v) [v at pos 2]:
#   K_0: (g,a,b,v) for g->a, g!=b,v -- "add before"
#   K_1: (a,w,b,v) for a->w, w->b, w!=v -- "add between a and b"
#   K_2: (a,b,w,v) for b->w, w->v, w!=a -- "add between b and v"
#   K_3: not applicable (v is already at end)

# For (v,a,b) [v at pos 0]:
#   K_0: not applicable (v is already at start)
#   K_1: (v,w,a,b) for v->w, w->a, w!=b -- "add between v and a"
#   K_2: (v,a,w,b) for a->w, w->b, w!=v -- "add between a and b"
#   K_3: (v,a,b,c) for b->c, c!=v,a -- "add after"

# Let me try sum of K_0(type2) + K_3(type0):
# This is exactly the front+back cone from before.

# Let me also try sum of ALL insertions.
w_all = np.zeros(len(paths3))

for (a, b), coeff in M.items():
    # For (a,b,v) with coefficient +coeff:
    c_abv = coeff

    # K_0: (g,a,b,v) for g->a, g not in {b,v}
    for g in range(n):
        if g != a and g != b and g != v and A[g][a] == 1:
            path = (g,a,b,v)
            if path in path3_idx:
                w_all[path3_idx[path]] += c_abv

    # K_1: (a,w,b,v) for a->w->b, w not in {v}
    for w in range(n):
        if w != a and w != b and w != v and A[a][w] == 1 and A[w][b] == 1:
            path = (a,w,b,v)
            if path in path3_idx:
                w_all[path3_idx[path]] += c_abv

    # K_2: (a,b,w,v) for b->w->v, w not in {a}
    # Wait: b->v is given (b in Q), so we need b->w->v? That means b->w AND w->v,
    # but v is NOT a successor of w in general. w->v means w is in Q? No, w->v means
    # A[w][v]=1, so v->w is false, meaning w is in Q.
    for w in range(n):
        if w != a and w != b and w != v and A[b][w] == 1 and A[w][v] == 1:
            path = (a,b,w,v)
            if path in path3_idx:
                w_all[path3_idx[path]] += c_abv

    # For (v,a,b) with coefficient -coeff:
    c_vab = -coeff

    # K_1: (v,w,a,b) for v->w->a, w not in {b}
    for w in range(n):
        if w != a and w != b and w != v and A[v][w] == 1 and A[w][a] == 1:
            path = (v,w,a,b)
            if path in path3_idx:
                w_all[path3_idx[path]] += c_vab

    # K_2: (v,a,w,b) for a->w->b, w not in {v}
    for w in range(n):
        if w != a and w != b and w != v and A[a][w] == 1 and A[w][b] == 1:
            path = (v,a,w,b)
            if path in path3_idx:
                w_all[path3_idx[path]] += c_vab

    # K_3: (v,a,b,c) for b->c, c not in {v,a}
    for c in range(n):
        if c != a and c != b and c != v and A[b][c] == 1:
            path = (v,a,b,c)
            if path in path3_idx:
                w_all[path3_idx[path]] += c_vab

d3_w_all = D3 @ w_all
err_all = np.max(np.abs(d3_w_all - z))
print(f"All-position cone d_3 error: {err_all:.2e}")

if err_all > 1e-6:
    # Discrepancy
    disc_all = d3_w_all - z
    print(f"  Discrepancy norm: {np.linalg.norm(disc_all):.4f}")
    print(f"  Nonzero discrepancy terms:")
    for j, path in enumerate(paths2):
        if abs(disc_all[j]) > 1e-10:
            print(f"    [{path}]: {disc_all[j]:+.4f}")

    # What factor would make it work?
    # d_3(w_all) = z + disc
    # If we scale: d_3(alpha * w_all) = alpha * z + alpha * disc
    # We want alpha * z + alpha * disc = z
    # So alpha * (z + disc) = z => alpha = z / (z + disc) -- not constant!

    # But maybe disc is proportional to z?
    if np.linalg.norm(z) > 1e-10:
        ratio = d3_w_all / (z + 1e-20)
        nonzero_ratios = [ratio[j] for j in range(len(z)) if abs(z[j]) > 1e-10]
        if nonzero_ratios:
            print(f"\n  d_3(w_all)/z ratios (where z!=0): {[f'{r:.4f}' for r in nonzero_ratios]}")
            if len(set(round(r, 4) for r in nonzero_ratios)) == 1:
                scale = nonzero_ratios[0]
                print(f"  CONSTANT RATIO = {scale}!")
                print(f"  => w = (1/{scale}) * w_all would work")
else:
    print(f"  ALL-POSITION CONE EXACTLY FILLS!")

# Try: w = (1/(n-2)) * [front + back cone]
# This would generalize if the cone operator has a uniform degree
for alpha in [1.0, 0.5, 1/3, 0.25, 1/(n-1), 1/(n-2), 2/(n-1), 2/(n-2)]:
    w_scaled = alpha * w_cone
    err_scaled = np.max(np.abs(D3 @ w_scaled - z))
    if err_scaled < 1e-6:
        print(f"\n  FOUND: alpha={alpha:.6f} ({alpha}) gives exact filling!")
        break


# ============================================================
# PART 4: Test front+back cone for ALL n=5 tournaments
# ============================================================
print(f"\n{'='*70}")
print("FRONT+BACK CONE TEST FOR ALL n=5")
print("=" * 70)

n = 5
success = 0
fail = 0
total = 0

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    path2_idx = {p: i for i, p in enumerate(paths2)}
    path3_idx = {p: i for i, p in enumerate(paths3)}
    D3 = build_full_boundary_matrix(paths3, paths2)

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            continue

        m = len(arcs_PQ)
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue

        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C

        if ker_dim == 0:
            continue

        _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)

        for ki in range(ker_basis.shape[0] if 'ker_basis' in dir() else Vt.shape[0] - rank_C):
            M_vec = Vt[rank_C + ki]

            # Build swap cycle z
            z = np.zeros(len(paths2))
            for j, (a, b) in enumerate(arcs_PQ):
                coeff = M_vec[j]
                if abs(coeff) < 1e-12: continue
                if (a,b,v) in path2_idx: z[path2_idx[(a,b,v)]] += coeff
                if (v,a,b) in path2_idx: z[path2_idx[(v,a,b)]] -= coeff

            if np.max(np.abs(z)) < 1e-12:
                continue

            total += 1

            # Build front+back cone
            w_cone = np.zeros(len(paths3))
            for j, (a, b) in enumerate(arcs_PQ):
                coeff = M_vec[j]
                if abs(coeff) < 1e-12: continue

                # Front cone for type 2: (g,a,b,v) for g->a
                for g in range(n):
                    if g != a and g != b and g != v and A[g][a] == 1:
                        path = (g,a,b,v)
                        if path in path3_idx:
                            w_cone[path3_idx[path]] += coeff

                # Back cone for type 0: (v,a,b,c) for b->c
                for c in range(n):
                    if c != a and c != b and c != v and A[b][c] == 1:
                        path = (v,a,b,c)
                        if path in path3_idx:
                            w_cone[path3_idx[path]] -= coeff  # type 0 has -coeff

            d3_w = D3 @ w_cone
            err = np.max(np.abs(d3_w - z))

            if err < 1e-6:
                success += 1
            else:
                fail += 1
                # Try scaling
                found_scale = False
                for alpha in [0.5, 1/3, 0.25, 1/(n-1), 1/(n-2)]:
                    err_s = np.max(np.abs(D3 @ (alpha * w_cone) - z))
                    if err_s < 1e-6:
                        success += 1
                        fail -= 1
                        found_scale = True
                        break
                if not found_scale and fail <= 3:
                    scores = tuple(sorted([sum(row) for row in A]))
                    print(f"  FAIL: T#{bits} v={v} scores={scores}: err={err:.2e}")
                    # Check ratio
                    nonzero = [(j, d3_w[j], z[j]) for j in range(len(z)) if abs(z[j]) > 1e-10]
                    if nonzero:
                        ratios = [dw/zv for j, dw, zv in nonzero]
                        print(f"    Ratios d3(w)/z: {[f'{r:.4f}' for r in ratios[:5]]}")

print(f"\nn=5: {total} swap cycles, {success} exact filling, {fail} failures")
print(f"  Front+back cone works? {'YES' if fail == 0 else 'NO'}")


print("\n\nDone.")
