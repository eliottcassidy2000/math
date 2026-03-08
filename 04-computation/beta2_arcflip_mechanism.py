#!/usr/bin/env python3
"""
beta2_arcflip_mechanism.py - Understand the algebraic mechanism of arc-flip β₂=0 preservation

For an arc flip (u→v) → (v→u), track EXACTLY which allowed paths are lost/gained
at each dimension, and how the Omega spaces and boundaries respond.

The key structural question: WHY does δΩ₃ ≥ δZ₂ always hold for tournaments?

Strategy: decompose allowed paths into:
  - Paths NOT using the flipped arc (unchanged)
  - Paths using u→v (lost)
  - Paths using v→u (gained)

Then analyze how the Omega space and boundary maps change.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_full_data(A, n, max_dim=4):
    allowed = {}
    for p in range(max_dim + 1):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break

    omega_basis = {}
    for p in range(max_dim):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p >= 1 and p-1 in allowed else [])
        omega_basis[p] = basis

    dims = {}
    for p in range(max_dim):
        dims[p] = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0

    # Compute Z_2 and im(d_3)
    if dims.get(2, 0) == 0:
        Z2 = 0
    else:
        bd2 = build_full_boundary_matrix(allowed[2], allowed[1] if 1 in allowed else [])
        bd2_om = bd2 @ omega_basis[2]
        S_v = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = int(np.sum(np.abs(S_v) > 1e-8))
        Z2 = dims[2] - rk2

    if dims.get(3, 0) == 0 or 3 not in allowed or not allowed[3]:
        im_d3 = 0
    else:
        bd3 = build_full_boundary_matrix(allowed[3], allowed[2] if 2 in allowed else [])
        bd3_om = bd3 @ omega_basis[3]
        S_v = np.linalg.svd(bd3_om, compute_uv=False)
        im_d3 = int(np.sum(np.abs(S_v) > 1e-8))

    return {
        'allowed': allowed, 'omega_basis': omega_basis, 'dims': dims,
        'Z2': Z2, 'im_d3': im_d3,
        'beta2': Z2 - im_d3,
        'surplus': dims.get(3, 0) - Z2
    }


def paths_using_arc(paths, u, v):
    """Count/list paths that use the arc u→v."""
    using = []
    for p in paths:
        for i in range(len(p) - 1):
            if p[i] == u and p[i+1] == v:
                using.append(p)
                break
    return using


n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"ARC-FLIP MECHANISM ANALYSIS AT n={n}")
print("=" * 70)

# For each surplus=0 tournament, analyze the exact path changes under each flip
surplus0_bits = []
all_data = {}

for bits in range(total):
    A = build_adj(n, bits)
    data = compute_full_data(A, n)
    all_data[bits] = data
    if data['surplus'] == 0:
        surplus0_bits.append(bits)

print(f"\n  Surplus=0 tournaments: {len(surplus0_bits)}")
print(f"\n  Analyzing path-level changes under arc flips from surplus=0...")

# Track: for each flip from surplus=0, how many allowed 3-paths and 4-paths
# use the flipped arc?
path_analysis = []

for bits in surplus0_bits[:20]:  # Detailed analysis of first 20
    A = build_adj(n, bits)
    data_before = all_data[bits]

    for u in range(n):
        for v in range(n):
            if u == v or A[u][v] == 0:
                continue

            # Flip u→v to v→u
            B = [row[:] for row in A]
            B[u][v] = 0
            B[v][u] = 1

            # Recompute bits for B
            bits_flip = 0
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if B[i][j] == 1:
                        bits_flip |= (1 << idx)
                    idx += 1

            data_after = all_data[bits_flip]

            # Count paths using the flipped arc at each dimension
            lost2 = paths_using_arc(data_before['allowed'].get(2, []), u, v)
            lost3 = paths_using_arc(data_before['allowed'].get(3, []), u, v)
            gained2 = paths_using_arc(data_after['allowed'].get(2, []), v, u)
            gained3 = paths_using_arc(data_after['allowed'].get(3, []), v, u)

            delta_A2 = len(data_after['allowed'].get(2, [])) - len(data_before['allowed'].get(2, []))
            delta_A3 = len(data_after['allowed'].get(3, [])) - len(data_before['allowed'].get(3, []))
            delta_O2 = data_after['dims'].get(2, 0) - data_before['dims'].get(2, 0)
            delta_O3 = data_after['dims'].get(3, 0) - data_before['dims'].get(3, 0)
            delta_Z2 = data_after['Z2'] - data_before['Z2']

            path_analysis.append({
                'bits': bits, 'u': u, 'v': v,
                'lost2': len(lost2), 'gained2': len(gained2), 'delta_A2': delta_A2,
                'lost3': len(lost3), 'gained3': len(gained3), 'delta_A3': delta_A3,
                'delta_O2': delta_O2, 'delta_O3': delta_O3, 'delta_Z2': delta_Z2,
                'surplus_after': data_after['surplus']
            })

# Summarize
print(f"\n  Total flips analyzed: {len(path_analysis)}")

# Joint distribution of (lost3, gained3)
lost_gained3 = Counter((d['lost3'], d['gained3']) for d in path_analysis)
print(f"\n  (lost_3paths, gained_3paths) distribution:")
for (l, g) in sorted(lost_gained3.keys()):
    print(f"    ({l}, {g}): {lost_gained3[(l,g)]}")

# Joint distribution of (lost2, gained2)
lost_gained2 = Counter((d['lost2'], d['gained2']) for d in path_analysis)
print(f"\n  (lost_2paths, gained_2paths) distribution:")
for (l, g) in sorted(lost_gained2.keys()):
    print(f"    ({l}, {g}): {lost_gained2[(l,g)]}")

# Joint distribution of (delta_A2, delta_A3)
delta_A = Counter((d['delta_A2'], d['delta_A3']) for d in path_analysis)
print(f"\n  (delta_|A_2|, delta_|A_3|) distribution:")
for (d2, d3) in sorted(delta_A.keys()):
    print(f"    ({d2:+d}, {d3:+d}): {delta_A[(d2,d3)]}")

# KEY: joint distribution of (delta_O3, delta_Z2, delta_O2)
joint_OZ = Counter((d['delta_O3'], d['delta_Z2'], d['delta_O2']) for d in path_analysis)
print(f"\n  (delta_O3, delta_Z2, delta_O2) distribution:")
for (dO3, dZ2, dO2) in sorted(joint_OZ.keys()):
    print(f"    ({dO3:+d}, {dZ2:+d}, {dO2:+d}): {joint_OZ[(dO3, dZ2, dO2)]}")

# Now do full analysis: ALL tournaments, track flipped-arc vertex degree info
print(f"\n{'='*70}")
print(f"FULL ANALYSIS: VERTEX DEGREE AND SURPLUS CHANGE")
print(f"{'='*70}")

# For each flip, record (out_degree(u), out_degree(v), surplus_before, surplus_after)
degree_surplus = []
for bits in range(total):
    A = build_adj(n, bits)
    data_before = all_data[bits]

    for u in range(n):
        for v in range(n):
            if u == v or A[u][v] == 0:
                continue

            B = [row[:] for row in A]
            B[u][v] = 0
            B[v][u] = 1
            bits_flip = 0
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if B[i][j] == 1:
                        bits_flip |= (1 << idx)
                    idx += 1

            data_after = all_data[bits_flip]

            out_u = sum(A[u])
            out_v = sum(A[v])

            degree_surplus.append({
                'out_u': out_u, 'out_v': out_v,
                'surplus_before': data_before['surplus'],
                'surplus_after': data_after['surplus'],
                'delta': data_after['surplus'] - data_before['surplus']
            })

# Minimum surplus by (out_u, out_v) pair
min_delta_by_deg = defaultdict(lambda: float('inf'))
count_by_deg = defaultdict(int)
for d in degree_surplus:
    key = (d['out_u'], d['out_v'])
    min_delta_by_deg[key] = min(min_delta_by_deg[key], d['delta'])
    count_by_deg[key] += 1

print(f"\n  Min delta surplus by (out_u, out_v):")
for key in sorted(min_delta_by_deg.keys()):
    print(f"    deg(u)={key[0]}, deg(v)={key[1]}: min_delta={min_delta_by_deg[key]}, count={count_by_deg[key]}")

# Minimum surplus after flip, by starting surplus
min_after_by_start = defaultdict(lambda: float('inf'))
for d in degree_surplus:
    min_after_by_start[d['surplus_before']] = min(
        min_after_by_start[d['surplus_before']], d['surplus_after'])

print(f"\n  Min surplus after flip, by starting surplus:")
for s in sorted(min_after_by_start.keys()):
    print(f"    surplus={s}: min_after={min_after_by_start[s]}")

# Count 3-cycles through (u,v) before and after
print(f"\n{'='*70}")
print(f"3-CYCLE ANALYSIS: FLIPS THROUGH (u,v)")
print(f"{'='*70}")

def count_3cycles_through_uv(A, n, u, v):
    """Count 3-cycles using arc u→v."""
    count = 0
    for w in range(n):
        if w == u or w == v:
            continue
        # u→v→w→u
        if A[u][v] and A[v][w] and A[w][u]:
            count += 1
        # w→u→v→w
        if A[w][u] and A[u][v] and A[v][w]:
            count += 1
    return count

cycle_data = []
for bits in surplus0_bits[:50]:
    A = build_adj(n, bits)
    for u in range(n):
        for v in range(n):
            if u == v or A[u][v] == 0:
                continue
            c3_through = count_3cycles_through_uv(A, n, u, v)

            B = [row[:] for row in A]
            B[u][v] = 0
            B[v][u] = 1
            bits_flip = 0
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if B[i][j] == 1:
                        bits_flip |= (1 << idx)
                    idx += 1

            c3_after_through = count_3cycles_through_uv(B, n, v, u)

            data_after = all_data[bits_flip]
            delta_surplus = data_after['surplus'] - all_data[bits]['surplus']

            cycle_data.append({
                'c3_before': c3_through,
                'c3_after': c3_after_through,
                'delta_surplus': delta_surplus
            })

# Joint distribution of (c3_through_before, c3_through_after, delta_surplus)
c3_joint = Counter((d['c3_before'], d['c3_after'], d['delta_surplus'])
                    for d in cycle_data)
print(f"\n  (c3_through_before, c3_through_after, delta_surplus):")
for key in sorted(c3_joint.keys()):
    print(f"    ({key[0]}, {key[1]}, {key[2]:+d}): {c3_joint[key]}")

print("\nDone.")
