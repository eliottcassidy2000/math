#!/usr/bin/env python3
"""
beta2_cone_construction.py - Explicit preimage construction for Z2 in tournaments

Strategy: For each Z2 element z, construct w in Om3 with d3(w) = z
using a CONE-LIKE construction specific to tournaments.

Key idea: Pick a vertex v. For a 2-path (a,b,c) in z, try to extend it
to a 3-path by prepending/appending v:
  (v,a,b,c) or (a,b,c,v)

For these to be in A3: need v->a->b->c or a->b->c->v (allowed paths).
For these to be in Om3: the junk faces must cancel.

Alternative: "Star" construction.
For each edge (u,w) in the boundary of z that doesn't cancel,
find a vertex connecting them to create a chain in Om3.

Let's first understand EXACTLY what Z2 looks like and then
try different constructions.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, time
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


# ============================================================
# For each Z2 element, analyze its support structure
# ============================================================
print("=" * 70)
print("Z2 SUPPORT STRUCTURE AND CONE CONSTRUCTION")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

# Track which construction works
cone_success = 0
cone_fail = 0
star_success = 0
vertex_extend_success = 0

for bits in range(total):
    A = build_adj(n, bits)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if d_om2 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2, a1)
    bd3 = build_full_boundary_matrix(a3, a2)

    # Compute Z2 basis
    bd2_om = bd2 @ om2
    U, S_vals, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    rk = sum(s > 1e-8 for s in S_vals)
    dim_z2 = d_om2 - rk

    if dim_z2 == 0:
        continue

    z2_om_basis = Vt[rk:].T
    z2_a2_basis = om2 @ z2_om_basis

    # For each Z2 basis element, try cone construction
    for col in range(dim_z2):
        z = z2_a2_basis[:, col]

        # Method 1: For each vertex v, try "cone over v"
        # w_v = sum_{(a,b,c) in z with v->a} alpha*(v,a,b,c)
        #      - sum_{(a,b,c) in z with c->v} alpha*(a,b,c,v)  (sign from orientation)
        found_cone = False
        for v in range(n):
            # Build candidate w by extending each path in z
            w = np.zeros(len(a3))
            a3_idx = {p: i for i, p in enumerate(a3)}

            valid = True
            for i in range(len(a2)):
                if abs(z[i]) < 1e-8:
                    continue
                a, b, c = a2[i]
                if v in (a, b, c):
                    continue  # skip paths involving v

                # Try prepending v: (v,a,b,c)
                if A[v][a]:
                    key = (v, a, b, c)
                    if key in a3_idx:
                        w[a3_idx[key]] += z[i]

            # Check if w is in Om3
            if d_om3 > 0 and np.linalg.norm(w) > 1e-8:
                # Project w onto Om3
                # w = om3 @ coords + w_perp
                coords = np.linalg.lstsq(om3, w, rcond=None)[0]
                w_proj = om3 @ coords
                if np.linalg.norm(w - w_proj) < 1e-6:
                    # w is in Om3! Check if d3(w) = z
                    bd3_w = bd3 @ w
                    if np.linalg.norm(bd3_w - z) < 1e-6:
                        found_cone = True
                        cone_success += 1
                        break
                    elif np.linalg.norm(bd3_w) > 1e-8:
                        # d3(w) != z but d3(w) != 0
                        # Maybe d3(w) is proportional to z?
                        pass

        if not found_cone:
            cone_fail += 1

print(f"\nCone prepend construction: {cone_success} success, {cone_fail} fail")


# ============================================================
# Better approach: Try ALL possible vertex extension methods
# ============================================================
print(f"\n{'='*70}")
print("VERTEX EXTENSION METHODS")
print("=" * 70)

methods_work = Counter()
n = 5

for bits in range(total):
    A = build_adj(n, bits)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if d_om2 == 0 or d_om3 == 0:
        continue

    bd2 = build_full_boundary_matrix(a2, a1)
    bd3 = build_full_boundary_matrix(a3, a2)

    bd2_om = bd2 @ om2
    U, S_vals, Vt = np.linalg.svd(bd2_om, full_matrices=True)
    rk = sum(s > 1e-8 for s in S_vals)
    dim_z2 = d_om2 - rk

    if dim_z2 == 0:
        continue

    z2_om_basis = Vt[rk:].T
    z2_a2_basis = om2 @ z2_om_basis

    # For each Z2 element, check if the GENERAL Om3 surjects
    bd3_om = bd3 @ om3
    for col in range(dim_z2):
        z = z2_a2_basis[:, col]
        result = np.linalg.lstsq(bd3_om, z, rcond=None)
        residual = np.linalg.norm(bd3_om @ result[0] - z)
        if residual < 1e-6:
            methods_work['general_om3'] += 1
        else:
            methods_work['general_om3_fail'] += 1

print(f"General Om3 surjectivity: {methods_work}")


# ============================================================
# KEY: What is the RANK DEFICIENCY structure?
# rk(d3|_Om3) vs dim(Z2) -- when are they exactly equal?
# ============================================================
print(f"\n{'='*70}")
print("RANK EQUATION ANALYSIS")
print("=" * 70)

# For each tournament, compute the exact dimensions
rank_data = []
for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    bd2 = build_full_boundary_matrix(a2, a1)
    bd3 = build_full_boundary_matrix(a3, a2)

    rk2 = 0
    if d_om2 > 0:
        S = np.linalg.svd(bd2 @ om2, compute_uv=False)
        rk2 = sum(s > 1e-8 for s in S)

    rk3 = 0
    if d_om3 > 0:
        S = np.linalg.svd(bd3 @ om3, compute_uv=False)
        rk3 = sum(s > 1e-8 for s in S)

    z2 = d_om2 - rk2
    b2 = z2 - rk3

    # Kernel of d3|_Om3
    ker3 = d_om3 - rk3  # this is dim(Z3) = beta3 (since Om4 usually empty at n=5)

    rank_data.append({
        'bits': bits, 'scores': scores,
        'A2': len(a2), 'A3': len(a3),
        'Om2': d_om2, 'Om3': d_om3,
        'rk2': rk2, 'rk3': rk3,
        'Z2': z2, 'Z3': ker3,
        'b2': b2,
        # Key ratio: how much of Om3 is used for filling Z2?
        'fill_ratio': rk3 / d_om3 if d_om3 > 0 else 0,
        'surplus': d_om3 - z2 if d_om3 > 0 else -z2,  # Om3 surplus over Z2
    })

# Show by score type
by_score = defaultdict(list)
for d in rank_data:
    by_score[d['scores']].append(d)

print(f"{'Scores':<20} {'Om2':>4} {'rk2':>4} {'Z2':>4} {'Om3':>4} {'rk3':>4} {'Z3':>4} {'surplus':>7} {'fill%':>6}")
for sc in sorted(by_score.keys()):
    entries = by_score[sc]
    # Check if all same
    vals = set((d['Om2'], d['rk2'], d['Z2'], d['Om3'], d['rk3'], d['Z3']) for d in entries)
    if len(vals) == 1:
        d = entries[0]
        print(f"{str(sc):<20} {d['Om2']:>4} {d['rk2']:>4} {d['Z2']:>4} {d['Om3']:>4} {d['rk3']:>4} {d['Z3']:>4} {d['surplus']:>7} {d['fill_ratio']*100:>5.1f}%")
    else:
        for v in sorted(vals):
            cnt = sum(1 for d in entries if (d['Om2'], d['rk2'], d['Z2'], d['Om3'], d['rk3'], d['Z3']) == v)
            print(f"{str(sc):<20} {v[0]:>4} {v[1]:>4} {v[2]:>4} {v[3]:>4} {v[4]:>4} {v[5]:>4}   ({cnt}x)")


# ============================================================
# CRITICAL OBSERVATION: Is dim(Om3) >= dim(Z2) always?
# ============================================================
print(f"\n{'='*70}")
print("IS dim(Om3) >= dim(Z2) ALWAYS?")
print("=" * 70)

always_ge = True
for d in rank_data:
    if d['Om3'] < d['Z2']:
        print(f"  FAIL at bits={d['bits']}: Om3={d['Om3']} < Z2={d['Z2']}")
        always_ge = False

if always_ge:
    print("  YES: dim(Om3) >= dim(Z2) for all n=5 tournaments")
    surplus_vals = Counter(d['surplus'] for d in rank_data)
    print(f"  Surplus distribution: {dict(sorted(surplus_vals.items()))}")


# ============================================================
# Now check n=6 (sample)
# ============================================================
print(f"\n{'='*70}")
print("DIMENSION CHECK AT n=6 (sample)")
print("=" * 70)

import random
random.seed(42)
n = 6
total6 = 1 << (n*(n-1)//2)
sample_size = 500

min_surplus = float('inf')
surplus_dist = Counter()
fill_dist = []

for _ in range(sample_size):
    bits = random.randint(0, total6 - 1)
    A = build_adj(n, bits)
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om2 = om2.shape[1] if om2.ndim == 2 else 0
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    bd2 = build_full_boundary_matrix(a2, a1)
    bd3 = build_full_boundary_matrix(a3, a2)

    rk2 = 0
    if d_om2 > 0:
        S = np.linalg.svd(bd2 @ om2, compute_uv=False)
        rk2 = sum(s > 1e-8 for s in S)

    rk3 = 0
    if d_om3 > 0:
        S = np.linalg.svd(bd3 @ om3, compute_uv=False)
        rk3 = sum(s > 1e-8 for s in S)

    z2 = d_om2 - rk2
    surplus = d_om3 - z2
    min_surplus = min(min_surplus, surplus)
    surplus_dist[surplus] += 1
    if d_om3 > 0:
        fill_dist.append(rk3 / d_om3)

print(f"n=6, {sample_size} samples:")
print(f"  Min surplus (Om3 - Z2): {min_surplus}")
print(f"  Surplus distribution: {dict(sorted(surplus_dist.items()))}")
if fill_dist:
    print(f"  Fill ratio (rk3/Om3): min={min(fill_dist):.3f}, max={max(fill_dist):.3f}, avg={sum(fill_dist)/len(fill_dist):.3f}")


# ============================================================
# DEEPER: What's the structure of ker(d3|_Om3)?
# At n=5, Z3 = ker(d3|_Om3) should be 0 (since beta3=0 at n=5)
# ============================================================
print(f"\n{'='*70}")
print("KERNEL OF d3|_Om3 AT n=5")
print("=" * 70)

n = 5
z3_dist = Counter()
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    d_om3 = om3.shape[1] if om3.ndim == 2 else 0

    if d_om3 == 0:
        z3_dist[0] += 1
        continue

    bd3 = build_full_boundary_matrix(a3, a2)
    S = np.linalg.svd(bd3 @ om3, compute_uv=False)
    rk3 = sum(s > 1e-8 for s in S)
    z3 = d_om3 - rk3
    z3_dist[z3] += 1

print(f"dim(Z3) distribution at n=5: {dict(sorted(z3_dist.items()))}")
print(f"beta3 = Z3 (since Om4 is usually empty at n=5)")

# At n=5: beta3=0 means Z3=0 means d3|_Om3 is INJECTIVE!
# So Om3 injects into Z2 via d3.
# beta2=0 means this injection is also SURJECTIVE onto Z2.
# So d3: Om3 -> Z2 is an ISOMORPHISM at n=5!

any_z3_nonzero = any(v > 0 for v in z3_dist.keys() if v in z3_dist)
if z3_dist.get(0, 0) == 1 << (n*(n-1)//2):
    print("\n  REMARKABLE: d3|_Om3 is INJECTIVE for ALL n=5 tournaments!")
    print("  Combined with beta2=0 (surjectivity): d3: Om3 -> Z2 is ISOMORPHISM")
    print("  This means dim(Om3) = dim(Z2) exactly!")

    # Verify
    exact_eq = sum(1 for d in rank_data if d['Om3'] == d['Z2'])
    print(f"  Verification: Om3 = Z2 in {exact_eq}/{len(rank_data)} cases")


print("\n\nDone.")
