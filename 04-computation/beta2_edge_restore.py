#!/usr/bin/env python3
"""
beta2_edge_restore.py - How does restoring a removed edge kill beta2?

When we remove edge u->v from tournament T, beta2 can become >0.
The "missing" edge creates a Z2 cycle that can't be filled.
Restoring the edge creates new Om3 elements that fill this cycle.

Strategy: For each edge removal that creates beta2>0:
1. Compute the "obstruction" = coker of d3 in T'
2. See exactly which Om3 elements are created when we restore u->v
3. Check that these new elements fill the obstruction

This should reveal the MECHANISM by which completeness forces beta2=0.

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
    build_full_boundary_matrix, path_betti_numbers
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
# Find ALL edge removals at n=5 that create beta2 > 0
# ============================================================
print("=" * 70)
print("EDGE REMOVAL -> RESTORATION ANALYSIS")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

removals_creating_beta2 = []

for bits in range(total):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            A_mod = [row[:] for row in A]
            A_mod[u][v] = 0
            betti = path_betti_numbers(A_mod, n, max_dim=2)
            if betti[2] > 0:
                removals_creating_beta2.append({
                    'bits': bits, 'scores': scores,
                    'u': u, 'v': v,
                    'du': sum(A[u]), 'dv': sum(A[v]),
                    'beta2': betti[2]
                })

print(f"Total removals creating beta2>0: {len(removals_creating_beta2)}")


# ============================================================
# For each such removal, analyze the obstruction in detail
# ============================================================
print(f"\n{'='*70}")
print("OBSTRUCTION ANALYSIS")
print("=" * 70)

analyzed = 0
for removal in removals_creating_beta2[:10]:  # First 10 cases
    bits = removal['bits']
    u, v = removal['u'], removal['v']

    # Full tournament T
    A_full = build_adj(n, bits)

    # Almost-tournament T' (edge u->v removed)
    A_mod = [row[:] for row in A_full]
    A_mod[u][v] = 0

    # Compute chain data for T'
    a1_mod = enumerate_allowed_paths(A_mod, n, 1)
    a2_mod = enumerate_allowed_paths(A_mod, n, 2)
    a3_mod = enumerate_allowed_paths(A_mod, n, 3)

    om2_mod = compute_omega_basis(A_mod, n, 2, a2_mod, a1_mod)
    om3_mod = compute_omega_basis(A_mod, n, 3, a3_mod, a2_mod)
    d_om2_mod = om2_mod.shape[1] if om2_mod.ndim == 2 else 0
    d_om3_mod = om3_mod.shape[1] if om3_mod.ndim == 2 else 0

    bd2_mod = build_full_boundary_matrix(a2_mod, a1_mod)
    bd3_mod = build_full_boundary_matrix(a3_mod, a2_mod)

    rk2_mod = 0
    if d_om2_mod > 0:
        S = np.linalg.svd(bd2_mod @ om2_mod, compute_uv=False)
        rk2_mod = sum(s > 1e-8 for s in S)

    rk3_mod = 0
    if d_om3_mod > 0:
        S = np.linalg.svd(bd3_mod @ om3_mod, compute_uv=False)
        rk3_mod = sum(s > 1e-8 for s in S)

    z2_mod = d_om2_mod - rk2_mod

    # Compute chain data for full T
    a1_full = enumerate_allowed_paths(A_full, n, 1)
    a2_full = enumerate_allowed_paths(A_full, n, 2)
    a3_full = enumerate_allowed_paths(A_full, n, 3)

    om2_full = compute_omega_basis(A_full, n, 2, a2_full, a1_full)
    om3_full = compute_omega_basis(A_full, n, 3, a3_full, a2_full)
    d_om2_full = om2_full.shape[1] if om2_full.ndim == 2 else 0
    d_om3_full = om3_full.shape[1] if om3_full.ndim == 2 else 0

    bd2_full = build_full_boundary_matrix(a2_full, a1_full)
    bd3_full = build_full_boundary_matrix(a3_full, a2_full)

    rk2_full = 0
    if d_om2_full > 0:
        S = np.linalg.svd(bd2_full @ om2_full, compute_uv=False)
        rk2_full = sum(s > 1e-8 for s in S)

    rk3_full = 0
    if d_om3_full > 0:
        S = np.linalg.svd(bd3_full @ om3_full, compute_uv=False)
        rk3_full = sum(s > 1e-8 for s in S)

    z2_full = d_om2_full - rk2_full

    # What 2-paths are GAINED when restoring u->v?
    a2_full_set = set(a2_full)
    a2_mod_set = set(a2_mod)
    new_2paths = a2_full_set - a2_mod_set
    lost_2paths = a2_mod_set - a2_full_set

    # What 3-paths are GAINED?
    a3_full_set = set(a3_full)
    a3_mod_set = set(a3_mod)
    new_3paths = a3_full_set - a3_mod_set
    lost_3paths = a3_mod_set - a3_full_set

    print(f"\nbits={bits}, remove {u}->{v} (d_u={removal['du']}, d_v={removal['dv']})")
    print(f"  T': Om2={d_om2_mod}, Z2={z2_mod}, Om3={d_om3_mod}, rk3={rk3_mod}")
    print(f"  T:  Om2={d_om2_full}, Z2={z2_full}, Om3={d_om3_full}, rk3={rk3_full}")
    print(f"  New 2-paths: {sorted(new_2paths)}")
    print(f"  Lost 2-paths: {sorted(lost_2paths)}")
    print(f"  New 3-paths ({len(new_3paths)}): {sorted(new_3paths)[:5]}...")
    print(f"  Lost 3-paths ({len(lost_3paths)})")
    print(f"  delta Om2={d_om2_full - d_om2_mod}, delta Z2={z2_full - z2_mod}")
    print(f"  delta Om3={d_om3_full - d_om3_mod}, delta rk3={rk3_full - rk3_mod}")
    print(f"  beta2(T')={z2_mod - rk3_mod}, beta2(T)={z2_full - rk3_full}")

    # The new 2-paths always contain u->v as a subpath
    new_2_using_uv = [p for p in new_2paths if (p[0]==u and p[1]==v) or (p[1]==u and p[2]==v)]
    print(f"  New 2-paths using u->v: {sorted(new_2_using_uv)}")

    # Does restoring u->v create new Om3 elements?
    # New 3-paths must pass through u->v
    new_3_using_uv = [p for p in new_3paths
                      if (p[0]==u and p[1]==v) or (p[1]==u and p[2]==v) or (p[2]==u and p[3]==v)]
    print(f"  New 3-paths using u->v ({len(new_3_using_uv)}): {sorted(new_3_using_uv)[:8]}")

    analyzed += 1

print(f"\n\nAnalyzed {analyzed} cases")


# ============================================================
# PATTERN: What types of edges, when removed, create beta2>0?
# ============================================================
print(f"\n{'='*70}")
print("PATTERN ANALYSIS")
print("=" * 70)

by_type = Counter()
for r in removals_creating_beta2:
    by_type[(r['scores'], r['du'], r['dv'])] += 1

for (sc, du, dv), cnt in sorted(by_type.items()):
    print(f"  scores={sc}, remove {du}->{dv}: {cnt} cases")


# ============================================================
# Check: which score types NEVER have beta2>0 under removal?
# ============================================================
print(f"\n{'='*70}")
print("SCORES WITH NO beta2-CREATING REMOVALS")
print("=" * 70)

all_scores = set()
scores_with_removal = set()
for bits in range(total):
    A = build_adj(n, bits)
    sc = tuple(sorted([sum(row) for row in A]))
    all_scores.add(sc)

for r in removals_creating_beta2:
    scores_with_removal.add(r['scores'])

safe_scores = all_scores - scores_with_removal
for sc in sorted(safe_scores):
    print(f"  {sc}: NEVER creates beta2>0 under edge removal")
for sc in sorted(scores_with_removal):
    print(f"  {sc}: CAN create beta2>0 under edge removal")


# ============================================================
# CRITICAL: Edge-type analysis
# Which edges are "essential" for beta2=0?
# ============================================================
print(f"\n{'='*70}")
print("ESSENTIAL EDGE STRUCTURE")
print("=" * 70)

# For a specific tournament, find ALL essential edges
for bits in [0]:  # Transitive tournament
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    essential = []
    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            A_mod = [row[:] for row in A]
            A_mod[u][v] = 0
            betti = path_betti_numbers(A_mod, n, max_dim=2)
            if betti[2] > 0:
                essential.append((u, v, sum(A[u]), sum(A[v])))

    print(f"\nbits={bits}, scores={scores}")
    if essential:
        for u, v, du, dv in essential:
            print(f"  Essential edge: {u}->{v} (d_u={du}, d_v={dv})")
    else:
        print(f"  NO essential edges (all removals preserve beta2=0)")

# Check a few more tournament types
for bits in [341, 10, 100, 512]:
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    essential = []
    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            A_mod = [row[:] for row in A]
            A_mod[u][v] = 0
            betti = path_betti_numbers(A_mod, n, max_dim=2)
            if betti[2] > 0:
                essential.append((u, v, sum(A[u]), sum(A[v])))

    print(f"\nbits={bits}, scores={scores}")
    if essential:
        for u, v, du, dv in essential:
            print(f"  Essential: {u}->{v} (d_u={du}, d_v={dv})")
    else:
        print(f"  No essential edges")


print("\n\nDone.")
