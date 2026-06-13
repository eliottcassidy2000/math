#!/usr/bin/env python3
"""
beta2_injectivity_analysis.py - Why is d_3 injective on Omega_3 at surplus=0?

At n=5, surplus=0 tournaments have dim(Omega_3) = dim(Z_2) = 3 and d_3 is injective.
All surplus=0 tournaments have Omega_3 = DT paths only (3 DT paths).

Key question: CAN two distinct DT 4-paths have the same boundary?
If not, then d_3 restricted to DT is always injective.

Boundary of DT path (a,b,c,d):
  d_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)

where all four faces are in A_2 (since DT has both shortcuts).

For two DT paths p1 = (a1,b1,c1,d1) and p2 = (a2,b2,c2,d2):
  d_3(p1) = d_3(p2) iff they share all four faces.
  Since (b1,c1,d1) = (b2,c2,d2) and (a1,b1,c1) = (a2,b2,c2) etc.
  This forces a1=a2, b1=b2, c1=c2, d1=d2. So d_3 is injective on DT basis vectors.

But d_3 injectivity on span(DT) means: if sum alpha_i * p_i has d_3(sum)=0,
then alpha_i = 0 for all i. This needs the boundaries to be linearly independent,
not just pairwise distinct.

Let's check: in what proportion of tournaments are the DT boundaries
linearly independent?

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

def count_DT(A, n):
    """Count DT 4-paths: (a,b,c,d) with a->b->c->d, a->c, b->d."""
    dt_paths = []
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                if not A[a][c]: continue  # DT condition: a->c
                for d in range(n):
                    if d == a or d == b or d == c or not A[c][d]: continue
                    if not A[b][d]: continue  # DT condition: b->d
                    dt_paths.append((a,b,c,d))
    return dt_paths

def compute_all(A, n):
    allowed = {}
    for p in range(5):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break

    omega_basis = {}
    for p in range(4):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p >= 1 and p-1 in allowed else [])
        omega_basis[p] = basis

    dim_O2 = omega_basis[2].shape[1] if omega_basis[2].ndim == 2 else 0
    dim_O3 = omega_basis[3].shape[1] if omega_basis[3].ndim == 2 else 0

    if dim_O2 == 0:
        Z2 = 0
    else:
        bd2 = build_full_boundary_matrix(allowed[2], allowed.get(1, []))
        bd2_om = bd2 @ omega_basis[2]
        S_v = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = int(np.sum(np.abs(S_v) > 1e-8))
        Z2 = dim_O2 - rk2

    surplus = dim_O3 - Z2
    return {
        'allowed': allowed, 'omega_basis': omega_basis,
        'O2': dim_O2, 'O3': dim_O3, 'Z2': Z2, 'surplus': surplus
    }


n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"DT AND OMEGA_3 RELATIONSHIP AT n={n}")
print("=" * 70)

# For all tournaments, compare |DT| vs dim(Omega_3) vs |A_3|
dt_vs_omega = defaultdict(int)
dt_details = []

for bits in range(total):
    A = build_adj(n, bits)
    data = compute_all(A, n)
    dt_paths = count_DT(A, n)
    n_dt = len(dt_paths)
    n_A3 = len(data['allowed'].get(3, []))

    dt_vs_omega[(n_dt, data['O3'], data['surplus'])] += 1
    dt_details.append({
        'bits': bits, 'n_dt': n_dt, 'n_A3': n_A3,
        'O3': data['O3'], 'Z2': data['Z2'], 'surplus': data['surplus']
    })

print(f"\n  (|DT|, dim(O3), surplus) distribution:")
for key in sorted(dt_vs_omega.keys()):
    print(f"    DT={key[0]}, O3={key[1]}, surplus={key[2]}: {dt_vs_omega[key]}")

# Key question: does |DT| = dim(Omega_3) always?
dt_eq_O3 = sum(1 for d in dt_details if d['n_dt'] == d['O3'])
print(f"\n  |DT| = dim(O3): {dt_eq_O3}/{total} ({100*dt_eq_O3/total:.1f}%)")

# When they differ, what's the relationship?
differ = [(d['n_dt'], d['O3'], d['surplus']) for d in dt_details if d['n_dt'] != d['O3']]
if differ:
    differ_dist = Counter(differ)
    print(f"  Cases where |DT| != dim(O3):")
    for key in sorted(differ_dist.keys()):
        print(f"    DT={key[0]}, O3={key[1]}, surplus={key[2]}: {differ_dist[key]}")

# Also check if |DT| >= dim(Z_2) always
dt_ge_Z2 = sum(1 for d in dt_details if d['n_dt'] >= d['Z2'])
print(f"\n  |DT| >= dim(Z_2): {dt_ge_Z2}/{total} ({100*dt_ge_Z2/total:.1f}%)")

# And |DT| vs |A_3|
print(f"\n  |DT|/|A_3| ratios:")
ratio_dist = Counter(d['n_dt'] / d['n_A3'] if d['n_A3'] > 0 else -1 for d in dt_details)
for r in sorted(ratio_dist.keys()):
    if ratio_dist[r] >= 10:
        print(f"    {r:.4f}: {ratio_dist[r]}")

# Now do n=6
print(f"\n{'='*70}")
print(f"n = 6")
print(f"{'='*70}")

n = 6
n_arcs = n*(n-1)//2
total = 1 << n_arcs

dt_vs_omega6 = defaultdict(int)
dt_eq_O3_6 = 0
dt_ge_Z2_6 = 0

for bits in range(total):
    if bits % 10000 == 0 and bits > 0:
        print(f"  ... {bits}/{total}")

    A = build_adj(n, bits)
    data = compute_all(A, n)
    dt_paths = count_DT(A, n)
    n_dt = len(dt_paths)

    if n_dt == data['O3']:
        dt_eq_O3_6 += 1
    if n_dt >= data['Z2']:
        dt_ge_Z2_6 += 1

    dt_vs_omega6[(n_dt, data['O3'], data['surplus'])] += 1

print(f"\n  |DT| = dim(O3): {dt_eq_O3_6}/{total} ({100*dt_eq_O3_6/total:.1f}%)")
print(f"  |DT| >= dim(Z_2): {dt_ge_Z2_6}/{total} ({100*dt_ge_Z2_6/total:.1f}%)")

# Show the distribution for small DT values
print(f"\n  (|DT|, dim(O3), surplus) for DT<=20:")
for key in sorted(dt_vs_omega6.keys()):
    if key[0] <= 20:
        print(f"    DT={key[0]}, O3={key[1]}, surplus={key[2]}: {dt_vs_omega6[key]}")

print("\nDone.")
