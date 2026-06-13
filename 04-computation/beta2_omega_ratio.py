#!/usr/bin/env python3
"""
beta2_omega_ratio.py - Does the (n-3) ratio hold for dim(Omega_p)?

We proved: delta_|A_3| = (n-3) * delta_|A_2|

Question: Does delta_dim(Omega_3) = (n-3) * delta_dim(Omega_2)?
If so, this would give a path to proving beta_2 = 0 via arc-flip induction.

Also check: delta_dim(Z_2) vs delta_dim(Omega_2).

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os, random
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

def compute_data(A, n):
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

    dims = {p: omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
            for p in range(4)}

    # Z_2
    if dims[2] == 0:
        Z2 = 0
    else:
        bd2 = build_full_boundary_matrix(allowed[2], allowed.get(1, []))
        bd2_om = bd2 @ omega_basis[2]
        S_v = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = int(np.sum(np.abs(S_v) > 1e-8))
        Z2 = dims[2] - rk2

    # im(d_3)
    if dims[3] == 0 or 3 not in allowed:
        im_d3 = 0
    else:
        bd3 = build_full_boundary_matrix(allowed[3], allowed.get(2, []))
        bd3_om = bd3 @ omega_basis[3]
        S_v = np.linalg.svd(bd3_om, compute_uv=False)
        im_d3 = int(np.sum(np.abs(S_v) > 1e-8))

    return {
        'A2': len(allowed.get(2, [])), 'A3': len(allowed.get(3, [])),
        'O2': dims[2], 'O3': dims[3],
        'Z2': Z2, 'im_d3': im_d3,
        'beta2': Z2 - im_d3,
        'surplus': dims[3] - Z2
    }


# Exhaustive n=5
n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"OMEGA RATIO ANALYSIS AT n={n}")
print("=" * 70)

all_data = {}
for bits in range(total):
    A = build_adj(n, bits)
    all_data[bits] = compute_data(A, n)

# Analyze all flips
ratios = []
for bits in range(total):
    A = build_adj(n, bits)
    d_before = all_data[bits]

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

            d_after = all_data[bits_flip]

            dA2 = d_after['A2'] - d_before['A2']
            dA3 = d_after['A3'] - d_before['A3']
            dO2 = d_after['O2'] - d_before['O2']
            dO3 = d_after['O3'] - d_before['O3']
            dZ2 = d_after['Z2'] - d_before['Z2']
            d_im = d_after['im_d3'] - d_before['im_d3']

            d_u = sum(A[u])
            d_v = sum(A[v])

            ratios.append({
                'dA2': dA2, 'dA3': dA3, 'dO2': dO2, 'dO3': dO3,
                'dZ2': dZ2, 'd_im': d_im,
                'd_u': d_u, 'd_v': d_v,
                'surplus_before': d_before['surplus'],
                'surplus_after': d_after['surplus']
            })

# Check if dO3 = (n-3) * dO2
dO3_vs_dO2 = Counter((d['dO3'], d['dO2']) for d in ratios)
print(f"\n  Joint (dO3, dO2) distribution:")
for (dO3, dO2) in sorted(dO3_vs_dO2.keys()):
    is_ratio = (dO3 == (n-3) * dO2)
    marker = " <-- ratio" if is_ratio else ""
    print(f"    dO3={dO3:+d}, dO2={dO2:+d}: {dO3_vs_dO2[(dO3, dO2)]}{marker}")

# Check ratio when dO2 != 0
omega_ratio = Counter(d['dO3'] / d['dO2'] if d['dO2'] != 0 else 'inf'
                       for d in ratios)
print(f"\n  dO3/dO2 ratio when dO2 != 0:")
for r in sorted(k for k in omega_ratio if k != 'inf'):
    print(f"    {r}: {omega_ratio[r]}")
if 'inf' in omega_ratio:
    print(f"    dO2=0: {omega_ratio['inf']}")

# Check dZ2 vs dO2
dZ2_vs_dO2 = Counter((d['dZ2'], d['dO2']) for d in ratios)
print(f"\n  Joint (dZ2, dO2):")
for (dZ2, dO2) in sorted(dZ2_vs_dO2.keys()):
    print(f"    dZ2={dZ2:+d}, dO2={dO2:+d}: {dZ2_vs_dO2[(dZ2, dO2)]}")

# What about dO3 as a function of (d_u, d_v)?
print(f"\n  dO3 by (d_u, d_v):")
by_deg = defaultdict(set)
for d in ratios:
    by_deg[(d['d_u'], d['d_v'])].add(d['dO3'])
for key in sorted(by_deg.keys()):
    vals = sorted(by_deg[key])
    print(f"    d_u={key[0]}, d_v={key[1]}: dO3 in {vals}")

# dO2 by (d_u, d_v)
print(f"\n  dO2 by (d_u, d_v):")
by_deg2 = defaultdict(set)
for d in ratios:
    by_deg2[(d['d_u'], d['d_v'])].add(d['dO2'])
for key in sorted(by_deg2.keys()):
    vals = sorted(by_deg2[key])
    print(f"    d_u={key[0]}, d_v={key[1]}: dO2 in {vals}")

# Is dO2 determined by (d_u, d_v)?
print(f"\n  Is dO2 determined by (d_u, d_v)?")
determined = all(len(v) == 1 for v in by_deg2.values())
print(f"    {determined}")

# Check what determines dO2 and dO3
# Try: dO2 function of (d_u - d_v)
print(f"\n  dO2 by (d_u - d_v):")
by_diff = defaultdict(set)
for d in ratios:
    by_diff[d['d_u'] - d['d_v']].add(d['dO2'])
for diff in sorted(by_diff.keys()):
    vals = sorted(by_diff[diff])
    print(f"    d_u-d_v={diff}: dO2 in {vals}")

print("\nDone.")
