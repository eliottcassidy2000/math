#!/usr/bin/env python3
"""
beta2_arcflip_proof.py - Analyze arc-flip mechanism for beta_2=0 preservation

Strategy: prove beta_2=0 for ALL tournaments via arc-flip induction.
Base: transitive tournament (beta_2=0 trivially).
Step: flipping any arc preserves beta_2=0.

Key quantities to track under arc flip (u,v) -> (v,u):
- dim(Omega_2) and dim(Z_2) = ker(d_2) in Omega_2
- dim(Omega_3) and rk(d_3) = image of d_3 in Z_2
- surplus = dim(Omega_3) - dim(Z_2) (must stay >= 0)

The change in allowed paths: paths using u->v are lost, paths using v->u are gained.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from itertools import combinations
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

def bits_from_adj(A, n):
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 1:
                bits |= (1 << idx)
            idx += 1
    return bits

def flip_arc(A, u, v, n):
    """Return new adj matrix with arc (u,v) flipped."""
    B = [row[:] for row in A]
    B[u][v] = 1 - B[u][v]
    B[v][u] = 1 - B[v][u]
    return B

def compute_chain_data(A, n, max_dim=4):
    """Compute detailed chain complex data."""
    allowed = {}
    for p in range(max_dim + 2):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break

    result = {'allowed': allowed}

    # Compute Omega bases
    omega_basis = {}
    for p in range(max_dim + 1):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p >= 1 and p-1 in allowed else [])
        omega_basis[p] = basis
    result['omega_dim'] = {p: omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
                           for p in range(max_dim + 1)}

    # Compute boundary matrices restricted to Omega
    bd_omega = {}
    for p in range(1, max_dim + 1):
        dim_p = result['omega_dim'][p]
        if dim_p == 0:
            continue
        bd = build_full_boundary_matrix(allowed[p], allowed[p-1] if p-1 in allowed else [])
        bd_omega[p] = bd @ omega_basis[p]

    # Compute ranks and kernels
    rk = {}
    ker = {}
    for p in range(1, max_dim + 1):
        if p not in bd_omega:
            rk[p] = 0
            ker[p] = result['omega_dim'][p]
        else:
            S_v = np.linalg.svd(bd_omega[p], compute_uv=False)
            rk[p] = int(np.sum(np.abs(S_v) > 1e-8))
            ker[p] = result['omega_dim'][p] - rk[p]

    result['rk'] = rk
    result['ker'] = ker

    # Z_2 = ker(d_2) in Omega_2
    result['Z2'] = ker.get(2, 0)
    # im(d_3) = rk(d_3|Omega_3)
    result['im_d3'] = rk.get(3, 0)
    # beta_2 = Z_2 - im(d_3)
    result['beta2'] = result['Z2'] - result['im_d3']
    # surplus
    result['surplus'] = result['omega_dim'].get(3, 0) - result['Z2']

    return result


print("=" * 70)
print("ARC-FLIP PRESERVATION OF BETA_2 = 0")
print("=" * 70)

for n in [5, 6]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    n_arcs = n * (n-1) // 2
    total = 1 << n_arcs

    # Track all flip effects
    delta_surplus = []  # (surplus_before, surplus_after, delta)
    min_surplus_after = float('inf')
    max_delta_negative = 0
    counterexamples = []

    count = 0
    for bits in range(total):
        if n == 6 and count % 5000 == 0 and count > 0:
            print(f"  ... {count}/{total}")
        count += 1

        A = build_adj(n, bits)
        data = compute_chain_data(A, n)

        if data['beta2'] != 0:
            print(f"  WARNING: beta_2 = {data['beta2']} for bits={bits}")
            continue

        # Try flipping each arc
        for u in range(n):
            for v in range(n):
                if u == v or A[u][v] == 0:
                    continue
                B = flip_arc(A, u, v, n)
                data_flip = compute_chain_data(B, n)

                s_before = data['surplus']
                s_after = data_flip['surplus']
                delta = s_after - s_before

                delta_surplus.append((s_before, s_after, delta))
                min_surplus_after = min(min_surplus_after, s_after)

                if data_flip['beta2'] != 0:
                    counterexamples.append((bits, u, v, data, data_flip))

    print(f"\n  Total flips analyzed: {len(delta_surplus)}")
    print(f"  Min surplus after flip: {min_surplus_after}")
    print(f"  Beta_2 preservation: {'YES' if not counterexamples else 'NO'}")

    if counterexamples:
        print(f"  COUNTEREXAMPLES: {len(counterexamples)}")
        for bits, u, v, d1, d2 in counterexamples[:3]:
            print(f"    bits={bits}, flip ({u},{v}): surplus {d1['surplus']}->{d2['surplus']}, beta2={d2['beta2']}")

    # Analyze delta distribution
    deltas = [d for _, _, d in delta_surplus]
    from collections import Counter
    delta_dist = Counter(deltas)
    print(f"\n  Delta surplus distribution:")
    for d in sorted(delta_dist.keys()):
        print(f"    delta={d}: {delta_dist[d]} ({100*delta_dist[d]/len(deltas):.1f}%)")

    # Key question: when surplus drops, what compensates?
    drops = [(sb, sa, d) for sb, sa, d in delta_surplus if d < 0]
    print(f"\n  Surplus drops: {len(drops)} ({100*len(drops)/len(delta_surplus):.1f}%)")
    if drops:
        max_drop = min(d for _, _, d in drops)
        print(f"  Max drop: {max_drop}")
        # Is drop always from surplus >= |drop|?
        safe = all(sb >= abs(d) for sb, sa, d in drops)
        print(f"  Always safe (surplus >= |drop|): {safe}")

    # Analyze by starting surplus
    surplus_groups = {}
    for sb, sa, d in delta_surplus:
        if sb not in surplus_groups:
            surplus_groups[sb] = []
        surplus_groups[sb].append((sa, d))

    print(f"\n  By starting surplus:")
    for sb in sorted(surplus_groups.keys()):
        group = surplus_groups[sb]
        min_sa = min(sa for sa, _ in group)
        max_drop = min(d for _, d in group) if any(d < 0 for _, d in group) else 0
        print(f"    surplus={sb}: {len(group)} flips, min_after={min_sa}, max_drop={max_drop}")

print("\nDone.")
