#!/usr/bin/env python3
"""
beta2_surplus_zero_stability.py - Why surplus=0 tournaments are stable under arc flips

At n=5, surplus=0 means dim(Omega_3) = dim(Z_2) exactly.
Key observation: max_drop=0 when starting from surplus=0.
This means flipping ANY arc either maintains or increases surplus.

Question: WHY? What structural property of surplus=0 tournaments
forces this stability?

Hypothesis: surplus=0 iff transitive triangle count = certain value,
and flipping from this state always creates MORE Omega_3 elements
relative to Z_2 elements.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
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

def count_3cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(n):
            if i == j: continue
            for k in range(n):
                if k == i or k == j: continue
                if A[i][j] and A[j][k] and A[k][i]:
                    count += 1
    return count // 3  # Each 3-cycle counted 3 times

def count_transitive_triples(A, n):
    count = 0
    for i in range(n):
        for j in range(n):
            if i == j: continue
            for k in range(n):
                if k == i or k == j: continue
                # i->j->k and i->k (transitive)
                if A[i][j] and A[j][k] and A[i][k]:
                    count += 1
    return count  # Each transitive triple counted once per ordering

def score_seq(A, n):
    return tuple(sorted(sum(row) for row in A))

def compute_surplus(A, n):
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

    # Z_2 = ker(d_2|Omega_2)
    if dim_O2 == 0:
        Z2 = 0
    elif 2 not in allowed or not allowed[2]:
        Z2 = 0
    else:
        bd2 = build_full_boundary_matrix(allowed[2], allowed[1] if 1 in allowed else [])
        bd2_om = bd2 @ omega_basis[2]
        S_v = np.linalg.svd(bd2_om, compute_uv=False)
        rk2 = int(np.sum(np.abs(S_v) > 1e-8))
        Z2 = dim_O2 - rk2

    surplus = dim_O3 - Z2
    return surplus, dim_O2, dim_O3, Z2


n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"SURPLUS=0 STABILITY ANALYSIS AT n={n}")
print("=" * 70)

# Classify all tournaments
surplus_data = {}
for bits in range(total):
    A = build_adj(n, bits)
    surplus, dim_O2, dim_O3, Z2 = compute_surplus(A, n)
    t3 = count_3cycles(A, n)
    tt = count_transitive_triples(A, n)
    score = score_seq(A, n)
    surplus_data[bits] = {
        'surplus': surplus, 'dim_O2': dim_O2, 'dim_O3': dim_O3,
        'Z2': Z2, 't3': t3, 'tt': tt, 'score': score
    }

# Group by surplus
from collections import Counter, defaultdict
surplus_groups = defaultdict(list)
for bits, data in surplus_data.items():
    surplus_groups[data['surplus']].append(bits)

print(f"\n  Surplus distribution:")
for s in sorted(surplus_groups.keys()):
    group = surplus_groups[s]
    t3_vals = Counter(surplus_data[b]['t3'] for b in group)
    score_vals = Counter(surplus_data[b]['score'] for b in group)
    print(f"  surplus={s}: {len(group)} tournaments")
    print(f"    t3 values: {dict(sorted(t3_vals.items()))}")
    print(f"    scores: {dict(sorted(score_vals.items()))}")
    if s == 0:
        for b in group[:5]:
            d = surplus_data[b]
            print(f"    Example bits={b}: O2={d['dim_O2']}, O3={d['dim_O3']}, Z2={d['Z2']}, t3={d['t3']}, TT={d['tt']}, score={d['score']}")

# For surplus=0: analyze what happens under flip
print(f"\n{'='*70}")
print(f"SURPLUS=0 FLIP ANALYSIS")
print(f"{'='*70}")

zero_bits = surplus_groups[0]
print(f"  Total surplus=0 tournaments: {len(zero_bits)}")

flip_effects = []  # (bits, u, v, surplus_before, surplus_after, delta_O3, delta_Z2)

for bits in zero_bits:
    A = build_adj(n, bits)
    d_before = surplus_data[bits]
    for u in range(n):
        for v in range(n):
            if u == v or A[u][v] == 0:
                continue
            # Flip u->v to v->u
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

            d_after = surplus_data[bits_flip]
            delta_O3 = d_after['dim_O3'] - d_before['dim_O3']
            delta_Z2 = d_after['Z2'] - d_before['Z2']
            flip_effects.append((bits, u, v,
                                d_before['surplus'], d_after['surplus'],
                                delta_O3, delta_Z2))

# Analyze
print(f"  Total flips from surplus=0: {len(flip_effects)}")

# surplus_after distribution
sa_dist = Counter(sa for _, _, _, sb, sa, _, _ in flip_effects)
print(f"  Surplus after flip: {dict(sorted(sa_dist.items()))}")

# delta_O3 and delta_Z2 distribution
dO3_dist = Counter(dO3 for _, _, _, _, _, dO3, _ in flip_effects)
dZ2_dist = Counter(dZ2 for _, _, _, _, _, _, dZ2 in flip_effects)
print(f"  Delta Omega_3: {dict(sorted(dO3_dist.items()))}")
print(f"  Delta Z_2:     {dict(sorted(dZ2_dist.items()))}")

# KEY: for surplus=0, we need delta_O3 >= delta_Z2 always
violations = [(bits, u, v, sb, sa, dO3, dZ2)
              for bits, u, v, sb, sa, dO3, dZ2 in flip_effects
              if dO3 < dZ2]
print(f"\n  Violations (delta_O3 < delta_Z2): {len(violations)}")
if violations:
    for bits, u, v, sb, sa, dO3, dZ2 in violations[:5]:
        print(f"    bits={bits}, flip ({u},{v}): dO3={dO3}, dZ2={dZ2}")

# When does Z2 increase?
z2_increase = [(bits, u, v, sb, sa, dO3, dZ2)
               for bits, u, v, sb, sa, dO3, dZ2 in flip_effects
               if dZ2 > 0]
print(f"\n  Z_2 increases: {len(z2_increase)} ({100*len(z2_increase)/len(flip_effects):.1f}%)")
if z2_increase:
    # What happens to O3 when Z2 increases?
    for bits, u, v, sb, sa, dO3, dZ2 in z2_increase[:10]:
        d = surplus_data[bits]
        print(f"    bits={bits}: t3={d['t3']}, flip ({u},{v}): dZ2=+{dZ2}, dO3={dO3:+d}")

# Joint distribution of (dO3, dZ2)
joint = Counter((dO3, dZ2) for _, _, _, _, _, dO3, dZ2 in flip_effects)
print(f"\n  Joint (delta_O3, delta_Z2):")
for (dO3, dZ2) in sorted(joint.keys()):
    print(f"    ({dO3:+d}, {dZ2:+d}): {joint[(dO3, dZ2)]}")

print("\nDone.")
