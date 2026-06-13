#!/usr/bin/env python3
"""
beta2_flip_algebra.py - Algebraic analysis of beta2 under arc flip

When we flip arc u->v to v->u in tournament T:
- Some 2-paths are gained, some lost
- Some 3-paths are gained, some lost
- Omega dimensions change
- Z2 dimension changes
- beta2 = dim(Z2) - rk(bd3|_Om3) remains 0

Key quantities under flip:
  delta_Om2 = dim(Om2') - dim(Om2)
  delta_Z2 = dim(Z2') - dim(Z2)
  delta_rk3 = rk(bd3'|_Om3') - rk(bd3|_Om3)

beta2=0 invariance means: delta_Z2 = delta_rk3 always.

We know from HYP-228: delta|A2| = 2(d_u - d_v - 1)
We know from HYP-227: delta|A3| = (n-3) * delta|A2|

Question: What are delta_Om2, delta_Z2, delta_Om3 in terms of
local graph properties?

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


def compute_chain_data(A, n):
    """Compute all chain complex data."""
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

    dd = sum(1 for p in a3 if A[p[0]][p[2]] and A[p[1]][p[3]])

    return {
        'A2': len(a2), 'A3': len(a3),
        'Om2': d_om2, 'Om3': d_om3,
        'rk2': rk2, 'rk3': rk3,
        'Z2': z2, 'b2': b2, 'DD': dd,
        'b1': (n*(n-1)//2 - (n-1)) - rk2,  # dim(Z1) - rk(bd2)
    }


# ============================================================
# Exhaustive arc-flip analysis at n=5
# ============================================================
print("=" * 70)
print("ARC-FLIP ALGEBRA: n=5")
print("=" * 70)

n = 5
total = 1 << (n*(n-1)//2)

# For each flip, compute deltas
delta_data = []

for bits in range(total):
    A = build_adj(n, bits)
    d_before = compute_chain_data(A, n)

    for u in range(n):
        for v in range(n):
            if u >= v or not A[u][v]:
                continue

            # Flip u->v to v->u
            A_flip = [row[:] for row in A]
            A_flip[u][v] = 0
            A_flip[v][u] = 1

            d_after = compute_chain_data(A_flip, n)

            du = sum(A[u])
            dv = sum(A[v])

            delta = {
                'bits': bits, 'u': u, 'v': v,
                'du': du, 'dv': dv,
                'dA2': d_after['A2'] - d_before['A2'],
                'dA3': d_after['A3'] - d_before['A3'],
                'dOm2': d_after['Om2'] - d_before['Om2'],
                'dOm3': d_after['Om3'] - d_before['Om3'],
                'drk2': d_after['rk2'] - d_before['rk2'],
                'drk3': d_after['rk3'] - d_before['rk3'],
                'dZ2': d_after['Z2'] - d_before['Z2'],
                'dDD': d_after['DD'] - d_before['DD'],
                'db1': d_after['b1'] - d_before['b1'],
            }
            delta_data.append(delta)

print(f"Total flips analyzed: {len(delta_data)}")

# Check: delta_Z2 = delta_rk3 always?
violations = sum(1 for d in delta_data if d['dZ2'] != d['drk3'])
print(f"delta_Z2 = delta_rk3 violations: {violations}/{len(delta_data)}")

# Analyze the relationship between delta values and local properties
print(f"\nBy (d_u, d_v):")
by_deg = defaultdict(list)
for d in delta_data:
    by_deg[(d['du'], d['dv'])].append(d)

for (du, dv), entries in sorted(by_deg.items()):
    dA2_vals = set(d['dA2'] for d in entries)
    dOm2_vals = set(d['dOm2'] for d in entries)
    dZ2_vals = set(d['dZ2'] for d in entries)
    drk3_vals = set(d['drk3'] for d in entries)
    dOm3_vals = set(d['dOm3'] for d in entries)
    db1_vals = set(d['db1'] for d in entries)

    print(f"  d_u={du}, d_v={dv} ({len(entries)} flips):")
    print(f"    dA2={sorted(dA2_vals)}, dOm2={sorted(dOm2_vals)}, dZ2={sorted(dZ2_vals)}")
    print(f"    dOm3={sorted(dOm3_vals)}, drk3={sorted(drk3_vals)}")
    print(f"    db1={sorted(db1_vals)}")


# ============================================================
# Key identity: dOm2 = dZ2 + drk2
# And: dZ2 = drk3 (beta2=0 invariance)
# So: dOm2 = drk3 + drk2
# ============================================================
print(f"\n{'='*70}")
print("IDENTITY: dOm2 = drk2 + drk3")
print("=" * 70)

check_ok = True
for d in delta_data:
    if d['dOm2'] != d['drk2'] + d['drk3']:
        check_ok = False
        print(f"  FAIL: dOm2={d['dOm2']}, drk2={d['drk2']}, drk3={d['drk3']}")
        break

if check_ok:
    print("  CONFIRMED: dOm2 = drk2 + drk3 for all flips")
    print("  (This is trivially true from dOm2 = dZ2 + drk2 and dZ2 = drk3)")


# ============================================================
# Is drk2 locally determined?
# ============================================================
print(f"\n{'='*70}")
print("IS drk2 LOCALLY DETERMINED?")
print("=" * 70)

# drk2 = delta(rk(bd2|_Om2)) = delta(dim(Z1) - beta1) = -delta(beta1)
# Since Z1 = C(n-1,2) is constant for tournaments.
# So drk2 = -db1.

print("drk2 should = -db1:")
for d in delta_data[:5]:
    print(f"  drk2={d['drk2']}, -db1={-d['db1']}, match={d['drk2'] == -d['db1']}")

# drk2 = -db1 is NOT locally determined (depends on global tournament structure)
# But drk3 = dZ2 is what we need.

# Let's check if db1 is locally determined
print(f"\ndb1 by (d_u, d_v):")
for (du, dv), entries in sorted(by_deg.items()):
    db1_vals = Counter(d['db1'] for d in entries)
    print(f"  d_u={du}, d_v={dv}: db1 values = {dict(sorted(db1_vals.items()))}")


# ============================================================
# Check: is dZ2 = dOm2 + db1?
# ============================================================
print(f"\n{'='*70}")
print("dZ2 = dOm2 + db1?")
print("=" * 70)

# dZ2 = d(Om2 - rk2) = dOm2 - drk2 = dOm2 + db1
# Since drk2 = -db1.
check = all(d['dZ2'] == d['dOm2'] + d['db1'] for d in delta_data)
print(f"  dZ2 = dOm2 + db1: {check}")
print(f"  Equivalently: delta(dim Z2) = delta(dim Om2) + delta(beta1)")

# This means beta2=0 invariance <=> drk3 = dOm2 + db1
# <=> delta(rk bd3) = delta(dim Om2) + delta(beta1)
# <=> delta(rk bd3) = delta(dim Om2 - rk bd2)

print(f"\n  So beta2=0 invariance under flip means:")
print(f"  delta(rk(bd3|_Om3)) = delta(dim(Om2)) + delta(beta1)")
print(f"  = delta(dim(Om2)) - delta(rk(bd2|_Om2))")


print("\n\nDone.")
