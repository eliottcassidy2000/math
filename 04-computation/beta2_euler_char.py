#!/usr/bin/env python3
"""
beta2_euler_char.py — Euler characteristic and rank analysis

If we can express rk(∂₂) and rk(∂₃) in terms of tournament invariants,
then β₂ = dim(Ω₂) - rk(∂₂) - rk(∂₃) gives a formula.

For tournaments on n vertices:
- dim(Ω₀) = n
- dim(Ω₁) = C(n,2) (all edges)
- rk(∂₁) = n - 1 (connected graph)
- dim(Ω₂), rk(∂₂), rk(∂₃) depend on tournament structure

Key hypothesis: rk(∂₂) = C(n,2) - n + 1 - δ where δ depends on c₃
And: Z₂ = dim(Ω₂) - rk(∂₂) = β₃ + rk(∂₃) (from β₂=0)

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def full_analysis(A, n):
    """Compute all Ω dimensions and boundary ranks."""
    result = {}
    oms = {}
    aps = {}
    
    for p in range(n):
        aps[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            oms[p] = np.eye(n)
        elif aps[p]:
            oms[p] = compute_omega_basis(A, n, p, aps[p], aps[p-1])
        else:
            oms[p] = np.zeros((0,0))
        result[f'dim_Om{p}'] = oms[p].shape[1] if oms[p].ndim == 2 and oms[p].shape[0] > 0 else 0
    
    for p in range(1, n):
        dp = result[f'dim_Om{p}']
        dpm1 = result[f'dim_Om{p-1}']
        if dp == 0 or dpm1 == 0:
            result[f'rk_bd{p}'] = 0
            continue
        bd = build_full_boundary_matrix(aps[p], aps[p-1])
        bd_om = bd @ oms[p]
        coords = np.linalg.lstsq(oms[p-1], bd_om, rcond=None)[0]
        result[f'rk_bd{p}'] = np.linalg.matrix_rank(coords, tol=1e-8)
    
    # Betti numbers
    for p in range(n):
        dp = result.get(f'dim_Om{p}', 0)
        rkp = result.get(f'rk_bd{p}', 0)
        rkp1 = result.get(f'rk_bd{p+1}', 0)
        result[f'beta{p}'] = dp - rkp - rkp1
    
    # 3-cycles
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[i][k] and A[k][j] and A[j][i]):
                    c3 += 1
    result['c3'] = c3
    
    # Score sequence
    result['scores'] = tuple(sorted(sum(A[i]) for i in range(n)))
    
    return result

# Analyze n=5
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

print("=" * 70)
print(f"COMPLETE RANK ANALYSIS AT n={n}")
print("=" * 70)

data = []
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    r = full_analysis(A, n)
    r['bits'] = bits
    data.append(r)

# Group by (c₃, dim Ω₂, rk∂₂, dim Ω₃, rk∂₃)
groups = Counter()
for d in data:
    key = (d['c3'], d['dim_Om2'], d['rk_bd2'], d['dim_Om3'], d['rk_bd3'])
    groups[key] += 1

print(f"\n(c₃, Ω₂, rk∂₂, Ω₃, rk∂₃): count  [Z₂ | B₂ | β₂ | β₃ | surplus]")
for key, cnt in sorted(groups.items()):
    c3, o2, r2, o3, r3 = key
    z2 = o2 - r2
    b2 = r3
    b2_actual = r3
    beta2 = z2 - b2
    beta3 = o3 - r3 - (0 if n < 5 else 0)  # rk(∂₄) for n=5
    # Need rk(∂₄)
    surplus = o3 - z2
    print(f"  c₃={c3}, Ω₂={o2}, rk∂₂={r2}, Ω₃={o3}, rk∂₃={r3}: {cnt}  "
          f"[Z₂={z2} | B₂={b2} | β₂={beta2} | ker∂₃={o3-r3} | surplus={surplus}]")

# Check the Euler characteristic
print(f"\n{'='*70}")
print("EULER CHARACTERISTIC")
print("=" * 70)

euler_dist = Counter()
for d in data:
    chi = sum((-1)**p * d.get(f'dim_Om{p}', 0) for p in range(n))
    beta_sum = sum((-1)**p * d.get(f'beta{p}', 0) for p in range(n))
    euler_dist[(chi, beta_sum)] += 1

print(f"(χ from dims, χ from betti): count")
for key, cnt in sorted(euler_dist.items()):
    print(f"  χ_dim={key[0]}, χ_betti={key[1]}: {cnt}")

# Formulas
print(f"\n{'='*70}")
print("FORMULA SEARCH")
print("=" * 70)

# Hypothesis: dim Ω₂ = C(n,3) - c₃ + f(c₃, n)
# At n=5: C(5,3) = 10
# Check dim Ω₂ - (10 - c₃)
print("\ndim Ω₂ - (C(n,3) - c₃) = ?")
delta_o2 = Counter()
for d in data:
    delta = d['dim_Om2'] - (10 - d['c3'])
    delta_o2[(d['c3'], delta)] += 1

for key, cnt in sorted(delta_o2.items()):
    print(f"  c₃={key[0]}: Ω₂ - TT_count = {key[1]} ({cnt} tournaments)")

# Hypothesis: rk(∂₂) = C(n,2) - n + 1 - f(c₃)
# C(5,2) - 4 = 6
print(f"\nrk(∂₂) as function of c₃:")
rk_by_c3 = {}
for d in data:
    c3 = d['c3']
    rk = d['rk_bd2']
    if c3 not in rk_by_c3:
        rk_by_c3[c3] = Counter()
    rk_by_c3[c3][rk] += 1

for c3 in sorted(rk_by_c3):
    print(f"  c₃={c3}: rk(∂₂) = {dict(sorted(rk_by_c3[c3].items()))}")

# Hypothesis: dim Ω₃ formula
print(f"\ndim Ω₃ as function of c₃ and score:")
o3_by_c3_score = {}
for d in data:
    key = (d['c3'], d['scores'])
    o3 = d['dim_Om3']
    if key not in o3_by_c3_score:
        o3_by_c3_score[key] = Counter()
    o3_by_c3_score[key][o3] += 1

for key in sorted(o3_by_c3_score):
    c3, scores = key
    print(f"  c₃={c3}, scores={scores}: Ω₃ = {dict(sorted(o3_by_c3_score[key].items()))}")

# KEY: Z₂ = rk(∂₃) always. What's the formula for rk(∂₃)?
print(f"\nrk(∂₃) = Z₂ as function of c₃ and score:")
rk3_by = {}
for d in data:
    key = (d['c3'], d['scores'])
    val = d['rk_bd3']
    if key not in rk3_by:
        rk3_by[key] = Counter()
    rk3_by[key][val] += 1

for key in sorted(rk3_by):
    c3, scores = key
    print(f"  c₃={c3}, scores={scores}: rk(∂₃) = {dict(sorted(rk3_by[key].items()))}")

# Also compute dim Ω₄ for full picture
print(f"\n{'='*70}")
print("FULL CHAIN COMPLEX AT n=5")
print("=" * 70)

full_groups = Counter()
for d in data:
    dims = tuple(d.get(f'dim_Om{p}', 0) for p in range(5))
    ranks = tuple(d.get(f'rk_bd{p}', 0) for p in range(1, 5))
    bettis = tuple(d.get(f'beta{p}', 0) for p in range(5))
    full_groups[(dims, ranks, bettis)] += 1

print(f"\n(Ω₀..Ω₄, rk∂₁..rk∂₄, β₀..β₄): count")
for key, cnt in sorted(full_groups.items()):
    dims, ranks, bettis = key
    print(f"  dims={list(dims)}, ranks={list(ranks)}, β={list(bettis)}: {cnt}")

print("\nDone.")
