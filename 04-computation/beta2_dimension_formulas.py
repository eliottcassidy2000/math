#!/usr/bin/env python3
"""Find universal formulas for dim(Ω_2), dim(Ω_3), dim(Z_2), dim(ker ∂_3).

Goal: Express these as functions of tournament invariants (n, t3, score seq, etc.)
and show algebraically that dim(Z_2) = rank(∂_3).

Also: analyze the STRUCTURE of ker(∂_3) — what 3-chains are boundaries that map to 0?
"""
import numpy as np
from itertools import combinations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    path_betti_numbers
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def tournament_invariants(A, n):
    """Compute key tournament invariants."""
    scores = sorted([sum(A[i]) for i in range(n)])
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: t3 += 1
                if A[j][i] and A[i][k] and A[k][j]: t3 += 1

    # Transitive triples
    tt = 0
    for a in range(n):
        for b in range(n):
            if b == a: continue
            if not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    tt += 1

    # Doubly-transitive 4-paths (DT)
    dt = 0
    for a in range(n):
        for b in range(n):
            if b == a or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b or not A[b][c]: continue
                for d in range(n):
                    if d == a or d == b or d == c or not A[c][d]: continue
                    if A[a][c] and A[b][d]:
                        dt += 1

    # Count 3-cycles (ordered)
    c3 = 0
    for a in range(n):
        for b in range(n):
            if b == a: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[a][b] and A[b][c] and A[c][a]:
                    c3 += 1

    return {
        'scores': tuple(scores),
        't3': t3,
        'tt': tt,
        'dt': dt,
        'c3': c3
    }

def full_chain_analysis(A, n):
    """Compute all chain complex dimensions and ranks."""
    a = {}
    for p in range(6):
        a[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

    om = {}
    for p in range(6):
        if p == 0:
            om[p] = np.eye(n)
        else:
            om[p] = compute_omega_basis(A, n, p, a[p], a[p-1])

    dims = {}
    ranks = {}
    kers = {}

    for p in range(6):
        if om[p].ndim == 2:
            dims[p] = om[p].shape[1]
        else:
            dims[p] = 0

    for p in range(1, 5):
        if dims[p] == 0:
            ranks[p] = 0
            kers[p] = 0
            continue
        bd = build_full_boundary_matrix(a[p], a[p-1])
        bd_om = bd @ om[p]
        if dims[p-1] > 0 and p > 0:
            # Express in Ω_{p-1} coords
            im_coords, _, _, _ = np.linalg.lstsq(om[p-1], bd_om, rcond=None)
            S = np.linalg.svd(im_coords, compute_uv=False)
        else:
            S = np.linalg.svd(bd_om, compute_uv=False)
        ranks[p] = sum(s > 1e-8 for s in S)
        kers[p] = dims[p] - ranks[p]

    return dims, ranks, kers, a, om

# ===== Analyze at n=5 =====
print("=" * 70)
print("DIMENSION FORMULAS FOR β_2 = 0")
print("=" * 70)

for n in [4, 5]:
    print(f"\n{'='*50}")
    print(f"n = {n}")
    print(f"{'='*50}")
    print(f"{'scores':>20} t3  tt   dt   c3  |  Ω2  Ω3  Ω4  |  Z2  B2  kerδ3 | β1  β2  β3")
    print("-" * 100)

    seen = set()
    for A in all_tournaments_gen(n):
        inv = tournament_invariants(A, n)
        dims, ranks, kers, a, om = full_chain_analysis(A, n)

        z2 = kers.get(2, 0)  # ker(∂_2)
        b2 = ranks.get(3, 0)  # im(∂_3) = rank(∂_3)
        beta1 = kers.get(1, 0) - ranks.get(2, 0) if 1 in kers and 2 in ranks else 0
        beta2 = z2 - b2
        beta3 = kers.get(3, 0) - ranks.get(4, 0) if 3 in kers and 4 in ranks else kers.get(3, 0)
        ker_d3 = dims.get(3, 0) - b2

        key = (inv['scores'], inv['t3'])
        if key in seen:
            continue
        seen.add(key)

        print(f"{str(inv['scores']):>20} {inv['t3']:3d} {inv['tt']:4d} {inv['dt']:4d} {inv['c3']:4d}  | "
              f"{dims.get(2,0):3d} {dims.get(3,0):3d} {dims.get(4,0):3d}  | "
              f"{z2:3d} {b2:3d} {ker_d3:5d}   | {beta1:2d}  {beta2:2d}  {beta3:2d}")

# ===== Look for formulas =====
print(f"\n\n{'='*70}")
print("FORMULA SEARCH")
print("=" * 70)

n = 5
data_rows = []
for A in all_tournaments_gen(n):
    inv = tournament_invariants(A, n)
    dims, ranks, kers, a, om = full_chain_analysis(A, n)
    z2 = kers.get(2, 0)
    b2 = ranks.get(3, 0)
    ker_d3 = dims.get(3, 0) - b2
    beta1 = max(0, kers.get(1, 0) - ranks.get(2, 0))
    data_rows.append({
        't3': inv['t3'],
        'tt': inv['tt'],
        'dt': inv['dt'],
        'c3': inv['c3'],
        'scores': inv['scores'],
        'dim_om1': dims.get(1, 0),
        'dim_om2': dims.get(2, 0),
        'dim_om3': dims.get(3, 0),
        'dim_om4': dims.get(4, 0),
        'z2': z2,
        'b2': b2,
        'ker_d3': ker_d3,
        'rank_d2': ranks.get(2, 0),
        'beta1': beta1
    })

# Check universal formulas
print(f"\nn={n}: Checking candidate formulas...")

# dim(Ω_1) = C(n,2) always for tournaments
om1_vals = set(d['dim_om1'] for d in data_rows)
print(f"  dim(Ω_1) values: {sorted(om1_vals)} (expect {n*(n-1)//2})")

# rank(∂_2) = ?
rd2_vals = set(d['rank_d2'] for d in data_rows)
print(f"  rank(∂_2) values: {sorted(rd2_vals)}")

# Check: rank(∂_2) = C(n,2) - n + 1 - β_1 ?
for d in data_rows:
    expected = n*(n-1)//2 - n + 1 - d['beta1']
    if d['rank_d2'] != expected:
        print(f"  MISMATCH: rank_d2={d['rank_d2']}, expected={expected}")
        break
else:
    print(f"  rank(∂_2) = C(n,2) - n + 1 - β_1: CONFIRMED")

# dim(Ω_2) vs tt (transitive triples)
print(f"\n  dim(Ω_2) vs tt:")
for d in data_rows[:20]:
    diff = d['dim_om2'] - d['tt']
    #print(f"    scores={d['scores']}, tt={d['tt']}, Ω_2={d['dim_om2']}, diff={diff}")

# Check: is dim(Ω_2) = tt + f(t3) for some function?
tt_to_om2 = defaultdict(set)
t3_to_om2 = defaultdict(set)
for d in data_rows:
    tt_to_om2[d['tt']].add(d['dim_om2'])
    t3_to_om2[d['t3']].add(d['dim_om2'])

print(f"  t3 → Ω_2: {dict((k, sorted(v)) for k, v in sorted(t3_to_om2.items()))}")

# Z_2 = ker(∂_2) = dim(Ω_2) - rank(∂_2)
# Is Z_2 a function of t3?
t3_to_z2 = defaultdict(set)
for d in data_rows:
    t3_to_z2[d['t3']].add(d['z2'])
print(f"  t3 → Z_2: {dict((k, sorted(v)) for k, v in sorted(t3_to_z2.items()))}")

# Is Z_2 a function of (t3, scores)?
key_to_z2 = defaultdict(set)
for d in data_rows:
    key_to_z2[(d['scores'], d['t3'])].add(d['z2'])
print(f"  (scores, t3) → Z_2: all unique? {all(len(v)==1 for v in key_to_z2.values())}")

# dim(Ω_3) as function of invariants
t3_to_om3 = defaultdict(set)
for d in data_rows:
    t3_to_om3[d['t3']].add(d['dim_om3'])
print(f"  t3 → Ω_3: {dict((k, sorted(v)) for k, v in sorted(t3_to_om3.items()))}")

# dt vs dim(Ω_3)
dt_to_om3 = defaultdict(set)
for d in data_rows:
    dt_to_om3[d['dt']].add(d['dim_om3'])
print(f"  dt → Ω_3: {dict((k, sorted(v)) for k, v in sorted(dt_to_om3.items()))}")

# Is Ω_3 determined by dt? Or by dt + something?
print(f"\n  Key relationship: Ω_3 - dt = ?")
for d in data_rows[:30]:
    if d['dt'] > 0 or d['dim_om3'] > 0:
        gap = d['dim_om3'] - d['dt']
        #print(f"    scores={d['scores']}, dt={d['dt']}, Ω_3={d['dim_om3']}, gap={gap}")

# Collect (dt, Ω_3) pairs
dt_om3_pairs = Counter()
for d in data_rows:
    dt_om3_pairs[(d['dt'], d['dim_om3'])] += 1
print(f"  (dt, Ω_3) distribution:")
for (dt, om3), count in sorted(dt_om3_pairs.items()):
    print(f"    dt={dt:3d}, Ω_3={om3:3d}, gap={om3-dt:3d}: ×{count}")

# ===== Euler characteristic =====
print(f"\n  Euler characteristic analysis:")
chi_dist = Counter()
for d in data_rows:
    chi = sum((-1)**p * d[f'dim_om{p}'] for p in [1,2,3,4] if f'dim_om{p}' in d)
    chi += 1  # dim_om0 = n but β_0 = 1...
    # Actually χ = Σ(-1)^p dim(Ω_p) = 1 - β_1 + 0 - β_3 + ...
    chi = n - d['dim_om1'] + d['dim_om2'] - d['dim_om3'] + d['dim_om4']
    chi_dist[chi] += 1
print(f"  χ distribution: {dict(sorted(chi_dist.items()))}")

# Check: χ = 1 - β_1?
for d in data_rows:
    chi = n - d['dim_om1'] + d['dim_om2'] - d['dim_om3'] + d['dim_om4']
    expected = 1 - d['beta1']
    if chi != expected:
        print(f"  χ ≠ 1-β_1: χ={chi}, β_1={d['beta1']}")
        break
else:
    print(f"  χ = 1 - β_1: CONFIRMED (consistent with β_2=β_3=0)")

# ===== Key formula: dim(Z_2) in terms of t3, tt, c3 =====
print(f"\n  Formula hunting for Z_2:")
for d in data_rows:
    # Test: Z_2 = t3?
    pass

t3_z2_check = all(d['z2'] == d['t3'] for d in data_rows)
print(f"  Z_2 = t3? {t3_z2_check}")

if not t3_z2_check:
    # Try other formulas
    for d in data_rows[:5]:
        print(f"    t3={d['t3']}, Z_2={d['z2']}, Ω_2={d['dim_om2']}, rank_d2={d['rank_d2']}")

    # Z_2 = Ω_2 - rank_d2 = Ω_2 - (C(n,2) - n + 1 - β_1)
    # So Z_2 = Ω_2 - C(n,2) + n - 1 + β_1
    print(f"\n  Z_2 = Ω_2 - C(n,2) + n - 1 + β_1:")
    for d in data_rows[:5]:
        cn2 = n*(n-1)//2
        computed = d['dim_om2'] - cn2 + n - 1 + d['beta1']
        print(f"    Ω_2={d['dim_om2']}, β_1={d['beta1']}, Z_2={d['z2']}, formula={computed}")

# ===== β_3 analysis =====
print(f"\n  β_3 values:")
b3_dist = Counter()
for d in data_rows:
    # β_3 = ker(∂_3) - rank(∂_4)
    # We need rank(∂_4) too... let me just check β_3 from betti
    b3_dist[d['ker_d3']] += 1
print(f"  ker(∂_3) distribution: {dict(sorted(b3_dist.items()))}")

print("\nDone.")
