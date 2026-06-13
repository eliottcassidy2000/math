#!/usr/bin/env python3
"""
trh_d2_cycle_predictor.py — opus-2026-03-13-S71

Find a combinatorial predicate for TRH d²=0.

Key observation: at n=5, t3=2 ⟺ d²≠0.
At n=6, (t3, c4) almost works but not perfectly.

Try: use ALL cycle counts (c3, c4, c5, c6) as predictor.
Or: use the full score+t3+c4 triple.
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict

def get_regular_paths(A, m):
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

def interior_boundary(paths_m, paths_m1):
    if not paths_m or not paths_m1:
        return np.zeros((len(paths_m1) if paths_m1 else 0,
                         len(paths_m) if paths_m else 0), dtype=int)
    idx = {p: i for i, p in enumerate(paths_m1)}
    m = len(paths_m[0]) - 1
    B = np.zeros((len(paths_m1), len(paths_m)), dtype=int)
    for j, path in enumerate(paths_m):
        for i in range(1, m):
            face = path[:i] + path[i+1:]
            if face in idx:
                B[idx[face], j] += (-1)**i
    return B

def check_d2(A):
    n = A.shape[0]
    all_paths = {m: get_regular_paths(A, m) for m in range(n)}
    for m in range(2, n):
        if all_paths[m] and all_paths[m-1] and all_paths[m-2]:
            B_m = interior_boundary(all_paths[m], all_paths[m-1])
            B_m1 = interior_boundary(all_paths[m-1], all_paths[m-2])
            d2 = B_m1 @ B_m
            if np.max(np.abs(d2)) > 0:
                return False
    return True

def count_kcycles(A, n, k):
    """Count directed k-cycles (divided by k for cyclic rotations)."""
    count = 0
    for combo in combinations(range(n), k):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[(i+1) % k]] for i in range(k)):
                count += 1
    return count // k

def canon_form(A):
    n = A.shape[0]
    best = None
    for perm in permutations(range(n)):
        enc = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1,n))
        if best is None or enc < best: best = enc
    return best

# ============================================================
print("="*70)
print("n=6: CYCLE-BASED PREDICTOR FOR d²=0")
print("="*70)

n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
types6 = {}
for bits in range(2**len(pairs)):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    cf = canon_form(A)
    if cf not in types6:
        types6[cf] = A

print(f"  {len(types6)} isomorphism types")

# Collect full invariant profiles
rows = []
for cf, A in types6.items():
    t3 = count_kcycles(A, n, 3)
    c4 = count_kcycles(A, n, 4)
    c5 = count_kcycles(A, n, 5)
    c6 = count_kcycles(A, n, 6)
    scores = tuple(sorted(int(sum(A[v][w] for w in range(n) if w != v)) for v in range(n)))
    d2ok = check_d2(A)
    rows.append({
        'cf': cf, 'scores': scores, 't3': t3, 'c4': c4, 'c5': c5, 'c6': c6,
        'd2ok': d2ok
    })

# Test: is (t3, c4) a perfect predictor?
t3c4_groups = defaultdict(list)
for r in rows:
    t3c4_groups[(r['t3'], r['c4'])].append(r['d2ok'])

mixed_t3c4 = [(k, v) for k, v in t3c4_groups.items() if len(set(v)) > 1]
print(f"\n  (t3, c4) as predictor: {len(mixed_t3c4)} mixed groups out of {len(t3c4_groups)}")
if mixed_t3c4:
    for (t3, c4), vals in sorted(mixed_t3c4):
        print(f"    (t3={t3}, c4={c4}): {sum(vals)} pass, {len(vals)-sum(vals)} fail")

# Test: is (t3, c4, c5) a perfect predictor?
t3c4c5_groups = defaultdict(list)
for r in rows:
    t3c4c5_groups[(r['t3'], r['c4'], r['c5'])].append(r['d2ok'])

mixed_t3c4c5 = [(k, v) for k, v in t3c4c5_groups.items() if len(set(v)) > 1]
print(f"  (t3, c4, c5) as predictor: {len(mixed_t3c4c5)} mixed groups")

# Test: is (t3, c4, c5, c6) a perfect predictor?
all_cycles_groups = defaultdict(list)
for r in rows:
    all_cycles_groups[(r['t3'], r['c4'], r['c5'], r['c6'])].append(r['d2ok'])

mixed_all = [(k, v) for k, v in all_cycles_groups.items() if len(set(v)) > 1]
print(f"  (t3, c4, c5, c6) as predictor: {len(mixed_all)} mixed groups")

# What DOES distinguish them? Try H (Hamiltonian path count)
def count_hp(A):
    n = A.shape[0]
    count = 0
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            count += 1
    return count

if mixed_all:
    print(f"\n  Adding H (Hamiltonian path count):")
    for r in rows:
        r['H'] = count_hp(types6[r['cf']])

    all_plus_H = defaultdict(list)
    for r in rows:
        all_plus_H[(r['t3'], r['c4'], r['c5'], r['c6'], r['H'])].append(r['d2ok'])
    mixed_H = [(k, v) for k, v in all_plus_H.items() if len(set(v)) > 1]
    print(f"  (t3,c4,c5,c6,H) as predictor: {len(mixed_H)} mixed groups")

# Actually, let's try: does the score sequence help?
score_groups = defaultdict(list)
for r in rows:
    score_groups[(r['scores'], r['t3'])].append(r['d2ok'])
mixed_score = [(k, v) for k, v in score_groups.items() if len(set(v)) > 1]
print(f"\n  (score, t3) as predictor: {len(mixed_score)} mixed groups")

# Full profile
full_groups = defaultdict(list)
for r in rows:
    full_groups[(r['scores'], r['t3'], r['c4'])].append(r['d2ok'])
mixed_full = [(k, v) for k, v in full_groups.items() if len(set(v)) > 1]
print(f"  (score, t3, c4) as predictor: {len(mixed_full)} mixed groups")

# ============================================================
print(f"\n{'='*70}")
print("FULL TABLE: ALL n=6 TYPES WITH d²=0 STATUS")
print("="*70)

print(f"\n  {'score':25s} {'t3':>3s} {'c4':>3s} {'c5':>3s} {'c6':>3s} {'H':>4s} {'d²=0':>5s}")
print(f"  {'-'*25} {'-'*3} {'-'*3} {'-'*3} {'-'*3} {'-'*4} {'-'*5}")

for r in rows:
    if 'H' not in r:
        r['H'] = count_hp(types6[r['cf']])

for r in sorted(rows, key=lambda x: (x['t3'], x['c4'], x['scores'])):
    d2_str = 'PASS' if r['d2ok'] else 'FAIL'
    print(f"  {str(r['scores']):25s} {r['t3']:3d} {r['c4']:3d} {r['c5']:3d} {r['c6']:3d} {r['H']:4d} {d2_str:>5s}")

# ============================================================
# Check: is the number of REGULAR paths itself a predictor?
print(f"\n{'='*70}")
print("REGULAR PATH COUNT AS PREDICTOR")
print("="*70)

for r in rows:
    A = types6[r['cf']]
    rp = [len(get_regular_paths(A, m)) for m in range(n)]
    r['reg_paths'] = tuple(rp)

rp_groups = defaultdict(list)
for r in rows:
    rp_groups[r['reg_paths']].append(r['d2ok'])

mixed_rp = [(k, v) for k, v in rp_groups.items() if len(set(v)) > 1]
print(f"  Regular path profile as predictor: {len(mixed_rp)} mixed groups out of {len(rp_groups)}")

# Check if unique
unique_rp = all(len(v) == 1 for v in rp_groups.values())
print(f"  Each regular path profile is unique? {unique_rp}")

# If not, check if (reg_paths, d2) is consistent
if mixed_rp:
    for rp, vals in sorted(mixed_rp)[:5]:
        print(f"    RP={list(rp)}: {sum(vals)} pass, {len(vals)-sum(vals)} fail")

print("\nDONE.")
