#!/usr/bin/env python3
"""
trh_d2_t3_deeper.py — opus-2026-03-13-S71

Deep dive into the TRH d²=0 condition.

At n=5: d²≠0 ⟺ t3=2 (exactly).
Question: Is t3 the key at n=6 too? Or is a more refined condition needed?

Also: d²=0 ⟺ d²=0 for T^op. Is this because t3(T) = t3(T^op)?
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

def count_3cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]+A[j][k]+A[k][i]==3 or A[j][i]+A[i][k]+A[k][j]==3:
                    count += 1
    return count

def canon_form(A):
    n = A.shape[0]
    best = None
    for perm in permutations(range(n)):
        enc = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1,n))
        if best is None or enc < best: best = enc
    return best

# ============================================================
print("="*70)
print("n=5: VERIFY t3=2 ⟺ d²≠0")
print("="*70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
types5 = {}
for bits in range(2**len(pairs)):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    cf = canon_form(A)
    if cf not in types5:
        types5[cf] = A

for cf, A in sorted(types5.items()):
    t3 = count_3cycles(A, n)
    d2ok = check_d2(A)
    scores = tuple(sorted(sum(A[v][w] for w in range(n) if w != v) for v in range(n)))
    if not d2ok or t3 == 2:
        print(f"  score={scores}, t3={t3}, d²=0? {d2ok}")

# ============================================================
print(f"\n{'='*70}")
print("n=6: ISOMORPHISM-CLASS ANALYSIS OF d²=0")
print("="*70)

n = 6
pairs6 = [(i,j) for i in range(n) for j in range(i+1,n)]
types6 = {}
for bits in range(2**len(pairs6)):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs6):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    cf = canon_form(A)
    if cf not in types6:
        types6[cf] = A

print(f"  {len(types6)} isomorphism types at n=6")

type_data = {}
for cf, A in types6.items():
    t3 = count_3cycles(A, n)
    d2ok = check_d2(A)
    scores = tuple(sorted(int(sum(A[v][w] for w in range(n) if w != v)) for v in range(n)))
    type_data[cf] = {'t3': t3, 'd2ok': d2ok, 'scores': scores}

# Count by t3 and d²=0
t3_d2 = defaultdict(lambda: {'pass': 0, 'fail': 0})
for td in type_data.values():
    if td['d2ok']:
        t3_d2[td['t3']]['pass'] += 1
    else:
        t3_d2[td['t3']]['fail'] += 1

print(f"\n  t3 vs d²=0 at n=6 (by isomorphism type):")
print(f"  {'t3':>3s}  {'pass':>5s}  {'fail':>5s}  {'total':>5s}")
for t3 in sorted(t3_d2.keys()):
    p, f = t3_d2[t3]['pass'], t3_d2[t3]['fail']
    print(f"  {t3:3d}  {p:5d}  {f:5d}  {p+f:5d}")

# Is t3 the sole predictor?
t3_only_fail = [t3 for t3, d in t3_d2.items() if d['pass'] == 0 and d['fail'] > 0]
t3_only_pass = [t3 for t3, d in t3_d2.items() if d['pass'] > 0 and d['fail'] == 0]
t3_mixed = [t3 for t3, d in t3_d2.items() if d['pass'] > 0 and d['fail'] > 0]

print(f"\n  t3 values where ALL types pass: {sorted(t3_only_pass)}")
print(f"  t3 values where ALL types fail: {sorted(t3_only_fail)}")
print(f"  t3 values with MIXED results:   {sorted(t3_mixed)}")

# For mixed t3 values, show what distinguishes pass from fail
if t3_mixed:
    print(f"\n  Detail for mixed t3 values:")
    for t3_val in sorted(t3_mixed):
        passing = [(cf, td) for cf, td in type_data.items()
                   if td['t3'] == t3_val and td['d2ok']]
        failing = [(cf, td) for cf, td in type_data.items()
                   if td['t3'] == t3_val and not td['d2ok']]
        print(f"\n  t3={t3_val}: {len(passing)} pass, {len(failing)} fail")
        print(f"    Passing scores: {sorted(set(td['scores'] for _, td in passing))}")
        print(f"    Failing scores: {sorted(set(td['scores'] for _, td in failing))}")

# ============================================================
print(f"\n{'='*70}")
print("n=6: COMPLEMENT SYMMETRY CHECK")
print("="*70)

comp_agree = 0
comp_disagree = 0
for cf, A in types6.items():
    # Complement: A^op[i][j] = 1 - A[i][j] for i≠j
    Aop = 1 - A - np.eye(n, dtype=int)
    cf_op = canon_form(Aop)
    if type_data[cf]['d2ok'] == type_data[cf_op]['d2ok']:
        comp_agree += 1
    else:
        comp_disagree += 1

print(f"  d²=0(T) ⟺ d²=0(T^op): {comp_agree}/{len(types6)} agree, {comp_disagree} disagree")

# Also: t3(T) = t3(T^op)?
t3_agree = 0
for cf, A in types6.items():
    Aop = 1 - A - np.eye(n, dtype=int)
    cf_op = canon_form(Aop)
    if type_data[cf]['t3'] == type_data[cf_op]['t3']:
        t3_agree += 1

print(f"  t3(T) = t3(T^op): {t3_agree}/{len(types6)}")

# ============================================================
print(f"\n{'='*70}")
print("n=6: ADDITIONAL INVARIANTS FOR MIXED t3")
print("="*70)

def count_kcycles(A, n, k):
    """Count directed k-cycles (up to cyclic rotation)."""
    count = 0
    for combo in combinations(range(n), k):
        for perm in permutations(combo):
            if all(A[perm[i]][perm[(i+1) % k]] for i in range(k)):
                count += 1
    return count // k  # divide by k rotations

for t3_val in sorted(t3_mixed):
    print(f"\n  t3={t3_val} types:")
    for cf, td in type_data.items():
        if td['t3'] != t3_val: continue
        A = types6[cf]
        # Count 4-cycles (directed)
        c4 = count_kcycles(A, n, 4)
        # Count transitive triples
        trans_triples = 0
        for i in range(n):
            for j in range(n):
                if i==j: continue
                for k in range(n):
                    if k==i or k==j: continue
                    if A[i][j] and A[j][k] and A[i][k]:
                        trans_triples += 1
        trans_triples //= 1  # each counted once (as ordered triple)

        print(f"    score={td['scores']}, t3={td['t3']}, c4={c4}, "
              f"trans={trans_triples}, d²=0? {td['d2ok']}")

print("\nDONE.")
