#!/usr/bin/env python3
"""
N=8 SAMPLING: Test H = 1 + 2(t3+t5+t7) + 4*(d33+d35) at n=8
opus-2026-03-14-S89b

At n=8: new possible terms:
- d35: disjoint (3-cycle, 5-cycle) pairs (need 8 vertices = n exactly)
- d55: disjoint (5,5) pairs need 10 > 8, IMPOSSIBLE
- d333: disjoint 3-cycle triples need 9 > 8, IMPOSSIBLE
- So only d33 and d35 at level 2, nothing at level 3.

Formula: H = 1 + 2(t3+t5+t7) + 4*(d33+d35)
"""

from itertools import combinations, permutations
import random, sys

n = 8
m = n*(n-1)//2  # 28

def compute_H(adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            val = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u): continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def tournament_adj(bits):
    adj = [[False]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            idx += 1
    return adj

def get_3cycle_triples(adj):
    triples = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    triples.append(frozenset([i,j,k]))
    return triples

def count_5cycles_on_subset(adj, verts):
    """Count directed 5-cycles on a specific 5-vertex subset."""
    count = 0
    v = verts
    # Fix v[0], permute rest
    for perm in permutations(v[1:]):
        ok = adj[v[0]][perm[0]]
        if not ok: continue
        for idx in range(3):
            if not adj[perm[idx]][perm[idx+1]]:
                ok = False
                break
        if ok and adj[perm[3]][v[0]]:
            count += 1
    return count

def count_5cycles(adj):
    total = 0
    for combo in combinations(range(n), 5):
        total += count_5cycles_on_subset(adj, list(combo))
    return total

def count_7cycles_on_subset(adj, verts):
    """Count directed 7-cycles on a specific 7-vertex subset."""
    count = 0
    v = verts
    for perm in permutations(v[1:]):
        ok = adj[v[0]][perm[0]]
        if not ok: continue
        for idx in range(5):
            if not adj[perm[idx]][perm[idx+1]]:
                ok = False
                break
        if ok and adj[perm[5]][v[0]]:
            count += 1
    return count

def count_7cycles(adj):
    total = 0
    for combo in combinations(range(n), 7):
        total += count_7cycles_on_subset(adj, list(combo))
    return total

def count_d33(triples):
    count = 0
    for i in range(len(triples)):
        for j in range(i+1, len(triples)):
            if len(triples[i] & triples[j]) == 0:
                count += 1
    return count

def count_d35(adj, triples):
    """Count disjoint (3-cycle, 5-cycle) pairs."""
    count = 0
    for tri in triples:
        remaining = sorted(set(range(n)) - tri)
        if len(remaining) < 5:
            continue
        for combo5 in combinations(remaining, 5):
            c5 = count_5cycles_on_subset(adj, list(combo5))
            count += c5
    return count

random.seed(54321)
sample_size = 5000
errors = 0
max_error = 0

print(f"N=8 SAMPLING: Testing H = 1 + 2(t3+t5+t7) + 4*(d33+d35)")
print(f"Sample size: {sample_size}")

for trial in range(sample_size):
    if trial % 500 == 0:
        print(f"  {trial}/{sample_size}...", file=sys.stderr, flush=True)

    bits = random.randint(0, (1 << m) - 1)
    adj = tournament_adj(bits)
    H = compute_H(adj)

    triples = get_3cycle_triples(adj)
    t3 = len(triples)
    t5 = count_5cycles(adj)
    t7 = count_7cycles(adj)
    d33 = count_d33(triples)
    d35 = count_d35(adj, triples)

    predicted = 1 + 2*(t3 + t5 + t7) + 4*(d33 + d35)

    if predicted != H:
        errors += 1
        err = abs(H - predicted)
        if err > max_error:
            max_error = err
        if errors <= 10:
            print(f"  MISMATCH #{errors}: H={H}, pred={predicted}, "
                  f"t3={t3}, t5={t5}, t7={t7}, d33={d33}, d35={d35}, "
                  f"diff={H-predicted}")

print(f"\nResults: {errors}/{sample_size} mismatches")
if errors == 0:
    print(f"*** ALL {sample_size} MATCH! Formula likely exact at n=8 ***")
else:
    print(f"Max error = {max_error}")
    # Check if residual is always even
    print("\nCollecting residual statistics...")
    random.seed(54321)
    resids = []
    for trial in range(min(sample_size, 1000)):
        bits = random.randint(0, (1 << m) - 1)
        adj = tournament_adj(bits)
        H = compute_H(adj)
        triples = get_3cycle_triples(adj)
        t3 = len(triples)
        t5 = count_5cycles(adj)
        t7 = count_7cycles(adj)
        d33 = count_d33(triples)
        d35 = count_d35(adj, triples)
        predicted = 1 + 2*(t3 + t5 + t7) + 4*(d33 + d35)
        if predicted != H:
            resids.append(H - predicted)

    if resids:
        from math import gcd
        from functools import reduce
        g = reduce(gcd, [abs(r) for r in resids])
        print(f"  Residuals GCD: {g}")
        print(f"  Residual values: {sorted(set(resids))[:20]}")
        print(f"  All even? {all(r % 2 == 0 for r in resids)}")
        print(f"  All divisible by 4? {all(r % 4 == 0 for r in resids)}")
        print(f"  All divisible by 8? {all(r % 8 == 0 for r in resids)}")
