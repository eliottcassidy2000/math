#!/usr/bin/env python3
"""
Large-sample verification of H = 1 + 2(t₃+t₅+t₇) + 4·d₃₃ at n=7
opus-2026-03-14-S89b
"""
from itertools import combinations, permutations
import random, sys

n = 7
m = n*(n-1)//2  # 21

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

# Precompute 5-subset permutations
perms5 = {}
for combo in combinations(range(n), 5):
    perms5[combo] = list(permutations(combo))

# Precompute 7-permutations (fix vertex 0)
perms6 = list(permutations(range(1, n)))

random.seed(12345)
sample_size = 50000
errors = 0
max_error = 0

for trial in range(sample_size):
    if trial % 5000 == 0:
        print(f"  {trial}/{sample_size}...", file=sys.stderr, flush=True)

    bits = random.randint(0, (1 << m) - 1)

    adj = [[False]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            idx += 1

    H = compute_H(adj)

    # t₃
    t3 = 0
    triples = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    t3 += 1
                    triples.append(frozenset([i,j,k]))

    # t₅
    t5 = 0
    for combo in combinations(range(n), 5):
        for perm in perms5[combo]:
            ok = True
            for idx2 in range(5):
                if not adj[perm[idx2]][perm[(idx2+1)%5]]:
                    ok = False
                    break
            if ok:
                t5 += 1
    t5 //= 5

    # t₇: directed Hamiltonian cycles (fix start at 0)
    t7 = 0
    for perm in perms6:
        if not adj[0][perm[0]]:
            continue
        ok = True
        for idx2 in range(len(perm)-1):
            if not adj[perm[idx2]][perm[idx2+1]]:
                ok = False
                break
        if ok and adj[perm[-1]][0]:
            t7 += 1
    # fixing start at 0 counts each directed cycle once

    # d₃₃
    d33 = 0
    for i in range(len(triples)):
        for j in range(i+1, len(triples)):
            if len(triples[i] & triples[j]) == 0:
                d33 += 1

    predicted = 1 + 2*(t3 + t5 + t7) + 4*d33

    if predicted != H:
        errors += 1
        err = abs(H - predicted)
        if err > max_error:
            max_error = err
        if errors <= 5:
            print(f"  MISMATCH #{errors}: H={H}, pred={predicted}, t3={t3}, t5={t5}, t7={t7}, d33={d33}, diff={H-predicted}")

print(f"\nResults: {errors}/{sample_size} mismatches")
if errors == 0:
    print(f"*** ALL {sample_size} SAMPLES MATCH! Formula H=1+2(t₃+t₅+t₇)+4d₃₃ very likely exact at n=7 ***")
else:
    print(f"Max error = {max_error}")
