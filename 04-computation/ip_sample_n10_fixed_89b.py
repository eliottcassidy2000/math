#!/usr/bin/env python3
"""
N=10 FIXED: d55 was double-counted! Each unordered (5,5) pair counted twice.
opus-2026-03-14-S89b

Fix: divide d55 by 2 (each pair {A,B} enumerated as both (A,B) and (B,A)).

Also test: does this fix make n=10 work perfectly?
"""

from itertools import combinations, permutations
import random, sys, time

def compute_H(n, adj):
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

def tournament_adj(n, bits):
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

def get_3cycle_triples(n, adj):
    triples = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    triples.append(frozenset([i,j,k]))
    return triples

def count_directed_kcycles_on_subset(adj, verts):
    k = len(verts)
    if k < 3: return 0
    count = 0
    v0 = verts[0]
    for perm in permutations(verts[1:]):
        ok = adj[v0][perm[0]]
        if not ok: continue
        for idx in range(len(perm)-1):
            if not adj[perm[idx]][perm[idx+1]]:
                ok = False
                break
        if ok and adj[perm[-1]][v0]:
            count += 1
    return count

def count_kcycles(n, adj, k):
    total = 0
    for combo in combinations(range(n), k):
        total += count_directed_kcycles_on_subset(adj, list(combo))
    return total

def count_d33(triples):
    count = 0
    for i in range(len(triples)):
        for j in range(i+1, len(triples)):
            if len(triples[i] & triples[j]) == 0:
                count += 1
    return count

def count_d35(n, adj, triples):
    count = 0
    for tri in triples:
        remaining = sorted(set(range(n)) - tri)
        for combo5 in combinations(remaining, 5):
            count += count_directed_kcycles_on_subset(adj, list(combo5))
    return count

def count_d37(n, adj, triples):
    count = 0
    for tri in triples:
        remaining = sorted(set(range(n)) - tri)
        if len(remaining) < 7: continue
        for combo7 in combinations(remaining, 7):
            count += count_directed_kcycles_on_subset(adj, list(combo7))
    return count

def count_d55(n, adj):
    """Disjoint (5,5) pairs — FIXED: use ordered pairs with first < second."""
    count = 0
    all_5subsets = list(combinations(range(n), 5))
    for idx_a in range(len(all_5subsets)):
        combo5a = all_5subsets[idx_a]
        c5a = count_directed_kcycles_on_subset(adj, list(combo5a))
        if c5a == 0: continue
        remaining = set(range(n)) - set(combo5a)
        for idx_b in range(idx_a+1, len(all_5subsets)):
            combo5b = all_5subsets[idx_b]
            if set(combo5b) != remaining and not set(combo5b).issubset(remaining):
                continue
            if set(combo5b).issubset(remaining):
                c5b = count_directed_kcycles_on_subset(adj, list(combo5b))
                count += c5a * c5b
    return count

def count_d333(triples):
    count = 0
    nt = len(triples)
    for i in range(nt):
        for j in range(i+1, nt):
            if len(triples[i] & triples[j]) > 0: continue
            for k in range(j+1, nt):
                if len(triples[i] & triples[k]) == 0 and \
                   len(triples[j] & triples[k]) == 0:
                    count += 1
    return count

n = 10
m = n*(n-1)//2  # 45
random.seed(11111)  # Same seed as before
sample_size = 200

print(f"N=10 FIXED SAMPLING: {sample_size} tournaments")
print(f"Formula: H = 1 + 2*Σt_odd + 4*(d33+d35+d37+d55) + 8*d333")
print(f"d55 now properly counted (unordered pairs)")
print()

errors = 0
max_error = 0
t0 = time.time()

for trial in range(sample_size):
    if trial % 20 == 0:
        print(f"  {trial}/{sample_size}...", file=sys.stderr, flush=True)

    bits = random.randint(0, (1 << m) - 1)
    adj = tournament_adj(n, bits)
    H = compute_H(n, adj)

    triples = get_3cycle_triples(n, adj)
    t3 = len(triples)
    t5 = count_kcycles(n, adj, 5)
    t7 = count_kcycles(n, adj, 7)
    t9 = count_kcycles(n, adj, 9)

    d33 = count_d33(triples)
    d35 = count_d35(n, adj, triples)
    d37 = count_d37(n, adj, triples)
    d55 = count_d55(n, adj)
    d333 = count_d333(triples)

    pred = 1 + 2*(t3+t5+t7+t9) + 4*(d33+d35+d37+d55) + 8*d333

    if pred != H:
        errors += 1
        err = abs(H - pred)
        if err > max_error:
            max_error = err
        if errors <= 5:
            print(f"  MISMATCH #{errors}: H={H}, pred={pred}, diff={H-pred}")
            print(f"    t3={t3} t5={t5} t7={t7} t9={t9}")
            print(f"    d33={d33} d35={d35} d37={d37} d55={d55} d333={d333}")

elapsed = time.time() - t0
print(f"\nn=10: {errors}/{sample_size} mismatches ({elapsed:.1f}s)")
if errors == 0:
    print("*** ALL MATCH at n=10! Independence Polynomial EXACT! ***")
else:
    print(f"Max error = {max_error}")
