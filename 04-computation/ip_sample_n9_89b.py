#!/usr/bin/env python3
"""
N=9 SAMPLING: Test the full Independence Polynomial formula at n=9
opus-2026-03-14-S89b

At n=9, new possible terms:
- t_9: directed 9-cycles (Hamiltonian cycles on 9-subsets... but at n=9, only one 9-subset)
- d333: three mutually disjoint 3-cycles (need 9 = n vertices exactly)
- d37: disjoint (3,7) pair needs 10 > 9, IMPOSSIBLE
- d55: disjoint (5,5) pair needs 10 > 9, IMPOSSIBLE
- d35: disjoint (3,5) pair needs 8 ≤ 9, POSSIBLE
- higher: d335 = disjoint (3,3,5) needs 11 > 9, IMPOSSIBLE

So at n=9:
  H = 1 + 2*(t3+t5+t7+t9) + 4*(d33+d35) + 8*d333

Level 0: 1
Level 1: 2*(t3+t5+t7+t9)
Level 2: 4*(d33+d35)
Level 3: 8*d333

This is a HUGE test — the first time level 3 (coefficient 8) appears!

Note: t9 = number of directed Hamiltonian 9-cycles at n=9 (fix start vertex,
check all 8! = 40320 permutations). This is expensive.

Strategy: small sample (500) but thorough computation.
"""

from itertools import combinations, permutations
import random, sys, time

n = 9
m = n*(n-1)//2  # 36

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

def count_directed_kcycles_on_subset(adj, verts):
    """Count directed k-cycles on a vertex subset (fix first vertex)."""
    k = len(verts)
    if k < 3:
        return 0
    count = 0
    v0 = verts[0]
    for perm in permutations(verts[1:]):
        ok = adj[v0][perm[0]]
        if not ok:
            continue
        for idx in range(len(perm)-1):
            if not adj[perm[idx]][perm[idx+1]]:
                ok = False
                break
        if ok and adj[perm[-1]][v0]:
            count += 1
    return count

def count_kcycles(adj, k):
    """Count all directed k-cycles in the tournament."""
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

def count_d35(adj, triples):
    """Count disjoint (3-cycle, 5-cycle) pairs.
    For each 3-cycle, count 5-cycles on disjoint vertices."""
    count = 0
    for tri in triples:
        remaining = sorted(set(range(n)) - tri)
        for combo5 in combinations(remaining, 5):
            c5 = count_directed_kcycles_on_subset(adj, list(combo5))
            count += c5
    return count

def count_d333(triples):
    """Count mutually disjoint 3-cycle triples."""
    count = 0
    nt = len(triples)
    for i in range(nt):
        for j in range(i+1, nt):
            if len(triples[i] & triples[j]) > 0:
                continue
            for k in range(j+1, nt):
                if len(triples[i] & triples[k]) == 0 and \
                   len(triples[j] & triples[k]) == 0:
                    count += 1
    return count

random.seed(99999)
sample_size = 500
errors_base = 0
errors_full = 0
max_error_base = 0
max_error_full = 0

print(f"N=9 SAMPLING")
print(f"Testing H = 1 + 2*(t3+t5+t7+t9) + 4*(d33+d35) + 8*d333")
print(f"Sample size: {sample_size}")
print()

t0 = time.time()

for trial in range(sample_size):
    if trial % 50 == 0:
        elapsed = time.time() - t0
        rate = trial / elapsed if elapsed > 0 else 0
        eta = (sample_size - trial) / rate if rate > 0 else 0
        print(f"  {trial}/{sample_size}... ({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)",
              file=sys.stderr, flush=True)

    bits = random.randint(0, (1 << m) - 1)
    adj = tournament_adj(bits)
    H = compute_H(adj)

    triples = get_3cycle_triples(adj)
    t3 = len(triples)
    t5 = count_kcycles(adj, 5)
    t7 = count_kcycles(adj, 7)

    # t9: Hamiltonian 9-cycles (fix vertex 0)
    t9 = count_directed_kcycles_on_subset(adj, list(range(n)))

    d33 = count_d33(triples)
    d35 = count_d35(adj, triples)
    d333 = count_d333(triples)

    # Base formula (without d333 or t9)
    pred_base = 1 + 2*(t3 + t5 + t7) + 4*(d33 + d35)
    # Full formula (with d333 and t9)
    pred_full = 1 + 2*(t3 + t5 + t7 + t9) + 4*(d33 + d35) + 8*d333

    if pred_base != H:
        errors_base += 1
        err = abs(H - pred_base)
        if err > max_error_base:
            max_error_base = err

    if pred_full != H:
        errors_full += 1
        err = abs(H - pred_full)
        if err > max_error_full:
            max_error_full = err
        if errors_full <= 5:
            print(f"  FULL MISMATCH #{errors_full}: H={H}, pred={pred_full}, "
                  f"t3={t3}, t5={t5}, t7={t7}, t9={t9}, "
                  f"d33={d33}, d35={d35}, d333={d333}, diff={H-pred_full}")

elapsed = time.time() - t0
print(f"\nDone in {elapsed:.1f}s")

print(f"\nBase formula (no t9, no d333): {errors_base}/{sample_size} mismatches, max error = {max_error_base}")
print(f"Full formula (with t9 and d333): {errors_full}/{sample_size} mismatches, max error = {max_error_full}")

if errors_full == 0:
    print(f"\n*** ALL {sample_size} MATCH! Independence polynomial formula EXACT at n=9! ***")
    print(f"*** First verification of 2^3 = 8 coefficient for disjoint triples! ***")
elif errors_base == 0:
    print(f"\n*** Base formula without d333/t9 already exact! ***")
else:
    print(f"\nNeither formula exact. Collecting residual data...")
