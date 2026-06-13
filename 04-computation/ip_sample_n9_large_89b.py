#!/usr/bin/env python3
"""
Larger n=9 sample (5000) and preliminary n=10 test (100)
opus-2026-03-14-S89b

n=9: H = 1 + 2*(t3+t5+t7+t9) + 4*(d33+d35) + 8*d333
n=10: NEW terms:
  - d55: disjoint (5,5) pair needs 10 = n, POSSIBLE
  - d37: disjoint (3,7) pair needs 10 = n, POSSIBLE
  - d335: disjoint (3,3,5) triple needs 11 > 10, IMPOSSIBLE
  - d3333: need 12 > 10, IMPOSSIBLE
  So at n=10: H = 1 + 2*Σt_{odd} + 4*(d33+d35+d55+d37) + 8*d333
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
    """Disjoint (3-cycle, 7-cycle) pairs."""
    count = 0
    for tri in triples:
        remaining = sorted(set(range(n)) - tri)
        if len(remaining) < 7: continue
        for combo7 in combinations(remaining, 7):
            count += count_directed_kcycles_on_subset(adj, list(combo7))
    return count

def count_d55(n, adj):
    """Disjoint (5-cycle, 5-cycle) pairs."""
    count = 0
    for combo5a in combinations(range(n), 5):
        c5a = count_directed_kcycles_on_subset(adj, list(combo5a))
        if c5a == 0: continue
        remaining = sorted(set(range(n)) - set(combo5a))
        if len(remaining) < 5: continue
        for combo5b in combinations(remaining, 5):
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

# ===== N=9 LARGE SAMPLE =====
print("="*60)
print("N=9 LARGE SAMPLE: 5000 tournaments")
print("="*60)

n = 9
m = n*(n-1)//2
random.seed(77777)
sample_size = 5000
errors = 0
t0 = time.time()

for trial in range(sample_size):
    if trial % 500 == 0:
        print(f"  {trial}/{sample_size}...", file=sys.stderr, flush=True)

    bits = random.randint(0, (1 << m) - 1)
    adj = tournament_adj(n, bits)
    H = compute_H(n, adj)

    triples = get_3cycle_triples(n, adj)
    t3 = len(triples)
    t5 = count_kcycles(n, adj, 5)
    t7 = count_kcycles(n, adj, 7)
    t9 = count_directed_kcycles_on_subset(adj, list(range(n)))
    d33 = count_d33(triples)
    d35 = count_d35(n, adj, triples)
    d333 = count_d333(triples)

    pred = 1 + 2*(t3+t5+t7+t9) + 4*(d33+d35) + 8*d333

    if pred != H:
        errors += 1
        if errors <= 3:
            print(f"  MISMATCH: H={H}, pred={pred}, diff={H-pred}")

print(f"\nn=9: {errors}/{sample_size} mismatches ({time.time()-t0:.1f}s)")
if errors == 0:
    print("*** ALL MATCH at n=9! ***")

# ===== N=10 SMALL SAMPLE =====
print("\n" + "="*60)
print("N=10 SMALL SAMPLE: 100 tournaments")
print("Testing H = 1 + 2*Σt + 4*(d33+d35+d37+d55) + 8*d333")
print("="*60)

n = 10
m = n*(n-1)//2  # 45
random.seed(11111)
sample_size = 100
errors = 0
t0 = time.time()

for trial in range(sample_size):
    if trial % 10 == 0:
        print(f"  {trial}/{sample_size}...", file=sys.stderr, flush=True)

    bits = random.randint(0, (1 << m) - 1)
    adj = tournament_adj(n, bits)
    H = compute_H(n, adj)

    triples = get_3cycle_triples(n, adj)
    t3 = len(triples)
    t5 = count_kcycles(n, adj, 5)
    t7 = count_kcycles(n, adj, 7)
    t9 = count_kcycles(n, adj, 9)  # 9-cycles on 9-subsets (C(10,9)=10 subsets)

    d33 = count_d33(triples)
    d35 = count_d35(n, adj, triples)
    d37 = count_d37(n, adj, triples)
    d55 = count_d55(n, adj)
    d333 = count_d333(triples)

    pred = 1 + 2*(t3+t5+t7+t9) + 4*(d33+d35+d37+d55) + 8*d333

    if pred != H:
        errors += 1
        if errors <= 5:
            print(f"  MISMATCH #{errors}: H={H}, pred={pred}, diff={H-pred}, "
                  f"t3={t3} t5={t5} t7={t7} t9={t9} d33={d33} d35={d35} d37={d37} d55={d55} d333={d333}")

elapsed = time.time() - t0
print(f"\nn=10: {errors}/{sample_size} mismatches ({elapsed:.1f}s)")
if errors == 0:
    print("*** ALL MATCH at n=10! ***")
else:
    print(f"Collecting residuals...")
