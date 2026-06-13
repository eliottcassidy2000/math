#!/usr/bin/env python3
"""
Mechanism of the H=21 gap: why sw(1) + sw(2) = 20 is never achievable.

INSIGHT: The gap arises from a parity-like constraint on (t_3, t_5):

At n=6:
- sw(1) = 2*t_3
- sw(2) = 2*t_5 + 4*p_{33}

For H=21: need sw(1) + sw(2) = 20 = 2*t_3 + 2*t_5 + 4*p_33
=> t_3 + t_5 + 2*p_33 = 10

The key structural constraint: t_3 and t_5 are not independent.
For a tournament on n=6 vertices, the number of cyclic triples t_3
constrains t_5 (5-cycles).

Formula: t_3 = C(n,3) - sum_v C(out(v),2) where out(v) is out-degree.
For n=6: t_3 = 20 - sum C(d_v, 2).

Let's understand the (t_3, t_5) constraint exhaustively.

opus-2026-03-07-S39
"""
from itertools import combinations, permutations
from collections import defaultdict


def make_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1
    return A


def count_cycles(n, A):
    """Count directed 3-cycles and 5-cycles."""
    edge_set = {(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]}

    t3 = 0
    three_cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if A[a][b] and A[b][c] and A[c][a]:
                    t3 += 1
                    three_cycles.append(frozenset({a,b,c}))
                elif A[a][c] and A[c][b] and A[b][a]:
                    t3 += 1
                    three_cycles.append(frozenset({a,b,c}))

    t5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all((p[i], p[(i+1)%5]) in edge_set for i in range(5)):
                t5 += 1
    t5 //= 5

    dp = sum(1 for i in range(len(three_cycles))
            for j in range(i+1, len(three_cycles))
            if not (three_cycles[i] & three_cycles[j]))

    return t3, t5, dp


# === n=6: achievable (t3, t5, p33) triples ===
print("=== n=6: achievable (t3, t5, p33) and H values ===")
n = 6
triples = defaultdict(list)
total = 1 << (n*(n-1)//2)

for bits in range(total):
    A = make_tournament(n, bits)
    t3, t5, dp = count_cycles(n, A)
    triples[(t3, t5, dp)].append(bits)

print(f"Achievable (t3, t5, p33) triples:")
for (t3, t5, dp) in sorted(triples.keys()):
    sw1 = 2*t3
    sw2 = 2*t5 + 4*dp
    H = 1 + sw1 + sw2
    target = t3 + t5 + 2*dp
    marker = " *** TARGET=10" if target == 10 else ""
    print(f"  ({t3}, {t5}, {dp}): H={H}, count={len(triples[(t3, t5, dp)])}, "
          f"t3+t5+2p33={target}{marker}")

# === Find WHY target=10 is impossible ===
print(f"\n=== Why t3 + t5 + 2*p33 = 10 is impossible at n=6 ===")
print(f"Need: (t3, t5, p33) with t3 + t5 + 2*p33 = 10")
print()

# List all (t3, target) pairs
print(f"Achievable (t3, t3+t5+2*p33) pairs:")
for (t3, t5, dp) in sorted(triples.keys()):
    target = t3 + t5 + 2*dp
    if 8 <= target <= 12:
        print(f"  t3={t3}: target = {t3} + {t5} + 2*{dp} = {target}")

# === The constraint: what t5 values are forced by t3? ===
print(f"\n=== t5 as function of t3 at n=6 ===")
t3_to_t5 = defaultdict(set)
for (t3, t5, dp) in triples.keys():
    t3_to_t5[t3].add(t5)

for t3 in sorted(t3_to_t5.keys()):
    print(f"  t3={t3}: achievable t5 = {sorted(t3_to_t5[t3])}")

# === Same analysis at n=7 ===
print(f"\n=== n=7: what (t3, t5) pairs achieve H near 21? ===")
n = 7

# Sample random tournaments
import random
random.seed(42)

h_near_21 = defaultdict(list)
t3_t5_by_H = defaultdict(set)

num_samples = 200000
for _ in range(num_samples):
    A = [[0]*7 for _ in range(7)]
    for i in range(7):
        for j in range(i+1, 7):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    # Quick H computation
    from functools import lru_cache

    adj = [0]*7
    for i in range(7):
        for j in range(7):
            if i != j and A[i][j]:
                adj[i] |= (1 << j)

    dp = [[0]*7 for _ in range(1<<7)]
    for v in range(7):
        dp[1<<v][v] = 1
    for S in range(1, 1<<7):
        for v in range(7):
            if not (S & (1<<v)): continue
            if dp[S][v] == 0: continue
            for u in range(7):
                if S & (1<<u): continue
                if adj[v] & (1<<u):
                    dp[S|(1<<u)][u] += dp[S][v]
    H = sum(dp[(1<<7)-1][v] for v in range(7))

    if 15 <= H <= 27:
        t3 = 0
        for a in range(7):
            for b in range(a+1, 7):
                for c in range(b+1, 7):
                    if A[a][b] and A[b][c] and A[c][a]:
                        t3 += 1
                    elif A[a][c] and A[c][b] and A[b][a]:
                        t3 += 1
        t3_t5_by_H[H].add(t3)

for h in range(15, 28):
    if h in t3_t5_by_H:
        print(f"  H={h}: t3 in {sorted(t3_t5_by_H[h])}")
    else:
        print(f"  H={h}: NOT FOUND (in {num_samples} samples)")
