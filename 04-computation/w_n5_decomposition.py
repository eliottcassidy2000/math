#!/usr/bin/env python3
"""
Decompose w_{n-5} into sigma contributions.

w_{n-5} = sigma_4 = sum_{|S|=4} sigma(S)

Position subsets S of size 4 from {0,...,n-2} have different "gap patterns":
  - (1,1,1): {i,i+1,i+2,i+3} — 4 consecutive positions
  - (1,1,g): {i,i+1,i+2,i+2+g} with g>=2
  - (1,g,1): {i,i+1,i+1+g,i+2+g} with g>=2
  - (g,1,1): {i,i+g,i+g+1,i+g+2} with g>=2
  - (1,g): {i,i+1,j,j+1} with j-i-1 >= 2 (two disjoint pairs)
  - etc.

For non-adjacent components: sigma = 0 (product of independent zero-mean factors).
So only CONNECTED gap patterns contribute.

Connected patterns of size 4 with gaps: (1,1,1), (2,1), (1,2), (3), (1,1,g) etc.
Actually "connected" here means the INDEX SET is consecutive or?

No - "connected" in position-adjacency means: positions form a connected
subgraph where i~j iff |i-j|=1. So a connected S of size 4 is {i,i+1,i+2,i+3}.

ANY other S of size 4 has at least two connected components, and if those
components are separated by gap >= 2, the sigma factorizes into a product of
independent sums, each with zero mean. Wait - not necessarily zero mean!

Actually for |S|=4, decomposed into connected components:
- If S = C1 union C2 with min gap >= 2:
  sigma(S) = sigma(C1) * sigma(C2) / (some normalization)?

No. sigma(S) = sum_P prod_{i in S} s_{p_i,p_{i+1}}.
If S = {a,a+1,c,c+1} with c-a-1 >= 2:
  sigma(S) = sum_P s_{p_a,p_{a+1}} * s_{p_{a+1},p_{a+2}} * s_{p_c,p_{c+1}} * s_{p_{c+1},p_{c+2}}

Wait, I confused myself. |S| = 4 means we pick 4 positions, and the PRODUCT is over
those 4 positions. So S = {i_1, i_2, i_3, i_4} and sigma(S) = sum_P prod_{i in S} s_{p_i,p_{i+1}}.

If S has two disconnected components (in the adjacency sense), say C1 = {a,a+1} and C2 = {c,c+1}
with c > a+2, then... the product factorizes into two independent parts? Not quite,
because the permutation is shared.

Actually, this is exactly opus's sigma pattern theory: sigma(S) depends on the
connected-component gap structure. And for S with disconnected components,
sigma(S) does NOT necessarily equal the product of individual sigmas.

Let me just compute at n=7 to see which patterns contribute.
kind-pasteur-2026-03-06-S25g
"""
from itertools import permutations, combinations
from math import factorial, comb
import random

def compute_sigma(A, n, positions):
    total = 0.0
    for p in permutations(range(n)):
        prod = 1.0
        for i in positions:
            prod *= (A[p[i]][p[i+1]] - 0.5)
        total += prod
    return total

n = 7
random.seed(42)

# Get all S of size 4 from {0,...,5} and classify by gap pattern
def gap_pattern(S):
    s = sorted(S)
    return tuple(s[i+1] - s[i] for i in range(len(s)-1))

all_S4 = list(combinations(range(n-1), 4))
patterns = {}
for S in all_S4:
    gp = gap_pattern(S)
    if gp not in patterns:
        patterns[gp] = []
    patterns[gp].append(S)

print(f"n={n}: |S|=4 position subsets from {{0,...,{n-2}}}")
print(f"Total: {len(all_S4)} subsets")
print(f"Patterns: {sorted(patterns.keys())}")
for gp, sets in sorted(patterns.items()):
    print(f"  {gp}: {len(sets)} subsets, e.g. {sets[0]}")

# Now compute sigma(S) for a few tournaments
print(f"\n{'='*60}")

for trial in range(5):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = sum(1 for v in combinations(range(n), 3)
             for p in permutations(v)
             if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3))) // 3

    print(f"\ntrial={trial}, t3={t3}")

    for gp, sets in sorted(patterns.items()):
        vals = [compute_sigma(A, n, S) for S in sets]
        avg = sum(vals) / len(vals)
        all_same = all(abs(v - vals[0]) < 0.01 for v in vals)
        total = sum(vals)
        print(f"  pattern {gp}: per_set={vals[0]:.4f}, all_same={'Y' if all_same else 'N'}, total={total:.4f}")

# Compute total w_{n-5} = w_2
print(f"\n{'='*60}")
print(f"w_{{n-5}} = w_2 decomposition at n={n}")

random.seed(42)
for trial in range(5):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = sum(1 for v in combinations(range(n), 3)
             for p in permutations(v)
             if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3))) // 3

    pattern_totals = {}
    for gp, sets in sorted(patterns.items()):
        total = sum(compute_sigma(A, n, S) for S in sets)
        pattern_totals[gp] = total

    w2_total = sum(pattern_totals.values())

    # Known formula: w2 = -60*t3 + 12*t5 + 24*a2 + 231
    print(f"\n  trial={trial}: t3={t3}, w2={w2_total:.4f}")
    for gp, total in sorted(pattern_totals.items()):
        print(f"    pattern {gp}: {total:.4f}")

print("\nDONE")
