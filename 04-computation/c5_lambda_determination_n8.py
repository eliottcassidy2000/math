"""
c5_lambda_determination_n8.py -- kind-pasteur-2026-03-13-S61
Does lambda determine c5_dir at n=8?
If yes, dc5=0 holds at n=8 (and likely all n).
If no, dc5=0 would need a different explanation.

n=8: C(8,2)=28 bits, 2^28 = 268M tournaments. Too many for exhaustive.
Sample instead.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def count_total_directed_k_cycles(A, n, k):
    total = 0
    for combo in combinations(range(n), k):
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = True
            for i in range(k):
                if A[path[i]][path[(i+1) % k]] != 1:
                    valid = False
                    break
            if valid:
                total += 1
    return total

n = 8
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"IS c5_dir A FUNCTION OF lambda? (n={n}, sampled)")
print("=" * 60)

np.random.seed(42)
lambda_groups = {}

for trial in range(20000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    key = tuple(L[i][j] for i in range(n) for j in range(n))
    c5 = count_total_directed_k_cycles(A, n, 5)
    if key not in lambda_groups:
        lambda_groups[key] = set()
    lambda_groups[key].add(c5)

    if trial % 5000 == 0 and trial > 0:
        amb = sum(1 for v in lambda_groups.values() if len(v) > 1)
        print(f"  ... {trial}/20000, groups={len(lambda_groups)}, ambiguous={amb}")

ambiguous = sum(1 for v in lambda_groups.values() if len(v) > 1)
print(f"\n  Total lambda groups: {len(lambda_groups)}")
print(f"  Ambiguous: {ambiguous}")
print(f"  c5 IS a function of lambda? {ambiguous == 0}")

if ambiguous > 0:
    shown = 0
    for key, vals in lambda_groups.items():
        if len(vals) > 1 and shown < 5:
            print(f"  Ambiguous: c5 values = {sorted(vals)}")
            shown += 1

# Also check c7 at n=8
print(f"\n{'='*60}")
print(f"IS c7_dir A FUNCTION OF lambda? (n={n}, sampled)")
print(f"{'='*60}")

np.random.seed(42)
lambda_groups_c7 = {}

for trial in range(10000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    key = tuple(L[i][j] for i in range(n) for j in range(n))
    c7 = count_total_directed_k_cycles(A, n, 7)
    if key not in lambda_groups_c7:
        lambda_groups_c7[key] = set()
    lambda_groups_c7[key].add(c7)

ambiguous_c7 = sum(1 for v in lambda_groups_c7.values() if len(v) > 1)
print(f"  Total lambda groups: {len(lambda_groups_c7)}")
print(f"  Ambiguous: {ambiguous_c7}")
print(f"  c7 IS a function of lambda? {ambiguous_c7 == 0}")

if ambiguous_c7 > 0:
    shown = 0
    for key, vals in lambda_groups_c7.items():
        if len(vals) > 1 and shown < 5:
            print(f"  Ambiguous: c7 values = {sorted(vals)}")
            shown += 1

print("\nDone.")
