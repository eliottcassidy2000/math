"""
dc5_from_lambda.py -- kind-pasteur-2026-03-13-S61
Can the total directed 5-cycle count be expressed as a function of lambda(u,v)?

If so, lambda preservation would TRIVIALLY imply dc5=0.

If NOT, dc5=0 requires a more subtle argument.

The total directed 3-cycle count IS a function of lambda:
  c3_dir = (1/3) * sum_v delta(v) where delta(v) = (1/2) * sum_u lambda(v,u)
  Actually: c3_dir = (1/3) * (number of c3 vertex sets)
  And |C| = (1/6) * sum lambda(u,v)

But for directed 5-cycles, the situation is more complex.
Let's check: do two tournaments with the same lambda graph
always have the same total directed 5-cycle count?
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

# At n=5: is the total directed 5-cycle count a function of lambda?
n = 5
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"IS c5_dir A FUNCTION OF lambda? (n={n})")
print("=" * 60)

# Group tournaments by lambda graph
lambda_groups = {}
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    key = tuple(L[i][j] for i in range(n) for j in range(n))
    c5 = count_total_directed_k_cycles(A, n, 5)
    if key not in lambda_groups:
        lambda_groups[key] = set()
    lambda_groups[key].add(c5)

ambiguous = 0
for key, c5_vals in lambda_groups.items():
    if len(c5_vals) > 1:
        ambiguous += 1
        if ambiguous <= 3:
            L_matrix = np.array(key).reshape(n, n)
            print(f"  Lambda group with MULTIPLE c5 values: {sorted(c5_vals)}")
            print(f"    Lambda: {[L_matrix[i][j] for i in range(n) for j in range(i+1,n)]}")

print(f"\n  Total lambda groups: {len(lambda_groups)}")
print(f"  Ambiguous groups (multiple c5): {ambiguous}")
print(f"  c5 IS a function of lambda? {ambiguous == 0}")

# Now check n=6
n = 6
total_bits = n * (n-1) // 2

print(f"\n{'='*60}")
print(f"IS c5_dir A FUNCTION OF lambda? (n={n})")
print(f"{'='*60}")

lambda_groups = {}
for bits in range(1 << total_bits):
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    key = tuple(L[i][j] for i in range(n) for j in range(n))
    c5 = count_total_directed_k_cycles(A, n, 5)
    if key not in lambda_groups:
        lambda_groups[key] = set()
    lambda_groups[key].add(c5)

    if bits % 10000 == 0 and bits > 0:
        pass  # Don't print for speed

ambiguous = 0
examples_shown = 0
for key, c5_vals in lambda_groups.items():
    if len(c5_vals) > 1:
        ambiguous += 1
        if examples_shown < 3:
            print(f"  Ambiguous: c5 values = {sorted(c5_vals)}")
            examples_shown += 1

print(f"\n  Total lambda groups: {len(lambda_groups)}")
print(f"  Ambiguous: {ambiguous}")
print(f"  c5 IS a function of lambda? {ambiguous == 0}")

# Now check n=7 (sampled)
n = 7
total_bits = n * (n-1) // 2

print(f"\n{'='*60}")
print(f"IS c5_dir A FUNCTION OF lambda? (n={n}, sampled)")
print(f"{'='*60}")

np.random.seed(42)
lambda_groups = {}
for trial in range(50000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    key = tuple(L[i][j] for i in range(n) for j in range(n))
    c5 = count_total_directed_k_cycles(A, n, 5)
    if key not in lambda_groups:
        lambda_groups[key] = set()
    lambda_groups[key].add(c5)

ambiguous = 0
for key, c5_vals in lambda_groups.items():
    if len(c5_vals) > 1:
        ambiguous += 1

print(f"  Total lambda groups: {len(lambda_groups)}")
print(f"  Ambiguous: {ambiguous}")
print(f"  c5 IS a function of lambda? {ambiguous == 0}")

if ambiguous > 0:
    # Show examples
    shown = 0
    for key, c5_vals in lambda_groups.items():
        if len(c5_vals) > 1 and shown < 3:
            print(f"  Example: c5 values = {sorted(c5_vals)}")
            shown += 1

# Also check c7 at n=7
print(f"\n{'='*60}")
print(f"IS c7_dir A FUNCTION OF lambda? (n={n})")
print(f"{'='*60}")

np.random.seed(42)
lambda_groups_c7 = {}
for trial in range(20000):
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
        if len(vals) > 1 and shown < 3:
            print(f"  Example: c7 values = {sorted(vals)}")
            shown += 1

print("\nDone.")
