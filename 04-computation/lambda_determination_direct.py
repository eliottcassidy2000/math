"""
lambda_determination_direct.py -- kind-pasteur-2026-03-13-S61
Direct verification using Vitali pairs (known lambda-preserving pairs).

For each Vitali pair (T, T'), lambda(T) = lambda(T').
Check: is c5_dir(T) = c5_dir(T') always?
       is c7_dir(T) = c7_dir(T') always?

We KNOW dc7 != 0 at n=8, so c7_dir can differ between Vitali pairs.
But we claim c5_dir is always the same.
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

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

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

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

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

# Direct test using Vitali pairs at n=7
for n in [7, 8]:
    total_bits = n * (n-1) // 2
    print("=" * 60)
    print(f"DIRECT LAMBDA DETERMINATION TEST (n={n})")
    print("=" * 60)

    np.random.seed(42)
    c3_diffs = []
    c5_diffs = []
    c7_diffs = []
    examples = 0

    for trial in range(3000 if n <= 8 else 1000):
        bits = np.random.randint(0, 1 << min(total_bits, 31))
        if total_bits > 31:
            bits |= int(np.random.randint(0, 1 << (total_bits - 31))) << 31
        A = bits_to_adj(bits, n)
        lam = lambda_graph(A, n)

        for subset in combinations(range(n), 4):
            ss = sub_scores(A, n, list(subset))
            if ss != (1, 1, 2, 2):
                continue
            B = reverse_subtournament(A, n, list(subset))
            if not np.array_equal(lam, lambda_graph(B, n)):
                continue

            c3_A = count_total_directed_k_cycles(A, n, 3)
            c3_B = count_total_directed_k_cycles(B, n, 3)
            c5_A = count_total_directed_k_cycles(A, n, 5)
            c5_B = count_total_directed_k_cycles(B, n, 5)
            c7_A = count_total_directed_k_cycles(A, n, 7)
            c7_B = count_total_directed_k_cycles(B, n, 7)

            c3_diffs.append(c3_B - c3_A)
            c5_diffs.append(c5_B - c5_A)
            c7_diffs.append(c7_B - c7_A)
            examples += 1
            break

        if trial % 500 == 0 and trial > 0:
            print(f"  ... {trial}, found {examples} Vitali pairs")

    print(f"\n  Total Vitali pairs: {examples}")
    print(f"  dc3 distribution: {dict(sorted(Counter(c3_diffs).items()))}")
    print(f"  dc5 distribution: {dict(sorted(Counter(c5_diffs).items()))}")
    print(f"  dc7 distribution: {dict(sorted(Counter(c7_diffs).items()))}")
    print(f"  dc3 = 0 always? {all(d == 0 for d in c3_diffs)}")
    print(f"  dc5 = 0 always? {all(d == 0 for d in c5_diffs)}")
    print(f"  dc7 = 0 always? {all(d == 0 for d in c7_diffs)}")

print("\nDone.")
