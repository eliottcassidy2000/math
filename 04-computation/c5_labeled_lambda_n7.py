#!/usr/bin/env python3
"""
Check: is c5_dir determined by the LABELED lambda graph at n=7?

The distinction:
- Lambda MULTISET: sorted tuple of all C(n,2) lambda values
- Lambda LABELED: the full matrix L[i][j] for all pairs (depends on vertex labeling)

kind-pasteur claims c5 is lambda-determined at n=7.
Our multiset test shows 35 ambiguous groups.
The resolution might be that c5 depends on the LABELED lambda, not just the multiset.

opus-2026-03-13-S71c
"""
import sys, time
import numpy as np
from itertools import combinations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

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
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

n = 7
tb = n*(n-1)//2
np.random.seed(42)

# Test labeled lambda
lam_groups_labeled = defaultdict(set)
# Test multiset lambda
lam_groups_multi = defaultdict(set)

t0 = time.time()
for trial in range(50000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    # Labeled key: upper triangle in row-major order
    labeled_key = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
    # Multiset key: sorted
    multi_key = tuple(sorted(labeled_key))

    # c5_dir: sum over all 5-subsets
    c5_total = 0
    for combo in combinations(range(n), 5):
        A_sub = A[np.ix_(list(combo), list(combo))]
        tr5_sub = int(np.trace(np.linalg.matrix_power(A_sub, 5)))
        c5_total += tr5_sub // 5

    lam_groups_labeled[labeled_key].add(c5_total)
    lam_groups_multi[multi_key].add(c5_total)

dt = time.time() - t0
print(f"n=7, 50000 samples, {dt:.1f}s")

ambig_labeled = sum(1 for v in lam_groups_labeled.values() if len(v) > 1)
ambig_multi = sum(1 for v in lam_groups_multi.values() if len(v) > 1)

print(f"\nLabeled lambda groups: {len(lam_groups_labeled)}, ambiguous: {ambig_labeled}")
print(f"Multiset lambda groups: {len(lam_groups_multi)}, ambiguous: {ambig_multi}")

if ambig_labeled > 0:
    print("\n  Labeled ambiguities:")
    count = 0
    for key, vals in sorted(lam_groups_labeled.items()):
        if len(vals) > 1:
            print(f"    labeled={key}: c5 = {sorted(vals)}")
            count += 1
            if count >= 5: break
else:
    print("\n  c5 IS determined by LABELED lambda! Kind-pasteur was right.")
    print("  c5 is NOT determined by lambda multiset alone — vertex identity matters.")

print(f"\nDone.")
