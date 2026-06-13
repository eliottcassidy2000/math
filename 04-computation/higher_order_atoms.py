"""
higher_order_atoms.py -- kind-pasteur-2026-03-13-S61

The simplex profile ambiguity comes in three types:
  - dc7 = 1: Connected by Vitali atom (4-vertex reversal preserving lambda)
  - dc7 = 2: Connected by ??? (5-vertex reversal? composition?)
  - dc7 = 3: Connected by ??? (6-vertex reversal? triple composition?)

This script investigates:
1. Do higher-order lambda-preserving reversals exist (5-vertex, 6-vertex)?
2. Are the dc7 = 2 cases related by composition of 2 Vitali atoms?
3. What is the EXACT mechanism for each dc7 gap?
4. Is there a "generalized Vitali atom" at each size?
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict

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

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def get_simplex_data(A, n):
    L = lambda_graph(A, n)
    A2 = A @ A
    pairs = []
    for u in range(n):
        for v in range(u+1, n):
            sig = n - 2 - int(A2[u][v]) - int(A2[v][u])
            lam = int(L[u][v])
            delta = n - 2 - sig - lam
            pairs.append((sig, lam, delta))
    return pairs

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print("HIGHER-ORDER LAMBDA-PRESERVING REVERSALS AT n=7")
print("=" * 60)

np.random.seed(42)

# 1. Search for k-vertex lambda-preserving reversals for k=4,5,6
print("\n--- Lambda-preserving reversal census ---")

np.random.seed(42)
reversal_stats = {k: {'found': 0, 'total': 0, 'dc7_dist': Counter()}
                  for k in [4, 5, 6]}

for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    c7_A = count_directed_k_cycles(A, n, 7)

    for k in [4, 5, 6]:
        for subset in combinations(range(n), k):
            reversal_stats[k]['total'] += 1
            # Reverse the sub-tournament
            B = A.copy()
            for i in subset:
                for j in subset:
                    if i != j:
                        B[i][j] = A[j][i]
            L2 = lambda_graph(B, n)
            if np.array_equal(L, L2):
                reversal_stats[k]['found'] += 1
                c7_B = count_directed_k_cycles(B, n, 7)
                dc7 = c7_B - c7_A
                reversal_stats[k]['dc7_dist'][dc7] += 1

for k in [4, 5, 6]:
    stats = reversal_stats[k]
    rate = 100 * stats['found'] / stats['total'] if stats['total'] > 0 else 0
    print(f"\n  k={k}: {stats['found']}/{stats['total']} lambda-preserving ({rate:.2f}%)")
    print(f"    dc7 distribution: {dict(sorted(stats['dc7_dist'].items()))}")

# 2. For the dc7=2 case: is it a composition of two Vitali atoms?
print(f"\n\n{'='*60}")
print("dc7=2 VIA COMPOSITION OF VITALI ATOMS")
print("=" * 60)

np.random.seed(42)
dc7_2_examples = []

for trial in range(3000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    c7_A = count_directed_k_cycles(A, n, 7)

    # Find all Vitali atoms (4-vertex, lambda-preserving)
    atoms = []
    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = A.copy()
        for i in subset:
            for j in subset:
                if i != j:
                    B[i][j] = A[j][i]
        if np.array_equal(L, lambda_graph(B, n)):
            c7_B = count_directed_k_cycles(B, n, 7)
            atoms.append({'subset': subset, 'B': B.copy(), 'dc7': c7_B - c7_A, 'c7_B': c7_B})

    # Try all pairs of atoms to see if composition gives dc7=2
    for i, a1 in enumerate(atoms):
        for j, a2 in enumerate(atoms):
            if i == j:
                continue
            # Apply a1 then a2
            B1 = a1['B']
            L_B1 = lambda_graph(B1, n)
            # Now apply a2's subset to B1
            C = B1.copy()
            for u in a2['subset']:
                for v in a2['subset']:
                    if u != v:
                        C[u][v] = B1[v][u]
            # Check if lambda still preserved
            L_C = lambda_graph(C, n)
            if np.array_equal(L, L_C):
                c7_C = count_directed_k_cycles(C, n, 7)
                dc7_total = c7_C - c7_A
                if abs(dc7_total) >= 2:
                    dc7_2_examples.append({
                        'c7_A': c7_A,
                        'dc7_1': a1['dc7'],
                        'dc7_2': c7_C - a1['c7_B'],
                        'dc7_total': dc7_total,
                        'atoms': (a1['subset'], a2['subset']),
                        'overlap': len(set(a1['subset']) & set(a2['subset'])),
                    })

    if len(dc7_2_examples) >= 20:
        break

print(f"Found {len(dc7_2_examples)} cases with |dc7_composition| >= 2")
for ex in dc7_2_examples[:10]:
    print(f"  Atoms {ex['atoms']}, overlap={ex['overlap']}: "
          f"dc7_1={ex['dc7_1']:+d}, dc7_2={ex['dc7_2']:+d}, total={ex['dc7_total']:+d}")

# 3. Check if dc7=2 cases from simplex ambiguity are related by
#    5-vertex lambda-preserving reversal
print(f"\n\n{'='*60}")
print("dc7=2 VIA 5-VERTEX REVERSAL")
print("=" * 60)

np.random.seed(42)
dc7_2_by_5 = []

for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    c7_A = count_directed_k_cycles(A, n, 7)

    for subset in combinations(range(n), 5):
        B = A.copy()
        for i in subset:
            for j in subset:
                if i != j:
                    B[i][j] = A[j][i]
        if np.array_equal(L, lambda_graph(B, n)):
            c7_B = count_directed_k_cycles(B, n, 7)
            dc7 = c7_B - c7_A
            if dc7 != 0:
                dc7_2_by_5.append({
                    'c7_A': c7_A, 'c7_B': c7_B, 'dc7': dc7,
                    'subset': subset,
                    'sub_scores': sub_scores(A, n, list(subset)),
                })

print(f"Found {len(dc7_2_by_5)} dc7!=0 from 5-vertex reversals")
dc7_5_dist = Counter(ex['dc7'] for ex in dc7_2_by_5)
print(f"dc7 distribution: {dict(sorted(dc7_5_dist.items()))}")
score_dist = Counter(ex['sub_scores'] for ex in dc7_2_by_5)
print(f"Sub-tournament scores: {dict(sorted(score_dist.items()))}")

for ex in dc7_2_by_5[:5]:
    print(f"  subset={ex['subset']}, scores={ex['sub_scores']}, "
          f"c7: {ex['c7_A']}->{ex['c7_B']} (dc7={ex['dc7']:+d})")

# 4. Check 6-vertex
print(f"\n\n{'='*60}")
print("6-VERTEX REVERSAL")
print("=" * 60)

np.random.seed(42)
dc7_6 = []

for trial in range(1000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    c7_A = count_directed_k_cycles(A, n, 7)

    for subset in combinations(range(n), 6):
        B = A.copy()
        for i in subset:
            for j in subset:
                if i != j:
                    B[i][j] = A[j][i]
        if np.array_equal(L, lambda_graph(B, n)):
            c7_B = count_directed_k_cycles(B, n, 7)
            dc7 = c7_B - c7_A
            dc7_6.append({
                'dc7': dc7,
                'sub_scores': sub_scores(A, n, list(subset)),
            })

print(f"Found {len(dc7_6)} 6-vertex lambda-preserving reversals")
dc7_6_dist = Counter(ex['dc7'] for ex in dc7_6)
print(f"dc7 distribution: {dict(sorted(dc7_6_dist.items()))}")
score6_dist = Counter(ex['sub_scores'] for ex in dc7_6)
print(f"Sub-tournament scores: {dict(sorted(score6_dist.items()))}")

# 5. Summary: the "atom hierarchy"
print(f"\n\n{'='*60}")
print("ATOM HIERARCHY SUMMARY")
print("=" * 60)
print("""
At n=7, lambda-preserving reversals exist at sizes:
  k=4: The Vitali atom. Score (1,1,2,2). |dc7| <= 1.
  k=5: ??? Score type ???. |dc7| range ???.
  k=6: ??? |dc7| range ???.

The simplex profile ambiguity has gaps 1, 2, 3.
These correspond to different "atom levels" in the hierarchy.
""")

print("\nDone.")
