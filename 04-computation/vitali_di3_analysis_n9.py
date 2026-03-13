"""
vitali_di3_analysis_n9.py -- kind-pasteur-2026-03-13-S61
Why is di3 = 0 at n=9?

Three mutually disjoint 3-cycles partition all 9 vertices into three
sets of 3. For di3 to change, the Vitali atom would need to change
the number of such triples.

But the Vitali atom only swaps c3 vertex sets at |V cap S| = 2.
The swap preserves the TOTAL number of c3 vertex sets.

KEY QUESTION: Does the swap also preserve the number of
mutually disjoint c3 triples? If so, di3 is structurally zero.

A disjoint c3 triple {V1, V2, V3} partitions {0,..,8}.
Each Vi has |Vi cap S| vertices from S and |Vi cap E| from ext.
Since |S|=4 and |E|=5: |V1 cap S| + |V2 cap S| + |V3 cap S| = 4.
So the possible S-overlap patterns are:
(0,0,4): impossible since |Vi|=3 but |S|=4
(0,1,3): one Vi fully ext, one with 1 S-vert, one with 3 S-verts
(0,2,2): one Vi fully ext, two with 2 S-verts each
(1,1,2): two with 1 S-vert, one with 2 S-verts
(1,2,1): same as (1,1,2) reordered
(2,2,0): same as (0,2,2) reordered

Wait: |Vi cap S| can be 0 only if all 3 vertices of Vi are in E.
Since |E|=5, a 3-subset of E exists. And |S|=4, so |Vi cap S| can be
0, 1, 2, or 3 (but not 4 since |Vi|=3).

Possible patterns (sorted): (0,1,3), (0,2,2), (1,1,2)
The (0,1,3) pattern has one Vi with 3 S-verts (all of S minus one).
The (0,2,2) pattern has one fully ext Vi and two with 2 S-verts.
The (1,1,2) pattern has no fully ext Vi.

The Vitali atom swaps c3 vertex sets at |Vi cap S| = 2.
So if a triple has a Vi with |Vi cap S| = 2, that Vi might be swapped.
If the swap preserves the triple, di3 is preserved.
If NOT, the triple could be created or destroyed.

Let's check computationally whether di3 is STRUCTURALLY zero or just rare.
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

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def find_c3_vertex_sets(A, n):
    sets = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            sets.append(frozenset(combo))
    return sets

def count_disjoint_triples(c3_sets):
    """Count mutually disjoint c3 triples."""
    n = len(c3_sets)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            if c3_sets[i] & c3_sets[j]:
                continue
            for k in range(j+1, n):
                if not (c3_sets[i] & c3_sets[k]) and not (c3_sets[j] & c3_sets[k]):
                    count += 1
    return count

n = 9
total_bits = n * (n - 1) // 2

print("=" * 70)
print("di3 ANALYSIS AT n=9: IS IT STRUCTURALLY ZERO?")
print("=" * 70)

np.random.seed(2026)
data = []

for trial in range(5000):
    bits_lo = int(np.random.randint(0, 1 << 31))
    bits_hi = int(np.random.randint(0, 1 << (total_bits - 31)))
    bits = bits_lo | (bits_hi << 31)
    A = bits_to_adj(bits, n)
    lam = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(lam, lambda_graph(B, n)):
            continue

        H_A = count_ham_paths(A, n)
        H_B = count_ham_paths(B, n)
        if H_A == H_B:
            continue

        # Count disjoint c3 triples before and after
        c3_A = find_c3_vertex_sets(A, n)
        c3_B = find_c3_vertex_sets(B, n)

        t3_A = count_disjoint_triples(c3_A)
        t3_B = count_disjoint_triples(c3_B)
        delta_t3 = t3_B - t3_A

        # Count c3 vertex set count
        delta_c3_vs = len(c3_B) - len(c3_A)

        # Disjoint c3 pairs
        def count_disjoint_pairs(c3_sets):
            count = 0
            for i in range(len(c3_sets)):
                for j in range(i+1, len(c3_sets)):
                    if not (c3_sets[i] & c3_sets[j]):
                        count += 1
            return count

        dp_A = count_disjoint_pairs(c3_A)
        dp_B = count_disjoint_pairs(c3_B)
        delta_dp = dp_B - dp_A

        # Lost/gained c3 and their S-overlap
        S = frozenset(subset)
        lost = set(c3_A) - set(c3_B)
        gained = set(c3_B) - set(c3_A)

        data.append({
            'delta_H': H_B - H_A,
            't3_A': t3_A, 't3_B': t3_B, 'delta_t3': delta_t3,
            'dp_A': dp_A, 'dp_B': dp_B, 'delta_dp': delta_dp,
            'c3_A': len(c3_A), 'delta_c3_vs': delta_c3_vs,
            'lost': len(lost), 'gained': len(gained),
        })
        break

    if trial % 1000 == 0 and trial > 0:
        print(f"  ... {trial}/5000, found {len(data)}")

print(f"\nTotal H-changing examples: {len(data)}")

# Analysis
print(f"\n  delta_t3 distribution: {dict(sorted(Counter(d['delta_t3'] for d in data).items()))}")
print(f"  delta_dp (c3-c3 disjoint pairs): {dict(sorted(Counter(d['delta_dp'] for d in data).items()))}")
print(f"  delta_c3_vs: {dict(sorted(Counter(d['delta_c3_vs'] for d in data).items()))}")

# Is di3 (the DIRECTED cycle version) also zero?
# Since c3 always has multiplicity 1, i3_directed = i3_vertex_set = t3
# So if delta_t3 = 0, then delta_i3 = 0 too.

t3_zero = all(d['delta_t3'] == 0 for d in data)
dp_zero = all(d['delta_dp'] == 0 for d in data)
c3_vs_zero = all(d['delta_c3_vs'] == 0 for d in data)

print(f"\n  delta_t3 = 0 always? {t3_zero}")
print(f"  delta_dp = 0 always? {dp_zero}")
print(f"  delta_c3_vs = 0 always? {c3_vs_zero}")

if not t3_zero:
    print(f"\n  di3 CAN change at n=9! Examples:")
    for d in data:
        if d['delta_t3'] != 0:
            print(f"    dH={d['delta_H']}, dt3={d['delta_t3']}, t3_A={d['t3_A']}, t3_B={d['t3_B']}")
else:
    print(f"\n  di3 = 0 STRUCTURALLY at n=9!")
    print(f"  This means the c3 swap preserves disjoint triple counts.")

# Show absolute values
print(f"\n  t3_A values: {dict(sorted(Counter(d['t3_A'] for d in data).items()))}")
print(f"  dp_A values (c3-c3 pairs): min={min(d['dp_A'] for d in data)}, max={max(d['dp_A'] for d in data)}")

# Cross-tab: di3 vs delta_H
if not t3_zero:
    print(f"\n  delta_H vs delta_t3:")
    for d in data:
        if d['delta_t3'] != 0:
            print(f"    dH={d['delta_H']}, dt3={d['delta_t3']}, ddp={d['delta_dp']}")

print("\nDone.")
