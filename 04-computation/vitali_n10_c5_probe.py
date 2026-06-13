"""
vitali_n10_c5_probe.py -- kind-pasteur-2026-03-13-S61
Does c5 become non-invariant at n=10?

At n<=9: dc5_dir = 0 ALWAYS (the Vitali atom preserves net c5 count).
At n=10: a 5-cycle V has |V ∩ S| in {0,1,2,3,4}.
  |V ∩ S| = 0: possible since |V|+|S|=9 <= 10. V fully external.
  |V ∩ S| = 1: possible. 1 arc touches S. No reversed arcs in cycle. INVARIANT.
  |V ∩ S| = 2: 1 reversed arc within S. CAN change direction count.
  |V ∩ S| = 3: 3 reversed arcs. CAN change.
  |V ∩ S| = 4: 6 reversed arcs. CAN change.

Key difference from n<=9: at n<=9, |V ∩ S| >= |V|+|S|-n = 9-n.
  n=8: |V ∩ S| >= 1. So |V ∩ S|=0 is impossible.
  n=9: |V ∩ S| >= 0. |V ∩ S|=0 is possible (V is 5 of 5 ext vertices).
  n=10: |V ∩ S| >= 0. Same, but now 6 ext vertices (more room).

The c5 NET preservation (dc5=0) is NOT just about arcs:
it requires that gained-minus-lost directed c5's cancel.
The question: does MORE room at n=10 allow this to break?

Also check: does di2 now include (c5,c5) pairs? Need 10 = n vertices.
And (c3,c7) pairs: 3+7=10 = n. YES, new pair types!

Test the full formula at n=10 (computationally expensive — use small sample).
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter
import time

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

def count_directed_cycles_by_set(A, n, length):
    """Return dict: vertex_set -> directed_cycle_count."""
    result = {}
    for combo in combinations(range(n), length):
        count = 0
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = True
            for k in range(length):
                if A[path[k]][path[(k+1) % length]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
        if count > 0:
            result[frozenset(combo)] = count
    return result

def count_disjoint_directed_pairs(cycle_dict):
    vsets = list(cycle_dict.keys())
    total = 0
    for i in range(len(vsets)):
        for j in range(i+1, len(vsets)):
            if len(vsets[i] & vsets[j]) == 0:
                total += cycle_dict[vsets[i]] * cycle_dict[vsets[j]]
    return total

n = 10
total_bits = n * (n - 1) // 2  # 45

print("=" * 70)
print(f"VITALI ATOM AT n={n}: c5 PHASE TRANSITION PROBE")
print(f"total_bits = {total_bits}")
print("=" * 70)

np.random.seed(2026)
data = []

t0 = time.time()
for trial in range(2000):
    # Generate random tournament bits (45 bits)
    bits = 0
    for chunk in range(0, total_bits, 30):
        chunk_size = min(30, total_bits - chunk)
        bits |= int(np.random.randint(0, 1 << chunk_size)) << chunk

    A = bits_to_adj(bits, n)
    lam = lambda_graph(A, n)

    found = False
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

        # Compute cycle data — c3, c5, c7, c9 (skip c10+ since not odd)
        # At n=10: odd cycles up to 9
        dc = {}
        all_cycles_A = {}
        all_cycles_B = {}

        for length in [3, 5, 7, 9]:
            c_A = count_directed_cycles_by_set(A, n, length)
            c_B = count_directed_cycles_by_set(B, n, length)
            dc[length] = sum(c_B.get(vs, 0) for vs in set(c_A) | set(c_B)) - sum(c_A.get(vs, 0) for vs in set(c_A) | set(c_B))

            # Actually just sum totals
            tot_A = sum(c_A.values())
            tot_B = sum(c_B.values())
            dc[length] = tot_B - tot_A

            for vs, m in c_A.items():
                all_cycles_A[vs] = m
            for vs, m in c_B.items():
                all_cycles_B[vs] = m

        # Disjoint pair count (directed, with multiplicity)
        i2_A = count_disjoint_directed_pairs(all_cycles_A)
        i2_B = count_disjoint_directed_pairs(all_cycles_B)
        di2 = i2_B - i2_A

        delta_H = H_B - H_A

        # Test formula: delta_H = 2*(dc7+dc9) + 4*di2
        predicted = 2 * (dc[7] + dc[9]) + 4 * di2

        data.append({
            'delta_H': delta_H,
            'dc3': dc[3], 'dc5': dc[5], 'dc7': dc[7], 'dc9': dc[9],
            'di2': di2,
            'predicted': predicted,
        })

        if dc[5] != 0:
            print(f"  *** c5 NON-ZERO at trial {trial}: dc5={dc[5]}, dH={delta_H}")

        found = True
        break

    if trial % 200 == 0 and trial > 0:
        elapsed = time.time() - t0
        print(f"  ... {trial}/2000, found {len(data)}, elapsed {elapsed:.0f}s")

elapsed = time.time() - t0
print(f"\nTotal: {len(data)} H-changing examples in {elapsed:.0f}s")

if not data:
    print("No examples found!")
else:
    matches = sum(1 for d in data if d['predicted'] == d['delta_H'])
    print(f"\nFormula delta_H = 2*(dc7+dc9) + 4*di2:")
    print(f"  Match rate: {matches}/{len(data)} ({100*matches/len(data):.1f}%)")

    print(f"\n  dc3 distribution: {dict(sorted(Counter(d['dc3'] for d in data).items()))}")
    print(f"  dc5 distribution: {dict(sorted(Counter(d['dc5'] for d in data).items()))}")
    print(f"  dc7 distribution: {dict(sorted(Counter(d['dc7'] for d in data).items()))}")
    print(f"  dc9 distribution: {dict(sorted(Counter(d['dc9'] for d in data).items()))}")
    print(f"  di2 distribution: {dict(sorted(Counter(d['di2'] for d in data).items()))}")
    print(f"  delta_H distribution: {dict(sorted(Counter(d['delta_H'] for d in data).items()))}")

    dc3_zero = all(d['dc3'] == 0 for d in data)
    dc5_zero = all(d['dc5'] == 0 for d in data)
    print(f"\n  dc3 = 0 always? {dc3_zero}")
    print(f"  dc5 = 0 always? {dc5_zero}")

    # Mismatches
    mismatches = [d for d in data if d['predicted'] != d['delta_H']]
    if mismatches:
        print(f"\n  MISMATCHES ({len(mismatches)}):")
        for d in mismatches[:10]:
            residual = d['delta_H'] - d['predicted']
            print(f"    dH={d['delta_H']}, pred={d['predicted']}, residual={residual}")
            print(f"      dc3={d['dc3']}, dc5={d['dc5']}, dc7={d['dc7']}, dc9={d['dc9']}, di2={d['di2']}")

print("\nDone.")
