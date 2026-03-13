"""
vitali_n9_prediction.py -- kind-pasteur-2026-03-13-S61
Predictions for Vitali atom at n=9 and first tests.

PROVEN HIERARCHY:
n<=6: delta_H = 0 (Vitali is gauge-trivial)
n=7:  delta_H = 2*dc7_dir (one channel: c7 count)
n=8:  delta_H = 2*dc7_dir + 4*di2 (two channels: c7 + disjoint pairs)

PREDICTIONS FOR n=9:
1. c9 channel opens: 9-cycles use all 9 vertices
   => delta_c9_dir contributes 2*dc9 (like c7 at n=7)

2. i2 grows: more room for disjoint pairs
   - (c3,c3): need 6 <= 9. More pairs possible.
   - (c3,c5): need 8 <= 9. Still tight but more room.
   - (c3,c7): need 10 > 9. Still impossible.
   - (c5,c5): need 10 > 9. Still impossible.

3. i3 channel opens: disjoint TRIPLES
   - (c3,c3,c3): need 9 = n. Tight! EXISTS.
   - Others: need > 9. Impossible.
   - Contributes 8*delta_i3

4. c5 may become non-invariant (phase transition?)

FORMULA PREDICTION:
delta_H = 2*(dc7_dir + dc9_dir) + 4*di2 + 8*di3

where di3 counts changes in vertex-disjoint directed cycle triples.

This script tests these predictions at n=9 (small sample due to cost).
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
    """Count vertex-disjoint directed cycle pairs (with multiplicity)."""
    vsets = list(cycle_dict.keys())
    total = 0
    for i in range(len(vsets)):
        for j in range(i+1, len(vsets)):
            if len(vsets[i] & vsets[j]) == 0:
                total += cycle_dict[vsets[i]] * cycle_dict[vsets[j]]
    return total

def count_disjoint_directed_triples(cycle_dict):
    """Count vertex-disjoint directed cycle triples (with multiplicity)."""
    vsets = list(cycle_dict.keys())
    total = 0
    nc = len(vsets)
    for i in range(nc):
        for j in range(i+1, nc):
            if vsets[i] & vsets[j]:
                continue
            for k in range(j+1, nc):
                if not (vsets[i] & vsets[k]) and not (vsets[j] & vsets[k]):
                    total += cycle_dict[vsets[i]] * cycle_dict[vsets[j]] * cycle_dict[vsets[k]]
    return total

n = 9
total_bits = n * (n - 1) // 2

print("=" * 70)
print(f"VITALI ATOM AT n={n}: TESTING PREDICTIONS")
print("=" * 70)

np.random.seed(2026)
data = []

t0 = time.time()
for trial in range(5000):  # need more trials at n=9
    bits = int(np.random.randint(0, 1 << 31)) | (int(np.random.randint(0, 1 << (total_bits - 31))) << 31)
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

        # Compute cycle data — include 9-cycles (C(9,9)=1 set, 8!=40320 perms, fast)
        all_cycles_A = {}
        all_cycles_B = {}
        for length in [3, 5, 7, 9]:  # include all odd lengths
            c_A = count_directed_cycles_by_set(A, n, length)
            c_B = count_directed_cycles_by_set(B, n, length)
            for vs, m in c_A.items():
                all_cycles_A[vs] = (length, m)
            for vs, m in c_B.items():
                all_cycles_B[vs] = (length, m)

        # Net directed cycle counts by length
        def total_dir_by_len(cycles):
            by_len = Counter()
            for vs, (l, m) in cycles.items():
                by_len[l] += m
            return by_len

        dir_A = total_dir_by_len(all_cycles_A)
        dir_B = total_dir_by_len(all_cycles_B)

        dc3 = dir_B.get(3, 0) - dir_A.get(3, 0)
        dc5 = dir_B.get(5, 0) - dir_A.get(5, 0)
        dc7 = dir_B.get(7, 0) - dir_A.get(7, 0)
        dc9 = dir_B.get(9, 0) - dir_A.get(9, 0)

        # Disjoint pair count (directed, with multiplicity)
        cycle_dict_A = {vs: m for vs, (l, m) in all_cycles_A.items()}
        cycle_dict_B = {vs: m for vs, (l, m) in all_cycles_B.items()}

        i2_A = count_disjoint_directed_pairs(cycle_dict_A)
        i2_B = count_disjoint_directed_pairs(cycle_dict_B)
        di2 = i2_B - i2_A

        # Disjoint triple count
        i3_A = count_disjoint_directed_triples(cycle_dict_A)
        i3_B = count_disjoint_directed_triples(cycle_dict_B)
        di3 = i3_B - i3_A

        delta_H = H_B - H_A

        # Test prediction
        predicted = 2 * (dc7 + dc9) + 4 * di2 + 8 * di3

        data.append({
            'delta_H': delta_H,
            'dc3': dc3, 'dc5': dc5, 'dc7': dc7, 'dc9': dc9,
            'di2': di2, 'di3': di3,
            'predicted': predicted,
            'i2_A': i2_A, 'i2_B': i2_B,
            'i3_A': i3_A, 'i3_B': i3_B,
        })

        if predicted != delta_H:
            print(f"  MISMATCH at trial {trial}: dH={delta_H}, pred={predicted}")
            print(f"    dc3={dc3}, dc5={dc5}, dc7={dc7}, dc9={dc9}")
            print(f"    di2={di2}, di3={di3}")

        found = True
        break

    if trial % 500 == 0 and trial > 0:
        elapsed = time.time() - t0
        print(f"  ... {trial}/200, found {len(data)}, elapsed {elapsed:.0f}s")

elapsed = time.time() - t0
print(f"\nTotal: {len(data)} H-changing examples in {elapsed:.0f}s")

if not data:
    print("No examples found!")
else:
    # Formula verification
    matches = sum(1 for d in data if d['predicted'] == d['delta_H'])
    print(f"\nFormula delta_H = 2*(dc7+dc9) + 4*di2 + 8*di3:")
    print(f"  Match rate: {matches}/{len(data)} ({100*matches/len(data):.1f}%)")

    # Channel analysis
    print(f"\n  dc3 distribution: {dict(sorted(Counter(d['dc3'] for d in data).items()))}")
    print(f"  dc5 distribution: {dict(sorted(Counter(d['dc5'] for d in data).items()))}")
    print(f"  dc7 distribution: {dict(sorted(Counter(d['dc7'] for d in data).items()))}")
    print(f"  dc9 distribution: {dict(sorted(Counter(d['dc9'] for d in data).items()))}")
    print(f"  di2 distribution: {dict(sorted(Counter(d['di2'] for d in data).items()))}")
    print(f"  di3 distribution: {dict(sorted(Counter(d['di3'] for d in data).items()))}")
    print(f"  delta_H distribution: {dict(sorted(Counter(d['delta_H'] for d in data).items()))}")

    # Is dc3 = 0 still? (c3 net preservation)
    dc3_zero = all(d['dc3'] == 0 for d in data)
    print(f"\n  dc3 = 0 always? {dc3_zero}")

    # Is dc5 = 0 still? (c5 net preservation)
    dc5_zero = all(d['dc5'] == 0 for d in data)
    print(f"  dc5 = 0 always? {dc5_zero}")

    # Does c9 channel now contribute?
    dc9_nonzero = [d for d in data if d['dc9'] != 0]
    print(f"  dc9 != 0: {len(dc9_nonzero)} / {len(data)}")

    # Does i3 channel now contribute?
    di3_nonzero = [d for d in data if d['di3'] != 0]
    print(f"  di3 != 0: {len(di3_nonzero)} / {len(data)}")

    # Show a few examples
    print(f"\n  Sample examples:")
    for d in data[:10]:
        match = "OK" if d['predicted'] == d['delta_H'] else "FAIL"
        print(f"    dH={d['delta_H']:>4}, dc7={d['dc7']:>3}, dc9={d['dc9']:>3}, di2={d['di2']:>3}, di3={d['di3']:>3}, pred={d['predicted']:>4} [{match}]")

print("\nDone.")
