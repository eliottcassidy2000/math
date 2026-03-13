"""
vitali_n10_pair_types.py -- kind-pasteur-2026-03-13-S61
At n=10: new disjoint pair types are geometrically possible:
  (c3,c7): 3+7=10=n  (FIRST TIME)
  (c5,c5): 5+5=10=n  (FIRST TIME)

The formula delta_H = 2*(dc7+dc9) + 4*di2 works at 100% (12/12).
But WHICH pair types drive di2?

Also: c11 impossible (11>10), so no new dc channel.
The question: do the new pair types contribute to di2?
If so, the mechanism is richer than at n=8-9 where only (c3,c5) contributed.
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

def count_disjoint_pairs_by_type(cycle_dicts_A, cycle_dicts_B, lengths):
    """Count disjoint directed pair changes broken down by pair type (len1, len2)."""
    results = {}
    for i, l1 in enumerate(lengths):
        for j, l2 in enumerate(lengths):
            if l1 > l2:
                continue
            if l1 + l2 > 10:  # can't be disjoint at n=10
                continue

            # Compute i2 for this pair type before and after
            vsets_A_1 = list(cycle_dicts_A.get(l1, {}).keys())
            vsets_A_2 = list(cycle_dicts_A.get(l2, {}).keys())
            vsets_B_1 = list(cycle_dicts_B.get(l1, {}).keys())
            vsets_B_2 = list(cycle_dicts_B.get(l2, {}).keys())

            def count_type(d1, d2, same_len):
                total = 0
                vs1 = list(d1.keys())
                vs2 = list(d2.keys())
                for a in range(len(vs1)):
                    start_b = a + 1 if same_len else 0
                    for b in range(start_b if same_len else 0, len(vs2)):
                        if same_len and b <= a:
                            continue
                        if len(vs1[a] & vs2[b]) == 0:
                            total += d1[vs1[a]] * d2[vs2[b]]
                return total

            same = (l1 == l2)
            d1_A = cycle_dicts_A.get(l1, {})
            d2_A = cycle_dicts_A.get(l2, {})
            d1_B = cycle_dicts_B.get(l1, {})
            d2_B = cycle_dicts_B.get(l2, {})

            i2_A = count_type(d1_A, d2_A, same)
            i2_B = count_type(d1_B, d2_B, same)

            results[(l1, l2)] = {'i2_A': i2_A, 'i2_B': i2_B, 'delta': i2_B - i2_A}

    return results

n = 10
total_bits = n * (n - 1) // 2

print("=" * 70)
print(f"VITALI ATOM AT n={n}: DISJOINT PAIR TYPE DECOMPOSITION")
print("=" * 70)

np.random.seed(2026)
data = []

t0 = time.time()
for trial in range(3000):
    bits = 0
    for chunk in range(0, total_bits, 30):
        chunk_size = min(30, total_bits - chunk)
        bits |= int(np.random.randint(0, 1 << chunk_size)) << chunk

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

        # Compute cycle dicts by length
        lengths = [3, 5, 7, 9]
        cd_A = {}
        cd_B = {}
        dc = {}
        for l in lengths:
            c_A = count_directed_cycles_by_set(A, n, l)
            c_B = count_directed_cycles_by_set(B, n, l)
            cd_A[l] = c_A
            cd_B[l] = c_B
            dc[l] = sum(c_B.values()) - sum(c_A.values())

        # Pair type decomposition
        pair_types = count_disjoint_pairs_by_type(cd_A, cd_B, lengths)

        delta_H = H_B - H_A
        total_di2 = sum(v['delta'] for v in pair_types.values())

        data.append({
            'delta_H': delta_H,
            'dc': dc,
            'pair_types': pair_types,
            'total_di2': total_di2,
        })
        break

    if trial % 500 == 0 and trial > 0:
        elapsed = time.time() - t0
        print(f"  ... {trial}/3000, found {len(data)}, elapsed {elapsed:.0f}s")

elapsed = time.time() - t0
print(f"\nTotal: {len(data)} H-changing examples in {elapsed:.0f}s")

if not data:
    print("No examples found!")
else:
    # Pair type analysis
    print(f"\n{'='*70}")
    print(f"DISJOINT PAIR TYPE CONTRIBUTIONS TO di2")
    print(f"{'='*70}")

    # Possible pair types at n=10: (3,3), (3,5), (3,7), (5,5)
    for pt in [(3,3), (3,5), (3,7), (5,5)]:
        deltas = [d['pair_types'].get(pt, {}).get('delta', 0) for d in data]
        nonzero = sum(1 for x in deltas if x != 0)
        print(f"\n  ({pt[0]},{pt[1]}) pairs:")
        print(f"    Nonzero: {nonzero}/{len(data)}")
        print(f"    Distribution: {dict(sorted(Counter(deltas).items()))}")
        if any(d['pair_types'].get(pt, {}).get('i2_A', 0) > 0 for d in data):
            i2_vals = [d['pair_types'].get(pt, {}).get('i2_A', 0) for d in data]
            print(f"    i2_A values: min={min(i2_vals)}, max={max(i2_vals)}, mean={np.mean(i2_vals):.1f}")

    # Total di2 check
    matches = sum(1 for d in data if 2*(d['dc'][7]+d['dc'][9]) + 4*d['total_di2'] == d['delta_H'])
    print(f"\n  Formula check: {matches}/{len(data)}")

    # Which pair types drive di2?
    print(f"\n  Pair type contribution to total |di2|:")
    for pt in [(3,3), (3,5), (3,7), (5,5)]:
        total_abs = sum(abs(d['pair_types'].get(pt, {}).get('delta', 0)) for d in data)
        total_abs_all = sum(abs(d['total_di2']) for d in data)
        pct = 100 * total_abs / total_abs_all if total_abs_all > 0 else 0
        print(f"    ({pt[0]},{pt[1]}): {total_abs} ({pct:.1f}%)")

print("\nDone.")
