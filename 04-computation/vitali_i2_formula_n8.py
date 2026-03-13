"""
vitali_i2_formula_n8.py -- kind-pasteur-2026-03-13-S61
VERIFY: delta_H = 2*delta_c7_dir + 4*delta_i2 at n=8

From directed_multiplicity_vitali.py:
- dc3_dir = dc5_dir = 0 always (net directed cycle changes cancel)
- residual = dH - 2*dc7_dir is always in {-4, 0, 4}
- This is EXACTLY the signature of delta_i2 (disjoint pair count)

The clique perturbation argument:
- c7 cycles use 7 of 8 vertices => universal conflict
  => only contributes 2*delta_c7_dir (alpha_1 only)
- c3/c5 reshuffling: net zero directed cycles, but WHICH cycles exist changes
  => this changes the NUMBER of vertex-disjoint directed cycle pairs
  => contributes 4*delta_i2

PREDICTION: delta_H = 2*delta_c7_dir + 4*delta_i2 EXACTLY for ALL examples
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

def enumerate_all_directed_cycles(A, n):
    """Enumerate all directed odd cycles."""
    cycles = []
    for k in range(3, n+1, 2):
        for combo in combinations(range(n), k):
            for perm in permutations(combo[1:]):
                path = (combo[0],) + perm
                valid = True
                for i in range(k):
                    if A[path[i]][path[(i+1) % k]] != 1:
                        valid = False
                        break
                if valid:
                    cycles.append(path)
    return cycles

def count_disjoint_pairs(cycles):
    """Count pairs of directed cycles with disjoint vertex sets."""
    vsets = [frozenset(c) for c in cycles]
    count = 0
    for i in range(len(vsets)):
        for j in range(i+1, len(vsets)):
            if len(vsets[i] & vsets[j]) == 0:
                count += 1
    return count

n = 8
total_bits = n * (n - 1) // 2

print("=" * 70)
print("VERIFY: delta_H = 2*delta_c7_dir + 4*delta_i2")
print("=" * 70)

verified = 0
failed = 0
data = []

for seed in range(20):
    np.random.seed(seed * 1000 + 42)
    for trial in range(500):
        bits = np.random.randint(0, 1 << total_bits)
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

            # Get ALL directed cycles
            dir_A = enumerate_all_directed_cycles(A, n)
            dir_B = enumerate_all_directed_cycles(B, n)

            # c7 directed cycle count (total, not per vertex set)
            c7_dir_A = sum(1 for c in dir_A if len(c) == 7)
            c7_dir_B = sum(1 for c in dir_B if len(c) == 7)
            delta_c7_dir = c7_dir_B - c7_dir_A

            # Disjoint directed cycle pairs
            i2_A = count_disjoint_pairs(dir_A)
            i2_B = count_disjoint_pairs(dir_B)
            delta_i2 = i2_B - i2_A

            delta_H = H_B - H_A
            predicted = 2 * delta_c7_dir + 4 * delta_i2

            if predicted == delta_H:
                verified += 1
            else:
                failed += 1
                print(f"  FAIL: dH={delta_H}, dc7_dir={delta_c7_dir}, di2={delta_i2}, predicted={predicted}")

            data.append({
                'delta_H': delta_H,
                'delta_c7_dir': delta_c7_dir,
                'delta_i2': delta_i2,
                'predicted': predicted,
                'i2_A': i2_A,
                'i2_B': i2_B,
                'total_dir_A': len(dir_A),
                'total_dir_B': len(dir_B),
            })

            break

    print(f"  seed {seed}: verified={verified}, failed={failed}")

print(f"\nRESULT: verified={verified}, failed={failed}")
print(f"Formula delta_H = 2*delta_c7_dir + 4*delta_i2: {'CONFIRMED' if failed == 0 else 'FAILS'}")

# Distribution analysis
print(f"\n  delta_c7_dir distribution: {dict(sorted(Counter(d['delta_c7_dir'] for d in data).items()))}")
print(f"  delta_i2 distribution: {dict(sorted(Counter(d['delta_i2'] for d in data).items()))}")
print(f"  delta_H distribution: {dict(sorted(Counter(d['delta_H'] for d in data).items()))}")

# Cross-tab
print(f"\n  Cross-tab (dc7_dir, di2) -> delta_H:")
cross = Counter()
for d in data:
    cross[(d['delta_c7_dir'], d['delta_i2'])] = cross.get((d['delta_c7_dir'], d['delta_i2']), [])
    if isinstance(cross[(d['delta_c7_dir'], d['delta_i2'])], list):
        cross[(d['delta_c7_dir'], d['delta_i2'])].append(d['delta_H'])
    else:
        cross[(d['delta_c7_dir'], d['delta_i2'])] += 1

# Fix: use defaultdict
cross2 = {}
for d in data:
    key = (d['delta_c7_dir'], d['delta_i2'])
    if key not in cross2:
        cross2[key] = []
    cross2[key].append(d['delta_H'])

for key in sorted(cross2.keys()):
    dhs = cross2[key]
    pred = 2*key[0] + 4*key[1]
    all_match = all(dh == pred for dh in dhs)
    print(f"  dc7_dir={key[0]:>2}, di2={key[1]:>2}: dH values={sorted(set(dhs))}, count={len(dhs)}, pred={pred}, match={all_match}")

# i2 absolute values
print(f"\n  i2_A range: {min(d['i2_A'] for d in data)}-{max(d['i2_A'] for d in data)}")
print(f"  i2_B range: {min(d['i2_B'] for d in data)}-{max(d['i2_B'] for d in data)}")

# Check higher-order terms: do we need i3?
print(f"\n  Are there examples with |delta_i2| > 1?")
large_di2 = [d for d in data if abs(d['delta_i2']) > 1]
if large_di2:
    print(f"  YES: {len(large_di2)} examples with |delta_i2| > 1")
    for d in large_di2[:5]:
        print(f"    dc7_dir={d['delta_c7_dir']}, di2={d['delta_i2']}, dH={d['delta_H']}")
else:
    print(f"  NO: delta_i2 always in {{-1, 0, 1}}")

# Compare with n=7 formula
print(f"\n{'='*70}")
print(f"COMPARISON WITH n=7")
print(f"{'='*70}")
print(f"""
At n=7: delta_H = 2*delta_c7_dir, delta_i2 = 0 always
  (all cycles use >= 3 vertices, disjoint pairs need 6 vertices,
   only 3-cycles can be disjoint, and disjoint 3-3 pairs are preserved)

At n=8: delta_H = 2*delta_c7_dir + 4*delta_i2
  - delta_c7_dir in {{-2,...,2}} (c7 uses 7 of 8 vertices)
  - delta_i2 in {{-1, 0, 1}} (from c3/c5 reshuffling)
  - Possible delta_H values: {{-6, -4, -2, 0, 2, 4, 6}}

The i2 channel opens at n=8 because:
1. Disjoint (3,5)-pairs need exactly 8 vertices = n
2. The Vitali atom reshuffles WHICH 3-cycles and 5-cycles exist
   (total counts preserved, but identities change)
3. This changes whether specific (3,5) pairs are disjoint
4. Each such change contributes +-4 to delta_H

At n=9, predict:
- delta_c9_dir channel opens (9-cycles on all 9 vertices)
- delta_i2 grows (more room for disjoint pairs)
- delta_i3 channel opens (disjoint triples, e.g., three 3-cycles)
- Formula: delta_H = 2*(dc7+dc9) + 4*di2 + 8*di3
""")

print("Done.")
