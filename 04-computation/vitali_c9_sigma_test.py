"""
vitali_c9_sigma_test.py -- kind-pasteur-2026-03-13-S61

Use Vitali pairs (known lambda-preserving pairs) to test:
1. Does sigma preservation imply dc9=0?
   (We know sigma is NOT preserved, so dc9 can be nonzero)
2. Is dc9 determined by delta_sigma?
3. What is the relationship between dc7 and dc9?

At n=9: c9 = Hamiltonian directed cycles = tr(A^9)/9.
"""

import numpy as np
from itertools import combinations
from collections import defaultdict, Counter

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

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

n = 9
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"VITALI PAIR ANALYSIS AT n={n}")
print("=" * 60)

np.random.seed(42)
examples = []

for trial in range(5000):
    bits = int(np.random.randint(0, 2**31)) | (int(np.random.randint(0, 2**5)) << 31)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        # Found a Vitali pair
        c5_A = count_directed_k_cycles(A, n, 5)
        c5_B = count_directed_k_cycles(B, n, 5)
        c7_A = count_directed_k_cycles(A, n, 7)
        c7_B = count_directed_k_cycles(B, n, 7)
        c9_A = count_directed_k_cycles(A, n, 9)
        c9_B = count_directed_k_cycles(B, n, 9)

        dc5 = c5_B - c5_A
        dc7 = c7_B - c7_A
        dc9 = c9_B - c9_A

        # Sigma change
        A2_A = A @ A
        A2_B = B @ B
        dsig_sum = 0
        dsig_abs = 0
        dsig_pairs = 0
        for u in range(n):
            for v in range(u+1, n):
                sig_A = n - 2 - int(A2_A[u][v]) - int(A2_A[v][u])
                sig_B = n - 2 - int(A2_B[u][v]) - int(A2_B[v][u])
                ds = sig_B - sig_A
                dsig_sum += ds
                dsig_abs += abs(ds)
                if ds != 0:
                    dsig_pairs += 1

        examples.append({
            'bits': bits,
            'subset': subset,
            'dc5': dc5,
            'dc7': dc7,
            'dc9': dc9,
            'dsig_sum': dsig_sum,
            'dsig_abs': dsig_abs,
            'dsig_pairs': dsig_pairs,
        })
        break

    if trial % 1000 == 0 and trial > 0:
        print(f"  {trial}: {len(examples)} Vitali pairs found")

print(f"\nTotal Vitali pairs: {len(examples)}")

# Filter to non-trivial (H-changing) pairs
nontrivial = [ex for ex in examples if ex['dc7'] != 0 or ex['dc9'] != 0]
print(f"Non-trivial (dc7 or dc9 nonzero): {len(nontrivial)}")

# dc5 check
dc5_dist = Counter(ex['dc5'] for ex in examples)
print(f"\ndc5 distribution: {dict(sorted(dc5_dist.items()))}")
print(f"dc5=0 always? {all(ex['dc5'] == 0 for ex in examples)}")

# dc7 and dc9 distributions
dc7_dist = Counter(ex['dc7'] for ex in examples)
dc9_dist = Counter(ex['dc9'] for ex in examples)
print(f"dc7 distribution: {dict(sorted(dc7_dist.items()))}")
print(f"dc9 distribution: {dict(sorted(dc9_dist.items()))}")

# Sigma change statistics
print(f"\nSigma change statistics:")
print(f"  dsig_sum always 0? {all(ex['dsig_sum'] == 0 for ex in examples)}")
dsig_abs_dist = Counter(ex['dsig_abs'] for ex in examples)
dsig_pairs_dist = Counter(ex['dsig_pairs'] for ex in examples)
print(f"  dsig_abs distribution: {dict(sorted(dsig_abs_dist.items()))}")
print(f"  dsig_pairs distribution: {dict(sorted(dsig_pairs_dist.items()))}")

# Relationship between dc7, dc9, and sigma changes
print(f"\n{'='*60}")
print("RELATIONSHIP: dc7, dc9, delta_sigma")
print(f"{'='*60}")

# Is dc7+dc9 determined by sigma change?
dc7_dc9_vs_dsig = defaultdict(list)
for ex in nontrivial:
    key = (ex['dsig_abs'], ex['dsig_pairs'])
    dc7_dc9_vs_dsig[key].append((ex['dc7'], ex['dc9']))

for key, vals in sorted(dc7_dc9_vs_dsig.items()):
    print(f"  dsig_abs={key[0]}, dsig_pairs={key[1]}: (dc7,dc9) values = {Counter(vals)}")

# Is there a linear relation dc9 = a*dc7 + b?
dc7_vals = np.array([ex['dc7'] for ex in nontrivial], dtype=float)
dc9_vals = np.array([ex['dc9'] for ex in nontrivial], dtype=float)

if len(nontrivial) > 2:
    from numpy.linalg import lstsq
    X = np.column_stack([np.ones(len(dc7_vals)), dc7_vals])
    coeffs, _, _, _ = lstsq(X, dc9_vals, rcond=None)
    err = np.max(np.abs(X @ coeffs - dc9_vals))
    print(f"\n  dc9 = {coeffs[0]:.4f} + {coeffs[1]:.4f} * dc7")
    print(f"  Max error: {err:.4f}")
    print(f"  Perfect linear relation? {err < 0.001}")

    # Joint distribution
    print(f"\n  Joint (dc7, dc9) distribution:")
    joint = Counter((ex['dc7'], ex['dc9']) for ex in nontrivial)
    for k in sorted(joint.keys()):
        print(f"    dc7={k[0]:+d}, dc9={k[1]:+d}: {joint[k]}")

# Scatter plot summary
print(f"\n  Full joint (dc7, dc9) over ALL examples:")
joint_all = Counter((ex['dc7'], ex['dc9']) for ex in examples)
for k in sorted(joint_all.keys()):
    print(f"    dc7={k[0]:+d}, dc9={k[1]:+d}: {joint_all[k]}")

print("\nDone.")
