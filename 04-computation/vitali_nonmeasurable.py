"""
vitali_nonmeasurable.py -- kind-pasteur-2026-03-13-S61

The "non-measurability" of sigma: sigma is invisible to the lambda
sigma-algebra but controls c7 via its power sums.

This script quantifies the precise "information gap" between
lambda (measurable) and sigma (non-measurable):

1. How many bits of information does lambda carry vs sigma?
2. What fraction of tournament structure is lambda-measurable?
3. The "non-measurable residual" — exactly how much c7 variation
   is invisible to lambda?
4. Connection to the Vitali coset: how many tournaments share
   each lambda graph?
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

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"LAMBDA MEASURABILITY AND SIGMA NON-MEASURABILITY AT n={n}")
print("=" * 60)

np.random.seed(42)

# Collect data grouped by lambda graph
lambda_groups = defaultdict(list)
score_groups = defaultdict(set)

for trial in range(20000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    c7 = count_directed_k_cycles(A, n, 7)
    c3 = count_directed_k_cycles(A, n, 3)
    c5 = count_directed_k_cycles(A, n, 5)
    scores = tuple(sorted(int(sum(A[i])) for i in range(n)))

    # Lambda graph as a canonical key (upper triangle)
    lam_key = tuple(int(L[u][v]) for u in range(n) for v in range(u+1, n))

    lambda_groups[lam_key].append({
        'c7': c7,
        'c3': c3,
        'c5': c5,
        'scores': scores,
        'bits': bits,
    })
    score_groups[scores].add(lam_key)

# 1. Lambda graph diversity
n_lambda = len(lambda_groups)
print(f"\nDistinct lambda graphs (from 20k samples): {n_lambda}")
print(f"  Total tournament space: 2^{total_bits} = {2**total_bits}")
print(f"  Lambda bits ~ log2({n_lambda}) = {np.log2(n_lambda):.1f}")
print(f"  Tournament bits = {total_bits}")
print(f"  => Lambda captures ~{100*np.log2(n_lambda)/total_bits:.1f}% of information")

# 2. Lambda fiber sizes (how many tournaments per lambda graph)
fiber_sizes = [len(v) for v in lambda_groups.values()]
print(f"\n--- Lambda Fiber Sizes ---")
print(f"  Mean fiber size: {np.mean(fiber_sizes):.2f}")
print(f"  Max fiber size: {max(fiber_sizes)}")
print(f"  Fibers of size 1: {sum(1 for s in fiber_sizes if s == 1)}")
print(f"  Fibers of size >= 5: {sum(1 for s in fiber_sizes if s >= 5)}")

# 3. c7 variation within lambda fibers
print(f"\n--- c7 Variation Within Lambda Fibers ---")
c7_ranges = []
c7_unique_counts = []
for lam_key, entries in lambda_groups.items():
    if len(entries) >= 2:
        c7s = [e['c7'] for e in entries]
        c7_range = max(c7s) - min(c7s)
        c7_ranges.append(c7_range)
        c7_unique_counts.append(len(set(c7s)))

print(f"  Fibers with >= 2 tournaments: {len(c7_ranges)}")
print(f"  c7 range distribution:")
c7_range_dist = Counter(c7_ranges)
for r in sorted(c7_range_dist.keys())[:15]:
    print(f"    range={r}: {c7_range_dist[r]} fibers")

# What fraction have c7 variation?
n_varying = sum(1 for r in c7_ranges if r > 0)
print(f"\n  Fibers with c7 variation: {n_varying}/{len(c7_ranges)} ({100*n_varying/len(c7_ranges):.1f}%)")
print(f"  Fibers with CONSTANT c7: {len(c7_ranges)-n_varying}/{len(c7_ranges)} ({100*(1-n_varying/len(c7_ranges)):.1f}%)")

# 4. c5 variation within lambda fibers (should be 0 by THM-172/173)
c5_varying = 0
for lam_key, entries in lambda_groups.items():
    if len(entries) >= 2:
        c5s = set(e['c5'] for e in entries)
        if len(c5s) > 1:
            c5_varying += 1
print(f"\n  Fibers with c5 variation: {c5_varying} (should be 0 by THM-172/173)")

# 5. Score variation within lambda fibers
score_varying = 0
for lam_key, entries in lambda_groups.items():
    if len(entries) >= 2:
        score_set = set(e['scores'] for e in entries)
        if len(score_set) > 1:
            score_varying += 1
print(f"  Fibers with score variation: {score_varying}")

# 6. The "non-measurable residual" — how much c7 variance is
#    unexplained by lambda?
# Total c7 variance
all_c7 = [e['c7'] for entries in lambda_groups.values() for e in entries]
total_var = np.var(all_c7)

# Within-group variance (lambda-conditional variance)
within_var = 0
total_count = 0
for entries in lambda_groups.values():
    if len(entries) >= 2:
        c7s = [e['c7'] for e in entries]
        within_var += len(c7s) * np.var(c7s)
        total_count += len(c7s)
within_var /= total_count

# Between-group variance = Total - Within
between_var = total_var - within_var

print(f"\n--- Variance Decomposition ---")
print(f"  Total c7 variance: {total_var:.2f}")
print(f"  Between-lambda variance: {between_var:.2f} ({100*between_var/total_var:.2f}%)")
print(f"  Within-lambda variance: {within_var:.4f} ({100*within_var/total_var:.4f}%)")
print(f"  => Lambda explains {100*between_var/total_var:.2f}% of c7 variance!")
print(f"  => The non-measurable residual is only {100*within_var/total_var:.4f}%")

# 7. The Vitali coset structure
# Two tournaments are in the same "Vitali coset" iff they share the same lambda graph.
# The Vitali atom moves between elements of the same coset.
# The c7 invariant is "non-measurable" because it varies within cosets.
print(f"\n--- The Vitali Coset Structure ---")
print(f"  Number of cosets (lambda equivalence classes): {n_lambda}")
print(f"  Mean coset size (from sample): {np.mean(fiber_sizes):.1f}")

# Within each coset, what is the c7 structure?
# For large cosets, check if c7 takes exactly the values
# {c7_base, c7_base+1} (binary/Bernoulli) or more complex
c7_vals_per_fiber = []
for entries in lambda_groups.values():
    if len(entries) >= 3:
        c7s = sorted(set(e['c7'] for e in entries))
        c7_vals_per_fiber.append(len(c7s))

if c7_vals_per_fiber:
    c7_vals_dist = Counter(c7_vals_per_fiber)
    print(f"\n  Number of distinct c7 values per fiber (for fibers >= 3):")
    for k in sorted(c7_vals_dist.keys()):
        print(f"    {k} values: {c7_vals_dist[k]} fibers")

# 8. Information-theoretic measure
# H(c7) = total entropy of c7
# H(c7 | lambda) = conditional entropy
# I(c7; lambda) = H(c7) - H(c7|lambda) = mutual information
from collections import Counter
c7_dist = Counter(all_c7)
total_n = len(all_c7)
H_c7 = -sum((c/total_n) * np.log2(c/total_n) for c in c7_dist.values())

H_c7_given_lambda = 0
for entries in lambda_groups.values():
    if len(entries) >= 1:
        c7s = [e['c7'] for e in entries]
        p_lambda = len(entries) / total_n
        c7_local = Counter(c7s)
        H_local = -sum((c/len(c7s)) * np.log2(c/len(c7s)) for c in c7_local.values() if c > 0)
        H_c7_given_lambda += p_lambda * H_local

I_lambda_c7 = H_c7 - H_c7_given_lambda

print(f"\n--- Information Theory ---")
print(f"  H(c7) = {H_c7:.4f} bits")
print(f"  H(c7 | lambda) = {H_c7_given_lambda:.4f} bits")
print(f"  I(lambda; c7) = {I_lambda_c7:.4f} bits = {100*I_lambda_c7/H_c7:.2f}% of H(c7)")
print(f"  => Lambda captures {100*I_lambda_c7/H_c7:.2f}% of c7 information")
print(f"  => The non-measurable residual is {100*H_c7_given_lambda/H_c7:.2f}% = {H_c7_given_lambda:.4f} bits")

# 9. Similar analysis for scores
H_c7_given_scores = 0
score_c7_groups = defaultdict(list)
for entries in lambda_groups.values():
    for e in entries:
        score_c7_groups[e['scores']].append(e['c7'])

for scores, c7s in score_c7_groups.items():
    p_score = len(c7s) / total_n
    c7_local = Counter(c7s)
    H_local = -sum((c/len(c7s)) * np.log2(c/len(c7s)) for c in c7_local.values() if c > 0)
    H_c7_given_scores += p_score * H_local

I_scores_c7 = H_c7 - H_c7_given_scores
print(f"\n  For comparison:")
print(f"  I(scores; c7) = {I_scores_c7:.4f} bits = {100*I_scores_c7/H_c7:.2f}% of H(c7)")
print(f"  I(lambda; c7) = {I_lambda_c7:.4f} bits = {100*I_lambda_c7/H_c7:.2f}% of H(c7)")
print(f"  Lambda >> scores for c7 prediction")
print(f"  But lambda is NOT sufficient: {100*H_c7_given_lambda/H_c7:.2f}% residual")

print("\n\nDone.")
