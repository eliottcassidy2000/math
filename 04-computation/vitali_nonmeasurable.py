"""
HISTORICAL SAMPLE; SCOPE-CORRECTED 2026-10-05.

The 20000-sample experiment concerns a finite coarse sigma-algebra generated
by labelled lambda fibres. It does not construct a classical Vitali set or
establish non-Lebesgue-measurability. Every observable is measurable on the
full finite tournament space.

Correction lineage: trace(A**7)/7 counts closed walks modulo cyclic starting
position here, NOT simple directed 7-cycles. All such quantities and output
labels below are now called walk7. For a strong four-vertex tournament it is2
although the simple7-cycle count is0. The within-fibre variance is normalized
by the full sample size, including singleton fibres with zero contribution.
See 05-knowledge/results/vitali_selector_certificate_20261005.md.

This is the historical seed42 sample, not an exhaustive fibre census.
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

def closed_walks_divided_by_length(A, n, k):
    """Exact integer quotient for the prime lengths3,5,7 used here.

    No loop or2-cycle exists in a tournament, so at lengths3 and5 every
    closed walk is simple. Length7 also admits a3-cycle/4-cycle concatenation.
    General lengths need not give an integral quotient and are not accepted.
    """
    if k not in (3, 5, 7):
        raise ValueError('this historical quotient uses only prime lengths3,5,7')
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"FINITE LAMBDA-FIBRE INFORMATION LOSS AT n={n}")
print("=" * 60)

np.random.seed(42)

# Collect data grouped by lambda graph
lambda_groups = defaultdict(list)
score_groups = defaultdict(set)

for trial in range(20000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    walk7 = closed_walks_divided_by_length(A, n, 7)
    c3 = closed_walks_divided_by_length(A, n, 3)
    c5 = closed_walks_divided_by_length(A, n, 5)
    scores = tuple(sorted(int(sum(A[i])) for i in range(n)))

    # Lambda graph as a canonical key (upper triangle)
    lam_key = tuple(int(L[u][v]) for u in range(n) for v in range(u+1, n))

    lambda_groups[lam_key].append({
        'walk7': walk7,
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
print(f"  Observed fibre-count log2({n_lambda}) = {np.log2(n_lambda):.1f}")
print(f"  Tournament bits = {total_bits}")
print(f"  => Observed log-count ratio {100*np.log2(n_lambda)/total_bits:.1f}% (not an entropy estimate)")

# 2. Lambda fiber sizes (how many tournaments per lambda graph)
fiber_sizes = [len(v) for v in lambda_groups.values()]
print(f"\n--- Lambda Fiber Sizes ---")
print(f"  Mean fiber size: {np.mean(fiber_sizes):.2f}")
print(f"  Max fiber size: {max(fiber_sizes)}")
print(f"  Fibers of size 1: {sum(1 for s in fiber_sizes if s == 1)}")
print(f"  Fibers of size >= 5: {sum(1 for s in fiber_sizes if s >= 5)}")

# 3. walk7 variation within lambda fibers
print(f"\n--- walk7 Variation Within Lambda Fibers ---")
walk7_ranges = []
walk7_unique_counts = []
for lam_key, entries in lambda_groups.items():
    if len(entries) >= 2:
        walk7s = [e['walk7'] for e in entries]
        walk7_range = max(walk7s) - min(walk7s)
        walk7_ranges.append(walk7_range)
        walk7_unique_counts.append(len(set(walk7s)))

print(f"  Fibers with >= 2 tournaments: {len(walk7_ranges)}")
print(f"  walk7 range distribution:")
walk7_range_dist = Counter(walk7_ranges)
for r in sorted(walk7_range_dist.keys())[:15]:
    print(f"    range={r}: {walk7_range_dist[r]} fibers")

# What fraction have walk7 variation?
n_varying = sum(1 for r in walk7_ranges if r > 0)
print(f"\n  Fibers with walk7 variation: {n_varying}/{len(walk7_ranges)} ({100*n_varying/len(walk7_ranges):.1f}%)")
print(f"  Fibers with CONSTANT walk7: {len(walk7_ranges)-n_varying}/{len(walk7_ranges)} ({100*(1-n_varying/len(walk7_ranges)):.1f}%)")

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

# 6. The "lambda-unresolved residual" — how much walk7 variance is
#    unexplained by lambda?
# Total walk7 variance
all_walk7 = [e['walk7'] for entries in lambda_groups.values() for e in entries]
total_var = np.var(all_walk7)

# Within-group variance (lambda-conditional variance)
within_var = 0
for entries in lambda_groups.values():
    if len(entries) >= 2:
        walk7s = [e['walk7'] for e in entries]
        within_var += len(walk7s) * np.var(walk7s)
within_var /= len(all_walk7)  # Singleton fibres contribute zero, but remain in the denominator.

# Between-group variance = Total - Within
between_var = total_var - within_var

print(f"\n--- Variance Decomposition ---")
print(f"  Total walk7 variance: {total_var:.2f}")
print(f"  Between-lambda variance: {between_var:.2f} ({100*between_var/total_var:.2f}%)")
print(f"  Within-lambda variance: {within_var:.4f} ({100*within_var/total_var:.4f}%)")
print(f"  => Lambda explains {100*between_var/total_var:.2f}% of walk7 variance!")
print(f"  => The lambda-unresolved residual is only {100*within_var/total_var:.4f}%")

# 7. The Vitali coset structure
# Two tournaments are in the same "Vitali coset" iff they share the same lambda graph.
# The Vitali atom moves between elements of the same coset.
# The walk7 invariant is "not lambda-determined" because it varies within cosets.
print(f"\n--- The Labelled Lambda Fibre Structure ---")
print(f"  Number of sampled lambda fibres: {n_lambda}")
print(f"  Mean fibre size (from sample): {np.mean(fiber_sizes):.1f}")

# Within each coset, what is the walk7 structure?
# For large cosets, check if walk7 takes exactly the values
# {walk7_base, walk7_base+1} (binary/Bernoulli) or more complex
walk7_vals_per_fiber = []
for entries in lambda_groups.values():
    if len(entries) >= 3:
        walk7s = sorted(set(e['walk7'] for e in entries))
        walk7_vals_per_fiber.append(len(walk7s))

if walk7_vals_per_fiber:
    walk7_vals_dist = Counter(walk7_vals_per_fiber)
    print(f"\n  Number of distinct walk7 values per fiber (for fibers >= 3):")
    for k in sorted(walk7_vals_dist.keys()):
        print(f"    {k} values: {walk7_vals_dist[k]} fibers")

# 8. Information-theoretic measure
# H(walk7) = total entropy of walk7
# H(walk7 | lambda) = conditional entropy
# I(walk7; lambda) = H(walk7) - H(walk7|lambda) = mutual information
from collections import Counter
walk7_dist = Counter(all_walk7)
total_n = len(all_walk7)
H_walk7 = -sum((c/total_n) * np.log2(c/total_n) for c in walk7_dist.values())

H_walk7_given_lambda = 0
for entries in lambda_groups.values():
    if len(entries) >= 1:
        walk7s = [e['walk7'] for e in entries]
        p_lambda = len(entries) / total_n
        walk7_local = Counter(walk7s)
        H_local = -sum((c/len(walk7s)) * np.log2(c/len(walk7s)) for c in walk7_local.values() if c > 0)
        H_walk7_given_lambda += p_lambda * H_local

I_lambda_walk7 = H_walk7 - H_walk7_given_lambda

print(f"\n--- Information Theory ---")
print(f"  H(walk7) = {H_walk7:.4f} bits")
print(f"  H(walk7 | lambda) = {H_walk7_given_lambda:.4f} bits")
print(f"  I(lambda; walk7) = {I_lambda_walk7:.4f} bits = {100*I_lambda_walk7/H_walk7:.2f}% of H(walk7)")
print(f"  => Lambda captures {100*I_lambda_walk7/H_walk7:.2f}% of walk7 information")
print(f"  => The lambda-unresolved residual is {100*H_walk7_given_lambda/H_walk7:.2f}% = {H_walk7_given_lambda:.4f} bits")

# 9. Similar analysis for scores
H_walk7_given_scores = 0
score_walk7_groups = defaultdict(list)
for entries in lambda_groups.values():
    for e in entries:
        score_walk7_groups[e['scores']].append(e['walk7'])

for scores, walk7s in score_walk7_groups.items():
    p_score = len(walk7s) / total_n
    walk7_local = Counter(walk7s)
    H_local = -sum((c/len(walk7s)) * np.log2(c/len(walk7s)) for c in walk7_local.values() if c > 0)
    H_walk7_given_scores += p_score * H_local

I_scores_walk7 = H_walk7 - H_walk7_given_scores
print(f"\n  For comparison:")
print(f"  I(scores; walk7) = {I_scores_walk7:.4f} bits = {100*I_scores_walk7/H_walk7:.2f}% of H(walk7)")
print(f"  I(lambda; walk7) = {I_lambda_walk7:.4f} bits = {100*I_lambda_walk7/H_walk7:.2f}% of H(walk7)")
print(f"  Lambda >> scores for walk7 prediction")
print(f"  But lambda is NOT sufficient: {100*H_walk7_given_lambda/H_walk7:.2f}% residual")

print("\n\nDone.")
