"""
measurability_n8.py -- kind-pasteur-2026-03-13-S61

At n=7, lambda captures 99.99% of c7 information.
At n=8, the Vitali atom is more active (dc7 ranges [-3,+3]).
Does the non-measurable residual grow?

Also: does the sigma power sum determination extend?
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

n = 8
total_bits = n * (n-1) // 2

print("=" * 60)
print(f"MEASURABILITY ANALYSIS AT n={n}")
print("=" * 60)

np.random.seed(42)

# Collect data grouped by lambda graph
lambda_groups = defaultdict(list)

for trial in range(5000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    c7 = count_directed_k_cycles(A, n, 7)
    c3 = count_directed_k_cycles(A, n, 3)

    # Build sigma matrix
    A2 = A @ A
    sig_mat = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            sig_mat[u][v] = sig_mat[v][u] = n - 2 - int(A2[u][v]) - int(A2[v][u])

    tr_s2 = int(np.trace(sig_mat @ sig_mat))
    tr_s3 = int(np.trace(sig_mat @ sig_mat @ sig_mat))
    tr_s4 = int(np.trace(np.linalg.matrix_power(sig_mat, 4)))

    # Lambda graph key
    lam_key = tuple(int(L[u][v]) for u in range(n) for v in range(u+1, n))

    scores = tuple(sorted(int(sum(A[i])) for i in range(n)))

    lambda_groups[lam_key].append({
        'c7': c7,
        'c3': c3,
        'scores': scores,
        'tr_s2': tr_s2,
        'tr_s3': tr_s3,
        'tr_s4': tr_s4,
    })

n_lambda = len(lambda_groups)
print(f"\nDistinct lambda graphs (from 5k samples): {n_lambda}")
print(f"  Lambda bits ~ log2({n_lambda}) = {np.log2(n_lambda):.1f}")
print(f"  Tournament bits = {total_bits}")
print(f"  => Lambda captures ~{100*np.log2(n_lambda)/total_bits:.1f}% of bits")

# c7 variation within lambda fibers
print(f"\n--- c7 Within Lambda Fibers ---")
fibers_multi = [(k, v) for k, v in lambda_groups.items() if len(v) >= 2]
print(f"  Multi-tournament fibers: {len(fibers_multi)}")

c7_ranges = []
for _, entries in fibers_multi:
    c7s = [e['c7'] for e in entries]
    c7_ranges.append(max(c7s) - min(c7s))

c7_range_dist = Counter(c7_ranges)
for r in sorted(c7_range_dist.keys())[:20]:
    print(f"    range={r}: {c7_range_dist[r]}")

n_varying = sum(1 for r in c7_ranges if r > 0)
print(f"  c7 varying fibers: {n_varying}/{len(fibers_multi)} ({100*n_varying/len(fibers_multi) if fibers_multi else 0:.1f}%)")

# Variance decomposition
all_c7 = [e['c7'] for entries in lambda_groups.values() for e in entries]
total_var = np.var(all_c7)

within_var = 0
total_count = 0
for entries in lambda_groups.values():
    if len(entries) >= 2:
        c7s = [e['c7'] for e in entries]
        within_var += len(c7s) * np.var(c7s)
        total_count += len(c7s)
if total_count > 0:
    within_var /= total_count

print(f"\n--- Variance Decomposition ---")
print(f"  Total c7 variance: {total_var:.2f}")
print(f"  Within-lambda variance: {within_var:.4f}")
if total_var > 0:
    print(f"  Lambda explains {100*(1-within_var/total_var):.4f}% of c7 variance")

# Information theory
H_c7 = 0
c7_dist = Counter(all_c7)
total_n = len(all_c7)
for c in c7_dist.values():
    p = c / total_n
    if p > 0:
        H_c7 -= p * np.log2(p)

H_c7_given_lambda = 0
for entries in lambda_groups.values():
    p_lambda = len(entries) / total_n
    c7s = [e['c7'] for e in entries]
    c7_local = Counter(c7s)
    H_local = 0
    for c in c7_local.values():
        p = c / len(c7s)
        if p > 0:
            H_local -= p * np.log2(p)
    H_c7_given_lambda += p_lambda * H_local

I_lambda_c7 = H_c7 - H_c7_given_lambda
print(f"\n--- Information Theory ---")
print(f"  H(c7) = {H_c7:.4f} bits")
print(f"  H(c7 | lambda) = {H_c7_given_lambda:.4f} bits")
print(f"  I(lambda; c7) = {I_lambda_c7:.4f} bits = {100*I_lambda_c7/H_c7:.2f}% of H(c7)")

# Sigma power sum determination at n=8
print(f"\n--- Sigma Power Sums at n={n} ---")
tr_c7 = defaultdict(set)
for entries in lambda_groups.values():
    for e in entries:
        tr_c7[(e['tr_s2'], e['tr_s3'], e['tr_s4'])].add(e['c7'])

ambig = sum(1 for c7s in tr_c7.values() if len(c7s) > 1)
print(f"  (tr_s2, tr_s3, tr_s4) -> c7: {len(tr_c7)} groups, {ambig} ambiguous")
if ambig > 0:
    for (t2, t3, t4), c7s in sorted(tr_c7.items()):
        if len(c7s) > 1:
            print(f"    ({t2}, {t3}, {t4}): c7 in {sorted(c7s)}")
            if ambig > 10 and list(tr_c7.values()).index(c7s) > 10:
                print(f"    ... ({ambig} total ambiguous)")
                break

# Also check (tr_s2, tr_s3)
tr23_c7 = defaultdict(set)
for entries in lambda_groups.values():
    for e in entries:
        tr23_c7[(e['tr_s2'], e['tr_s3'])].add(e['c7'])
ambig23 = sum(1 for c7s in tr23_c7.values() if len(c7s) > 1)
print(f"  (tr_s2, tr_s3) -> c7: {len(tr23_c7)} groups, {ambig23} ambiguous")

# Score prediction
H_c7_given_scores = 0
score_c7 = defaultdict(list)
for entries in lambda_groups.values():
    for e in entries:
        score_c7[e['scores']].append(e['c7'])
for scores, c7s in score_c7.items():
    p_score = len(c7s) / total_n
    c7_local = Counter(c7s)
    H_local = 0
    for c in c7_local.values():
        p = c / len(c7s)
        if p > 0:
            H_local -= p * np.log2(p)
    H_c7_given_scores += p_score * H_local
I_scores_c7 = H_c7 - H_c7_given_scores
print(f"\n  I(scores; c7) = {I_scores_c7:.4f} bits = {100*I_scores_c7/H_c7:.2f}%")

print("\n\nDone.")
