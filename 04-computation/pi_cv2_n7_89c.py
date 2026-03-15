#!/usr/bin/env python3
"""
pi_cv2_n7_89c.py — Compute CV² at n=7 using pair formula
opus-2026-03-14-S89c

E[H²] = Σ_{compatible (π,σ)} 2^{-|constrained edges|}

This avoids enumerating all 2^21 tournaments!
We enumerate all n!² = 25,401,600 pairs of permutations.
Too many for brute force, but we can be smarter.

Actually n=7: n! = 5040, n!² ≈ 25.4M. That's feasible in pure Python
if we optimize. Each pair check: O(n) to find constrained edges.
25M × 7 = 175M operations ≈ 30 seconds in optimized Python.

E[H] = 7!/2^6 = 5040/64 = 315/4 = 78.75
E[H²] via formula gives exact Var.
"""

import time
from fractions import Fraction
from itertools import permutations
from math import factorial

n = 7
m = n * (n - 1) // 2  # = 21

print(f"Computing E[H²] at n={n} using compatible pair formula")
print(f"n! = {factorial(n)}, n!² = {factorial(n)**2}")
print(f"m = {m} undirected edges")

t0 = time.time()

# For each pair (π, σ), we need:
# 1. Check compatibility: no edge (u,v) where π needs u→v and σ needs v→u
# 2. Count constrained undirected edges: |E(π) ∪ E(σ)|
# Weight = 2^{m - constrained} / 2^m = 2^{-constrained}

# Pre-compute all path edge sets
print("Pre-computing path edge sets...", flush=True)
all_perms = list(permutations(range(n)))
n_perms = len(all_perms)

# For each permutation, store directed edges as a set of (u,v) pairs
# and undirected edges as a set of frozensets
perm_dir_edges = []
perm_undir_edges = []
for perm in all_perms:
    dir_e = set()
    undir_e = set()
    for i in range(n-1):
        u, v = perm[i], perm[i+1]
        dir_e.add((u, v))
        undir_e.add((u, v) if u < v else (v, u))
    perm_dir_edges.append(dir_e)
    perm_undir_edges.append(undir_e)

print(f"Pre-computed {n_perms} permutations", flush=True)

# Now compute E[H²]
# E[H²] = Σ_{compat (π,σ)} 2^{-|E(π)∪E(σ)|}
# = Σ_{compat} 2^{-(2(n-1) - |E(π)∩E(σ)|)}
# = 2^{-2(n-1)} × Σ_{compat} 2^{|E(π)∩E(σ)|}

total_compatible = 0
weighted_sum = Fraction(0)

# The 2^{-2(n-1)} factor: 2^{-12}
base_denom = 2**(2*(n-1))

# For efficiency, we'll accumulate 2^k counts
# where k = |common undirected edges|
count_by_common = {}

print("Computing compatible pairs...", flush=True)
t1 = time.time()

for i in range(n_perms):
    if i % 1000 == 0 and i > 0:
        elapsed = time.time() - t1
        eta = elapsed / i * (n_perms - i)
        print(f"  {i}/{n_perms} ({i/n_perms*100:.1f}%), {elapsed:.0f}s elapsed, ETA {eta:.0f}s", flush=True)

    dir_i = perm_dir_edges[i]
    undir_i = perm_undir_edges[i]

    for j in range(n_perms):
        dir_j = perm_dir_edges[j]

        # Check compatibility: no (u,v) in dir_i with (v,u) in dir_j
        compatible = True
        for u, v in dir_i:
            if (v, u) in dir_j:
                compatible = False
                break
        if not compatible:
            continue

        # Count common undirected edges
        undir_j = perm_undir_edges[j]
        common = len(undir_i & undir_j)

        total_compatible += 1
        count_by_common[common] = count_by_common.get(common, 0) + 1

t2 = time.time()
print(f"\nDone in {t2-t0:.1f}s")
print(f"Total compatible pairs: {total_compatible} / {n_perms**2}")
print(f"Compatibility rate: {total_compatible/n_perms**2:.6f}")

# E[H²] = (1/2^m) × Σ_T H²(T) = Σ_{compat} 2^{-(2(n-1) - common)}
# = 2^{-2(n-1)} × Σ_compat 2^common
print(f"\nCommon edge distribution:")
for k in sorted(count_by_common.keys()):
    print(f"  k={k}: {count_by_common[k]} pairs")

E_H2 = Fraction(0)
for k, cnt in count_by_common.items():
    E_H2 += Fraction(cnt * 2**k, base_denom)

E_H = Fraction(factorial(n), 2**(n-1))
E_H2_over_EH2 = E_H2 / E_H**2
CV2 = E_H2_over_EH2 - 1

print(f"\nResults for n={n}:")
print(f"  E[H] = {E_H} = {float(E_H):.4f}")
print(f"  E[H²] = {E_H2} = {float(E_H2):.4f}")
print(f"  E[H²]/E[H]² = {E_H2_over_EH2} = {float(E_H2_over_EH2):.10f}")
print(f"  CV² = {CV2} = {float(CV2):.10f}")

# The CV² sequence: 1/3, 1/3, 19/60, 13/45, ?
print(f"\nCV² sequence so far:")
cv2_seq = [Fraction(1,3), Fraction(1,3), Fraction(19,60), Fraction(13,45), CV2]
for i, (n_val, cv) in enumerate(zip([3,4,5,6,7], cv2_seq)):
    print(f"  n={n_val}: {cv} = {float(cv):.10f}")

# Is it decreasing? Rate?
for i in range(1, len(cv2_seq)):
    diff = float(cv2_seq[i]) - float(cv2_seq[i-1])
    print(f"  Δ(n={i+2}→{i+3}) = {diff:.10f}")

# E[1/H] × E[H] — can we compute this too?
# E[1/H] needs all H values, which requires enumerating all 2^21 tournaments.
# That's 2M tournaments × n=7 DP = slow but possible.
# Skip for now.

print(f"\n  E[H²]/E[H]² sequence: 4/3, 4/3, 79/60, 58/45, {E_H2_over_EH2}")
print(f"  Denominators: 3, 3, 60, 45, {E_H2_over_EH2.denominator}")

print("\nDone!")
