#!/usr/bin/env python3
"""
beta2_omega_poincare.py - Poincare polynomial of Omega chain complex

For each tournament T, compute P_T(t) = sum dim(Omega_p) * t^p.

Key observations from prior analysis:
- Transitive: P(t) = ((1+t)^n - 1) / t (binomial coefficients)
- Regular n=5: P(t) = 5 + 10t + 10t^2 + 10t^3 + 5t^4 (palindromic!)
- chi = P(-1) = 1 - beta_1

Questions:
1. When is P(t) palindromic?
2. Is there a formula for P(t) in terms of tournament invariants?
3. Does palindromicity correlate with regularity or cycle structure?
4. What happens at n=6,7?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def omega_dims(A, n, max_p=None):
    """Compute dim(Omega_p) for all p."""
    if max_p is None:
        max_p = n - 1

    allowed = {}
    for p in range(-1, max_p + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    dims = []
    for p in range(max_p + 1):
        omega_p = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        dim = omega_p.shape[1] if omega_p.ndim == 2 else 0
        dims.append(dim)

    return dims


def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]):
                    c3 += 1
    return c3


def is_palindrome(lst):
    return lst == lst[::-1]


# ============================================================
# PART 1: Omega Poincare polynomials at n=4,5
# ============================================================
print("=" * 70)
print("OMEGA POINCARE POLYNOMIALS")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 1 << m
    t0 = time.time()

    profile_count = Counter()
    profile_scores = defaultdict(set)
    profile_c3 = defaultdict(set)
    profile_palin = defaultdict(int)

    for bits in range(total):
        A = build_adj(n, bits)
        scores = tuple(sorted([sum(row) for row in A]))
        c3 = count_3cycles(A, n)
        dims = omega_dims(A, n)
        profile = tuple(dims)
        profile_count[profile] += 1
        profile_scores[profile].add(scores)
        profile_c3[profile].add(c3)
        if is_palindrome(dims):
            profile_palin[profile] += 1

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments ({elapsed:.1f}s)")
    print(f"  Distinct Omega profiles: {len(profile_count)}")
    for prof in sorted(profile_count.keys()):
        cnt = profile_count[prof]
        scores = sorted(profile_scores[prof])
        c3s = sorted(profile_c3[prof])
        palin = is_palindrome(list(prof))
        chi = sum((-1)**p * prof[p] for p in range(len(prof)))
        print(f"  {list(prof)}: {cnt} tournaments, c3={c3s}, chi={chi}, palindromic={palin}")


# ============================================================
# PART 2: Palindromicity and regularity
# ============================================================
print(f"\n{'='*70}")
print("PALINDROMICITY ANALYSIS")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 1 << m
    palin_count = 0
    palin_scores = Counter()
    non_palin_scores = Counter()

    for bits in range(total):
        A = build_adj(n, bits)
        scores = tuple(sorted([sum(row) for row in A]))
        dims = omega_dims(A, n)
        if is_palindrome(dims):
            palin_count += 1
            palin_scores[scores] += 1
        else:
            non_palin_scores[scores] += 1

    print(f"\nn={n}: palindromic Omega profiles: {palin_count}/{total}")
    print(f"  Palindromic by score: {dict(sorted(palin_scores.items()))}")
    print(f"  Non-palindromic by score: {dict(sorted(non_palin_scores.items()))}")


# ============================================================
# PART 3: dim(A_p) vs dim(Omega_p) — what's the Omega rank deficiency?
# ============================================================
print(f"\n{'='*70}")
print("OMEGA RANK DEFICIENCY: dim(A_p) - dim(Omega_p)")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 1 << m

    deficiency_data = defaultdict(lambda: defaultdict(set))

    for bits in range(total):
        A = build_adj(n, bits)
        c3 = count_3cycles(A, n)

        allowed = {}
        for p in range(-1, n + 1):
            if p < 0:
                allowed[p] = []
            else:
                allowed[p] = enumerate_allowed_paths(A, n, p)

        for p in range(n):
            omega_p = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
            dim_omega = omega_p.shape[1] if omega_p.ndim == 2 else 0
            dim_A = len(allowed[p])
            deficiency = dim_A - dim_omega
            deficiency_data[p][(c3, deficiency)].add(bits)

    print(f"\nn={n}:")
    for p in sorted(deficiency_data.keys()):
        print(f"  p={p}:")
        for (c3, deficiency) in sorted(deficiency_data[p].keys()):
            count = len(deficiency_data[p][(c3, deficiency)])
            print(f"    c3={c3}, dim(A_{p})-dim(Omega_{p})={deficiency}: {count} tournaments")


# ============================================================
# PART 4: dim(Omega_2) formula
# ============================================================
print(f"\n{'='*70}")
print("dim(Omega_2) FORMULA")
print("=" * 70)

# dim(A_2) = C(n,3) + 2*c3 (each transitive triple gives 1 path, each 3-cycle gives 3)
# dim(Omega_2) = dim(A_2) - rank(Omega constraints)
# What determines the Omega constraint rank?

n = 5
m = n*(n-1)//2

# Collect: c3, dim(A2), dim(Omega2), constraint rank
formula_data = []
for bits in range(1 << m):
    A = build_adj(n, bits)
    c3 = count_3cycles(A, n)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    dim_A2 = len(paths2)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_Omega2 = omega2.shape[1] if omega2.ndim == 2 else 0
    constraint_rank = dim_A2 - dim_Omega2

    formula_data.append({
        'bits': bits, 'c3': c3,
        'dim_A2': dim_A2, 'dim_O2': dim_Omega2,
        'constraint_rank': constraint_rank
    })

# Check: dim(A_2) = C(5,3) + 2*c3?
print(f"\nn=5: dim(A_2) = C(5,3) + 2*c3?")
for row in formula_data[:10]:
    expected = 10 + 2*row['c3']
    match = "Y" if row['dim_A2'] == expected else "N"
    print(f"  T#{row['bits']}: c3={row['c3']}, dim(A2)={row['dim_A2']}, C(5,3)+2*c3={expected} [{match}]")

all_A2_match = all(row['dim_A2'] == 10 + 2*row['c3'] for row in formula_data)
print(f"  ALL match? {all_A2_match}")

# Now look at constraint_rank vs c3
print(f"\n  Constraint rank by c3:")
c3_cr = defaultdict(set)
for row in formula_data:
    c3_cr[row['c3']].add(row['constraint_rank'])
for c3 in sorted(c3_cr.keys()):
    print(f"    c3={c3}: constraint_rank in {sorted(c3_cr[c3])}")

# Check: constraint_rank = 2*c3?
print(f"\n  constraint_rank = 2*c3?")
match_2c3 = all(row['constraint_rank'] == 2*row['c3'] for row in formula_data)
print(f"    ALL match? {match_2c3}")
# If not, check 3*c3
match_3c3 = all(row['constraint_rank'] == 3*row['c3'] for row in formula_data)
print(f"    constraint_rank = 3*c3? {match_3c3}")

# What IS the constraint rank?
# For c3=0: constraint_rank=0 (no Omega constraints needed for 2-paths)
# For c3=1: constraint_rank=?
# The non-allowed 1-pairs (a,c) with c->a give constraints

# Count non-allowed pairs with mediators
print(f"\n  Non-allowed pairs with mediators:")
for c3_val in sorted(c3_cr.keys()):
    # Sample one
    for row in formula_data:
        if row['c3'] == c3_val:
            bits = row['bits']
            A = build_adj(n, bits)
            na_with_med = 0
            na_total = 0
            for a in range(n):
                for c in range(n):
                    if a != c and A[c][a] == 1:
                        na_total += 1
                        meds = [b for b in range(n) if b != a and b != c and A[a][b] == 1 and A[b][c] == 1]
                        if len(meds) > 0:
                            na_with_med += 1
            print(f"    c3={c3_val}: na_total={na_total}, na_with_mediators={na_with_med}, constraint_rank={row['constraint_rank']}")
            break


# ============================================================
# PART 5: Sample n=6,7 Omega profiles
# ============================================================
print(f"\n{'='*70}")
print("OMEGA PROFILES AT n=6 (SAMPLE)")
print("=" * 70)

import random
random.seed(42)

n = 6
total = 1 << (n*(n-1)//2)
sample_size = 2000

profile_count_6 = Counter()
for trial in range(sample_size):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    dims = omega_dims(A, n)
    profile = tuple(dims)
    profile_count_6[profile] += 1

print(f"\nn=6: {sample_size} sampled tournaments")
print(f"  Distinct profiles: {len(profile_count_6)}")
for prof in sorted(profile_count_6.keys(), key=lambda x: -profile_count_6[x])[:20]:
    cnt = profile_count_6[prof]
    palin = is_palindrome(list(prof))
    chi = sum((-1)**p * prof[p] for p in range(len(prof)))
    print(f"  {list(prof)}: {cnt} ({100*cnt/sample_size:.1f}%), chi={chi}, palin={palin}")


# n=7
print(f"\n{'='*70}")
print("OMEGA PROFILES AT n=7 (SAMPLE)")
print("=" * 70)

n = 7
total = 1 << (n*(n-1)//2)
sample_size = 500

profile_count_7 = Counter()
for trial in range(sample_size):
    bits = random.randint(0, total - 1)
    A = build_adj(n, bits)
    dims = omega_dims(A, n)
    profile = tuple(dims)
    profile_count_7[profile] += 1

print(f"\nn=7: {sample_size} sampled tournaments")
print(f"  Distinct profiles: {len(profile_count_7)}")
for prof in sorted(profile_count_7.keys(), key=lambda x: -profile_count_7[x])[:20]:
    cnt = profile_count_7[prof]
    palin = is_palindrome(list(prof))
    chi = sum((-1)**p * prof[p] for p in range(len(prof)))
    print(f"  {list(prof)}: {cnt} ({100*cnt/sample_size:.1f}%), chi={chi}, palin={palin}")


print("\n\nDone.")
