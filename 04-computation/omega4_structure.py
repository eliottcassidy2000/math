#!/usr/bin/env python3
"""
omega4_structure.py — opus-2026-03-13-S70

Investigate the structure of Omega_4 (regular 4-paths).
At n=5: Omega_4 is the regular HP count.
At n=6: Omega_4 uses 5 vertices from 6.

Key questions:
1. What local invariant determines Omega_4 on 5-vertex subsets?
2. Does it follow the pattern: Omega_m = C(n,m+1) - correction?
3. What's the correction term for m=4?

Known pattern:
  Omega_2 = C(n,3) - t_3  [Walsh degree 2]
  Omega_3 = C(n,4) - Σ_v[c3(N_in) + c3(N_out)]  [Walsh degree 4]
  Omega_4 = ???  [Walsh degree 6 expected]
"""

import numpy as np
from math import comb
from itertools import combinations, permutations
from collections import Counter, defaultdict

def adj_matrix(n, T_bits):
    A = np.zeros((n,n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (T_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_regular_m_paths(A, m):
    """Count regular m-paths in tournament A."""
    n = A.shape[0]
    count = 0
    def dfs(path, depth):
        nonlocal count
        if depth == m:
            count += 1
            return
        last = path[-1]
        for v in range(n):
            if v not in path and A[last][v] == 1:
                if depth >= 1 and A[path[-2]][v] != 1:
                    continue
                path.append(v)
                dfs(path, depth + 1)
                path.pop()
    for start in range(n):
        dfs([start], 0)
    return count

def omega4_on_subset(A, vertices):
    """Count regular 4-paths using exactly these 5 vertices."""
    count = 0
    for perm in permutations(vertices):
        v0, v1, v2, v3, v4 = perm
        # Path: v0→v1→v2→v3→v4
        # Regularity: v0→v2, v1→v3, v2→v4
        if (A[v0][v1] and A[v1][v2] and A[v2][v3] and A[v3][v4] and
            A[v0][v2] and A[v1][v3] and A[v2][v4]):
            count += 1
    return count

def count_3cycles(A_sub, vertices):
    """Count 3-cycles in subtournament on given vertices."""
    count = 0
    vlist = list(vertices)
    for i in range(len(vlist)):
        for j in range(i+1, len(vlist)):
            for k in range(j+1, len(vlist)):
                a, b, c = vlist[i], vlist[j], vlist[k]
                if (A_sub[a][b] and A_sub[b][c] and A_sub[c][a]) or \
                   (A_sub[a][c] and A_sub[c][b] and A_sub[b][a]):
                    count += 1
    return count

def score_sequence(A, vertices):
    """Get sorted score sequence of subtournament."""
    scores = []
    for v in vertices:
        s = sum(A[v][w] for w in vertices if w != v)
        scores.append(s)
    return tuple(sorted(scores))

# ============================================================
# n = 5: Omega_4 is the regular HP count
# ============================================================
print("="*70)
print("OMEGA_4 STRUCTURE: n = 5 (regular HPs)")
print("="*70)

n = 5
num_pairs = n*(n-1)//2

omega4_vals = {}
omega3_vals = {}
omega2_vals = {}
t3_vals = {}

for bits in range(2**num_pairs):
    A = adj_matrix(n, bits)
    omega4_vals[bits] = count_regular_m_paths(A, 4)
    omega3_vals[bits] = count_regular_m_paths(A, 3)
    omega2_vals[bits] = count_regular_m_paths(A, 2)
    t3_vals[bits] = count_3cycles(A, range(n))

print(f"  Omega_4 values: {sorted(set(omega4_vals.values()))}")
print(f"  Omega_4 distribution: {dict(sorted(Counter(omega4_vals.values()).items()))}")

# At n=5, Omega_4 uses all 5 vertices — it's the regular HP count
# The whole tournament determines it, no "local" structure on subsets
# But at n=6, it will decompose over 5-vertex subsets

# ============================================================
# n = 6: Omega_4 per 5-vertex subset analysis
# ============================================================
print(f"\n{'='*70}")
print("OMEGA_4 STRUCTURE: n = 6")
print("="*70)

n = 6
num_pairs = n*(n-1)//2
total = 2**num_pairs
print(f"  Total tournaments: {total}")

# Compute Omega_4 for all tournaments
print("  Computing Omega_4 for all tournaments...")
omega4_6 = {}
for bits in range(total):
    A = adj_matrix(n, bits)
    omega4_6[bits] = count_regular_m_paths(A, 4)
    if bits % 5000 == 0 and bits > 0:
        print(f"    {bits}/{total}...")

print(f"  Omega_4 values at n=6: {sorted(set(omega4_6.values()))}")
print(f"  Distribution: {dict(sorted(Counter(omega4_6.values()).items()))}")

# Check: does Omega_4 decompose as sum over 5-vertex subsets?
print("\n  Verifying Omega_4 = sum over 5-vertex subsets...")
mismatches = 0
for bits in range(min(500, total)):
    A = adj_matrix(n, bits)
    total_from_subsets = 0
    for combo in combinations(range(n), 5):
        total_from_subsets += omega4_on_subset(A, list(combo))
    if total_from_subsets != omega4_6[bits]:
        mismatches += 1
        if mismatches <= 3:
            print(f"    MISMATCH at T={bits}: Omega_4={omega4_6[bits]}, sum={total_from_subsets}")
if mismatches == 0:
    print(f"    CONFIRMED: Omega_4 = sum over C(6,5)=6 five-vertex subsets (500 checked)")
else:
    print(f"    {mismatches}/500 mismatches!")

# What does omega4_local depend on for 5-vertex subtournaments?
print("\n  omega4_local by 5-vertex tournament type (score sequence):")
per_sub_by_type = defaultdict(list)

for bits in range(total):
    A = adj_matrix(n, bits)
    for combo in combinations(range(n), 5):
        val = omega4_on_subset(A, list(combo))
        scores = score_sequence(A, combo)
        t3 = count_3cycles(A, combo)
        per_sub_by_type[(scores, t3)].append(val)

for (scores, t3), vals in sorted(per_sub_by_type.items()):
    unique_vals = sorted(set(vals))
    print(f"    scores={scores}, t3={t3}: omega4_local ∈ {unique_vals}, counts={Counter(vals)}")

# Key question: is omega4_local determined by score sequence + t3?
# Or does it need finer invariants?
print("\n  Is omega4_local a function of (score_seq, t3)?")
determined = True
for (scores, t3), vals in per_sub_by_type.items():
    if len(set(vals)) > 1:
        determined = False
        print(f"    NOT determined at scores={scores}, t3={t3}: vals={sorted(set(vals))}")

if determined:
    print("    YES — omega4_local is determined by (score_seq, t3)")

# Also check: is it determined by score sequence alone?
per_sub_by_score = defaultdict(list)
for (scores, t3), vals in per_sub_by_type.items():
    per_sub_by_score[scores].extend(vals)

print("\n  Is omega4_local a function of score_seq alone?")
for scores, vals in sorted(per_sub_by_score.items()):
    unique_vals = sorted(set(vals))
    if len(unique_vals) > 1:
        print(f"    NOT determined at scores={scores}: vals={unique_vals}")

# Correlation of Omega_4 with known invariants
print("\n  Correlations:")
o4_list = [omega4_6[b] for b in range(total)]
t3_list = []
for bits in range(total):
    A = adj_matrix(n, bits)
    t3_list.append(count_3cycles(A, range(n)))

print(f"    corr(Omega_4, t3) = {np.corrcoef(o4_list, t3_list)[0,1]:.6f}")

# Check the formula pattern: Omega_4 = C(n,5) - correction
C_n_5 = comb(n, 5)
print(f"\n  C({n},5) = {C_n_5}")
corrections = [C_n_5 - v for v in o4_list]
print(f"  Correction = C(n,5) - Omega_4: values = {sorted(set(corrections))}")

print("\nDONE.")
