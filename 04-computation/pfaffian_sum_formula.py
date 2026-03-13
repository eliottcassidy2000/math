#!/usr/bin/env python3
"""
pfaffian_sum_formula.py -- kind-pasteur-2026-03-13-S61

PROVED: det(I + 2A(T)) = (sum_i (-1)^i Pf(S_ii))^2

where S = A - A^T is the skew-adjacency matrix.

At n=5: |Pf_sum| takes values {1,3,5,7,9} and maps bijectively to c5:
  |Pf_sum|=1: c5=0
  |Pf_sum|=3: c5=0 (but different score class)
  |Pf_sum|=5: c5=1
  ...wait, this isn't quite right. Let me check the exact relationship.

At n=7 (the ambiguous pair):
  H=109: |Pf_sum|=13, c7=8  -> 2*c7-3 = 13? 2*8-3=13 YES!
  H=111: |Pf_sum|=19, c7=9  -> 2*c7+1 = 19? 2*9+1=19 YES!

These are DIFFERENT formulas for different lambda classes.
So Pf_sum is NOT a simple function of c7 alone.

This script investigates the exact relationship between Pf_sum, H, c_n, and other
tournament invariants.

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k <= 2:
        return 0
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def pfaffian(M):
    n = len(M)
    if n == 0:
        return 1
    if n == 1:
        return 0
    if n == 2:
        return M[0][1]
    result = 0
    for j in range(1, n):
        if M[0][j] == 0:
            continue
        indices = [i for i in range(n) if i != 0 and i != j]
        sub = [[M[indices[a]][indices[b]] for b in range(len(indices))] for a in range(len(indices))]
        result += ((-1) ** (j + 1)) * M[0][j] * pfaffian(sub)
    return result


def pfaffian_sum(A, n):
    """Compute the Pfaffian sum = sum_i (-1)^i * Pf(S_ii)."""
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
    total = 0
    for i in range(n):
        remaining = [j for j in range(n) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
        pf = pfaffian(S_del)
        total += ((-1) ** i) * pf
    return total


# ========================================================================
# EXHAUSTIVE n=5: Pf_sum vs all invariants
# ========================================================================
print("=" * 70)
print("n=5 EXHAUSTIVE: Pfaffian sum vs tournament invariants")
print("=" * 70)

n5 = 5
data_5 = []

for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    H = count_ham_paths(A, n5)
    scores = tuple(sorted(sum(A[v]) for v in range(n5)))
    c3 = sum(count_directed_ham_cycles_on_subset(A, list(sub))
             for sub in combinations(range(n5), 3))
    c5 = count_directed_ham_cycles_on_subset(A, list(range(n5)))
    ps = pfaffian_sum(A, n5)
    data_5.append({'bits': bits, 'H': H, 'scores': scores, 'c3': c3, 'c5': c5, 'ps': ps})

# Group by (scores, c3, c5) and check Pf_sum
by_class = defaultdict(list)
for d in data_5:
    by_class[(d['scores'], d['c3'], d['c5'])].append(d)

print(f"\n  Classes by (scores, c3, c5):")
for key in sorted(by_class.keys()):
    group = by_class[key]
    ps_vals = sorted(set(d['ps'] for d in group))
    H_vals = sorted(set(d['H'] for d in group))
    print(f"    scores={key[0]}, c3={key[1]}, c5={key[2]}: H={H_vals}, |Pf_sum|={sorted(set(abs(d['ps']) for d in group))}, Pf_sum={ps_vals}")

# KEY: Is |Pf_sum| = 2*c_n + 1?
print(f"\n  Testing |Pf_sum| = 2*c5 + 1:")
for d in data_5:
    if abs(d['ps']) != 2 * d['c5'] + 1:
        print(f"    FAIL: bits={d['bits']}, H={d['H']}, c5={d['c5']}, Pf_sum={d['ps']}, |Pf_sum|={abs(d['ps'])}, 2*c5+1={2*d['c5']+1}")
        break
else:
    print(f"    VERIFIED for all {len(data_5)} tournaments!")


# ========================================================================
# EXHAUSTIVE n=3: Pf_sum vs c3
# ========================================================================
print(f"\n{'='*70}")
print("n=3 EXHAUSTIVE")
print("=" * 70)

n3 = 3
for bits in range(1 << 3):
    A = binary_to_tournament(bits, n3)
    H = count_ham_paths(A, n3)
    c3 = count_directed_ham_cycles_on_subset(A, list(range(n3)))
    ps = pfaffian_sum(A, n3)
    print(f"  bits={bits}: H={H}, c3={c3}, Pf_sum={ps}, |Pf_sum|={abs(ps)}, 2*c3+1={2*c3+1}")


# ========================================================================
# EXHAUSTIVE n=4: Pf(S) (even n, full Pfaffian)
# ========================================================================
print(f"\n{'='*70}")
print("n=4 EXHAUSTIVE: Pf(S) vs tournament invariants")
print("=" * 70)

n4 = 4
for bits in range(1 << 6):
    A = binary_to_tournament(bits, n4)
    S = [[A[i][j] - A[j][i] for j in range(n4)] for i in range(n4)]
    pf = pfaffian(S)
    H = count_ham_paths(A, n4)
    c3 = sum(count_directed_ham_cycles_on_subset(A, list(sub))
             for sub in combinations(range(n4), 3))
    scores = tuple(sorted(sum(A[v]) for v in range(n4)))
    if abs(pf) > 0 or H == 1:  # show interesting ones
        print(f"  bits={bits}: H={H}, c3={c3}, scores={scores}, Pf(S)={pf}")


# ========================================================================
# n=7 SAMPLE: Pf_sum vs c7
# ========================================================================
print(f"\n{'='*70}")
print("n=7 SAMPLE: Pfaffian sum vs c7")
print("=" * 70)

n7 = 7
import random
random.seed(42)

data_7 = []
sample_bits = random.sample(range(1 << 21), 200)

for bits in sample_bits:
    A = binary_to_tournament(bits, n7)
    H = count_ham_paths(A, n7)
    c7 = count_directed_ham_cycles_on_subset(A, list(range(n7)))
    ps = pfaffian_sum(A, n7)
    scores = tuple(sorted(sum(A[v]) for v in range(n7)))
    data_7.append({'bits': bits, 'H': H, 'c7': c7, 'ps': ps, 'scores': scores})

# Check |Pf_sum| = 2*c7 + 1?
for d in data_7:
    if abs(d['ps']) != 2 * d['c7'] + 1:
        print(f"  COUNTEREXAMPLE: bits={d['bits']}, H={d['H']}, c7={d['c7']}, |Pf_sum|={abs(d['ps'])}, 2*c7+1={2*d['c7']+1}")

# Show distribution
by_c7 = defaultdict(list)
for d in data_7:
    by_c7[d['c7']].append(d)

print(f"\n  c7 -> |Pf_sum| distribution:")
for c7_val in sorted(by_c7.keys()):
    ps_vals = sorted(set(abs(d['ps']) for d in by_c7[c7_val]))
    count = len(by_c7[c7_val])
    print(f"    c7={c7_val}: |Pf_sum|={ps_vals}, count={count}")


# Check: does |Pf_sum| depend on more than just c7?
# For each c7, is |Pf_sum| unique?
print(f"\n  c7 -> unique |Pf_sum|?")
ambig_count = 0
for c7_val in sorted(by_c7.keys()):
    ps_abs_vals = set(abs(d['ps']) for d in by_c7[c7_val])
    if len(ps_abs_vals) > 1:
        ambig_count += 1
        print(f"    c7={c7_val}: MULTIPLE |Pf_sum| = {sorted(ps_abs_vals)}")
print(f"  Ambiguous: {ambig_count} c7 values")


# ========================================================================
# THE DEEP FORMULA
# ========================================================================
print(f"\n{'='*70}")
print("FORMULA SEARCH: |Pf_sum| as function of (H, c3, c5, c7, ...)")
print("=" * 70)

# At n=5 we found |Pf_sum| = 2*c5 + 1 exactly.
# Is there a similar formula at n=7?

# Check: |Pf_sum| = a*c7 + b*c5 + c*c3 + d?
# We need c3, c5, c7 for the sample

for d in data_7[:10]:
    bits = d['bits']
    A = binary_to_tournament(bits, n7)
    c3 = sum(count_directed_ham_cycles_on_subset(A, list(sub))
             for sub in combinations(range(n7), 3))
    c5 = sum(count_directed_ham_cycles_on_subset(A, list(sub))
             for sub in combinations(range(n7), 5))
    c7 = d['c7']
    ps = d['ps']
    H = d['H']

    # From OCF: H = 1 + 2*(c3+c5+c7) + 4*alpha_2
    alpha_1 = c3 + c5 + c7
    alpha_2 = (H - 1 - 2*alpha_1) // 4

    print(f"  bits={bits}: H={H}, c3={c3}, c5={c5}, c7={c7}, alpha2={alpha_2}, Pf_sum={ps}")
    print(f"    2*(c3+c5+c7)+1 = {2*(c3+c5+c7)+1}")
    print(f"    H - 4*alpha2 = {H - 4*alpha_2}")
    print(f"    |Pf_sum| = {abs(ps)}")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
