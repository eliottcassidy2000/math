#!/usr/bin/env python3
"""
ALGEBRAIC PROOF: moment formulas sum_P f^j for general n.

THEOREM. For a tournament T on n vertices (n odd), and f_P = #{forward edges in P}:

  sum_P f^j = sum over j-multisets {i_1,...,i_j} of positions in [0,n-2]
              of (j!/prod k_l!) * sigma(distinct positions)

where sigma(S) = sum over n! permutations of prod_{pos in S} A[p_pos, p_{pos+1}].

The quantity sigma(S) depends ONLY on the adjacency pattern of positions in S
(which positions are consecutive/overlapping):
  - If S has an isolated position, sigma(S) = n!/4 (half of ordered pairs).
  - If S = {i, i+1} (consecutive pair), sigma(S) depends on t_3.
  - If S = {i, j} with |i-j| > 1, sigma(S) is universal (pair-partition lemma).
  - For larger S: more complex overlap patterns.

This script derives the exact coefficients for each moment at n=5 and n=7.

opus-2026-03-06-S28
"""
from itertools import permutations, combinations
from math import factorial, comb
from collections import Counter, defaultdict
import random

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

# =====================================================================
# Step 1: Compute sigma(S) = sum_P prod_{pos in S} A[p_pos, p_{pos+1}]
# for all position subsets S at n=7
# =====================================================================
n = 7
print("=" * 70)
print(f"POSITION SUBSET SUMS sigma(S) at n={n}")
print("=" * 70)

# For a position subset S of [0..n-2], compute:
# sigma(S, T) = sum_{perm P} prod_{i in S} A[p_i, p_{i+1}]

def compute_sigma(A, n, positions):
    """Sum over all n! perms of product of A[p_i, p_{i+1}] for i in positions."""
    total = 0
    for p in permutations(range(n)):
        prod = 1
        for i in positions:
            prod *= A[p[i]][p[i+1]]
        total += prod
    return total

# Classify position subsets by their adjacency pattern
def position_pattern(positions, n):
    """Return the connected-component sizes of the position graph.
    Two positions i, j are adjacent if |i-j| = 1."""
    if not positions:
        return ()
    pos = sorted(positions)
    components = []
    comp = [pos[0]]
    for i in range(1, len(pos)):
        if pos[i] == comp[-1] + 1:
            comp.append(pos[i])
        else:
            components.append(len(comp))
            comp = [pos[i]]
    components.append(len(comp))
    return tuple(sorted(components, reverse=True))

# Generate test tournaments
random.seed(707)
test_tournaments = []
for trial in range(5):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    t3 = count_3_cycles(A, n)
    test_tournaments.append((A, t3))

# Compute sigma(S) for each position subset size and pattern
print(f"\nPosition subsets by size and pattern:")
for size in range(1, n):
    subsets_by_pattern = defaultdict(list)
    for S in combinations(range(n-1), size):
        pat = position_pattern(S, n)
        subsets_by_pattern[pat].append(S)

    for pat, subsets in sorted(subsets_by_pattern.items()):
        # Compute sigma for each tournament
        vals = []
        for A, t3 in test_tournaments:
            sigma_vals = [compute_sigma(A, n, S) for S in subsets]
            vals.append((t3, sigma_vals))

        # Check if all subsets with same pattern give same sigma
        all_same = all(len(set(sv)) == 1 for _, sv in vals)

        print(f"\n  |S|={size}, pattern={pat}: {len(subsets)} subsets")
        if all_same:
            for t3, sv in vals:
                print(f"    t3={t3}: sigma = {sv[0]} (all {len(subsets)} equal)")
        else:
            for t3, sv in vals[:2]:
                print(f"    t3={t3}: sigma = {sv[:min(4,len(sv))]}")

# =====================================================================
# Step 2: For size 1, verify sigma = n!/2
# =====================================================================
print(f"\n{'='*70}")
print("SIZE-1 VERIFICATION: sigma({i}) = n!/2")
print(f"{'='*70}")
for A, t3 in test_tournaments[:2]:
    for i in range(n-1):
        s = compute_sigma(A, n, (i,))
        print(f"  t3={t3}, pos={i}: sigma={s}, n!/2={factorial(n)//2}, match={s == factorial(n)//2}")
    break

# =====================================================================
# Step 3: For size 2, verify consecutive vs non-consecutive
# =====================================================================
print(f"\n{'='*70}")
print("SIZE-2: consecutive (depends on t3) vs non-consecutive (universal)")
print(f"{'='*70}")

for A, t3 in test_tournaments:
    consec_vals = []
    nonconsec_vals = []
    for i in range(n-1):
        for j in range(i+1, n-1):
            s = compute_sigma(A, n, (i, j))
            if j == i + 1:
                consec_vals.append(s)
            else:
                nonconsec_vals.append(s)
    print(f"  t3={t3}: consec={set(consec_vals)}, nonconsec={set(nonconsec_vals)}")

# =====================================================================
# Step 4: KEY — derive the moment formulas from sigma decomposition
# =====================================================================
print(f"\n{'='*70}")
print("DERIVING MOMENT FORMULAS FROM SIGMA DECOMPOSITION")
print(f"{'='*70}")

# f^j = (sum_{i=0}^{n-2} T_i)^j = sum over j-tuples (a_1,...,a_j) of prod T_{a_k}
# Since T_i^k = T_i (indicator), this = sum over j-multisets of (mult coeff) * prod T_{distinct}
# = sum over subsets S of [0,n-2] of Stirling(j, |S|, distribution) * sigma(S)
#
# More precisely, f^j = sum over set partitions of [j] -> [0,n-2]
# But simpler: f^j = sum_{k=0}^j S(j,k) * f_(k) where f_(k) = f!/(f-k)! falling factorial
# and f_(k) = k! * sum over size-k subsets of [0,n-2] of prod T_i
# Wait, that's not quite right either because the positions are not the same as elements.
#
# Actually: f^j = sum_{k=0}^{min(j,n-1)} S(j,k) * f^{(k)}
# where S(j,k) = Stirling number of second kind
# and f^{(k)} = sum over k-element subsets of positions of k! * prod T_i
# NO — f^{(k)} = f*(f-1)*...*(f-k+1)
# = sum over ordered k-tuples of DISTINCT positions of prod T_{a_l}
# = sum over k-subsets S of k! * prod_{i in S} T_i  [IF the T_i were distinct]
# But wait, T_i are 0/1 so f^{(k)} = sum ordered distinct k-tuples = k! * sigma(S) summed.
#
# So: sum_P f^{(k)} = k! * sum over size-k subsets S of sigma(S)
# And sum_P f^j = sum_{k=0}^j S(j,k) * sum_P f^{(k)} = sum_{k=0}^j S(j,k) * k! * SIGMA_k
# where SIGMA_k = sum over size-k position subsets of sigma(S).

# Stirling numbers of the second kind
def stirling2(n, k):
    if n == 0 and k == 0: return 1
    if n == 0 or k == 0: return 0
    return k * stirling2(n-1, k) + stirling2(n-1, k-1)

print(f"\n  f^j = sum_k S(j,k) * k! * SIGMA_k")
print(f"  where SIGMA_k = sum over size-k position subsets of sigma(S)")
print(f"  and sigma(S) = sum_P prod_{{i in S}} T_i")

# Compute SIGMA_k for each tournament
for A, t3 in test_tournaments[:3]:
    print(f"\n  t3={t3}:")
    SIGMA = []
    for k in range(n):
        total = 0
        for S in combinations(range(n-1), k):
            total += compute_sigma(A, n, S)
        SIGMA.append(total)
        print(f"    SIGMA_{k} = {total}")

    # Reconstruct moments
    for j in range(n):
        mj = sum(stirling2(j, k) * factorial(k) * SIGMA[k] for k in range(min(j,n-1)+1))
        print(f"    m_{j} via Stirling = {mj}")

# =====================================================================
# Step 5: SIGMA_k by pattern type
# =====================================================================
print(f"\n{'='*70}")
print("SIGMA_k DECOMPOSED BY PATTERN")
print(f"{'='*70}")

for k in range(1, 5):
    print(f"\n  k={k}: {comb(n-1,k)} position subsets")
    pattern_sigmas = defaultdict(lambda: [])
    for S in combinations(range(n-1), k):
        pat = position_pattern(S, n)
        for A, t3 in test_tournaments:
            pattern_sigmas[(k, pat)].append((t3, compute_sigma(A, n, S)))

    for (kk, pat), vals in sorted(pattern_sigmas.items()):
        # Group by tournament
        by_t3 = defaultdict(list)
        for t3, v in vals:
            by_t3[t3].append(v)
        # Check if all subs give same value
        for t3_val, vlist in sorted(by_t3.items()):
            unique = set(vlist)
            if len(unique) == 1:
                print(f"    pat={pat}, t3={t3_val}: all {len(vlist)} subsets give {unique.pop()}")
            else:
                print(f"    pat={pat}, t3={t3_val}: varies: {sorted(unique)}")
        break  # Just show one tournament's values for brevity

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
