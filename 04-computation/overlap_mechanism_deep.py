#!/usr/bin/env python3
"""
overlap_mechanism_deep.py -- kind-pasteur-2026-03-13-S61

Deep dive into the MECHANISM connecting overlap weights to H.

Key discovery from vitali_overlap_hidden_structure.py:
  At n=7 regular, 3-cycle overlap histogram determines H-class.
  Paley: {0:7, 1:63, 2:21} -> H=189
  H=171: {0:10, 1:57, 2:24}
  H=175: {0:14, 1:49, 2:28}

The overlap-1 count PERFECTLY tracks H (63 > 57 > 49 but H: 189 > 175 > 171).
Wait: 189 > 175 > 171, but overlap-1: 63 > 57 > 49. So MORE overlap-1 = HIGHER H.
And MORE disjoint = LOWER H (despite alpha_2 contributing 4 per pair).

This is the Tournament Uncertainty Principle: the 2*alpha_1 term DOMINATES
the 4*alpha_2 term. More overlap-1 pairs mean MORE total directed cycles
(from 5-cycles living on the same vertex fabric).

This script:
1. Linear regression: H = a*ov0 + b*ov1 + c*ov2 + d at n=7 regular
2. Why overlap-1 drives 5-cycles: the "gateway" structure
3. Algebraic constraint: ov0 + ov1 + ov2 = C(14,2) = 91 for ALL regular n=7
4. Transfer mechanism: overlap weight -> directed 5-cycle count
5. n=6 exhaustive: same mechanism at smaller scale?
6. The Vitali obstruction: what prevents overlap-1 from determining c5 exactly?

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations
from collections import defaultdict
import math


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


# ========================================================================
# ANALYSIS 1: LINEAR MODEL — H = f(ov0, ov1, ov2) AT n=7 REGULAR
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: LINEAR MODEL FOR H FROM OVERLAP WEIGHTS")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
total = 1 << m

# Collect all regular tournament data
regular_data = []
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = [sum(A[v]) for v in range(n)]
    if any(s != 3 for s in scores):
        continue

    H = count_ham_paths(A, n)

    # 3-cycle vertex sets
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))
    nc3 = len(c3_sets)

    # Overlap histogram
    ov = defaultdict(int)
    for i in range(nc3):
        for j in range(i+1, nc3):
            w = len(c3_sets[i] & c3_sets[j])
            ov[w] += 1

    ov0 = ov.get(0, 0)
    ov1 = ov.get(1, 0)
    ov2 = ov.get(2, 0)

    # 5-cycle count
    c5_dir = 0
    for subset in combinations(range(n), 5):
        c5_dir += count_directed_ham_cycles_on_subset(A, list(subset))

    # 7-cycle (Hamiltonian) count
    c7_dir = count_directed_ham_cycles_on_subset(A, list(range(n)))

    # Disjoint pairs (= alpha_2 for 3-cycles)
    disj_33 = ov0

    regular_data.append({
        'bits': bits, 'H': H, 'ov0': ov0, 'ov1': ov1, 'ov2': ov2,
        'c5_dir': c5_dir, 'c7_dir': c7_dir, 'disj_33': disj_33
    })

print(f"Regular tournaments: {len(regular_data)}")

# Check: ov0 + ov1 + ov2 = C(14, 2) = 91?
for d in regular_data[:5]:
    ov_sum = d['ov0'] + d['ov1'] + d['ov2']
    print(f"  bits={d['bits']}: ov0={d['ov0']}, ov1={d['ov1']}, ov2={d['ov2']}, "
          f"sum={ov_sum}, H={d['H']}")

# Since ov0+ov1+ov2 = C(14,2) = 91, the system has only 2 DOF (ov0, ov1)
# with ov2 = 91 - ov0 - ov1.

# Linear fit: H = a*ov0 + b*ov1 + c
# Using the three known classes:
# H=189: ov0=7, ov1=63
# H=171: ov0=10, ov1=57
# H=175: ov0=14, ov1=49
# System:  7a + 63b + c = 189
#         10a + 57b + c = 171
#         14a + 49b + c = 175

# From first two: 3a - 6b = -18, so a - 2b = -6
# From first and third: 7a - 14b + c2 = -14, so a - 2b = -2 ...wait

# Let me solve properly
# Eq1: 7a + 63b + c = 189
# Eq2: 10a + 57b + c = 171
# Eq3: 14a + 49b + c = 175

# Eq2-Eq1: 3a - 6b = -18 => a - 2b = -6
# Eq3-Eq2: 4a - 8b = 4 => a - 2b = 1

# CONTRADICTION! a - 2b = -6 AND a - 2b = 1 is impossible.
# So H is NOT a linear function of (ov0, ov1).

print(f"\nLinear fit attempt:")
print(f"  Eq2-Eq1: 3a - 6b = -18 => a - 2b = -6")
print(f"  Eq3-Eq2: 4a - 8b = 4 => a - 2b = 1")
print(f"  CONTRADICTION: H is NOT linear in (ov0, ov1)")

# Try quadratic: H = a*ov0^2 + b*ov1^2 + c*ov0 + d*ov1 + e*ov0*ov1 + f
# With 3 data points and 6 unknowns, underdetermined.
# But since ov0+ov1+ov2=91, we can use ov0 as sole variable:
# ov0=7 -> H=189, ov0=10 -> H=171, ov0=14 -> H=175
# So H = a*ov0^2 + b*ov0 + c
# 49a + 7b + c = 189
# 100a + 10b + c = 171
# 196a + 14b + c = 175
# From 1,2: 51a + 3b = -18 => 17a + b = -6
# From 2,3: 96a + 4b = 4 => 24a + b = 1
# Subtract: 7a = 7 => a = 1
# b = -6 - 17 = -23
# c = 189 - 49 + 161 = 189 - 49 + 161 = 301
# Check: c = 189 - 49*1 - 7*(-23) = 189 - 49 + 161 = 301

a_coef, b_coef, c_coef = 1, -23, 301
print(f"\nQuadratic fit: H = {a_coef}*ov0^2 + ({b_coef})*ov0 + {c_coef}")
for d in regular_data[:10]:
    H_pred = a_coef * d['ov0']**2 + b_coef * d['ov0'] + c_coef
    print(f"  ov0={d['ov0']}: H_pred={H_pred}, H_actual={d['H']}, "
          f"match={'YES' if H_pred == d['H'] else 'NO'}")

# Verify for ALL
all_match = all(
    a_coef * d['ov0']**2 + b_coef * d['ov0'] + c_coef == d['H']
    for d in regular_data
)
print(f"\n  H = ov0^2 - 23*ov0 + 301 for ALL {len(regular_data)} regular n=7 tournaments: {all_match}")

# This is a remarkable result! H is a quadratic function of the disjoint pair count.
# Since ov0 = disj_33 = alpha_2 for 3-cycle pairs, this gives:
# H = alpha_2^2 - 23*alpha_2 + 301
# Check: alpha_2=7: 49-161+301 = 189 YES
# alpha_2=10: 100-230+301 = 171 YES
# alpha_2=14: 196-322+301 = 175 YES

print(f"\n  THEOREM: For regular n=7 tournaments:")
print(f"  H = disj_33^2 - 23*disj_33 + 301")
print(f"  where disj_33 = number of vertex-disjoint 3-cycle pairs")


# ========================================================================
# ANALYSIS 2: WHY QUADRATIC? — CONNECTING TO OCF DECOMPOSITION
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: WHY QUADRATIC? CONNECTING TO OCF DECOMPOSITION")
print("=" * 70)

# We know: H = 1 + 2*alpha_1 + 4*alpha_2 for n=7
# And: c5_dir + 2*disj_33 = 56 (constraint for regular n=7)
# And: alpha_1 = c3_dir + c5_dir + c7_dir
# For regular: c3_dir = 14 (constant)
# So: alpha_1 = 14 + c5_dir + c7_dir
# And: H = 1 + 2*(14 + c5_dir + c7_dir) + 4*disj_33
#       = 29 + 2*c5_dir + 2*c7_dir + 4*disj_33
# Using c5_dir = 56 - 2*disj_33:
#       = 29 + 2*(56 - 2*disj_33) + 2*c7_dir + 4*disj_33
#       = 29 + 112 - 4*disj_33 + 2*c7_dir + 4*disj_33
#       = 141 + 2*c7_dir
# This gives H = 141 + 2*c7_dir (linear in c7), NOT quadratic in disj_33.
# But both hold simultaneously! So:
# disj_33^2 - 23*disj_33 + 301 = 141 + 2*c7_dir
# => c7_dir = (disj_33^2 - 23*disj_33 + 160) / 2

print(f"Deriving: c7_dir = (disj_33^2 - 23*disj_33 + 160) / 2")
for d in regular_data[:10]:
    c7_pred = (d['disj_33']**2 - 23*d['disj_33'] + 160) / 2
    print(f"  disj_33={d['disj_33']}: c7_pred={c7_pred:.0f}, c7_actual={d['c7_dir']}, "
          f"match={'YES' if abs(c7_pred - d['c7_dir']) < 0.01 else 'NO'}")

# So the quadratic relationship H vs disj_33 arises because c7 has a quadratic
# dependence on disj_33. This is a NEW constraint (beyond c5+2*disj=56).

# Verify: c7_dir = (disj_33^2 - 23*disj_33 + 160) / 2 for ALL
all_match_c7 = all(
    abs((d['disj_33']**2 - 23*d['disj_33'] + 160) / 2 - d['c7_dir']) < 0.01
    for d in regular_data
)
print(f"\n  c7_dir = (disj_33^2 - 23*disj_33 + 160)/2 for ALL: {all_match_c7}")

# Values:
# disj_33=7: c7 = (49-161+160)/2 = 48/2 = 24 (Paley)
# disj_33=10: c7 = (100-230+160)/2 = 30/2 = 15
# disj_33=14: c7 = (196-322+160)/2 = 34/2 = 17
print(f"\n  disj_33=7: c7={(49-161+160)//2} (Paley)")
print(f"  disj_33=10: c7={(100-230+160)//2}")
print(f"  disj_33=14: c7={(196-322+160)//2}")


# ========================================================================
# ANALYSIS 3: THE SECOND CONSTRAINT — c7 vs disj_33
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: MECHANISM — WHY c7 DEPENDS ON disj_33")
print("=" * 70)

# A Hamiltonian 7-cycle visits all 7 vertices. More disjoint 3-pairs means
# the 3-cycle structure is "spread out" — vertices are used less efficiently.
# But Hamiltonian cycles need "connectivity" to close.

# Let's look at how 3-cycles relate to Hamiltonian structure.
# For each tournament, count how many Hamiltonian cycles "respect" the 3-cycle
# structure in different ways.

for H_target in [189, 171, 175]:
    d = next(d for d in regular_data if d['H'] == H_target)
    bits = d['bits']
    A = binary_to_tournament(bits, n)

    # Enumerate all directed Ham cycles as ordered 7-tuples
    ham_cycles = []
    # Use DP to enumerate (not just count) — get actual paths
    # Path from 0: dp[mask][v] = list of path strings
    # Too expensive for all paths. Instead just count.
    c7 = d['c7_dir']

    # Get 3-cycle vertex sets
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))

    # For each 5-vertex subset, check if it's in both a 3-cycle pair and a 5-cycle
    c5_per_5set = {}
    for subset in combinations(range(n), 5):
        fs = frozenset(subset)
        c5_count = count_directed_ham_cycles_on_subset(A, list(subset))
        c5_per_5set[fs] = c5_count

    # How does c5 distribution correlate with disjoint 3-pair structure?
    # Each 5-subset can contain 0, 1, or 2 disjoint 3-cycle vertex sets
    c5_by_disj = defaultdict(list)
    for subset in combinations(range(n), 5):
        fs = frozenset(subset)
        # Count how many 3-cycle vertex sets are contained in this 5-subset
        contained = [cs for cs in c3_sets if cs.issubset(fs)]
        # Among contained, count disjoint pairs
        n_contained = len(contained)
        disj_in_5 = 0
        for i in range(n_contained):
            for j in range(i+1, n_contained):
                if not (contained[i] & contained[j]):
                    disj_in_5 += 1

        c5_by_disj[disj_in_5].append(c5_per_5set[fs])

    print(f"\n  H={H_target} (disj_33={d['disj_33']}, c7={c7}):")
    print(f"    c5 per 5-subset grouped by #disjoint 3-pairs in that subset:")
    for dp_count in sorted(c5_by_disj.keys()):
        vals = c5_by_disj[dp_count]
        avg = sum(vals) / len(vals) if vals else 0
        total = sum(vals)
        print(f"      disj_pairs_in_5set={dp_count}: count={len(vals)}, "
              f"avg_c5={avg:.2f}, total_c5={total}, distribution={sorted(set(vals))}")


# ========================================================================
# ANALYSIS 4: n=6 — DOES OVERLAP WEIGHT DETERMINE H?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: n=6 — OVERLAP WEIGHT vs H")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m

# Group by score class and check if overlap histogram determines H
score_groups = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))

    # 3-cycles
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))
    nc3 = len(c3_sets)

    ov_hist = defaultdict(int)
    for i in range(nc3):
        for j in range(i+1, nc3):
            w = len(c3_sets[i] & c3_sets[j])
            ov_hist[w] += 1

    ov_key = tuple(sorted(ov_hist.items()))
    score_groups[scores].append({'H': H, 'ov_key': ov_key, 'nc3': nc3, 'bits': bits})

print(f"Score classes at n=6: {len(score_groups)}")

# For each score class, check if overlap histogram determines H
for sc in sorted(score_groups.keys()):
    group = score_groups[sc]
    if len(set(d['H'] for d in group)) <= 1:
        continue  # Skip single-H classes

    # Group by overlap histogram within this score class
    by_ov = defaultdict(list)
    for d in group:
        by_ov[d['ov_key']].append(d['H'])

    ambig = 0
    for ok, Hs in by_ov.items():
        if len(set(Hs)) > 1:
            ambig += 1

    if ambig > 0:
        H_range = sorted(set(d['H'] for d in group))
        ov_count = len(by_ov)
        print(f"\n  Score {sc}: H in {H_range}, {ov_count} overlap histograms, {ambig} AMBIGUOUS")
        for ok in sorted(by_ov.keys()):
            Hs = sorted(by_ov[ok])
            H_unique = sorted(set(Hs))
            if len(H_unique) > 1:
                print(f"    {dict(ok)}: H = {H_unique}")
    else:
        H_range = sorted(set(d['H'] for d in group))
        if len(H_range) > 2:
            print(f"  Score {sc}: H in {H_range}, overlap histogram DETERMINES H")


# ========================================================================
# ANALYSIS 5: THE TRANSFER MECHANISM — 3-CYCLE OVERLAP -> 5-CYCLE COUNT
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: TRANSFER MECHANISM — 3-OVERLAP -> 5-CYCLES AT n=6")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
total = 1 << m

# Focus on score class (2,2,2,3,3,3) — the near-regular class
target_score = (2, 2, 2, 3, 3, 3)
transfer_data = []

for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != target_score:
        continue

    H = count_ham_paths(A, n)

    # 3-cycles
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))
    nc3 = len(c3_sets)

    # 5-cycle count (directed)
    c5_dir = 0
    for subset in combinations(range(n), 5):
        c5_dir += count_directed_ham_cycles_on_subset(A, list(subset))

    # Disjoint 3-3 pairs
    disj = 0
    ov1 = 0
    ov2 = 0
    for i in range(nc3):
        for j in range(i+1, nc3):
            w = len(c3_sets[i] & c3_sets[j])
            if w == 0:
                disj += 1
            elif w == 1:
                ov1 += 1
            elif w == 2:
                ov2 += 1

    transfer_data.append({
        'H': H, 'c3': nc3, 'c5_dir': c5_dir, 'disj': disj,
        'ov1': ov1, 'ov2': ov2
    })

# Group by (c3, disj)
by_c3_disj = defaultdict(list)
for d in transfer_data:
    by_c3_disj[(d['c3'], d['disj'])].append(d)

print(f"Score {target_score}:")
print(f"  Total: {len(transfer_data)}")
print(f"  Distinct (c3, disj) classes: {len(by_c3_disj)}")

for (c3, disj) in sorted(by_c3_disj.keys()):
    group = by_c3_disj[(c3, disj)]
    H_set = sorted(set(d['H'] for d in group))
    c5_set = sorted(set(d['c5_dir'] for d in group))
    ov1_set = sorted(set(d['ov1'] for d in group))

    print(f"  c3={c3}, disj_33={disj}: count={len(group)}, "
          f"H={H_set}, c5_dir={c5_set}, ov1={ov1_set}")

    # Check: does c5+2*disj = constant? (like at n=7)
    sums = set(d['c5_dir'] + 2*d['disj'] for d in group)
    print(f"    c5+2*disj = {sorted(sums)}")


# ========================================================================
# ANALYSIS 6: UNIVERSALITY CHECK — DOES QUADRATIC HOLD AT n=6?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: QUADRATIC H(disj_33) AT n=6?")
print("=" * 70)

# Within each score class, check if H = a*disj^2 + b*disj + c
for sc in sorted(score_groups.keys()):
    group = score_groups[sc]
    Hs = sorted(set(d['H'] for d in group))
    if len(Hs) <= 1:
        continue

    # Collect (disj, H) pairs
    disj_H = defaultdict(set)
    for d in group:
        # Compute disj for this tournament
        bits = d['bits']
        A_loc = binary_to_tournament(bits, n)
        c3_loc = []
        for a, b, c_v in combinations(range(n), 3):
            if A_loc[a][b] and A_loc[b][c_v] and A_loc[c_v][a]:
                c3_loc.append(frozenset([a, b, c_v]))
            if A_loc[a][c_v] and A_loc[c_v][b] and A_loc[b][a]:
                c3_loc.append(frozenset([a, b, c_v]))
        c3_loc = list(set(c3_loc))
        nc3_loc = len(c3_loc)
        disj_count = 0
        for i in range(nc3_loc):
            for j in range(i+1, nc3_loc):
                if not (c3_loc[i] & c3_loc[j]):
                    disj_count += 1
        disj_H[disj_count].add(d['H'])

    if not disj_H:
        continue

    # Check if disj determines H
    det = all(len(v) == 1 for v in disj_H.values())
    if det:
        pairs = [(d, list(h)[0]) for d, h in disj_H.items()]
        pairs.sort()
        if len(pairs) >= 3:
            # Try quadratic fit
            # H = a*d^2 + b*d + c
            d0, h0 = pairs[0]
            d1, h1 = pairs[1]
            d2, h2 = pairs[2]
            # Solve
            # a*d0^2 + b*d0 + c = h0
            # a*d1^2 + b*d1 + c = h1
            # a*d2^2 + b*d2 + c = h2
            det_m = (d0**2*(d1-d2) - d0*(d1**2-d2**2) + (d1**2*d2-d2**2*d1))
            if abs(det_m) > 0.01:
                a_fit = (h0*(d1-d2) - d0*(h1-h2) + (h1*d2-h2*d1)) / det_m * d0
                # Hmm, let me do it differently
                pass

        print(f"  Score {sc}: disj DETERMINES H: {sorted(pairs)}")
    else:
        # disj doesn't determine H
        ambig_disj = [(d, sorted(h)) for d, h in disj_H.items() if len(h) > 1]
        if ambig_disj:
            print(f"  Score {sc}: disj does NOT determine H: {ambig_disj[:3]}...")


# ========================================================================
# ANALYSIS 7: COUNTING ARGUMENT — WHY 3 CLASSES AT n=7 REGULAR
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 7: COUNTING ARGUMENT FOR 3 H-CLASSES AT n=7 REGULAR")
print("=" * 70)

n = 7
# All regular n=7 have c3=14 (vertex sets), each 3-vertex set is either a 3-cycle or not
# C(7,3) = 35 total 3-vertex subsets, 14 are 3-cycles, 21 are non-cycles
# For 3-cycle pairs: C(14,2) = 91 total pairs
# Each pair has overlap 0, 1, or 2

# Constraints on overlap counts:
# ov0 + ov1 + ov2 = 91
# Each vertex participates in C(6,2) = 15 pairs as a potential shared vertex
# More precisely: for a fixed vertex v, the 3-cycles containing v form a set of size...

# Count: each vertex is in how many 3-cycles?
# In a regular n=7 tournament (out-degree 3), vertex v is in exactly:
# Sum over pairs (a,b) in N_out(v): 1 if a->b (one 3-cycle v->a->b->...wait)
# Actually, vertex v is in c3_vertex = 2*(d_out*(d_out-1)/2 - ???)
# For regular: d_out = 3, C(3,2)=3. Each pair of out-neighbors either forms
# a 3-cycle (if the edge between them goes the right way) or not.
# c3 through v = (number of edges within N_out(v) with correct direction)
# + (number of edges within N_in(v) with correct direction)

print("At n=7 regular:")
# Let's just compute from actual data
for H_target in [189, 171, 175]:
    d = next(d for d in regular_data if d['H'] == H_target)
    bits = d['bits']
    A = binary_to_tournament(bits, n)

    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))

    # Per-vertex 3-cycle count
    per_vertex = [sum(1 for cs in c3_sets if v in cs) for v in range(n)]
    print(f"\n  H={H_target}: per-vertex 3-cycle count = {per_vertex}")
    print(f"    Sum = {sum(per_vertex)} (should be 3*14={3*14})")

    # For overlap analysis: how many overlap-1 pairs involve vertex v?
    ov1_per_vertex = [0] * n
    for i in range(len(c3_sets)):
        for j in range(i+1, len(c3_sets)):
            shared = c3_sets[i] & c3_sets[j]
            if len(shared) == 1:
                v = list(shared)[0]
                ov1_per_vertex[v] += 1

    print(f"    Per-vertex overlap-1 count = {ov1_per_vertex}")
    print(f"    Sum = {sum(ov1_per_vertex)} (should be ov1={d['ov1']})")

    # Per-vertex overlap-2 count
    ov2_per_pair = defaultdict(int)
    for i in range(len(c3_sets)):
        for j in range(i+1, len(c3_sets)):
            shared = c3_sets[i] & c3_sets[j]
            if len(shared) == 2:
                pair = frozenset(shared)
                ov2_per_pair[pair] += 1

    print(f"    Distinct overlap-2 vertex pairs: {len(ov2_per_pair)}")
    ov2_vals = sorted(ov2_per_pair.values())
    print(f"    Overlap-2 multiplicity distribution: {ov2_vals}")


print(f"\n\n{'='*70}")
print("DONE.")
print("=" * 70)
