#!/usr/bin/env python3
"""
vitali_deep_structure.py -- kind-pasteur-2026-03-13-S61

THE DEEP VITALI STRUCTURE: Hidden dimensions in tournaments

THM-163 revealed that H has effective dimension ~1.8 DOF on tournament space.
This means the "tournament manifold" (the image of H on {-1,+1}^m) is
essentially 2-dimensional, with:
  - DOF 1: Score variance (captures ~97% of H variation)
  - DOF 2: Cycle disjointness structure (captures ~3% via H_4)

The Vitali connection deepens:
  - The "measurable" part of tournament structure = score sequence
  - The "non-measurable" part = cycle interactions within a score class
  - The quotient by score equivalence is FINITE (like Vitali mod Q on [0,1])
  - But the LIMIT n->inf makes the quotient increasingly complex

This script explores:
1. The 2D "H-manifold" embedding for n=5,6
2. How many DOF survive in the Fourier expansion at each n
3. The Vitali coset structure: score classes as analogues of Q-cosets
4. The tournament "uncertainty principle": score regularity vs cycle control
5. Connection to overlap weight via the {0,1,2} structure

Author: kind-pasteur-2026-03-13-S61
"""

import math
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
    if n <= 1:
        return 1
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


def compute_H2(A, n):
    c2 = math.factorial(n - 2) / (2 ** (n - 2))
    scores = [sum(A[v]) for v in range(n)]
    half = (n - 1) / 2
    Z = [-2 * (s - half)**2 + half for s in scores]
    return c2 * sum(Z)


def score_variance(A, n):
    scores = [sum(A[v]) for v in range(n)]
    half = (n - 1) / 2
    return sum((s - half)**2 for s in scores) / n


# ========================================================================
# ANALYSIS 1: THE 2D H-MANIFOLD
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: THE 2D H-MANIFOLD")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    total = 1 << m
    EH = math.factorial(n) / (2 ** (n - 1))

    # Compute (score_variance, H_4, H) for all tournaments
    data = []
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        var_s = score_variance(A, n)
        H2 = compute_H2(A, n)
        H4 = H - EH - H2
        data.append((var_s, H4, H))

    # Group by (var_s, H4) to see if this DETERMINES H
    manifold = defaultdict(set)
    for var_s, h4, H in data:
        manifold[(round(var_s, 4), round(h4, 4))].add(H)

    # Check: is (var_s, H4) a unique predictor of H?
    all_unique = all(len(v) == 1 for v in manifold.values())

    print(f"\nn={n}: H-manifold analysis")
    print(f"  Total tournaments: {total}")
    print(f"  (var_s, H_4) determines H uniquely? {all_unique}")
    print(f"  Number of (var_s, H_4) pairs: {len(manifold)}")
    print(f"  Number of distinct H values: {len(set(H for _, _, H in data))}")

    if all_unique:
        print(f"  H is a FUNCTION of (score_variance, H_4) -- confirmed!")
        print(f"  The H-manifold is EXACTLY 2-dimensional.")

    # Since H = EH + H2(var_s) + H4, and H2 is determined by var_s,
    # H is determined by (var_s, H4). This is TRIVIALLY TRUE by construction!
    # The real question: is H4 determined by anything simpler?

    # How many distinct H4 values exist for each score class?
    by_var = defaultdict(set)
    for var_s, h4, H in data:
        by_var[round(var_s, 4)].add(round(h4, 4))

    print(f"\n  Score variance | #H4 values | H4 range")
    print(f"  ---------------+------------+---------")
    for var_s in sorted(by_var.keys()):
        h4_vals = sorted(by_var[var_s])
        if len(h4_vals) == 1:
            print(f"  {var_s:>15.4f} | {len(h4_vals):>10d} | {h4_vals[0]:.1f}")
        else:
            print(f"  {var_s:>15.4f} | {len(h4_vals):>10d} | {h4_vals[0]:.1f} to {h4_vals[-1]:.1f}")


# ========================================================================
# ANALYSIS 2: THE TOURNAMENT UNCERTAINTY PRINCIPLE
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: THE TOURNAMENT UNCERTAINTY PRINCIPLE")
print("=" * 70)

print("""
CLAIM: There is a TOURNAMENT UNCERTAINTY PRINCIPLE:
  Score regularity (small Var_s) and cycle control (specific H_4) are
  partially incompatible. You can optimize one at the cost of the other.

  - Regular tournaments (Var_s = 0) have the MOST favorable H_2,
    but H_4 is FIXED (= 0 at n=5, may vary at larger n).
  - Irregular tournaments sacrifice H_2 but gain more control over H_4.

This is analogous to the Heisenberg uncertainty principle: position
(score exactness) and momentum (cycle structure) cannot both be
precisely controlled.

It's also the Vitali analogy: trying to "choose" both a specific
score sequence AND a specific cycle structure from each equivalence
class is like trying to construct a measurable Vitali set.
""")

# Test: at each score variance level, what's the range of H_4?
for n in [5, 6]:
    m = n * (n - 1) // 2
    total = 1 << m
    EH = math.factorial(n) / (2 ** (n - 1))

    by_var = defaultdict(list)
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        var_s = round(score_variance(A, n), 4)
        H2 = compute_H2(A, n)
        H4 = round(H - EH - H2, 4)
        by_var[var_s].append(H4)

    print(f"\nn={n}: Score variance vs H_4 spread")
    print(f"  {'Var_s':>8s} | {'#tours':>6s} | {'H4_range':>12s} | {'H4_spread':>10s} | {'H4_std':>8s}")
    print(f"  {'':->8s}-+-{'':->6s}-+-{'':->12s}-+-{'':->10s}-+-{'':->8s}")

    for var_s in sorted(by_var.keys()):
        h4s = by_var[var_s]
        spread = max(h4s) - min(h4s)
        std_h4 = (sum((h - sum(h4s)/len(h4s))**2 for h in h4s) / len(h4s)) ** 0.5
        range_str = f"{min(h4s):.1f} to {max(h4s):.1f}" if spread > 0 else f"{h4s[0]:.1f}"
        print(f"  {var_s:>8.4f} | {len(h4s):>6d} | {range_str:>12s} | {spread:>10.1f} | {std_h4:>8.4f}")


# ========================================================================
# ANALYSIS 3: THE {0,1,2} INTERPRETATION REVISITED
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: THE {0,1,2} STRUCTURE IN OVERLAP WEIGHTS")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m
EH = math.factorial(n) / (2 ** (n - 1))

print(f"""
At n={n}: The three levels {0, 1, 2} appear naturally in three places:

1. SCORE VALUES: s_v ranges from 0 to n-1 = {n-1}
   The "ternary" restriction to {0, 1, 2} captures n=3 exactly.

2. OVERLAP WEIGHTS: Two 3-cycles share 0, 1, or 2 vertices.
   W(C_i, C_j) in {{0, 1, 2}} (never 3, since |C|=3).

3. FOURIER DEGREES: Only degrees {0, 2, 4} appear (mapped to {0, 1, 2}).
   The degree-2k term captures k-th order interactions.

The VITALI STRUCTURE maps these three views:
  Score level (fine)   ->  Fourier degree 2 (coarse)
  Overlap weight (fine) -> Fourier degree 4 (finer)
  Full tournament (finest) -> All degrees (complete)

Each coarsening is a "Vitali quotient" — grouping by equivalence class.
""")

# Show the overlap weight distribution for all n=5 tournaments
c3_overlap_by_H = defaultdict(lambda: defaultdict(int))
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)

    # 3-cycle vertex sets
    c3_sets = []
    for a, b, c in combinations(range(n), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset([a, b, c]))
    c3_sets = list(set(c3_sets))

    # Overlap distribution
    for i in range(len(c3_sets)):
        for j in range(i + 1, len(c3_sets)):
            ov = len(c3_sets[i] & c3_sets[j])
            c3_overlap_by_H[H][ov] += 1

print(f"\n  3-Cycle overlap distribution by H value:")
print(f"  {'H':>5s} | {'ov=0':>8s} | {'ov=1':>8s} | {'ov=2':>8s} | {'total':>8s}")
print(f"  {'':->5s}-+-{'':->8s}-+-{'':->8s}-+-{'':->8s}-+-{'':->8s}")
for H in sorted(c3_overlap_by_H.keys()):
    ov_dist = c3_overlap_by_H[H]
    total_pairs = sum(ov_dist.values())
    # Average over tournaments with this H
    # Actually we accumulated across tournaments, need to normalize
    # Count tournaments with this H
    n_tours = sum(1 for bits in range(total)
                  if count_ham_paths(binary_to_tournament(bits, n), n) == H)
    print(f"  {H:>5d} | {ov_dist.get(0,0)/n_tours:>8.2f} | {ov_dist.get(1,0)/n_tours:>8.2f} | "
          f"{ov_dist.get(2,0)/n_tours:>8.2f} | {total_pairs/n_tours:>8.2f}")


# ========================================================================
# ANALYSIS 4: THE HIDDEN THIRD DIMENSION — WHAT H_4 REALLY MEASURES
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: WHAT H_4 REALLY MEASURES")
print("=" * 70)

# At n=5, only score class (1,2,2,2,3) has H variation.
# H_4 = {-1, 0, 3} for this class.
# What distinguishes H=11 from H=13 from H=15?

n = 5
m = n * (n - 1) // 2
total = 1 << m
EH = math.factorial(n) / (2 ** (n - 1))

print(f"\nn={n}: Score class (1,2,2,2,3)")

for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != (1, 2, 2, 2, 3):
        continue

    H = count_ham_paths(A, n)
    H2 = compute_H2(A, n)
    H4 = round(H - EH - H2, 4)

    # Count directed cycles by length
    c3_dir = 0
    c5_dir = 0
    c3_sets = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            c3_dir += 1
            c3_sets.append(frozenset([a, b, c]))
        if A[a][c] and A[c][b] and A[b][a]:
            c3_dir += 1
            c3_sets.append(frozenset([a, b, c]))
    c3_sets = list(set(c3_sets))

    # 5-cycle count
    sub_A = [[A[i][j] for j in range(n)] for i in range(n)]
    dp = {}
    dp[(1, 0)] = 1
    for mask_val in range(1, 1 << n):
        for v in range(n):
            if not (mask_val & (1 << v)):
                continue
            key = (mask_val, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask_val & (1 << w):
                    continue
                if sub_A[v][w]:
                    nk = (mask_val | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << n) - 1
    for v in range(1, n):
        if (full, v) in dp and sub_A[v][0]:
            c5_dir += dp[(full, v)]

    # Common-neighbor matrix AA^T
    AA_T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                AA_T[i][j] = sum(A[i][k] * A[j][k] for k in range(n))

    # AA^T off-diagonal values
    off_diag = []
    for i in range(n):
        for j in range(i+1, n):
            off_diag.append(AA_T[i][j])
    aat_var = sum((x - sum(off_diag)/len(off_diag))**2 for x in off_diag) / len(off_diag)

    # Total directed cycles
    total_dir = c3_dir + c5_dir

    # Only print first few of each H value
    if H == 11:
        label = "FIRST"
    elif H == 13:
        label = "FIRST"
    elif H == 15:
        label = "FIRST"
    else:
        label = ""

    if label:
        actual_scores = [sum(A[v]) for v in range(n)]
        print(f"\n  H={H} (H_4={H4}): bits={bits}")
        print(f"    Scores (by vertex): {actual_scores}")
        print(f"    c3_dir={c3_dir}, c5_dir={c5_dir}, total_dir={total_dir}")
        print(f"    c3_sets={len(c3_sets)}")
        print(f"    AA^T off-diag: {off_diag}, var={aat_var:.4f}")

        # Adjacency description
        adj_desc = []
        for i in range(n):
            beats = [j for j in range(n) if A[i][j]]
            adj_desc.append(f"    v{i}(s={actual_scores[i]}): beats {beats}")
        for line in adj_desc:
            print(line)


# Summary
print(f"\n  Summary across all (1,2,2,2,3) tournaments:")
by_H = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != (1, 2, 2, 2, 3):
        continue

    H = count_ham_paths(A, n)

    c3_dir = 0
    c5_dir = 0
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            c3_dir += 1
        if A[a][c] and A[c][b] and A[b][a]:
            c3_dir += 1

    dp = {}
    dp[(1, 0)] = 1
    for mask_val in range(1, 1 << n):
        for v in range(n):
            if not (mask_val & (1 << v)):
                continue
            key = (mask_val, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask_val & (1 << w):
                    continue
                if A[v][w]:
                    nk = (mask_val | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    c5_dir = 0
    full = (1 << n) - 1
    for v in range(1, n):
        if (full, v) in dp and A[v][0]:
            c5_dir += dp[(full, v)]

    total_dir = c3_dir + c5_dir

    off_diag = []
    for i in range(n):
        for j in range(i+1, n):
            off_diag.append(sum(A[i][k] * A[j][k] for k in range(n)))
    aat_var = sum((x - sum(off_diag)/len(off_diag))**2 for x in off_diag) / len(off_diag)

    by_H[H].append({
        'c3_dir': c3_dir, 'c5_dir': c5_dir, 'total_dir': total_dir,
        'aat_var': aat_var
    })

print(f"\n  {'H':>5s} | {'count':>6s} | {'c3_dir':>7s} | {'c5_dir':>7s} | {'total':>6s} | {'AA^T var':>10s}")
print(f"  {'':->5s}-+-{'':->6s}-+-{'':->7s}-+-{'':->7s}-+-{'':->6s}-+-{'':->10s}")
for H in sorted(by_H.keys()):
    group = by_H[H]
    c3 = set(d['c3_dir'] for d in group)
    c5 = set(d['c5_dir'] for d in group)
    td = set(d['total_dir'] for d in group)
    av = set(round(d['aat_var'], 4) for d in group)
    print(f"  {H:>5d} | {len(group):>6d} | {str(sorted(c3)):>7s} | {str(sorted(c5)):>7s} | "
          f"{str(sorted(td)):>6s} | {str(sorted(av)):>10s}")


print("\n" + "=" * 70)
print("DONE.")
