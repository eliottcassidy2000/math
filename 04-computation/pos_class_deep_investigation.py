#!/usr/bin/env python3
"""
pos_class_deep_investigation.py — kind-pasteur-2026-03-21-S16b

DEEP INVESTIGATION: The (1,2,2,2,3) score class at n=5.

This is the ONLY ambiguous score class at n=5 — it contains tournaments with
H = 11, 13, or 15. The ENTIRE OCR residual (1/19 of Var(H)) comes from this class.

The score sequence (1,2,2,2,3) is SELF-COMPLEMENTARY: s_i + s_{n-1-i} = n-1 = 4.
  s_0 + s_4 = 1 + 3 = 4 ✓
  s_1 + s_3 = 2 + 2 = 4 ✓
  s_2 = 2 = (n-1)/2 ✓

This means tournaments with this score CAN be self-complementary (T ≅ T^op).

KEY FACTS FROM THE REPO:
1. H = 9 + 2*c5_dir within this class (THM-164)
2. c5_dir ∈ {1, 2, 3} (directed Hamiltonian cycle count)
3. The "cone over C₃" is in this class: {1,2,3} cycle, 4 beats all, all beat 0
4. 40/64 of the n=5 H-maximizers have this score (non-hereditary, MISTAKE-010)
5. This class has the LOWEST nonzero S₂ = 2 (weakest Cartan entanglement)

QUESTIONS TO INVESTIGATE:
A. What is the EXACT structure of the 280 tournaments with score (1,2,2,2,3)?
   How many have each H value? What is the automorphism group distribution?
B. Which tournaments are self-complementary? Do SC ones achieve max H?
C. Why is c5_dir ∈ {1,2,3} and not {0}? Is c5_dir > 0 forced by the score?
D. What is the conflict graph Ω for each H subclass?
E. How does this class GENERALIZE to n=6, n=7? What is the "POS class" at each n?
F. Does the POS class always contribute the dominant share of the OCR residual?
G. Is there a self-complementary analog of the OCR formula?

Author: kind-pasteur-2026-03-21-S16b
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import comb, factorial, gcd
from fractions import Fraction

def build_tournament(n, bits):
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

def ham_paths_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    key = (mask | (1 << w), w)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def is_self_complementary(A, n):
    """Check if T ≅ T^op (exists permutation sigma with A[sigma(i)][sigma(j)] = 1-A[i][j] for i≠j)."""
    from itertools import permutations as perms
    for sigma in perms(range(n)):
        is_anti = True
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if A[sigma[i]][sigma[j]] != (1 - A[i][j]):
                    is_anti = False
                    break
            if not is_anti:
                break
        if is_anti:
            return True
    return False

def count_directed_cycles(A, n, length):
    """Count directed cycles of given length (fixing start vertex)."""
    count = 0
    for verts in combinations(range(n), length):
        sub = [[A[verts[a]][verts[b]] for b in range(length)] for a in range(length)]
        for perm in permutations(range(1, length)):
            cycle = (0,) + perm
            if all(sub[cycle[k]][cycle[(k+1) % length]] for k in range(length)):
                count += 1
    return count

def automorphism_group_size(A, n):
    """Count automorphisms of tournament A."""
    count = 0
    for sigma in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(i+1, n):
                if A[sigma[i]][sigma[j]] != A[i][j]:
                    is_aut = False
                    break
            if not is_aut:
                break
        if is_aut:
            count += 1
    return count

print("=" * 72)
print("  THE (1,2,2,2,3) CLASS: DEEP INVESTIGATION")
print("  kind-pasteur-2026-03-21-S16b")
print("=" * 72)

n = 5
m = n*(n-1)//2  # 10
total = 1 << m   # 1024

# ============================================================
# PART 1: Complete classification
# ============================================================

print("\n" + "=" * 72)
print("  PART 1: Complete classification of (1,2,2,2,3) class")
print("=" * 72)

pos_class = []  # Tournaments with score (1,2,2,2,3)
all_data = []

for bits in range(total):
    A = build_tournament(n, bits)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    H = ham_paths_dp(A, n)
    all_data.append({'bits': bits, 'A': A, 'scores': scores, 'H': H})

    if scores == (1, 2, 2, 2, 3):
        # Compute additional invariants
        c3 = comb(n, 3) - sum(s*(s-1)//2 for s in [sum(A[i]) for i in range(n)])
        c5_dir = count_directed_cycles(A, n, 5)

        pos_class.append({
            'bits': bits, 'A': A, 'H': H, 'c3': c3, 'c5_dir': c5_dir,
        })

print(f"  Total (1,2,2,2,3) tournaments: {len(pos_class)}")

# Distribution by H
H_dist = Counter(d['H'] for d in pos_class)
c5_dist = Counter(d['c5_dir'] for d in pos_class)
print(f"  H distribution: {dict(sorted(H_dist.items()))}")
print(f"  c5_dir distribution: {dict(sorted(c5_dist.items()))}")
print(f"  H = 9 + 2*c5_dir check: {all(d['H'] == 9 + 2*d['c5_dir'] for d in pos_class)}")
print(f"  c3 values: {set(d['c3'] for d in pos_class)}")

# ============================================================
# PART 2: Self-complementary analysis
# ============================================================

print("\n" + "=" * 72)
print("  PART 2: Self-complementary tournaments in (1,2,2,2,3)")
print("=" * 72)

sc_count = 0
sc_by_H = Counter()
nsc_by_H = Counter()

for d in pos_class:
    sc = is_self_complementary(d['A'], n)
    d['sc'] = sc
    if sc:
        sc_count += 1
        sc_by_H[d['H']] += 1
    else:
        nsc_by_H[d['H']] += 1

print(f"  Self-complementary: {sc_count} / {len(pos_class)}")
print(f"  SC by H: {dict(sorted(sc_by_H.items()))}")
print(f"  NSC by H: {dict(sorted(nsc_by_H.items()))}")

# ============================================================
# PART 3: Automorphism groups
# ============================================================

print("\n" + "=" * 72)
print("  PART 3: Automorphism groups")
print("=" * 72)

# Sample a few from each H class
for h_val in sorted(H_dist.keys()):
    examples = [d for d in pos_class if d['H'] == h_val][:5]
    aut_sizes = []
    for d in examples:
        aut = automorphism_group_size(d['A'], n)
        aut_sizes.append(aut)
    print(f"  H={h_val}: |Aut| = {aut_sizes} (first {len(examples)} examples)")

# ============================================================
# PART 4: The "cone over C₃" — the canonical example
# ============================================================

print("\n" + "=" * 72)
print("  PART 4: The 'cone over C₃' structure")
print("=" * 72)

# Cone over C₃: vertices 1,2,3 form cycle, vertex 4 beats 1,2,3, vertices 1,2,3 beat 0
# In our encoding: need to find which bits give this
# Let's construct it explicitly:
# Vertex 0: loses to 1,2,3; loses/beats 4
# Vertex 4: beats 0,1,2,3
# Vertices 1,2,3: cycle 1→2→3→1, beat 0, lose to 4
# Scores: 0 has 0 or 1 wins. 4 has 4 wins. Actually:
# 0 beats nobody among 1,2,3,4? Then score(0)=0.
# Wait, (1,2,2,2,3) has min score 1, not 0.

# Let me think. Score 1: beats exactly 1 vertex. Score 3: beats exactly 3.
# For a "cone": source beats all, sink loses to all except one.
# With scores (1,2,2,2,3):
# - vertex with score 3 beats 3 others
# - vertex with score 1 beats 1 other

# Let me find a concrete example with |Aut|=3 (the Z/3 symmetry)
for d in pos_class:
    aut = automorphism_group_size(d['A'], n)
    if aut == 3:
        A = d['A']
        actual_scores = [sum(A[i]) for i in range(n)]
        print(f"  Found |Aut|=3 tournament: bits={d['bits']}, H={d['H']}, c5_dir={d['c5_dir']}")
        print(f"    Scores (per vertex): {actual_scores}")
        print(f"    Adjacency:")
        for i in range(n):
            row = [A[i][j] for j in range(n)]
            print(f"      {i}: {row}")

        # Which vertex has which score?
        for i in range(n):
            print(f"    Vertex {i}: score={actual_scores[i]}, beats {[j for j in range(n) if A[i][j]]}")
        print()
        break

# ============================================================
# PART 5: WHY is c5_dir > 0 forced? (is it?)
# ============================================================

print("\n" + "=" * 72)
print("  PART 5: Is c5_dir = 0 possible with score (1,2,2,2,3)?")
print("=" * 72)

# c5_dir = 0 means NO directed Hamiltonian cycles in the 5-vertex tournament
# Is this compatible with score (1,2,2,2,3)?
# The transitive tournament has c5_dir = 0, scores (0,1,2,3,4) — NOT (1,2,2,2,3).
# For (1,2,2,2,3): we must have c3 = 4 (Rao's formula: 10 - sum C(si,2) = 10 - (0+1+1+1+3) = 4)
# A tournament with 4 directed 3-cycles and no 5-cycle: is this possible?

c5_zero = [d for d in pos_class if d['c5_dir'] == 0]
print(f"  Tournaments with score (1,2,2,2,3) and c5_dir=0: {len(c5_zero)}")
if c5_zero:
    print(f"    H values: {set(d['H'] for d in c5_zero)}")
else:
    print(f"  c5_dir = 0 is IMPOSSIBLE with this score sequence!")
    print(f"  Minimum c5_dir = {min(d['c5_dir'] for d in pos_class)}")
    print(f"  This is a STRUCTURAL CONSTRAINT: 4 directed 3-cycles in 5 vertices force at least 1 directed 5-cycle.")

# ============================================================
# PART 6: What determines c5_dir within the class?
# ============================================================

print("\n" + "=" * 72)
print("  PART 6: What determines c5_dir within (1,2,2,2,3)?")
print("=" * 72)

# Within this class, all tournaments have same scores and same c3.
# c5_dir varies. What LOCAL structure distinguishes them?

for h_val in sorted(H_dist.keys()):
    examples = [d for d in pos_class if d['H'] == h_val]
    # Compute common-neighbor statistics
    common_out = []
    common_in = []
    for d in examples[:10]:
        A = d['A']
        # For each pair (i,j), count common out-neighbors
        co = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(n):
                    if k != i and k != j:
                        if A[i][k] and A[j][k]:
                            co += 1
        common_out.append(co)
    print(f"  H={h_val} (c5_dir={h_val//2-4}): common_out = {common_out[:5]}")

# ============================================================
# PART 7: Contribution to OCR residual
# ============================================================

print("\n" + "=" * 72)
print("  PART 7: Exact contribution of (1,2,2,2,3) to OCR residual")
print("=" * 72)

# Var(H) decomposition: total = between-class + within-class
# OCR residual = within-class / total

H_all = np.array([d['H'] for d in all_data], dtype=float)
total_var = np.var(H_all)

# Between-class: variance of class means
score_groups = defaultdict(list)
for d in all_data:
    score_groups[d['scores']].append(d['H'])

between_var = 0
within_var = 0
for s, Hs in score_groups.items():
    Hs = np.array(Hs, dtype=float)
    class_mean = np.mean(Hs)
    between_var += len(Hs) * (class_mean - np.mean(H_all))**2
    within_var += np.sum((Hs - class_mean)**2)
between_var /= total
within_var /= total

print(f"  Total Var(H) = {total_var:.6f}")
print(f"  Between-class = {between_var:.6f} ({100*between_var/total_var:.2f}%)")
print(f"  Within-class = {within_var:.6f} ({100*within_var/total_var:.2f}%)")

# Now: what fraction of within-class variance comes from (1,2,2,2,3)?
pos_Hs = np.array([d['H'] for d in pos_class], dtype=float)
pos_within = np.sum((pos_Hs - np.mean(pos_Hs))**2) / total

print(f"\n  (1,2,2,2,3) within-class contribution:")
print(f"    Var contribution = {pos_within:.6f}")
print(f"    Fraction of total within-class = {pos_within/within_var:.6f}")
print(f"    Fraction of total Var(H) = {pos_within/total_var:.6f}")

# Does this EXACTLY equal the OCR residual?
# OCR = 18/19, so 1 - OCR = 1/19 ≈ 0.05263
# within_var / total_var should equal 1/19 minus a correction for imperfect linear fit
print(f"    1/19 = {1/19:.6f}")
print(f"    within_var/total_var = {within_var/total_var:.6f}")
# Note: within-class / total ≠ 1-R^2 exactly.
# R^2 is 1 - SS_res/SS_tot where SS_res involves LINEAR regression residuals.
# Within-class variance is a DIFFERENT decomposition (ANOVA).
# For one-hot encoding of score classes, these ARE the same!
# Because the best linear predictor within a group is the group mean.
# But S₂ is NOT a one-hot encoding — it collapses some classes.

# Let's check: does S₂ UNIQUELY determine score class at n=5?
s2_to_scores = defaultdict(set)
for s, Hs in score_groups.items():
    s2_val = sum((si - 2)**2 for si in s)
    s2_to_scores[s2_val].add(s)

print(f"\n  S₂ → score class mapping at n=5:")
for s2_val in sorted(s2_to_scores.keys()):
    classes = s2_to_scores[s2_val]
    print(f"    S₂ = {s2_val}: {classes}")

# If multiple score classes share the same S₂, the linear regression on S₂
# can't distinguish them, and the ANOVA within-class variance ≠ regression residual.

# ============================================================
# PART 8: The POS property — self-complementary scores
# ============================================================

print("\n" + "=" * 72)
print("  PART 8: Self-complementary score sequences at each n")
print("=" * 72)

# A score sequence is self-complementary iff s_i + s_{n-1-i} = n-1 for all i.
# This is the "point of symmetry" (POS) condition.

for n_check in [3, 4, 5, 6, 7]:
    # Find all score sequences at this n
    m_check = n_check * (n_check - 1) // 2
    total_check = 1 << m_check

    all_scores = set()
    sc_scores = set()

    if total_check <= 2100000:
        for bits in range(total_check):
            A = build_tournament(n_check, bits)
            scores = tuple(sorted([sum(A[i]) for i in range(n_check)]))
            all_scores.add(scores)

        for s in all_scores:
            is_sc = all(s[i] + s[n_check - 1 - i] == n_check - 1 for i in range(n_check))
            if is_sc:
                sc_scores.add(s)

        print(f"\n  n={n_check}: {len(all_scores)} score sequences, {len(sc_scores)} self-complementary")
        for s in sorted(sc_scores):
            # How many tournaments?
            count = sum(1 for bits in range(total_check)
                        if tuple(sorted([sum(build_tournament(n_check, bits)[i])
                                         for i in range(n_check)])) == s)
            # How many H values?
            Hs = set()
            for bits in range(total_check):
                A = build_tournament(n_check, bits)
                if tuple(sorted([sum(A[i]) for i in range(n_check)])) == s:
                    Hs.add(ham_paths_dp(A, n_check))
            print(f"    {s}: {count} tours, H ∈ {sorted(Hs)}, ambiguous={len(Hs)>1}")
    else:
        print(f"\n  n={n_check}: (sampled — too large for exhaustive)")

# ============================================================
# PART 9: The KEY QUESTION — Why does the SC score generate the residual?
# ============================================================

print("\n" + "=" * 72)
print("  PART 9: Why does the SC score class generate the OCR residual?")
print("=" * 72)

print("""
  THESIS: The self-complementary score class has EXTRA degrees of freedom
  that non-SC classes lack. Specifically:

  1. In a non-SC score class, the score sequence is ASYMMETRIC, which
     CONSTRAINS the cycle structure. High-score vertices dominate,
     low-score vertices are dominated, and this forces a specific c5_dir.

  2. In the SC class (1,2,2,2,3), the three middle vertices all have
     score 2 = (n-1)/2. They are INTERCHANGEABLE in terms of out-degree.
     The cycle structure among them is NOT determined by scores.

  3. The number of "free" cycle choices is exactly the number of ways
     to orient the 3-cycle on the middle vertices: either 1→2→3→1 or
     1→3→2→1. This gives TWO fundamentally different structures within
     the same score class.

  4. c5_dir = {1, 2, 3} corresponds to different orientations of the
     internal structure. The spread |c5_max - c5_min| = 2 = one flip.

  PREDICTION: At n=6, the SC score classes should contribute most of the
  OCR residual. At n=7, the regular class (3,3,3,3,3,3,3) — the most
  symmetric SC class — should have the largest within-class spread.
""")

# ============================================================
# PART 10: The overlap structure — why 1/19?
# ============================================================

print("=" * 72)
print("  PART 10: Why 1/19? The arithmetic of the residual")
print("=" * 72)

# Exact computation:
# Var(H) = 285/16 at n=5
# Within-class Var from (1,2,2,2,3) = ?
# Let's compute exactly

N = total  # 1024
pos_count = len(pos_class)  # 280
pos_H = [d['H'] for d in pos_class]
pos_mean_H = Fraction(sum(pos_H), pos_count)
pos_var_H_exact = sum(Fraction(h) - pos_mean_H for h in pos_H)  # Wait, need SS

# sum(H) over POS class
sum_H_pos = sum(pos_H)
sum_H2_pos = sum(h**2 for h in pos_H)
mean_H_pos = Fraction(sum_H_pos, pos_count)
var_H_pos = Fraction(sum_H2_pos, pos_count) - mean_H_pos**2

# Within-class SS = pos_count * var_H_pos
SS_within_pos = pos_count * var_H_pos

# Total SS = N * Var(H) = N * 285/16
SS_total = N * Fraction(285, 16)

# Fraction of residual from POS class
frac = SS_within_pos / SS_total

# Also: total within-class SS = sum over all classes of class_size * class_var
SS_within_total = Fraction(0)
for s, Hs_list in score_groups.items():
    if len(set(Hs_list)) > 1:
        sz = len(Hs_list)
        s_h = sum(Hs_list)
        s_h2 = sum(h**2 for h in Hs_list)
        class_var = Fraction(s_h2, sz) - Fraction(s_h, sz)**2
        SS_within_total += sz * class_var

print(f"  (1,2,2,2,3) class: {pos_count} tournaments")
print(f"    sum(H) = {sum_H_pos}, mean(H) = {mean_H_pos}")
print(f"    Var(H within) = {var_H_pos}")
print(f"    SS_within_pos = {SS_within_pos}")
print(f"    SS_total = {SS_total}")
print(f"    SS_within_pos / SS_total = {frac} = {float(frac):.10f}")
print(f"    Total SS_within (all classes) = {SS_within_total}")
print(f"    SS_within_total / SS_total = {SS_within_total / SS_total}")
print(f"    Note: 1 - OCR = 1/19 = {Fraction(1, 19)}")
print()

# R^2 from ANOVA (one-hot score class) vs R^2 from S₂ regression
# At n=5, S₂ uniquely identifies 8 of 9 score classes.
# Two classes share S₂=8: (0,1,3,3,3) and (0,2,2,2,4) and (1,1,1,3,4).
# Wait, let me recheck.

print("  S₂ collision check at n=5:")
for s2_val in sorted(s2_to_scores.keys()):
    classes = s2_to_scores[s2_val]
    if len(classes) > 1:
        print(f"    S₂={s2_val}: COLLISION — {classes}")
        for s in classes:
            Hs = score_groups[s]
            print(f"      {s}: count={len(Hs)}, H={set(Hs)}")
    # If colliding classes have DIFFERENT H values, R^2(S₂) < R^2(score_class)
    # and the difference accounts for part of the residual.

# CRITICAL: R^2(S₂, H) = 18/19.
# R^2(ANOVA score class, H) = 1 - SS_within/SS_total
R2_anova = 1 - SS_within_total / SS_total
print(f"\n  R²(ANOVA score class) = {R2_anova}")
print(f"  R²(S₂ regression) = 18/19 = {Fraction(18, 19)}")
print(f"  Equal? {R2_anova == Fraction(18, 19)}")

# If these are equal, then:
# 1. S₂ captures ALL between-class variation (no info lost from collisions)
# 2. The residual 1/19 comes ENTIRELY from within-class variation
# 3. And the ONLY within-class variation is in (1,2,2,2,3)

if R2_anova == Fraction(18, 19):
    print(f"\n  *** S₂ regression = ANOVA EXACTLY ***")
    print(f"  This means: despite S₂ collisions, the colliding classes")
    print(f"  have the SAME mean H. So S₂ loses no information.")
    # Check
    for s2_val in sorted(s2_to_scores.keys()):
        classes = s2_to_scores[s2_val]
        if len(classes) > 1:
            means = {}
            for s in classes:
                means[s] = Fraction(sum(score_groups[s]), len(score_groups[s]))
            print(f"    S₂={s2_val}: class means = {means}")
            all_same = len(set(means.values())) == 1
            print(f"      All same mean? {all_same}")
