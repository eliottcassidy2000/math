#!/usr/bin/env python3
"""
pos_class_fast.py — kind-pasteur-2026-03-21-S16b (fast version)

Focus on the critical questions about the (1,2,2,2,3) class.
Skip slow operations (automorphism groups, self-complementary check for all).
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import comb
from fractions import Fraction
import sys

# Flush output immediately
sys.stdout.reconfigure(line_buffering=True)

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

n = 5
m = n*(n-1)//2
total = 1 << m

print("=" * 72)
print("  POS CLASS (1,2,2,2,3) — FAST VERSION")
print("=" * 72)

# Build all n=5 tournaments
all_data = []
pos_class = []

for bits in range(total):
    A = build_tournament(n, bits)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    H = ham_paths_dp(A, n)
    all_data.append({'bits': bits, 'scores': scores, 'H': H})
    if scores == (1, 2, 2, 2, 3):
        # c5_dir
        c5d = 0
        A_flat = A
        for perm in permutations(range(1, 5)):
            cycle = (0,) + perm
            if all(A_flat[cycle[k]][cycle[(k+1)%5]] for k in range(5)):
                c5d += 1
        c3 = comb(5, 3) - sum(s*(s-1)//2 for s in [sum(A[i]) for i in range(5)])
        pos_class.append({'bits': bits, 'H': H, 'c3': c3, 'c5d': c5d, 'A': A})

print(f"\nTotal (1,2,2,2,3) tournaments: {len(pos_class)}")
H_dist = Counter(d['H'] for d in pos_class)
c5_dist = Counter(d['c5d'] for d in pos_class)
print(f"H distribution: {dict(sorted(H_dist.items()))}")
print(f"c5_dir distribution: {dict(sorted(c5_dist.items()))}")
print(f"H = 9 + 2*c5_dir: {all(d['H'] == 9 + 2*d['c5d'] for d in pos_class)}")
print(f"c3 = {set(d['c3'] for d in pos_class)} (all same)")

# ============================================================
# S2 collision check
# ============================================================
print("\n" + "=" * 72)
print("  S2 COLLISION CHECK")
print("=" * 72)

score_groups = defaultdict(list)
for d in all_data:
    score_groups[d['scores']].append(d['H'])

s2_to_scores = defaultdict(set)
s2_to_mean_H = defaultdict(list)
for s in score_groups:
    s2 = sum((si - 2)**2 for si in s)
    s2_to_scores[s2].add(s)
    s2_to_mean_H[s2].append((s, Fraction(sum(score_groups[s]), len(score_groups[s]))))

print(f"\nS2 -> score class mapping:")
for s2 in sorted(s2_to_scores.keys()):
    classes = s2_to_scores[s2]
    means = [(s, m) for s, m in s2_to_mean_H[s2]]
    if len(classes) > 1:
        print(f"  S2={s2}: COLLISION {classes}")
        for s, mean in means:
            print(f"    {s}: mean H = {mean} = {float(mean):.4f}, count = {len(score_groups[s])}")
        same_mean = len(set(m for _, m in means)) == 1
        print(f"    Same mean? {same_mean}")
    else:
        s = list(classes)[0]
        print(f"  S2={s2}: {s}, mean H = {float(means[0][1]):.4f}, count = {len(score_groups[s])}")

# ============================================================
# ANOVA R^2 vs regression R^2
# ============================================================
print("\n" + "=" * 72)
print("  ANOVA R^2 vs S2 REGRESSION R^2")
print("=" * 72)

N = total
H_all = [d['H'] for d in all_data]
sum_H = sum(H_all)
sum_H2 = sum(h**2 for h in H_all)
grand_mean = Fraction(sum_H, N)
SS_total = Fraction(sum_H2, 1) - N * grand_mean**2

# ANOVA SS_between
SS_between = Fraction(0)
for s, Hs in score_groups.items():
    class_mean = Fraction(sum(Hs), len(Hs))
    SS_between += len(Hs) * (class_mean - grand_mean)**2

SS_within = SS_total - SS_between
R2_anova = SS_between / SS_total

print(f"  SS_total = {SS_total}")
print(f"  SS_between (score classes) = {SS_between}")
print(f"  SS_within = {SS_within}")
print(f"  R^2(ANOVA) = {R2_anova} = {float(R2_anova):.10f}")
print(f"  R^2(S2 regression) = 18/19 = {float(Fraction(18,19)):.10f}")
print(f"  Equal? {R2_anova == Fraction(18, 19)}")
print()

# The within-class SS comes ENTIRELY from (1,2,2,2,3)
pos_Hs = [d['H'] for d in pos_class]
pos_mean = Fraction(sum(pos_Hs), len(pos_Hs))
SS_within_pos = sum(Fraction(h) - pos_mean for h in pos_Hs)  # Hmm, need squares
SS_within_pos = sum((Fraction(h) - pos_mean)**2 for h in pos_Hs)

print(f"  SS_within from (1,2,2,2,3): {SS_within_pos}")
print(f"  SS_within total: {SS_within}")
print(f"  Fraction: {SS_within_pos / SS_within}")
print(f"  The (1,2,2,2,3) class accounts for {float(SS_within_pos / SS_within)*100:.1f}% of the residual")

# ============================================================
# THE KEY IDENTITY: Why exactly 1/19?
# ============================================================
print("\n" + "=" * 72)
print("  WHY 1/19?")
print("=" * 72)

print(f"\n  1/19 = SS_within / SS_total = {SS_within} / {SS_total}")
print(f"  = {SS_within} / {SS_total}")
print(f"  Simplified: {SS_within / SS_total}")

# Factor everything
print(f"\n  SS_total = {SS_total}")
print(f"  SS_within = {SS_within}")
print(f"  Ratio = {Fraction(SS_within, SS_total)}")

# What is SS_within exactly?
# 280 tournaments, H in {11, 13, 15}
# H counts: 11 -> ?, 13 -> ?, 15 -> ?
# Mean H in class = ?
print(f"\n  (1,2,2,2,3) details:")
print(f"    Count: {len(pos_class)}")
print(f"    H distribution: {dict(sorted(H_dist.items()))}")
print(f"    Mean H = {pos_mean} = {float(pos_mean):.6f}")

# SS_within = n_11*(11-mean)^2 + n_13*(13-mean)^2 + n_15*(15-mean)^2
n_11 = H_dist[11]
n_13 = H_dist[13]
n_15 = H_dist[15]
print(f"    n_11={n_11}, n_13={n_13}, n_15={n_15}")
print(f"    SS = {n_11}*(11-{float(pos_mean):.2f})^2 + {n_13}*(13-{float(pos_mean):.2f})^2 + {n_15}*(15-{float(pos_mean):.2f})^2")

ss_check = n_11 * (Fraction(11) - pos_mean)**2 + n_13 * (Fraction(13) - pos_mean)**2 + n_15 * (Fraction(15) - pos_mean)**2
print(f"    = {ss_check}")
print(f"    Match: {ss_check == SS_within_pos}")

# ============================================================
# WHAT DETERMINES H WITHIN THE CLASS?
# ============================================================
print("\n" + "=" * 72)
print("  WHAT DETERMINES c5_dir?")
print("=" * 72)

# Group by c5_dir and look at adjacency patterns
for c5_val in sorted(c5_dist.keys()):
    examples = [d for d in pos_class if d['c5d'] == c5_val]
    print(f"\n  c5_dir = {c5_val} (H = {9 + 2*c5_val}): {len(examples)} tournaments")

    # Look at per-vertex out-sets
    # In score (1,2,2,2,3): one vertex has score 1, one has score 3, three have score 2
    # The arrangement of arcs among the three score-2 vertices determines the cycle structure

    # Compute: does the score-1 vertex beat the score-3 vertex?
    beat_pattern = Counter()
    for d in examples:
        A = d['A']
        actual_scores = [sum(A[i]) for i in range(5)]
        v1 = actual_scores.index(1)  # Could be ambiguous if multiple have score 1
        # Find score-1 and score-3 vertices
        score_1_verts = [i for i in range(5) if actual_scores[i] == 1]
        score_3_verts = [i for i in range(5) if actual_scores[i] == 3]
        score_2_verts = [i for i in range(5) if actual_scores[i] == 2]

        # Does score-1 beat score-3?
        if len(score_1_verts) == 1 and len(score_3_verts) == 1:
            v1 = score_1_verts[0]
            v3 = score_3_verts[0]
            beats = A[v1][v3]
            beat_pattern[beats] += 1

    print(f"    Score-1 beats Score-3: {dict(beat_pattern)}")

# ============================================================
# THE SELF-COMPLEMENTARY STRUCTURE
# ============================================================
print("\n" + "=" * 72)
print("  SELF-COMPLEMENTARY SCORE: s_i + s_{n-1-i} = n-1")
print("=" * 72)

print(f"\n  Score (1,2,2,2,3): self-complementary check:")
s = (1, 2, 2, 2, 3)
for i in range(5):
    print(f"    s[{i}] + s[{4-i}] = {s[i]} + {s[4-i]} = {s[i] + s[4-i]}", end="")
    print(f" {'= 4 ✓' if s[i] + s[4-i] == 4 else '≠ 4 ✗'}")

print(f"\n  Other SC scores at n=5:")
all_scores_5 = set()
for d in all_data:
    all_scores_5.add(d['scores'])

for score in sorted(all_scores_5):
    is_sc = all(score[i] + score[4-i] == 4 for i in range(5))
    if is_sc:
        Hs = score_groups[score]
        n_ambig = len(set(Hs))
        print(f"    {score}: SC=True, count={len(Hs)}, distinct H={n_ambig}, H={sorted(set(Hs))}")

print(f"\n  Non-SC scores at n=5:")
for score in sorted(all_scores_5):
    is_sc = all(score[i] + score[4-i] == 4 for i in range(5))
    if not is_sc:
        Hs = score_groups[score]
        n_ambig = len(set(Hs))
        print(f"    {score}: SC=False, count={len(Hs)}, distinct H={n_ambig}, H={sorted(set(Hs))}")

# ============================================================
# THE GENERALIZATION TO n=6
# ============================================================
print("\n" + "=" * 72)
print("  GENERALIZATION: SC SCORES AND AMBIGUITY AT n=6")
print("=" * 72)

n6 = 6
m6 = n6*(n6-1)//2
total6 = 1 << m6

score_groups_6 = defaultdict(list)
for bits in range(total6):
    A = build_tournament(n6, bits)
    scores = tuple(sorted([sum(A[i]) for i in range(n6)]))
    H = ham_paths_dp(A, n6)
    score_groups_6[scores].append(H)

print(f"\n  n=6: {len(score_groups_6)} score classes")
print(f"\n  SC vs non-SC ambiguity:")

sc_ambig = 0
sc_total = 0
nsc_ambig = 0
nsc_total = 0

for score in sorted(score_groups_6.keys()):
    is_sc = all(score[i] + score[5-i] == 5 for i in range(6))
    Hs = score_groups_6[score]
    n_distinct = len(set(Hs))
    ambig = n_distinct > 1

    if is_sc:
        sc_total += 1
        if ambig:
            sc_ambig += 1
            print(f"    SC  {score}: {n_distinct} H values, spread {max(Hs)-min(Hs)}")
    else:
        nsc_total += 1
        if ambig:
            nsc_ambig += 1
            print(f"    NSC {score}: {n_distinct} H values, spread {max(Hs)-min(Hs)}")

print(f"\n  SC classes: {sc_ambig}/{sc_total} ambiguous")
print(f"  NSC classes: {nsc_ambig}/{nsc_total} ambiguous")

# Compute SS contributions
SS_total_6 = Fraction(0)
grand_mean_6 = Fraction(sum(sum(Hs) for Hs in score_groups_6.values()), total6)
for Hs in score_groups_6.values():
    for h in Hs:
        SS_total_6 += (Fraction(h) - grand_mean_6)**2

SS_within_sc = Fraction(0)
SS_within_nsc = Fraction(0)
for score, Hs in score_groups_6.items():
    is_sc = all(score[i] + score[5-i] == 5 for i in range(6))
    class_mean = Fraction(sum(Hs), len(Hs))
    ss = sum((Fraction(h) - class_mean)**2 for h in Hs)
    if is_sc:
        SS_within_sc += ss
    else:
        SS_within_nsc += ss

SS_within_6 = SS_within_sc + SS_within_nsc

print(f"\n  Variance decomposition:")
print(f"    SS_total = {SS_total_6}")
print(f"    SS_within (total) = {SS_within_6}")
print(f"    SS_within (SC classes) = {SS_within_sc}")
print(f"    SS_within (NSC classes) = {SS_within_nsc}")
print(f"    SC fraction of residual = {float(SS_within_sc / SS_within_6):.4f}")
print(f"    NSC fraction of residual = {float(SS_within_nsc / SS_within_6):.4f}")
print(f"    R^2(ANOVA) = {float(1 - SS_within_6 / SS_total_6):.10f}")
print(f"    Expected = 12/13 = {float(Fraction(12,13)):.10f}")
