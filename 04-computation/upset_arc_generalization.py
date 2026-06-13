#!/usr/bin/env python3
"""
upset_arc_generalization.py — kind-pasteur-2026-03-21-S16b

THE "UPSET ARC" PRINCIPLE at n=5: When the score-1 vertex beats the score-3
vertex (an "upset"), c5_dir is maximized and H increases.

QUESTION: Does this generalize? At n=6, within an ambiguous score class,
does the presence/absence of "upset arcs" (low-score beating high-score)
predict H?

BROADER QUESTION: Is H determined by the score sequence PLUS the number
of "upset arcs" (arcs going "against" the score ordering)?

An "upset arc" for score sequence s is an arc i→j where s_i < s_j.
The total number of upset arcs is related to Kendall's tau distance
between the tournament and the "score-sorted" transitive tournament.
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import comb
from fractions import Fraction
import sys

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

def count_upset_arcs(A, n):
    """Count arcs i→j where score(i) < score(j)."""
    scores = [sum(A[i]) for i in range(n)]
    upsets = 0
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j] and scores[i] < scores[j]:
                upsets += 1
    return upsets

def count_same_score_arcs(A, n):
    """Among vertices with the SAME score, count arcs (these are "free")."""
    scores = [sum(A[i]) for i in range(n)]
    count = 0
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j] and scores[i] == scores[j]:
                count += 1
    return count

print("=" * 72)
print("  THE UPSET ARC PRINCIPLE")
print("  kind-pasteur-2026-03-21-S16b")
print("=" * 72)

# ============================================================
# n=5: Validate
# ============================================================
print("\n  n=5: upset arcs within ambiguous class (1,2,2,2,3)")

n = 5
m = n*(n-1)//2
total = 1 << m

score_groups = defaultdict(list)
for bits in range(total):
    A = build_tournament(n, bits)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    H = ham_paths_dp(A, n)
    upsets = count_upset_arcs(A, n)
    same_score = count_same_score_arcs(A, n)
    score_groups[scores].append({'H': H, 'upsets': upsets, 'same': same_score})

target = (1, 2, 2, 2, 3)
entries = score_groups[target]
H_by_upsets = defaultdict(list)
for e in entries:
    H_by_upsets[e['upsets']].append(e['H'])

print(f"\n  Upset arc distribution within (1,2,2,2,3):")
for u in sorted(H_by_upsets.keys()):
    Hs = H_by_upsets[u]
    print(f"    upsets={u}: count={len(Hs)}, H values={sorted(set(Hs))}")

# Same-score arcs
H_by_same = defaultdict(list)
for e in entries:
    H_by_same[e['same']].append(e['H'])

print(f"\n  Same-score arc distribution:")
for s in sorted(H_by_same.keys()):
    Hs = H_by_same[s]
    print(f"    same_score_arcs={s}: count={len(Hs)}, H values={sorted(set(Hs))}")

# Upset + same together
H_by_both = defaultdict(list)
for e in entries:
    H_by_both[(e['upsets'], e['same'])].append(e['H'])

print(f"\n  (upsets, same_score) -> H:")
for (u, s) in sorted(H_by_both.keys()):
    Hs = H_by_both[(u, s)]
    n_distinct = len(set(Hs))
    print(f"    ({u}, {s}): count={len(Hs)}, H={sorted(set(Hs))}, ambig={n_distinct > 1}")

# ============================================================
# n=5: R^2 of upset count predicting H (across all tournaments)
# ============================================================

print(f"\n  R^2(upsets, H) across all n=5 tournaments:")
all_upsets = []
all_H = []
for s, entries_s in score_groups.items():
    for e in entries_s:
        all_upsets.append(e['upsets'])
        all_H.append(e['H'])

all_upsets = np.array(all_upsets, dtype=float)
all_H = np.array(all_H, dtype=float)
r2_upsets = np.corrcoef(all_upsets, all_H)[0, 1]**2
print(f"  R^2(upsets, H) = {r2_upsets:.6f}")
print(f"  Compare R^2(S2, H) = {18/19:.6f}")

# ============================================================
# n=6: Same analysis
# ============================================================

print("\n" + "=" * 72)
print("  n=6: Upset arcs and H")
print("=" * 72)

n = 6
m = n*(n-1)//2
total = 1 << m

score_groups_6 = defaultdict(list)
for bits in range(total):
    A = build_tournament(n, bits)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    H = ham_paths_dp(A, n)
    upsets = count_upset_arcs(A, n)
    same = count_same_score_arcs(A, n)
    score_groups_6[scores].append({'H': H, 'upsets': upsets, 'same': same})

# Check ambiguous classes
ambiguous_6 = {s: entries for s, entries in score_groups_6.items()
               if len(set(e['H'] for e in entries)) > 1}

print(f"\n  Ambiguous score classes at n=6: {len(ambiguous_6)}")

# For the worst class (1,2,2,3,3,4) — 6 H values
worst = (1, 2, 2, 3, 3, 4)
if worst in ambiguous_6:
    entries = ambiguous_6[worst]
    print(f"\n  Class {worst}: {len(entries)} tournaments, H in {sorted(set(e['H'] for e in entries))}")

    # Does upset count predict H within this class?
    H_by_u = defaultdict(list)
    for e in entries:
        H_by_u[e['upsets']].append(e['H'])

    print(f"  Upsets -> H:")
    for u in sorted(H_by_u.keys()):
        Hs = H_by_u[u]
        print(f"    upsets={u}: count={len(Hs)}, H={sorted(set(Hs))}")

    # R^2 within this class
    up_arr = np.array([e['upsets'] for e in entries], dtype=float)
    h_arr = np.array([e['H'] for e in entries], dtype=float)
    if np.std(up_arr) > 0:
        r2_within = np.corrcoef(up_arr, h_arr)[0, 1]**2
        print(f"\n  R^2(upsets, H) within {worst} = {r2_within:.6f}")

# R^2 across all n=6
all_up = np.array([e['upsets'] for entries in score_groups_6.values() for e in entries], dtype=float)
all_h = np.array([e['H'] for entries in score_groups_6.values() for e in entries], dtype=float)
all_s2 = []
for s, entries in score_groups_6.items():
    s2 = sum((si - 2.5)**2 for si in s)
    for _ in entries:
        all_s2.append(s2)
all_s2 = np.array(all_s2, dtype=float)

r2_up_6 = np.corrcoef(all_up, all_h)[0, 1]**2
r2_s2_6 = np.corrcoef(all_s2, all_h)[0, 1]**2

# Joint: S2 + upsets
X = np.column_stack([np.ones(len(all_h)), all_s2, all_up])
coefs = np.linalg.lstsq(X, all_h, rcond=None)[0]
pred = X @ coefs
ss_res = np.sum((all_h - pred)**2)
ss_tot = np.sum((all_h - np.mean(all_h))**2)
r2_joint = 1 - ss_res / ss_tot

print(f"\n  Global R^2 comparisons at n=6:")
print(f"    R^2(S2, H)         = {r2_s2_6:.6f}")
print(f"    R^2(upsets, H)     = {r2_up_6:.6f}")
print(f"    R^2(S2+upsets, H)  = {r2_joint:.6f}")

# Also: same-score arcs
all_same = np.array([e['same'] for entries in score_groups_6.values() for e in entries], dtype=float)
X2 = np.column_stack([np.ones(len(all_h)), all_s2, all_same])
coefs2 = np.linalg.lstsq(X2, all_h, rcond=None)[0]
pred2 = X2 @ coefs2
r2_joint2 = 1 - np.sum((all_h - pred2)**2) / ss_tot

print(f"    R^2(S2+same_arcs, H) = {r2_joint2:.6f}")

# Triple
X3 = np.column_stack([np.ones(len(all_h)), all_s2, all_up, all_same])
coefs3 = np.linalg.lstsq(X3, all_h, rcond=None)[0]
pred3 = X3 @ coefs3
r2_triple = 1 - np.sum((all_h - pred3)**2) / ss_tot
print(f"    R^2(S2+upsets+same, H) = {r2_triple:.6f}")

# Ideal: ANOVA (score class)
# Already know R^2(ANOVA) = 12/13 for S2 regression, but ANOVA gives more
# Let's compute ANOVA with score+upsets as group variable
score_upset_groups = defaultdict(list)
for s, entries in score_groups_6.items():
    for e in entries:
        score_upset_groups[(s, e['upsets'])].append(e['H'])

SS_between_su = Fraction(0)
grand_mean = Fraction(sum(all_h), len(all_h))
for key, Hs in score_upset_groups.items():
    class_mean = Fraction(sum(Hs), len(Hs))
    SS_between_su += len(Hs) * (class_mean - grand_mean)**2

SS_total_6 = sum((Fraction(h) - grand_mean)**2 for h in all_h)
R2_anova_su = float(SS_between_su / SS_total_6)
print(f"    R^2(ANOVA score+upsets) = {R2_anova_su:.6f}")

# Score + same arcs
score_same_groups = defaultdict(list)
for s, entries in score_groups_6.items():
    for e in entries:
        score_same_groups[(s, e['same'])].append(e['H'])

SS_between_ss = Fraction(0)
for key, Hs in score_same_groups.items():
    class_mean = Fraction(sum(Hs), len(Hs))
    SS_between_ss += len(Hs) * (class_mean - grand_mean)**2

R2_anova_ss = float(SS_between_ss / SS_total_6)
print(f"    R^2(ANOVA score+same)   = {R2_anova_ss:.6f}")
