#!/usr/bin/env python3
"""
Self-Converse Maximizer: Towards a Proof
==========================================
Key observations from paired_class_H_analysis.py:

1. Global H-maximizer is ALWAYS SC at n=3,4,5,6,7 (exhaustive/sampling)
2. SC mean H >> NSC mean H at all tested n
3. Score regularity (low variance) strongly correlates with high H
4. SC tournaments tend toward regular scores

Can we PROVE that max H is achieved by an SC tournament?

Approach: Use the averaging argument.
  H(T) + H(T^op) = 2*H(T) (since H(T) = H(T^op) by path reversal).
  For NSC: T ≇ T^op, so T and T^op are in different classes.
  The average H over {T, T^op} = H(T).
  No averaging gain — the symmetry argument doesn't help directly.

Better approach: Use the independence polynomial.
  I(Omega(T), x) = I(Omega(T^op), x) (Theorem 1 of THM-022).
  So H(T) = I(Omega(T), 2) is the SAME for T and T^op.
  This means H is constant on T,T^op pairs.

Key question: Among all tournaments with score sequence s, which has max H?
Is it always SC? Does the SC structure create more odd cycles?

Also: Connection to the FLIP operation from the tiling model.
Grid-symmetric (blueself) tilings correspond to specific SC tournaments.
Do they maximize H within their score sequence class?

kind-pasteur-2026-03-06-S18
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, hamiltonian_path_count,
                             opposite_tournament, find_odd_cycles)
from itertools import permutations
from collections import defaultdict

def canonical_form(T):
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def is_self_converse(T):
    Top = opposite_tournament(T)
    return canonical_form(T) == canonical_form(Top)

def score_sequence(T):
    n = len(T)
    return tuple(sorted([sum(T[i]) for i in range(n)], reverse=True))

# ══════════════════════════════════════════════════════════════
# Within each score sequence class: who maximizes H?
# ══════════════════════════════════════════════════════════════
print("=" * 70)
print("SC MAXIMIZER WITHIN SCORE SEQUENCE CLASSES")
print("=" * 70)

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    # Group by (score_sequence, isomorphism_class)
    by_scores = defaultdict(list)  # score_seq -> list of (H, is_SC, canonical)

    seen_classes = {}  # canonical -> (H, SC)
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        ss = score_sequence(T)
        cf = canonical_form(T)
        if cf not in seen_classes:
            sc = is_self_converse(T)
            seen_classes[cf] = (h, sc, ss)

    # Group classes by score sequence
    score_classes = defaultdict(list)
    for cf, (h, sc, ss) in seen_classes.items():
        score_classes[ss].append((h, sc, cf))

    print(f"\nn={n}:")
    for ss in sorted(score_classes.keys()):
        classes = sorted(score_classes[ss], key=lambda x: -x[0])
        max_h = classes[0][0]
        max_sc = classes[0][1]
        sc_count = sum(1 for _, sc, _ in classes if sc)
        nsc_count = len(classes) - sc_count

        if len(classes) > 1:
            h_vals = [(h, 'SC' if sc else 'NSC') for h, sc, _ in classes]
            print(f"  scores={ss}: {len(classes)} classes, {sc_count} SC, {nsc_count} NSC")
            print(f"    H values: {h_vals}")
            print(f"    MAX H={max_h} is {'SC' if max_sc else 'NSC'}")
        else:
            print(f"  scores={ss}: 1 class, {'SC' if classes[0][1] else 'NSC'}, H={classes[0][0]}")

# ══════════════════════════════════════════════════════════════
# Cycle count comparison: SC vs NSC with same score sequence
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("CYCLE COUNT: SC vs NSC WITHIN SCORE CLASSES (n=6)")
print("=" * 70)

n = 6
m = n * (n - 1) // 2

# For each tournament, count odd cycles by length
def cycle_counts(T):
    cycles = find_odd_cycles(T)
    counts = defaultdict(int)
    for c in cycles:
        counts[len(c)] += 1
    return dict(counts)

score_data = defaultdict(list)
seen = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    cf = canonical_form(T)
    if cf in seen:
        continue
    seen[cf] = True
    h = hamiltonian_path_count(T)
    ss = score_sequence(T)
    sc = is_self_converse(T)
    cc = cycle_counts(T)
    score_data[ss].append((h, sc, cc))

for ss in sorted(score_data.keys()):
    entries = score_data[ss]
    if len(entries) > 1:
        print(f"\nscores={ss}:")
        for h, sc, cc in sorted(entries, key=lambda x: -x[0]):
            label = 'SC ' if sc else 'NSC'
            c3 = cc.get(3, 0)
            c5 = cc.get(5, 0)
            total = sum(cc.values())
            print(f"  {label} H={h:3d}, c3={c3:2d}, c5={c5:2d}, total_cycles={total:3d}")

# ══════════════════════════════════════════════════════════════
# Is there a simple inequality: SC + most regular => max H?
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("PROOF ATTEMPT: H(T) <= H(T_regular) for regular score seq")
print("=" * 70)

# For the most regular score sequence at n=6: (3,3,3,2,2,2)
# There are multiple SC classes. Do they ALL have H > all NSC with same scores?
regular_classes = [(h, sc, cc) for h, sc, cc in score_data[(3,3,3,2,2,2)]]
print(f"\nScore (3,3,3,2,2,2): {len(regular_classes)} classes")
sc_h = sorted([h for h, sc, _ in regular_classes if sc], reverse=True)
nsc_h = sorted([h for h, sc, _ in regular_classes if not sc], reverse=True)
print(f"  SC H values: {sc_h}")
print(f"  NSC H values: {nsc_h}")
if nsc_h:
    print(f"  All SC > all NSC? {min(sc_h) > max(nsc_h)}")
else:
    print(f"  (No NSC classes with this score sequence)")

# For score (4,3,3,2,2,1)
mixed_classes = [(h, sc, cc) for h, sc, cc in score_data[(4,3,3,2,2,1)]]
print(f"\nScore (4,3,3,2,2,1): {len(mixed_classes)} classes")
sc_h2 = sorted([h for h, sc, _ in mixed_classes if sc], reverse=True)
nsc_h2 = sorted([h for h, sc, _ in mixed_classes if not sc], reverse=True)
print(f"  SC H values: {sc_h2}")
print(f"  NSC H values: {nsc_h2}")
if sc_h2 and nsc_h2:
    print(f"  All SC > all NSC? {min(sc_h2) > max(nsc_h2)}")

# ══════════════════════════════════════════════════════════════
# KEY TEST: For each score sequence, is max H always SC?
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("FOR EACH SCORE SEQUENCE: IS MAX H ALWAYS SC?")
print("=" * 70)

for n in [4, 5, 6]:
    m = n * (n - 1) // 2
    score_h = defaultdict(lambda: {'sc_max': 0, 'nsc_max': 0, 'sc_min': float('inf'), 'nsc_min': float('inf')})

    seen2 = {}
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        cf = canonical_form(T)
        if cf in seen2:
            continue
        seen2[cf] = True
        h = hamiltonian_path_count(T)
        ss = score_sequence(T)
        sc = is_self_converse(T)

        if sc:
            score_h[ss]['sc_max'] = max(score_h[ss]['sc_max'], h)
            score_h[ss]['sc_min'] = min(score_h[ss]['sc_min'], h)
        else:
            score_h[ss]['nsc_max'] = max(score_h[ss]['nsc_max'], h)
            score_h[ss]['nsc_min'] = min(score_h[ss]['nsc_min'], h)

    print(f"\nn={n}:")
    all_sc_max = True
    for ss in sorted(score_h.keys()):
        d = score_h[ss]
        has_both = d['sc_max'] > 0 and d['nsc_max'] > 0
        if has_both:
            overall_max = max(d['sc_max'], d['nsc_max'])
            max_is_sc = d['sc_max'] >= d['nsc_max']
            if not max_is_sc:
                all_sc_max = False
            print(f"  {ss}: SC max={d['sc_max']}, NSC max={d['nsc_max']} -> max is {'SC' if max_is_sc else 'NSC!!!'}")
    print(f"  MAX H always SC within score class? {all_sc_max}")

print("\nDone.")
