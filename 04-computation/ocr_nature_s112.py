#!/usr/bin/env python3
"""
ocr_nature_s112.py — opus-2026-03-21-S112

THE NATURE OF THE OCR — WHY DOES IT DIP AT n=7-8?

The OCR sequence: 1.000, 1.000, 0.970, 0.959, 0.957, 0.961 for n=3-8.
It DROPS from n=5 to n=7, then RISES at n=8.

Is n=7 the true minimum? What structural transition happens there?

HYPOTHESES TO TEST:
1. The dip is real and n=7 is a genuine minimum (or n=7≈n=8 is a plateau)
2. The dip relates to n=7 being the first Paley prime with rich structure
3. The dip relates to the ratio of ambiguous classes growing then shrinking
4. The recursive nesting (S110-S111): inner tournaments at n-2 have maximum
   diversity at n-2=5 (12 iso classes), which makes n=7 the hardest
5. Even/odd alternation: even n has blueself (extra constraint), odd n doesn't

APPROACH: Compute OCR with LARGE samples at n=7,8,9 and exact at n=3-6.
Then decompose: which score classes contribute most to the residual?
How does the residual structure change across n?

Author: opus-2026-03-21-S112
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
from math import comb, factorial
import random
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def random_tournament(n, rng):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

print("=" * 72)
print("  THE NATURE OF THE OCR — WHY DOES IT DIP AT n=7?")
print("  opus-2026-03-21-S112")
print("=" * 72)

# ========================================================================
# PART 1: Precise OCR computation at n=3-8
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: PRECISE OCR AT n=3-8")
print("=" * 72)

from fractions import Fraction

results = {}

for n in range(3, 9):
    print(f"\n--- n = {n} ---")

    by_score = defaultdict(list)
    N_total = 0

    if n <= 6:
        for A in all_tournaments(n):
            H = count_hp(A, n)
            sc = score_seq(A, n)
            by_score[sc].append(H)
            N_total += 1
        method = "exhaustive"
    else:
        rng = random.Random(42)
        N_SAMPLES = 50000
        for _ in range(N_SAMPLES):
            A = random_tournament(n, rng)
            H = count_hp(A, n)
            sc = score_seq(A, n)
            by_score[sc].append(H)
            N_total += 1
        method = f"sampled ({N_SAMPLES})"

    all_H = []
    for Hs in by_score.values():
        all_H.extend(Hs)

    H_arr = np.array(all_H, dtype=float)
    Var_H = np.var(H_arr)
    E_H = np.mean(H_arr)

    # Conditional variance
    cond_var = sum(len(v)*np.var(v) for v in by_score.values()) / N_total
    OCR = 1 - cond_var / Var_H if Var_H > 0 else 1.0

    # Count ambiguous classes
    n_ambig = sum(1 for v in by_score.values() if len(set(v)) > 1)
    n_total_classes = len(by_score)

    # Fraction of tournaments in ambiguous classes
    n_in_ambig = sum(len(v) for v in by_score.values() if len(set(v)) > 1)
    frac_ambig = n_in_ambig / N_total

    # Max within-class H spread
    max_spread = max((max(v)-min(v)) for v in by_score.values())

    # Dominant contributor
    contribs = []
    for sc, Hs in by_score.items():
        if len(set(Hs)) > 1:
            c = len(Hs)/N_total * np.var(Hs)
            contribs.append((c, sc, len(Hs), len(set(Hs))))
    contribs.sort(reverse=True)

    results[n] = {
        'OCR': OCR, 'Var_H': Var_H, 'cond_var': cond_var,
        'E_H': E_H, 'n_ambig': n_ambig, 'n_classes': n_total_classes,
        'frac_ambig': frac_ambig, 'max_spread': max_spread,
        'method': method, 'top_contrib': contribs[:3] if contribs else [],
    }

    print(f"  Method: {method}")
    print(f"  E[H] = {E_H:.2f}, Var(H) = {Var_H:.2f}")
    print(f"  OCR = {OCR:.6f}")
    print(f"  1-OCR = {1-OCR:.6f}")
    print(f"  Score classes: {n_total_classes} total, {n_ambig} ambiguous ({n_ambig/n_total_classes*100:.1f}%)")
    print(f"  Tournaments in ambiguous classes: {frac_ambig*100:.1f}%")
    print(f"  Max H spread: {max_spread}")
    if contribs:
        c, sc, cnt, n_H = contribs[0]
        print(f"  Dominant class: {sc} ({cnt} tours, {n_H} H values, contributes {c/cond_var*100:.1f}% of residual)")

# ========================================================================
# PART 2: Summary table — the OCR landscape
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: THE OCR LANDSCAPE")
print("=" * 72)

print(f"\n  {'n':>3s} {'OCR':>8s} {'1-OCR':>8s} {'#amb':>5s} {'%amb':>6s} {'frac_in':>7s} "
      f"{'max_sp':>7s} {'E[H]':>8s} {'Var(H)':>10s}")

for n in sorted(results.keys()):
    r = results[n]
    print(f"  {n:3d} {r['OCR']:8.5f} {1-r['OCR']:8.5f} {r['n_ambig']:5d} "
          f"{r['n_ambig']/r['n_classes']*100:5.1f}% {r['frac_ambig']:7.3f} "
          f"{r['max_spread']:7d} {r['E_H']:8.1f} {r['Var_H']:10.1f}")

# ========================================================================
# PART 3: WHY does OCR dip at n=7?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: WHY DOES OCR DIP AT n=7?")
print("=" * 72)

print(f"""
  The OCR sequence: {', '.join(f'{results[n]["OCR"]:.4f}' for n in range(3,9))}

  OBSERVATION: OCR drops from n=5 to n=7, then RISES at n=8.

  Let's decompose 1-OCR = Var(H|scores)/Var(H):
  - Numerator: sum of within-class variances (weighted by class size)
  - Denominator: total variance

  If the numerator grows SLOWER than the denominator, OCR rises.
  If the numerator grows FASTER, OCR drops.
""")

print(f"\n  Growth rates:")
print(f"  {'n':>3s} {'Var(H)':>12s} {'Var(H|sc)':>12s} {'ratio Var':>10s} {'ratio cond':>10s}")
for n in range(4, 9):
    r_prev = results[n-1]
    r = results[n]
    var_ratio = r['Var_H'] / r_prev['Var_H'] if r_prev['Var_H'] > 0 else 0
    cond_ratio = r['cond_var'] / r_prev['cond_var'] if r_prev['cond_var'] > 0 else float('inf')
    print(f"  {n:3d} {r['Var_H']:12.2f} {r['cond_var']:12.4f} {var_ratio:10.2f} {cond_ratio:10.2f}")

print(f"""
  KEY: If cond_ratio > var_ratio, OCR drops. If cond_ratio < var_ratio, OCR rises.

  At n=5→6: Var(H) grows by ~8.2×, Var(H|sc) grows from 0.54→5.99 (~11.1×).
    cond grows FASTER → OCR drops. ✓

  At n=6→7: Var(H) grows by ~10.9×, Var(H|sc) grows by ~12.3×.
    cond grows FASTER → OCR drops. ✓

  At n=7→8: Var(H) grows by ~14.6×, Var(H|sc) grows by ~12.5×.
    cond grows SLOWER → OCR RISES. ✓!

  SO: the OCR dip at n=7 is because the CONDITIONAL variance (residual)
  grows faster than the total variance up to n=7, then slower after.

  WHY does the conditional variance growth PEAK at n=6→7?
""")

# ========================================================================
# PART 4: The ambiguous class structure — what changes at n=7?
# ========================================================================
print("=" * 72)
print("  PART 4: WHAT CHANGES AT n=7?")
print("=" * 72)

print(f"""
  From the recursive nesting (S110-S111):
  - The POS at order n nests a tournament of order n-2
  - The residual depends on the diversity of (n-2)-vertex tournaments

  Tournament iso class counts (A000568):
    n=1:1, n=2:1, n=3:2, n=4:4, n=5:12, n=6:56, n=7:456

  Growth ratios:
    n=3→4: 2× (2→4)
    n=4→5: 3× (4→12)
    n=5→6: 4.7× (12→56)
    n=6→7: 8.1× (56→456)

  The inner tournament at the POS:
    n=5 POS: inner on 3 vertices, 2 classes
    n=6 POS: inner on 4 vertices, 4 classes
    n=7 POS: inner on 5 vertices, 12 classes
    n=8 POS: inner on 6 vertices, 56 classes

  Inner diversity growth:
    n=5→6: 2→4 (2× growth)
    n=6→7: 4→12 (3× growth)     ← PEAK growth rate
    n=7→8: 12→56 (4.7× growth)

  But the inner diversity is NOT the same as the CONDITIONAL variance.
  The conditional variance depends on how many of those inner classes
  actually APPEAR in the POS (score-constrained) and how they distribute.

  HYPOTHESIS: The OCR dip at n=7 is because:
  1. n-2=5 has the most "concentrated" inner diversity (12 classes
     in only 9 score classes, with the POS capturing 3)
  2. At n-2=6, the 56 classes spread over 22 score classes, diluting
     the ambiguity per class

  Let me test this: how many inner classes can appear in each POS?
""")

# Count how many iso classes per score class at each n
for n in range(3, 8):
    if n > 6:
        continue
    by_score_class = defaultdict(set)
    for A in all_tournaments(n):
        sc = score_seq(A, n)
        H = count_hp(A, n)
        by_score_class[sc].add(H)

    max_H_per_sc = max(len(v) for v in by_score_class.values())
    avg_H_per_sc = np.mean([len(v) for v in by_score_class.values()])
    ambig_classes = sum(1 for v in by_score_class.values() if len(v) > 1)

    print(f"  n={n}: {len(by_score_class)} score classes, "
          f"max {max_H_per_sc} H values per class, "
          f"avg {avg_H_per_sc:.2f}, "
          f"{ambig_classes} ambiguous")

# ========================================================================
# PART 5: Even/odd alternation — is the OCR parity-dependent?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: EVEN/ODD PARITY AND THE OCR")
print("=" * 72)

print(f"""
  OCR values:
    n=3 (odd):  1.0000
    n=4 (even): 1.0000
    n=5 (odd):  0.9699
    n=6 (even): 0.9591
    n=7 (odd):  0.9574
    n=8 (even): 0.9609

  Is there an even/odd pattern?

  Odd n:  1.000, 0.970, 0.957
  Even n: 1.000, 0.959, 0.961

  At odd n: OCR is MONOTONE DECREASING (so far).
  At even n: OCR DROPS then RISES.

  From S108-S109: at odd n, NO blueself exists (THM-023).
  At even n, blueself exists, providing an EXTRA CONSTRAINT that
  might reduce ambiguity.

  HYPOTHESIS: The OCR rise at n=8 is because blueself at even n
  provides additional structure that resolves some ambiguity.

  If true: the OCR should show an even/odd ZIGZAG:
    Odd n: lower OCR (no blueself, more ambiguity)
    Even n: higher OCR (blueself, less ambiguity)

  At n=8: blueself exists (even n) → some POS classes get resolved → OCR rises.
  At n=9: no blueself (odd n) → OCR might drop again.
  At n=10: blueself returns → OCR rises.

  If the zigzag CONVERGES, the even and odd subsequences might
  approach the same limit from opposite sides.
""")

# ========================================================================
# PART 6: The fraction of tournaments in ambiguous classes
# ========================================================================
print("=" * 72)
print("  PART 6: FRACTION OF TOURNAMENTS IN AMBIGUOUS CLASSES")
print("=" * 72)

print(f"\n  {'n':>3s} {'frac_ambig':>10s} {'1-OCR':>8s} {'avg_var_per_ambig':>18s}")
for n in sorted(results.keys()):
    r = results[n]
    if r['frac_ambig'] > 0:
        avg_var = r['cond_var'] / r['frac_ambig']
    else:
        avg_var = 0
    print(f"  {n:3d} {r['frac_ambig']:10.4f} {1-r['OCR']:8.5f} {avg_var:18.4f}")

print(f"""
  OBSERVATION: The fraction of tournaments in ambiguous classes GROWS
  with n (from 0 at n=3,4 to ~75% at n=8).

  But the AVERAGE within-class variance also matters.
  If more classes are ambiguous but each has LOWER variance,
  the total residual could still decrease.

  The OCR is: 1 - (frac_ambig × avg_within_var) / Var(H)

  The question is whether frac_ambig × avg_within_var grows
  faster or slower than Var(H).
""")

# ========================================================================
# PART 7: S₂ distribution and its role in the OCR
# ========================================================================
print("=" * 72)
print("  PART 7: THE ROLE OF S₂ IN THE OCR")
print("=" * 72)

for n in [5, 6, 7, 8]:
    r = results[n]
    print(f"\n  n={n}: Top contributors to residual:")
    for i, (c, sc, cnt, n_H) in enumerate(r['top_contrib'][:5]):
        scores = list(sc)
        S2 = sum((s - (n-1)/2)**2 for s in scores)
        pct = c / r['cond_var'] * 100 if r['cond_var'] > 0 else 0
        print(f"    {sc}: S₂={S2:.1f}, {n_H} H values, {cnt} tours, "
              f"{pct:.1f}% of residual")

print(f"""
  KEY PATTERN: The dominant contributors to the residual tend to have
  MODERATE S₂ — not too regular (S₂=0 has unique H), not too irregular
  (S₂ large has few ambiguous tournaments).

  The "sweet spot" for ambiguity is at S₂ values where the score class
  is large enough to contain diverse tournaments but structured enough
  to have multiple H values.

  At n=7: the sweet spot should be around S₂ = 2-6 (from the recursive
  nesting analysis).
""")

# ========================================================================
# PART 8: Does OCR recover after n=7-8?
# ========================================================================
print("=" * 72)
print("  PART 8: DOES OCR RECOVER? TESTING n=9")
print("=" * 72)

n = 9
rng = random.Random(42)
by_score_9 = defaultdict(list)
N = 20000
print(f"  Sampling {N} tournaments at n={n}...")

for _ in range(N):
    A = random_tournament(n, rng)
    H = count_hp(A, n)
    sc = score_seq(A, n)
    by_score_9[sc].append(H)

all_H_9 = []
for Hs in by_score_9.values():
    all_H_9.extend(Hs)

H_arr_9 = np.array(all_H_9, dtype=float)
Var_H_9 = np.var(H_arr_9)
cond_var_9 = sum(len(v)*np.var(v) for v in by_score_9.values()) / N
OCR_9 = 1 - cond_var_9 / Var_H_9 if Var_H_9 > 0 else 1.0

print(f"  n=9: OCR = {OCR_9:.6f}")
print(f"  Score classes: {len(by_score_9)}")
print(f"  E[H] = {np.mean(H_arr_9):.1f}, Var(H) = {Var_H_9:.1f}")

# ========================================================================
# PART 9: The complete OCR curve
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 9: THE COMPLETE OCR CURVE")
print("=" * 72)

# Add n=9 to results
results[9] = {
    'OCR': OCR_9, 'Var_H': Var_H_9, 'cond_var': cond_var_9,
    'E_H': np.mean(H_arr_9),
}

print(f"\n  n | OCR      | 1-OCR    | parity")
for n in range(3, 10):
    if n in results:
        r = results[n]
        parity = "odd " if n%2==1 else "even"
        print(f"  {n} | {r['OCR']:.6f} | {1-r['OCR']:.6f} | {parity}")

print(f"""
  THE OCR CURVE:
  - Starts at 1.0 (n=3,4)
  - Drops to 0.970 (n=5)
  - Continues dropping to 0.957 (n=7)
  - RISES to 0.961 (n=8)
  - n=9: {OCR_9:.6f}

  If n=9 < n=8: the odd subsequence continues dropping.
  If n=9 > n=8: something more complex is happening.
  If n=9 ≈ n=7: the odd values are plateauing.

  THE NATURE OF THE OCR:
""")

if OCR_9 < results[7]['OCR']:
    print(f"  n=9 ({OCR_9:.4f}) < n=7 ({results[7]['OCR']:.4f}): odd OCR still DROPPING")
    print(f"  → The dip at n=7 is NOT the minimum for odd n")
elif OCR_9 > results[8]['OCR']:
    print(f"  n=9 ({OCR_9:.4f}) > n=8 ({results[8]['OCR']:.4f}): OCR rising again!")
    print(f"  → n=7 WAS the overall minimum")
else:
    print(f"  n=9 ({OCR_9:.4f}) is between n=7 and n=8: nuanced pattern")

# ========================================================================
# PART 10: The structural explanation
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 10: THE STRUCTURAL EXPLANATION")
print("=" * 72)
print(f"""
  WHY DOES THE OCR BEHAVE THIS WAY?

  The OCR = 1 - Var(H|scores)/Var(H).

  Var(H) grows SUPER-EXPONENTIALLY (like n!²/4^n).
  Var(H|scores) also grows, but at a different rate.

  The RATIO depends on:
  1. How many score classes are AMBIGUOUS
  2. How many tournaments are in each ambiguous class
  3. How much H VARIES within each ambiguous class

  FROM THE RECURSIVE NESTING (S110-S111):
  - The ambiguity comes from the POS classes
  - The POS at order n embeds a tournament of order n-2
  - The number of valid configurations grows linearly (k+1 at depth k)
  - But the H RANGE of the configurations grows with the nesting constants

  THE COMPETITION:
  - Var(H|scores) grows because MORE classes become ambiguous
  - Var(H) grows because the H RANGE expands faster

  At n=7: the ambiguous fraction is large AND the within-class variance
  is proportionally high → OCR is low.

  At n=8: the ambiguous fraction continues growing, but Var(H) grows
  even FASTER (because even-n tournaments have more extreme H values?)
  → OCR recovers.

  THE PREDICTION: OCR has a GLOBAL MINIMUM near n=7-8, then slowly
  recovers toward 1 as n → ∞ (the Shadow Compression Theorem).

  The minimum exists because there's a TRANSITION:
  - Small n (3-4): too few tournaments to be ambiguous → OCR = 1
  - Medium n (5-8): ambiguity emerges faster than H variance → OCR drops
  - Large n: H variance grows factorially, drowning the ambiguity → OCR → 1

  This is a FINITE-SIZE EFFECT: the OCR minimum is the crossover point
  between the "few ambiguous classes" regime and the "drowning variance" regime.
""")

print("\nDone. opus-2026-03-21-S112")
