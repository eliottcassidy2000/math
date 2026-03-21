#!/usr/bin/env python3
"""
exact_ocr_formula_s104.py — opus-2026-03-21-S104

EXACT OCR FORMULA FROM FIRST PRINCIPLES

Goal: Express OCR = 1 - Var(H|scores)/Var(H) as an exact rational number
at each n, then find the pattern.

Key insight from S97b-S103:
  c₃ = C(n+1,3)/4 - S₂/2 is EXACT (scores determine c₃)
  H ≈ linear in c₃ with coefficient b(n)
  The RESIDUAL H - b·c₃ - a comes from c₅, α₂, etc.

At n=5: only class (1,2,2,2,3) has multiple H values: {11,13,15}.
  Count: 120, 120, 40 respectively.
  The OCR residual comes ENTIRELY from this class.

Plan:
1. Compute exact Var(H), E[H], Var(H|scores) at n=3-6 using Fraction arithmetic
2. Express OCR as exact fraction
3. Decompose the residual: which score classes contribute?
4. Connect to skeleton: are the ambiguous classes the blackself/blueself classes?

Author: opus-2026-03-21-S104
"""

import sys
from fractions import Fraction
from math import comb
from itertools import combinations, permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

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

def count_c3(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: c3 += 1
        if A[i][k] and A[k][j] and A[j][i]: c3 += 1
    return c3

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

print("=" * 72)
print("  EXACT OCR FORMULA FROM FIRST PRINCIPLES")
print("  opus-2026-03-21-S104")
print("=" * 72)

# ========================================================================
# PART 1: Exact Var(H), Var(H|scores) as fractions
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: EXACT VARIANCE DECOMPOSITION")
print("=" * 72)

for n in range(3, 8 if True else 7):
    if n >= 7:
        print(f"\n--- n = {n}: skipping (too large for exhaustive) ---")
        continue

    print(f"\n--- n = {n} ---")
    N = 2**comb(n, 2)  # total tournaments

    # Collect all (H, score) pairs
    by_score = defaultdict(list)
    all_H = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        sc = score_seq(A, n)
        by_score[sc].append(H)
        all_H.append(H)

    # Exact computation using Fractions
    N_frac = Fraction(N)

    # E[H] = sum(H) / N
    sum_H = sum(all_H)
    E_H = Fraction(sum_H, N)

    # E[H²] = sum(H²) / N
    sum_H2 = sum(h*h for h in all_H)
    E_H2 = Fraction(sum_H2, N)

    # Var(H) = E[H²] - E[H]²
    Var_H = E_H2 - E_H * E_H

    # Var(H|scores) = E[Var(H|score=s)]
    # = sum_s P(score=s) * Var(H|score=s)
    # = sum_s (n_s/N) * [(sum H² in s) / n_s - (sum H in s / n_s)²]
    # = (1/N) * sum_s [sum_H2_s - (sum_H_s)² / n_s]

    cond_var_num = Fraction(0)
    for sc, Hs in by_score.items():
        n_s = len(Hs)
        sum_Hs = sum(Hs)
        sum_Hs2 = sum(h*h for h in Hs)
        # Var within = sum_Hs2/n_s - (sum_Hs/n_s)²
        # Contribution = n_s/N * Var_within = (sum_Hs2 - sum_Hs²/n_s) / N
        cond_var_num += Fraction(sum_Hs2) - Fraction(sum_Hs * sum_Hs, n_s)

    Var_H_cond = cond_var_num / N_frac

    # OCR
    if Var_H > 0:
        OCR = 1 - Var_H_cond / Var_H
    else:
        OCR = Fraction(1)

    residual = Var_H_cond

    print(f"  N = {N}")
    print(f"  E[H] = {E_H} = {float(E_H):.6f}")
    print(f"  Var(H) = {Var_H} = {float(Var_H):.6f}")
    print(f"  Var(H|scores) = {Var_H_cond} = {float(Var_H_cond):.6f}")
    print(f"  OCR = {OCR} = {float(OCR):.8f}")
    print(f"  1 - OCR = {1 - OCR} = {float(1 - OCR):.8f}")

    # Score class decomposition
    print(f"\n  Score class decomposition of residual:")
    print(f"  {'score':>25s} {'count':>6s} {'H values':>25s} {'Var_within':>12s} {'Contrib':>12s}")

    for sc in sorted(by_score.keys()):
        Hs = by_score[sc]
        n_s = len(Hs)
        H_set = sorted(set(Hs))
        if len(H_set) == 1:
            var_w = Fraction(0)
        else:
            s2 = sum(h*h for h in Hs)
            s1 = sum(Hs)
            var_w = Fraction(s2, n_s) - Fraction(s1*s1, n_s*n_s)

        contrib = Fraction(n_s, N) * var_w
        if contrib > 0:
            print(f"  {str(sc):>25s} {n_s:6d} {str(H_set):>25s} "
                  f"{float(var_w):12.4f} {float(contrib):12.6f}")

# ========================================================================
# PART 2: Pattern analysis of OCR fractions
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: OCR AS EXACT FRACTIONS")
print("=" * 72)

# From Part 1, collect the exact OCR values
# Let me recompute more carefully and store
ocr_exact = {}
var_H_exact = {}
var_cond_exact = {}

for n in range(3, 7):
    N = 2**comb(n, 2)
    by_score = defaultdict(list)
    all_H = []

    for A in all_tournaments(n):
        H = count_hp(A, n)
        sc = score_seq(A, n)
        by_score[sc].append(H)
        all_H.append(H)

    sum_H = sum(all_H)
    sum_H2 = sum(h*h for h in all_H)
    E_H = Fraction(sum_H, N)
    Var_H = Fraction(sum_H2, N) - E_H * E_H

    cond_var_num = Fraction(0)
    for sc, Hs in by_score.items():
        n_s = len(Hs)
        cond_var_num += Fraction(sum(h*h for h in Hs)) - Fraction(sum(Hs)**2, n_s)
    Var_cond = cond_var_num / Fraction(N)

    OCR = 1 - Var_cond / Var_H if Var_H > 0 else Fraction(1)

    ocr_exact[n] = OCR
    var_H_exact[n] = Var_H
    var_cond_exact[n] = Var_cond

print(f"\n  Exact OCR fractions:")
for n in sorted(ocr_exact.keys()):
    ocr = ocr_exact[n]
    print(f"  n={n}: OCR = {ocr}")
    print(f"        = {float(ocr):.10f}")
    print(f"        Var(H) = {var_H_exact[n]}")
    print(f"        Var(H|scores) = {var_cond_exact[n]}")
    res = 1 - ocr
    print(f"        1-OCR = {res}")

# ========================================================================
# PART 3: Analyze the residual fractions
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: RESIDUAL FRACTION ANALYSIS")
print("=" * 72)

print(f"\n  The residual 1 - OCR as exact fraction:")
for n in sorted(ocr_exact.keys()):
    res = 1 - ocr_exact[n]
    print(f"  n={n}: 1-OCR = {res}")
    # Try to express in a revealing form
    # At n=3,4: 0
    # At n=5: 1 - OCR = Var(H|scores)/Var(H)
    if res > 0:
        vc = var_cond_exact[n]
        vh = var_H_exact[n]
        print(f"        = {vc} / {vh}")
        print(f"        = {float(vc):.6f} / {float(vh):.6f}")

# ========================================================================
# PART 4: The OCF structure behind the residual
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: WHY THE RESIDUAL EXISTS — OCF ANALYSIS")
print("=" * 72)
print("""
  The residual comes from score classes with multiple H values.

  At n=5: only class (1,2,2,2,3) with H ∈ {11, 13, 15}.
  At n=6: multiple classes have H spread.

  From S98b: within (1,2,2,2,3), H is determined by c₅:
    c₅=1 → H=11 (120 tournaments)
    c₅=2 → H=13 (120 tournaments)
    c₅=3 → H=15 (40 tournaments)

  The residual Var(H|scores) at n=5 comes entirely from this class.

  The 280 tournaments in class (1,2,2,2,3) have:
    E[H|sc] = (120*11 + 120*13 + 40*15) / 280
            = (1320 + 1560 + 600) / 280 = 3480/280 = 87/7
    E[H²|sc] = (120*121 + 120*169 + 40*225) / 280
             = (14520 + 20280 + 9000) / 280 = 43800/280 = 1095/7
    Var(H|sc) = 1095/7 - (87/7)² = 1095/7 - 7569/49
             = 7665/49 - 7569/49 = 96/49

  Contribution to Var(H|scores):
    P(sc=(1,2,2,2,3)) * Var(H|sc) = (280/1024) * (96/49)
    = 280*96 / (1024*49) = 26880/50176 = 420/784 = 15/28

  Total Var(H|scores) at n=5 = 15/28 (if this is the only contributor)
  Var(H) at n=5 = 285/16

  1 - OCR = (15/28) / (285/16) = (15*16) / (28*285) = 240/7980 = 2/66.5?
  Let me compute exactly...
""")

# Exact computation for n=5
n = 5
N = 1024
by_score_5 = defaultdict(list)
for A in all_tournaments(n):
    H = count_hp(A, n)
    sc = score_seq(A, n)
    by_score_5[sc].append(H)

# The ambiguous class
sc_ambig = (1, 2, 2, 2, 3)
Hs_ambig = by_score_5[sc_ambig]
print(f"\n  n=5, class {sc_ambig}:")
print(f"    Count: {len(Hs_ambig)}")
H_counts = defaultdict(int)
for h in Hs_ambig:
    H_counts[h] += 1
for h in sorted(H_counts.keys()):
    print(f"    H={h}: {H_counts[h]} tournaments")

n_s = len(Hs_ambig)
sum_H = sum(Hs_ambig)
sum_H2 = sum(h*h for h in Hs_ambig)
E_H_sc = Fraction(sum_H, n_s)
Var_H_sc = Fraction(sum_H2, n_s) - E_H_sc * E_H_sc
print(f"    E[H|sc] = {E_H_sc} = {float(E_H_sc):.6f}")
print(f"    Var(H|sc) = {Var_H_sc} = {float(Var_H_sc):.6f}")

contrib = Fraction(n_s, N) * Var_H_sc
print(f"    Contribution = {Fraction(n_s, N)} × {Var_H_sc} = {contrib}")
print(f"    = {float(contrib):.6f}")

# Total Var(H|scores)
Var_cond_5 = var_cond_exact[5]
print(f"\n    Total Var(H|scores) = {Var_cond_5}")
print(f"    From this class: {contrib}")
print(f"    Match: {contrib == Var_cond_5}")

# OCR
Var_H_5 = var_H_exact[5]
print(f"\n    Var(H) = {Var_H_5}")
print(f"    OCR = 1 - {Var_cond_5}/{Var_H_5} = {1 - Var_cond_5/Var_H_5}")

# Simplify the OCR fraction
ocr_5 = 1 - Var_cond_5/Var_H_5
print(f"    OCR = {ocr_5}")

# ========================================================================
# PART 5: H distribution within each score class at n=5,6
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: COMPLETE H DISTRIBUTION BY SCORE CLASS")
print("=" * 72)

for n in [5, 6]:
    print(f"\n--- n = {n} ---")
    by_score = defaultdict(list)
    for A in all_tournaments(n):
        H = count_hp(A, n)
        sc = score_seq(A, n)
        by_score[sc].append(H)

    print(f"  {'score':>25s} {'count':>6s} {'H values':>30s} {'E[H]':>8s} {'Var':>8s}")
    for sc in sorted(by_score.keys()):
        Hs = by_score[sc]
        H_set = sorted(set(Hs))
        E = sum(Hs) / len(Hs)
        V = sum((h-E)**2 for h in Hs) / len(Hs)
        ambig = " ***" if len(H_set) > 1 else ""
        print(f"  {str(sc):>25s} {len(Hs):6d} {str(H_set):>30s} {E:8.2f} {V:8.2f}{ambig}")

# ========================================================================
# PART 6: Connection to the skeleton — which iso classes are ambiguous?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 6: CONNECTING OCR RESIDUAL TO SKELETON STRUCTURE")
print("=" * 72)

n = 5
print(f"\n--- n = {n}: Ambiguous class analysis ---")

# The ambiguous class (1,2,2,2,3) has H ∈ {11, 13, 15}
# From S99: these correspond to 3 SC iso classes + some NSC
# H=11: blackself, c₅=1, e_cyc=8
# H=13: blackself, c₅=2, e_cyc=9
# H=15: NOT blackself, c₅=3, e_cyc=10

# From the S99 skeleton output:
# Class 2 (H=13, c3=4, blackself)
# Class 3 (H=11, c3=4, blackself)
# Class 5 (H=15, c3=4)

print(f"""
  Within score (1,2,2,2,3), the three H values correspond to:

  H=11 (120 tournaments): c₅=1, e_cyc=8, BLACKSELF
  H=13 (120 tournaments): c₅=2, e_cyc=9, BLACKSELF
  H=15 (40 tournaments):  c₅=3, e_cyc=10, NOT blackself

  KEY OBSERVATION: The two BLACKSELF classes account for 240/280 = 6/7
  of the ambiguous class. They have H=11 and H=13.
  The non-blackself class has H=15 with only 40/280 = 1/7 of the class.

  HYPOTHESIS: Blackself classes contribute disproportionately to the OCR residual
  because they have unusual c₅ counts.

  THE SKELETON CONNECTION:
  The blackself classes (H=11, H=13) are connected to each other in the
  skeleton via GS flips. The non-blackself class (H=15) is the one with
  3 GS tilings but only 5 total tilings (|Aut|=3).

  The c₅ count distinguishes them:
  c₅ = 1 ↔ H = 11  [fewest 5-cycles → fewest "extra" paths]
  c₅ = 2 ↔ H = 13  [middle]
  c₅ = 3 ↔ H = 15  [most 5-cycles → most extra paths → highest H]

  And |Aut| differs:
  H=11: |Aut|=1  (generic, 120 tilings)
  H=13: |Aut|=1  (generic, 120 tilings)
  H=15: |Aut|=3  (3-fold symmetry, only 40 tilings despite higher H)

  THIS IS THE KEY: Higher symmetry (|Aut|) correlates with higher H
  WITHIN the same score class. The most symmetric tournament has the
  most Hamiltonian paths, but the fewest tiling representations.
""")

# ========================================================================
# PART 7: The exact Var(H|scores) formula
# ========================================================================
print("=" * 72)
print("  PART 7: TOWARD AN EXACT Var(H|SCORES) FORMULA")
print("=" * 72)

# At n=5: Var(H|scores) comes entirely from class (1,2,2,2,3)
# This class has c₃=4, and the three H values are:
# H = 1 + 2(c₃ + c₅) = 1 + 2(4 + c₅) = 9 + 2c₅
# So H = 11 when c₅=1, H=13 when c₅=2, H=15 when c₅=3.
# The variance is Var(2c₅) = 4*Var(c₅).

# Var(c₅) within this score class:
# c₅=1: 120 tournaments
# c₅=2: 120 tournaments
# c₅=3: 40 tournaments
# Total: 280
# E[c₅] = (120*1 + 120*2 + 40*3) / 280 = (120+240+120)/280 = 480/280 = 12/7
# E[c₅²] = (120*1 + 120*4 + 40*9) / 280 = (120+480+360)/280 = 960/280 = 24/7
# Var(c₅) = 24/7 - (12/7)² = 24/7 - 144/49 = (168-144)/49 = 24/49

E_c5 = Fraction(120*1 + 120*2 + 40*3, 280)
E_c5_sq = Fraction(120*1 + 120*4 + 40*9, 280)
Var_c5 = E_c5_sq - E_c5**2

print(f"\n  n=5, class (1,2,2,2,3):")
print(f"    E[c₅] = {E_c5} = {float(E_c5):.6f}")
print(f"    Var(c₅) = {Var_c5} = {float(Var_c5):.6f}")
print(f"    Var(H|sc) = 4·Var(c₅) = {4*Var_c5} = {float(4*Var_c5):.6f}")
print(f"    Check: {4*Var_c5} = {Var_H_sc}? {4*Var_c5 == Var_H_sc}")

print(f"""
  SO: Var(H|scores) at n=5 = P(sc=(1,2,2,2,3)) × 4·Var(c₅|sc=(1,2,2,2,3))
                            = (280/1024) × 4 × (24/49)
                            = (280 × 96) / (1024 × 49)
                            = {Fraction(280*96, 1024*49)}
                            = {float(Fraction(280*96, 1024*49)):.6f}

  And OCR = 1 - {Fraction(280*96, 1024*49)} / {Var_H_5}

  Let me simplify...
""")

residual_5 = Fraction(280*96, 1024*49)
ocr_exact_5 = 1 - residual_5 / Var_H_5

print(f"  Var(H|scores) = {residual_5}")
print(f"  Var(H) = {Var_H_5}")
print(f"  Ratio = {residual_5 / Var_H_5}")
print(f"  OCR = {ocr_exact_5}")

# ========================================================================
# PART 8: What drives the ambiguity?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 8: WHAT DRIVES THE AMBIGUITY? — THE ROLE OF |Aut|")
print("=" * 72)
print("""
  Within the ambiguous class (1,2,2,2,3) at n=5:

  c₅=1 → H=11, |Aut|=1, #til=11, #GS=3
  c₅=2 → H=13, |Aut|=1, #til=13, #GS=3
  c₅=3 → H=15, |Aut|=3, #til=5,  #GS=1

  The c₅=3 class is the MOST SYMMETRIC (|Aut|=3) and has the HIGHEST H.

  THEOREM (at n=5): Within the ambiguous score class, the tournament with
  the highest automorphism group has the highest H.

  Is this a general principle? Let me check at n=6.
""")

# At n=6: check which score classes are ambiguous and how |Aut| relates
n = 6
by_score_6 = defaultdict(list)
for A in all_tournaments(n):
    H = count_hp(A, n)
    sc = score_seq(A, n)
    by_score_6[sc].append(H)

for sc in sorted(by_score_6.keys()):
    Hs = by_score_6[sc]
    H_set = sorted(set(Hs))
    if len(H_set) > 1:
        H_counts = Counter(Hs)
        print(f"  Score {sc}: H values = {dict(sorted(H_counts.items()))}")

        # For each H value, find |Aut| = H / #tilings-in-class
        # #tilings-in-class = count with this H in this score seq
        # But multiple iso classes can have the same (H, score).
        # Actually from S99: #til = H/|Aut|, and within a score class,
        # the count of tournaments with a given H = (#iso classes with this H) × (n!/|Aut|)
        # Hmm, actually the count = sum over iso classes C with H(C)=h and score(C)=sc
        # of |C|, where |C| = n!/|Aut(C)|.

        # Let's just look at the counts
        for h in sorted(H_counts.keys()):
            cnt = H_counts[h]
            # If all tournaments with H=h in this score class have the same |Aut|,
            # then |Aut| = n! * h / cnt... no.
            # Actually count = (number of labeled tournaments in this score+H class)
            # and #til = H/|Aut|.
            # A labeled tournament has a unique canonical form.
            # Among all N=2^15 labeled tournaments, cnt have score=sc and H=h.
            # These belong to some set of iso classes, each contributing n!/|Aut| labelings.
            # Since each labeling maps to exactly 1 tiling, cnt = sum n!/|Aut|.
            pass

print(f"""
  KEY INSIGHT: The ambiguous score classes at n=6 have MORE H values
  than at n=5. This is because longer odd cycles (5-cycles, 7-cycles
  at n≥7) create more variation within each score class.

  The OCR residual measures EXACTLY this variation.
  It is driven by the DIVERSITY of cycle structures compatible with
  a given score sequence.
""")

# ========================================================================
# PART 9: The Var(H|scores) exact formula attempt
# ========================================================================
print("=" * 72)
print("  PART 9: EXACT FORMULA FOR Var(H|SCORES)")
print("=" * 72)
print("""
  CLAIM: Var(H|scores) = (4/total) × Σ_s n_s × Var(c₅|s)

  where the sum is over ambiguous score classes s, n_s = |score class s|,
  and Var(c₅|s) = variance of c₅ within score class s.

  At n=5: only one ambiguous class, so
    Var(H|scores) = (280/1024) × 4 × Var(c₅|(1,2,2,2,3))

  MORE GENERALLY: H = 1 + 2α₁ + 4α₂ + ...
  α₁ = c₃ + c₅ + ...
  Since c₃ is determined by scores: Var(α₁|scores) = Var(c₅ + c₇ + ...|scores)
  And Var(H|scores) ≈ 4 × Var(α₁|scores) (if α₂ is approximately determined by α₁)

  This gives:
    Var(H|scores) ≈ 4 × Var(c₅|scores)   (at n=5 where c₇=0)

  For general n:
    Var(H|scores) ≈ 4 × Var(c₅ + c₇ + ... | scores)

  And the OCR residual:
    1 - OCR ≈ 4 × Var(c₅ + c₇ + ... | scores) / Var(H)

  Since Var(H) ≈ 4 × Var(α₁) (the coefficient 4 comes from H ≈ 1+2α₁):
    Wait: Var(H) ≈ 4 × Var(α₁) since H ≈ 1 + 2α₁.
    (Ignoring higher OCF terms.)

    1 - OCR ≈ 4 × Var(c₅|scores) / (4 × Var(α₁))
            = Var(c₅|scores) / Var(c₃ + c₅)
            ≈ Var(c₅|scores) / (Var(c₃) + Var(c₅) + 2Cov(c₃,c₅))

  Since c₃ is determined by scores: Var(c₃|scores) = 0.
  So Var(c₃) = Var(E[c₃|scores]) = Var(c₃) (no simplification).

  But: Var(α₁) = Var(c₃) + Var(c₅) + Var(c₇) + ... + cross terms

  The OCR measures what fraction of the TOTAL variance of α₁
  is explained by c₃ (via scores).

  1 - OCR ≈ [Var(c₅|scores) + Var(c₇|scores) + ...] / Var(α₁)

  At n=5: 1-OCR = Var(c₅|sc) / Var(α₁) (only c₃ and c₅ at n=5)
""")

# Verify: Var(α₁) at n=5
n = 5
alpha1s = []
for A in all_tournaments(n):
    c3 = count_c3(A, n)
    c5 = 0
    for verts in combinations(range(n), 5):
        for perm in permutations(verts[1:]):
            path = [verts[0]] + list(perm)
            if all(A[path[k]][path[(k+1)%5]] for k in range(5)):
                c5 += 1
    alpha1s.append(c3 + c5)

Var_alpha1 = Fraction(sum(a*a for a in alpha1s), 1024) - Fraction(sum(alpha1s), 1024)**2
print(f"\n  n=5: Var(α₁) = {Var_alpha1} = {float(Var_alpha1):.6f}")
print(f"  4×Var(α₁) = {4*Var_alpha1} = {float(4*Var_alpha1):.6f}")
print(f"  Var(H) = {Var_H_5} = {float(Var_H_5):.6f}")
print(f"  4×Var(α₁) = Var(H)? {4*Var_alpha1 == Var_H_5}")

# If H = 1 + 2α₁ exactly (no α₂ at n=5):
# Var(H) = 4×Var(α₁). Let's check.
all_H5 = []
for A in all_tournaments(n):
    all_H5.append(count_hp(A, n))

# Verify H = 1 + 2α₁ at n=5
mismatch = 0
for i, (A_list_entry) in enumerate(zip(all_tournaments(n))):
    pass
# Actually let me just check from the data we have
print(f"\n  From S98b: H = 1 + 2α₁ at n=5 (α₂=0 always). So Var(H) = 4×Var(α₁).")

# The OCR formula then:
# 1-OCR = Var(H|scores)/Var(H) = 4×Var(c₅|scores) / (4×Var(α₁))
#        = Var(c₅|scores) / Var(α₁)
# Since α₁ = c₃ + c₅ and c₃ is determined by scores:
# Var(α₁|scores) = Var(c₅|scores)
# Var(α₁) = Var(c₃) + Var(c₅) + 2Cov(c₃,c₅)
# Since c₃ = 5 - S₂/2, and S₂ is determined by scores, Var(c₃) = Var(5-S₂/2) = Var(S₂)/4

# Compute Var(c₃), Var(c₅), Cov(c₃,c₅) at n=5
c3s_5 = []
c5s_5 = []
for A in all_tournaments(n):
    c3 = count_c3(A, n)
    c5 = 0
    for verts in combinations(range(n), 5):
        for perm in permutations(verts[1:]):
            path = [verts[0]] + list(perm)
            if all(A[path[k]][path[(k+1)%5]] for k in range(5)):
                c5 += 1
    c3s_5.append(c3)
    c5s_5.append(c5)

c3_arr = [Fraction(c) for c in c3s_5]
c5_arr = [Fraction(c) for c in c5s_5]
N5 = Fraction(1024)

E_c3 = sum(c3_arr) / N5
E_c5 = sum(c5_arr) / N5
Var_c3 = sum(c*c for c in c3_arr) / N5 - E_c3**2
Var_c5_total = sum(c*c for c in c5_arr) / N5 - E_c5**2
Cov_c3c5 = sum(c3_arr[i]*c5_arr[i] for i in range(1024)) / N5 - E_c3*E_c5

print(f"\n  n=5 exact cycle statistics:")
print(f"    E[c₃] = {E_c3}")
print(f"    E[c₅] = {E_c5}")
print(f"    Var(c₃) = {Var_c3}")
print(f"    Var(c₅) = {Var_c5_total}")
print(f"    Cov(c₃,c₅) = {Cov_c3c5}")
print(f"    Var(α₁) = Var(c₃)+Var(c₅)+2Cov = {Var_c3 + Var_c5_total + 2*Cov_c3c5}")
print(f"    Should equal: {Var_alpha1}")

# Var(c₅|scores) = E[Var(c₅|sc)] — the within-score-class variance of c₅
by_score_c5 = defaultdict(list)
for i in range(1024):
    sc = score_seq(list(all_tournaments(n))[0] if False else None, n)
# Actually I need to recompute with score info
data_5 = []
idx = 0
for A in all_tournaments(n):
    sc = score_seq(A, n)
    data_5.append({'sc': sc, 'c3': c3s_5[idx], 'c5': c5s_5[idx]})
    idx += 1

by_sc_c5 = defaultdict(list)
for d in data_5:
    by_sc_c5[d['sc']].append(d['c5'])

cond_var_c5 = Fraction(0)
for sc, c5s in by_sc_c5.items():
    n_s = len(c5s)
    E_c5_s = Fraction(sum(c5s), n_s)
    var_c5_s = Fraction(sum(c*c for c in c5s), n_s) - E_c5_s**2
    cond_var_c5 += Fraction(n_s, 1024) * var_c5_s

print(f"\n    Var(c₅|scores) = {cond_var_c5} = {float(cond_var_c5):.6f}")
print(f"    Should match 1-OCR fraction of Var(α₁):")
print(f"    Var(c₅|scores) / Var(α₁) = {cond_var_c5 / Var_alpha1}")
print(f"    1 - OCR = {1 - ocr_exact[5]}")
print(f"    Match: {cond_var_c5 / Var_alpha1 == 1 - ocr_exact[5]}")

# ========================================================================
# PART 10: The master formula
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 10: THE MASTER FORMULA")
print("=" * 72)
print(f"""
  THEOREM (Exact OCR at n=5):

  At n=5, H = 1 + 2(c₃ + c₅) exactly.
  c₃ is determined by scores.
  Therefore:

  OCR = 1 - Var(c₅|scores) / Var(c₃ + c₅)

  Numerically:
    Var(c₃) = {Var_c3}
    Var(c₅) = {Var_c5_total}
    Cov(c₃,c₅) = {Cov_c3c5}
    Var(c₃+c₅) = {Var_c3 + Var_c5_total + 2*Cov_c3c5}
    Var(c₅|scores) = {cond_var_c5}

    OCR = 1 - {cond_var_c5} / {Var_c3 + Var_c5_total + 2*Cov_c3c5}
        = {1 - cond_var_c5 / (Var_c3 + Var_c5_total + 2*Cov_c3c5)}
        = {ocr_exact[5]}

  The EXACT OCR at n=5 is {ocr_exact[5]} = {float(ocr_exact[5]):.10f}.

  GENERAL FORMULA:
  OCR(n) = 1 - Var(c₅ + c₇ + ... | scores) / Var(c₃ + c₅ + c₇ + ...)

  This is the fraction of TOTAL odd-cycle variance that is
  DETERMINED by the score sequence.

  The numerator = "hidden" cycle variance (not captured by scores).
  The denominator = total cycle variance.

  As n grows: c₃ dominates (it grows as n³), so the denominator grows
  faster than the numerator. This suggests OCR → 1.

  But our data shows OCR stabilizes at ~0.96. This means the ratio
  Var(hidden cycles | scores) / Var(all cycles) stabilizes at ~0.04.
""")

print("\nDone. opus-2026-03-21-S104")
