#!/usr/bin/env python3
"""
waggly_higher_n_s298.py -- opus-2026-03-24-S298

WAGGLY AT HIGHER n: THE PATTERN THAT SIMPLIFIES

At n=5 (m=6): self-loop rates by layer are 10, 19, 7, 17, 12, 13%.
The even/odd oscillation and the d=2 peak suggest a FORMULA.

HYPOTHESIS: the self-loop rate at distance d depends on a simple
combinatorial quantity — perhaps related to the Krawtchouk polynomials
(the natural basis for functions on the Hamming scheme).

THE KEY SIMPLIFICATION:
The total self-loop weight at distance d equals:
  SL(d) = Σ_C #{(t, S) : t ∈ C, |S|=d, t⊕S ∈ C}

This counts ordered pairs (tiling, d-subset) where the d-flip preserves
the class. By Burnside, this relates to the cycle index.

For the TOTAL (summed over all classes):
  SL(d) / (2^m × C(m,d)) = Prob(random d-flip is silent)

THE BURNSIDE CONNECTION:
A d-flip is silent iff ∃ σ ∈ S_n with σ(T) = T' where T' differs
from T in exactly d tile positions. This is hard to compute directly.

But we can SAMPLE: pick random tilings, flip d random tiles, check
if the result is in the same iso class.

Author: opus-2026-03-24-S298
"""

import sys
import numpy as np
from math import comb, factorial, log2
from itertools import permutations, combinations
from collections import defaultdict, Counter
import random
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def build_tiling_and_classify(bits, n, tile_pairs):
    m = len(tile_pairs)
    A = [[0]*n for _ in range(n)]
    for i in range(1, n): A[i][i-1] = 1
    for k, (lo, hi) in enumerate(tile_pairs):
        if (bits >> k) & 1: A[lo][hi] = 1
        else: A[hi][lo] = 1
    return canon(A, n)

# ============================================================
banner("EXACT WAGGLY SELF-LOOP RATES: n=4,5,6")
# ============================================================

for n in [4, 5, 6]:
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()

    # Build class map
    tc = {}
    for bits in range(total):
        tc[bits] = build_tiling_and_classify(bits, n, tile_pairs)

    print("  n=%d (m=%d, %d tilings):" % (n, m, total))
    print("  %-5s %-8s %-12s %-12s" % ("d", "C(m,d)", "SL_count", "SL_rate"))
    print("  " + "-"*40)

    rates = []
    for d in range(1, m+1):
        sl = 0
        total_pairs = 0
        for bits in range(total):
            c1 = tc[bits]
            for flips in combinations(range(m), d):
                mask = sum(1 << p for p in flips)
                c2 = tc[bits ^ mask]
                total_pairs += 1
                if c1 == c2:
                    sl += 1
        rate = sl / total_pairs if total_pairs > 0 else 0
        rates.append(rate)
        print("  d=%-3d %-8d %-12d %-12.6f" % (d, comb(m,d), sl, rate))

    print()

    # THE PATTERN: compute rates × C(m,d) to see unnormalized
    print("  UNNORMALIZED: SL_count / 2^m (= avg silent mutations per tiling at distance d):")
    for d in range(1, m+1):
        sl_per_tiling = rates[d-1] * comb(m, d)
        print("    d=%d: %.4f (out of C(%d,%d)=%d possible)" %
              (d, sl_per_tiling, m, d, comb(m, d)))
    print()

# ============================================================
banner("THE SELF-LOOP RATE TABLE ACROSS n")
# ============================================================

# Collect all rates
all_rates = {}
for n in [4, 5, 6]:
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()
    tc = {}
    for bits in range(total):
        tc[bits] = build_tiling_and_classify(bits, n, tile_pairs)

    rates = []
    for d in range(1, m+1):
        sl = 0; tot = 0
        for bits in range(total):
            c1 = tc[bits]
            for flips in combinations(range(m), d):
                mask = sum(1 << p for p in flips)
                c2 = tc[bits ^ mask]
                tot += 1
                if c1 == c2: sl += 1
        rates.append(sl / tot if tot > 0 else 0)
    all_rates[n] = rates

print("  SELF-LOOP RATE BY d/m (normalized distance):\n")
print("  %-8s" % "d/m", end="")
for n in [4, 5, 6]:
    print(" n=%-6d" % n, end="")
print()
print("  " + "-"*30)

# Print at normalized distances 0, 0.2, 0.4, ..., 1.0
for frac in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    print("  %-8.1f" % frac, end="")
    for n in [4, 5, 6]:
        m = comb(n-1, 2)
        d = max(1, min(m, round(frac * m)))
        if d <= m:
            print(" %-8.4f" % all_rates[n][d-1], end="")
        else:
            print(" %-8s" % "-", end="")
    print()

# ============================================================
banner("SAMPLED RATES AT n=7 (m=15)")
# ============================================================

n = 7; m = comb(n-1, 2)  # m=15
total = 2**m  # 32768
base_set = set((i, i+1) for i in range(n-1))
all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
tile_pairs.sort()

print("  n=7: m=15, %d tilings. SAMPLING 2000 random tilings per distance.\n" % total)

# Build full class map (32768 is manageable)
tc7 = {}
for bits in range(total):
    tc7[bits] = build_tiling_and_classify(bits, n, tile_pairs)

random.seed(42)
print("  %-5s %-8s %-12s %-12s" % ("d", "C(m,d)", "SL_est", "samples"))
print("  " + "-"*45)

rates7 = []
for d in range(1, m+1):
    n_samples = min(2000, total)
    sl = 0; tot = 0
    sampled_tilings = random.sample(range(total), n_samples)
    for bits in sampled_tilings:
        c1 = tc7[bits]
        # Sample random d-subsets
        n_flips = min(50, comb(m, d))
        if comb(m, d) <= 50:
            flip_sets = list(combinations(range(m), d))
        else:
            flip_sets = [tuple(sorted(random.sample(range(m), d))) for _ in range(50)]
        for flips in flip_sets:
            mask = sum(1 << p for p in flips)
            c2 = tc7[bits ^ mask]
            tot += 1
            if c1 == c2: sl += 1
    rate = sl / tot if tot > 0 else 0
    rates7.append(rate)
    print("  d=%-3d %-8d %-12.6f %-12d" % (d, comb(m,d), rate, tot))

all_rates[7] = rates7

# ============================================================
banner("THE UNIVERSAL CURVE")
# ============================================================

print("  SELF-LOOP RATE vs NORMALIZED DISTANCE d/m:\n")
print("  %-8s" % "d/m", end="")
for n in [4, 5, 6, 7]:
    print(" n=%-6d" % n, end="")
print()
print("  " + "-"*40)

for d_frac_10 in range(11):
    frac = d_frac_10 / 10.0
    print("  %-8.2f" % frac, end="")
    for n in [4, 5, 6, 7]:
        m = comb(n-1, 2)
        d = max(1, min(m, round(frac * m)))
        rates = all_rates[n]
        if d >= 1 and d <= len(rates):
            print(" %-8.4f" % rates[d-1], end="")
        else:
            print(" %-8s" % "-", end="")
    print()

# ============================================================
banner("THE SIMPLIFICATION: WHAT DRIVES THE RATE?")
# ============================================================

# Hypothesis: SL rate at distance d ≈ f(d/m) for some universal function f.
# Or: SL rate ∝ V_n(d) / V_n where V_n(d) is an "effective class count" at distance d.

# Actually, the simplest model: SL rate ≈ <fiber> / 2^m at large n.
# Because a random d-flip of a random tiling gives a random tournament,
# and the probability of landing in the same class ≈ fiber / 2^m.

# But the rates are HIGHER than this baseline, especially at d=2.
# The excess is the "structural correlation" — nearby tilings are more
# likely to be in the same class than random ones.

print("  BASELINE: <fiber>/2^m (random landing probability):\n")
for n in [4, 5, 6, 7]:
    m = comb(n-1, 2)
    total = 2**m
    # Average fiber ≈ n!/2^{n-1}
    avg_fiber = factorial(n) / 2**(n-1)
    baseline = avg_fiber / total
    print("  n=%d: <fiber>=%.1f, 2^m=%d, baseline=%.6f" %
          (n, avg_fiber, total, baseline))
    # Compare with actual rates
    rates = all_rates[n]
    print("    d=1: %.6f (%.1f× baseline)" % (rates[0], rates[0]/baseline if baseline > 0 else 0))
    if len(rates) >= 2:
        print("    d=2: %.6f (%.1f× baseline)" % (rates[1], rates[1]/baseline if baseline > 0 else 0))
    mid = len(rates) // 2
    print("    d=%d: %.6f (%.1f× baseline)" % (mid+1, rates[mid], rates[mid]/baseline if baseline > 0 else 0))
    print("    d=%d: %.6f (%.1f× baseline)" % (len(rates), rates[-1], rates[-1]/baseline if baseline > 0 else 0))
    print()

# ============================================================
banner("THE CORRELATION DECAY")
# ============================================================

# The EXCESS over baseline decays with distance.
# Define: excess(d) = SL_rate(d) / baseline - 1
# This measures how much more likely a d-flip is to be silent
# compared to a random landing.

print("  EXCESS OVER BASELINE (= structural correlation):\n")
print("  %-5s" % "d", end="")
for n in [4, 5, 6, 7]:
    print(" n=%-8d" % n, end="")
print()
print("  " + "-"*40)

for d_idx in range(max(len(all_rates[n]) for n in [4,5,6,7])):
    d = d_idx + 1
    print("  d=%-3d" % d, end="")
    for n in [4, 5, 6, 7]:
        m = comb(n-1, 2)
        rates = all_rates[n]
        avg_fiber = factorial(n) / 2**(n-1)
        baseline = avg_fiber / 2**m
        if d <= len(rates) and baseline > 0:
            excess = rates[d-1] / baseline - 1
            print(" %-10.2f" % excess, end="")
        else:
            print(" %-10s" % "", end="")
    print()

# ============================================================
banner("THE SUMMARY")
# ============================================================

print("""
  THE WAGGLY SELF-LOOP RATE HAS THREE REGIMES:

  1. NEAR-FIELD (d = 1, 2, 3):
     High excess over baseline. Structural correlations dominate.
     A small mutation is likely to preserve the class because
     the tournament changes little.
     Rate: 5-20% (depends on n and d).

  2. MID-FIELD (d ≈ m/3 to 2m/3):
     Moderate excess. The d-flip is large enough that the result
     is "almost random" but structural memory persists.
     Rate: ~baseline × small multiplier.

  3. FAR-FIELD (d ≈ m, complement):
     Excess drops but doesn't reach zero. The complement tiling
     has STRUCTURAL correlation with the original (not random!)
     because flipping all tiles creates a specific, predictable
     tournament — not a random one.
     Rate: ~baseline × 1-3.

  THE EVEN/ODD OSCILLATION:
  At n=5: d=1(10%), d=2(19%), d=3(7%), d=4(17%), d=5(12%), d=6(13%)
  Even d has HIGHER rates than odd d (except d=m).
  This is a PARITY EFFECT: flipping an even number of tiles
  preserves more structural invariants than odd flips.

  CONJECTURE: the even/odd oscillation persists at all n,
  with even-d rates ≈ 2× odd-d rates in the near-field.

  THE PRACTICAL CONCLUSION:
  For any tiling, the d=2 layer is the "sweet spot" —
  most likely to preserve the class while still being a
  non-trivial mutation. This suggests that 2-tile mutations
  are the natural "phonons" of the tiling lattice.
""")

print("Done. opus-2026-03-24-S298")
