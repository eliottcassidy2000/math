#!/usr/bin/env python3
"""
n8_sampling_s114.py — Sampling-based OCR estimation at n=8
opus-2026-03-21-S114

n=8: 2^28 = 268,435,456 tournaments. Full enumeration ~19 hours.
Instead: SAMPLE ~500K random tournaments, compute H and scores,
estimate OCR from the sample.

Also compute:
  - E[H] exactly (= n!/2^{n-1} = 40320/128 = 315)
  - E[H] from sample (should match)
  - Var(H) from sample
  - Score class statistics from sample
  - Which score classes are ambiguous
  - Top contributors to the OCR residual
"""

import numpy as np
from math import comb, factorial
from collections import defaultdict, Counter
from fractions import Fraction
import time
import random

print("=" * 76)
print("  SAMPLING-BASED OCR AT n=8")
print("  opus-2026-03-21-S114")
print("=" * 76)
print()

n = 8
m = comb(n, 2)  # 28
N_total = 2**m  # 268,435,456

print(f"  n = {n}, C(n,2) = {m}, total tournaments = {N_total:,}")
print(f"  E[H] = n!/2^{{n-1}} = {factorial(n)}/{2**(n-1)} = {factorial(n) // 2**(n-1)}")
print()

# =============================================================================
# FAST H COMPUTATION
# =============================================================================

def H_dp_fast(adj_flat, n):
    """Bitmask DP for Hamiltonian path count."""
    nn = n
    dp = [0] * ((1 << nn) * nn)
    for v in range(nn):
        dp[(1 << v) * nn + v] = 1
    for S in range(1, 1 << nn):
        for v in range(nn):
            if not (S & (1 << v)):
                continue
            val = dp[S * nn + v]
            if val == 0:
                continue
            for u in range(nn):
                if S & (1 << u):
                    continue
                if adj_flat[v * nn + u]:
                    dp[(S | (1 << u)) * nn + u] += val
    full = (1 << nn) - 1
    return sum(dp[full * nn + v] for v in range(nn))

# =============================================================================
# SAMPLING
# =============================================================================

N_samples = 500000
print(f"  Sampling {N_samples:,} random tournaments...")
t0 = time.time()

random.seed(42)

# Data collection
score_data = defaultdict(lambda: {'count': 0, 'H_sum': 0, 'H2_sum': 0,
                                   'H_min': 10**9, 'H_max': 0, 'H_vals': Counter()})
all_H = []

for sample_idx in range(N_samples):
    # Generate random tournament
    adj_flat = [0] * (n * n)
    scores = [0] * n
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                adj_flat[i * n + j] = 1
                scores[i] += 1
            else:
                adj_flat[j * n + i] = 1
                scores[j] += 1
            idx += 1

    sc = tuple(sorted(scores))
    H = H_dp_fast(adj_flat, n)

    d = score_data[sc]
    d['count'] += 1
    d['H_sum'] += H
    d['H2_sum'] += H * H
    d['H_min'] = min(d['H_min'], H)
    d['H_max'] = max(d['H_max'], H)
    d['H_vals'][H] += 1
    all_H.append(H)

    if (sample_idx + 1) % 100000 == 0:
        elapsed = time.time() - t0
        print(f"    {sample_idx + 1:,} samples in {elapsed:.1f}s ({(sample_idx+1)/elapsed:.0f}/s)")

elapsed = time.time() - t0
print(f"  Done: {N_samples:,} samples in {elapsed:.1f}s")
print()

# =============================================================================
# GLOBAL STATISTICS
# =============================================================================

print("=" * 76)
print("  GLOBAL STATISTICS")
print("=" * 76)
print()

all_H_arr = np.array(all_H, dtype=float)
E_H_sample = np.mean(all_H_arr)
E_H_exact = factorial(n) / 2**(n-1)
Var_H_sample = np.var(all_H_arr)

print(f"  E[H] exact  = {E_H_exact}")
print(f"  E[H] sample = {E_H_sample:.4f}")
print(f"  Var(H) sample = {Var_H_sample:.4f}")
print(f"  H range: [{min(all_H)}, {max(all_H)}]")
print()

# Count score classes observed
n_scores_observed = len(score_data)
print(f"  Score classes observed: {n_scores_observed} (expected 167 for n=8)")
print()

# =============================================================================
# WITHIN-CLASS VARIANCE (OCR ESTIMATION)
# =============================================================================

print("=" * 76)
print("  OCR ESTIMATION FROM SAMPLE")
print("=" * 76)
print()

# Compute within-class variance
Var_within = 0.0
total_weight = 0.0

for sc, d in score_data.items():
    ns = d['count']
    if ns < 2:
        continue
    E_H_s = d['H_sum'] / ns
    E_H2_s = d['H2_sum'] / ns
    Var_s = E_H2_s - E_H_s ** 2
    weight = ns / N_samples
    Var_within += weight * Var_s
    total_weight += weight

OCR_est = 1 - Var_within / Var_H_sample if Var_H_sample > 0 else 1.0

print(f"  Var(H|scores) ≈ {Var_within:.4f}")
print(f"  Var(H) ≈ {Var_H_sample:.4f}")
print(f"  OCR ≈ {OCR_est:.8f}")
print(f"  1 - OCR ≈ {1 - OCR_est:.8f}")
print()

# Comparison with previous n
print(f"  OCR SEQUENCE (exact for n≤7, estimated for n=8):")
ocr_seq = [(3, 1.0), (4, 1.0), (5, 0.969925), (6, 0.959056), (7, 0.958107)]
for nn, ocr in ocr_seq:
    print(f"    n={nn}: OCR = {ocr:.6f}")
print(f"    n=8: OCR ≈ {OCR_est:.6f} (estimated from {N_samples:,} samples)")
print()

# (1-OCR)*n trend
print(f"  (1-OCR)×n trend:")
for nn, ocr in ocr_seq:
    print(f"    n={nn}: (1-OCR)×n = {(1-ocr)*nn:.4f}")
print(f"    n=8: (1-OCR)×n ≈ {(1-OCR_est)*8:.4f}")
print()

# =============================================================================
# SCORE CLASS TABLE
# =============================================================================

print("=" * 76)
print("  SCORE CLASSES — SORTED BY CONTRIBUTION TO RESIDUAL")
print("=" * 76)
print()

# Compute S₂ and contribution for each class
contributions = []
for sc, d in score_data.items():
    ns = d['count']
    if ns < 5:
        continue
    S2 = sum((s - (n - 1) / 2) ** 2 for s in sc)
    E_H_s = d['H_sum'] / ns
    Var_s = d['H2_sum'] / ns - E_H_s ** 2
    n_H = len(d['H_vals'])
    weight = ns / N_samples
    contrib = weight * Var_s

    # Self-complementary?
    comp = tuple(sorted(n - 1 - s for s in sc))
    is_sc = (sc == comp)

    contributions.append({
        'sc': sc, 'count': ns, 'S2': S2, 'E_H': E_H_s, 'Var': Var_s,
        'n_H': n_H, 'H_min': d['H_min'], 'H_max': d['H_max'],
        'contrib': contrib, 'is_sc': is_sc, 'weight': weight,
    })

contributions.sort(key=lambda x: -x['contrib'])

# Print top contributors
print(f"  TOP 20 CONTRIBUTORS TO RESIDUAL:")
print()
header = (f"  {'Rank':>4} | {'Score':>30} | {'Count':>6} | {'%':>5} | {'S₂':>5} | {'SC':>2} | "
          f"{'#H':>3} | {'H range':>12} | {'Var(H|s)':>10} | {'Contrib%':>8}")
print(header)
print("  " + "-" * (len(header) - 2))

cumulative = 0
for rank, c in enumerate(contributions[:20], 1):
    pct = c['weight'] * 100
    H_range = f"[{c['H_min']},{c['H_max']}]"
    cum_contrib = cumulative + c['contrib']
    cum_pct = cum_contrib / Var_within * 100 if Var_within > 0 else 0
    cumulative = cum_contrib

    print(f"  {rank:4d} | {str(c['sc']):>30} | {c['count']:6d} | {pct:4.1f}% | {c['S2']:5.1f} | "
          f"{'✓' if c['is_sc'] else ' ':>2} | {c['n_H']:3d} | {H_range:>12} | {c['Var']:10.2f} | {cum_pct:7.1f}%")

print()
print(f"  Top 1:  {contributions[0]['contrib']/Var_within*100:.1f}% of residual")
print(f"  Top 5:  {sum(c['contrib'] for c in contributions[:5])/Var_within*100:.1f}%")
print(f"  Top 10: {sum(c['contrib'] for c in contributions[:10])/Var_within*100:.1f}%")
print(f"  Top 20: {sum(c['contrib'] for c in contributions[:20])/Var_within*100:.1f}%")
print()

# =============================================================================
# SCORE CLASS SUMMARY
# =============================================================================

print("=" * 76)
print("  SCORE CLASS SUMMARY")
print("=" * 76)
print()

n_ambiguous = sum(1 for c in contributions if c['n_H'] > 1)
n_total_classes = len(contributions)
print(f"  Total score classes observed: {n_total_classes}")
print(f"  Ambiguous (multiple H values): {n_ambiguous} ({n_ambiguous/n_total_classes*100:.0f}%)")
print()

# Size distribution
sizes = sorted([c['count'] for c in contributions], reverse=True)
print(f"  Class sizes: largest = {sizes[0]:,}, smallest = {sizes[-1]:,}, median = {sizes[len(sizes)//2]:,}")
print()

# SC analysis
sc_classes = [c for c in contributions if c['is_sc']]
non_sc_classes = [c for c in contributions if not c['is_sc']]
if sc_classes and non_sc_classes:
    avg_sc = sum(c['count'] for c in sc_classes) / len(sc_classes)
    avg_non_sc = sum(c['count'] for c in non_sc_classes) / len(non_sc_classes)
    print(f"  Self-complementary: {len(sc_classes)} classes, avg size {avg_sc:.0f}")
    print(f"  Non-self-comp: {len(non_sc_classes)} classes, avg size {avg_non_sc:.0f}")
    print(f"  SC ratio: {avg_sc/avg_non_sc:.2f}×")
print()

# =============================================================================
# THE PoS PREDICTION CHECK
# =============================================================================

print("=" * 76)
print("  PoS PREDICTION CHECK AT n=8")
print("=" * 76)
print()

print("""
  From S107/S110, the predicted dominant PoS at n=8 should have:
  - Almost-source/sink structure
  - Self-complementary
  - Near-regular (low S₂)

  The near-regular SC scores at n=8 (even n, no regular exists):
  Closest to regular: all scores near (n-1)/2 = 3.5
  → scores like (3,3,3,3,4,4,4,4) or (2,3,3,4,4,4,4,4)...
  Wait: n=8 is even, so scores are half-integers from 3.5.
  Possible near-regular: (3,3,3,3,4,4,4,4) with S₂ = 8×0.25 = 2.

  Let me check which score class has the most residual contribution.
""")

if contributions:
    top = contributions[0]
    print(f"  ACTUAL dominant class: {top['sc']}")
    print(f"    S₂ = {top['S2']}")
    print(f"    SC = {top['is_sc']}")
    print(f"    Count = {top['count']:,}")
    print(f"    #H values = {top['n_H']}")
    print(f"    H range = [{top['H_min']}, {top['H_max']}]")
    print(f"    Var(H|s) = {top['Var']:.2f}")
    print(f"    Contributes {top['contrib']/Var_within*100:.1f}% of residual")

print()

# =============================================================================
# EXACT QUANTITIES
# =============================================================================

print("=" * 76)
print("  EXACT AND ESTIMATED QUANTITIES")
print("=" * 76)
print()

print(f"  EXACT:")
print(f"    E[H] = {factorial(n)}/{2**(n-1)} = {factorial(n) // 2**(n-1)}")
print(f"    #score sequences = 167 (OEIS A000571)")
print(f"    #iso classes = 6880 (OEIS A000568)")
print()

print(f"  ESTIMATED (from {N_samples:,} samples):")
print(f"    E[H] ≈ {E_H_sample:.2f} (exact: {E_H_exact})")
print(f"    Var(H) ≈ {Var_H_sample:.2f}")
print(f"    Var(H|scores) ≈ {Var_within:.4f}")
print(f"    OCR ≈ {OCR_est:.6f}")
print(f"    1-OCR ≈ {1-OCR_est:.6f}")
print(f"    (1-OCR)×n ≈ {(1-OCR_est)*n:.4f}")
print()

# E[H] sequence (exact)
print(f"  E[H] SEQUENCE (exact):")
for nn in range(3, 9):
    eh = factorial(nn) // 2**(nn-1)
    print(f"    n={nn}: E[H] = {eh}")
print()

# Var(H) sequence
print(f"  Var(H) SEQUENCE:")
var_seq = {3: 0.75, 4: 3.0, 5: 17.8125, 6: 146.25, 7: 1611.91}
for nn in sorted(var_seq.keys()):
    print(f"    n={nn}: Var(H) = {var_seq[nn]:.4f}")
print(f"    n=8: Var(H) ≈ {Var_H_sample:.2f}")
if 7 in var_seq:
    ratio = Var_H_sample / var_seq[7]
    print(f"    Ratio Var(8)/Var(7) ≈ {ratio:.2f}")
print()

print("opus-2026-03-21-S114: SAMPLING-BASED OCR AT n=8")
