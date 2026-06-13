#!/usr/bin/env python3
"""
n7_integration_s20bz.py -- kind-pasteur-2026-03-22-S20bz

INTEGRATING THE n=7 RESULTS FROM OPUS S212.

E(G_7) = 4086 CONFIRMED.
Edge sequence: 1, 5, 30, 290, 4086.

Update all formulas and predictions with the new data point.

Author: kind-pasteur-2026-03-22-S20bz
"""
import sys
from math import comb, factorial, sqrt, pi, log
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def fiber_frac(n):
    k = n - 2
    result = Fraction(1)
    for j in range(k):
        result *= Fraction(1 + 2*j, 2*(j+1))
    return result

def A000568(n):
    """Tournament iso class count from Burnside."""
    known = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536, 10:9733056}
    return known.get(n, None)

print("=" * 70)
print("  N=7 INTEGRATION: UPDATING ALL FORMULAS")
print("=" * 70)

# ================================================================
# 1. THE COMPLETE EDGE SEQUENCE
# ================================================================
print(f"\n{'='*70}")
print(f"  1. THE COMPLETE EDGE SEQUENCE + PREDICTIONS")
print(f"{'='*70}\n")

edges = {3: 1, 4: 5, 5: 30, 6: 290, 7: 4086}

print(f"  {'n':>3s} {'V':>8s} {'E':>8s} {'m':>4s} {'f':>8s} {'E_pred':>8s} {'c(n)':>8s} {'eps':>8s}")
for n in range(3, 8):
    V = A000568(n)
    E = edges[n]
    m = comb(n, 2)
    f = float(fiber_frac(n))
    E_pred = 0.5 * V * m * (1 - f)
    c = E / E_pred
    eps = 1 - c
    print(f"  {n:>3d} {V:>8d} {E:>8d} {m:>4d} {f:>8.4f} {E_pred:>8.0f} {c:>8.4f} {eps:>8.4f}")

# The correction sequence: 0.667, 0.667, 0.727, 0.950, ???
# At n=7: c = 4086 / (0.5 * 456 * 21 * (1 - 0.2461)) = 4086 / 3610 = 1.132?!
# Wait: c > 1? Let me recompute.
n = 7
V = 456
E = 4086
m = 21
f = float(fiber_frac(7))
E_pred = 0.5 * V * m * (1 - f)
c = E / E_pred
print(f"\n  CORRECTION AT n=7: c = {c:.4f}")
print(f"  E_pred = {E_pred:.0f}, E_actual = {E}")
print(f"  c > 1! The approximation UNDERESTIMATES at n=7!")

# ================================================================
# 2. WHY c > 1 AT n=7?
# ================================================================
print(f"\n{'='*70}")
print(f"  2. WHY THE FORMULA UNDERESTIMATES AT n=7")
print(f"{'='*70}\n")

# The formula E ~ V*m*(1-f)/2 assumes:
# - Each class has degree ~ m*(1-f) (= cross-arcs per class)
# - No duplicate neighbors (each cross-arc goes to a different class)
# - All classes have |Aut|=1

# At n=7: avg_degree = 2*4086/456 = 17.92
# Predicted avg_degree = m*(1-f) = 21*(1-0.2461) = 15.83
# Actual > predicted! This means classes have MORE distinct neighbors
# than the naive cross-arc count suggests.

# This can happen if:
# - Some |Aut|>1 classes have HIGHER degree than expected
# - The Aut-orbit structure creates FEWER self-loops than predicted
# - The self-loop fraction per class is LOWER than f(n)

print(f"  Average degree at n=7: {2*E/V:.2f}")
print(f"  Predicted avg degree: {m*(1-f):.2f}")
print(f"  Ratio: {2*E/(V*m*(1-f)):.4f}")

# The issue: the self-loop fraction f(n) is the AVERAGE over all
# labeled tournaments. But the self-loop fraction per ISO CLASS
# may be different (weighted differently).

# Opus S212 reports: T_7 (transition orbits) = 8912.
# avg_deg = 8912/456 = wait, no: avg_deg = 2*E/V = 8172/456 = 17.92.
# T_7/V = 8912/456 = 19.54. And T_7/(2*E) = 8912/8172 = 1.09.

# From opus: T_n/E_n = 4.0, 3.2, 2.93, 2.43, 2.18
# This is the "orbit-to-edge" ratio: on average each edge is
# represented by T/E arc-orbit types.

# The self-loop count at n=7: 201 self-loop classes (from opus)
# Wait: is that 201 self-loop EDGES or classes with self-loops?
# From the output: "Self-loops: 201" -- this is likely the number
# of iso class pairs that are self-loops (i.e., 201 arcs that
# produce a tournament in the same class).

# ================================================================
# 3. THE UPDATED SEQUENCE TABLE
# ================================================================
print(f"\n{'='*70}")
print(f"  3. THE COMPLETE UPDATED TABLE")
print(f"{'='*70}\n")

print(f"  {'Quantity':>25s}", end="")
for n in range(3, 8):
    print(f"  {'n='+str(n):>8s}", end="")
print()

rows = {
    'Vertices': [2, 4, 12, 56, 456],
    'Edges': [1, 5, 30, 290, 4086],
    'Width': [1, 2, 3, 6, '?'],
    'Level edges': [0, 0, 1, 15, '?'],
    'Self-loops': ['?', '?', '?', '?', 201],
    'avg_degree': [1.0, 2.5, 5.0, 10.36, 17.92],
    'Density': [1.000, 0.833, 0.455, 0.188, 0.039],
    'H-levels': [2, 3, 7, 19, '?'],
    'SC classes': [2, 2, 8, 12, '?'],
    'Chains': [1, 3, 99, 292510, '?'],
}

for name, vals in rows.items():
    print(f"  {name:>25s}", end="")
    for v in vals:
        print(f"  {str(v):>8s}", end="")
    print()

# ================================================================
# 4. VERIFY WIDTH PREDICTION FOR n=7
# ================================================================
print(f"\n{'='*70}")
print(f"  4. WIDTH PREDICTION: C(5,2) = 10 at n=7")
print(f"{'='*70}\n")

# Width = C(n-2, floor((n-2)/2)) = C(5, 2) = 10.
# This predicts the maximum number of iso classes at the same H value.
# We need to check this against the actual n=7 data.
# Opus didn't report the width directly, but the H-levels are important.

print(f"  Predicted width at n=7: C(5, 2) = {comb(5, 2)}")
print(f"  This means: at most 10 iso classes share the same H value.")
print(f"  Need to verify from the n=7 computation data.")

# ================================================================
# 5. IMPROVED EDGE FORMULA
# ================================================================
print(f"\n{'='*70}")
print(f"  5. IMPROVED EDGE FORMULA")
print(f"{'='*70}\n")

# The simple formula E ~ V*m*(1-f)/2 is:
# n=3: 67%, n=4: 67%, n=5: 73%, n=6: 95%, n=7: 88%
# The correction factor is NOT monotonically increasing!
# c = 0.667, 0.667, 0.727, 0.950, 1.132

# Actually c > 1 at n=7 means the formula UNDER-predicts.
# This is unexpected. Let me use opus's T_n formula instead.

# Opus: T_n = Burnside-style count of transition orbits.
# T_n = (1/n!) * sum_{sigma} Fix(sigma) * C(f(sigma), 2)
# where f(sigma) = Fix(sigma) for tournaments.

# T_n and E_n: T_n/E_n = 4.0, 3.2, 2.93, 2.43, 2.18
# This ratio is DECREASING toward 2. If T_n/E_n -> 2:
# E_n ~ T_n/2 asymptotically.

# T_n values: 4, 16, 88, 704, 8912
# Check: T_n/V_n = 2.0, 4.0, 7.33, 12.57, 19.54
# Compare with m = 3, 6, 10, 15, 21
# T_n/V_n < m always. The deficit = m - T_n/V_n = 1.0, 2.0, 2.67, 2.43, 1.46

# Actually T_n is the number of DIRECTED arc-class transitions
# (ordered pairs: class C with specific arc orbit type).
# E = number of UNORDERED pairs of distinct classes connected by a flip.
# T_n includes self-loops and multi-edges.

# T_n = 2*E + self_loop_orbits + multi_edge_excess
# From opus: T - 2E = 2, 6, 28, 124, 740
# This is the self-loop + multi-edge overhead.

print(f"  OPUS'S TRANSITION ORBIT FORMULA:")
T = {3: 4, 4: 16, 5: 88, 6: 704, 7: 8912}
for n in range(3, 8):
    V = A000568(n)
    E = edges[n]
    t = T[n]
    overhead = t - 2*E
    print(f"  n={n}: T={t}, 2E={2*E}, overhead(T-2E)={overhead}, T/(2E)={t/(2*E):.4f}")

# Overhead sequence: 2, 6, 28, 124, 740
# Ratios: 3, 4.67, 4.43, 5.97
# Is overhead ~ something?
overhead_seq = [T[n] - 2*edges[n] for n in range(3, 8)]
print(f"\n  Overhead sequence: {overhead_seq}")
print(f"  Overhead/V: {[overhead_seq[i]/A000568(i+3) for i in range(5)]}")

# ================================================================
# 6. THE REFINED PREDICTION FOR n=8
# ================================================================
print(f"\n{'='*70}")
print(f"  6. PREDICTIONS FOR n=8")
print(f"{'='*70}\n")

# Using T_n/E_n ~ 2.18 at n=7, extrapolate T_n/E_n at n=8.
# The ratio decreases roughly as: 4.0, 3.2, 2.93, 2.43, 2.18
# Differences: -0.8, -0.27, -0.50, -0.25
# The ratio seems to converge to 2, approaching from above.

# At n=8: T_8/E_8 ~ 2.10 (extrapolation)
# If we had T_8, we could compute E_8 = T_8 / 2.10.
# T_n grows very fast. From the data:
# T_n = 4, 16, 88, 704, 8912
# Ratios: 4, 5.5, 8.0, 12.66
# These ratios are growing. T_n might grow as V_n * m / something.

# T_n/V_n = 2.0, 4.0, 7.33, 12.57, 19.54
# Compare m = 3, 6, 10, 15, 21
# T_n/V_n ~ m - small_correction
# Correction = 1.0, 2.0, 2.67, 2.43, 1.46 (NOT monotone!)

# At n=8: m=28. If T/V ~ m-1.2 = 26.8. V=6880.
# T_8 ~ 6880 * 26.8 = 184,384.
# E_8 ~ T_8/2.10 = 87,802.
# Previously I predicted E_8 ~ 74,592 from V*m*(1-f)/2.

print(f"  PREDICTION FOR n=8:")
n8_V = 6880
n8_m = 28
n8_f = float(fiber_frac(8))

# Method 1: V*m*(1-f)/2
E8_m1 = 0.5 * n8_V * n8_m * (1 - n8_f)
print(f"  Method 1 (V*m*(1-f)/2): {E8_m1:.0f}")

# Method 2: T/ratio extrapolation
n8_TV = n8_m - 1.2  # extrapolated T/V
n8_T = n8_V * n8_TV
n8_ratio = 2.10  # extrapolated T/E
E8_m2 = n8_T / n8_ratio
print(f"  Method 2 (T/ratio extrap): T~{n8_T:.0f}, E~{E8_m2:.0f}")

# Method 3: average of methods
E8_m3 = (E8_m1 + E8_m2) / 2
print(f"  Method 3 (average): {E8_m3:.0f}")
print(f"  Range: [{min(E8_m1, E8_m2):.0f}, {max(E8_m1, E8_m2):.0f}]")
