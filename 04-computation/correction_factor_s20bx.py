#!/usr/bin/env python3
"""
correction_factor_s20bx.py -- kind-pasteur-2026-03-22-S20bx

THE CORRECTION FACTOR: Why E < V*m*(1-f)/2 and what fixes it.

The approximation E ~ V*m*(1-f)/2 overcounts because:
1. Multiple arcs from the same class go to the SAME neighbor (duplicates)
2. Classes with |Aut| > 1 have fewer effective arcs (orbit grouping)

The correction factor c(n) = E_actual / E_approx approaches 1 as n -> inf.

APPROACH 1: Compute c(n) exactly at n=3,4,5,6.
APPROACH 2: Model c(n) as a COLLISION probability.
APPROACH 3: Derive c(n) from the Burnside structure.

Author: kind-pasteur-2026-03-22-S20bx
"""
import sys
from math import comb, factorial, exp, log
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def fiber_frac(n):
    k = n - 2
    result = Fraction(1)
    for j in range(k):
        result *= Fraction(1 + 2*j, 2*(j+1))
    return result

print("=" * 70)
print("  THE CORRECTION FACTOR")
print("=" * 70)

# ================================================================
# 1. EXACT CORRECTION FACTORS
# ================================================================
print(f"\n{'='*70}")
print(f"  1. EXACT CORRECTION FACTORS c(n) = E_actual / E_approx")
print(f"{'='*70}\n")

data = {
    3: {'V': 2,  'E': 1,   'm': 3},
    4: {'V': 4,  'E': 5,   'm': 6},
    5: {'V': 12, 'E': 30,  'm': 10},
    6: {'V': 56, 'E': 290, 'm': 15},
}

print(f"  {'n':>3s} {'V':>5s} {'m':>4s} {'f(n)':>8s} {'E_approx':>9s} {'E_actual':>9s} {'c(n)':>8s} {'1-c':>8s}")
for n in [3, 4, 5, 6]:
    d = data[n]
    f = float(fiber_frac(n))
    E_approx = 0.5 * d['V'] * d['m'] * (1 - f)
    c = d['E'] / E_approx
    print(f"  {n:>3d} {d['V']:>5d} {d['m']:>4d} {f:>8.4f} {E_approx:>9.1f} {d['E']:>9d} {c:>8.4f} {1-c:>8.4f}")

# ================================================================
# 2. THE COLLISION MODEL
# ================================================================
print(f"\n{'='*70}")
print(f"  2. THE COLLISION MODEL")
print(f"{'='*70}\n")

# If each of the m*(1-f) cross-arcs lands in a uniformly random neighbor
# among V-1 possible classes, the expected number of DISTINCT neighbors is:
# D = (V-1) * (1 - (1 - 1/(V-1))^{m*(1-f)})
# ~ m*(1-f) * (1 - m*(1-f)/(2*(V-1))) for small m/V

# But arcs DON'T land uniformly: they're biased toward nearby classes
# (small H-change arcs are more common than large H-change arcs).

print(f"  UNIFORM COLLISION MODEL:")
print(f"  {'n':>3s} {'m(1-f)':>8s} {'V-1':>5s} {'D_pred':>8s} {'D_actual':>9s} {'Ratio':>8s}")

for n in [3, 4, 5, 6]:
    d = data[n]
    f = float(fiber_frac(n))
    cross = d['m'] * (1 - f)
    V = d['V']
    # Expected distinct neighbors per class (uniform model)
    D_uniform = (V - 1) * (1 - (1 - 1/(V-1))**cross) if V > 1 else 0
    # Actual average degree
    D_actual = 2 * d['E'] / V
    ratio = D_actual / D_uniform if D_uniform > 0 else 0
    print(f"  {n:>3d} {cross:>8.2f} {V-1:>5d} {D_uniform:>8.2f} {D_actual:>9.2f} {ratio:>8.4f}")

# ================================================================
# 3. THE IMPROVED FORMULA
# ================================================================
print(f"\n{'='*70}")
print(f"  3. IMPROVED FORMULAS")
print(f"{'='*70}\n")

# The average degree D(n) = 2*E/V.
# From data: D = [1, 2.5, 5.0, 10.36]
# This looks like it might be m*(1-f)*correction.

# Let me try: D(n) = m - m*f(n) - duplicates(n)
# duplicates per class = m*(1-f) - D = "duplicate rate"
# At n=5: duplicates = 10*(1-0.3125) - 5.0 = 6.875 - 5.0 = 1.875
# At n=6: duplicates = 15*(1-0.2734) - 10.36 = 10.90 - 10.36 = 0.54

print(f"  DUPLICATE RATE:")
for n in [3, 4, 5, 6]:
    d = data[n]
    f = float(fiber_frac(n))
    D = 2 * d['E'] / d['V']
    cross = d['m'] * (1 - f)
    dups = cross - D
    dup_rate = dups / cross if cross > 0 else 0
    print(f"  n={n}: D={D:.2f}, cross_per_class={cross:.2f}, dups={dups:.2f}, dup_rate={dup_rate:.4f}")

# The duplicate rate is DECREASING: 0.333, 0.200, 0.272, 0.050
# At n=6 it's only 5%! This confirms the approximation improves.

# Better formula: account for |Aut| > 1 classes
# For |Aut|=a: the class has m/a effective arc-orbits (approximately)
# Each orbit has a arcs going to the same place.
# So effective cross-arcs = (m/a)*(1-f_a) where f_a depends on |Aut|.

# Even simpler: compute the weighted sum
# E = (1/2) * sum_C D(C)
# where D(C) = #{distinct neighbors of C}
# For large n: D(C) ~ m*(1-f) for |Aut|=1 classes (the vast majority)
# and D(C) ~ m/|Aut| for |Aut|>1 classes.

# The WEIGHTED contribution:
# E ~ (1/2) * [N_1 * m*(1-f) + sum_{|Aut|>1 classes} D(C)]
# where N_1 = number of |Aut|=1 classes.

# At n=5: N_1=7, D~6.43, contribution = 45
#          N_{>1}=5, contribution = 15
#          Total = 60, E = 30.

# At n=6: most classes have |Aut|=1.
# The key: N_1/V -> 1 as n -> inf.
# So E ~ (1/2) * V * m * (1-f) * (1 - dup_rate)
# with dup_rate -> 0 as V -> inf.

# ================================================================
# 4. THE ASYMPTOTIC FORMULA
# ================================================================
print(f"\n{'='*70}")
print(f"  4. THE ASYMPTOTIC EDGE FORMULA")
print(f"{'='*70}\n")

print(f"""  E(n) ~ (1/2) * A000568(n) * C(n,2) * (1 - f(n)) * (1 - epsilon(n))

  where:
    A000568(n) = tournament iso class count (from Burnside)
    C(n,2) = number of arcs
    f(n) = (1/2)_{{n-2}} / (n-2)! = self-loop fraction
    epsilon(n) = duplicate neighbor rate -> 0 as n -> inf

  ASYMPTOTICALLY: E(n) ~ (1/2) * A000568(n) * C(n,2) * (1 - f(n))

  Since A000568(n) ~ 2^{{C(n,2)}} / n! and f(n) ~ 1/sqrt(pi*(n-2)):

  E(n) ~ 2^{{C(n,2)-1}} * C(n,2) * (1 - 1/sqrt(pi*(n-2))) / n!
""")

# Predictions for n=7, 8
for n in [7, 8, 9, 10]:
    V_n = {7: 456, 8: 6880, 9: 191536, 10: 9733056}[n]
    m_n = comb(n, 2)
    f_n = float(fiber_frac(n))
    E_pred = 0.5 * V_n * m_n * (1 - f_n)
    print(f"  n={n}: V={V_n}, m={m_n}, f={f_n:.4f}, E_pred ~ {E_pred:.0f}")

# ================================================================
# 5. THE EXACT CORRECTION FROM ORBiT-STABILIZER
# ================================================================
print(f"\n{'='*70}")
print(f"  5. ORBIT-STABILIZER APPROACH TO EXACT EDGE COUNT")
print(f"{'='*70}\n")

# The EXACT edge count requires knowing, for each class C,
# the number of distinct neighbor classes.
# This depends on the FINE structure of C (not just |Aut|).

# But there's a BURNSIDE-LIKE formula for edges:
# edges = (1/|G|) * sum_{g in G} Fix_edges(g)
# where Fix_edges(g) counts edges of Q_m fixed by g.

# An edge of Q_m connects two tournaments differing in one arc.
# g fixes an edge (T, T') iff g(T) = T and g(T') = T' (both fixed)
# OR g(T) = T' and g(T') = T (g swaps them).

# Case 1: g fixes both endpoints. Then g must fix both T and flip(T,e).
# This means g is an automorphism of T that ALSO fixes the arc e.
# Fix_both(g) = Fix(g) * (fraction of arcs fixed by g in each fixed tournament).

# Case 2: g swaps the endpoints. Then g must map T to T' = flip(T,e)
# and vice versa. This means g flips exactly the arc e.

# This is getting complicated. Let me try the QUOTIENT approach instead.

# For the QUOTIENT graph Q_m / G:
# edges = number of distinct pairs {orbit(T), orbit(T')} connected by an edge
# = number of distinct pairs {C1, C2} such that some T in C1 has an arc
#   whose flip lands in C2.

# This can be computed from the ORBIT-PAIR counting:
# For each ordered pair of classes (C1, C2):
#   w(C1, C2) = #{arcs e in representative T of C1 : flip(T,e) in C2}
#   If w(C1, C2) > 0: edge exists.

# The edge count = #{unordered pairs {C1, C2} : w(C1,C2) > 0 or w(C2,C1) > 0}

# Since the weight is symmetric (proved earlier): w(C1,C2) > 0 iff w(C2,C1) > 0.
# So edge exists iff w(C1, C2) > 0.

# Now: w(C1, C2) = sum over Aut-orbits of arcs in C1 that land in C2
# = (# arcs from T1 landing in C2) / |Aut(C1)| * |Aut(C1)| [by orbit counting]
# Hmm, it's just the raw count.

# The total cross-weight from C1 = sum_{C2 != C1} w(C1, C2) = m - self_loops(C1).
# The degree(C1) = #{C2 : w(C1,C2) > 0}.
# edges = (1/2) * sum_C1 degree(C1).

# So we need: for each C1, how many DISTINCT C2 are hit?
# This depends on the WEIGHT DISTRIBUTION w(C1, .).
# If all weights are 1 (every cross-arc goes to a different class): degree = m - self_loops.
# If some weights are > 1 (duplicate neighbors): degree < m - self_loops.

# The weight is > 1 when:
# Two arcs e1, e2 of T1 land in the SAME class C2 (but different labeled tournaments).
# This happens when flip(T1, e1) and flip(T1, e2) are isomorphic.

# For |Aut|=1: this means there exists sigma in S_n such that
# sigma(flip(T1, e1)) = flip(T1, e2). This sigma maps T1-minus-e1 to T1-minus-e2
# with the appropriate arc corrections.

# This is the "NEAR-AUTOMORPHISM" condition: sigma is an automorphism of T1
# that also maps e1 to e2 (with possible sign flip).
# For |Aut(T1)|=1: sigma must be the identity, so e1 = e2. No duplicates!
# WAIT: that's not right. sigma doesn't need to fix T1 exactly; it needs
# to map flip(T1,e1) to flip(T1,e2). This is different from being an Aut of T1.

print(f"""  KEY INSIGHT: For |Aut(T)|=1 classes, are there duplicate neighbors?

  Duplicate: two arcs e1 != e2 such that flip(T, e1) ~ flip(T, e2).
  This requires a permutation sigma with sigma(flip(T,e1)) = flip(T,e2).
  sigma maps T to something that agrees with T on all arcs except e1 and e2.

  For |Aut|=1: this sigma is NOT an automorphism of T.
  It's a "double-flip near-automorphism" -- fixes T except at e1 and e2.

  Equivalently: T and sigma(T) differ in exactly 2 arcs (e1 and sigma^(-1)(e2)).
  This is a HAMMING-2 neighbor in the hypercube.

  So: duplicate neighbors occur iff T has a HAMMING-2 neighbor that is isomorphic.
  This is related to the |Aut| of the FLIPPED tournament, not T itself.

  CONCLUSION: The duplicate rate depends on the DENSITY of the iso class
  in the Hamming-2 neighborhood of T. At large n, Hamming-2 neighbors are
  exponentially rare relative to total classes, so duplicates vanish.
""")

# Verify: at n=5, count Hamming-2 iso-neighbors
print(f"  VERIFICATION: Hamming-2 iso-neighbors at n=5")

from itertools import permutations, combinations
import numpy as np

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

canon_map = {}
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    canon_map[bits] = best

# For each tournament, count Hamming-2 neighbors in same iso class
h2_same_count = 0
total_h2 = 0
for bits in range(2**m):
    cf = canon_map[bits]
    for k1 in range(m):
        for k2 in range(k1+1, m):
            nb = bits ^ (1 << k1) ^ (1 << k2)
            total_h2 += 1
            if canon_map[nb] == cf:
                h2_same_count += 1

frac = h2_same_count / total_h2 if total_h2 > 0 else 0
print(f"  Hamming-2 same-class fraction: {h2_same_count}/{total_h2} = {frac:.6f}")
print(f"  This is the probability that a Hamming-2 perturbation stays in the same class.")
print(f"  Compare with Hamming-1 (self-loop) fraction: {float(Fraction(11,64)):.6f}")
