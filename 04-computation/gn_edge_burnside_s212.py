#!/usr/bin/env python3
"""
gn_edge_burnside_s212.py — Burnside approach to |E(G_n)|
opus-2026-03-22-S212

KEY INSIGHT: the transition orbit count T_n = (1/n!) Σ_σ Fix_T(σ) × Fix_arcs(σ)
counts orbits of S_n on (tournament, arc-position) pairs. Each orbit is either:
  - a self-loop orbit (flip stays in same class)
  - a cross-class orbit (flip changes class)

Multiple cross-class orbits can map to the SAME edge of G_n.
The number of cross-class orbits mapping to edge (i,j) =
  number of distinct Aut-orbit types of arc positions producing this transition.

For n=3,4: exactly 2 cross-class orbits per edge!
For n=5: 72/30 = 2.4 orbits per edge (average).

QUESTION: is there a formula for the self-loop orbit count S_n?
If so, |E| is bounded by T_n - S_n ≥ |E|.

Also: can we compute S_n via Burnside?
S_n = (1/n!) Σ_σ |{(T,e) self-loop : σ(T,e) = (T,e)}|
"""

import sys
from math import comb, factorial, gcd
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  BURNSIDE APPROACH TO |E(G_n)|")
print("  opus-2026-03-22-S212")
print("=" * 80)

def partitions(n, max_part=None):
    if max_part is None:
        max_part = n
    if n == 0:
        yield []
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield [first] + rest

def cycle_type_count(n, cycle_type):
    ct = Counter(cycle_type)
    result = factorial(n)
    for length, count in ct.items():
        result //= (length ** count) * factorial(count)
    return result

def fix_tournaments(cycle_type):
    for c in cycle_type:
        if c % 2 == 0:
            return 0
    k = len(cycle_type)
    exp = sum((c - 1) // 2 for c in cycle_type)
    for i in range(k):
        for j in range(i + 1, k):
            exp += gcd(cycle_type[i], cycle_type[j])
    return 2 ** exp

def fix_arcs(cycle_type, n):
    """Number of arc positions fixed by σ with given cycle type.
    {u,v} is fixed iff both u,v are fixed points, or they form a 2-cycle."""
    f = cycle_type.count(1)
    t = cycle_type.count(2)
    return comb(f, 2) + t

def fix_arcs_detailed(cycle_type, n):
    """For each fixed arc position, count how many are:
    (a) both endpoints fixed: C(f, 2)
    (b) endpoints forming a 2-cycle: t"""
    f = cycle_type.count(1)
    t = cycle_type.count(2)
    return comb(f, 2), t

# ============================================================================
# PART 1: Transition orbit counts
# ============================================================================

print("\n  PART 1: Transition orbit counts T_n = (1/n!) Σ Fix_T × Fix_arcs")
print("-" * 60)

for n in range(3, 12):
    m = comb(n, 2)
    nfact = factorial(n)

    T_sum = 0
    V_sum = 0
    # Also decompose by arc-fix type
    type_a_sum = 0  # both endpoints fixed
    type_b_sum = 0  # endpoints form 2-cycle

    for part in partitions(n):
        ct = cycle_type_count(n, part)
        ft = fix_tournaments(part)
        fa = fix_arcs(part, n)
        fa_a, fa_b = fix_arcs_detailed(part, n)

        T_sum += ct * ft * fa
        V_sum += ct * ft
        type_a_sum += ct * ft * fa_a
        type_b_sum += ct * ft * fa_b

    T_n = T_sum // nfact
    V_n = V_sum // nfact
    Ta = type_a_sum // nfact
    Tb = type_b_sum // nfact

    print(f"  n={n:2d}: |V|={V_n:>8d}, T={T_n:>10d}, Ta={Ta:>8d}, Tb={Tb:>8d}, "
          f"T/|V|={T_n/V_n:.4f}, m={m}")

# ============================================================================
# PART 2: Exact self-loop orbit counts from data
# ============================================================================

print("\n  PART 2: Self-loop analysis from known data")
print("-" * 60)

# Known: transition orbits and edge counts
T_known = {3: 4, 4: 16, 5: 88, 6: 704}
E_known = {3: 1, 4: 5, 5: 30, 6: 290}

# Self-loop weight (sum over all classes)
SL_weight = {3: 12, 4: 144, 5: 1760, 6: 37920}

for n in sorted(T_known):
    m = comb(n, 2)
    nfact = factorial(n)

    T = T_known[n]
    E = E_known[n]
    SL_w = SL_weight[n]

    # Total (tournament, arc) pairs = 2^m × m
    total_pairs = (2**m) * m

    # Total orbits on (tournament, arc) pairs = T
    # These split into self-loop orbits and cross-class orbits
    # Cross-class orbits ≥ E (since each edge has ≥ 1 orbit)

    # The cross-class pairs = total - self-loop pairs = 2^m × m - SL_w
    cross_pairs = total_pairs - SL_w

    # Self-loop orbits: number of S_n-orbits on self-loop (T, e) pairs
    # The average orbit size is n! / avg_stab_size
    # For most tournaments (|Aut|=1), orbit size = n!
    # So self-loop orbits ≈ SL_w / n!

    SL_orbits_est = SL_w / nfact
    cross_orbits = T - round(SL_orbits_est)

    # Better: compute self-loop orbits exactly
    # SL orbits = (1/n!) Σ_σ |{self-loop (T,e) : σ(T)=T, σ(e)=e, [flip(T,e)]=[T]}|
    # This is hard because it requires knowing which flips are self-loops.

    # But we can compute: SL_w / n! gives a lower bound on SL orbits
    # (equals SL orbits when all self-loop tournaments have |Aut|=1)

    print(f"  n={n}: T={T}, |E|={E}, SL_w={SL_w}")
    print(f"    SL_orbits ≈ SL_w/n! = {SL_orbits_est:.2f}")
    print(f"    cross_orbits ≈ T - SL_orbits = {T - SL_orbits_est:.2f}")
    print(f"    cross_orbits / |E| = {(T - SL_orbits_est) / E:.4f}")
    print(f"    (avg orbits per edge)")

# ============================================================================
# PART 3: Look for cross_orbits / |E| pattern
# ============================================================================

print("\n  PART 3: Cross-orbits per edge")
print("-" * 60)

for n in sorted(T_known):
    nfact = factorial(n)
    T = T_known[n]
    E = E_known[n]
    SL_w = SL_weight[n]

    SL_orb = SL_w / nfact
    cross_orb = T - SL_orb
    ratio = cross_orb / E

    print(f"  n={n}: cross_orbits/|E| = {ratio:.6f}")

# n=3: (4 - 12/6) / 1 = (4-2)/1 = 2.0
# n=4: (16 - 144/24) / 5 = (16-6)/5 = 10/5 = 2.0
# n=5: (88 - 1760/120) / 30 = (88 - 14.67) / 30 = 73.33/30 = 2.444
# n=6: (704 - 37920/720) / 290 = (704 - 52.67) / 290 = 651.33/290 = 2.246

# If this ratio approached a constant c, we'd have |E| ≈ (T - SL_w/n!) / c

# ============================================================================
# PART 4: The Burnside-computable upper bound
# ============================================================================

print("\n  PART 4: Predicted edge counts using the transition formula")
print("-" * 60)

# If cross_orbits/|E| ≈ 2 (as suggested by n=3,4,6):
# |E| ≈ (T - SL_w/n!) / 2

# But we don't know SL_w for n≥7. Can we estimate SL_w?
# SL_w / total_pairs = self-loop fraction
# Self-loop fraction: 0.5, 0.375, 0.172, 0.077 for n=3,4,5,6
# These are decreasing. Let's see if there's a formula.

sl_frac = {}
for n in sorted(SL_weight):
    m = comb(n, 2)
    sl_frac[n] = SL_weight[n] / ((2**m) * m)

print(f"  Self-loop fractions: {sl_frac}")

# Sequence: 0.5, 0.375, 0.172, 0.077
# Ratios: 0.75, 0.459, 0.449
# Nearly halving each time for n≥5

# Can we compute self-loop weight via Burnside?
# SL_w = Σ_T Σ_e [flip(T,e) in [T]]
# = Σ_[T] |[T]| × |{e : flip(T,e) ~ T}|

# For a tournament T, an arc flip produces an isomorphic tournament iff
# ∃ σ ∈ S_n such that σ(flip(T,e)) = T.
# This means σ is an "almost-automorphism" that acts like an automorphism
# everywhere except at position e.

# This is related to the "partial symmetry" or "local symmetry defect" of T at e.

# ============================================================================
# PART 5: Extended transition orbit table
# ============================================================================

print("\n  PART 5: Extended transition orbit table (for formula chasing)")
print("-" * 60)

T_vals = {}
V_vals = {}
for n in range(3, 14):
    m = comb(n, 2)
    nfact = factorial(n)

    T_sum = 0
    V_sum = 0

    for part in partitions(n):
        ct = cycle_type_count(n, part)
        ft = fix_tournaments(part)
        fa = fix_arcs(part, n)
        T_sum += ct * ft * fa
        V_sum += ct * ft

    T_vals[n] = T_sum // nfact
    V_vals[n] = V_sum // nfact

print(f"  {'n':>3} {'|V|':>10} {'T_n':>12} {'T/|V|':>8} {'m':>4} {'T/m':>10}")
for n in sorted(T_vals):
    m = comb(n, 2)
    T = T_vals[n]
    V = V_vals[n]
    print(f"  {n:3d} {V:10d} {T:12d} {T/V:8.4f} {m:4d} {T/m:10.2f}")

# ============================================================================
# PART 6: T_n via Burnside for large n
# ============================================================================

print("\n  PART 6: T_n sequence and ratios")
print("-" * 60)

T_list = [T_vals[n] for n in sorted(T_vals)]
print(f"  T_n: {T_list}")
for i in range(1, len(T_list)):
    print(f"  T_{i+3}/T_{i+2} = {T_list[i]/T_list[i-1]:.4f}")

# ============================================================================
# PART 7: Can we derive |E| from T_n?
# ============================================================================

print("\n  PART 7: If |E| ≈ T/(c_n) for some sequence c_n...")
print("-" * 60)

for n in sorted(T_known):
    T = T_known[n]
    E = E_known[n]
    c = T / E
    print(f"  n={n}: T/|E| = {c:.4f}")

# c values: 4.0, 3.2, 2.933, 2.428
# Decreasing toward 2?
# This makes sense: as n grows, self-loops become rarer (self-loop fraction → 0),
# so T ≈ cross_orbits, and cross_orbits/|E| should approach some constant.

# If the average edge has weight ≈ n! (which is approximately true from Part 5 earlier),
# and each orbit has size n!/k for some k,
# then orbits per edge ≈ weight / (n!/1) ≈ avg_weight/n! ≈ 1.
# But we saw orbits per edge ≈ 2. This suggests average orbit size ≈ n!/2.

# Actually: an edge (i,j) has weight w(i,j). The number of orbits contributing
# to this edge is w(i,j) / (average orbit size).
# Average orbit size for a cross-class transition = n! / avg(|Stab|).
# For |Aut(T)|=1, |Stab(T,e)| ≤ 1, so orbit size = n!.
# But the edge weight counts transitions from BOTH sides: w(i,j) includes
# transitions from i→j and j→i. Each direction has its own orbits.
# So the number of orbits = (w_ij→ / n!_i) + (w_ji→ / n!_j)
# where n!_i = n!/|Stab| for transition from class i.

# For generic classes: w(i,j) = 2 × |C_i| × k(i,j) / total... actually this is wrong.
# w(i,j) = |C_i| × k(i,j) where k(i,j) = # arcs of T_i whose flip gives class j.
# And w(j,i) = |C_j| × k(j,i). And w(i,j) = w(j,i) (symmetric).

# The orbits from the i-side: the k(i,j) arcs of T_i that flip to class j
# fall into Aut(T_i)-orbits. Each orbit gives one S_n-orbit of transitions.
# So # orbits from i-side = # Aut(T_i)-orbits among the k(i,j) arcs.

# For |Aut|=1: every arc is its own orbit. So orbits from i-side = k(i,j).
# Total orbits for edge (i,j) = k(i,j) + k(j,i).
# And w(i,j) = |C_i| × k(i,j) = |C_j| × k(j,i).
# For generic classes: |C_i| = |C_j| = n!, so k(i,j) = k(j,i).
# Total orbits = 2 × k(i,j).

# Since w(i,j) = n! × k(i,j), and average edge weight ≈ n! (from data),
# k(i,j) ≈ 1 for most edges. So most edges have k(i,j) = 1.
# And orbits per edge = 2 × 1 = 2. This matches n=3,4!

# For n=5,6: some edges have k(i,j) > 1, pushing the average above 2.

# This is KEY: most edges correspond to a UNIQUE arc position type.
# The edges with multiple arc types are the exception.

print(f"\n  PREDICTION: most edges have exactly 1 arc-position type per side")
print(f"  → cross_orbits ≈ 2 × |E| for generic cases")
print(f"  → |E| ≈ (T_n - SL_orbits) / 2 as lower bound")

# Let's predict |E| for n=7:
T7 = T_vals[7]
# Estimate SL_orbits at n=7:
# Self-loop fraction decreases: 0.5, 0.375, 0.172, 0.077
# Roughly SL_frac ≈ 0.5^{something}
# Let's extrapolate: SL_frac(7) ≈ 0.03 (rough)
m7 = comb(7, 2)
total_pairs_7 = (2**m7) * m7
SL_w_est_7 = total_pairs_7 * 0.03
SL_orb_est_7 = SL_w_est_7 / factorial(7)
cross_orb_7 = T7 - SL_orb_est_7

print(f"\n  n=7 prediction:")
print(f"    T_7 = {T7}")
print(f"    SL_w_est ≈ {SL_w_est_7:.0f}")
print(f"    SL_orb_est ≈ {SL_orb_est_7:.1f}")
print(f"    cross_orb ≈ {cross_orb_7:.0f}")
print(f"    |E|_est (c=2) ≈ {cross_orb_7/2:.0f}")
print(f"    |E|_est (c=2.3) ≈ {cross_orb_7/2.3:.0f}")
print(f"    |E|_est (c=2.2) ≈ {cross_orb_7/2.2:.0f}")

# Also for n=8
T8 = T_vals[8]
m8 = comb(8, 2)
total_pairs_8 = (2**m8) * m8
SL_w_est_8 = total_pairs_8 * 0.01  # even smaller fraction
SL_orb_est_8 = SL_w_est_8 / factorial(8)

print(f"\n  n=8 prediction:")
print(f"    T_8 = {T8}")
print(f"    |E|_est (c=2.0) ≈ {(T8 - SL_orb_est_8)/2:.0f}")
print(f"    |E|_est (c=2.2) ≈ {(T8 - SL_orb_est_8)/2.2:.0f}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
