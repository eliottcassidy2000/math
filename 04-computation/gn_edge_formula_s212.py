#!/usr/bin/env python3
"""
gn_edge_formula_s212.py — Chase the edge count formula for G_n
opus-2026-03-22-S212

Known data: |E(G_n)| = 1, 5, 30, 290 for n=3,4,5,6.

Approach:
1. Compute ALL quantities derivable from Burnside theory
2. Look for relationships between |E|, |V|, self-loops, weights, etc.
3. Compute the Burnside-style "transition orbit" count and compare
4. Look for generating function patterns
5. Try to express |E| in terms of known sequences

The edge count is the KEYSTONE of the G_n description.
"""

import sys
from math import comb, factorial, gcd
from itertools import permutations
from collections import defaultdict, Counter
from functools import reduce
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  CHASING THE EDGE COUNT FORMULA FOR G_n")
print("  opus-2026-03-22-S212")
print("=" * 80)

# ============================================================================
# Known data
# ============================================================================

# From exact computation
V = {3: 2, 4: 4, 5: 12, 6: 56}  # A000568
E = {3: 1, 4: 5, 5: 30, 6: 290}
SL = {3: 1, 4: 2, 5: 7, 6: 30}  # self-loop classes
SC = {3: 2, 4: 2, 5: 8, 6: 12}  # SC classes

# Edge decomposition
E_blue = {3: 1, 4: 1, 5: 14, 6: 200}
E_black = {3: 0, 4: 4, 5: 16, 6: 90}
E_scsc = {3: 1, 4: 1, 5: 12, 6: 13}
E_nsns = {3: 0, 4: 0, 5: 2, 6: 187}
E_scns = {3: 0, 4: 4, 5: 16, 6: 90}

# Sum of degrees = 2|E|
sum_deg = {n: 2 * E[n] for n in E}

# ============================================================================
# Part 1: Burnside formula quantities
# ============================================================================

print("\n  PART 1: Burnside-derived quantities")
print("-" * 60)

def partitions(n, max_part=None):
    """Generate all partitions of n."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield []
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield [first] + rest

def cycle_type_count(n, cycle_type):
    """Number of permutations in S_n with given cycle type."""
    # cycle_type is a list of cycle lengths (partition of n)
    ct = Counter(cycle_type)
    result = factorial(n)
    for length, count in ct.items():
        result //= (length ** count) * factorial(count)
    return result

def fix_tournaments(cycle_type):
    """Number of tournaments fixed by a permutation with given cycle type.
    Non-zero only when ALL cycle lengths are odd."""
    for c in cycle_type:
        if c % 2 == 0:
            return 0

    # For all-odd cycles: Fix = 2^exponent
    # exponent = sum((c_i - 1)/2) + sum_{i<j} gcd(c_i, c_j)
    k = len(cycle_type)
    exp = sum((c - 1) // 2 for c in cycle_type)
    for i in range(k):
        for j in range(i + 1, k):
            exp += gcd(cycle_type[i], cycle_type[j])
    return 2 ** exp

def fix_arcs(cycle_type, n):
    """Number of arc positions {u,v} fixed by permutation with given cycle type.
    An arc {u,v} is fixed iff σ(u)=u,σ(v)=v or σ(u)=v,σ(v)=u.
    = C(f,2) + t where f = fixed points, t = 2-cycles."""
    f = cycle_type.count(1)  # fixed points
    t = cycle_type.count(2)  # 2-cycles
    return comb(f, 2) + t

for n in range(3, 10):
    m = comb(n, 2)
    nfact = factorial(n)

    # Burnside vertex count
    vertex_sum = 0
    # Burnside transition orbit count
    transition_sum = 0
    # Also: compute Fix_tournaments for only-odd-cycle types
    all_odd_sum = 0

    for part in partitions(n):
        ct_count = cycle_type_count(n, part)
        fix_t = fix_tournaments(part)
        fix_a = fix_arcs(part, n)

        vertex_sum += ct_count * fix_t
        transition_sum += ct_count * fix_t * fix_a

        if all(c % 2 == 1 for c in part):
            all_odd_sum += ct_count * fix_t

    V_n = vertex_sum // nfact
    transition_orbits = transition_sum // nfact

    print(f"  n={n}: |V|={V_n}, m={m}, "
          f"transition_orbits={transition_orbits}, "
          f"total_transitions={2**m * m}")

    if n <= 6:
        print(f"         KNOWN: |E|={E[n]}, 2|E|={2*E[n]}, "
              f"self-loop classes={SL[n]}")
        print(f"         transition_orbits - |E| = {transition_orbits - E[n]}")
        print(f"         These are: self-loop orbits + multi-weight edge orbits")

# ============================================================================
# Part 2: Relationship hunting
# ============================================================================

print("\n  PART 2: Relationship hunting")
print("-" * 60)

for n in range(3, 7):
    m = comb(n, 2)
    nfact = factorial(n)
    v = V[n]
    e = E[n]
    sl = SL[n]

    total_weight = 2**m * m
    total_self_weight_exact = total_weight - sum_deg[n]  # wrong: sum_deg counts cross-class only

    avg_deg = 2 * e / v
    density = 2 * e / (v * (v - 1))

    print(f"\n  n={n}: |V|={v}, |E|={e}, m={m}")
    print(f"    2^m = {2**m}")
    print(f"    2^m × m = {total_weight}")
    print(f"    2|E| = {2*e}")
    print(f"    avg_deg = {avg_deg:.4f}")
    print(f"    avg_deg / m = {avg_deg/m:.4f}")
    print(f"    density = {density:.4f}")
    print(f"    |E| / |V| = {e/v:.4f}")
    print(f"    |E| / |V|^2 = {e/v**2:.6f}")
    print(f"    |E| / m = {e/m:.4f}")
    print(f"    |E| / (|V| × m) = {e/(v*m):.6f}")
    print(f"    |E| × n! = {e * nfact}")
    print(f"    |E| × n! / 2^m = {e * nfact / 2**m:.4f}")
    print(f"    |E| × n!^2 / 2^m = {e * nfact**2 / 2**m:.4f}")

# ============================================================================
# Part 3: Try specific formulas
# ============================================================================

print("\n  PART 3: Testing specific formula candidates")
print("-" * 60)

for n in range(3, 7):
    m = comb(n, 2)
    v = V[n]
    e = E[n]
    nfact = factorial(n)

    candidates = {
        'V*(V-1)/2': v * (v-1) // 2,
        'V*m/2': v * m / 2,
        'V*(m-1)/2': v * (m-1) / 2,
        '(2^m * m - S) / (2*n!)': None,  # need self-loop total
        'V^2 * m / (2 * n!)': v**2 * m / (2 * nfact),
        '2^m * m / (n! * 2)': 2**m * m / (nfact * 2),
        '2^(m-1) * m / n!': 2**(m-1) * m / nfact,
    }

    print(f"\n  n={n}: actual |E|={e}")
    for name, val in candidates.items():
        if val is not None:
            print(f"    {name:30s} = {val:12.4f}  {'✓' if abs(val - e) < 0.001 else ''}")

# ============================================================================
# Part 4: OEIS lookup prep
# ============================================================================

print("\n  PART 4: Sequences for OEIS lookup")
print("-" * 60)

print(f"  Edge count: {[E[n] for n in sorted(E)]}")
print(f"  Double edge count: {[2*E[n] for n in sorted(E)]}")
print(f"  Edge × n!: {[E[n]*factorial(n) for n in sorted(E)]}")
print(f"  Edge × n!^2: {[E[n]*factorial(n)**2 for n in sorted(E)]}")
print(f"  Edge × 2^C(n,2): {[E[n]*(2**comb(n,2)) for n in sorted(E)]}")

# Also compute some derived quantities
for n in sorted(E):
    m = comb(n, 2)
    print(f"  n={n}: |E|×n! = {E[n]*factorial(n)}, "
          f"|E|×2^m = {E[n]*(2**m)}, "
          f"|E|/(V×m) = {E[n]/(V[n]*m):.6f}")

# ============================================================================
# Part 5: The total weight approach
# ============================================================================

print("\n  PART 5: Total weight decomposition")
print("-" * 60)

# For each n, we know:
# total_weight = 2^m × m (all tournament-arc pairs)
# = self_loop_weight + cross_weight
# cross_weight = Σ_{edges} w(i,j)  (bidirectional: w(i,j) = w(j,i))
# = 2 × Σ_{edges} w_one_dir(i,j)

# The number of distinct edge weights gives a clue about the edge count.
# If all edges had the same weight w, then:
# |E| × 2w = cross_weight = total - self_loop_weight

# From the data:
self_loop_total = {
    3: 12,     # class 0: 12
    4: 144,    # classes 0,3: 72 each
    5: 1760,
    6: None    # need to sum from full data
}

# Sum self-loops from n=6 data (from S211 output):
# class 0: 3600, 1: 480, 2: 240, 3: 240, 4: 480, 5: 2880, 6: 2160, 7: 2880,
# 8: 1440, 11: 1440, 13: 1440, 14: 1440, 19: 720, 20: 720, 23: 2160,
# 25: 720, 26: 720, 29: 2160, 30: 2160, 31: 720, 32: 720, 35: 720,
# 41: 720, 42: 240, 44: 240, 46: 2160, 47: 720, 50: 720, 51: 2160, 54: 720
sl_n6 = (3600 + 480 + 240 + 240 + 480 + 2880 + 2160 + 2880 +
         1440 + 1440 + 1440 + 1440 + 720 + 720 + 2160 +
         720 + 720 + 2160 + 2160 + 720 + 720 + 720 +
         720 + 240 + 240 + 2160 + 720 + 720 + 2160 + 720)
self_loop_total[6] = sl_n6

# From n=5 data:
sl_n5 = 480 + 40 + 40 + 360 + 360 + 240 + 240
self_loop_total[5] = sl_n5

print(f"  Self-loop weight:")
for n in sorted(self_loop_total):
    if self_loop_total[n] is not None:
        m = comb(n, 2)
        tw = 2**m * m
        cw = tw - self_loop_total[n]
        avg_edge_weight = cw / (2 * E[n]) if n in E else 0
        print(f"    n={n}: total={tw}, self={self_loop_total[n]}, "
              f"cross={cw}, avg_edge_wt={avg_edge_weight:.2f}")

# ============================================================================
# Part 6: The average degree formula
# ============================================================================

print("\n  PART 6: Average degree analysis")
print("-" * 60)

for n in range(3, 7):
    m = comb(n, 2)
    v = V[n]
    e = E[n]
    nfact = factorial(n)
    avg_deg = 2 * e / v

    print(f"  n={n}: avg_deg = {avg_deg:.4f}")
    print(f"    m = {m}")
    print(f"    avg_deg / (m-1) = {avg_deg/(m-1):.4f}")
    print(f"    avg_deg / m = {avg_deg/m:.4f}")
    print(f"    1 - avg_deg/m = {1 - avg_deg/m:.4f}")
    # The "missing fraction" = 1 - avg_deg/m tells us how many arcs
    # produce same-class neighbors (on average)

    # Self-loop fraction of total
    if self_loop_total[n] is not None:
        tw = 2**m * m
        sf = self_loop_total[n] / tw
        print(f"    self-loop fraction = {sf:.4f}")
        # avg arcs producing self-loop per tournament = m * sf
        print(f"    avg self-loop arcs = {m * sf:.2f}")
        # avg arcs producing distinct classes = avg_deg? No, that's wrong.
        # A tournament has m arcs. Each produces a neighbor class.
        # Some arcs produce the same class (self-loop arcs) → self_loop_frac
        # Remaining arcs produce neighbor classes (possibly repeated)
        # deg(class_i) = # distinct neighbor classes

# ============================================================================
# Part 7: The transition orbit formula (Burnside for edges)
# ============================================================================

print("\n  PART 7: Burnside transition orbit count vs edge count")
print("-" * 60)

for n in range(3, 10):
    m = comb(n, 2)
    nfact = factorial(n)

    # Compute: T = (1/n!) Σ_{σ} Fix_T(σ) × Fix_arcs(σ)
    # This counts orbits of S_n on (tournament, arc-position) pairs

    total_sum = 0
    for part in partitions(n):
        ct = cycle_type_count(n, part)
        ft = fix_tournaments(part)
        fa = fix_arcs(part, n)
        total_sum += ct * ft * fa

    T_n = total_sum // nfact

    # Also compute just the "all-arcs" Burnside for the arc-flip graph
    # orbits of S_n on tournament-arc pairs
    # = orbits on vertices of Q_m (= |V(G_n)|) × orbits on edges of Q_m?
    # No, it's orbits on (vertex, dimension) pairs of Q_m.

    if n <= 6:
        print(f"  n={n}: transition_orbits = {T_n}, "
              f"|E| = {E[n]}, |V| = {V[n]}, "
              f"T - |E| - |V| = {T_n - E[n] - V[n]}")
        # T counts all orbits including self-loops.
        # Self-loop orbits + cross-class orbits = T
        # But multiple orbits can give the same edge... hmm
    else:
        print(f"  n={n}: transition_orbits = {T_n}, |V| = {V.get(n, '?')}")

# ============================================================================
# Part 8: Exact formula attempt via the structure theorem
# ============================================================================

print("\n  PART 8: Structure theorem approach")
print("-" * 60)

# Key identity: for each iso class i with representative T_i,
# m arcs of T_i produce some set of neighbor classes.
# deg(i) = |{neighbor classes}| (excluding self)
#
# For |Aut(T_i)| = 1 (trivial automorphism):
#   The m arcs are all "different" (no automorphism relates them)
#   But some arcs can still produce the same neighbor class
#   deg(i) ≤ m
#
# For |Aut(T_i)| > 1:
#   Arcs in the same Aut-orbit MUST produce the same neighbor
#   Number of Aut-orbits on arcs ≤ m / |Aut|
#   deg(i) ≤ # Aut-orbits on arcs

# At n=5: m=10
# Generic (|Aut|=1) classes: 4 classes, degs = 6,6,7,7. Avg = 6.5
# SC classes with |Aut|=1: 4 classes, degs = 6,7,6,6. Avg = 6.25
# |Aut|=3 classes: 3 classes, degs = 3,4,3. Avg = 3.33
# |Aut|=5 class: 1 class, deg = 2

# The ratio deg/m:
# |Aut|=1: 6.5/10 = 0.65
# |Aut|=3: 3.33/10 = 0.33
# |Aut|=5: 2/10 = 0.20

# At n=6: m=15
# |Aut|=1: most have degs 10-14, avg ≈ 12.4, ratio ≈ 0.83
# |Aut|=3: degs 5-7, avg ≈ 5.9, ratio ≈ 0.39
# |Aut|=5: deg 3, ratio ≈ 0.20
# |Aut|=9: deg 3, ratio ≈ 0.20

# Pattern: ratio ≈ 1 - 1/|Aut| for |Aut|=1,3 but not for 5,9

# Let me compute the Burnside orbit count for arcs per class
for n in range(3, 7):
    m = comb(n, 2)
    print(f"\n  n={n}:")
    print(f"    |Aut|  deg  orbits  deg/m  deg/orbits")

    # Compute for each possible |Aut| size
    # For a class with |Aut|=a, the number of Aut-orbits on arcs depends
    # on the specific automorphism group structure, not just |Aut|.
    # But as an estimate: orbits ≈ m / a (if Aut acts freely) or ≈ m (if Aut is trivial)

    # Actually we can compute exactly from the data at n=5,6
    # Using |Aut| and the known degree from S211

# Let me just print the key ratios
print("\n  RATIO TABLE: deg(i) / C(n,2) for each n")
known_degs = {
    3: [(1, 1, 1), (3, 1, 3)],
    4: [(1, 3, 1), (3, 2, 3), (3, 2, 3), (1, 3, 5)],
    5: [(1, 6, 1), (3, 3, 3), (3, 4, 3), (3, 3, 3),
        (1, 7, 5), (1, 7, 5),
        (1, 6, 9), (1, 7, 9), (1, 6, 11), (1, 6, 13),
        (3, 3, 15), (5, 2, 15)],
}
# Format: (|Aut|, deg, H)

for n in sorted(known_degs):
    m = comb(n, 2)
    print(f"\n  n={n}, m={m}:")
    for aut, deg, h in known_degs[n]:
        orbits_est = m  # for |Aut|=1
        if aut > 1:
            orbits_est = m // aut + (m % aut > 0)  # rough upper bound
        print(f"    H={h:3d} |Aut|={aut} deg={deg:2d} deg/m={deg/m:.3f} "
              f"m/|Aut|={m/aut:.1f}")

# ============================================================================
# Part 9: Generating function / recurrence search
# ============================================================================

print("\n  PART 9: Generating function / recurrence search")
print("-" * 60)

e_seq = [1, 5, 30, 290]  # for n=3,4,5,6

# Test: a(n) = c * a(n-1) + d * a(n-2)?
# 5 = c*1 + d*?  (no a(1))
# Let's test a(n)/a(n-1):
ratios = [e_seq[i+1]/e_seq[i] for i in range(len(e_seq)-1)]
print(f"  Ratios a(n+1)/a(n): {ratios}")
# [5.0, 6.0, 9.667]

# Test: a(n) = f(n) * a(n-1)?
# f(4) = 5, f(5) = 6, f(6) = 290/30 = 9.667
# Not a clean function of n

# Test second differences
diffs = [e_seq[i+1] - e_seq[i] for i in range(len(e_seq)-1)]
print(f"  First differences: {diffs}")
# [4, 25, 260]
diffs2 = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
print(f"  Second differences: {diffs2}")
# [21, 235]

# Test: |E| = polynomial in n?
# n=3: 1, n=4: 5, n=5: 30, n=6: 290
# Degree 1: 1 = 3a+b, 5 = 4a+b → a=4, b=-11 → 5*5-11=14≠30. No.
# Degree 2: 1=9a+3b+c, 5=16a+4b+c, 30=25a+5b+c
#   5-1=7a+b → 7a+b=4
#   30-5=9a+b → 9a+b=25
#   → 2a=21, a=10.5. Not integer. No clean degree-2 polynomial.

# Test: |E| = a * V^2 + b * V + c?
# V=2: 1 = 4a+2b+c
# V=4: 5 = 16a+4b+c
# V=12: 30 = 144a+12b+c
# V=56: 290 = 3136a+56b+c
# System: 12a+2b=4, 128a+8b=25, 2992a+42b=256
# From first: 6a+b=2, b=2-6a
# Second: 128a+8(2-6a)=25 → 128a+16-48a=25 → 80a=9 → a=9/80=0.1125
# Not clean.

# Test: |E| ~ |V|^alpha * m^beta?
import math
print(f"\n  Regression: log(|E|) vs log(|V|), log(m)")
for i, n in enumerate([3,4,5,6]):
    le = math.log(E[n])
    lv = math.log(V[n])
    lm = math.log(comb(n,2))
    print(f"    n={n}: log|E|={le:.4f}, log|V|={lv:.4f}, log(m)={lm:.4f}")

# ============================================================================
# Part 10: Sum of degree formula
# ============================================================================

print("\n  PART 10: Total degree = sum over classes of deg(i)")
print("-" * 60)

# We know: for each class i, deg(i) ≤ m - self_arcs(i)/|C_i|
# where self_arcs(i) = number of arc flips in class i that stay in class i.

# And: Σ_i deg(i) × |C_i| ≤ 2^m × m - Σ_i self_arcs(i)
# This is because each tournament contributes m neighbors, and deg(i) counts
# DISTINCT neighbor classes per class, not per tournament.

# The weighted sum: Σ_i deg(i) × |C_i| is NOT 2^m × m - self_loop.
# Instead: Σ_i deg(i) is just 2|E|.

# But we can relate the weighted degree:
# w_deg(i) = sum of weights of edges incident to i = Σ_{j≠i} w(i,j)
# = |C_i| × m - self_arcs(i)
# And Σ_i w_deg(i) = 2 × (total cross-class weight) = 2 × (2^m × m - total_self_loop)

for n in range(3, 7):
    m = comb(n, 2)
    v = V[n]
    e = E[n]
    tw = 2**m * m
    sl_w = self_loop_total[n] if self_loop_total[n] is not None else 0
    cross_w = tw - sl_w

    # avg weight per edge (each edge counted once from each direction = weight/2)
    avg_w = cross_w / (2 * e) if e > 0 else 0

    # avg tournament-neighbors per class
    avg_w_deg = cross_w / v

    print(f"  n={n}: total_wt={tw}, self_wt={sl_w}, cross_wt={cross_w}")
    print(f"    2|E|={2*e}, avg_deg={2*e/v:.2f}")
    print(f"    cross_wt / (2|E|) = avg_edge_weight = {avg_w:.2f}")
    print(f"    cross_wt / |V| = avg_weighted_deg = {avg_w_deg:.2f}")
    print(f"    avg_edge_weight / (n!/1) = {avg_w / factorial(n):.4f}")
    print(f"    avg_edge_weight / n! = {avg_w / factorial(n):.4f}")

# ============================================================================
# Part 11: RATIO patterns
# ============================================================================

print("\n  PART 11: Master ratio table")
print("-" * 60)

print(f"  {'n':>3} {'|V|':>6} {'|E|':>6} {'m':>4} {'n!':>6} {'2^m':>8} "
      f"{'|E|/|V|':>8} {'density':>8} {'|E|×n!':>10} {'|E|×n!/2^m':>12}")

for n in range(3, 7):
    m = comb(n, 2)
    v = V[n]
    e = E[n]
    nf = factorial(n)
    two_m = 2**m
    density = 2*e/(v*(v-1)) if v > 1 else 0

    print(f"  {n:3d} {v:6d} {e:6d} {m:4d} {nf:6d} {two_m:8d} "
          f"{e/v:8.3f} {density:8.5f} {e*nf:10d} {e*nf/two_m:12.4f}")

# ============================================================================
# Part 12: The critical quantity: |E| × n! / 2^m
# ============================================================================

print("\n  PART 12: |E| × n! / 2^m")
print("-" * 60)

for n in range(3, 7):
    m = comb(n, 2)
    nf = factorial(n)
    val = E[n] * nf / (2**m)
    print(f"  n={n}: |E|×n!/2^m = {val:.4f}")
    # n=3: 1×6/8 = 0.75
    # n=4: 5×24/64 = 1.875
    # n=5: 30×120/1024 = 3.515625
    # n=6: 290×720/32768 = 6.3720703125

# These are: 3/4, 15/8, 225/64, ... ?
# 0.75 = 3/4
# 1.875 = 15/8
# 3.515625 = 225/64 = (15/8)^2 / (15/8) ... no = 225/64
# 6.3720703125 = ?

# Let me check: 6.3720703125 * 32768 = 208800 = 290 × 720. Yes.
# 208800 = 2^5 × 3^2 × 5^2 × 29
# Hmm, 29 is a prime. That's ugly.

# Let me try: 3/4, 15/8, 225/64 = (15)^2/64
# 3/4 = 3/4
# 15/8 = 15/8
# 225/64: 30×120/1024 = 3600/1024 = 225/64. Yes!
# 225 = 15^2. 64 = 8^2. So 225/64 = (15/8)^2.
# Pattern: 3/4, 15/8, (15/8)^2?
# Check: (3/4) × (15/8)^0 = 3/4 = 0.75 ✓
# (3/4) × (15/8)^1 = 45/32 = 1.40625 ≠ 1.875. No.

# Different pattern: ratios
# 1.875/0.75 = 2.5
# 3.515625/1.875 = 1.875
# 6.372.../3.515625 = 1.8125

# Not clean.

# Let me try |E| * n!^2 / 2^m
for n in range(3, 7):
    m = comb(n, 2)
    nf = factorial(n)
    val = E[n] * nf**2 / (2**m)
    print(f"  n={n}: |E|×n!²/2^m = {val:.2f}")

# n=3: 1×36/8 = 4.5
# n=4: 5×576/64 = 45
# n=5: 30×14400/1024 = 421.875
# n=6: 290×518400/32768 = 4587.890625

# 4.5, 45, 421.875, 4587.890625
# ratios: 10, 9.375, 10.886... Not clean.

print(f"\n  PART 12b: |E| * 2 / (|V| * (|V|-1))  [density]")
for n in range(3, 7):
    v = V[n]
    e = E[n]
    d = 2*e / (v*(v-1))
    print(f"  n={n}: density = {d:.6f}, 1/density = {1/d:.4f}")
    # n=3: 1.0, n=4: 0.833, n=5: 0.4545, n=6: 0.1883
    # 1/density: 1, 1.2, 2.2, 5.31

# ============================================================================
# Part 13: Self-loop weight Burnside formula
# ============================================================================

print("\n  PART 13: Self-loop weight via Burnside")
print("-" * 60)

# A self-loop at class [T] occurs when flipping an arc e gives a tournament
# isomorphic to T (in the same class). The self-loop weight for class [T] is:
# sl([T]) = Σ_{T' in [T]} |{e : flip(T',e) in [T]}|
#         = |[T]| × |{e : flip(T,e) in [T]}|   (by symmetry)

# By Burnside: the self-loop weight across all classes is
# Σ_i sl(i) = total_self_loop_weight

# Can we compute this via Burnside?
# An arc flip (T, e) is a self-loop iff [flip(T,e)] = [T],
# i.e., ∃ σ ∈ S_n such that σ(flip(T,e)) = T.

# This means T and flip(T,e) are isomorphic. Let σ be the isomorphism.
# Then σ(T)(a,b) = T(a,b) for all (a,b) ≠ e, and σ(T)(e) = 1-T(e).
# So σ is an "almost-automorphism" that flips exactly one arc.

# This is hard to count via Burnside directly.

# But the total self-loop weight can be related to:
# total_self_loop = Σ_{T} |{e : [flip(T,e)] = [T]}|
# = Σ_{T} |{e : ∃σ∈S_n, σ(flip(T,e)) = T}|

# For a tournament T with |Aut(T)| = a:
# |{e : flip(T,e) in [T]}| depends on the specific tournament.

# Actually, by orbit-stabilizer:
# |{T' in [T] : d(T, T') = 1}| = |[T]| × (avg arcs that stay in class)
# where the average is over all T in [T].

# The total self-loop weight is:
# Σ_classes Σ_{T in class} |{arcs e : flip(T,e) in class}|
# = Σ_T |{e : [flip(T,e)] = [T]}|

# This is the number of (T, e) pairs where the flip preserves the iso class.
# By Burnside, this is n! times the number of orbits of S_n on
# {(T, e) : [flip(T,e)] = [T]}.

# Hmm, but this is circular.

# Let me just compute the self-loop weight sequence and check OEIS
print(f"  Self-loop weight sequence: {[self_loop_total[n] for n in [3,4,5,6]]}")
# 12, 144, 1760, ?
# 12 = 2^2 × 3
# 144 = 2^4 × 3^2 = 12^2
# 1760 = 2^5 × 5 × 11
# Hmm, no clean pattern

print(f"\n  Self-loop fraction:")
for n in [3,4,5,6]:
    m = comb(n, 2)
    tw = 2**m * m
    sf = self_loop_total[n] / tw
    print(f"  n={n}: self_loop/total = {self_loop_total[n]}/{tw} = {sf:.6f}")

print(f"\n{'='*80}")
print(f"  DONE — awaiting n=7 result to extend all sequences")
print(f"{'='*80}")
