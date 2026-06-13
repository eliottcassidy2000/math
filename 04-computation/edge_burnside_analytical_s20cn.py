#!/usr/bin/env python3
"""
edge_burnside_analytical_s20cn.py -- kind-pasteur-2026-03-22-S20cn

ANALYTICAL COMPUTATION of T_n, D_n, and overhead structure.

T_n = (1/n!) Σ_{σ all-odd} Fix(σ) × C(f(σ), 2)
V_n = (1/n!) Σ_{σ all-odd} Fix(σ)
D_n = V_n × m - T_n

where:
  Fix(σ) = 2^{p(s)}, p(s) = Σ_{i<j} gcd(c_i,c_j) + Σ_i (c_i-1)/2
  f(σ) = number of fixed points (1-cycles)
  C(f,2) = number of fixed arcs

Key insight: ONLY odd-cycle permutations contribute.

Author: kind-pasteur-2026-03-22-S20cn
"""
import sys
from math import comb, factorial, gcd, log2
from collections import Counter
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def partitions_odd(n, max_part=None):
    """Generate all partitions of n into odd parts."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield []
        return
    start = min(max_part, n)
    if start % 2 == 0:
        start -= 1  # force odd start
    for k in range(start, 0, -2):  # odd parts only: k, k-2, k-4, ...
        for rest in partitions_odd(n - k, k):
            yield [k] + rest

def cycle_type_count(parts, n):
    """Number of permutations in S_n with given cycle type.
    = n! / (c1^a1 * a1! * c2^a2 * a2! * ...)
    where ai = multiplicity of ci."""
    mult = Counter(parts)
    denom = 1
    for c, a in mult.items():
        denom *= (c ** a) * factorial(a)
    return factorial(n) // denom

def p_value(parts):
    """p(s) = Σ_{i<j} gcd(c_i, c_j) + Σ_i (c_i - 1)/2"""
    r = len(parts)
    val = 0
    for i in range(r):
        for j in range(i+1, r):
            val += gcd(parts[i], parts[j])
    for c in parts:
        val += (c - 1) // 2
    return val

def fixed_points(parts):
    """Number of 1-cycles in the partition."""
    return parts.count(1)

print("=" * 70)
print("  ANALYTICAL T_n, D_n COMPUTATION VIA BURNSIDE")
print("  kind-pasteur-2026-03-22-S20cn")
print("=" * 70)

# Known edge counts
known_E = {3: 1, 4: 5, 5: 30, 6: 290, 7: 4086}
# Known overhead decomposition from computation
known_SL = {3: 2, 4: 6, 5: 16, 6: 58}
known_MW = {3: 0, 4: 0, 5: 12, 6: 66}

print(f"\n  {'n':>3s} {'m':>4s} {'V_n':>12s} {'T_n':>15s} {'D_n':>12s} {'D/V':>8s} {'T/2':>15s} {'|E|':>8s} {'OH':>6s} {'OH/V':>6s}")

for n in range(3, 21):
    m = comb(n, 2)
    nfact = factorial(n)

    V_sum = 0
    T_sum = 0

    for parts in partitions_odd(n):
        ct = cycle_type_count(parts, n)
        p = p_value(parts)
        f = fixed_points(parts)
        fix = 2 ** p  # = 1 << p if p < 63

        V_sum += ct * fix
        T_sum += ct * fix * comb(f, 2)

    V_n = V_sum // nfact
    T_n = T_sum // nfact
    D_n = V_n * m - T_n

    E_n = known_E.get(n, '?')
    if isinstance(E_n, int):
        OH = T_n - 2 * E_n
        oh_str = f"{OH:>6d}"
        oh_v = f"{OH/V_n:.3f}"
    else:
        oh_str = "?"
        oh_v = "?"

    if n <= 13:
        print(f"  {n:>3d} {m:>4d} {V_n:>12d} {T_n:>15d} {D_n:>12d} {D_n/V_n:>8.4f} {T_n//2:>15d} {str(E_n):>8s} {oh_str:>6s} {str(oh_v):>6s}")
    else:
        print(f"  {n:>3d} {m:>4d} {V_n:>12.4e} {T_n:>15.6e} {D_n:>12.4e} {D_n/V_n:>8.4f}")

# Detailed breakdown by cycle type
print(f"\n{'='*70}")
print(f"  CYCLE TYPE BREAKDOWN")
print(f"{'='*70}")

for n in range(3, 10):
    m = comb(n, 2)
    nfact = factorial(n)
    print(f"\n  n={n} (m={m}, n!={nfact}):")
    print(f"    {'cycle type':>20s} {'ct':>8s} {'p(s)':>6s} {'Fix':>10s} {'f':>3s} {'C(f,2)':>7s} {'V contrib':>12s} {'T contrib':>12s}")

    total_V = 0
    total_T = 0
    for parts in sorted(partitions_odd(n)):
        ct = cycle_type_count(parts, n)
        p = p_value(parts)
        f = fixed_points(parts)
        fix = 2 ** p

        V_contrib = ct * fix / nfact
        T_contrib = ct * fix * comb(f, 2) / nfact

        total_V += V_contrib
        total_T += T_contrib

        print(f"    {str(parts):>20s} {ct:>8d} {p:>6d} {fix:>10d} {f:>3d} {comb(f,2):>7d} {V_contrib:>12.4f} {T_contrib:>12.4f}")

    E = known_E.get(n)
    if E:
        print(f"    V={total_V:.0f}, T={total_T:.0f}, D={total_V*m - total_T:.0f}, E={E}, OH={total_T - 2*E:.0f}")

# Edge approximation analysis
print(f"\n{'='*70}")
print(f"  EDGE COUNT APPROXIMATION ANALYSIS")
print(f"{'='*70}")
print(f"\n  T_n/2 vs |E_n| (T/2 is upper bound for E):")
print(f"  {'n':>3s} {'T/2':>10s} {'E':>8s} {'T/2-E':>8s} {'(T/2-E)/E':>10s} {'T/(2E)':>8s}")
for n in sorted(known_E.keys()):
    m = comb(n, 2)
    nfact = factorial(n)
    V_sum = T_sum = 0
    for parts in partitions_odd(n):
        ct = cycle_type_count(parts, n)
        p = p_value(parts)
        f = fixed_points(parts)
        fix = 2 ** p
        V_sum += ct * fix
        T_sum += ct * fix * comb(f, 2)
    T_n = T_sum // nfact
    E = known_E[n]
    print(f"  {n:>3d} {T_n//2:>10d} {E:>8d} {T_n//2-E:>8d} {(T_n/2-E)/E:>10.4f} {T_n/(2*E):>8.4f}")

# Predict E_n from T_n and overhead trend
print(f"\n  OVERHEAD TRENDS:")
print(f"  {'n':>3s} {'OH':>6s} {'SL':>5s} {'MW':>5s} {'SL/OH':>7s} {'MW/OH':>7s}")
for n in sorted(known_SL.keys()):
    SL = known_SL[n]
    MW = known_MW[n]
    OH = SL + MW
    print(f"  {n:>3d} {OH:>6d} {SL:>5d} {MW:>5d} {SL/OH:>7.3f} {MW/OH:>7.3f}")

# Predict using E ~ T/2 - OH/2
print(f"\n  EDGE PREDICTIONS (from T_n and overhead extrapolation):")
for n in range(8, 14):
    m = comb(n, 2)
    nfact = factorial(n)
    V_sum = T_sum = 0
    for parts in partitions_odd(n):
        ct = cycle_type_count(parts, n)
        p = p_value(parts)
        f = fixed_points(parts)
        fix = 2 ** p
        V_sum += ct * fix
        T_sum += ct * fix * comb(f, 2)
    V_n = V_sum // nfact
    T_n = T_sum // nfact
    D_n = V_n * m - T_n
    E_approx_upper = T_n // 2
    # Extrapolate OH from trend
    # At n=7: OH/V ~ 1.62, at n=6: 2.21, at n=5: 2.33
    # OH/V decreases as ~ O(1/n) maybe
    # More precisely: OH ~ SL + MW, both ~ V * correction
    # At large n: OH/V -> 0, so E -> T/2
    print(f"  n={n}: V={V_n}, T={T_n}, T/2={T_n//2}, D={D_n}, D/V={D_n/V_n:.4f}")

print(f"\n{'='*70}")
print(f"  T_n SEQUENCE (new, not in OEIS)")
print(f"{'='*70}")
for n in range(3, 16):
    m = comb(n, 2)
    nfact = factorial(n)
    T_sum = 0
    for parts in partitions_odd(n):
        ct = cycle_type_count(parts, n)
        p = p_value(parts)
        f = fixed_points(parts)
        fix = 2 ** p
        T_sum += ct * fix * comb(f, 2)
    T_n = T_sum // nfact
    print(f"  T({n}) = {T_n}")

# Key insight: T_n / (V_n * m) -> 1 as n -> inf
# So E_n ~ V_n * m / 2 asymptotically
# The EXACT formula: E_n = (T_n - OH_n) / 2
# where OH_n = self-loop orbits + multi-weight correction
# T_n is EXACTLY computable from Burnside.
# OH_n is the remaining unknown.
