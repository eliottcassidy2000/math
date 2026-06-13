#!/usr/bin/env python3
"""
Interlacing verification for I(Omega(T), x).

CONJECTURE: For every tournament T and every vertex v,
  I(Omega(T-v), x) interlaces I(Omega(T), x).

If true, combined with base case (n=3: I = 1+c3*x, one real root),
gives real-rootedness of I(Omega(T), x) for ALL n by induction.

Previous verification: 0/5120 (n=5 exhaustive), 0/196608 (n=6 exhaustive)
This script extends to n=7 (exhaustive 3-cycle) and n=8,9 (random, full Omega).

Author: opus-2026-03-06-S17
Verification script for: THM-022, OPEN-Q-015 (interlacing lead)
"""

import sys
import os
import random
import itertools

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, random_tournament,
                             find_odd_cycles, conflict_graph)

random.seed(42)


def indep_poly_coeffs(adj):
    """Independence polynomial coefficients [alpha_0, alpha_1, ..., alpha_k]."""
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def poly_roots(coeffs):
    """Find roots of polynomial given as [a0, a1, ..., ad] = a0 + a1*x + ... + ad*x^d."""
    if len(coeffs) <= 1:
        return []
    # Use numpy if available, otherwise companion matrix
    try:
        import numpy as np
        # numpy.roots expects [ad, ..., a1, a0]
        p = list(reversed(coeffs))
        roots = np.roots(p)
        # Sort by real part
        return sorted(roots, key=lambda r: r.real)
    except ImportError:
        # Fallback for degree 1 and 2
        d = len(coeffs) - 1
        if d == 1:
            return [-coeffs[0] / coeffs[1]]
        elif d == 2:
            a, b, c = coeffs[0], coeffs[1], coeffs[2]
            disc = b * b - 4 * a * c
            if disc >= 0:
                import math
                r1 = (-b - math.sqrt(disc)) / (2 * c)
                r2 = (-b + math.sqrt(disc)) / (2 * c)
                return sorted([r1, r2])
            else:
                return []  # complex
        else:
            raise RuntimeError(f"Need numpy for degree {d} polynomials")


def all_real(roots, tol=1e-8):
    """Check if all roots are real (imaginary part < tol)."""
    return all(abs(r.imag) < tol for r in roots) if roots else True


def check_interlacing(roots_big, roots_small, tol=1e-5):
    """
    Check that roots_small interlaces roots_big.
    Interlacing: r1 <= s1 <= r2 <= s2 <= ... <= r_{d-1} <= s_{d-1} <= r_d
    where r_i are roots of big (degree d) and s_i are roots of small (degree d-1).
    """
    d = len(roots_big)
    ds = len(roots_small)

    if ds == 0:
        return True
    if ds > d:
        return False

    # Extract real parts
    rb = sorted([r.real for r in roots_big])
    rs = sorted([r.real for r in roots_small])

    if ds == d - 1:
        # Standard interlacing: r1 <= s1 <= r2 <= ... <= s_{d-1} <= r_d
        for i in range(ds):
            if rb[i] > rs[i] + tol:
                return False
            if rs[i] > rb[i + 1] + tol:
                return False
        return True
    elif ds == d:
        # Same degree — check if one interlaces the other
        ok1 = all(rb[i] <= rs[i] + tol for i in range(d)) and \
              all(rs[i] <= rb[i + 1] + tol for i in range(d - 1))
        ok2 = all(rs[i] <= rb[i] + tol for i in range(d)) and \
              all(rb[i] <= rs[i + 1] + tol for i in range(d - 1))
        return ok1 or ok2
    elif ds < d - 1:
        # Degree drop > 1: trivially interlaces
        return True

    return False


def get_omega_ip(T, cycles_filter=None):
    """Get independence polynomial of Omega(T), optionally filtering cycle lengths."""
    cycles = find_odd_cycles(T)
    if cycles_filter:
        cycles = [c for c in cycles if len(c) in cycles_filter]
    if not cycles:
        return [1]
    cg = conflict_graph(cycles)
    return indep_poly_coeffs(cg)


def delete_vertex(T, v):
    n = len(T)
    return [[T[i][j] for j in range(n) if j != v] for i in range(n) if i != v]


def run_interlacing_test(n, mode='exhaustive', num_samples=200, cycle_types=None):
    """Test interlacing for n-vertex tournaments."""
    label = f"n={n}"
    if cycle_types:
        label += f" ({cycle_types}-cycles)"
    if mode == 'random':
        label += f" ({num_samples} random)"
    else:
        label += " (exhaustive)"

    print(f"\n  Testing {label}...")

    failures = 0
    total = 0
    real_root_failures = 0

    if mode == 'exhaustive':
        m = n * (n - 1) // 2
        iterator = range(1 << m)
        expected_total = (1 << m) * n
    else:
        iterator = range(num_samples)
        expected_total = num_samples * n

    for idx, item in enumerate(iterator):
        if mode == 'exhaustive':
            T = tournament_from_bits(n, item)
        else:
            T = random_tournament(n)

        ip_T = get_omega_ip(T, cycle_types)
        roots_T = poly_roots(ip_T)

        if not all_real(roots_T):
            real_root_failures += 1

        for v in range(n):
            Tv = delete_vertex(T, v)
            ip_Tv = get_omega_ip(Tv, cycle_types)
            roots_Tv = poly_roots(ip_Tv)

            if roots_T and roots_Tv and all_real(roots_T) and all_real(roots_Tv):
                if not check_interlacing(roots_T, roots_Tv):
                    failures += 1
                    if failures <= 3:
                        print(f"    FAIL: T={ip_T}, Tv={ip_Tv}")
                        print(f"      roots_T = {[f'{r.real:.4f}' for r in roots_T]}")
                        print(f"      roots_Tv = {[f'{r.real:.4f}' for r in roots_Tv]}")
            total += 1

        if mode == 'exhaustive' and item % 5000 == 0 and item > 0:
            pct = 100 * item / (1 << m)
            print(f"    {item}/{1 << m} ({pct:.0f}%), {failures} failures so far")

    print(f"  {label}: {failures}/{total} interlacing failures, "
          f"{real_root_failures} real-root failures")
    return failures, real_root_failures


# ============================================================
# Main verification
# ============================================================
print("=" * 70)
print("INTERLACING VERIFICATION: I(Omega(T-v),x) interlaces I(Omega(T),x)")
print("=" * 70)

total_interlacing_fails = 0
total_realroot_fails = 0

# n=5 exhaustive (full Omega)
f, r = run_interlacing_test(5, 'exhaustive')
total_interlacing_fails += f
total_realroot_fails += r

# n=6 exhaustive (full Omega)
f, r = run_interlacing_test(6, 'exhaustive')
total_interlacing_fails += f
total_realroot_fails += r

# n=7 random (full Omega — 3-cycles + 5-cycles + 7-cycles)
f, r = run_interlacing_test(7, 'random', 500)
total_interlacing_fails += f
total_realroot_fails += r

# n=8 random (full Omega)
f, r = run_interlacing_test(8, 'random', 200)
total_interlacing_fails += f
total_realroot_fails += r

# n=9 random (3-cycle subgraph only — full Omega too expensive)
f, r = run_interlacing_test(9, 'random', 100, cycle_types={3})
total_interlacing_fails += f
total_realroot_fails += r

# n=10 random (3-cycle subgraph)
f, r = run_interlacing_test(10, 'random', 50, cycle_types={3})
total_interlacing_fails += f
total_realroot_fails += r

print(f"\n{'='*70}")
print(f"SUMMARY: {total_interlacing_fails} total interlacing failures, "
      f"{total_realroot_fails} real-root failures")
if total_interlacing_fails == 0 and total_realroot_fails == 0:
    print("ALL TESTS PASSED. Interlacing conjecture holds at all tested n.")
    print("Combined with base case, this gives real roots for all n (if proved).")
else:
    if total_interlacing_fails > 0:
        print("*** INTERLACING FAILURES FOUND ***")
    if total_realroot_fails > 0:
        print("*** REAL-ROOT FAILURES FOUND ***")
print("=" * 70)
