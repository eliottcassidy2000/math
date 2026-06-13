#!/usr/bin/env python3
"""
D2_determines_H.py -- kind-pasteur-2026-03-13-S60

KEY QUESTION: For circulant tournaments on Z_p, does the D^2 multiset
(= unordered eigenvalue moduli) determine H?

If yes: this is a deep structural theorem about permanent and spectrum.
If no: there exist orientations with the same eigenvalue moduli but different H.

The "S10 cancellation at p=11" would then be a trivial consequence of
having exactly n_types = n_params (4 orbit types, 4 linear params).

This script:
1. For each orientation, compute D^2 multiset AND H
2. Group by D^2 multiset
3. Check if each group has unique H
4. Count orbit types and verify the dimension argument

ALSO: test for NON-LINEAR moment dependence -- does H = polynomial(S4,...)?
"""

import cmath
import numpy as np
from collections import defaultdict
import time


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def compute_H_heldkarp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def compute_D2_multiset(p, S):
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D2_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D2_vals.append(lam.imag ** 2)
    return tuple(sorted(round(d, 8) for d in D2_vals))


def compute_moments(p, S, max_moment=None):
    m = (p - 1) // 2
    if max_moment is None:
        max_moment = 2 * m
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)
    moments = {}
    for k in range(2, max_moment + 1, 2):
        moments[k] = sum(d**k for d in D_vals)
    return moments


def test_D2_determines_H(p):
    """Core test: does D^2 multiset determine H?"""
    m = (p - 1) // 2
    N = 1 << m
    limit = min(N, 128)

    print(f"\n{'='*70}")
    print(f"D^2 DETERMINES H? at p={p}, m={m}, testing {limit} orientations")
    print(f"{'='*70}")

    t0 = time.time()
    by_D2 = defaultdict(list)

    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        D2 = compute_D2_multiset(p, S)
        A = build_adj(p, S)
        H = compute_H_heldkarp(A, p)

        by_D2[D2].append({'bits': bits, 'H': H, 'S': S})

    t1 = time.time()
    print(f"  Computed in {t1-t0:.1f}s")

    n_types = len(by_D2)
    print(f"  {n_types} unique D^2 multisets")

    # Check uniqueness
    all_unique = True
    for D2, entries in by_D2.items():
        H_values = set(e['H'] for e in entries)
        if len(H_values) > 1:
            all_unique = False
            print(f"  *** FAILURE: D^2={D2} has H values {H_values} ***")

    if all_unique:
        print(f"  *** D^2 multiset DETERMINES H at p={p} ***")

    # Show the orbit types
    print(f"\n  Orbit types:")
    for D2 in sorted(by_D2.keys()):
        entries = by_D2[D2]
        H = entries[0]['H']
        print(f"    D^2={[round(d,4) for d in D2]}: H={H}, count={len(entries)}")

    # Dimension analysis
    print(f"\n  Dimension analysis:")
    print(f"    m = {m} (number of eigenvalue pairs)")
    print(f"    m-1 = {m-1} (dimension of D^2 hyperplane, since S2=const)")
    print(f"    n_types = {n_types} (number of unique D^2 multisets)")
    print(f"    Moments needed for LINEAR fit: n_types - 1 = {n_types - 1}")
    print(f"    Available moment parameters: m-1 = {m-1}")

    if n_types - 1 <= m - 1:
        print(f"    FIT IS TRIVIALLY EXACT (n_types-1 <= m-1)")
        print(f"    Any {n_types-1} independent moments + const will fit H exactly")
    else:
        print(f"    FIT IS GENUINELY OVERCONSTRAINED (n_types-1 > m-1)")
        print(f"    Cancellations would be a genuine algebraic identity")

    return by_D2, n_types


def count_orbit_types_no_H(p, limit=None):
    """Count orbit types without computing H (fast, works for larger p)."""
    m = (p - 1) // 2
    N = 1 << m
    if limit is None:
        limit = min(N, 256)

    D2_set = set()
    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        D2 = compute_D2_multiset(p, S)
        D2_set.add(D2)

    return len(D2_set)


def orbit_count_pattern():
    """Determine the number of orbit types for various primes."""
    print(f"\n{'='*70}")
    print(f"ORBIT TYPE COUNTS")
    print(f"{'='*70}")

    print(f"\n  {'p':>4} {'m':>3} {'2^m':>6} {'n_types':>8} {'m-1':>4} {'p mod 4':>7} "
          f"{'Trivial?':>8}")
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
        m = (p - 1) // 2
        N = 1 << m
        limit = min(N, 512)
        n_types = count_orbit_types_no_H(p, limit)
        trivial = "YES" if n_types - 1 <= m - 1 else "NO"
        print(f"  {p:4d} {m:3d} {N:6d} {n_types:8d} {m-1:4d} {p%4:7d} {trivial:>8}")

    # For larger primes, we can only sample
    print(f"\n  Sampling larger primes (128 orientations):")
    for p in [37, 41, 43, 47, 53, 59, 61]:
        m = (p - 1) // 2
        n_types = count_orbit_types_no_H(p, 128)
        print(f"    p={p}, m={m}: {n_types} types (of 128 sampled)")


def nonlinear_moment_test(p):
    """Test if H is a NONLINEAR function of few moments."""
    m = (p - 1) // 2
    N = 1 << m
    limit = min(N, 128)

    print(f"\n{'='*70}")
    print(f"NONLINEAR MOMENT TESTS at p={p}")
    print(f"{'='*70}")

    data = []
    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        A = build_adj(p, S)
        H = compute_H_heldkarp(A, p)
        moments = compute_moments(p, S)
        data.append({'H': H, 'moments': moments})

    H_arr = np.array([d['H'] for d in data], dtype=float)

    # Test polynomial degrees
    moment_names = [f'S{2*j}' for j in range(2, m+1)]

    # Linear fit with k moments
    for k in range(1, min(len(moment_names), 6) + 1):
        cols = [np.array([d['moments'][2*(j+2)] for d in data], dtype=float)
                for j in range(k)]
        X = np.column_stack(cols + [np.ones(len(data))])
        _, _, _, _ = np.linalg.lstsq(X, H_arr, rcond=None)
        pred = X @ np.linalg.lstsq(X, H_arr, rcond=None)[0]
        err = np.max(np.abs(H_arr - pred))
        print(f"  Linear({','.join(moment_names[:k])}): max_err = {err:.2f}")

    # Quadratic fit with fewer moments
    if len(data) > 10:
        S4 = np.array([d['moments'][4] for d in data], dtype=float)
        S6 = np.array([d['moments'][6] for d in data], dtype=float)

        # Degree 2 in S4 alone
        X_q1 = np.column_stack([S4, S4**2, np.ones(len(data))])
        pred_q1 = X_q1 @ np.linalg.lstsq(X_q1, H_arr, rcond=None)[0]
        err_q1 = np.max(np.abs(H_arr - pred_q1))
        print(f"  Quadratic(S4): max_err = {err_q1:.2f}")

        # Degree 2 in (S4, S6)
        X_q2 = np.column_stack([S4, S6, S4**2, S4*S6, S6**2, np.ones(len(data))])
        pred_q2 = X_q2 @ np.linalg.lstsq(X_q2, H_arr, rcond=None)[0]
        err_q2 = np.max(np.abs(H_arr - pred_q2))
        print(f"  Quadratic(S4,S6): max_err = {err_q2:.2f}")

        if 8 in data[0]['moments']:
            S8 = np.array([d['moments'][8] for d in data], dtype=float)
            X_q3 = np.column_stack([S4, S6, S8, S4**2, S4*S6, S6**2, np.ones(len(data))])
            pred_q3 = X_q3 @ np.linalg.lstsq(X_q3, H_arr, rcond=None)[0]
            err_q3 = np.max(np.abs(H_arr - pred_q3))
            print(f"  Quadratic(S4,S6,S8): max_err = {err_q3:.2f}")

        # Degree 3 in S4 alone
        X_c1 = np.column_stack([S4, S4**2, S4**3, np.ones(len(data))])
        pred_c1 = X_c1 @ np.linalg.lstsq(X_c1, H_arr, rcond=None)[0]
        err_c1 = np.max(np.abs(H_arr - pred_c1))
        print(f"  Cubic(S4): max_err = {err_c1:.2f}")

        # Degree 3 in (S4, S6)
        X_c2 = np.column_stack([S4, S6, S4**2, S4*S6, S6**2,
                                 S4**3, S4**2*S6, S4*S6**2, S6**3,
                                 np.ones(len(data))])
        pred_c2 = X_c2 @ np.linalg.lstsq(X_c2, H_arr, rcond=None)[0]
        err_c2 = np.max(np.abs(H_arr - pred_c2))
        print(f"  Cubic(S4,S6): max_err = {err_c2:.2f}")


# ================================================================
# MAIN
# ================================================================
print("=" * 70)
print("D^2 MULTISET -> H DETERMINATION TEST")
print("=" * 70)

for p in [7, 11, 13]:
    test_D2_determines_H(p)

orbit_count_pattern()

for p in [11, 13]:
    nonlinear_moment_test(p)

print("\nDONE.")
