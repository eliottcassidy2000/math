#!/usr/bin/env python3
"""
H_moment_p13_p17.py -- kind-pasteur-2026-03-13-S60

At p=11: 4 unique moment orbits, H = f(S4,S6,S8) trivially (4 equations, 4 unknowns).
At p=13: 5 unique energy orbits => 5 > 4 unknowns, so H = f(S4,S6,S8) is a genuine test!
At p=17: 9 unique energy orbits => 9 >> 4 unknowns, even more stringent.

This script computes H via Held-Karp at p=13 (feasible: 2^13 * 169 ~ 1.4M operations)
and tests how many eigenvalue moments are needed to determine H.

Key prediction from theory:
  c_k needs floor((k-3)/2) moment parameters
  alpha_1 = sum c_k needs (p-3)/2 parameters
  alpha_2 needs 1 parameter (S4)
  If alpha_3 introduces no new parameters, then H needs at most (p-3)/2 parameters
  But we saw S10 gets cancelled at p=11, suggesting fewer suffice.
"""

import cmath
import numpy as np
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s)%p] = 1
    return A


def compute_H_heldkarp(A, n):
    """Compute H(T) = number of Hamiltonian paths using Held-Karp."""
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


def compute_eigenvalue_moments(p, S, max_moment=14):
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)
    moments = {}
    for k in range(2, max_moment + 1, 2):
        moments[k] = sum(d**k for d in D_vals)
    return moments


def additive_energy(S, p):
    S_set = set(S)
    energy = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    energy += 1
    return energy


def analyze_prime(p):
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"H vs MOMENTS at p={p}, m={m}, N={N}")
    print(f"{'='*70}")

    data = []
    limit = min(N, 128)  # cap for very large p

    import time
    t0 = time.time()

    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        H = compute_H_heldkarp(A, p)
        moments = compute_eigenvalue_moments(p, S, max_moment=min(p + 1, 14))
        E = additive_energy(S, p)

        data.append({'bits': bits, 'H': H, 'E': E, **{f'S{k}': v for k, v in moments.items()}})

    t1 = time.time()
    print(f"  Computed {len(data)} orientations in {t1-t0:.1f}s")

    # Group by energy to find unique types
    by_E = defaultdict(list)
    for d in data:
        by_E[d['E']].append(d)

    n_types = len(by_E)
    print(f"  {n_types} unique energy levels")

    # Show all unique types
    print(f"\n  {'E':>5} {'H':>10} {'S4':>10} {'S6':>12} {'S8':>14} {'S10':>14} {'count':>6}")
    for E in sorted(by_E.keys()):
        d = by_E[E][0]
        H_set = set(dd['H'] for dd in by_E[E])
        print(f"  {E:5d} {d['H']:10d} {d['S4']:10.2f} {d.get('S6',0):12.2f} "
              f"{d.get('S8',0):14.2f} {d.get('S10',0):14.2f} {len(by_E[E]):6d}", end='')
        if len(H_set) > 1:
            print(f"  *** MULTIPLE H VALUES: {H_set} ***", end='')
        print()

    # Check if E determines H
    E_determines_H = all(len(set(dd['H'] for dd in by_E[E])) == 1 for E in by_E)
    print(f"\n  E determines H: {E_determines_H}")

    # Test H = f(moments) with increasing number of moments
    H_arr = np.array([d['H'] for d in data], dtype=float)

    moment_names = ['S4', 'S6', 'S8', 'S10', 'S12']
    for n_mom in range(1, min(len(moment_names), n_types) + 1):
        cols = [np.array([d[moment_names[i]] for d in data], dtype=float)
                for i in range(n_mom)]
        X = np.column_stack(cols + [np.ones(len(data))])
        coeffs, _, _, _ = np.linalg.lstsq(X, H_arr, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(H_arr - pred))
        mom_str = ','.join(moment_names[:n_mom])
        print(f"  H = f({mom_str}): max_err = {max_err:.4f} "
              f"{'*** EXACT ***' if max_err < 0.01 else ''}")

    # Also test: H = f(E) only
    E_arr = np.array([d['E'] for d in data], dtype=float)
    X_E = np.column_stack([E_arr, np.ones(len(data))])
    coeffs_E, _, _, _ = np.linalg.lstsq(X_E, H_arr, rcond=None)
    max_err_E = np.max(np.abs(H_arr - X_E @ coeffs_E))
    print(f"\n  H = f(E): max_err = {max_err_E:.4f}")

    # H = f(E, S6)
    if 'S6' in data[0]:
        S6_arr = np.array([d['S6'] for d in data], dtype=float)
        X_ES = np.column_stack([E_arr, S6_arr, np.ones(len(data))])
        coeffs_ES, _, _, _ = np.linalg.lstsq(X_ES, H_arr, rcond=None)
        max_err_ES = np.max(np.abs(H_arr - X_ES @ coeffs_ES))
        print(f"  H = f(E,S6): max_err = {max_err_ES:.4f}")

    if 'S8' in data[0]:
        S8_arr = np.array([d['S8'] for d in data], dtype=float)
        X_ESS = np.column_stack([E_arr, S6_arr, S8_arr, np.ones(len(data))])
        coeffs_ESS, _, _, _ = np.linalg.lstsq(X_ESS, H_arr, rcond=None)
        max_err_ESS = np.max(np.abs(H_arr - X_ESS @ coeffs_ESS))
        print(f"  H = f(E,S6,S8): max_err = {max_err_ESS:.4f}")

    if 'S10' in data[0]:
        S10_arr = np.array([d['S10'] for d in data], dtype=float)
        X_ESSS = np.column_stack([E_arr, S6_arr, S8_arr, S10_arr, np.ones(len(data))])
        coeffs_ESSS, _, _, _ = np.linalg.lstsq(X_ESSS, H_arr, rcond=None)
        max_err_ESSS = np.max(np.abs(H_arr - X_ESSS @ coeffs_ESSS))
        print(f"  H = f(E,S6,S8,S10): max_err = {max_err_ESSS:.4f}")


# ================================================================
# MAIN
# ================================================================
print("=" * 70)
print("H vs EIGENVALUE MOMENTS: Testing at p=7,11,13")
print("=" * 70)

for p in [7, 11, 13]:
    analyze_prime(p)

print("\nDONE.")
