#!/usr/bin/env python3
"""
ck_energy_relationship.py -- kind-pasteur-2026-03-13-S60

Since disj3 = (p/4)*E + C (THM-156), and E is related to sum D_t^4,
the question is: how do HIGHER cycle counts c_k relate to E?

From the eigenvalue expansion:
  c_k = tr(A^k)/k = [m^k + 2*sum Re(z_t^k)] / k
  Re(z^k) = sum_{j even} C(k,j)(-1/2)^{k-j}(-1)^{j/2} D^j

The D^2 coefficient of Re(z^k) is C(k,2)*(-1/2)^{k-2}*(-1) = -k(k-1)/2^{k-1}
The D^4 coefficient is C(k,4)*(-1/2)^{k-4}*1 = k(k-1)(k-2)(k-3)/(24*2^{k-4})

So c_k depends on sum D^2 (constant), sum D^4 (= linear in E), sum D^6, etc.
For k=3: only D^2 term => c_3 constant (known)
For k=5: D^2 and D^4 => c_5 linear in E (known, THM-155/156)
For k=7: D^2, D^4, D^6 => c_7 depends on E AND sum D^6
For k=9: D^2, D^4, D^6, D^8 => c_9 depends on E, sum D^6, sum D^8

This script:
1. Verifies this structure explicitly
2. Finds the MINIMAL set of eigenvalue moments that determine each c_k
3. Tests if c_7 is a function of (E, one extra parameter)
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


def compute_eigenvalue_moments(p, S, max_moment=12):
    """Compute S_k = sum_{t=1}^{m} D_t^k for k=2,4,...,max_moment."""
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)

    moments = {}
    for k in range(2, max_moment + 1, 2):
        moments[k] = sum(d**k for d in D_vals)
    return moments, D_vals


def compute_cycle_counts(A, p, max_k=None):
    if max_k is None:
        max_k = p
    A_np = np.array(A, dtype=np.float64)
    counts = {}
    for k in range(3, max_k + 1, 2):
        Ak = np.linalg.matrix_power(A_np, k)
        counts[k] = int(round(np.trace(Ak))) // k
    return counts


def eigenvalue_dependency_analysis(p):
    """For each c_k, determine which moments S_{2j} = sum D^{2j} it depends on."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"EIGENVALUE DEPENDENCY at p={p}, m={m}")
    print(f"{'='*70}")

    # Collect data for all orientations
    data = []
    limit = min(N, 128)
    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        ck = compute_cycle_counts(A, p, max_k=min(p, 13))
        moments, _ = compute_eigenvalue_moments(p, S, max_moment=12)
        E = additive_energy(S, p)

        row = {'bits': bits, 'E': E}
        row.update({f'c{k}': v for k, v in ck.items()})
        row.update({f'S{k}': v for k, v in moments.items()})
        data.append(row)

    # For each c_k, try to express as function of moments
    max_k = min(p, 13)
    for k in range(3, max_k + 1, 2):
        ck_arr = np.array([d[f'c{k}'] for d in data], dtype=float)
        if np.std(ck_arr) < 1e-6:
            print(f"\n  c_{k}: CONSTANT = {int(ck_arr[0])}")
            continue

        # Try with increasing numbers of moments
        print(f"\n  c_{k} dependency analysis:")

        # 1-variable: S4 only (= linear in E)
        S4_arr = np.array([d['S4'] for d in data], dtype=float)
        if np.std(S4_arr) > 1e-6:
            X1 = np.column_stack([S4_arr, np.ones(len(data))])
            coeffs, _, _, _ = np.linalg.lstsq(X1, ck_arr, rcond=None)
            pred = X1 @ coeffs
            err1 = np.max(np.abs(ck_arr - pred))
            print(f"    c_{k} = f(S4):         max_err = {err1:.4f} {'EXACT' if err1 < 0.01 else ''}")

        # 2-variable: S4, S6
        if f'S6' in data[0]:
            S6_arr = np.array([d['S6'] for d in data], dtype=float)
            X2 = np.column_stack([S4_arr, S6_arr, np.ones(len(data))])
            coeffs2, _, _, _ = np.linalg.lstsq(X2, ck_arr, rcond=None)
            pred2 = X2 @ coeffs2
            err2 = np.max(np.abs(ck_arr - pred2))
            print(f"    c_{k} = f(S4,S6):      max_err = {err2:.4f} {'EXACT' if err2 < 0.01 else ''}")

        # 3-variable: S4, S6, S8
        if f'S8' in data[0]:
            S8_arr = np.array([d['S8'] for d in data], dtype=float)
            X3 = np.column_stack([S4_arr, S6_arr, S8_arr, np.ones(len(data))])
            coeffs3, _, _, _ = np.linalg.lstsq(X3, ck_arr, rcond=None)
            pred3 = X3 @ coeffs3
            err3 = np.max(np.abs(ck_arr - pred3))
            print(f"    c_{k} = f(S4,S6,S8):   max_err = {err3:.4f} {'EXACT' if err3 < 0.01 else ''}")

        # 4-variable: S4, S6, S8, S10
        if f'S10' in data[0]:
            S10_arr = np.array([d['S10'] for d in data], dtype=float)
            X4 = np.column_stack([S4_arr, S6_arr, S8_arr, S10_arr, np.ones(len(data))])
            coeffs4, _, _, _ = np.linalg.lstsq(X4, ck_arr, rcond=None)
            pred4 = X4 @ coeffs4
            err4 = np.max(np.abs(ck_arr - pred4))
            print(f"    c_{k} = f(S4..S10):    max_err = {err4:.4f} {'EXACT' if err4 < 0.01 else ''}")

        # 5-variable: S4, S6, S8, S10, S12
        if f'S12' in data[0]:
            S12_arr = np.array([d['S12'] for d in data], dtype=float)
            X5 = np.column_stack([S4_arr, S6_arr, S8_arr, S10_arr, S12_arr, np.ones(len(data))])
            coeffs5, _, _, _ = np.linalg.lstsq(X5, ck_arr, rcond=None)
            pred5 = X5 @ coeffs5
            err5 = np.max(np.abs(ck_arr - pred5))
            print(f"    c_{k} = f(S4..S12):    max_err = {err5:.4f} {'EXACT' if err5 < 0.01 else ''}")

    # The number of independent parameters needed is floor((k-3)/2) + 1
    print(f"\n  THEORY: c_k depends on S4, S6, ..., S_{{k-1}}")
    print(f"  Number of parameters: floor((k-3)/2) for k odd")
    print(f"  k=3: 0 parameters (constant)")
    print(f"  k=5: 1 parameter (S4 = linear in E)")
    print(f"  k=7: 2 parameters (S4, S6)")
    print(f"  k=9: 3 parameters (S4, S6, S8)")
    print(f"  k=11: 4 parameters (S4, S6, S8, S10)")
    print(f"  k=13: 5 parameters (S4, S6, S8, S10, S12)")


def H_moment_dependency(p):
    """Since H = sum 2^j * alpha_j and alpha_j depends on cycle interaction structure,
    test how many moments determine H."""
    m = (p - 1) // 2
    N = 1 << m

    if p > 11:
        print(f"\n  [p={p} too large for H computation]")
        return

    print(f"\n{'='*70}")
    print(f"H vs EIGENVALUE MOMENTS at p={p}")
    print(f"{'='*70}")

    from itertools import combinations

    data = []
    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        # Compute H via Held-Karp
        n = p
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
        H = sum(dp.get((full, v), 0) for v in range(n))

        moments, _ = compute_eigenvalue_moments(p, S, max_moment=10)
        E = additive_energy(S, p)

        data.append({'bits': bits, 'H': H, 'E': E, **{f'S{k}': v for k, v in moments.items()}})

    H_arr = np.array([d['H'] for d in data], dtype=float)

    # Test H = f(S4)
    S4_arr = np.array([d['S4'] for d in data], dtype=float)
    X1 = np.column_stack([S4_arr, np.ones(N)])
    coeffs, _, _, _ = np.linalg.lstsq(X1, H_arr, rcond=None)
    err1 = np.max(np.abs(H_arr - X1 @ coeffs))
    print(f"  H = f(S4): max_err = {err1:.2f}")

    # H = f(S4, S6)
    S6_arr = np.array([d['S6'] for d in data], dtype=float)
    X2 = np.column_stack([S4_arr, S6_arr, np.ones(N)])
    coeffs2, _, _, _ = np.linalg.lstsq(X2, H_arr, rcond=None)
    err2 = np.max(np.abs(H_arr - X2 @ coeffs2))
    print(f"  H = f(S4,S6): max_err = {err2:.2f}")

    # H = f(S4, S6, S8)
    S8_arr = np.array([d['S8'] for d in data], dtype=float)
    X3 = np.column_stack([S4_arr, S6_arr, S8_arr, np.ones(N)])
    coeffs3, _, _, _ = np.linalg.lstsq(X3, H_arr, rcond=None)
    err3 = np.max(np.abs(H_arr - X3 @ coeffs3))
    print(f"  H = f(S4,S6,S8): max_err = {err3:.2f}")

    # H = f(S4, S6, S8, S10)
    S10_arr = np.array([d['S10'] for d in data], dtype=float)
    X4 = np.column_stack([S4_arr, S6_arr, S8_arr, S10_arr, np.ones(N)])
    coeffs4, _, _, _ = np.linalg.lstsq(X4, H_arr, rcond=None)
    err4 = np.max(np.abs(H_arr - X4 @ coeffs4))
    print(f"  H = f(S4,S6,S8,S10): max_err = {err4:.2f}")

    if err4 < 0.01:
        print(f"  *** H is an EXACT linear function of (S4,S6,S8,S10) ***")
        print(f"  Coefficients: {coeffs4}")

    # Also try with energy E instead of S4 (since E = S4/p + const)
    E_arr = np.array([d['E'] for d in data], dtype=float)
    X_E = np.column_stack([E_arr, S6_arr, S8_arr, S10_arr, np.ones(N)])
    coeffs_E, _, _, _ = np.linalg.lstsq(X_E, H_arr, rcond=None)
    err_E = np.max(np.abs(H_arr - X_E @ coeffs_E))
    print(f"\n  H = f(E,S6,S8,S10): max_err = {err_E:.2f}")
    if err_E < 0.01:
        print(f"  *** H = {coeffs_E[0]:.4f}*E + {coeffs_E[1]:.4f}*S6 + {coeffs_E[2]:.4f}*S8 + {coeffs_E[3]:.4f}*S10 + {coeffs_E[4]:.2f} ***")

    # Show unique data
    print(f"\n  Unique (E, S6, H) triples:")
    seen = {}
    for d in data:
        key = (d['E'], round(d['S6'], 2))
        if key not in seen:
            seen[key] = d['H']
    for key in sorted(seen.keys()):
        print(f"    E={key[0]:>4}, S6={key[1]:>10.2f}, H={seen[key]:>8}")


# ================================================================
# MAIN
# ================================================================

print("=" * 70)
print("CYCLE COUNT - ENERGY MOMENT DEPENDENCY")
print("=" * 70)

for p in [7, 11, 13]:
    eigenvalue_dependency_analysis(p)

for p in [7, 11]:
    H_moment_dependency(p)

print("\nDONE.")
