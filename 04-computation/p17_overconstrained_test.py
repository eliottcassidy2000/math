#!/usr/bin/env python3
"""
p17_overconstrained_test.py -- kind-pasteur-2026-03-13-S60

At p=17 (m=8), there are 16 unique D^2 multisets (orbit types).
With m-1 = 7 independent moments (S4,...,S16), a linear fit has 8 params.
16 data points with 8 params => GENUINELY OVERCONSTRAINED.

If H = linear(S4,...,S16) is still exact, it would be a real algebraic
identity, not a trivial consequence of dimension counting.

BUT: H computation at p=17 is O(2^17 * 17^2) per orientation, too slow
in Python (~10 min each, 128 orientations = ~20 hours).

ALTERNATIVE APPROACH: Test cycle counts c_k as linear functions of moments.
Since c_k = (1/k)[m^k + 2*sum Re(z^k)] and Re(z^k) depends on S_2,...,S_k,
this is algebraically exact. But does the OVERCONSTRAINED linear system
for H (which combines alpha_j contributions) remain consistent?

Strategy:
1. Compute cycle counts c_k for all orientations (fast via matrix power)
2. Test c_k = f(S4,...,S_{k-1}) -- should be exact by theory
3. Test alpha_1 = sum c_k as function of moments -- should need (p-3)/2 moments
4. Compute alpha_2 = disj3 as function of S4 (from THM-156)
5. The question becomes: can we predict H from alpha_1, alpha_2 alone?
"""

import cmath
import numpy as np
from collections import defaultdict
import time


def build_adj(p, S):
    S_set = set(S)
    A = np.zeros((p, p), dtype=np.float64)
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


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
    return moments, D_vals


def compute_D2_multiset(p, S):
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D2_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D2_vals.append(lam.imag ** 2)
    return tuple(sorted(round(d, 8) for d in D2_vals))


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


def analyze_cycle_moments(p, max_k=None):
    """Analyze cycle count - moment dependency at prime p."""
    m = (p - 1) // 2
    N = 1 << m
    limit = min(N, 256)
    if max_k is None:
        max_k = min(p, 17)

    print(f"\n{'='*70}")
    print(f"CYCLE-MOMENT ANALYSIS at p={p}, m={m}, {limit} orientations")
    print(f"{'='*70}")

    t0 = time.time()
    data = []
    by_D2 = defaultdict(list)

    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        moments, D_vals = compute_moments(p, S, max_moment=2*m)
        D2 = compute_D2_multiset(p, S)
        E = additive_energy(S, p)

        # Cycle counts via matrix power
        ck = {}
        for k in range(3, max_k + 1, 2):
            Ak = np.linalg.matrix_power(A, k)
            ck[k] = int(round(np.trace(Ak))) // k

        alpha_1 = sum(ck.values())

        row = {
            'bits': bits, 'S': S, 'E': E, 'D2': D2,
            'moments': moments, 'ck': ck, 'alpha_1': alpha_1
        }
        data.append(row)
        by_D2[D2].append(row)

    t1 = time.time()
    n_types = len(by_D2)
    print(f"  Computed in {t1-t0:.1f}s")
    print(f"  {n_types} unique D^2 multisets")

    # Verify: c_k is exact linear function of moments
    print(f"\n  --- CYCLE COUNT VERIFICATION ---")
    for k in range(3, max_k + 1, 2):
        ck_arr = np.array([d['ck'][k] for d in data], dtype=float)
        if np.std(ck_arr) < 1e-6:
            print(f"  c_{k}: CONSTANT = {int(ck_arr[0])}")
            continue

        # Theory: c_k depends on S4,...,S_{k-1} = floor((k-3)/2) moments
        n_params = (k - 3) // 2
        if n_params == 0:
            print(f"  c_{k}: should be constant (0 params)")
            continue

        moment_indices = [2*j for j in range(2, 2+n_params)]
        cols = [np.array([d['moments'][mi] for d in data], dtype=float)
                for mi in moment_indices]
        X = np.column_stack(cols + [np.ones(limit)])
        coeffs, _, _, _ = np.linalg.lstsq(X, ck_arr, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(ck_arr - pred))
        mom_names = ','.join(f'S{mi}' for mi in moment_indices)
        print(f"  c_{k} = f({mom_names}): max_err = {max_err:.4f} "
              f"{'EXACT' if max_err < 0.01 else ''}")

    # Alpha_1 = sum c_k as function of moments
    print(f"\n  --- ALPHA_1 MOMENT DEPENDENCY ---")
    a1_arr = np.array([d['alpha_1'] for d in data], dtype=float)

    moment_names = [f'S{2*j}' for j in range(2, m+1)]
    for n_mom in range(1, min(len(moment_names), 10) + 1):
        cols = [np.array([d['moments'][2*(j+2)] for d in data], dtype=float)
                for j in range(n_mom)]
        X = np.column_stack(cols + [np.ones(limit)])
        coeffs, _, _, _ = np.linalg.lstsq(X, a1_arr, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(a1_arr - pred))
        mom_str = ','.join(moment_names[:n_mom])
        exact = max_err < 0.01
        print(f"  alpha_1 = f({mom_str}): max_err = {max_err:.4f} "
              f"{'EXACT' if exact else ''}")
        if exact:
            break

    # Check: does D^2 multiset determine alpha_1?
    print(f"\n  --- D^2 MULTISET -> alpha_1 ---")
    a1_unique = True
    for D2, entries in by_D2.items():
        a1_vals = set(e['alpha_1'] for e in entries)
        if len(a1_vals) > 1:
            a1_unique = False
            print(f"  FAILURE: D^2 multiset has alpha_1 values {a1_vals}")
    if a1_unique:
        print(f"  D^2 multiset determines alpha_1: YES")

    # Check: does D^2 determine each c_k?
    for k in range(3, max_k + 1, 2):
        ck_unique = True
        for D2, entries in by_D2.items():
            ck_vals = set(e['ck'][k] for e in entries)
            if len(ck_vals) > 1:
                ck_unique = False
                break
        if not ck_unique:
            print(f"  D^2 -> c_{k}: FAILS")
        else:
            if np.std([data[0]['ck'][k] for data in by_D2.values()]) > 0:
                print(f"  D^2 -> c_{k}: YES (non-trivially)")

    # Show orbit types with cycle counts
    print(f"\n  --- ORBIT TYPES ---")
    header = f"  {'Type':>4} {'count':>5} {'E':>5}"
    for k in range(3, min(max_k+1, 18), 2):
        header += f" {'c'+str(k):>8}"
    header += f" {'alpha_1':>12}"
    print(header)

    for i, D2 in enumerate(sorted(by_D2.keys())):
        entries = by_D2[D2]
        d = entries[0]
        line = f"  {i:4d} {len(entries):5d} {d['E']:5d}"
        for k in range(3, min(max_k+1, 18), 2):
            line += f" {d['ck'][k]:8d}"
        line += f" {d['alpha_1']:12d}"
        print(line)

    # Key test: alpha_1 as function of (S4,...,S_{p-3}) -- the OVERCONSTRAINED test
    if n_types > m:
        print(f"\n  --- OVERCONSTRAINED TEST ---")
        print(f"  n_types = {n_types}, m-1 = {m-1}")
        print(f"  {n_types} data points, {m-1+1} = {m} linear params")
        print(f"  OVERCONSTRAINED by {n_types - m} equations")

        # Can alpha_1 = linear(S4,...,S_{p-3})?
        # Build the reduced system using one representative per D^2 type
        a1_rep = []
        moment_rep = []
        for D2 in sorted(by_D2.keys()):
            d = by_D2[D2][0]
            a1_rep.append(d['alpha_1'])
            moment_rep.append([d['moments'][2*j] for j in range(2, m+1)])

        a1_rep = np.array(a1_rep, dtype=float)
        X_rep = np.column_stack([np.array(moment_rep), np.ones(n_types)])
        coeffs_rep, _, _, _ = np.linalg.lstsq(X_rep, a1_rep, rcond=None)
        pred_rep = X_rep @ coeffs_rep
        max_err_rep = np.max(np.abs(a1_rep - pred_rep))
        print(f"  alpha_1 = linear(S4,...,S{2*m}): max_err = {max_err_rep:.4f}")

        if max_err_rep < 0.01:
            print(f"  *** GENUINE LINEAR IDENTITY! ***")
            print(f"  Coefficients: {coeffs_rep}")
        else:
            print(f"  Linear identity FAILS -- need nonlinear terms")

            # Try quadratic
            mom_arr = np.array(moment_rep)
            quad_cols = [mom_arr]
            # Add quadratic terms
            for i in range(mom_arr.shape[1]):
                for j in range(i, mom_arr.shape[1]):
                    quad_cols.append((mom_arr[:, i] * mom_arr[:, j]).reshape(-1, 1))
            X_quad = np.column_stack(quad_cols + [np.ones(n_types)])
            n_quad_params = X_quad.shape[1]
            print(f"  Quadratic params: {n_quad_params}")

            if n_quad_params <= n_types:
                coeffs_q, _, _, _ = np.linalg.lstsq(X_quad, a1_rep, rcond=None)
                pred_q = X_quad @ coeffs_q
                max_err_q = np.max(np.abs(a1_rep - pred_q))
                print(f"  alpha_1 = quadratic(S4,...,S{2*m}): max_err = {max_err_q:.4f}")
                if max_err_q < 0.01:
                    print(f"  *** QUADRATIC IDENTITY! ***")

    return data, by_D2


# ================================================================
# MAIN
# ================================================================
print("=" * 70)
print("OVERCONSTRAINED MOMENT TESTS (cycle counts)")
print("=" * 70)

for p in [13, 17, 19]:
    analyze_cycle_moments(p)

print("\nDONE.")
