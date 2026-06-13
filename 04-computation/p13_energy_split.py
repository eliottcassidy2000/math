#!/usr/bin/env python3
"""
p13_energy_split.py -- kind-pasteur-2026-03-13-S60

At p=13, E=118 has TWO distinct H values: {3683797, 3704857}.
This script investigates what distinguishes them:
1. Their eigenvalue moment profiles (S4,...,S12)
2. Their cycle count profiles (c3,...,c13)
3. Their Fourier spectra
4. Whether there's a simple invariant separating them

Key question: is the split related to p=1 mod 4 (no Paley tournament)?
"""

import cmath
import numpy as np
from collections import defaultdict


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


def compute_eigenvalue_moments(p, S, max_moment=14):
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    lam_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)
        lam_vals.append(lam)
    moments = {}
    for k in range(2, max_moment + 1, 2):
        moments[k] = sum(d**k for d in D_vals)
    return moments, D_vals, lam_vals


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


def compute_cycle_counts(A, p, max_k=None):
    if max_k is None:
        max_k = p
    A_np = np.array(A, dtype=np.float64)
    counts = {}
    for k in range(3, max_k + 1, 2):
        Ak = np.linalg.matrix_power(A_np, k)
        counts[k] = int(round(np.trace(Ak))) // k
    return counts


def fourier_spectrum(S, p):
    omega = cmath.exp(2j * cmath.pi / p)
    spec = []
    for t in range(p):
        val = sum(omega ** (t * s) for s in S)
        spec.append(abs(val)**2)
    return spec


def analyze_p13_split():
    p = 13
    m = 6
    N = 64

    print("=" * 70)
    print(f"DETAILED ANALYSIS OF E=118 SPLIT AT p={p}")
    print("=" * 70)

    # Collect ALL data
    data = []
    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        H = compute_H_heldkarp(A, p)
        E = additive_energy(S, p)
        moments, D_vals, lam_vals = compute_eigenvalue_moments(p, S, max_moment=12)
        ck = compute_cycle_counts(A, p, max_k=p)
        fourier = fourier_spectrum(S, p)

        data.append({
            'bits': bits, 'S': S, 'H': H, 'E': E,
            'moments': moments, 'D_vals': D_vals, 'ck': ck,
            'fourier': fourier, 'lam_vals': lam_vals
        })

    # Find all E=118 entries
    e118 = [d for d in data if d['E'] == 118]
    print(f"\nE=118 entries: {len(e118)}")

    # Group by H
    by_H = defaultdict(list)
    for d in e118:
        by_H[d['H']].append(d)

    for H_val in sorted(by_H.keys()):
        group = by_H[H_val]
        print(f"\n  H = {H_val} ({len(group)} orientations):")
        for d in group:
            print(f"    bits={d['bits']:>3d}, S={d['S']}")

        d0 = group[0]
        print(f"    Moments: S4={d0['moments'][4]:.4f}, S6={d0['moments'][6]:.4f}, "
              f"S8={d0['moments'][8]:.4f}")
        print(f"    S10={d0['moments'][10]:.4f}, S12={d0['moments'][12]:.4f}")
        print(f"    D^2 vals: {[round(d**2, 4) for d in d0['D_vals']]}")
        print(f"    Cycle counts: ", end="")
        for k in sorted(d0['ck']):
            print(f"c{k}={d0['ck'][k]}", end=" ")
        print()

    # Compare the two groups
    if len(by_H) == 2:
        H_vals = sorted(by_H.keys())
        d_lo = by_H[H_vals[0]][0]
        d_hi = by_H[H_vals[1]][0]

        print(f"\n  --- COMPARISON (H={H_vals[0]} vs H={H_vals[1]}) ---")

        # Moment differences
        for k in sorted(d_lo['moments']):
            diff = d_hi['moments'][k] - d_lo['moments'][k]
            print(f"    S{k}: {d_lo['moments'][k]:.6f} vs {d_hi['moments'][k]:.6f}, "
                  f"diff = {diff:.6f}")

        # Cycle count differences
        print(f"\n    Cycle count differences:")
        for k in sorted(d_lo['ck']):
            diff = d_hi['ck'][k] - d_lo['ck'][k]
            if diff != 0:
                print(f"      c{k}: {d_lo['ck'][k]} vs {d_hi['ck'][k]}, diff = {diff}")
            else:
                print(f"      c{k}: {d_lo['ck'][k]} (same)")

        # D^2 spectrum
        d2_lo = sorted([d**2 for d in d_lo['D_vals']])
        d2_hi = sorted([d**2 for d in d_hi['D_vals']])
        print(f"\n    Sorted D^2 spectrum:")
        print(f"      Lo: {[round(x, 6) for x in d2_lo]}")
        print(f"      Hi: {[round(x, 6) for x in d2_hi]}")

        # Are D^2 multisets the same?
        same_d2 = (d2_lo == d2_hi)
        print(f"    Same D^2 multiset: {same_d2}")

        # Fourier spectrum comparison
        f_lo = sorted(d_lo['fourier'])
        f_hi = sorted(d_hi['fourier'])
        print(f"\n    Sorted Fourier |S_hat|^2 spectrum:")
        print(f"      Lo: {[round(x, 4) for x in f_lo]}")
        print(f"      Hi: {[round(x, 4) for x in f_hi]}")
        same_fourier = all(abs(a-b) < 1e-6 for a, b in zip(f_lo, f_hi))
        print(f"    Same Fourier spectrum: {same_fourier}")

        # alpha_1 (sum of cycle counts)
        a1_lo = sum(d_lo['ck'].values())
        a1_hi = sum(d_hi['ck'].values())
        print(f"\n    alpha_1 (total cycles): {a1_lo} vs {a1_hi}, diff = {a1_hi - a1_lo}")

        # H difference analysis
        delta_H = H_vals[1] - H_vals[0]
        print(f"\n    Delta H = {delta_H}")
        print(f"    Delta H / 2 = {delta_H / 2}")
        print(f"    Delta H / 4 = {delta_H / 4}")
        print(f"    Delta H / 8 = {delta_H / 8}")
        print(f"    H difference in terms of OCF: 2*da1 + 4*da2 + 8*da3 + ...")

    # Full table: all energy levels with moment profiles
    print(f"\n{'='*70}")
    print(f"ALL ORBIT TYPES AT p={p}")
    print(f"{'='*70}")

    by_E = defaultdict(list)
    for d in data:
        by_E[d['E']].append(d)

    print(f"\n  {'E':>4} {'H':>10} {'c3':>5} {'c5':>5} {'c7':>5} {'c9':>6} "
          f"{'c11':>6} {'c13':>8} {'alpha1':>8} {'count':>5}")
    for E in sorted(by_E.keys()):
        # Group by H within this energy
        h_groups = defaultdict(list)
        for d in by_E[E]:
            h_groups[d['H']].append(d)

        for H_val in sorted(h_groups.keys()):
            d = h_groups[H_val][0]
            ck = d['ck']
            a1 = sum(ck.values())
            count = len(h_groups[H_val])
            print(f"  {E:4d} {H_val:10d} {ck.get(3,0):5d} {ck.get(5,0):5d} {ck.get(7,0):5d} "
                  f"{ck.get(9,0):6d} {ck.get(11,0):6d} {ck.get(13,0):8d} {a1:8d} {count:5d}")

    # Test: does cycle count profile determine H?
    print(f"\n  --- Does cycle profile (c3,...,c13) determine H? ---")
    by_ck = defaultdict(set)
    for d in data:
        ck_key = tuple(d['ck'].get(k, 0) for k in range(3, p+1, 2))
        by_ck[ck_key].add(d['H'])

    all_unique = all(len(v) == 1 for v in by_ck.values())
    print(f"  Cycle profile determines H: {all_unique}")
    if not all_unique:
        for ck_key, h_set in by_ck.items():
            if len(h_set) > 1:
                print(f"    AMBIGUOUS: ck={ck_key}, H values={h_set}")

    # Test: does (c3, c5, c7) determine H?
    by_ck3 = defaultdict(set)
    for d in data:
        ck_key = tuple(d['ck'].get(k, 0) for k in [3, 5, 7])
        by_ck3[ck_key].add(d['H'])
    all_unique_3 = all(len(v) == 1 for v in by_ck3.values())
    print(f"  (c3,c5,c7) determines H: {all_unique_3}")

    # Quadratic test: H = f(S4, S6, S8, S4^2, S4*S6, S6^2)?
    print(f"\n  --- Quadratic moment fits ---")
    H_arr = np.array([d['H'] for d in data], dtype=float)
    S4_arr = np.array([d['moments'][4] for d in data], dtype=float)
    S6_arr = np.array([d['moments'][6] for d in data], dtype=float)
    S8_arr = np.array([d['moments'][8] for d in data], dtype=float)

    # Linear: S4, S6, S8
    X_lin = np.column_stack([S4_arr, S6_arr, S8_arr, np.ones(len(data))])
    c_lin, _, _, _ = np.linalg.lstsq(X_lin, H_arr, rcond=None)
    err_lin = np.max(np.abs(H_arr - X_lin @ c_lin))
    print(f"  Linear H = f(S4,S6,S8): max_err = {err_lin:.2f}")

    # Quadratic: S4, S6, S8, S4^2, S4*S6, S6^2
    X_quad = np.column_stack([S4_arr, S6_arr, S8_arr,
                              S4_arr**2, S4_arr*S6_arr, S6_arr**2,
                              np.ones(len(data))])
    c_quad, _, _, _ = np.linalg.lstsq(X_quad, H_arr, rcond=None)
    err_quad = np.max(np.abs(H_arr - X_quad @ c_quad))
    print(f"  Quadratic H = f(S4,S6,S8,S4^2,S4*S6,S6^2): max_err = {err_quad:.2f}")

    # Quadratic with S10
    S10_arr = np.array([d['moments'][10] for d in data], dtype=float)
    X_quad2 = np.column_stack([S4_arr, S6_arr, S8_arr, S10_arr,
                               S4_arr**2, S4_arr*S6_arr,
                               np.ones(len(data))])
    c_quad2, _, _, _ = np.linalg.lstsq(X_quad2, H_arr, rcond=None)
    err_quad2 = np.max(np.abs(H_arr - X_quad2 @ c_quad2))
    print(f"  H = f(S4,S6,S8,S10,S4^2,S4*S6): max_err = {err_quad2:.2f}")

    # p=1 mod 4 analysis
    print(f"\n{'='*70}")
    print(f"p mod 4 ANALYSIS")
    print(f"{'='*70}")
    print(f"  p={p} = 1 mod 4")
    print(f"  -1 is a QR mod {p}: {pow(p-1, (p-1)//2, p) == 1}")
    print(f"  QR set: {sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)}")
    print(f"  Since -1 is QR, QR set is symmetric under negation")
    print(f"  => QR is NOT a valid connection set (S cap (-S) != empty)")
    print(f"  => No Paley tournament at p=13")

    # Check which orientations are at E=118
    print(f"\n  E=118 orientations and their S sets:")
    for d in e118:
        S_sorted = sorted(d['S'])
        neg_S = sorted([(p - s) % p for s in d['S']])
        overlap = set(d['S']) & set(neg_S)
        print(f"    bits={d['bits']:>3d}, S={S_sorted}, H={d['H']}")


# Also check p=17 and p=19 to see if cancellation pattern holds
def check_moment_counts():
    """How many moments does H need at each prime?"""
    print(f"\n{'='*70}")
    print(f"MOMENT COUNT BY PRIME (how many S_k moments to determine H)")
    print(f"{'='*70}")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        N = 1 << m

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
            H = compute_H_heldkarp(A, p)
            moments, _, _ = compute_eigenvalue_moments(p, S, max_moment=2*m)

            data.append({'H': H, 'moments': moments})

        H_arr = np.array([d['H'] for d in data], dtype=float)

        # Count unique H values
        n_unique_H = len(set(d['H'] for d in data))
        print(f"\n  p={p}, m={m}, {limit} orientations, {n_unique_H} unique H values")

        # Test increasing moment counts
        moment_names = [f'S{2*j}' for j in range(2, m+1)]
        for n_mom in range(1, len(moment_names)+1):
            cols = [np.array([d['moments'][2*j] for d in data], dtype=float)
                    for j in range(2, 2+n_mom)]
            X = np.column_stack(cols + [np.ones(len(data))])
            coeffs, _, _, _ = np.linalg.lstsq(X, H_arr, rcond=None)
            pred = X @ coeffs
            max_err = np.max(np.abs(H_arr - pred))
            mom_str = ','.join(moment_names[:n_mom])
            exact = max_err < 0.01
            print(f"    H = f({mom_str}): max_err = {max_err:.4f} {'EXACT' if exact else ''}")
            if exact:
                print(f"    => H needs {n_mom} moments at p={p} "
                      f"(m-1={m-1}, p mod 4 = {p % 4})")
                break


if __name__ == '__main__':
    analyze_p13_split()
    check_moment_counts()
    print("\nDONE.")
