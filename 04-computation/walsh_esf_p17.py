#!/usr/bin/env python3
"""
walsh_esf_p17.py -- Test H = linear(e_j) at p=17

At p=17, m=8, there are 256 circulant tournament orientations.
We need to:
1. Compute H for all 256 (expensive: ~6 min total)
2. Compute Q_k = |S_hat(k)|^2 for all orientations (fast)
3. Group by unique (e_2,...,e_8) profiles
4. Test if H is a LINEAR function of (e_2,...,e_8) [7 variables]
   Need > 8 unique profiles for this to be overdetermined.

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
import time
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def main():
    p = 17
    m = (p - 1) // 2  # = 8
    QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    omega = cmath.exp(2j * cmath.pi / p)

    print("=" * 70)
    print(f"ESF ANALYSIS AT p={p}, m={m}")
    print("=" * 70)

    pairs = [(s, p - s) for s in range(1, m + 1)]

    # Step 1: Compute H for all 256 orientations
    print(f"\nComputing H for all {1 << m} orientations...")
    t0 = time.time()

    H_dict = {}
    Q_dict = {}

    for bits in range(1 << m):
        sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
        S = []
        for i, (a, b) in enumerate(pairs):
            S.append(a if sigma[i] == 1 else b)
        S = sorted(S)
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_dict[sigma] = H

        # Compute Q_k
        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S)
            Q_vals.append(abs(val)**2)
        Q_dict[sigma] = Q_vals

        if bits % 32 == 0:
            elapsed = time.time() - t0
            print(f"  {bits}/{1<<m} done ({elapsed:.1f}s)")

    t1 = time.time()
    print(f"All H computed in {t1-t0:.1f}s")

    # Step 2: Elementary symmetric functions
    print(f"\nComputing elementary symmetric functions...")

    all_esf = {}
    for sigma in H_dict:
        Q = Q_dict[sigma]
        e_vals = [1.0]
        for j in range(1, m + 1):
            ej = sum(math.prod(Q[i] for i in combo)
                    for combo in combinations(range(m), j))
            e_vals.append(ej)
        all_esf[sigma] = e_vals

    # Check integrality
    max_frac = 0
    for sigma in H_dict:
        e_vals = all_esf[sigma]
        for j in range(1, m + 1):
            frac = abs(e_vals[j] - round(e_vals[j]))
            max_frac = max(max_frac, frac)
    print(f"Max fractional part: {max_frac:.2e}")
    print(f"All e_j integers: {'YES' if max_frac < 1e-4 else 'NO'}")

    # Step 3: Group by e_j profile
    esf_groups = defaultdict(list)
    for sigma in H_dict:
        H = H_dict[sigma]
        e_profile = tuple(round(all_esf[sigma][j]) for j in range(2, m + 1))
        esf_groups[e_profile].append((sigma, H))

    n_unique = len(esf_groups)
    n_vars = m - 1  # = 7
    print(f"\nNumber of unique (e_2,...,e_{m}) profiles: {n_unique}")
    print(f"Number of variables (e_2,...,e_{m}): {n_vars}")
    print(f"System is {'overdetermined' if n_unique > n_vars + 1 else 'underdetermined'}")

    # Print table
    print(f"\n{'e_2':>8} {'e_3':>8} {'e_4':>10} {'e_5':>10} {'e_6':>10} "
          f"{'e_7':>10} {'e_8':>10} {'H':>14} {'count':>6}")
    print("-" * 100)
    for e_profile, items in sorted(esf_groups.items(), key=lambda x: -x[1][0][1]):
        H = items[0][1]
        e_str = " ".join(f"{e_profile[j]:>10}" for j in range(min(len(e_profile), 7)))
        print(f"    {e_str} {H:>14} {len(items):>6}")

    # Step 4: Test linearity
    print(f"\n{'='*70}")
    print(f"LINEARITY TEST: H = c_0 + sum c_j * e_{{j+1}}")
    print(f"{'='*70}")

    profiles_list = list(esf_groups.items())
    A_mat = []
    b_vec = []
    for e_profile, items in profiles_list:
        A_mat.append([1.0] + list(e_profile))
        b_vec.append(float(items[0][1]))

    n_eq = len(A_mat)
    n_col = n_vars + 1  # = 8

    print(f"Equations: {n_eq}, Variables: {n_col}")

    # Solve via least squares
    XTX = [[0.0]*n_col for _ in range(n_col)]
    XTb = [0.0]*n_col
    for r in range(n_eq):
        for i in range(n_col):
            for j in range(n_col):
                XTX[i][j] += A_mat[r][i] * A_mat[r][j]
            XTb[i] += A_mat[r][i] * b_vec[r]

    # Gaussian elimination
    aug = [XTX[i][:] + [XTb[i]] for i in range(n_col)]
    for col in range(n_col):
        max_row = max(range(col, n_col), key=lambda r: abs(aug[r][col]))
        aug[col], aug[max_row] = aug[max_row], aug[col]
        pivot = aug[col][col]
        if abs(pivot) < 1e-20:
            print(f"Singular at col {col}")
            break
        for j in range(col, n_col + 1):
            aug[col][j] /= pivot
        for r in range(n_col):
            if r == col:
                continue
            factor = aug[r][col]
            for j in range(col, n_col + 1):
                aug[r][j] -= factor * aug[col][j]

    coeffs = [aug[i][n_col] for i in range(n_col)]

    # Check residuals
    residuals = []
    for r in range(n_eq):
        H_pred = sum(coeffs[i] * A_mat[r][i] for i in range(n_col))
        residuals.append(b_vec[r] - H_pred)

    max_resid = max(abs(r) for r in residuals)
    rms_resid = (sum(r**2 for r in residuals) / n_eq) ** 0.5
    print(f"\nMax residual = {max_resid:.6f}")
    print(f"RMS residual = {rms_resid:.6f}")
    print(f"Is H linear in e_j? {'YES' if max_resid < 0.01 else 'NO'}")

    if max_resid < 0.01:
        print(f"\nCoefficients:")
        for j in range(n_col):
            if j == 0:
                print(f"  c_0 = {coeffs[j]:.10f}")
            else:
                print(f"  c(e_{j+1}) = {coeffs[j]:.10f}")

        # Rational reconstruction
        print(f"\nRational reconstruction:")
        for j in range(n_col):
            found = False
            for denom in range(1, 1000):
                numer = round(coeffs[j] * denom)
                if abs(coeffs[j] - numer / denom) < 1e-4:
                    print(f"  c_{j} = {numer}/{denom}")
                    found = True
                    break
            if not found:
                print(f"  c_{j} = {coeffs[j]:.10f} (not recognized)")
    else:
        print(f"\nLinear model fails. Testing quadratic...")

        # Show largest residuals
        print(f"\nTop 5 residuals:")
        indexed = sorted(enumerate(residuals), key=lambda x: -abs(x[1]))[:5]
        for idx, resid in indexed:
            e_profile = profiles_list[idx][0]
            H = profiles_list[idx][1][0][1]
            print(f"  e_profile={e_profile[:3]}..., H={H}, resid={resid:>+.2f}")

        # Try quadratic
        quad_features = []
        for e_profile, items in profiles_list:
            row = [1.0] + list(e_profile)
            for i in range(len(e_profile)):
                for j in range(i, len(e_profile)):
                    row.append(e_profile[i] * e_profile[j])
            quad_features.append(row)

        n_qcol = len(quad_features[0])
        print(f"\nQuadratic: {n_qcol} features, {n_eq} data points")

        if n_qcol <= n_eq:
            XTX = [[0.0]*n_qcol for _ in range(n_qcol)]
            XTb = [0.0]*n_qcol
            for r in range(n_eq):
                for i in range(n_qcol):
                    for j in range(n_qcol):
                        XTX[i][j] += quad_features[r][i] * quad_features[r][j]
                    XTb[i] += quad_features[r][i] * b_vec[r]

            aug = [XTX[i][:] + [XTb[i]] for i in range(n_qcol)]
            for col in range(n_qcol):
                max_row = max(range(col, n_qcol), key=lambda r: abs(aug[r][col]))
                aug[col], aug[max_row] = aug[max_row], aug[col]
                pivot = aug[col][col]
                if abs(pivot) < 1e-20:
                    break
                for j in range(col, n_qcol + 1):
                    aug[col][j] /= pivot
                for r in range(n_qcol):
                    if r == col:
                        continue
                    factor = aug[r][col]
                    for j in range(col, n_qcol + 1):
                        aug[r][j] -= factor * aug[col][j]

            qcoeffs = [aug[i][n_qcol] for i in range(n_qcol)]
            qresid = []
            for r in range(n_eq):
                H_pred = sum(qcoeffs[i] * quad_features[r][i] for i in range(n_qcol))
                qresid.append(b_vec[r] - H_pred)

            max_qresid = max(abs(r) for r in qresid)
            print(f"Quadratic max residual = {max_qresid:.6f}")
            print(f"Is H quadratic in e_j? {'YES' if max_qresid < 0.01 else 'NO'}")

    # Step 5: Product analysis
    print(f"\n{'='*70}")
    print(f"PRODUCT AND SPECIAL VALUES")
    print(f"{'='*70}")

    # Product of Q_k for each profile
    print(f"\nProducts (e_m = prod Q_k):")
    e_m_vals = set()
    for sigma in H_dict:
        e_m = round(all_esf[sigma][m])
        e_m_vals.add(e_m)
    for val in sorted(e_m_vals):
        count = sum(1 for sigma in H_dict if round(all_esf[sigma][m]) == val)
        print(f"  e_{m} = {val}, count = {count}")

    # H/p values
    print(f"\nH/p values:")
    H_vals = sorted(set(H_dict.values()), reverse=True)
    for H in H_vals:
        print(f"  H = {H}, H/p = {H/p:.2f}")

    # Interval and Paley
    sigma_int = tuple([1] * m)
    sigma_paley = tuple(1 if (i+1) in QR else -1 for i in range(m))
    H_int = H_dict[sigma_int]
    H_paley = H_dict.get(sigma_paley, "N/A")

    print(f"\nInterval: sigma={sigma_int}, H={H_int}")
    print(f"Paley: sigma={sigma_paley}, H={H_paley}")
    print(f"Interval is max: {H_int == max(H_dict.values())}")


if __name__ == '__main__':
    main()
