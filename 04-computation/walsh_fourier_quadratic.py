#!/usr/bin/env python3
"""
walsh_fourier_quadratic.py -- The quadratic form structure of H(sigma)

KEY DISCOVERY:
For circulant tournament on Z_p with chord pairs (g_i, p-g_i), g_i = i+1:

  S_hat(k; sigma) = sum_i f_i(k, sigma_i)

where f_i(k, +1) = omega^{k*g_i} and f_i(k, -1) = omega^{-k*g_i}.

Decomposition:
  f_i(k, sigma_i) = cos(2*pi*k*g_i/p) + i*sigma_i*sin(2*pi*k*g_i/p)

Therefore:
  S_hat(k; sigma) = C(k) + i * D(k; sigma)

where:
  C(k) = sum_i cos(2*pi*k*g_i/p)              [sigma-INDEPENDENT]
  D(k; sigma) = sum_i sigma_i * sin(2*pi*k*g_i/p)   [sigma-LINEAR]

And:
  |S_hat(k; sigma)|^2 = C(k)^2 + D(k; sigma)^2

Since D is linear in sigma, |S_hat|^2 is a QUADRATIC FORM in sigma.

CONSEQUENCES:
1. H depends on sigma through |S_hat(k)|^2 values (transfer matrix structure)
2. Each |S_hat(k)|^2 is a degree-2 polynomial in sigma
3. H = F(Q_1(sigma), ..., Q_m(sigma)) where Q_k are quadratic forms
4. Degree-2n Walsh terms come from n-th order interactions among Q_k
5. Walsh odd-degree vanishing is AUTOMATIC (each Q_k is even-degree)

This script:
1. Verifies the C/D decomposition
2. Expresses each |S_hat(k)|^2 as a quadratic form in sigma
3. Tests if H = linear function of {|S_hat(k)|^2} (which would give exact degree-2 Walsh)
4. Tests higher-order dependence (which generates degree-4+ Walsh terms)

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
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
    print("=" * 70)
    print("WALSH-FOURIER QUADRATIC FORM ANALYSIS")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}, p mod 4 = {p%4}")
        print(f"{'='*70}")

        # ====== 1. C(k) and sin coefficients ======
        print(f"\n  1. FOURIER DECOMPOSITION STRUCTURE")
        print(f"     S_hat(k) = C(k) + i*D(k;sigma), C sigma-indep, D linear in sigma")

        # g_i = i+1 for i=0,...,m-1
        gaps = [i + 1 for i in range(m)]

        # C(k) = sum cos(2*pi*k*g_i/p)
        C_k = []
        # sin_coeff[k][i] = sin(2*pi*k*g_i/p)
        sin_coeffs = []

        for k in range(1, m + 1):
            c = sum(math.cos(2 * math.pi * k * g / p) for g in gaps)
            C_k.append(c)
            sines = [math.sin(2 * math.pi * k * g / p) for g in gaps]
            sin_coeffs.append(sines)

        print(f"\n    C(k) values:")
        for k in range(m):
            print(f"      C({k+1}) = {C_k[k]:>+.6f}")

        print(f"\n    Sine coefficient matrix sin(2*pi*k*g_i/p):")
        header = "k\\i  " + " ".join(f"   i={i}" for i in range(m))
        print(f"      {header}")
        for k in range(m):
            row = f"  k={k+1}  " + " ".join(f"{sin_coeffs[k][i]:>+.4f}" for i in range(m))
            print(f"      {row}")

        # ====== 2. Verify |S_hat(k)|^2 = C^2 + D^2 ======
        print(f"\n  2. VERIFY |S_hat(k)|^2 = C(k)^2 + D(k;sigma)^2")

        pairs = [(s, p - s) for s in range(1, m + 1)]
        H_dict = {}
        shat_sq_dict = {}

        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = []
            for i, (a, b) in enumerate(pairs):
                S.append(a if sigma[i] == 1 else b)
            S = sorted(S)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            H_dict[sigma] = H

            # Compute |S_hat(k)|^2 directly
            shat_sq = []
            for k in range(1, m + 1):
                val = sum(omega ** (k * s) for s in S)
                shat_sq.append(abs(val)**2)
            shat_sq_dict[sigma] = shat_sq

            # Compute via C + iD
            shat_sq_formula = []
            for ki in range(m):
                D = sum(sigma[i] * sin_coeffs[ki][i] for i in range(m))
                sq = C_k[ki]**2 + D**2
                shat_sq_formula.append(sq)

            # Verify
            max_err = max(abs(shat_sq[ki] - shat_sq_formula[ki])
                         for ki in range(m))
            if bits < 4:  # show first few
                print(f"    sigma={sigma}: max error = {max_err:.1e}")

        all_errors = []
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            shat_sq = shat_sq_dict[sigma]
            for ki in range(m):
                D = sum(sigma[i] * sin_coeffs[ki][i] for i in range(m))
                sq = C_k[ki]**2 + D**2
                all_errors.append(abs(shat_sq[ki] - sq))
        print(f"    Max error over all sigma: {max(all_errors):.2e}")

        # ====== 3. Express Q_k(sigma) = |S_hat(k)|^2 as explicit quadratic ======
        print(f"\n  3. QUADRATIC FORM COEFFICIENTS")
        print(f"     Q_k(sigma) = C(k)^2 + (sum_i s_ki * sigma_i)^2")
        print(f"     = C(k)^2 + sum_i s_ki^2 + sum_{{i<j}} 2*s_ki*s_kj * sigma_i*sigma_j")
        print(f"     (constant + pure quadratic: NO linear terms)")

        Q_const = []  # constant part of Q_k: C(k)^2 + sum s_ki^2
        Q_coeff = []  # Q_coeff[k][(i,j)] = 2*s_ki*s_kj for i<j

        for ki in range(m):
            const = C_k[ki]**2 + sum(sin_coeffs[ki][i]**2 for i in range(m))
            coeffs = {}
            for i in range(m):
                for j in range(i+1, m):
                    coeffs[(i,j)] = 2 * sin_coeffs[ki][i] * sin_coeffs[ki][j]
            Q_const.append(const)
            Q_coeff.append(coeffs)

        print(f"\n    Quadratic form constants (sigma-independent part):")
        for ki in range(m):
            print(f"      Q_{ki+1} constant = {Q_const[ki]:.6f}")

        # The sum over all k gives the total quadratic form
        total_const = sum(Q_const)
        total_coeff = {}
        for i in range(m):
            for j in range(i+1, m):
                total_coeff[(i,j)] = sum(Q_coeff[ki][(i,j)] for ki in range(m))

        print(f"\n    Total quadratic form (sum of all Q_k):")
        print(f"      Constant = {total_const:.6f}")
        for (i,j), c in sorted(total_coeff.items()):
            if abs(c) > 1e-10:
                print(f"      sigma_{i}*sigma_{j}: {c:>+.6f}")
        if all(abs(c) < 1e-10 for c in total_coeff.values()):
            print(f"      ALL cross terms = 0 (orthogonal decomposition!)")

        # ====== 4. Test H as function of {Q_k} ======
        print(f"\n  4. IS H A LINEAR FUNCTION OF {{Q_k}}?")
        print(f"     Testing H(sigma) = a_0 + sum_k a_k * Q_k(sigma)")

        # Build matrix equation: for each sigma, compute Q values
        n_sigma = 1 << m
        Q_matrix = []
        H_vector = []
        for bits in range(n_sigma):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            Q_vals = []
            for ki in range(m):
                D = sum(sigma[i] * sin_coeffs[ki][i] for i in range(m))
                Q_vals.append(C_k[ki]**2 + D**2)
            Q_matrix.append([1.0] + Q_vals)  # intercept + Q_1, ..., Q_m
            H_vector.append(H_dict[sigma])

        # Solve least squares using normal equations
        # X^T X a = X^T H
        n_params = m + 1
        XTX = [[0.0]*n_params for _ in range(n_params)]
        XTH = [0.0]*n_params
        for row_idx in range(n_sigma):
            for i in range(n_params):
                for j in range(n_params):
                    XTX[i][j] += Q_matrix[row_idx][i] * Q_matrix[row_idx][j]
                XTH[i] += Q_matrix[row_idx][i] * H_vector[row_idx]

        # Solve via Gaussian elimination
        aug = [XTX[i][:] + [XTH[i]] for i in range(n_params)]
        for col in range(n_params):
            # Find pivot
            max_row = max(range(col, n_params), key=lambda r: abs(aug[r][col]))
            aug[col], aug[max_row] = aug[max_row], aug[col]
            pivot = aug[col][col]
            if abs(pivot) < 1e-15:
                print(f"    Singular matrix at column {col}")
                continue
            for j in range(col, n_params + 1):
                aug[col][j] /= pivot
            for r in range(n_params):
                if r == col:
                    continue
                factor = aug[r][col]
                for j in range(col, n_params + 1):
                    aug[r][j] -= factor * aug[col][j]

        coeffs_linear = [aug[i][n_params] for i in range(n_params)]

        # Check residuals
        residuals = []
        for row_idx in range(n_sigma):
            H_pred = sum(coeffs_linear[i] * Q_matrix[row_idx][i] for i in range(n_params))
            residuals.append(H_vector[row_idx] - H_pred)

        max_resid = max(abs(r) for r in residuals)
        rms_resid = (sum(r**2 for r in residuals) / n_sigma) ** 0.5
        print(f"    Max residual = {max_resid:.6f}")
        print(f"    RMS residual = {rms_resid:.6f}")
        print(f"    Is H linear in Q_k? {'YES' if max_resid < 0.01 else 'NO'}")

        if max_resid < 0.01:
            print(f"    Coefficients:")
            print(f"      a_0 = {coeffs_linear[0]:.6f}")
            for ki in range(m):
                print(f"      a_{ki+1} = {coeffs_linear[ki+1]:.6f}")
        else:
            print(f"    Non-linear dependence detected!")
            print(f"    Top 5 residuals:")
            sigma_list = []
            for bits in range(n_sigma):
                sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
                sigma_list.append(sigma)
            indexed = sorted(enumerate(residuals), key=lambda x: -abs(x[1]))[:5]
            for idx, resid in indexed:
                print(f"      sigma={sigma_list[idx]}: H={H_vector[idx]}, "
                      f"H_pred={H_vector[idx]-resid:.1f}, residual={resid:>+.2f}")

        # ====== 5. Test H as quadratic function of {Q_k} ======
        print(f"\n  5. IS H A QUADRATIC FUNCTION OF {{Q_k}}?")
        print(f"     Testing H = a_0 + sum a_k Q_k + sum a_kl Q_k Q_l")

        # Add quadratic terms Q_k * Q_l (k <= l)
        Q_matrix_quad = []
        for bits in range(n_sigma):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            Q_vals = []
            for ki in range(m):
                D = sum(sigma[i] * sin_coeffs[ki][i] for i in range(m))
                Q_vals.append(C_k[ki]**2 + D**2)

            row = [1.0] + Q_vals  # constant + linear
            for ki in range(m):
                for kj in range(ki, m):
                    row.append(Q_vals[ki] * Q_vals[kj])
            Q_matrix_quad.append(row)

        n_params_quad = len(Q_matrix_quad[0])
        print(f"    Number of parameters: {n_params_quad}")

        # Check if we have enough equations
        if n_params_quad > n_sigma:
            print(f"    Underdetermined ({n_params_quad} params > {n_sigma} equations)")
        else:
            # Solve via least squares
            XTX = [[0.0]*n_params_quad for _ in range(n_params_quad)]
            XTH = [0.0]*n_params_quad
            for row_idx in range(n_sigma):
                for i in range(n_params_quad):
                    for j in range(n_params_quad):
                        XTX[i][j] += Q_matrix_quad[row_idx][i] * Q_matrix_quad[row_idx][j]
                    XTH[i] += Q_matrix_quad[row_idx][i] * H_vector[row_idx]

            aug = [XTX[i][:] + [XTH[i]] for i in range(n_params_quad)]
            for col in range(n_params_quad):
                max_row = max(range(col, n_params_quad), key=lambda r: abs(aug[r][col]))
                aug[col], aug[max_row] = aug[max_row], aug[col]
                pivot = aug[col][col]
                if abs(pivot) < 1e-15:
                    continue
                for j in range(col, n_params_quad + 1):
                    aug[col][j] /= pivot
                for r in range(n_params_quad):
                    if r == col:
                        continue
                    factor = aug[r][col]
                    for j in range(col, n_params_quad + 1):
                        aug[r][j] -= factor * aug[col][j]

            coeffs_quad = [aug[i][n_params_quad] for i in range(n_params_quad)]

            residuals_q = []
            for row_idx in range(n_sigma):
                H_pred = sum(coeffs_quad[i] * Q_matrix_quad[row_idx][i]
                            for i in range(n_params_quad))
                residuals_q.append(H_vector[row_idx] - H_pred)

            max_resid_q = max(abs(r) for r in residuals_q)
            rms_resid_q = (sum(r**2 for r in residuals_q) / n_sigma) ** 0.5
            print(f"    Max residual = {max_resid_q:.6f}")
            print(f"    RMS residual = {rms_resid_q:.6f}")
            print(f"    Is H quadratic in Q_k? {'YES' if max_resid_q < 0.01 else 'NO'}")

        # ====== 6. CONNECTION TO WALSH DECOMPOSITION ======
        print(f"\n  6. WALSH DECOMPOSITION FROM QUADRATIC FORMS")
        print(f"     Each Q_k contributes degree-2 Walsh coefficients via s_ki*s_kj terms")

        # The Walsh decomposition of Q_k(sigma):
        # Q_k = const + sum_{i<j} 2*s_ki*s_kj * sigma_i*sigma_j
        # The Walsh coefficient of Q_k at subset {i,j} is 2*s_ki*s_kj.

        # If H = a_0 + sum a_k * Q_k, then:
        # h_hat[{i,j}] = sum_k a_k * 2 * s_ki * s_kj

        # Compute the predicted Walsh degree-2 coefficients
        # Full Walsh decomposition for comparison
        n = 1 << m
        h_hat = {}
        for bits in range(n):
            S_idx = tuple(i for i in range(m) if bits & (1 << i))
            coeff = 0
            for sigma_bits in range(n):
                sigma = tuple((1 if sigma_bits & (1 << i) else -1) for i in range(m))
                H = H_dict[sigma]
                prod = 1
                for i in S_idx:
                    prod *= sigma[i]
                coeff += H * prod
            h_hat[S_idx] = coeff / n

        print(f"\n    Walsh degree-2 vs quadratic form prediction:")
        for i in range(m):
            for j in range(i+1, m):
                actual = h_hat.get((i,j), 0)
                if abs(actual) < 0.01:
                    continue

                # Predicted from linear model: h_pred = sum_k a_k * 2*s_ki*s_kj
                if max_resid < 0.01:  # linear model works
                    pred = sum(coeffs_linear[ki+1] * 2 * sin_coeffs[ki][i] * sin_coeffs[ki][j]
                              for ki in range(m))
                    print(f"      h_hat[{{{i},{j}}}] = {actual:>+.4f}, "
                          f"predicted = {pred:>+.4f}, "
                          f"error = {abs(actual-pred):.4f}")
                else:
                    print(f"      h_hat[{{{i},{j}}}] = {actual:>+.4f} (linear model fails)")

        # ====== 7. ORTHOGONALITY of {Q_k} ======
        print(f"\n  7. ORTHOGONALITY OF QUADRATIC FORMS")
        # <Q_k, Q_l> = (1/2^m) sum_sigma Q_k(sigma) Q_l(sigma)
        # If the Q_k are "orthogonal" in some sense, this simplifies the analysis.

        Q_arrays = []
        for ki in range(m):
            Q_vals = []
            for bits in range(n_sigma):
                sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
                D = sum(sigma[i] * sin_coeffs[ki][i] for i in range(m))
                Q_vals.append(C_k[ki]**2 + D**2)
            Q_arrays.append(Q_vals)

        # Compute Gram matrix <Q_k, Q_l> (after centering)
        Q_means = [sum(Q_arrays[ki]) / n_sigma for ki in range(m)]
        print(f"\n    Q_k means: {[f'{q:.4f}' for q in Q_means]}")

        gram = [[0.0]*m for _ in range(m)]
        for ki in range(m):
            for kj in range(ki, m):
                val = sum((Q_arrays[ki][s] - Q_means[ki]) * (Q_arrays[kj][s] - Q_means[kj])
                         for s in range(n_sigma)) / n_sigma
                gram[ki][kj] = val
                gram[kj][ki] = val

        print(f"\n    Gram matrix (centered):")
        for ki in range(m):
            row = " ".join(f"{gram[ki][kj]:>12.2f}" for kj in range(m))
            print(f"      {row}")

        # Check: are Q_k paired? (Q_k and Q_{p-k} should be related by complex conjugation)
        print(f"\n    Q_k pairing structure (k and p-k):")
        for ki in range(m):
            k = ki + 1
            k_conj = p - k
            if k_conj > m:
                k_conj_idx = p - k_conj - 1  # canonical
            else:
                k_conj_idx = k_conj - 1

            # Check if Q_k == Q_{p-k}
            diff = sum(abs(Q_arrays[ki][s] - Q_arrays[k_conj_idx][s])
                      for s in range(n_sigma))
            if diff < 1e-6:
                print(f"      Q_{k} == Q_{k_conj_idx+1} (paired)")
            else:
                print(f"      Q_{k} != Q_{k_conj_idx+1} (diff = {diff:.2f})")

        # ====== 8. CHIRALITY: Does H depend on sign of S_hat(k)? ======
        print(f"\n  8. CHIRALITY: DOES H DEPEND ON PHASE OF S_hat(k)?")
        # Group sigmas by their sorted {Q_k} multiset
        q_groups = defaultdict(list)
        for bits in range(n_sigma):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            Q_vals = tuple(round(Q_arrays[ki][bits], 6) for ki in range(m))
            q_sorted = tuple(sorted(Q_vals))
            q_groups[q_sorted].append((sigma, H_dict[sigma]))

        n_chiral = 0
        for q_sorted, items in q_groups.items():
            H_set = set(H for _, H in items)
            if len(H_set) > 1:
                n_chiral += 1
                if n_chiral <= 3:
                    print(f"    CHIRAL: Q_sorted={q_sorted[:3]}..., "
                          f"H values = {sorted(H_set)}")

        print(f"    Total chiral Q-multisets: {n_chiral} / {len(q_groups)}")
        if n_chiral == 0:
            print(f"    => H is determined by {'{'}Q_k{'}'} multiset (no chirality)")


if __name__ == '__main__':
    main()
