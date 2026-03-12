"""
esym_p13_analysis.py — Elementary symmetric polynomial decomposition at p=13.

AT p=13 (mod 4 = 1):
  - No Paley tournament (QR_13 not a valid tournament set)
  - The cyclic interval (Dirichlet kernel) maximizes H
  - Spectral landscape reverses: concentrated > flat
  - THM-134 (opus) proved local max at p=7,11 — what happens at p=13?

QUESTIONS:
  1. Does H = c_0 + sum c_k * e_k(x) still hold?
  2. What are the signs of c_k?
  3. Is the Hessian POSITIVE definite at the "center" (making it a LOCAL MIN)?
  4. Is the interval tournament a local max?

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import cmath
import math
from itertools import combinations
from collections import defaultdict

sys.path.insert(0, '04-computation')


def all_circulant_tournaments(n):
    pairs, used = [], set()
    for a in range(1, n):
        if a not in used:
            b = n - a
            if a == b: return []
            pairs.append((a, b)); used.add(a); used.add(b)
    results = []
    for bits in range(2 ** len(pairs)):
        S = [a if (bits >> i) & 1 else b for i, (a, b) in enumerate(pairs)]
        results.append(tuple(sorted(S)))
    return results


def ham_count_dp(n, S):
    S_set = set(S)
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i)%n in S_set: adj[i][j] = True
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if dp[mask][v]==0 or not(mask&(1<<v)): continue
            for w in range(n):
                if mask&(1<<w): continue
                if adj[v][w]: dp[mask|(1<<w)][w] += dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))


def circulant_eigenvalues(n, S):
    omega = cmath.exp(2j*cmath.pi/n)
    return [sum(omega**(k*s) for s in S) for k in range(n)]


def elementary_symmetric(xs, k):
    """Compute e_k(x_1, ..., x_m) = sum over k-subsets of product."""
    if k == 0: return 1
    if k > len(xs): return 0
    total = 0
    for subset in combinations(range(len(xs)), k):
        prod = 1
        for i in subset:
            prod *= xs[i]
        total += prod
    return total


def solve_linear_system(A, b):
    """Solve Ax = b by Gauss elimination. Returns x or None."""
    n = len(b)
    # Augmented matrix
    M = [row[:] + [bi] for row, bi in zip(A, b)]

    for col in range(n):
        # Pivot
        max_row = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[max_row][col]) < 1e-12:
            return None
        M[col], M[max_row] = M[max_row], M[col]
        for row in range(col+1, n):
            f = M[row][col] / M[col][col]
            for j in range(col, n+1):
                M[row][j] -= f * M[col][j]

    # Back substitution
    x = [0]*n
    for i in range(n-1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j]*x[j] for j in range(i+1, n))) / M[i][i]
    return x


def main():
    print("=" * 70)
    print("ELEMENTARY SYMMETRIC POLYNOMIAL ANALYSIS — p = 13")
    print("=" * 70)

    for p in [7, 11, 13]:
        print(f"\n{'=' * 60}")
        print(f"p = {p} (mod 4 = {p % 4}), m = {(p-1)//2}")
        print(f"{'=' * 60}")

        m = (p-1)//2
        all_S = all_circulant_tournaments(p)

        # Compute H and spectral data
        data = []
        for S in all_S:
            H = ham_count_dp(p, S)
            eigs = circulant_eigenvalues(p, S)
            y2 = [eigs[k].imag**2 for k in range(1, m+1)]
            # Elementary symmetric polynomials e_k(y^2)
            esyms = [elementary_symmetric(y2, k) for k in range(m+1)]
            data.append({'S': S, 'H': H, 'y2': y2, 'esyms': esyms})

        data.sort(key=lambda d: d['H'], reverse=True)
        H_vals = sorted(set(d['H'] for d in data), reverse=True)
        n_H = len(H_vals)

        print(f"  {len(all_S)} tournaments, {n_H} distinct H values")

        # Check e_1 universality
        e1_vals = set(round(d['esyms'][1], 8) for d in data)
        print(f"  e_1 universal: {len(e1_vals) == 1} (values: {e1_vals})")

        # Print table of H vs e_k
        header = f"  {'H':>12}" + "".join(f"  {'e_'+str(k):>12}" for k in range(1, min(m+1, 7)))
        print(f"\n{header}")
        for H_val in H_vals:
            d = next(x for x in data if x['H'] == H_val)
            row = f"  {H_val:>12}"
            for k in range(1, min(m+1, 7)):
                row += f"  {d['esyms'][k]:>12.4f}"
            count = sum(1 for x in data if x['H'] == H_val)
            # Mark special
            if H_val == H_vals[0]:
                row += "  MAX"
            if H_val == H_vals[-1]:
                row += "  MIN"
            row += f"  (x{count})"
            print(row)

        # Fit H = c_0 + c_2*e_2 + ... + c_m*e_m
        # (e_1 is constant, absorbed into c_0)
        # Need n_H >= m to have enough data points

        # Get one representative per H class
        representatives = []
        for H_val in H_vals:
            d = next(x for x in data if x['H'] == H_val)
            representatives.append(d)

        # Try fitting with increasing number of e_k
        for n_terms in range(1, min(n_H+1, m)):
            # Use e_2, ..., e_{n_terms+1}
            k_start = 2
            k_end = k_start + n_terms

            # Build system: for each H class, H = c_0 + sum c_k * e_k
            # Variables: c_0, c_{k_start}, ..., c_{k_end-1}
            n_vars = 1 + n_terms
            if n_H < n_vars:
                continue

            # Least squares if n_H > n_vars
            # Build design matrix
            X = []
            y_vec = []
            for d in representatives:
                row = [1.0]  # constant term
                for k in range(k_start, k_end):
                    row.append(d['esyms'][k])
                X.append(row)
                y_vec.append(float(d['H']))

            # Normal equations: X^T X c = X^T y
            XtX = [[sum(X[i][r]*X[i][c2] for i in range(n_H)) for c2 in range(n_vars)] for r in range(n_vars)]
            Xty = [sum(X[i][r]*y_vec[i] for i in range(n_H)) for r in range(n_vars)]

            coeffs = solve_linear_system(XtX, Xty)
            if coeffs is None:
                print(f"\n  Fit with e_{k_start},...,e_{k_end-1}: SINGULAR")
                continue

            # Compute residuals
            residuals = []
            for d in representatives:
                pred = coeffs[0]
                for j, k in enumerate(range(k_start, k_end)):
                    pred += coeffs[j+1] * d['esyms'][k]
                residuals.append(d['H'] - pred)

            max_res = max(abs(r) for r in residuals)
            ss_res = sum(r**2 for r in residuals)
            mean_H = sum(d['H'] for d in representatives) / n_H
            ss_tot = sum((d['H'] - mean_H)**2 for d in representatives)
            r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 1

            formula = f"H = {coeffs[0]:.4f}"
            signs = []
            for j, k in enumerate(range(k_start, k_end)):
                sign = "+" if coeffs[j+1] >= 0 else "-"
                formula += f" {sign}{abs(coeffs[j+1]):.6f}*e_{k}"
                signs.append(sign)

            print(f"\n  Fit: {formula}")
            print(f"  R^2 = {r2:.10f}, max |residual| = {max_res:.6f}")

            if max_res < 0.01:
                print(f"  *** EXACT FIT ***")
                pos = sum(1 for s in signs if s == "+")
                neg = sum(1 for s in signs if s == "-")
                print(f"  Signs: {pos} positive, {neg} negative")

                # Hessian at "center" (all x_k = p/4)
                v = p / 4.0
                # Hessian of e_k at uniform point, restricted to sum=0:
                # Hess(e_k)|_{sum=0} = -C(m-2, k-2) * v^{k-2} * I
                lambda_H = 0
                for j, k in enumerate(range(k_start, k_end)):
                    if k >= 2:
                        binom = math.comb(m-2, k-2)
                        contrib = coeffs[j+1] * binom * v**(k-2)
                        lambda_H -= contrib
                        print(f"    Hessian contrib from e_{k}: "
                              f"c={coeffs[j+1]:.4f}, C({m-2},{k-2})={binom}, "
                              f"v^{k-2}={v**(k-2):.4f}, term={-contrib:.4f}")

                print(f"  lambda_H = {lambda_H:.6f}")
                if lambda_H < 0:
                    print(f"  Hessian NEGATIVE DEFINITE at center => center is LOCAL MAX")
                else:
                    print(f"  Hessian POSITIVE DEFINITE at center => center is LOCAL MIN!")
                    print(f"  *** LANDSCAPE REVERSAL CONFIRMED ***")

                break  # Found exact fit

        # Check if center (uniform) is closest to max or min
        center_y2 = p / 4.0
        for d in data:
            d['dist_center'] = sum((y - center_y2)**2 for y in d['y2'])

        max_tournament = data[0]
        min_tournament = data[-1]
        print(f"\n  MAX H={max_tournament['H']}: dist from center = {max_tournament['dist_center']:.4f}")
        print(f"  MIN H={min_tournament['H']}: dist from center = {min_tournament['dist_center']:.4f}")

        # Verify: all e_k are maximized at uniform for p=3 mod 4 (Schur)
        if p % 4 == 3:
            print(f"\n  Schur concavity check: e_k(Paley) > e_k(others)?")
            paley_d = next(d for d in data if d['dist_center'] < 1e-6)
            for k in range(2, min(m+1, 7)):
                paley_ek = paley_d['esyms'][k]
                max_other = max(d['esyms'][k] for d in data if d['dist_center'] > 1e-6)
                print(f"    e_{k}: Paley={paley_ek:.4f}, max other={max_other:.4f}, "
                      f"Paley > all? {paley_ek > max_other + 1e-8}")


if __name__ == '__main__':
    main()
