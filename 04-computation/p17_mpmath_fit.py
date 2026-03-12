"""
p17_mpmath_fit.py — High-precision e_k fitting at p=17 using mpmath.

The float sin() values have ~15 digits of precision, which is insufficient
when computing e_k of 8 variables (products of 8 numbers, each ~O(1),
give values that lose precision at high k).

Strategy: use mpmath with 100+ decimal digits for y_k^2, e_k, and the
linear algebra. If the fit is truly exact, the residuals should drop to
machine epsilon of the mpmath precision.

Also test: is H affine in e_k, or does it need QUADRATIC terms?

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import math
import time
from collections import defaultdict
from itertools import combinations

try:
    import mpmath
    mpmath.mp.dps = 100  # 100 decimal digits
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False
    print("mpmath not available, falling back to standard float")

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
    adj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S_set:
                adj[i] |= (1 << j)
    full_mask = (1 << n) - 1
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        if bin(mask).count('1') >= n:
            continue
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp.get((mask, v), 0)
            if cnt == 0:
                continue
            candidates = adj[v] & ~mask
            w = 0
            while candidates:
                if candidates & 1:
                    dp[(mask | (1 << w), w)] += cnt
                candidates >>= 1
                w += 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))


def multiplicative_orbit(p, S):
    S_set = frozenset(S)
    orbit = {S_set}
    for a in range(2, p):
        if math.gcd(a, p) != 1:
            continue
        orbit.add(frozenset((a * s) % p for s in S))
    return orbit


def find_orbit_representatives(p):
    all_S = all_circulant_tournaments(p)
    seen = set()
    reps = []
    for S in all_S:
        S_frozen = frozenset(S)
        if S_frozen in seen:
            continue
        orbit = multiplicative_orbit(p, S)
        reps.append(tuple(sorted(S)))
        for member in orbit:
            seen.add(member)
    return reps


def y_squared_mpmath(p, S):
    """Compute y_k^2 with mpmath high precision."""
    m = (p - 1) // 2
    pi = mpmath.pi
    y2 = []
    for k in range(1, m + 1):
        val = sum(mpmath.sin(2 * pi * k * s / p) for s in S)
        y2.append(val * val)
    return y2


def elementary_symmetric_mpmath(xs, k):
    """Compute e_k with mpmath precision."""
    if k == 0: return mpmath.mpf(1)
    if k > len(xs): return mpmath.mpf(0)
    total = mpmath.mpf(0)
    for subset in combinations(range(len(xs)), k):
        prod = mpmath.mpf(1)
        for i in subset:
            prod *= xs[i]
        total += prod
    return total


def solve_mpmath(A, b):
    """Solve Ax=b with mpmath matrix operations."""
    n = len(b)
    M = mpmath.matrix(n, n + 1)
    for r in range(n):
        for c in range(n):
            M[r, c] = A[r][c]
        M[r, n] = b[r]

    # Gaussian elimination with partial pivoting
    for col in range(n):
        max_row = col
        for r in range(col + 1, n):
            if abs(M[r, col]) > abs(M[max_row, col]):
                max_row = r
        if abs(M[max_row, col]) < mpmath.mpf(10) ** (-80):
            return None
        # Swap
        for j in range(n + 1):
            M[col, j], M[max_row, j] = M[max_row, j], M[col, j]
        # Eliminate
        for row in range(col + 1, n):
            f = M[row, col] / M[col, col]
            for j in range(col, n + 1):
                M[row, j] -= f * M[col, j]

    # Back substitution
    x = [mpmath.mpf(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i, n] - sum(M[i, j] * x[j] for j in range(i + 1, n))) / M[i, i]
    return x


def main():
    if not HAS_MPMATH:
        print("ERROR: mpmath required")
        return

    p = 17
    m = (p - 1) // 2  # = 8
    print(f"{'=' * 70}")
    print(f"p = {p}, m = {m}: MPMATH HIGH-PRECISION e_k FIT")
    print(f"mpmath precision: {mpmath.mp.dps} decimal digits")
    print(f"{'=' * 70}")

    # Step 1: Compute H and high-precision y^2 for orbit reps
    reps = find_orbit_representatives(p)
    print(f"\n  {len(reps)} orbit representatives")

    data = []
    for i, S in enumerate(reps):
        t0 = time.time()
        H = ham_count_dp(p, S)
        elapsed = time.time() - t0
        y2 = y_squared_mpmath(p, S)
        esyms = [elementary_symmetric_mpmath(y2, k) for k in range(m + 1)]
        orbit_size = len(multiplicative_orbit(p, S))
        data.append({'S': S, 'H': H, 'y2': y2, 'esyms': esyms, 'orbit': orbit_size})
        print(f"  [{i+1}/{len(reps)}] H={H}, time={elapsed:.1f}s")

    data.sort(key=lambda d: d['H'], reverse=True)
    H_vals = sorted(set(d['H'] for d in data), reverse=True)
    n_H = len(H_vals)
    print(f"\n  {n_H} distinct H values")

    # Get one rep per H value
    reps_per_H = [next(d for d in data if d['H'] == H_val) for H_val in H_vals]

    # Step 2: e_1 check
    e1_vals = set(mpmath.nstr(d['esyms'][1], 20) for d in data)
    print(f"  e_1 values: {e1_vals} (should be single value = {(p-1)/2})")

    # Step 3: AFFINE FIT H = c_0 + sum c_k * e_k
    print(f"\n  === AFFINE FIT (e_2 through e_8) ===")
    for n_terms in range(1, m):
        k_start = 2
        k_end = k_start + n_terms
        n_vars = 1 + n_terms

        if n_H < n_vars:
            continue

        # Build design matrix
        X = []
        y_vec = []
        for d in reps_per_H:
            row = [mpmath.mpf(1)] + [d['esyms'][k] for k in range(k_start, k_end)]
            X.append(row)
            y_vec.append(mpmath.mpf(d['H']))

        # Normal equations: X^T X c = X^T y
        XtX = [[sum(X[i][r] * X[i][c] for i in range(n_H))
                for c in range(n_vars)] for r in range(n_vars)]
        Xty = [sum(X[i][r] * y_vec[i] for i in range(n_H)) for r in range(n_vars)]

        coeffs = solve_mpmath(XtX, Xty)
        if coeffs is None:
            print(f"  n_terms={n_terms}: SINGULAR")
            continue

        # Compute residuals
        residuals = []
        for d in reps_per_H:
            pred = coeffs[0]
            for j, k in enumerate(range(k_start, k_end)):
                pred += coeffs[j + 1] * d['esyms'][k]
            residuals.append(mpmath.mpf(d['H']) - pred)

        max_res = max(abs(r) for r in residuals)
        ss_res = sum(r**2 for r in residuals)
        mean_H = sum(mpmath.mpf(d['H']) for d in reps_per_H) / n_H
        ss_tot = sum((mpmath.mpf(d['H']) - mean_H)**2 for d in reps_per_H)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1

        print(f"  n_terms={n_terms} (e_2..e_{k_end-1}): R²={mpmath.nstr(r2, 15)}, "
              f"max|res|={mpmath.nstr(max_res, 10)}")

        if max_res < mpmath.mpf('0.001'):
            print(f"\n  *** EXACT FIT with {n_terms} terms! ***")
            for j in range(n_vars):
                k = 0 if j == 0 else k_start + j - 1
                print(f"    c_{k} = {mpmath.nstr(coeffs[j], 20)}")

            # Hessian eigenvalue
            v = mpmath.mpf(p) / 4
            lambda_H = mpmath.mpf(0)
            for j, k in enumerate(range(k_start, k_end)):
                if k >= 2:
                    binom = math.comb(m - 2, k - 2)
                    lambda_H -= coeffs[j + 1] * binom * v**(k - 2)
            print(f"\n  lambda_H = {mpmath.nstr(lambda_H, 15)}")
            if lambda_H > 0:
                print(f"  *** POSITIVE — center is LOCAL MIN ***")
            else:
                print(f"  *** NEGATIVE — center is LOCAL MAX ***")
            break

    # Step 4: If affine didn't work, try QUADRATIC terms
    print(f"\n  === QUADRATIC FIT ===")
    print(f"  Testing H = c_0 + sum c_k*e_k + sum c_kl*e_k*e_l")

    # Use e_2, e_3, e_4 + quadratic terms
    # Terms: 1, e_2, e_3, e_4, e_2^2, e_2*e_3, e_2*e_4, e_3^2, e_3*e_4, e_4^2
    # = 10 terms (need at least 10 data points)
    terms_labels = ['1', 'e2', 'e3', 'e4', 'e2²', 'e2·e3', 'e2·e4', 'e3²', 'e3·e4', 'e4²']

    X_quad = []
    y_quad = []
    for d in reps_per_H:
        e2, e3, e4 = d['esyms'][2], d['esyms'][3], d['esyms'][4]
        row = [mpmath.mpf(1), e2, e3, e4, e2*e2, e2*e3, e2*e4, e3*e3, e3*e4, e4*e4]
        X_quad.append(row)
        y_quad.append(mpmath.mpf(d['H']))

    n_terms_q = len(terms_labels)
    XtX_q = [[sum(X_quad[i][r] * X_quad[i][c] for i in range(n_H))
              for c in range(n_terms_q)] for r in range(n_terms_q)]
    Xty_q = [sum(X_quad[i][r] * y_quad[i] for i in range(n_H)) for r in range(n_terms_q)]

    coeffs_q = solve_mpmath(XtX_q, Xty_q)
    if coeffs_q:
        residuals_q = []
        for d in reps_per_H:
            e2, e3, e4 = d['esyms'][2], d['esyms'][3], d['esyms'][4]
            row = [mpmath.mpf(1), e2, e3, e4, e2*e2, e2*e3, e2*e4, e3*e3, e3*e4, e4*e4]
            pred = sum(coeffs_q[j] * row[j] for j in range(n_terms_q))
            residuals_q.append(mpmath.mpf(d['H']) - pred)

        max_res_q = max(abs(r) for r in residuals_q)
        ss_res_q = sum(r**2 for r in residuals_q)
        r2_q = 1 - ss_res_q / ss_tot if ss_tot > 0 else 1
        print(f"  Quad(e2,e3,e4): R²={mpmath.nstr(r2_q, 15)}, max|res|={mpmath.nstr(max_res_q, 10)}")

        if max_res_q < mpmath.mpf('0.001'):
            print(f"  *** EXACT QUADRATIC FIT ***")
            for j in range(n_terms_q):
                print(f"    {terms_labels[j]:>8}: {mpmath.nstr(coeffs_q[j], 15)}")

    # Step 5: Try all e_k + e_2^2 (quadratic in e_2 only)
    print(f"\n  === MIXED: AFFINE e_2..e_8 + QUADRATIC e_2^2 ===")
    X_mix = []
    for d in reps_per_H:
        row = [mpmath.mpf(1)]
        for k in range(2, m + 1):
            row.append(d['esyms'][k])
        row.append(d['esyms'][2] ** 2)  # e_2^2
        X_mix.append(row)

    n_vars_mix = 1 + (m - 1) + 1  # = 9
    XtX_m = [[sum(X_mix[i][r] * X_mix[i][c] for i in range(n_H))
              for c in range(n_vars_mix)] for r in range(n_vars_mix)]
    Xty_m = [sum(X_mix[i][r] * mpmath.mpf(reps_per_H[i]['H']) for i in range(n_H))
             for r in range(n_vars_mix)]

    coeffs_m = solve_mpmath(XtX_m, Xty_m)
    if coeffs_m:
        residuals_m = []
        for d in reps_per_H:
            row = [mpmath.mpf(1)]
            for k in range(2, m + 1):
                row.append(d['esyms'][k])
            row.append(d['esyms'][2] ** 2)
            pred = sum(coeffs_m[j] * row[j] for j in range(n_vars_mix))
            residuals_m.append(mpmath.mpf(d['H']) - pred)

        max_res_m = max(abs(r) for r in residuals_m)
        ss_res_m = sum(r**2 for r in residuals_m)
        r2_m = 1 - ss_res_m / ss_tot if ss_tot > 0 else 1
        print(f"  R²={mpmath.nstr(r2_m, 15)}, max|res|={mpmath.nstr(max_res_m, 10)}")

    # Step 6: Verify smaller primes with same method
    print(f"\n  === VERIFICATION: p=7, 11, 13 with mpmath ===")
    for p_test in [7, 11, 13]:
        m_test = (p_test - 1) // 2
        all_S_test = all_circulant_tournaments(p_test)
        data_test = []
        for S in all_S_test:
            H = ham_count_dp(p_test, S)
            y2 = y_squared_mpmath(p_test, S)
            esyms = [elementary_symmetric_mpmath(y2, k) for k in range(m_test + 1)]
            data_test.append({'H': H, 'esyms': esyms})

        data_test.sort(key=lambda d: d['H'], reverse=True)
        H_vals_test = sorted(set(d['H'] for d in data_test), reverse=True)
        n_H_test = len(H_vals_test)
        reps_test = [next(d for d in data_test if d['H'] == H_val) for H_val in H_vals_test]

        # Fit affine
        for n_terms_test in range(1, min(n_H_test, m_test)):
            k_start_t = 2
            k_end_t = k_start_t + n_terms_test
            n_vars_t = 1 + n_terms_test
            if n_H_test < n_vars_t:
                continue

            X_t = [[mpmath.mpf(1)] + [d['esyms'][k] for k in range(k_start_t, k_end_t)]
                   for d in reps_test]
            y_t = [mpmath.mpf(d['H']) for d in reps_test]

            XtX_t = [[sum(X_t[i][r] * X_t[i][c] for i in range(n_H_test))
                      for c in range(n_vars_t)] for r in range(n_vars_t)]
            Xty_t = [sum(X_t[i][r] * y_t[i] for i in range(n_H_test)) for r in range(n_vars_t)]

            coeffs_t = solve_mpmath(XtX_t, Xty_t)
            if coeffs_t is None:
                continue

            residuals_t = []
            for d in reps_test:
                pred = coeffs_t[0]
                for j, k in enumerate(range(k_start_t, k_end_t)):
                    pred += coeffs_t[j + 1] * d['esyms'][k]
                residuals_t.append(mpmath.mpf(d['H']) - pred)

            max_res_t = max(abs(r) for r in residuals_t)

            if max_res_t < mpmath.mpf('0.001'):
                print(f"  p={p_test}: EXACT with {n_terms_test} terms (e_2..e_{k_end_t-1}), "
                      f"max|res|={mpmath.nstr(max_res_t, 5)}")
                break
            elif n_terms_test == min(n_H_test, m_test) - 1:
                print(f"  p={p_test}: NOT exact with {n_terms_test} terms, "
                      f"max|res|={mpmath.nstr(max_res_t, 10)}")


if __name__ == '__main__':
    main()
