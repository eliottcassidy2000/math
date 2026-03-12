"""
interpolation_artifact.py — Verify that affine H(e_k) at p<=13 is interpolation.

KEY INSIGHT: At p=7,11,13, the number of Z_p* orbits equals the number
of affine parameters (m-1+1). So the "exact fit" is just solving a square
system — ANY set of m linearly independent features would give R^2=1.

TEST: Fit H to RANDOM linear combinations of y^2 values. If R^2=1.0
even with random features (as long as we use the right number), this
confirms the interpolation artifact.

Then: investigate what DOES characterize H-maximization robustly,
regardless of the polynomial basis.

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import math
import random
from collections import defaultdict
from itertools import combinations

import mpmath
mpmath.mp.dps = 50

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


def solve_mpmath(A, b):
    n = len(b)
    M = mpmath.matrix(n, n + 1)
    for r in range(n):
        for c in range(n):
            M[r, c] = A[r][c]
        M[r, n] = b[r]
    for col in range(n):
        max_row = col
        for r in range(col + 1, n):
            if abs(M[r, col]) > abs(M[max_row, col]):
                max_row = r
        if abs(M[max_row, col]) < mpmath.mpf(10) ** (-40):
            return None
        for j in range(n + 1):
            M[col, j], M[max_row, j] = M[max_row, j], M[col, j]
        for row in range(col + 1, n):
            f = M[row, col] / M[col, col]
            for j in range(col, n + 1):
                M[row, j] -= f * M[col, j]
    x = [mpmath.mpf(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i, n] - sum(M[i, j] * x[j] for j in range(i + 1, n))) / M[i, i]
    return x


def main():
    print("=" * 70)
    print("INTERPOLATION ARTIFACT TEST")
    print("=" * 70)

    random.seed(42)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        all_S = all_circulant_tournaments(p)

        # Compute H and y^2 for all connection sets
        data = []
        for S in all_S:
            H = ham_count_dp(p, S)
            pi = mpmath.pi
            y2 = []
            for k in range(1, m + 1):
                val = sum(mpmath.sin(2 * pi * k * s / p) for s in S)
                y2.append(val * val)
            data.append({'S': S, 'H': H, 'y2': y2})

        # Get orbit representatives
        H_vals = sorted(set(d['H'] for d in data), reverse=True)
        n_H = len(H_vals)
        reps = [next(d for d in data if d['H'] == hv) for hv in H_vals]

        print(f"\n  p={p}: m={m}, {len(all_S)} tours, {n_H} distinct H values")

        # Test 1: Fit H to RANDOM features of y^2
        # Use n_H-1 random linear combinations of y_k^2 as features
        n_feats = n_H - 1  # So total variables = n_H, exactly determined

        successes = 0
        n_trials = 20
        for trial in range(n_trials):
            # Generate random feature vectors
            rand_coeffs = [[random.gauss(0, 1) for _ in range(m)]
                          for _ in range(n_feats)]

            # Compute features for each representative
            X = []
            y_vec = []
            for d in reps:
                row = [mpmath.mpf(1)]  # constant
                for fc in rand_coeffs:
                    feat = sum(mpmath.mpf(fc[j]) * d['y2'][j] for j in range(m))
                    row.append(feat)
                X.append(row)
                y_vec.append(mpmath.mpf(d['H']))

            n_vars = len(X[0])
            # If exactly determined (n_H == n_vars), solve directly
            if n_H == n_vars:
                A = [[X[i][j] for j in range(n_vars)] for i in range(n_H)]
                b = y_vec
                coeffs = solve_mpmath(A, b)
                if coeffs is not None:
                    # Check residuals
                    residuals = []
                    for i in range(n_H):
                        pred = sum(coeffs[j] * X[i][j] for j in range(n_vars))
                        residuals.append(abs(y_vec[i] - pred))
                    max_res = float(max(residuals))
                    if max_res < 1e-30:
                        successes += 1

        print(f"    Random features ({n_feats}+1 params, {n_H} points): "
              f"{successes}/{n_trials} exact (all should be exact)")

        if successes == n_trials:
            print(f"    *** CONFIRMED: Any {n_feats} features give exact fit ***")
            print(f"    The 'affine in e_k' property is just interpolation!")

        # Test 2: Use FEWER features than n_H-1
        # If H is truly degree-1 in e_k, then using e_2 alone should give
        # high R^2 even when overdetermined
        for n_test in range(1, n_H):
            # Use first n_test random features
            if n_test >= n_H:
                break
            rand_coeffs_test = [[random.gauss(0, 1) for _ in range(m)]
                               for _ in range(n_test)]

            X = []
            for d in reps:
                row = [mpmath.mpf(1)]
                for fc in rand_coeffs_test:
                    feat = sum(mpmath.mpf(fc[j]) * d['y2'][j] for j in range(m))
                    row.append(feat)
                X.append(row)
            y_vec = [mpmath.mpf(d['H']) for d in reps]

            n_vars = len(X[0])
            # Least squares via normal equations
            XtX = [[sum(X[i][r] * X[i][c] for i in range(n_H))
                    for c in range(n_vars)] for r in range(n_vars)]
            Xty = [sum(X[i][r] * y_vec[i] for i in range(n_H))
                   for r in range(n_vars)]
            coeffs = solve_mpmath(XtX, Xty)
            if coeffs is None:
                continue

            residuals = [y_vec[i] - sum(coeffs[j] * X[i][j] for j in range(n_vars))
                        for i in range(n_H)]
            max_res = float(max(abs(r) for r in residuals))
            ss_res = float(sum(r**2 for r in residuals))
            mean_H = float(sum(y_vec)) / n_H
            ss_tot = sum((float(y) - mean_H)**2 for y in y_vec)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1

            if n_test <= 3 or max_res < 1e-10 or n_test == n_H - 1:
                print(f"    {n_test} random feats + const ({n_test+1} params, "
                      f"{n_H} data): R^2={r2:.6f}, max|res|={max_res:.4f}")

    # ================================================================
    # Now the REAL question: what's a ROBUST characterization of H-max?
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"ROBUST CHARACTERIZATION: What maximizes H?")
    print(f"{'=' * 60}")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        all_S = all_circulant_tournaments(p)

        data = []
        for S in all_S:
            H = ham_count_dp(p, S)
            pi = mpmath.pi
            y2 = []
            for k in range(1, m + 1):
                val = sum(mpmath.sin(2 * pi * k * s / p) for s in S)
                y2.append(val * val)

            # Spectral statistics
            y2_float = [float(y) for y in y2]
            mean_y2 = sum(y2_float) / m
            var_y2 = sum((y - mean_y2)**2 for y in y2_float) / m
            max_y2 = max(y2_float)
            min_y2 = min(y2_float)
            spread = max_y2 - min_y2

            # How "flat" is the spectrum?
            # Flat = all y_k^2 equal = p/4
            dist_flat = sum((y - p/4)**2 for y in y2_float)

            # QR check
            qr = set(pow(a, 2, p) % p for a in range(1, p)) - {0}
            is_paley = set(S) == qr
            is_interval = set(S) == set(range(p - m, p))

            data.append({
                'S': S, 'H': H, 'y2': y2_float,
                'var_y2': var_y2, 'spread': spread,
                'dist_flat': dist_flat, 'is_paley': is_paley,
                'is_interval': is_interval
            })

        data.sort(key=lambda d: d['H'], reverse=True)
        H_vals = sorted(set(d['H'] for d in data), reverse=True)
        reps = [next(d for d in data if d['H'] == hv) for hv in H_vals]

        print(f"\n  p={p} (mod 4 = {p%4}), m={m}:")
        print(f"  {'H':>10} {'var_y2':>10} {'spread':>8} {'dist_flat':>10} {'paley':>6} {'interval':>8}")
        for d in reps:
            print(f"  {d['H']:>10} {d['var_y2']:>10.4f} {d['spread']:>8.4f} "
                  f"{d['dist_flat']:>10.4f} {str(d['is_paley']):>6} "
                  f"{str(d['is_interval']):>8}")

        # Correlation: H vs var_y2
        H_list = [float(d['H']) for d in reps]
        var_list = [d['var_y2'] for d in reps]
        n = len(reps)
        mean_H = sum(H_list) / n
        mean_var = sum(var_list) / n
        cov = sum((H_list[i] - mean_H) * (var_list[i] - mean_var) for i in range(n))
        std_H = (sum((h - mean_H)**2 for h in H_list))**0.5
        std_var = (sum((v - mean_var)**2 for v in var_list))**0.5
        corr = cov / (std_H * std_var) if std_H * std_var > 0 else 0
        print(f"  Corr(H, var_y2) = {corr:.4f}")
        if corr < -0.5:
            print(f"  => H-maximizer tends to have LOW spectral variance (FLAT)")
        elif corr > 0.5:
            print(f"  => H-maximizer tends to have HIGH spectral variance (CONCENTRATED)")

        # Who is the maximizer?
        max_d = reps[0]
        if max_d['is_paley']:
            print(f"  MAXIMIZER: Paley (flat spectrum, dist_flat={max_d['dist_flat']:.6f})")
        elif max_d['is_interval']:
            print(f"  MAXIMIZER: Cyclic interval (concentrated, dist_flat={max_d['dist_flat']:.4f})")
        else:
            print(f"  MAXIMIZER: S={max_d['S']}, dist_flat={max_d['dist_flat']:.4f}")

    # ================================================================
    # MONOTONICITY: Is H monotonically decreasing in dist_flat
    #               for p=3 mod 4 (Paley primes)?
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"MONOTONICITY: H vs dist_flat")
    print(f"{'=' * 60}")

    for p in [7, 11]:
        m = (p - 1) // 2
        all_S = all_circulant_tournaments(p)

        data = []
        for S in all_S:
            H = ham_count_dp(p, S)
            pi = mpmath.pi
            y2 = [float(sum(mpmath.sin(2 * pi * k * s / p) for s in S) ** 2)
                  for k in range(1, m + 1)]
            dist_flat = sum((y - p/4)**2 for y in y2)
            data.append({'S': S, 'H': H, 'dist_flat': dist_flat})

        # Sort by dist_flat
        data.sort(key=lambda d: d['dist_flat'])

        print(f"\n  p={p}: Sorted by dist_flat:")
        print(f"  {'dist_flat':>12} {'H':>10}")
        prev_H = None
        monotone = True
        for d in data:
            marker = ""
            if prev_H is not None and d['H'] > prev_H:
                monotone = False
                marker = " <-- VIOLATION"
            prev_H = d['H']
            h = d['H']
            df = d['dist_flat']
            # Only print unique (H, dist_flat) pairs
            print(f"  {df:>12.6f} {h:>10}{marker}")

        print(f"  Monotone decreasing: {monotone}")


if __name__ == '__main__':
    main()
