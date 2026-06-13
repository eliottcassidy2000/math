"""
p17_exact_fit.py — Fix the p=17 e_k fitting using EXACT arithmetic.

The p17_test.py fitting gave R²=0.865 with 7 terms, which is unexpected
since p=7,11,13 all gave exact fits. Hypothesis: floating-point precision
loss in elementary_symmetric() and the linear solver.

Strategy:
  1. Compute eigenvalues symbolically: lambda_k = sum_{s in S} omega^{ks}
     For circulant tournaments, y_k^2 = Im(lambda_k)^2 are algebraic numbers.
  2. Use Fraction arithmetic for the linear system.
  3. Also test: does the number of distinct (e_2,...,e_m) tuples = n_H?
     If not, the e_k basis doesn't separate H-classes and we need more variables.

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import cmath
import math
import time
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
from decimal import Decimal, getcontext

# High precision
getcontext().prec = 50

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
    """Held-Karp DP for Hamiltonian path count."""
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
        popcount = bin(mask).count('1')
        if popcount >= n:
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
        new_S = frozenset((a * s) % p for s in S)
        orbit.add(new_S)
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


def y_squared_exact(p, S):
    """Compute y_k^2 = Im(lambda_k)^2 using high-precision trig.

    lambda_k = sum_{s in S} omega^{ks}, omega = exp(2pi i/p)
    Im(lambda_k) = sum_{s in S} sin(2*pi*k*s/p)
    y_k^2 = (sum sin(2*pi*k*s/p))^2

    For EXACT computation: use the identity that for circulant tournaments,
    y_k^2 values are algebraic. But let's first check if the issue is
    just the elementary_symmetric computation precision.
    """
    m = (p - 1) // 2
    y2 = []
    for k in range(1, m + 1):
        val = sum(math.sin(2 * math.pi * k * s / p) for s in S)
        y2.append(val * val)
    return y2


def elementary_symmetric_decimal(xs, k):
    """Compute e_k using Decimal for higher precision."""
    if k == 0: return Decimal(1)
    if k > len(xs): return Decimal(0)
    xs_d = [Decimal(str(x)) for x in xs]
    total = Decimal(0)
    for subset in combinations(range(len(xs_d)), k):
        prod = Decimal(1)
        for i in subset:
            prod *= xs_d[i]
        total += prod
    return total


def solve_exact(A, b):
    """Solve Ax=b using Fraction arithmetic (exact)."""
    n = len(b)
    M = [[Fraction(A[r][c]) for c in range(len(A[r]))] + [Fraction(b[r])] for r in range(n)]

    for col in range(n):
        # Find pivot
        max_row = col
        for r in range(col + 1, n):
            if abs(M[r][col]) > abs(M[max_row][col]):
                max_row = r
        if M[max_row][col] == 0:
            return None
        M[col], M[max_row] = M[max_row], M[col]

        for row in range(col + 1, n):
            f = M[row][col] / M[col][col]
            for j in range(col, n + 1):
                M[row][j] -= f * M[col][j]

    x = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]
    return x


def main():
    p = 17
    m = (p - 1) // 2  # = 8
    print(f"{'=' * 70}")
    print(f"p = {p}, m = {m}: EXACT ARITHMETIC e_k FIT")
    print(f"{'=' * 70}")

    # Step 1: Get orbit representatives and H values
    reps = find_orbit_representatives(p)
    print(f"\n  {len(reps)} orbit representatives")

    data = []
    for i, S in enumerate(reps):
        t0 = time.time()
        H = ham_count_dp(p, S)
        elapsed = time.time() - t0
        y2 = y_squared_exact(p, S)

        # Compute e_k with high precision
        esyms_float = [float(elementary_symmetric_decimal(y2, k)) for k in range(m + 1)]

        orbit_size = len(multiplicative_orbit(p, S))
        data.append({'S': S, 'H': H, 'y2': y2, 'esyms': esyms_float, 'orbit': orbit_size})
        print(f"  [{i+1}/{len(reps)}] H={H}, orbit={orbit_size}, time={elapsed:.1f}s")

    # Step 2: Check how many DISTINCT e_k tuples there are
    data.sort(key=lambda d: d['H'], reverse=True)
    H_vals = sorted(set(d['H'] for d in data), reverse=True)
    n_H = len(H_vals)
    print(f"\n  {n_H} distinct H values")

    # Check if e_2 alone distinguishes H classes
    e2_to_H = defaultdict(set)
    for d in data:
        e2_to_H[round(d['esyms'][2], 6)].add(d['H'])

    n_e2 = len(e2_to_H)
    collisions_e2 = sum(1 for v in e2_to_H.values() if len(v) > 1)
    print(f"  Distinct e_2 values: {n_e2}, collisions: {collisions_e2}")

    # Check if (e_2, e_3) distinguishes
    e23_to_H = defaultdict(set)
    for d in data:
        key = (round(d['esyms'][2], 4), round(d['esyms'][3], 4))
        e23_to_H[key].add(d['H'])

    n_e23 = len(e23_to_H)
    collisions_e23 = sum(1 for v in e23_to_H.values() if len(v) > 1)
    print(f"  Distinct (e_2, e_3) pairs: {n_e23}, collisions: {collisions_e23}")

    # Check full tuple
    ek_to_H = defaultdict(set)
    for d in data:
        key = tuple(round(d['esyms'][k], 4) for k in range(2, m + 1))
        ek_to_H[key].add(d['H'])

    n_ek = len(ek_to_H)
    collisions_ek = sum(1 for v in ek_to_H.values() if len(v) > 1)
    print(f"  Distinct (e_2,...,e_{m}) tuples: {n_ek}, collisions: {collisions_ek}")

    # Print detailed table
    print(f"\n  {'H':>15} {'e_2':>12} {'e_3':>14} {'e_4':>14} {'e_5':>14}")
    reps_per_H = []
    for H_val in H_vals:
        d = next(x for x in data if x['H'] == H_val)
        reps_per_H.append(d)
        print(f"  {H_val:>15} {d['esyms'][2]:>12.4f} {d['esyms'][3]:>14.4f} "
              f"{d['esyms'][4]:>14.4f} {d['esyms'][5]:>14.4f}")

    # Step 3: Check if different H values share the same e_k tuple
    print(f"\n  === COLLISION ANALYSIS ===")
    for key, h_set in sorted(ek_to_H.items()):
        if len(h_set) > 1:
            print(f"  COLLISION: e_k tuple maps to H values: {sorted(h_set)}")
            # Find the actual data points
            for d in data:
                key2 = tuple(round(d['esyms'][k], 4) for k in range(2, m + 1))
                if key2 == key:
                    print(f"    S={list(d['S'])}, H={d['H']}")
                    print(f"    y^2 = {[round(y, 6) for y in d['y2']]}")

    if collisions_ek == 0:
        print(f"  No collisions — e_k tuple uniquely determines H")
        print(f"  The fitting failure must be NUMERICAL PRECISION")
    else:
        print(f"  {collisions_ek} collisions — e_k basis INSUFFICIENT at p={p}!")
        print(f"  Need additional polynomial terms (power sums, mixed terms)")

    # Step 4: Exact fit using Fraction arithmetic
    print(f"\n  === EXACT FIT (Fraction arithmetic) ===")

    # Convert e_k values to high-precision Decimals for the linear system
    # Use ALL m-1 terms: e_2, ..., e_m
    n_terms = min(n_H - 1, m - 1)  # max terms we can fit
    k_start = 2
    k_end = k_start + n_terms
    n_vars = 1 + n_terms

    print(f"  Using {n_terms} terms (e_{k_start}..e_{k_end-1}), {n_vars} variables")
    print(f"  {n_H} data points ({'exact' if n_H == n_vars else 'overdetermined'})")

    # If overdetermined, use least squares with exact arithmetic
    # Build X matrix and y vector with high precision
    X = []
    y_vec = []
    for d in reps_per_H:
        row = [1.0] + [d['esyms'][k] for k in range(k_start, k_end)]
        X.append(row)
        y_vec.append(float(d['H']))

    # Normal equations with exact Fraction arithmetic
    XtX = [[sum(Fraction(X[i][r]).limit_denominator(10**15) *
                Fraction(X[i][c]).limit_denominator(10**15)
                for i in range(n_H)) for c in range(n_vars)] for r in range(n_vars)]
    Xty = [sum(Fraction(X[i][r]).limit_denominator(10**15) *
               Fraction(y_vec[i]).limit_denominator(10**15)
               for i in range(n_H)) for r in range(n_vars)]

    coeffs = solve_exact(XtX, Xty)
    if coeffs is None:
        print(f"  SINGULAR system!")
    else:
        # Compute residuals
        residuals = []
        for d in reps_per_H:
            pred = float(coeffs[0])
            for j, k in enumerate(range(k_start, k_end)):
                pred += float(coeffs[j + 1]) * d['esyms'][k]
            residuals.append(d['H'] - pred)

        max_res = max(abs(r) for r in residuals)
        ss_res = sum(r**2 for r in residuals)
        mean_H = sum(d['H'] for d in reps_per_H) / n_H
        ss_tot = sum((d['H'] - mean_H)**2 for d in reps_per_H)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 1

        print(f"  R² = {r2:.12f}")
        print(f"  max |residual| = {max_res:.6f}")

        print(f"\n  Coefficients:")
        for j in range(n_vars):
            k = 0 if j == 0 else k_start + j - 1
            print(f"    c_{k} = {float(coeffs[j]):.6f} = {coeffs[j]}")

    # Step 5: Try power sum basis instead
    print(f"\n  === POWER SUM BASIS ===")
    # p_k = sum y_i^{2k} (power sums of y^2)
    for d in data:
        d['psums'] = [sum(y**(2*k) for y in [yi**0.5 for yi in d['y2']]) if k > 0
                      else float(m) for k in range(m + 1)]
        # Actually: p_k(y^2) = sum y_i^{2k}
        d['psums'] = [sum(yi**k for yi in d['y2']) for k in range(m + 1)]

    # Check separation
    ps_to_H = defaultdict(set)
    for d in data:
        key = tuple(round(d['psums'][k], 4) for k in range(2, m + 1))
        ps_to_H[key].add(d['H'])

    n_ps = len(ps_to_H)
    collisions_ps = sum(1 for v in ps_to_H.values() if len(v) > 1)
    print(f"  Distinct power sum tuples: {n_ps}, collisions: {collisions_ps}")

    # Step 6: Check individual y_k^2 values — are there degeneracies?
    print(f"\n  === y^2 VALUE ANALYSIS ===")
    # How many distinct SORTED y^2 tuples are there?
    sorted_y2_to_H = defaultdict(set)
    for d in data:
        key = tuple(sorted([round(y, 6) for y in d['y2']]))
        sorted_y2_to_H[key].add(d['H'])

    n_sorted = len(sorted_y2_to_H)
    collisions_sorted = sum(1 for v in sorted_y2_to_H.values() if len(v) > 1)
    print(f"  Distinct sorted y^2 tuples: {n_sorted}, collisions: {collisions_sorted}")

    if collisions_sorted > 0:
        print(f"  *** CRITICAL: Same y^2 multiset gives DIFFERENT H! ***")
        print(f"  This means H is NOT a symmetric function of y^2 at p=17!")
        for key, h_set in sorted(sorted_y2_to_H.items()):
            if len(h_set) > 1:
                print(f"\n  Collision: sorted y^2 = {key}")
                print(f"  H values: {sorted(h_set)}")
                for d in data:
                    key2 = tuple(sorted([round(y, 6) for y in d['y2']]))
                    if key2 == key:
                        print(f"    S={list(d['S'])}: H={d['H']}")
                        print(f"    y^2 = {[round(y, 6) for y in d['y2']]}")
    else:
        print(f"  No collisions — y^2 multiset uniquely determines H")
        print(f"  H IS a symmetric function of y^2 (as expected)")
        print(f"  The fitting issue is purely numerical")

    # Step 7: Use UNSORTED y^2 (order matters via conjugate pairing)
    print(f"\n  === ORDERED y^2 ANALYSIS ===")
    # y_k^2 for k=1,...,m. Note y_k = y_{p-k} for circulant, so only m independent
    for d in data[:3]:
        print(f"  S={list(d['S'])}:")
        for k in range(1, m + 1):
            print(f"    y_{k}^2 = {d['y2'][k-1]:.8f}")


if __name__ == '__main__':
    main()
