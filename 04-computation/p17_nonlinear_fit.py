"""
p17_nonlinear_fit.py — Find the correct polynomial form for H at p=17.

PROVEN: H is NOT affine in e_k(y^2) at p=17 (with 100-digit mpmath).
This script finds what polynomial expression H IS.

Key question: what DEGREE and which MONOMIALS in e_k are needed?

With 16 data points, we can fit up to 15 free parameters.
Strategy:
  1. Start with power sums S_k = sum y_i^{2k} (equivalent to e_k via Newton)
  2. Try degree-2 monomials: S_j * S_k
  3. Stepwise selection: add the most informative term at each step

Also: investigate whether H depends on NON-SYMMETRIC functions of y^2,
which would indicate dependence on the eigenvalue ORDERING (cyclic structure).

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import math
import time
from collections import defaultdict
from itertools import combinations

import mpmath
mpmath.mp.dps = 50  # 50 digits is plenty (verified exact at p≤13 with 100)

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


def compute_data(p):
    """Compute H, y^2, and derived invariants for orbit representatives."""
    m = (p - 1) // 2
    reps = find_orbit_representatives(p)
    data = []
    for S in reps:
        H = ham_count_dp(p, S)

        # High-precision y_k^2
        pi = mpmath.pi
        y2 = []
        for k in range(1, m + 1):
            val = sum(mpmath.sin(2 * pi * k * s / p) for s in S)
            y2.append(val * val)

        # Elementary symmetric polynomials
        esyms = {}
        esyms[0] = mpmath.mpf(1)
        for kk in range(1, m + 1):
            total = mpmath.mpf(0)
            for subset in combinations(range(m), kk):
                prod = mpmath.mpf(1)
                for i in subset:
                    prod *= y2[i]
                total += prod
            esyms[kk] = total

        # Power sums
        psums = {}
        for kk in range(1, m + 1):
            psums[kk] = sum(y ** kk for y in y2)

        data.append({
            'S': S, 'H': H, 'y2': y2, 'esyms': esyms, 'psums': psums
        })

    data.sort(key=lambda d: d['H'], reverse=True)
    return data


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


def try_fit(data, feature_funcs, feature_names, label=""):
    """Try fitting H to given features. Returns (R^2, max_res, coeffs)."""
    n = len(data)
    H_vals = sorted(set(d['H'] for d in data), reverse=True)
    n_H = len(H_vals)
    reps = [next(d for d in data if d['H'] == hv) for hv in H_vals]

    n_feats = len(feature_funcs)
    n_vars = 1 + n_feats

    if n_H < n_vars:
        return None, None, None

    X = []
    y_vec = []
    for d in reps:
        row = [mpmath.mpf(1)] + [f(d) for f in feature_funcs]
        X.append(row)
        y_vec.append(mpmath.mpf(d['H']))

    XtX = [[sum(X[i][r] * X[i][c] for i in range(n_H))
            for c in range(n_vars)] for r in range(n_vars)]
    Xty = [sum(X[i][r] * y_vec[i] for i in range(n_H)) for r in range(n_vars)]

    coeffs = solve_mpmath(XtX, Xty)
    if coeffs is None:
        return None, None, None

    residuals = []
    for d in reps:
        row = [mpmath.mpf(1)] + [f(d) for f in feature_funcs]
        pred = sum(coeffs[j] * row[j] for j in range(n_vars))
        residuals.append(mpmath.mpf(d['H']) - pred)

    max_res = max(abs(r) for r in residuals)
    mean_H = sum(y_vec) / n_H
    ss_tot = sum((y - mean_H)**2 for y in y_vec)
    ss_res = sum(r**2 for r in residuals)
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else 1.0

    return r2, float(max_res), coeffs


def main():
    p = 17
    m = (p - 1) // 2
    print(f"{'=' * 70}")
    print(f"p = {p}, m = {m}: NONLINEAR POLYNOMIAL SEARCH FOR H")
    print(f"{'=' * 70}")

    t0 = time.time()
    data = compute_data(p)
    print(f"  Computed {len(data)} orbit reps in {time.time()-t0:.1f}s")

    H_vals = sorted(set(d['H'] for d in data), reverse=True)
    n_H = len(H_vals)
    print(f"  {n_H} distinct H values")

    # ================================================================
    # SECTION 1: Power sum basis (equivalent to e_k for degree 1)
    # ================================================================
    print(f"\n  === POWER SUMS S_k = sum y_i^(2k) ===")
    for n_terms in range(1, min(m, n_H)):
        # Skip S_1 since it's constant (Parseval); use S_2, ..., S_{n_terms+1}
        feats = [lambda d, k=k: d['psums'][k] for k in range(2, 2 + n_terms)]
        names = [f"S_{k}" for k in range(2, 2 + n_terms)]
        r2, max_res, _ = try_fit(data, feats, names)
        if r2 is not None:
            print(f"  S_2..S_{1+n_terms}: R²={r2:.8f}, max|res|={max_res:.2f}")

    # ================================================================
    # SECTION 2: Products of y^2 values — NON-SYMMETRIC basis
    # ================================================================
    print(f"\n  === INDIVIDUAL y_k^2 (NON-SYMMETRIC) ===")
    # Check if H depends on the ORDERING of y^2 (not just symmetric functions)
    # Use y_1^2, ..., y_m^2 directly (m=8 features, 9 variables total)
    feats_ind = [lambda d, k=k: d['y2'][k] for k in range(m)]
    names_ind = [f"y_{k+1}^2" for k in range(m)]
    r2_ind, max_res_ind, coeffs_ind = try_fit(data, feats_ind, names_ind)
    if r2_ind is not None:
        print(f"  y_1^2..y_8^2 (affine): R²={r2_ind:.8f}, max|res|={max_res_ind:.2f}")
    else:
        print(f"  y_1^2..y_8^2 (affine): SINGULAR")

    # Compare with symmetric (e_k)
    feats_ek = [lambda d, k=k: d['esyms'][k] for k in range(2, m + 1)]
    names_ek = [f"e_{k}" for k in range(2, m + 1)]
    r2_ek, max_res_ek, _ = try_fit(data, feats_ek, names_ek)
    print(f"  e_2..e_8 (affine): R²={r2_ek:.8f}, max|res|={max_res_ek:.2f}")

    if r2_ind is None or r2_ek is None:
        print(f"  Cannot compare (one fit failed)")
    elif abs(r2_ind - r2_ek) > 0.01:
        print(f"  *** NON-SYMMETRIC terms help! H depends on eigenvalue ORDERING ***")
    else:
        print(f"  Same R² — H is symmetric in y^2 (as expected)")

    # ================================================================
    # SECTION 3: Degree-2 monomials — which products matter?
    # ================================================================
    print(f"\n  === DEGREE-2 MONOMIAL SEARCH ===")

    # Start with best affine fit (e_2..e_8), then add product terms
    base_feats = [lambda d, k=k: d['esyms'][k] for k in range(2, m + 1)]
    base_names = [f"e_{k}" for k in range(2, m + 1)]

    # Try adding each possible product e_j * e_k
    best_r2 = r2_ek
    best_prod = None
    for j in range(2, m + 1):
        for k in range(j, m + 1):
            feat = lambda d, j=j, k=k: d['esyms'][j] * d['esyms'][k]
            test_feats = base_feats + [feat]
            r2, max_res, _ = try_fit(data, test_feats, base_names + [f"e_{j}*e_{k}"])
            if r2 is not None and r2 > best_r2 + 0.001:
                if best_prod is None or r2 > best_r2:
                    best_r2 = r2
                    best_prod = (j, k)
                    print(f"  +e_{j}*e_{k}: R²={r2:.8f}, max|res|={max_res:.2f}")

    if best_prod:
        print(f"  Best product term: e_{best_prod[0]}*e_{best_prod[1]}")
    else:
        print(f"  No single product term improves R² significantly")

    # ================================================================
    # SECTION 4: Greedy stepwise selection of degree-2 terms
    # ================================================================
    print(f"\n  === GREEDY STEPWISE (adding best term at each step) ===")

    # Available terms: all products e_j * e_k for j,k >= 2, j <= k
    all_products = [(j, k) for j in range(2, m + 1) for k in range(j, m + 1)]

    current_feats = list(base_feats)
    current_names = list(base_names)
    used_products = set()

    for step in range(min(8, n_H - len(current_feats) - 2)):  # leave room
        best_r2_step = -1
        best_term = None
        best_feat = None
        for j, k in all_products:
            if (j, k) in used_products:
                continue
            feat = lambda d, j=j, k=k: d['esyms'][j] * d['esyms'][k]
            test_feats = current_feats + [feat]
            r2, max_res, _ = try_fit(data, test_feats,
                                      current_names + [f"e_{j}*e_{k}"])
            if r2 is not None and r2 > best_r2_step:
                best_r2_step = r2
                best_term = (j, k)
                best_max_res = max_res

        if best_term is None:
            break

        j, k = best_term
        feat = lambda d, j=j, k=k: d['esyms'][j] * d['esyms'][k]
        current_feats.append(feat)
        current_names.append(f"e_{j}*e_{k}")
        used_products.add(best_term)

        print(f"  Step {step+1}: add e_{j}*e_{k}, "
              f"R²={best_r2_step:.10f}, max|res|={best_max_res:.2f}")

        if best_max_res < 0.01:
            print(f"  *** EXACT FIT with {len(current_names)} total terms! ***")
            # Get final coefficients
            _, _, final_coeffs = try_fit(data, current_feats, current_names)
            print(f"\n  Formula:")
            for i, name in enumerate(current_names):
                c = float(final_coeffs[i + 1])
                print(f"    {name}: c = {c:.6f}")
            print(f"    const: {float(final_coeffs[0]):.6f}")
            break

    # ================================================================
    # SECTION 5: Trace-based invariants
    # ================================================================
    print(f"\n  === TRACE-BASED INVARIANTS ===")
    # tr(A^k) = sum lambda_j^k. For circulant, lambda_k = sum_{s in S} omega^{ks}
    # The trace is sum_{k=0}^{n-1} lambda_k^power
    for d in data:
        S = d['S']
        traces = {}
        for power in range(1, 12):
            tr = mpmath.mpf(0)
            pi = mpmath.pi
            for k in range(p):
                lam = sum(mpmath.exp(2j * pi * k * s / p) for s in S)
                tr += lam ** power
            traces[power] = tr
        d['traces'] = traces

    # H vs traces
    for n_terms in [2, 3, 4, 5, 6, 7]:
        # Use traces of A^3, A^4, ..., A^{n_terms+2}
        feats = [lambda d, k=k: d['traces'][k].real for k in range(3, 3 + n_terms)]
        names = [f"tr(A^{k})" for k in range(3, 3 + n_terms)]
        r2, max_res, _ = try_fit(data, feats, names)
        if r2 is not None:
            print(f"  tr(A^3)..tr(A^{2+n_terms}): R²={r2:.8f}, max|res|={max_res:.2f}")

    # ================================================================
    # SECTION 6: Cycle counts directly
    # ================================================================
    print(f"\n  === CYCLE COUNTS c_k ===")
    # c_k = number of directed k-cycles in tournament
    # For circulant: c_k = tr(A^k)/k for k <= threshold (no non-simple walk corrections)
    # At p=17, c_3 = tr(A^3)/3 is constant for all regular circulants
    # c_5 = tr(A^5)/5 (exact for tournaments)
    for d in data:
        d['c3'] = d['traces'][3].real / 3
        d['c5'] = d['traces'][5].real / 5
        d['c7'] = d['traces'][7].real / 7  # might need correction

    # Check c3 constancy
    c3_vals = set(round(float(d['c3']), 2) for d in data)
    print(f"  c_3 values: {c3_vals} (should be single value)")

    c5_vals = sorted(set(round(float(d['c5']), 2) for d in data))
    print(f"  c_5 values ({len(c5_vals)}): {c5_vals[:5]}...")

    # Fit H to c_5
    feats_c5 = [lambda d: d['c5'].real]
    r2_c5, max_res_c5, _ = try_fit(data, feats_c5, ['c5'])
    print(f"  H vs c_5: R²={r2_c5:.8f}")

    # H vs (c_5, c_7)
    feats_c57 = [lambda d: d['c5'].real, lambda d: d['c7'].real]
    r2_c57, max_res_c57, _ = try_fit(data, feats_c57, ['c5', 'c7'])
    print(f"  H vs (c_5, c_7): R²={r2_c57:.8f}")

    # H vs (c_5, c_7, c_5^2)
    feats_c572 = [lambda d: d['c5'].real, lambda d: d['c7'].real,
                  lambda d: d['c5'].real ** 2]
    r2_c572, _, _ = try_fit(data, feats_c572, ['c5', 'c7', 'c5²'])
    print(f"  H vs (c_5, c_7, c_5²): R²={r2_c572:.8f}")

    # ================================================================
    # SECTION 7: Compare with p=13 to see if degree-2 IS needed there
    # ================================================================
    print(f"\n  === CROSS-CHECK: p=13 ===")
    data13 = compute_data(13)
    m13 = 6

    feats13 = [lambda d, k=k: d['esyms'][k] for k in range(2, m13 + 1)]
    r2_13, max_res_13, _ = try_fit(data13, feats13, [f"e_{k}" for k in range(2, m13 + 1)])
    print(f"  p=13 affine (e_2..e_6): R²={r2_13:.10f}, max|res|={max_res_13}")

    # Does p=13 need degree 2? Add ALL products
    feats13_ext = list(feats13)
    names13_ext = [f"e_{k}" for k in range(2, m13 + 1)]
    for j in range(2, m13 + 1):
        for k in range(j, m13 + 1):
            feats13_ext.append(lambda d, j=j, k=k: d['esyms'][j] * d['esyms'][k])
            names13_ext.append(f"e_{j}*e_{k}")

    r2_13q, max_res_13q, _ = try_fit(data13, feats13_ext, names13_ext)
    if r2_13q is not None:
        print(f"  p=13 affine+quadratic: R²={r2_13q:.10f}, max|res|={max_res_13q}")


if __name__ == '__main__':
    main()
