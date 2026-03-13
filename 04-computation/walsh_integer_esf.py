#!/usr/bin/env python3
"""
walsh_integer_esf.py -- Integer elementary symmetric functions of Q_k

DISCOVERY: The elementary symmetric functions e_j(Q_1,...,Q_m) of the
Fourier power spectrum Q_k = |S_hat(k)|^2 are ALL INTEGERS for every
circulant tournament orientation.

Moreover:
  - e_1 = m(m+1)/2 = sum Q_k is CONSTANT (Parseval)
  - For Interval: prod Q_k = e_m = 1 EXACTLY
  - For Paley: prod Q_k = ((p+1)/4)^m (Q_k constant = (p+1)/4)

This script:
1. Verifies integrality of e_j at all tested primes
2. Finds the minimal polynomial of Q_k values
3. Determines the exact algebraic structure
4. Tests if H is a polynomial in e_2,...,e_m (since e_1 is constant)
5. Explores the connection between e_m (= prod Q_k) and the tournament

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
from fractions import Fraction
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
    print("INTEGER ELEMENTARY SYMMETRIC FUNCTIONS OF Q_k")
    print("=" * 70)

    for p in [5, 7, 11, 13]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"{'='*70}")

        # Compute H and Q_k for all orientations
        pairs = [(s, p - s) for s in range(1, m + 1)]
        H_dict = {}
        Q_dict = {}  # sigma -> Q_k values (k=1..m)

        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = []
            for i, (a, b) in enumerate(pairs):
                S.append(a if sigma[i] == 1 else b)
            S = sorted(S)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            H_dict[sigma] = H

            # Compute Q_k = |S_hat(k)|^2 for k=1,...,m
            Q_vals = []
            for k in range(1, m + 1):
                val = sum(omega ** (k * s) for s in S)
                Q_vals.append(abs(val)**2)
            Q_dict[sigma] = Q_vals

        # ====== 1. ELEMENTARY SYMMETRIC FUNCTIONS ======
        print(f"\n  1. ELEMENTARY SYMMETRIC FUNCTIONS")

        all_esf = {}  # sigma -> list of e_j values
        for sigma in H_dict:
            Q = Q_dict[sigma]
            e_vals = [1.0]  # e_0 = 1
            for j in range(1, m + 1):
                ej = sum(math.prod(Q[i] for i in combo)
                        for combo in combinations(range(m), j))
                e_vals.append(ej)
            all_esf[sigma] = e_vals

        # Check integrality
        print(f"\n    Integrality check:")
        max_frac = 0
        for sigma in H_dict:
            e_vals = all_esf[sigma]
            for j in range(1, m + 1):
                frac = abs(e_vals[j] - round(e_vals[j]))
                max_frac = max(max_frac, frac)

        print(f"    Max fractional part of any e_j: {max_frac:.2e}")
        print(f"    All e_j are integers: {'YES' if max_frac < 1e-6 else 'NO'}")

        # ====== 2. TABLE OF (e_j, H) VALUES ======
        print(f"\n  2. TABLE OF (e_j) -> H")
        # Group by e_j profile
        esf_groups = defaultdict(list)
        for sigma in H_dict:
            H = H_dict[sigma]
            e_vals = tuple(round(all_esf[sigma][j]) for j in range(1, m + 1))
            esf_groups[e_vals].append((sigma, H))

        print(f"\n    {'e_1':>6} {'e_2':>8} {'e_3':>10} {'e_4':>10} {'e_5':>10} "
              f"{'e_m':>10} {'H':>12} {'count':>6}")
        print("    " + "-" * 80)
        for e_vals, items in sorted(esf_groups.items(), key=lambda x: -x[1][0][1]):
            H = items[0][1]
            e_str = " ".join(f"{e_vals[j]:>10}" for j in range(min(m, 5)))
            if m > 5:
                e_str += " ..."
            e_m_val = e_vals[m - 1]
            print(f"    {e_str} {e_m_val:>10} {H:>12} {len(items):>6}")

        # ====== 3. PROD Q_k ANALYSIS ======
        print(f"\n  3. PRODUCT OF Q_k (= e_m)")
        # For Interval: compute analytically
        # |S_hat(k)|^2 for S={1,...,m} is:
        # |sum_{j=1}^m omega^{jk}|^2 = |(omega^{(m+1)k} - omega^k)/(omega^k - 1)|^2
        # = sin^2(pi*k*(m+1)/p) / sin^2(pi*k/p)
        # Since m+1 = (p+1)/2:
        # = sin^2(pi*k*(p+1)/(2p)) / sin^2(pi*k/p)

        print(f"\n    Interval |S_hat(k)|^2 formula:")
        for k in range(1, m + 1):
            num = math.sin(math.pi * k * (m + 1) / p) ** 2
            den = math.sin(math.pi * k / p) ** 2
            Q_formula = num / den
            sigma_int = tuple([1] * m)
            Q_actual = Q_dict[sigma_int][k - 1]
            print(f"      k={k}: sin^2({k*(m+1)}/{p}*pi)/sin^2({k}/{p}*pi) "
                  f"= {Q_formula:.8f} [actual: {Q_actual:.8f}]")

        # Product of sin^2(k*(m+1)/p * pi) / sin^2(k/p * pi) for k=1..m
        prod_num = math.prod(math.sin(math.pi * k * (m + 1) / p) ** 2
                            for k in range(1, m + 1))
        prod_den = math.prod(math.sin(math.pi * k / p) ** 2
                            for k in range(1, m + 1))
        prod_Q_interval = prod_num / prod_den

        print(f"\n    prod Q_k (Interval) = {prod_Q_interval:.10f}")
        print(f"    Should be 1: {'YES' if abs(prod_Q_interval - 1) < 1e-8 else 'NO'}")

        # Known identity: prod_{k=1}^{n-1} sin(k*pi/n) = n / 2^{n-1}
        # So prod_{k=1}^{m} sin(k*pi/p) is related to this.
        # prod_{k=1}^{p-1} sin(k*pi/p) = p / 2^{p-1}
        # By symmetry sin(k*pi/p) = sin((p-k)*pi/p), so
        # prod_{k=1}^{m} sin(k*pi/p)^2 = p / 2^{p-1}

        prod_sin_1 = math.prod(math.sin(math.pi * k / p) for k in range(1, m + 1))
        prod_sin_m1 = math.prod(math.sin(math.pi * k * (m + 1) / p)
                                for k in range(1, m + 1))
        print(f"\n    prod sin(k*pi/p), k=1..m: {prod_sin_1:.10f}")
        print(f"    sqrt(p) / 2^m: {math.sqrt(p) / (2**m):.10f}")
        print(f"    prod sin(k*(m+1)*pi/p), k=1..m: {prod_sin_m1:.10f}")
        print(f"    These are equal: {abs(prod_sin_1 - prod_sin_m1) < 1e-8}")
        # If prod_sin_m1 = prod_sin_1, then prod Q_k = 1. Let's verify.

        # Why? m+1 = (p+1)/2. The set {k*(m+1) mod p : k=1,...,m}
        # Since gcd(m+1, p) = gcd((p+1)/2, p) = 1 (p prime, p odd),
        # k*(m+1) runs over all nonzero residues mod p as k runs.
        # But we want k=1,...,m. Then k*(m+1) mod p runs over some subset.
        mult_set = sorted([(k * (m + 1)) % p for k in range(1, m + 1)])
        print(f"\n    {{k*(m+1) mod p : k=1,...,m}} = {mult_set}")
        print(f"    {{1,...,m}} = {list(range(1, m+1))}")
        is_same_set = (set(mult_set) == set(range(1, m + 1)) or
                       set(mult_set) == set(range(m + 1, p)))

        # Actually: k*(m+1) mod p for k=1,...,m
        # Since m+1 = (p+1)/2, k*(p+1)/2 mod p = k/2 mod p (if 2 has inverse)
        # 2^{-1} mod p = (p+1)/2. So k*(p+1)/2 = k * 2^{-1} mod p.
        # {k*2^{-1} mod p : k=1,...,m} = {2^{-1}, 2*2^{-1},..., m*2^{-1}} mod p
        # Since 2 is invertible, this is a permutation of some m elements of Z_p*.
        # But sin(x*pi/p) = sin((p-x)*pi/p), so we need to check if
        # {k*(m+1) mod p} gives the same multiset of sin values as {k}.

        # The set {k*(m+1) mod p : k=1,...,m} as canonical:
        canonical_mult = sorted(min(x, p - x) for x in mult_set)
        canonical_orig = list(range(1, m + 1))
        print(f"    Canonical mult set: {canonical_mult}")
        print(f"    Canonical orig set: {canonical_orig}")
        print(f"    Same multiset: {canonical_mult == canonical_orig}")

        if canonical_mult == canonical_orig:
            print(f"    => prod sin(k*(m+1)*pi/p) = prod sin(k*pi/p)")
            print(f"    => prod Q_k (Interval) = 1 [PROVED]")

        # ====== 4. PALEY PRODUCT ======
        print(f"\n    Paley product: ((p+1)/4)^m = {((p+1)/4)**m:.0f}")
        sigma_paley = tuple(1 if (i+1) in QR else -1 for i in range(m))
        if sigma_paley in Q_dict:
            prod_Q_paley = math.prod(Q_dict[sigma_paley])
            print(f"    Actual prod Q_k (Paley) = {prod_Q_paley:.6f}")

        # ====== 5. ALL DISTINCT e_m VALUES ======
        print(f"\n  4. DISTINCT e_m = prod Q_k VALUES")
        e_m_vals = set()
        for sigma in H_dict:
            e_m = round(all_esf[sigma][m])
            e_m_vals.add(e_m)

        for val in sorted(e_m_vals):
            count = sum(1 for sigma in H_dict if round(all_esf[sigma][m]) == val)
            # Factor
            factors = []
            rem = abs(val)
            for f in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
                while rem > 0 and rem % f == 0:
                    factors.append(f)
                    rem //= f
            if rem > 1:
                factors.append(rem)
            fac_str = " * ".join(str(f) for f in factors) if factors else "1"
            print(f"    e_m = {val} = {fac_str}, count = {count}")

        # ====== 6. POLYNOMIAL FIT H = f(e_2, e_3, ..., e_m) ======
        print(f"\n  5. POLYNOMIAL FIT H = polynomial(e_2, ..., e_m)")
        # e_1 is constant, so drop it.

        # Get unique (e_2,...,e_m) profiles
        unique_profiles = {}
        for sigma in H_dict:
            H = H_dict[sigma]
            e_profile = tuple(round(all_esf[sigma][j]) for j in range(2, m + 1))
            if e_profile not in unique_profiles:
                unique_profiles[e_profile] = H

        n_unique = len(unique_profiles)
        n_vars = m - 1  # e_2, ..., e_m
        print(f"    {n_unique} unique profiles, {n_vars} variables")

        # Try linear: H = a_0 + sum a_j * e_{j+1}
        profiles_list = list(unique_profiles.items())
        if n_unique >= n_vars + 1:
            # Build matrix
            A_mat = []
            b_vec = []
            for e_profile, H in profiles_list:
                A_mat.append([1.0] + list(e_profile))
                b_vec.append(float(H))

            n_eq = len(A_mat)
            n_col = n_vars + 1

            # Solve via least squares (normal equations)
            XTX = [[0.0]*n_col for _ in range(n_col)]
            XTb = [0.0]*n_col
            for r in range(n_eq):
                for i in range(n_col):
                    for j in range(n_col):
                        XTX[i][j] += A_mat[r][i] * A_mat[r][j]
                    XTb[i] += A_mat[r][i] * b_vec[r]

            # Gaussian elimination
            aug = [XTX[i][:] + [XTb[i]] for i in range(n_col)]
            for col in range(min(n_col, n_eq)):
                max_row = max(range(col, n_col), key=lambda r: abs(aug[r][col]))
                aug[col], aug[max_row] = aug[max_row], aug[col]
                pivot = aug[col][col]
                if abs(pivot) < 1e-20:
                    print(f"    Singular at col {col}")
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

            residuals = []
            for r in range(n_eq):
                H_pred = sum(coeffs[i] * A_mat[r][i] for i in range(n_col))
                residuals.append(b_vec[r] - H_pred)

            max_resid = max(abs(r) for r in residuals)
            print(f"    Linear model: max residual = {max_resid:.6f}")
            print(f"    Perfect fit: {'YES' if max_resid < 0.01 else 'NO'}")

            if max_resid < 0.01:
                print(f"    H = {coeffs[0]:.4f}", end="")
                for j in range(n_vars):
                    if abs(coeffs[j+1]) > 1e-10:
                        print(f" + ({coeffs[j+1]:.4f})*e_{j+2}", end="")
                print()

                # Try to identify rational coefficients
                print(f"\n    Rational reconstruction:")
                for j in range(n_col):
                    # Try fractions with small denominators
                    found = False
                    for denom in range(1, 200):
                        numer = round(coeffs[j] * denom)
                        if abs(coeffs[j] - numer / denom) < 1e-6:
                            print(f"      c_{j} = {numer}/{denom}")
                            found = True
                            break
                    if not found:
                        print(f"      c_{j} = {coeffs[j]:.10f} (not recognized)")
            else:
                # Try quadratic
                if n_unique >= n_col * (n_col + 1) // 2:
                    print(f"    Trying quadratic model...")
                    quad_features = []
                    for e_profile, H in profiles_list:
                        row = [1.0] + list(e_profile)
                        for i in range(n_vars):
                            for j in range(i, n_vars):
                                row.append(e_profile[i] * e_profile[j])
                        quad_features.append(row)

                    n_qcol = len(quad_features[0])
                    print(f"    {n_qcol} features")

                    if n_unique >= n_qcol:
                        XTX = [[0.0]*n_qcol for _ in range(n_qcol)]
                        XTb = [0.0]*n_qcol
                        for r in range(n_eq):
                            for i in range(n_qcol):
                                for j in range(n_qcol):
                                    XTX[i][j] += quad_features[r][i] * quad_features[r][j]
                                XTb[i] += quad_features[r][i] * b_vec[r]

                        aug = [XTX[i][:] + [XTb[i]] for i in range(n_qcol)]
                        for col in range(n_qcol):
                            max_row = max(range(col, n_qcol),
                                         key=lambda r: abs(aug[r][col]))
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
                            H_pred = sum(qcoeffs[i] * quad_features[r][i]
                                        for i in range(n_qcol))
                            qresid.append(b_vec[r] - H_pred)

                        max_qresid = max(abs(r) for r in qresid)
                        print(f"    Quadratic: max residual = {max_qresid:.6f}")
                        print(f"    Perfect fit: {'YES' if max_qresid < 0.01 else 'NO'}")

        # ====== 7. CHARACTERISTIC POLYNOMIAL OF Q_k ======
        print(f"\n  6. CHARACTERISTIC POLYNOMIAL OF Q_k")
        # The Q_k values for each orientation are roots of:
        # t^m - e_1 t^{m-1} + e_2 t^{m-2} - ... + (-1)^m e_m = 0
        # Since e_j are integers, Q_k are algebraic integers!

        for e_vals, items in sorted(esf_groups.items(), key=lambda x: -x[1][0][1]):
            H = items[0][1]
            poly_str = f"t^{m}"
            for j in range(1, m + 1):
                sign = "+" if (j % 2 == 0) else "-"
                poly_str += f" {sign} {abs(e_vals[j-1])}*t^{m-j}"
            Q_vals = Q_dict[items[0][0]]
            Q_sorted = sorted(Q_vals, reverse=True)
            q_str = ", ".join(f"{q:.4f}" for q in Q_sorted)
            print(f"    H={H}: {poly_str} = 0")
            print(f"           roots: [{q_str}]")


if __name__ == '__main__':
    main()
