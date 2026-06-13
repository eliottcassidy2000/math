#!/usr/bin/env python3
"""
p17_exact_arithmetic.py -- kind-pasteur-2026-03-13-S60

The floating-point cycle count computations at p=17 have errors ~24 due to
double precision limits (numbers ~10^14). This script uses:
1. Integer matrix multiplication for EXACT c_k values
2. mpmath high-precision for eigenvalue moments
3. Exact rational arithmetic for the linear fit

KEY TEST: With exact c_k values, does alpha_1 = linear(moments) hold as a
GENUINE identity when overconstrained (16 types, 8 params at p=17)?
"""

import cmath
from fractions import Fraction
from collections import defaultdict
import time


def build_adj_int(p, S):
    """Build adjacency matrix as list of lists of ints."""
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def mat_mul_int(A, B, n):
    """Integer matrix multiplication."""
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if A[i][k] == 0:
                continue
            for j in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C


def mat_pow_int(A, n, k):
    """Integer matrix power A^k."""
    # Identity matrix
    result = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    base = [row[:] for row in A]
    while k > 0:
        if k % 2 == 1:
            result = mat_mul_int(result, base, n)
        base = mat_mul_int(base, base, n)
        k //= 2
    return result


def trace_int(A, n):
    """Trace of integer matrix."""
    return sum(A[i][i] for i in range(n))


def compute_ck_exact(A, p, max_k=None):
    """Compute exact cycle counts using integer arithmetic."""
    if max_k is None:
        max_k = p
    counts = {}
    for k in range(3, max_k + 1, 2):
        Ak = mat_pow_int(A, p, k)
        tr = trace_int(Ak, p)
        counts[k] = tr // k
    return counts


def compute_D2_multiset(p, S):
    """Compute D^2 multiset using high-precision arithmetic."""
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D2_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D2_vals.append(lam.imag ** 2)
    return tuple(sorted(round(d, 8) for d in D2_vals))


def compute_moments_exact(p, S, max_moment=None):
    """Compute moments using high precision. Return as Fraction approximations."""
    m = (p - 1) // 2
    if max_moment is None:
        max_moment = 2 * m

    # Use the eigenvalue expansion: c_k = (1/k)[m^k + 2*sum Re(z_t^k)]
    # and Re(z^k) = sum_r C(k,2r)*(-1/2)^{k-2r}*(-1)^r * D_t^{2r}
    # So S_{2r} = sum_t D_t^{2r} enters as: c_k += (2/k)*C(k,2r)*(-1/2)^{k-2r}*(-1)^r * S_{2r}

    # For exact moments, we need exact D_t values. Use the formula:
    # D_t = Im(sum_{s in S} omega^{st}) = sum_{s in S} sin(2*pi*s*t/p)

    # Instead: compute S_{2r} from cycle counts!
    # From c_k = (1/k)[m^k + 2*sum_r coeff_r * S_{2r}], we can INVERT to get S_{2r}
    # from exact c_k values.

    # Actually, let's just compute the moments from the c_k inversion.
    # c_3: depends on S_2 only (constant) => S_2 = mp/4 (known)
    # c_5: depends on S_2, S_4 => S_4 = ...
    # c_7: depends on S_2, S_4, S_6 => S_6 = ...
    # etc.

    # For now, compute moments using standard complex arithmetic (double precision)
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)

    moments = {}
    for k in range(2, max_moment + 1, 2):
        moments[k] = sum(d**k for d in D_vals)
    return moments


def invert_moments_from_ck(p, ck):
    """Invert the eigenvalue expansion to get S_{2r} from exact c_k values.

    c_k = (1/k)[m^k + 2*sum_{r=0}^{floor(k/2)} C(k,2r)*(-1/2)^{k-2r}*(-1)^r * S_{2r}]
    where S_0 = m.

    Rearranging:
    k*c_k - m^k = 2*sum_{r=0}^{floor(k/2)} coeff_r * S_{2r}

    For r=0: coeff_0 = (-1/2)^k, S_0 = m => known contribution

    So: sum_{r=1}^{floor(k/2)} coeff_r * S_{2r} = (k*c_k - m^k - 2*(-1/2)^k*m) / 2

    Use this iteratively: c_3 gives S_2 (which we know is mp/4),
    c_5 gives S_4, c_7 gives S_6, etc.
    """
    from math import comb

    m = (p - 1) // 2
    S = {0: Fraction(m)}  # S_0 = m

    # S_2 = mp/4 (universal constant)
    S[2] = Fraction(m * p, 4)

    for k in range(3, p + 1, 2):
        if k not in ck:
            break
        # k*c_k = m^k + 2*sum_{r=0}^{floor(k/2)} C(k,2r)*(-1/2)^{k-2r}*(-1)^r * S_{2r}
        lhs = Fraction(k * ck[k]) - Fraction(m)**k

        # Subtract known contributions (r=0,...,(k-3)/2)
        known_sum = Fraction(0)
        for r in range(0, (k-1)//2):
            coeff = Fraction(comb(k, 2*r)) * Fraction(-1, 2)**(k - 2*r) * (-1)**r
            if 2*r in S:
                known_sum += coeff * S[2*r]

        # The new term is r = (k-1)/2 (the highest r in the sum)
        r_new = (k - 1) // 2
        coeff_new = Fraction(comb(k, 2*r_new)) * Fraction(-1, 2)**(k - 2*r_new) * (-1)**r_new

        # lhs = 2*(known_sum + coeff_new * S_{k-1})
        # => S_{k-1} = (lhs/2 - known_sum) / coeff_new
        S_new = (lhs / 2 - known_sum) / coeff_new
        S[2*r_new] = S_new

    return S


def analyze_p17_exact():
    """Exact arithmetic analysis at p=17."""
    p = 17
    m = 8
    N = 256

    print(f"{'='*70}")
    print(f"EXACT ARITHMETIC at p={p}, m={m}")
    print(f"{'='*70}")

    t0 = time.time()
    data = []
    by_D2 = defaultdict(list)

    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj_int(p, S)

        # Exact cycle counts
        ck = compute_ck_exact(A, p)

        # Exact moments from ck inversion
        moments = invert_moments_from_ck(p, ck)

        D2 = compute_D2_multiset(p, S)
        alpha_1 = sum(ck.values())

        row = {
            'bits': bits, 'S': S, 'D2': D2,
            'ck': ck, 'moments_exact': moments,
            'alpha_1': alpha_1
        }
        data.append(row)
        by_D2[D2].append(row)

        if bits % 32 == 0:
            elapsed = time.time() - t0
            print(f"  Progress: {bits}/{N} ({elapsed:.0f}s)", flush=True)

    t1 = time.time()
    n_types = len(by_D2)
    print(f"\n  Computed {N} orientations in {t1-t0:.1f}s")
    print(f"  {n_types} unique D^2 multisets")

    # Verify exact c_k linear dependency
    print(f"\n  --- EXACT c_k VERIFICATION ---")
    for k in range(5, p + 1, 2):
        ck_vals = [d['ck'][k] for d in data]
        if len(set(ck_vals)) == 1:
            print(f"  c_{k}: CONSTANT = {ck_vals[0]}")
            continue

        # c_k should be exact linear in S_4,...,S_{k-1}
        n_params = (k - 3) // 2
        # Build exact linear system using Fraction
        # Use one representative per D^2 type
        types = sorted(by_D2.keys())
        X_rows = []
        y_vals = []
        for D2 in types:
            d = by_D2[D2][0]
            row = [d['moments_exact'].get(2*j, Fraction(0)) for j in range(2, 2+n_params)]
            row.append(Fraction(1))
            X_rows.append(row)
            y_vals.append(Fraction(d['ck'][k]))

        # Check if system is consistent using exact arithmetic
        # If n_types > n_params+1, it's overconstrained
        if len(types) <= n_params + 1:
            print(f"  c_{k}: {len(types)} types, {n_params+1} params (trivially exact)")
        else:
            # Check by solving with first n_params+1 equations, verify rest
            from fractions import Fraction as F
            n_eq = len(types)
            n_var = n_params + 1
            aug = [X_rows[i][:] + [y_vals[i]] for i in range(n_eq)]

            # Forward elimination on first n_var rows
            for col in range(n_var):
                pivot = None
                for row in range(col, min(n_var, n_eq)):
                    if aug[row][col] != 0:
                        pivot = row
                        break
                if pivot is None:
                    print(f"  c_{k}: SINGULAR at col {col}")
                    break
                aug[col], aug[pivot] = aug[pivot], aug[col]
                for row in range(col+1, n_eq):
                    if aug[row][col] != 0:
                        factor = aug[row][col] / aug[col][col]
                        for j in range(n_var + 1):
                            aug[row][j] -= factor * aug[col][j]

            # Back substitution
            coeffs = [F(0)] * n_var
            for col in range(n_var-1, -1, -1):
                s = aug[col][n_var]
                for j in range(col+1, n_var):
                    s -= aug[col][j] * coeffs[j]
                coeffs[col] = s / aug[col][col]

            # Verify remaining equations
            max_err = F(0)
            for i in range(n_var, n_eq):
                pred = sum(coeffs[j] * X_rows[i][j] for j in range(n_var))
                err = abs(pred - y_vals[i])
                max_err = max(max_err, err)

            print(f"  c_{k}: {n_eq} types, {n_var} params, "
                  f"max_err = {float(max_err):.1f} "
                  f"{'EXACT' if max_err == 0 else ''}")

    # Alpha_1 exact test
    print(f"\n  --- EXACT alpha_1 LINEAR TEST ---")
    types = sorted(by_D2.keys())
    n_eq = len(types)

    # All m-1 moments: S_4, S_6, ..., S_{2m}
    n_var = m  # m-1 moments + constant
    X_rows = []
    y_vals = []
    for D2 in types:
        d = by_D2[D2][0]
        row = [d['moments_exact'].get(2*j, Fraction(0)) for j in range(2, m+1)]
        row.append(Fraction(1))
        X_rows.append(row)
        y_vals.append(Fraction(d['alpha_1']))

    aug = [X_rows[i][:] + [y_vals[i]] for i in range(n_eq)]

    # Forward elimination
    for col in range(n_var):
        pivot = None
        for row in range(col, n_eq):
            if aug[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            print(f"  SINGULAR at col {col}")
            break
        aug[col], aug[pivot] = aug[pivot], aug[col]
        for row in range(col+1, n_eq):
            if aug[row][col] != 0:
                factor = aug[row][col] / aug[col][col]
                for j in range(n_var + 1):
                    aug[row][j] -= factor * aug[col][j]

    # Back substitution
    coeffs = [Fraction(0)] * n_var
    for col in range(n_var-1, -1, -1):
        s = aug[col][n_var]
        for j in range(col+1, n_var):
            s -= aug[col][j] * coeffs[j]
        coeffs[col] = s / aug[col][col]

    # Verify ALL equations
    max_err = Fraction(0)
    for i in range(n_eq):
        pred = sum(coeffs[j] * X_rows[i][j] for j in range(n_var))
        err = abs(pred - y_vals[i])
        max_err = max(max_err, err)

    print(f"  alpha_1 = linear(S4,...,S{2*m}): {n_eq} types, {n_var} params")
    print(f"  Max error = {float(max_err):.4f}")
    if max_err == 0:
        print(f"  *** GENUINE LINEAR IDENTITY CONFIRMED! ***")
        print(f"  Coefficients:")
        for j in range(m-1):
            print(f"    S{2*(j+2)}: {coeffs[j]} = {float(coeffs[j]):.8f}")
        print(f"    const: {coeffs[-1]} = {float(coeffs[-1]):.2f}")
    else:
        print(f"  Linear identity FAILS with exact arithmetic")
        print(f"  First few residuals:")
        for i in range(min(5, n_eq)):
            pred = sum(coeffs[j] * X_rows[i][j] for j in range(n_var))
            err = pred - y_vals[i]
            print(f"    Type {i}: alpha_1={y_vals[i]}, pred={pred}, err={err}")

    return data, by_D2


if __name__ == '__main__':
    analyze_p17_exact()
    print("\nDONE.")
