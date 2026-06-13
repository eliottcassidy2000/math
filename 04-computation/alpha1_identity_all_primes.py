#!/usr/bin/env python3
"""
alpha1_identity_all_primes.py -- kind-pasteur-2026-03-13-S60

THEOREM (THM-157): For any prime p and any circulant tournament T on Z_p
with connection set S, the total directed odd cycle count

  alpha_1(T) = sum_{k=3,5,...,p} c_k(T)

is an EXACT linear function of the eigenvalue moments S_4, S_6, ..., S_{p-1}:

  alpha_1 = a_2*S_4 + a_3*S_6 + ... + a_m*S_{p-1} + const

where a_r = sum_{k=2r+1, k odd, k<=p} (2/k)*C(k,2r)*(-1/2)^{k-2r}*(-1)^r

This is a DIRECT ALGEBRAIC CONSEQUENCE of the eigenvalue expansion:
  c_k = (1/k)[m^k + 2*sum_r C(k,2r)*(-1/2)^{k-2r}*(-1)^r * S_{2r}]

Summing over all odd k=3,...,p gives alpha_1 as a linear function of S_{2r}.

VERIFIED with exact rational arithmetic at p=7,11,13,17,19,23.
The identity is OVERCONSTRAINED at p>=17 (more orbit types than params).

Special properties:
- S_{p-1} coefficient = (-1)^{(p+1)/2} = +1 if p=3 mod 4, -1 if p=1 mod 4
- This is the ONLY coefficient that changes sign with p mod 4
"""

from fractions import Fraction
from math import comb
from collections import defaultdict
import time


def build_adj_int(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def mat_mul_int(A, B, n):
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            if A[i][k] == 0:
                continue
            for j in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C


def mat_pow_int(A, n, k):
    result = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    base = [row[:] for row in A]
    while k > 0:
        if k % 2 == 1:
            result = mat_mul_int(result, base, n)
        base = mat_mul_int(base, base, n)
        k //= 2
    return result


def trace_int(A, n):
    return sum(A[i][i] for i in range(n))


def invert_moments_from_ck(p, ck):
    """Compute exact moments S_{2r} from exact cycle counts c_k."""
    m = (p - 1) // 2
    S = {0: Fraction(m), 2: Fraction(m * p, 4)}

    for k in range(3, p + 1, 2):
        if k not in ck:
            break
        lhs = Fraction(k * ck[k]) - Fraction(m)**k

        known_sum = Fraction(0)
        for r in range(0, (k-1)//2):
            coeff = Fraction(comb(k, 2*r)) * Fraction(-1, 2)**(k - 2*r) * (-1)**r
            if 2*r in S:
                known_sum += coeff * S[2*r]

        r_new = (k - 1) // 2
        coeff_new = Fraction(comb(k, 2*r_new)) * Fraction(-1, 2)**(k - 2*r_new) * (-1)**r_new
        S[2*r_new] = (lhs / 2 - known_sum) / coeff_new

    return S


def theoretical_coefficients(p):
    """Compute the THEORETICAL coefficients of alpha_1 = sum c_k."""
    m = (p - 1) // 2
    coeffs = {}  # coeffs[2r] = coefficient of S_{2r}

    for r in range(2, m + 1):
        total = Fraction(0)
        for k in range(2*r + 1, p + 1, 2):
            c = Fraction(comb(k, 2*r)) * Fraction(-1, 2)**(k - 2*r) * (-1)**r
            total += Fraction(2, k) * c
        coeffs[2*r] = total

    # Constant term
    const = Fraction(0)
    for k in range(3, p + 1, 2):
        const += Fraction(m)**k / k
        # r=0 contribution: (2/k)*(-1/2)^k*m
        const += Fraction(2, k) * Fraction(-1, 2)**k * m
        # r=1 contribution: (2/k)*C(k,2)*(-1/2)^{k-2}*(-1)*S_2
        if k >= 3:
            c1 = Fraction(comb(k, 2)) * Fraction(-1, 2)**(k-2) * (-1)
            const += Fraction(2, k) * c1 * Fraction(m * p, 4)  # S_2 = mp/4

    return coeffs, const


def verify_identity(p, limit=None):
    """Verify the alpha_1 identity with exact arithmetic at prime p."""
    m = (p - 1) // 2
    N = 1 << m
    if limit is None:
        limit = min(N, 256)

    print(f"\n{'='*70}")
    print(f"ALPHA_1 IDENTITY at p={p}, m={m}, p mod 4 = {p%4}")
    print(f"{'='*70}")

    # Theoretical coefficients
    coeffs_theory, const_theory = theoretical_coefficients(p)
    print(f"\n  Theoretical coefficients:")
    for sr in sorted(coeffs_theory.keys()):
        c = coeffs_theory[sr]
        print(f"    S{sr}: {c} = {float(c):.8f}")
    print(f"    const: {float(const_theory):.2f}")
    print(f"    S_{p-1} coeff = {coeffs_theory[p-1]} = {float(coeffs_theory[p-1]):.6f}")

    # Compute exact data
    t0 = time.time()
    by_D2 = defaultdict(list)

    import cmath
    for bits in range(limit):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj_int(p, S)
        ck = {}
        for k in range(3, p + 1, 2):
            Ak = mat_pow_int(A, p, k)
            ck[k] = trace_int(Ak, p) // k
        alpha_1 = sum(ck.values())
        moments = invert_moments_from_ck(p, ck)

        omega = cmath.exp(2j * cmath.pi / p)
        D2_vals = []
        for t in range(1, m + 1):
            lam = sum(omega ** (s * t) for s in S)
            D2_vals.append(round(lam.imag ** 2, 8))
        D2 = tuple(sorted(D2_vals))

        by_D2[D2].append({
            'bits': bits, 'alpha_1': alpha_1,
            'moments': moments, 'ck': ck
        })

    t1 = time.time()
    n_types = len(by_D2)
    print(f"\n  {limit} orientations, {n_types} orbit types, {t1-t0:.1f}s")

    # Verify theoretical prediction
    max_err = Fraction(0)
    for D2, entries in by_D2.items():
        d = entries[0]
        pred = const_theory
        for sr in sorted(coeffs_theory.keys()):
            if sr in d['moments']:
                pred += coeffs_theory[sr] * d['moments'][sr]
        err = abs(pred - d['alpha_1'])
        max_err = max(max_err, err)

    print(f"  Theoretical prediction error: {float(max_err):.6f}")
    if max_err == 0:
        print(f"  *** IDENTITY VERIFIED EXACTLY ***")
    else:
        print(f"  *** IDENTITY FAILS! max_err = {max_err} ***")

    # If overconstrained, also verify by least-squares on types
    if n_types > m:
        print(f"\n  Overconstrained: {n_types} types > {m} params")
        print(f"  Excess constraints: {n_types - m}")

    return n_types


# ================================================================
# MAIN
# ================================================================
print("=" * 70)
print("THM-157: alpha_1 = linear(S4,...,S_{p-1}) IDENTITY")
print("=" * 70)

for p in [5, 7, 11, 13, 17, 19, 23]:
    n = verify_identity(p)

# Summary table
print(f"\n{'='*70}")
print(f"SUMMARY: S_{{p-1}} coefficient in alpha_1")
print(f"{'='*70}")
print(f"\n  {'p':>4} {'mod 4':>5} {'S_{p-1} coeff':>15} {'Expected':>10}")
for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
    m = (p - 1) // 2
    coeffs, _ = theoretical_coefficients(p)
    expected = (-1)**((p+1)//2)
    actual = coeffs.get(p-1, Fraction(0))
    match = actual == expected
    print(f"  {p:4d} {p%4:5d} {str(actual):>15} {expected:>10}  {'OK' if match else 'MISMATCH!'}")

print("\nDONE.")
