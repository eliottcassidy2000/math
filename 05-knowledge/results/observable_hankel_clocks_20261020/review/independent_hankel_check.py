#!/usr/bin/env python3
"""Independent construction/check of the observable Hankel precision packet.

This imports no root checker. SymPy supplies symbolic polynomials and Smith
forms; the determinant is independently computed with integer Bareiss steps.
"""

import hashlib
import json
import math
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def a(n):
    if n % 3 == 0:
        return -8*n**4 + 16*n**5 + 32*n**6 - 64*n**7
    if n % 3 == 1:
        return 8*n**4 + 16*n**5 - 96*n**6 - 64*n**7
    return -n + 12*n**2 - 60*n**3 + 160*n**4 - 240*n**5 + 192*n**6 - 64*n**7


def valuation(z, p):
    require(z != 0, "No zero Smith invariant expected")
    z, v = abs(int(z)), 0
    while z % p == 0:
        z //= p
        v += 1
    return v


def bareiss(matrix):
    b = [list(map(int, row)) for row in matrix]
    n, previous, sign = len(b), 1, 1
    for k in range(n-1):
        if b[k][k] == 0:
            row = next((j for j in range(k+1, n) if b[j][k]), None)
            if row is None:
                return 0
            b[k], b[row] = b[row], b[k]
            sign = -sign
        pivot = b[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = pivot*b[i][j] - b[i][k]*b[k][j]
                require(numerator % previous == 0, "Bareiss exact division")
                b[i][j] = numerator // previous
        for i in range(k+1, n):
            b[i][k] = 0
        previous = pivot
    return sign*b[-1][-1]


def rank_mod(matrix, p):
    b = [[int(z) % p for z in row] for row in matrix]
    r, rows, cols = 0, len(b), len(b[0])
    for col in range(cols):
        pivot = next((j for j in range(r, rows) if b[j][col]), None)
        if pivot is None:
            continue
        b[r], b[pivot] = b[pivot], b[r]
        inv = pow(b[r][col], -1, p)
        b[r] = [(z*inv) % p for z in b[r]]
        for i in range(rows):
            if i != r:
                coefficient = b[i][col]
                b[i] = [(u-coefficient*v) % p for u, v in zip(b[i], b[r])]
        r += 1
        if r == rows:
            break
    return r


def main():
    n, x = sp.symbols("n x")
    polynomials = [
        -8*n**4+16*n**5+32*n**6-64*n**7,
        8*n**4+16*n**5-96*n**6-64*n**7,
        -n+12*n**2-60*n**3+160*n**4-240*n**5+192*n**6-64*n**7,
    ]
    chi = sp.Poly((x-1)**8 * (x*x+x+1)**7, x)
    coeffs = [int(chi.nth(j)) for j in range(23)]
    for residue in range(3):
        residual = sp.expand(sum(coeffs[j]*polynomials[(residue+j) % 3].subs(n, n+j)
                                 for j in range(23)))
        require(residual == 0, "Global symbolic recurrence in each residue lane")

    matrix = [[a(i+j) for j in range(22)] for i in range(22)]
    H = sp.Matrix(matrix)
    determinant = bareiss(matrix)
    expected_det = -(2**206)*(3**191)*(5**22)*(7**8)*(61**7)
    require(determinant == expected_det, "Independent normalized Hankel determinant")

    S = sp.zeros(22)
    for i in range(21):
        S[i, i+1] = 1
    for j in range(22):
        S[21, j] = -coeffs[j]
    C = S.T
    require(H*C == S*H, "Hankel-companion intertwining")
    require(int(S.det()) == 1, "Unimodular integer shift")

    D = smith_normal_form(H, domain=ZZ)
    diagonal = [abs(int(D[i, i])) for i in range(22)]
    require(all(diagonal[i+1] % diagonal[i] == 0 for i in range(21)), "Smith divisibility chain")
    require(math.prod(diagonal) == abs(determinant), "Smith determinant product")
    expected = {
        2: [0]*6+[6]*6+[14]*2+[15]*4+[20]*2+[21]*2,
        3: [0]*3+[2]*3+[4]*3+[8]*3+[10]*3+[12]*3+[19,20,22,22],
        5: [0]*14+[2]*4+[3]*2+[4]*2,
        7: [0]*21+[8],
        61: [0]*21+[7],
    }
    valuations, ranks = {}, {}
    for p, values in expected.items():
        got = [valuation(d, p) for d in diagonal]
        require(got == values, "Published prime valuation list")
        valuations[p] = got
        ranks[p] = rank_mod(matrix, p)
        require(ranks[p] == got.count(0), "Independent modular rank")
    for p in [11, 13, 17, 19]:
        require(rank_mod(matrix, p) == 22, "Good-prime full rank control")

    prefixes = {}
    for r in [6, 12, 14, 18, 20, 22]:
        dr = smith_normal_form(H[:, :r], domain=ZZ)
        vals = [valuation(dr[i, i], 2) for i in range(r)]
        require(vals == expected[2][:r], "Prefix realizes the active Smith directions")
        prefixes[r] = vals
    staircase = []
    for k in range(1, 25):
        active = sum(v < k for v in expected[2])
        full_log2_size = sum(max(k-v, 0) for v in expected[2])
        require(active in prefixes, "Active layer has a checked prefix")
        prefix_log2_size = sum(max(k-v, 0) for v in prefixes[active])
        require(prefix_log2_size == full_log2_size, "Prefix image equals full image")
        staircase.append({"precision": k, "minimal_monic_recurrence_order": active,
                          "log2_observable_size": full_log2_size})

    period_tests = []
    for k in range(1, 13):
        m, period = 2**k, 3*2**k
        values = [a(i) % m for i in range(2*period)]
        require(values[:period] == values[period:], "Claimed period finite window")
        for divisor in sp.divisors(period):
            if divisor < period:
                require(any(values[j] != values[j+divisor] for j in range(period)),
                        "Every proper divisor period fails")
        period_tests.append({"precision": k, "least_period": period})

    H6 = H[:6, :6]
    b6 = sp.Matrix([a(i+6) for i in range(6)])
    six_register = []
    for modulus, target in [(64, [7,48,16,2,16,48]), (128, [7,112,16,66,80,48])]:
        candidate = [int(z) % modulus for z in H6.inv_mod(modulus)*b6]
        require(candidate == target, "Unique six-register candidate")
        residuals = [(a(i+6)-sum(candidate[j]*a(i+j) for j in range(6))) % modulus
                     for i in range(22)]
        if modulus == 64:
            require(not any(residuals), "First 22 residuals vanish; recurrence propagates")
        else:
            require(residuals[6] != 0, "Precision-seven obstruction at index six")
        six_register.append({"modulus": modulus, "coefficients": candidate,
                             "first_22_residuals": residuals})

    report = {
        "status": "PASS: independently constructed exact Hankel data and typed completion audit",
        "global_recurrence_symbolic_lanes": 3,
        "rational_minimal_order_by_nonzero_hankel_determinant": 22,
        "characteristic_polynomial_coefficients_constant_first": coeffs,
        "hankel_determinant": determinant,
        "hankel_companion_intertwining": True,
        "integer_shift_determinant": 1,
        "smith_diagonal_absolute": diagonal,
        "prime_valuations": valuations,
        "independent_modular_ranks": ranks,
        "prefix_smith_2valuations": prefixes,
        "minimal_monic_recurrence_order_staircase_checks": staircase,
        "finite_period_checks": period_tests,
        "six_register_candidates": six_register,
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
