#!/usr/bin/env python3
"""Exact companion for THM-2601.

The input is the thirteen septimal Bockstein factors printed by THM-2585.
All calculations are over F_13.  The script verifies an explicit scalar
sheet coordinate, its nonlinear translation law, all six owner-conjugate
coordinates, a full-rank affine-covariance obstruction, and the sharp CRT
component collisions.
"""

from collections import Counter
from itertools import combinations


P = 13
Q7 = 7

Y_RAW = [
    (0, 9, 9, 0, 0, 4, 4),
    (5, 9, 7, 11, 11, 4, 9),
    (3, 11, 7, 10, 10, 11, 1),
    (9, 11, 4, 2, 10, 5, 10),
    (11, 5, 4, 10, 12, 11, 8),
    (11, 1, 11, 2, 8, 1, 1),
    (6, 9, 12, 8, 4, 4, 7),
    (7, 6, 9, 9, 5, 1, 4),
    (2, 12, 12, 5, 11, 2, 12),
    (2, 5, 2, 1, 3, 9, 8),
    (4, 3, 8, 3, 11, 9, 2),
    (10, 12, 2, 3, 3, 6, 2),
    (8, 4, 9, 2, 2, 6, 4),
]

LAMBDA = (1, 11, 8, 1, 0, 1)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right)) % P


def reduce_phi7(coefficients):
    """Reduce modulo Phi_7=1+z+...+z^6 into degrees 0,...,5."""
    out = [value % P for value in coefficients]
    if len(out) < 6:
        out += [0] * (6 - len(out))
    for degree in range(len(out) - 1, 5, -1):
        value = out[degree] % P
        if not value:
            continue
        # z^degree*Phi_7 shifted by degree-6 is zero.
        for j in range(7):
            out[degree - 6 + j] = (out[degree - 6 + j] - value) % P
    return tuple(out[:6])


def sigma(vector, kappa):
    expanded = [0] * (5 * kappa + 1)
    for degree, value in enumerate(vector):
        expanded[degree * kappa] = value
    return reduce_phi7(expanded)


def poly_eval(poly, x):
    ans = 0
    power = 1
    for coefficient in poly:
        ans = (ans + coefficient * power) % P
        power = power * x % P
    return ans


def interpolate(values):
    """Degree-at-most-12 polynomial f with f(x)=values[x] on F_13."""
    out = [0] * P
    for x, target in enumerate(values):
        basis = [1]
        denominator = 1
        for y in range(P):
            if y == x:
                continue
            nxt = [0] * (len(basis) + 1)
            for j, value in enumerate(basis):
                nxt[j] = (nxt[j] - y * value) % P
                nxt[j + 1] = (nxt[j + 1] + value) % P
            basis = nxt
            denominator = denominator * (x - y) % P
        scale = target * pow(denominator, -1, P) % P
        for j, value in enumerate(basis):
            out[j] = (out[j] + scale * value) % P
    return tuple(out)


def determinant(matrix):
    a = [list(row) for row in matrix]
    n = len(a)
    det = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col] % P), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        det = det * a[col][col] % P
        inv = pow(a[col][col] % P, -1, P)
        a[col] = [value * inv % P for value in a[col]]
        for row in range(col + 1, n):
            factor = a[row][col] % P
            if factor:
                a[row] = [(x - factor * y) % P
                          for x, y in zip(a[row], a[col])]
    return det % P


def rank(matrix):
    a = [list(row) for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    pivot_row = 0
    for col in range(cols):
        pivot = next((r for r in range(pivot_row, rows) if a[r][col] % P), None)
        if pivot is None:
            continue
        a[pivot_row], a[pivot] = a[pivot], a[pivot_row]
        inv = pow(a[pivot_row][col] % P, -1, P)
        a[pivot_row] = [value * inv % P for value in a[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = a[row][col] % P
            if factor:
                a[row] = [(x - factor * y) % P
                          for x, y in zip(a[row], a[pivot_row])]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def mul_quadratic(a, b, middle):
    # t^2=-middle*t-1 modulo t^2+middle*t+1.
    a0, a1 = a
    b0, b1 = b
    return ((a0 * b0 - a1 * b1) % P,
            (a0 * b1 + a1 * b0 - middle * a1 * b1) % P)


def eval_quadratic(vector, middle):
    ans = (0, 0)
    power = (1, 0)
    for value in vector:
        ans = ((ans[0] + value * power[0]) % P,
               (ans[1] + value * power[1]) % P)
        power = mul_quadratic(power, (0, 1), middle)
    return ans


def main():
    factors = [reduce_phi7(row) for row in Y_RAW]
    require(factors == [tuple((row[i] - row[6]) % P for i in range(6))
                        for row in Y_RAW], "Phi_7 reduction changed")

    sheet_values = tuple(dot(LAMBDA, factor) for factor in factors)
    expected_values = (9, 3, 10, 1, 8, 0, 7, 2, 12, 11, 4, 6, 5)
    require(sheet_values == expected_values, "linear sheet coordinate changed")
    require(set(sheet_values) == set(range(P)), "sheet coordinate is not bijective")

    q_by_value = [None] * P
    for q, value in enumerate(sheet_values):
        q_by_value[value] = q
    q_by_value = tuple(q_by_value)
    inverse_poly = interpolate(q_by_value)
    expected_inverse = (5, 4, 3, 0, 7, 8, 12, 6, 11, 6, 0, 6, 0)
    require(inverse_poly == expected_inverse, "inverse polynomial changed")
    require(all(poly_eval(inverse_poly, value) == q
                for q, value in enumerate(sheet_values)),
            "inverse polynomial failed")

    successor_values = tuple(
        sheet_values[(q_by_value[value] + 1) % P] for value in range(P)
    )
    expected_successor_values = (7, 8, 12, 10, 6, 9, 5, 2, 0, 3, 1, 4, 11)
    require(successor_values == expected_successor_values,
            "successor permutation changed")
    successor_poly = interpolate(successor_values)
    expected_successor = (7, 4, 0, 1, 6, 4, 6, 1, 7, 6, 3, 2, 0)
    require(successor_poly == expected_successor, "successor polynomial changed")
    require(all(poly_eval(successor_poly, value) == successor_values[value]
                for value in range(P)), "successor polynomial failed")
    require(not any(all(successor_values[x] == (a * x + b) % P
                        for x in range(P))
                    for a in range(P) for b in range(P)),
            "successor unexpectedly became affine")
    for start in range(P):
        value = start
        orbit = []
        for _ in range(P):
            orbit.append(value)
            value = successor_values[value]
        require(value == start and len(set(orbit)) == P,
                "successor is not one 13-cycle")

    # Every owner automorphism has a conjugate scalar coordinate recovering
    # the same sheet value.
    owner_lambdas = []
    for kappa in range(1, Q7):
        inverse_kappa = pow(kappa, -1, Q7)
        conjugate = tuple(dot(LAMBDA, sigma(tuple(int(i == j) for i in range(6)),
                                                  inverse_kappa))
                          for j in range(6))
        owner_lambdas.append(conjugate)
        require(tuple(dot(conjugate, sigma(factor, kappa)) for factor in factors)
                == sheet_values, "owner-conjugate coordinate failed")

    # If mu(Y_q)=alpha*q+beta for all q, then
    # (mu_0,...,mu_5,alpha,beta) lies in the kernel of this 13x8 matrix.
    affine_matrix = [factor + ((-q) % P, (-1) % P)
                     for q, factor in enumerate(factors)]
    affine_rank = rank(affine_matrix)
    require(affine_rank == 8, "affine-covariance obstruction lost full rank")
    minor_rows = None
    minor_det = None
    for rows in combinations(range(P), 8):
        det = determinant([affine_matrix[row] for row in rows])
        if det:
            minor_rows, minor_det = rows, det
            break
    require(minor_rows == (0, 1, 2, 3, 4, 5, 6, 7) and minor_det == 4,
            "affine no-go minor certificate changed")

    # CRT hostile: two quadratic components separate every q, while the
    # third has precisely two double fibres.
    crt_data = {}
    for middle in (3, 5, 6):
        fibres = {}
        for q, factor in enumerate(factors):
            fibres.setdefault(eval_quadratic(factor, middle), []).append(q)
        crt_data[middle] = tuple(sorted(tuple(v) for v in fibres.values()))
    expected_crt = {
        3: tuple((q,) for q in range(P)),
        5: tuple((q,) for q in range(P)),
        6: ((0,), (1, 9), (2,), (3,), (4, 12), (5,), (6,), (7,), (8,),
            (10,), (11,)),
    }
    require(crt_data == expected_crt, "CRT collision boundary changed")

    # Sharp additive boundary.  A separating scalar coordinate is faithful on
    # singleton sheets but cannot be a Boolean-section charge detector.
    subset_sums = Counter()
    for mask in range(1 << P):
        subset_sums[sum(sheet_values[q] for q in range(P) if mask >> q & 1) % P] += 1
    require(subset_sums == Counter({0: 632, **{value: 630 for value in range(1, P)}}),
            "Boolean scalar subset-sum histogram changed")
    require(sheet_values[5] == 0, "singleton scalar-zero hostile changed")

    print("== exact THM-2601 linear Bockstein sheet-coordinate probe ==")
    print("lambda power-basis row:", LAMBDA)
    print("lambda(Y_q), q=0..12:", sheet_values)
    print("inverse q=P(t) coefficients low-to-high:", inverse_poly, "degree 11")
    print("successor t(q)->t(q+1) values:", successor_values)
    print("successor polynomial coefficients low-to-high:", successor_poly,
          "degree 11")
    print("successor orbit: one 13-cycle; affine scalar realizations: 0/169")
    print("owner-conjugate lambda rows:", tuple(owner_lambdas))
    print("affine-covariant linear system: rank", affine_rank,
          "certificate rows", minor_rows, "det", minor_det)
    print("CRT component fibres:", crt_data)
    print("Boolean scalar subset sums: zero=632 (631 nonempty), each nonzero=630;")
    print("singleton scalar-zero hostile: q=5, while Y_5 is a THM-2585 unit")
    print("exact checks: separator/inverse/successor/owner/affine/CRT/subsets PASS")


if __name__ == "__main__":
    main()
