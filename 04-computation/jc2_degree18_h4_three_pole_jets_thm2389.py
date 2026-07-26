#!/usr/bin/env python3
"""Exact companion for THM-2389.

Dependency-free checks for the H4 three-pole Hermite-jet synchronization.
Every failed check raises explicitly, including under ``python -O``.
"""

from fractions import Fraction as Q
from itertools import permutations
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def add(left, right, limit=None):
    size = max(len(left), len(right))
    if limit is not None:
        size = min(size, limit)
    out = [Q(0) for _ in range(size)]
    for index in range(size):
        if index < len(left):
            out[index] += left[index]
        if index < len(right):
            out[index] += right[index]
    return trim(out)


def scale(poly, scalar, limit=None):
    size = len(poly) if limit is None else min(len(poly), limit)
    return trim([scalar * poly[index] for index in range(size)])


def multiply(left, right, limit=None):
    size = len(left) + len(right) - 1
    if limit is not None:
        size = min(size, limit)
    out = [Q(0) for _ in range(size)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            if i + j >= size:
                break
            out[i + j] += left_value * right_value
    return trim(out)


def power(poly, exponent, limit=None):
    out = [Q(1)]
    base = list(poly)
    value = exponent
    while value:
        if value & 1:
            out = multiply(out, base, limit)
        value //= 2
        if value:
            base = multiply(base, base, limit)
    return out


def evaluate(poly, value):
    total = Q(0)
    for coefficient in reversed(poly):
        total = total * value + coefficient
    return total


def shifted(poly, shift):
    """Return coefficients of poly(x+shift)."""
    out = [Q(0) for _ in poly]
    for degree, coefficient in enumerate(poly):
        for index in range(degree + 1):
            out[index] += coefficient * Q(comb(degree, index)) * shift ** (
                degree - index
            )
    return trim(out)


def matrix_rank(matrix):
    rows = [list(map(Q, row)) for row in matrix]
    if not rows:
        return 0
    row_count = len(rows)
    column_count = len(rows[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        pivot_value = rows[pivot_row][column]
        rows[pivot_row] = [value / pivot_value for value in rows[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                rows[row][entry] - factor * rows[pivot_row][entry]
                for entry in range(column_count)
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def matrix_times_vector(matrix, vector):
    return [
        sum((Q(entry) * Q(value) for entry, value in zip(row, vector)), Q(0))
        for row in matrix
    ]


ALPHA = Q(16, 964467)
BETA = Q(64, 703096443)
F_POLY = [Q(704, 14348907), Q(80, 19683), Q(0), Q(1)]


def f(value):
    return evaluate(F_POLY, value)


def f_prime(value):
    return 3 * value * value + Q(80, 19683)


def check_slope_cubic():
    require(ALPHA * 245 == Q(80, 19683), "alpha slope coefficient")
    require(BETA * 539 == Q(704, 14348907), "beta slope coefficient")

    l_poly = [Q(1127), Q(-138915), Q(1607445), Q(-26040609)]
    transformed = scale(
        shifted(l_poly, Q(5, 243)), Q(-1, 26040609)
    )
    require(transformed == F_POLY, "depressed L-to-f translation")

    discriminant = -4 * F_POLY[1] ** 3 - 27 * F_POLY[0] ** 2
    require(
        discriminant == Q(-94208, 282429536481),
        "slope cubic discriminant",
    )

    g2_zero = Q(-8, 243)
    g4_zero = Q(4, 243)
    require(f(g2_zero) == Q(-64, 531441), "g2 hostile value")
    require(f(g4_zero) == Q(64, 531441), "g4 hostile value")

    require(
        ALPHA * 1890 == Q(160 * 243, 1240029),
        "g2 linear coefficient",
    )
    require(
        BETA * 11340 == Q(160 * 8, 1240029),
        "g2 constant coefficient",
    )
    require(BETA * 183708 == Q(256, 15309), "h3 pivot")
    require(
        ALPHA * 122472 == Q(128 * 243, 15309),
        "g4 linear coefficient",
    )
    require(
        -BETA * 367416 == Q(-128 * 4, 15309),
        "g4 constant coefficient",
    )
    require(BETA * 2480058 == Q(128, 567), "h5 pivot")
    return discriminant


def hermite_coefficients(a, slopes, kappa):
    a0, a1, a2, a3 = map(Q, a)
    s0, s1, sinf = map(Q, slopes)
    a_at_one = a0 + a1 + a2 + a3
    derivative_at_one = a1 + 2 * a2 + 3 * a3
    c0 = s0 * a0 * a0
    c1 = 2 * s0 * a0 * a1
    c6 = sinf * a3 * a3
    c5 = 2 * sinf * a3 * a2
    sigma = s1 * a_at_one * a_at_one - (c0 + c1 + c5 + c6)
    tau = (
        2 * s1 * a_at_one * derivative_at_one
        - (c1 + 5 * c5 + 6 * c6)
    )
    return [
        c0,
        c1,
        3 * sigma - tau + kappa,
        tau - 2 * sigma - 2 * kappa,
        kappa,
        c5,
        c6,
    ]


def check_hermite_family():
    constraint_matrix = [
        [1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1],
        [0, 1, 2, 3, 4, 5, 6],
        [0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1],
    ]
    kernel = [0, 0, 1, -2, 1, 0, 0]
    require(matrix_rank(constraint_matrix) == 6, "Hermite constraint rank")
    require(
        matrix_times_vector(constraint_matrix, kernel) == [Q(0)] * 6,
        "Hermite kernel B2^2",
    )

    rational_controls = 0
    a_samples = [
        (2, 3, -1, 5),
        (-3, 4, 2, 7),
        (5, -2, 6, 3),
        (7, 1, -4, 2),
    ]
    slope_samples = [
        (Q(1, 2), Q(-2, 3), Q(5, 7)),
        (Q(-4, 5), Q(3, 8), Q(9, 11)),
    ]
    for a in a_samples:
        a0, a1, a2, a3 = map(Q, a)
        a_at_one = sum(a, Q(0))
        derivative_at_one = a1 + 2 * a2 + 3 * a3
        require(a0 * a3 * a_at_one != 0, "rational A pole units")
        for slopes in slope_samples:
            for kappa in (Q(-2), Q(0), Q(7, 3)):
                c = hermite_coefficients(a, slopes, kappa)
                s0, s1, sinf = slopes
                require(c[0] == s0 * a0 * a0, "Hermite at zero")
                require(c[1] == 2 * s0 * a0 * a1, "Hermite derivative zero")
                require(c[6] == sinf * a3 * a3, "Hermite at infinity")
                require(c[5] == 2 * sinf * a3 * a2, "Hermite infinity jet")
                require(
                    sum(c, Q(0)) == s1 * a_at_one * a_at_one,
                    "Hermite at one",
                )
                require(
                    sum((Q(index) * c[index] for index in range(7)), Q(0))
                    == 2 * s1 * a_at_one * derivative_at_one,
                    "Hermite derivative one",
                )
                rational_controls += 1

    prime = 59
    alpha_mod = (16 * pow(964467, -1, prime)) % prime
    beta_mod = (64 * pow(703096443, -1, prime)) % prime
    roots = [
        value
        for value in range(prime)
        if (
            value**3
            + alpha_mod * 245 * value
            + beta_mod * 539
        )
        % prime
        == 0
    ]
    require(roots == [5, 11, 43], "first split-prime slope roots")
    split_controls = 0
    a = (2, 3, 5, 7)
    for slopes in permutations(roots):
        for kappa in (0, 1, 7):
            c = hermite_coefficients(a, slopes, kappa)
            c_mod = [
                (value.numerator * pow(value.denominator, -1, prime)) % prime
                for value in c
            ]
            a0, a1, a2, a3 = a
            s0, s1, sinf = slopes
            a_at_one = sum(a) % prime
            derivative_at_one = (a1 + 2 * a2 + 3 * a3) % prime
            require(c_mod[0] == s0 * a0 * a0 % prime, "split zero value")
            require(c_mod[1] == 2 * s0 * a0 * a1 % prime, "split zero jet")
            require(c_mod[6] == sinf * a3 * a3 % prime, "split infinity")
            require(c_mod[5] == 2 * sinf * a3 * a2 % prime, "split infinity jet")
            require(sum(c_mod) % prime == s1 * a_at_one**2 % prime, "split one")
            require(
                sum(index * c_mod[index] for index in range(7)) % prime
                == 2 * s1 * a_at_one * derivative_at_one % prime,
                "split one jet",
            )
            split_controls += 1

    return rational_controls, split_controls


def coefficient(poly, index):
    return poly[index] if index < len(poly) else Q(0)


def check_jet_reconstruction():
    controls = 0
    samples = [
        (Q(2, 5), Q(3, 7), Q(-1, 3), Q(5, 4), Q(2, 9), Q(-7, 5)),
        (Q(-3, 8), Q(4, 9), Q(5, 6), Q(-2, 7), Q(11, 10), Q(3, 5)),
        (Q(7, 11), Q(-5, 12), Q(13, 9), Q(1, 6), Q(-4, 7), Q(8, 13)),
    ]
    parameter_samples = [
        (Q(2), Q(-3), Q(5), Q(7)),
        (Q(-1, 2), Q(4, 3), Q(-5, 6), Q(9, 5)),
        (Q(7, 4), Q(2, 9), Q(11, 8), Q(-13, 7)),
    ]

    for jet in samples:
        s, k, ell, m, n, o = jet
        r = [s, Q(0), k, ell, m, n, o]
        for bp, cp, dp, wp in parameter_samples:
            p_series = [
                Q(245),
                Q(0),
                1890 * bp,
                Q(0),
                -24300 * bp * bp + 122472 * dp,
                Q(0),
                Q(0),
            ]
            q_series = [
                Q(539),
                Q(0),
                11340 * bp,
                183708 * cp,
                72900 * bp * bp - 367416 * dp,
                2361960 * bp * cp + 2480058 * wp,
                Q(0),
            ]
            phi = add(
                add(power(r, 3, 7), scale(multiply(p_series, r, 7), ALPHA, 7), 7),
                scale(q_series, BETA, 7),
                7,
            )

            derivative = f_prime(s)
            expected = [
                f(s),
                Q(0),
                derivative * k
                + (1890 * ALPHA * s + 11340 * BETA) * bp,
                derivative * ell + 183708 * BETA * cp,
                derivative * m
                + 3 * s * k * k
                + ALPHA
                * (
                    1890 * bp * k
                    + (-24300 * bp * bp + 122472 * dp) * s
                )
                + BETA * (72900 * bp * bp - 367416 * dp),
                derivative * n
                + 6 * s * k * ell
                + 1890 * ALPHA * bp * ell
                + BETA * (2361960 * bp * cp + 2480058 * wp),
                derivative * o
                + 3 * s * (2 * k * m + ell * ell)
                + k**3
                + ALPHA
                * (
                    1890 * bp * m
                    + (-24300 * bp * bp + 122472 * dp) * k
                ),
            ]
            require(
                [coefficient(phi, index) for index in range(7)] == expected,
                "jet coefficient expansion",
            )
            controls += 1
    return controls


def check_global_residual_kernel():
    matrix = []

    # Order at least six at t=0.
    for degree in range(6):
        row = [Q(0)] * 19
        row[degree] = Q(1)
        matrix.append(row)

    # Order at least six at t=1: coefficients of u^j in N(1+u).
    for jet_order in range(6):
        row = [Q(0)] * 19
        for degree in range(jet_order, 19):
            row[degree] = Q(comb(degree, jet_order))
        matrix.append(row)

    # Order at least six at infinity after division by A^6:
    # the degree-18 through degree-13 coefficients vanish.
    for degree in range(13, 19):
        row = [Q(0)] * 19
        row[degree] = Q(1)
        matrix.append(row)

    require(len(matrix) == 18, "global condition count")
    rank = matrix_rank(matrix)
    require(rank == 18, "global residual rank")

    b6 = [Q(0)] * 19
    for index in range(7):
        b6[6 + index] = Q(comb(6, index) * (-1) ** (6 - index))
    require(
        matrix_times_vector(matrix, b6) == [Q(0)] * 18,
        "global kernel B2^6",
    )

    # Sharp hostile: every lambda*B2^6 passes all order <=5 data.
    for scalar in (Q(-5, 7), Q(0), Q(11, 3)):
        hostile = [scalar * value for value in b6]
        require(
            matrix_times_vector(matrix, hostile) == [Q(0)] * 18,
            "lambda hostile",
        )
        require(hostile[6] == scalar, "sixth-order lambda lock")
    return rank


def main():
    discriminant = check_slope_cubic()
    rational_controls, split_controls = check_hermite_family()
    jet_controls = check_jet_reconstruction()
    global_rank = check_global_residual_kernel()

    print("THM-2389 H4 THREE-POLE HERMITE-JET EXACT COMPANION")
    print(
        "slope cubic: "
        f"disc={discriminant}; L-shift=5/243; "
        "hostile pivots f(-8/243)=-64/531441, f(4/243)=64/531441"
    )
    print(
        "Hermite numerator: "
        f"rank=6; kernel=B2^2; rational controls={rational_controls}; "
        f"split F59 controls={split_controls}; roots=(5,11,43)"
    )
    print(
        "pole reconstruction: "
        f"orders 0..6 verified on {jet_controls} exact rational controls; "
        "all four parameter pivots nonzero"
    )
    print(
        "global residual: "
        f"conditions=18; rank={global_rank}; kernel=B2^6; "
        "one sixth-order lambda lock is sharp"
    )
    print(
        "VERDICT: three synchronized order-five pole jets leave exactly "
        "lambda*B2^6; one order-six equation is the full spectral identity"
    )


if __name__ == "__main__":
    main()
