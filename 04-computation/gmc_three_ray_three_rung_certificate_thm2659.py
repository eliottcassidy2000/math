#!/usr/bin/env python3
"""Exact referee for THM-2659's free three-ray three-rung certificate."""

from itertools import product
from math import comb, factorial, gcd


CHARGES = (24, 9, 4, -6)
RAYS = (
    (1, 0, 0, 4),
    (0, 2, 0, 3),
    (0, 0, 3, 2),
)
ELL = 5
W_EXPONENTS = (0, 15, 20, 30)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first,) + tail


def channel_weight(channel):
    total = sum(channel)
    value = factorial(total)
    for entry in channel:
        value //= factorial(entry)
    return value


def determinant3(matrix):
    require(len(matrix) == 3 and all(len(row) == 3 for row in matrix),
            "3 by 3 determinant shape")
    a, b, c = matrix
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def add_channels(coefficients):
    return tuple(
        sum(coefficients[index] * RAYS[index][coordinate]
            for index in range(3))
        for coordinate in range(4)
    )


def rung(level):
    result = {}
    for exponents in compositions(level, 3):
        channel = add_channels(exponents)
        result[exponents] = channel_weight(channel)
    return result


def substitute_x(poly, a=-2, b=-2):
    """Substitute x=a*y+b*z in a homogeneous Q[x,y,z] dictionary."""
    result = {}
    for (x_degree, y_degree, z_degree), coefficient in poly.items():
        for y_from_x in range(x_degree + 1):
            exponent = (y_degree + y_from_x,
                        z_degree + x_degree - y_from_x)
            term = (coefficient * comb(x_degree, y_from_x)
                    * a ** y_from_x * b ** (x_degree - y_from_x))
            result[exponent] = result.get(exponent, 0) + term
    return {key: value for key, value in result.items() if value}


def poly_add(left, right):
    result = dict(left)
    for degree, coefficient in right.items():
        result[degree] = result.get(degree, 0) + coefficient
        if result[degree] == 0:
            del result[degree]
    return result


def poly_scale(poly, scalar):
    return {degree: scalar * coefficient for degree, coefficient in poly.items()
            if scalar * coefficient}


def poly_multiply(left, right):
    result = {}
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            result[degree] = (result.get(degree, 0)
                              + left_coefficient * right_coefficient)
    return {degree: coefficient for degree, coefficient in result.items()
            if coefficient}


def bareiss_determinant(matrix):
    work = [row[:] for row in matrix]
    size = len(work)
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next((row for row in range(pivot_index + 1, size)
                         if work[row][pivot_index] != 0), None)
            require(swap is not None, "singular Sylvester pivot")
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (work[row][column] * pivot
                             - work[row][pivot_index] * work[pivot_index][column])
                require(numerator % previous == 0, "nonintegral Bareiss division")
                work[row][column] = numerator // previous
        previous = pivot
        for row in range(pivot_index + 1, size):
            work[row][pivot_index] = 0
    return sign * work[-1][-1]


def resultant_quadratic_cubic(q_coefficients, c_coefficients):
    require(len(q_coefficients) == 3 and len(c_coefficients) == 4,
            "resultant degree")
    matrix = []
    for shift in range(3):
        row = [0] * 5
        row[shift:shift + 3] = q_coefficients
        matrix.append(row)
    for shift in range(2):
        row = [0] * 5
        row[shift:shift + 4] = c_coefficients
        matrix.append(row)
    return bareiss_determinant(matrix)


def main():
    require(all(sum(charge * entry for charge, entry in zip(CHARGES, ray)) == 0
                for ray in RAYS), "ray balance")
    require(all(sum(ray) == ELL for ray in RAYS), "equal ray mass")
    require(all(gcd(*ray) == 1 for ray in RAYS), "ray primitivity")

    maximal_minors = tuple(
        abs(determinant3([
            [RAYS[column][coordinate] for column in range(3)]
            for coordinate in range(4) if coordinate != omitted
        ]))
        for omitted in range(4)
    )
    require(tuple(sorted(maximal_minors)) == (4, 6, 9, 24),
            "ray maximal minors")
    require(gcd(*maximal_minors) == 1, "balanced lattice index")

    channel_checks = 0
    hilbert = []
    for mass in range(31):
        actual = set()
        for prefix in product(range(mass + 1), repeat=3):
            fourth = mass - sum(prefix)
            if fourth < 0:
                continue
            channel = prefix + (fourth,)
            if sum(charge * entry for charge, entry in zip(CHARGES, channel)) == 0:
                actual.add(channel)
        if mass % ELL:
            expected = set()
        else:
            level = mass // ELL
            expected = {add_channels(exponents)
                        for exponents in compositions(level, 3)}
        require(actual == expected, f"channel semigroup at mass {mass}")
        expected_count = 0 if mass % ELL else comb(mass // ELL + 2, 2)
        require(len(actual) == expected_count, f"Hilbert count at mass {mass}")
        channel_checks += len(actual)
        if actual:
            hilbert.append((mass, len(actual)))

    t5 = rung(1)
    t10 = rung(2)
    t15 = rung(3)
    expected_t5 = {(1, 0, 0): 5, (0, 1, 0): 10, (0, 0, 1): 10}
    expected_t10 = {
        (2, 0, 0): 45, (1, 1, 0): 360, (1, 0, 1): 840,
        (0, 2, 0): 210, (0, 1, 1): 2520, (0, 0, 2): 210,
    }
    expected_t15 = {
        (3, 0, 0): 455, (2, 1, 0): 8190, (2, 0, 1): 30030,
        (1, 2, 0): 15015, (1, 1, 1): 300300, (1, 0, 2): 45045,
        (0, 3, 0): 5005, (0, 2, 1): 225225,
        (0, 1, 2): 180180, (0, 0, 3): 5005,
    }
    require(t5 == expected_t5 and t10 == expected_t10 and t15 == expected_t15,
            "first three rung coefficients")

    q_homogeneous = {(2, 0): -330, (1, 1): 480, (0, 2): -1290}
    c_homogeneous = {
        (3, 0): 4095, (2, 1): -230685,
        (1, 2): -248430, (0, 3): 31395,
    }
    require(substitute_x(t10) == q_homogeneous, "quadratic reduction")
    require(substitute_x(t15) == c_homogeneous, "cubic reduction")

    q = {2: 11, 1: -16, 0: 43}
    c = {3: 3, 2: -169, 1: -182, 0: 23}
    u = {2: -157251, 1: 8616505, 0: 23433513}
    v = {1: 576587, 0: 48544}
    bezout = poly_add(poly_multiply(u, q), poly_multiply(v, c))
    require(bezout == {0: 1008757571}, "Bezout identity")
    resultant = resultant_quadratic_cubic((11, -16, 43), (3, -169, -182, 23))
    require(resultant == 1008757571, "quadratic-cubic resultant")
    q_discriminant = (-16) ** 2 - 4 * 11 * 43
    q_at_zero = 43
    q_at_minus_one = 11 + 16 + 43
    require(q_discriminant == -1636 and q_at_zero != 0 and q_at_minus_one != 0,
            "two-rung torus hostile")

    wick_checks = 0
    for mass in range(16):
        for channel in compositions(mass, 4):
            charge = sum(value * exponent for value, exponent in zip(CHARGES, channel))
            z_degree = 24 * mass
            w_degree = sum(value * exponent
                           for value, exponent in zip(W_EXPONENTS, channel))
            require((z_degree == w_degree) == (charge == 0),
                    "horizontal Wick balance")
            wick_checks += 1

    require(q_homogeneous[(2, 0)] != 0, "z=0 projective boundary")
    print("THM2659 GMC FREE THREE-RAY THREE-RUNG EXACT REFEREE")
    print("charges=(24,9,4,-6) rays=((1,0,0,4),(0,2,0,3),(0,0,3,2)) mass=5 index=1")
    print(f"channel_checks={channel_checks} hilbert_through_30={tuple(hilbert)}")
    print("T5=5*(x+2y+2z)")
    print("T10=15*(3x^2+24xy+56xz+14y^2+168yz+14z^2)")
    print("T15=455*(x^3+18x^2y+66x^2z+33xy^2+660xyz+99xz^2+11y^3+495y^2z+396yz^2+11z^3)")
    print("reduced_Q=11t^2-16t+43 reduced_C=3t^3-169t^2-182t+23")
    print(f"resultant_Q_C={resultant} bezout_rhs={bezout[0]}")
    print("two_rung_torus_hostile=Q_disc:-1636,Q(0):43,Q(-1):70")
    print(f"horizontal_support_W=(0,15,20,30) wick_checks={wick_checks} common_factor=(24m)!")
    print("PASS")


if __name__ == "__main__":
    main()
