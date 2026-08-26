#!/usr/bin/env python3
"""Independent exact audit for THM-4139.

This replay imports no primary-audit code.  It uses first-return tests rather
than a functional-graph census, a direct residue-tube linearization at 43,
and elementary polynomial arithmetic for the squaring-map dynatomic factor.
"""

from collections import Counter
from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def matrix_mul_2(left, right):
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def qmap(x: Fraction) -> Fraction:
    return x * x - Fraction(29, 16)


def tmap(u: int, modulus: int) -> int:
    return ((u * u - 29) * pow(4, -1, modulus)) % modulus


def first_return(value: int, modulus: int, limit: int, map_kind: str) -> int | None:
    current = value
    for step in range(1, limit + 1):
        if map_kind == "quadratic":
            current = tmap(current, modulus)
        elif map_kind == "doubling":
            current = 2 * current % modulus
        else:
            raise RuntimeError(f"unknown map kind: {map_kind}")
        if current == value:
            return step
    return None


def rational_and_ap_audit() -> None:
    numerators = (-7, -5, -3, -1, 1, 3, 5, 7)
    expected_images = (5, -1, -5, -7, -7, -5, -1, 5)
    actual_images = []
    for numerator in numerators:
        image = qmap(Fraction(numerator, 4))
        require(image.denominator == 4, "quarter-integer graph left denominator 4")
        actual_images.append(image.numerator)
    require(tuple(actual_images) == expected_images, "rational graph mismatch")

    # Direct AP solution: orient the cycle L=m-d -> H=m+d -> M=m -> L.
    m = Fraction(-1, 4)
    d = Fraction(3, 2)
    c = Fraction(-29, 16)
    points = (m - d, m + d, m)
    require(
        tuple(x * x + c for x in points) == (points[1], points[2], points[0]),
        "direct AP cycle solution fails",
    )
    require(-4 * m * d == d, "first AP difference equation fails")
    require(d * (2 * m + d) == d, "second AP difference equation fails")
    require(m - d - m * m == c, "AP parameter recovery fails")

    # The three marked parameter representatives and the SL_2 lift.
    def R(t: Fraction) -> Fraction:
        return -1 / (t + 1)

    orbit = [Fraction(1)]
    orbit.extend((R(orbit[-1]), R(R(orbit[-1]))))
    require(orbit == [Fraction(1), Fraction(-1, 2), Fraction(-2)], "R-orbit fails")

    ordered_cycle = [Fraction(-7, 4), Fraction(5, 4), Fraction(-1, 4)]

    def M(x: Fraction) -> Fraction:
        return (4 * x - 13) / (16 * x + 12)

    require(
        [M(x) for x in ordered_cycle] == ordered_cycle[1:] + ordered_cycle[:1],
        "unique point-coordinate interpolant fails",
    )
    B = ((Fraction(1, 4), Fraction(-13, 16)), (1, Fraction(3, 4)))
    require(B[0][0] * B[1][1] - B[0][1] * B[1][0] == 1, "det B != 1")
    require(
        matrix_mul_2(matrix_mul_2(B, B), B) == ((-1, 0), (0, -1)),
        "independent B^3=-I audit fails",
    )

    print("RATIONAL/AP: eight quarter-integers; sole cycle length 3; Q-period 6 absent")
    print("AP SOLUTION: m=-1/4, d=3/2, c=-29/16; marked t-orbit=(1,-1/2,-2)")


def horizontal_fibre_audit() -> None:
    a = Fraction(-48)
    x_cycle = [Fraction(-7, 4), Fraction(5, 4), Fraction(-1, 4)]
    parameters = [4 * x + 1 for x in x_cycle]
    require(parameters == [-6, 6, 0], "independent AP normalization fails")

    points = [
        (r * r + a, r * (r * r + 3 * a / 2))
        for r in parameters
    ]
    require(points == [(-12, 216), (-12, -216), (-48, 0)], "zero-fibre points fail")
    for U, V in points:
        residual = V * V - U**3 + 3 * a * a * U / 4 + a**3 / 4
        require(residual == 0, "point misses q=0 target fibre")

    def N(r: Fraction) -> Fraction:
        return 2 * (r - 6) / (r + 2)

    require(
        [N(r) for r in parameters] == parameters[1:] + parameters[:1],
        "normalized projective interpolant fails",
    )
    C = ((Fraction(1, 2), -3), (Fraction(1, 4), Fraction(1, 2)))
    require(
        matrix_mul_2(matrix_mul_2(C, C), C) == ((-1, 0), (0, -1)),
        "independent C^3=-I audit fails",
    )
    # Order-three N can preserve a two-point set only pointwise.  Its fixed
    # equation is r^2=-12, while the nodal conductor pair has r^2=-3a/2=72.
    require(Fraction(-12) != -3 * a / 2, "N unexpectedly descends through the node")

    rho_squared = -a**3 / 2
    section_u = a / 2 + 16 * rho_squared / (9 * a * a)
    section_v_over_rho = -1 - 64 * rho_squared / (27 * a**3)
    require((section_u, section_v_over_rho) == (Fraction(56, 3), Fraction(5, 27)), "3P q=0 specialization fails")
    require(all(section_u != U for U, _ in points), "3P section meets normalized cycle points")

    print("ZERO FIBRE: r=(-6,6,0) -> ((-12,216),(-12,-216),(-48,0)) on a=-48,q=0")
    print("CARRIER SEPARATION: 3P has U=56/3 at q=0; proposed section identification is false")


def finite_quadratic_audit() -> None:
    return_counts_63 = Counter(
        first_return(value, 63, 6, "quadratic") for value in range(63)
    )
    require(
        return_counts_63 == Counter({None: 54, 3: 9}),
        "first-return census for the quadratic map mod 63 fails",
    )

    p = 43
    modulus = p * p
    base_cycle = (-7, 5, -1)
    tube = {
        (base + p * digit) % modulus
        for base in base_cycle
        for digit in range(p)
    }
    require(len(tube) == 3 * p, "43^2 residue disks overlap")
    return_counts_tube = Counter(
        first_return(value, modulus, 6, "quadratic") for value in tube
    )
    require(
        return_counts_tube == Counter({6: 126, 3: 3}),
        "first-return census in the 43^2 tube fails",
    )
    return_counts_full = Counter(
        first_return(value, modulus, 6, "quadratic") for value in range(modulus)
    )
    require(
        return_counts_full == Counter({None: 1713, 6: 126, 5: 5, 3: 3, 2: 2}),
        "full first-return census modulo 43^2 fails",
    )

    inv2 = pow(2, -1, p)
    for index, base in enumerate(base_cycle):
        target = base_cycle[(index + 1) % 3]
        for digit in range(p):
            source = (base + p * digit) % modulus
            predicted = (target + p * ((base * inv2 * digit) % p)) % modulus
            require(tmap(source, modulus) == predicted, "one-step 43-adic linearization fails")

            after_three = source
            for _ in range(3):
                after_three = tmap(after_three, modulus)
            reflected = (base - p * digit) % modulus
            require(after_three == reflected, "T^3 does not negate the lift digit")

    print("FINITE QUADRATIC: mod63 point returns={3:9}; no six")
    print("FINITE QUADRATIC: 43^2 has exactly 126 period-6 points, all in the tube -> 21 cycles")


def trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_mul(left: list[int], right: list[int]) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            product[i + j] += a * b
    return trim(product)


def poly_div_exact(dividend: list[int], divisor: list[int]) -> list[int]:
    remainder = dividend[:]
    quotient = [0] * max(1, len(dividend) - len(divisor) + 1)
    require(divisor[-1] == 1, "independent polynomial division expects monic divisor")
    while len(remainder) >= len(divisor) and remainder != [0]:
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
        trim(remainder)
    require(remainder == [0], "polynomial division left a remainder")
    return trim(quotient)


def x_power_minus_one(n: int) -> list[int]:
    poly = [-1] + [0] * n
    poly[n] = 1
    return poly


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def cyclotomic(n: int, cache: dict[int, list[int]]) -> list[int]:
    if n not in cache:
        value = x_power_minus_one(n)
        for divisor in divisors(n):
            if divisor < n:
                value = poly_div_exact(value, cyclotomic(divisor, cache))
        cache[n] = value
    return cache[n]


def multiplicative_order_2(modulus: int) -> int:
    require(gcd(2, modulus) == 1, "order requested at an even modulus")
    value = 1
    for order in range(1, modulus + 1):
        value = 2 * value % modulus
        if value == 1:
            return order
    raise RuntimeError("order search failed")


def mersenne_squaring_audit() -> None:
    require(2**6 - 1 == 63 == 3**2 * 7, "63 factorization fails")
    require((2**2 - 1) % 3 == 0 and (2**3 - 1) % 7 == 0, "old-prime proof fails")

    point_returns = Counter(
        first_return(exponent, 63, 6, "doubling") for exponent in range(63)
    )
    require(
        point_returns == Counter({6: 54, 3: 6, 2: 2, 1: 1}),
        "doubling point-return census fails",
    )

    additive_order_counts: Counter[tuple[int, int]] = Counter()
    for exponent in range(63):
        conductor = 63 // gcd(exponent, 63)
        orbit_length = 1 if conductor == 1 else multiplicative_order_2(conductor)
        additive_order_counts[(conductor, orbit_length)] += 1
    require(
        additive_order_counts[(9, 6)] == 6
        and additive_order_counts[(21, 6)] == 12
        and additive_order_counts[(63, 6)] == 36,
        "conductor decomposition of the 54 exact-six points fails",
    )

    cache: dict[int, list[int]] = {1: [-1, 1]}
    exact_factor = [1]
    for conductor in (9, 21, 63):
        exact_factor = poly_mul(exact_factor, cyclotomic(conductor, cache))
    cross_left = poly_mul(x_power_minus_one(63), x_power_minus_one(1))
    cross_right = poly_mul(
        poly_mul(x_power_minus_one(7), x_power_minus_one(3)), exact_factor
    )
    require(cross_left == cross_right, "dynatomic/cyclotomic factor identity fails")
    require(len(exact_factor) - 1 == 54, "exact squaring factor has wrong degree")

    print("SQUARING: point returns={1:1,2:2,3:6,6:54} -> nine six-cycles")
    print("SQUARING: exact factor=Phi_9 Phi_21 Phi_63, degree 54; conductors 9/21/63")


if __name__ == "__main__":
    rational_and_ap_audit()
    horizontal_fibre_audit()
    finite_quadratic_audit()
    mersenne_squaring_audit()
    print("INDEPENDENT_AUDIT=PASS")
