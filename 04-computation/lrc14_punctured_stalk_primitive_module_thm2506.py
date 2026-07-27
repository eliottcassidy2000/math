#!/usr/bin/env python3
"""Exact primitive-module companion for THM-2506.

The script reconstructs a punctured 91-stalk cover, reduces all 72 primitive
Fourier transforms exactly modulo Phi_91, checks the explicit factorization,
verifies convolution by the unit-mask collision kernel, and exhausts every
affine pushforward to Z/13 and Z/169.  It uses only integer/Fraction arithmetic.
"""

from fractions import Fraction as F
from itertools import product


P = 13
Q = 7
N = P * Q
CENTRES = (84, 71, 33, 47, 58)
BLOCKER_SOURCE = 5


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def divisors(integer: int) -> list[int]:
    return [candidate for candidate in range(1, integer + 1) if integer % candidate == 0]


def poly_div_exact(numerator: list[int], denominator: list[int]) -> list[int]:
    """Exact division of low-to-high monic integer polynomials."""

    remainder = numerator[:]
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    while len(remainder) >= len(denominator):
        coefficient = remainder[-1]
        shift = len(remainder) - len(denominator)
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            remainder[index + shift] -= coefficient * value
        while remainder and remainder[-1] == 0:
            remainder.pop()
    require(not remainder, "exact cyclotomic division")
    return quotient


def cyclotomic_polynomial(integer: int, cache: dict[int, list[int]]) -> list[int]:
    if integer in cache:
        return cache[integer]
    polynomial = [-1] + [0] * (integer - 1) + [1]
    for divisor in divisors(integer)[:-1]:
        polynomial = poly_div_exact(polynomial, cyclotomic_polynomial(divisor, cache))
    cache[integer] = polynomial
    return polynomial


PHI91 = cyclotomic_polynomial(N, {1: [-1, 1]})
DEGREE = len(PHI91) - 1


def reduce_zeta91(raw: list[int | F]) -> tuple[F, ...]:
    """Reduce a polynomial of degree below 91 modulo Phi_91."""

    coefficients = [F(value) for value in raw] + [F(0)] * (N - len(raw))
    for exponent in range(N - 1, DEGREE - 1, -1):
        leading = coefficients[exponent]
        if leading:
            shift = exponent - DEGREE
            for index, value in enumerate(PHI91):
                coefficients[index + shift] -= leading * value
    require(all(value == 0 for value in coefficients[DEGREE:]), "cyclotomic remainder degree")
    return tuple(coefficients[:DEGREE])


def add_cyclic(left: dict[int, int], right: dict[int, int]) -> dict[int, int]:
    answer = left.copy()
    for exponent, value in right.items():
        answer[exponent % N] = answer.get(exponent % N, 0) + value
    return {exponent: value for exponent, value in answer.items() if value}


def multiply_cyclic(left: dict[int, int], right: dict[int, int]) -> dict[int, int]:
    answer: dict[int, int] = {}
    for left_exponent, left_value in left.items():
        for right_exponent, right_value in right.items():
            exponent = (left_exponent + right_exponent) % N
            answer[exponent] = answer.get(exponent, 0) + left_value * right_value
    return {exponent: value for exponent, value in answer.items() if value}


def dictionary_remainder(polynomial: dict[int, int]) -> tuple[F, ...]:
    raw = [0] * N
    for exponent, value in polynomial.items():
        raw[exponent % N] += value
    return reduce_zeta91(raw)


def crt_site(row: int, column: int) -> int:
    return next(site for site in range(N) if site % P == row and site % Q == column)


def centred_ap(centre: int, step: int = 1) -> set[int]:
    return {(centre + (offset - 6) * step) % N for offset in range(P)}


def stalk_defect() -> tuple[list[list[int]], set[int]]:
    guard = set(range(26))
    ordinary = [centred_ap(centre) for centre in CENTRES]
    blocker = {site for site in range(N) if site % Q == BLOCKER_SOURCE}
    cover = guard | blocker
    for progression in ordinary:
        cover |= progression
    require(cover == set(range(N)), "explicit guard/ordinary/blocker cover")

    defect = [[0 for _ in range(Q)] for _ in range(P)]
    for row, column in product(range(P), range(Q)):
        site = crt_site(row, column)
        multiplicity = int(site in guard) + sum(int(site in progression) for progression in ordinary)
        defect[row][column] = 1 - multiplicity
    require(all(sum(row) == 0 for row in defect), "row-sum law")
    require(
        {
            row: tuple(value for value in defect[row])
            for row in range(P)
            if any(defect[row])
        }
        == {
            0: (0, 0, 0, -1, 0, 1, 0),
            1: (0, 0, 0, 0, -1, 1, 0),
        },
        "two-row defect profile",
    )
    return defect, cover


def transform(defect: list[list[int]], alpha: int, beta: int) -> tuple[F, ...]:
    raw = [0] * N
    for row, column in product(range(P), range(Q)):
        exponent = (-Q * alpha * row - P * beta * column) % N
        raw[exponent] += defect[row][column]
    return reduce_zeta91(raw)


def factor_remainder(alpha: int, beta: int) -> tuple[F, ...]:
    """xi^3 (xi-1) (xi+1+omega xi), omega=zeta13^-a, xi=zeta7^-b."""

    omega = (-Q * alpha) % N
    xi = (-P * beta) % N
    factor = {3 * xi % N: 1}
    factor = multiply_cyclic(factor, {xi: 1, 0: -1})
    factor = multiply_cyclic(factor, {xi: 1, 0: 1, (omega + xi) % N: 1})
    return dictionary_remainder(factor)


def cyclic_convolution(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    answer = [[0 for _ in range(Q)] for _ in range(P)]
    for h, r, a, b in product(range(P), range(Q), range(P), range(Q)):
        answer[h][r] += left[a][b] * right[(h - a) % P][(r - b) % Q]
    return answer


def affine_pushforward(defect: list[list[int]], modulus: int, slope: int, offset: int) -> list[int]:
    """Every hom C13 x C7 -> Z/13^a kills the C7 coordinate."""

    require(modulus in (P, P * P), "tested 13-primary targets")
    scale = modulus // P
    answer = [0] * modulus
    for row, column in product(range(P), range(Q)):
        target = (scale * slope * row + offset) % modulus
        answer[target] += defect[row][column]
    return answer


def main() -> None:
    defect, cover = stalk_defect()
    require(DEGREE == 72, "Phi_91 degree")
    require(sum(abs(value) for row in defect for value in row) == 4, "explicit l1 norm")

    primitive = {}
    for alpha, beta in product(range(1, P), range(1, Q)):
        value = transform(defect, alpha, beta)
        require(any(value), f"primitive colour {(alpha, beta)}")
        require(value == factor_remainder(alpha, beta), f"factorization {(alpha, beta)}")
        primitive[(alpha, beta)] = value
    require(len(primitive) == 72, "primitive colour census")

    # The unit mask is the sharp squarefree THM-2474 collision kernel. Its
    # primitive transform is mu(91)=1, so convolution preserves every value.
    unit_mask = [
        [int(row != 0 and column != 0) for column in range(Q)]
        for row in range(P)
    ]
    convolved = cyclic_convolution(defect, unit_mask)
    for alpha, beta in primitive:
        require(
            transform(convolved, alpha, beta) == primitive[(alpha, beta)],
            "unit-mask primitive convolution automorphism",
        )

    pushforward_count = 0
    for modulus in (P, P * P):
        for slope in range(P):
            for offset in range(modulus):
                require(
                    not any(affine_pushforward(defect, modulus, slope, offset)),
                    f"affine pushforward to Z/{modulus}",
                )
                pushforward_count += 1
    require(pushforward_count == 2366, "affine pushforward census")

    floor_one = F(2, 7)
    floor_two = F(3, 7)
    require(floor_one / 72 == F(1, 252), "old k=1 union floor")
    require(floor_two / 72 == F(1, 168), "old k=2 union floor")

    print("THM-2506 PUNCTURED-STALK PRIMITIVE-MODULE EXACT COMPANION")
    print(f"cover_size={len(cover)};centres=" + ",".join(map(str, CENTRES)) + ";blocker_source=5")
    print("nonzero_defect_rows=0:e5-e3,1:e5-e4;l1=4;inherited_universal_l1_cap=18")
    print(f"phi91_degree={DEGREE};primitive_colours={len(primitive)};factorization=PASS")
    print("fixed_mode_floors=2/7,3/7;old_union_floors=1/252,1/168")
    print("unit_mask_convolution_primitive_colours=72;multiplier=mu(91)=1")
    print(f"affine_13_primary_pushforwards={pushforward_count};nonzero_pushforwards=0")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
