#!/usr/bin/env python3
"""Exact companion for THM-2860's Euler-tangent cubic-multipole census.

The companion has two jobs.

1. Verify the intrinsic algebra behind the Euler-transversality invariant

       Delta_A(E) = det[c, conjugate(c), A c, A conjugate(c)]

   and the Segre/Mobius description of its zero locus.
2. For every support 0 <= a0 < a1 < a2 < a3 <= 30, parameterize every
   Euler-tangent real chord by E=span(T, A T), construct the two exact
   degree-seven cubic-divisibility remainders, and prove that they are
   coprime by reduction modulo the good prime 1,000,033.

The finite census excludes Euler-tangent *cubic* Maxwell lines in the stated
box.  It says nothing about quartic-harmonic vanishing, shared transverse
lines, or scalar access to the Euler/Laguerre multiplier.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import combinations, product
from math import comb

import sympy as sp


PRIME = 1_000_033
MAX_EXPONENT = 30
MAX_MOMENT_DEGREE = 3 * MAX_EXPONENT


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------------------
# Small polynomial arithmetic over F_p, with coefficients in ascending order.
# ---------------------------------------------------------------------------


def trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] % PRIME == 0:
        poly.pop()
    return [coefficient % PRIME for coefficient in poly]


def add(first: list[int], second: list[int]) -> list[int]:
    answer = [0] * max(len(first), len(second))
    for index, coefficient in enumerate(first):
        answer[index] += coefficient
    for index, coefficient in enumerate(second):
        answer[index] += coefficient
    return trim(answer)


def scale(poly: list[int], scalar: int) -> list[int]:
    return trim([(scalar * coefficient) % PRIME for coefficient in poly])


def multiply(first: list[int], second: list[int]) -> list[int]:
    answer = [0] * (len(first) + len(second) - 1)
    for first_index, first_coefficient in enumerate(first):
        for second_index, second_coefficient in enumerate(second):
            answer[first_index + second_index] = (
                answer[first_index + second_index]
                + first_coefficient * second_coefficient
            ) % PRIME
    return trim(answer)


def polynomial_remainder(dividend: list[int], divisor: list[int]) -> list[int]:
    require(divisor != [0], "polynomial division by zero")
    remainder = trim(dividend[:])
    divisor = trim(divisor[:])
    inverse_lead = pow(divisor[-1], PRIME - 2, PRIME)
    while len(remainder) >= len(divisor) and remainder != [0]:
        shift = len(remainder) - len(divisor)
        quotient_term = remainder[-1] * inverse_lead % PRIME
        for index, coefficient in enumerate(divisor):
            remainder[shift + index] = (
                remainder[shift + index] - quotient_term * coefficient
            ) % PRIME
        remainder = trim(remainder)
    return remainder


def polynomial_gcd_degree(first: list[int], second: list[int]) -> int:
    first = trim(first[:])
    second = trim(second[:])
    while second != [0]:
        first, second = second, polynomial_remainder(first, second)
    return len(first) - 1


def has_full_septic_degrees(first: list[int], second: list[int]) -> bool:
    """Whether both affine forms retain their projective degree seven."""

    return len(trim(first[:])) - 1 == 7 and len(trim(second[:])) - 1 == 7


# Factorials are invertible because every degree used below is less than PRIME.
FACTORIALS = [1] * (MAX_MOMENT_DEGREE + 1)
for factorial_index in range(1, len(FACTORIALS)):
    FACTORIALS[factorial_index] = (
        FACTORIALS[factorial_index - 1] * factorial_index
    ) % PRIME

INVERSE_FACTORIALS = [
    pow(FACTORIALS[index], PRIME - 2, PRIME)
    for index in range(MAX_EXPONENT + 1)
]


def factorial_tensor_entry(
    support: tuple[int, int, int, int],
    indices: tuple[int, ...],
) -> int:
    """L(product f_(a_i)) modulo PRIME for normalized monomials f_a."""

    coefficient = FACTORIALS[sum(support[index] for index in indices)]
    for index in indices:
        coefficient = coefficient * INVERSE_FACTORIALS[support[index]] % PRIME
    return coefficient


def contract(
    support: tuple[int, int, int, int],
    vectors: tuple[tuple[list[int], ...], ...],
) -> list[int]:
    """Contract the factorial moment tensor against polynomial vectors."""

    answer = [0]
    order = len(vectors)
    for indices in product(range(4), repeat=order):
        term = [1]
        for vector_index, coordinate_index in enumerate(indices):
            term = multiply(term, vectors[vector_index][coordinate_index])
        answer = add(
            answer,
            scale(term, factorial_tensor_entry(support, indices)),
        )
    return answer


def tangent_pencil(
    support: tuple[int, int, int, int],
) -> tuple[tuple[list[int], ...], tuple[list[int], ...]]:
    """Return T(z) and A T(z), with T in ker(sum) intersect ker(sum a_j *)."""

    a0, a1, a2, a3 = support
    first_gap = a1 - a0
    tangent = (
        [a2 - a1, a3 - a1],
        [-(a2 - a0), -(a3 - a0)],
        [first_gap, 0],
        [0, first_gap],
    )
    euler_tangent = tuple(
        scale(tangent[index], support[index]) for index in range(4)
    )
    return tangent, euler_tangent


def tangent_remainders(
    support: tuple[int, int, int, int],
) -> tuple[list[int], list[int]]:
    """The two division-free conditions for Q to divide C on span(T,AT)."""

    tangent, euler_tangent = tangent_pencil(support)
    g11 = contract(support, (tangent, tangent))
    g12 = contract(support, (tangent, euler_tangent))
    g22 = contract(support, (euler_tangent, euler_tangent))
    t111 = contract(support, (tangent, tangent, tangent))
    t112 = contract(support, (tangent, tangent, euler_tangent))
    t122 = contract(support, (tangent, euler_tangent, euler_tangent))
    t222 = contract(
        support,
        (euler_tangent, euler_tangent, euler_tangent),
    )

    first = add(
        add(
            scale(multiply(multiply(t112, g11), g22), 3),
            scale(multiply(t222, multiply(g11, g11)), -1),
        ),
        scale(multiply(multiply(t111, g12), g22), -2),
    )
    second = add(
        add(
            scale(multiply(multiply(t122, g11), g22), 3),
            scale(multiply(multiply(t222, g12), g11), -2),
        ),
        scale(multiply(t111, multiply(g22, g22)), -1),
    )
    return first, second


# ---------------------------------------------------------------------------
# Independent exact constructor over Q for selected hostile/positive controls.
# ---------------------------------------------------------------------------


def rational_factorial_readout(poly: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    expanded = sp.Poly(sp.expand(poly), variable)
    return sp.expand(
        sum(
            coefficient * sp.factorial(exponent[0])
            for exponent, coefficient in expanded.terms()
        )
    )


def rational_tangent_remainders(
    support: tuple[int, int, int, int],
) -> tuple[sp.Poly, sp.Poly]:
    variable, parameter = sp.symbols("s parameter")
    a0, a1, a2, a3 = support
    first_gap = a1 - a0
    coefficients = (
        (a2 - a1) + (a3 - a1) * parameter,
        -(a2 - a0) - (a3 - a0) * parameter,
        sp.Integer(first_gap),
        sp.Integer(first_gap) * parameter,
    )
    tangent = sum(
        coefficients[index]
        * variable ** support[index]
        / sp.factorial(support[index])
        for index in range(4)
    )
    euler_tangent = sum(
        support[index]
        * coefficients[index]
        * variable ** support[index]
        / sp.factorial(support[index])
        for index in range(4)
    )

    def readout(*factors: sp.Expr) -> sp.Expr:
        return rational_factorial_readout(sp.prod(factors), variable)

    g11 = readout(tangent, tangent)
    g12 = readout(tangent, euler_tangent)
    g22 = readout(euler_tangent, euler_tangent)
    t111 = readout(tangent, tangent, tangent)
    t112 = readout(tangent, tangent, euler_tangent)
    t122 = readout(tangent, euler_tangent, euler_tangent)
    t222 = readout(euler_tangent, euler_tangent, euler_tangent)
    first = sp.Poly(
        sp.expand(
            3 * t112 * g11 * g22
            - t222 * g11**2
            - 2 * t111 * g12 * g22
        ),
        parameter,
        domain=sp.QQ,
    )
    second = sp.Poly(
        sp.expand(
            3 * t122 * g11 * g22
            - 2 * t222 * g12 * g11
            - t111 * g22**2
        ),
        parameter,
        domain=sp.QQ,
    )
    return first, second


def reduce_rational_poly(poly: sp.Poly) -> list[int]:
    answer = [0] * (poly.degree() + 1)
    for (degree,), coefficient in poly.terms():
        numerator, denominator = sp.fraction(coefficient)
        answer[degree] = (
            int(numerator) % PRIME
        ) * pow(int(denominator) % PRIME, PRIME - 2, PRIME) % PRIME
    return trim(answer)


def symbolic_geometry_controls() -> int:
    """Verify determinant, gauge, Mobius, and degenerate-face identities."""

    support_symbols = sp.symbols("a0:4", real=True)
    real_parts = sp.symbols("u0:4", real=True)
    imaginary_parts = sp.symbols("v0:4", real=True)
    multiplier_real, multiplier_imag = sp.symbols(
        "multiplier_real multiplier_imag",
        real=True,
    )

    vector = sp.Matrix(
        [
            real_parts[index] + sp.I * imaginary_parts[index]
            for index in range(4)
        ]
    )
    conjugate_vector = sp.conjugate(vector)
    euler = sp.diag(*support_symbols)
    complex_delta = sp.det(
        sp.Matrix.hstack(
            vector,
            conjugate_vector,
            euler * vector,
            euler * conjugate_vector,
        )
    )
    real_delta = -4 * sp.det(
        sp.Matrix.hstack(
            sp.Matrix(real_parts),
            sp.Matrix(imaginary_parts),
            euler * sp.Matrix(real_parts),
            euler * sp.Matrix(imaginary_parts),
        )
    )
    require(
        sp.expand(complex_delta - real_delta) == 0,
        "complex/real Euler determinant identity",
    )

    multiplier = multiplier_real + sp.I * multiplier_imag
    scaled_delta = sp.det(
        sp.Matrix.hstack(
            multiplier * vector,
            sp.conjugate(multiplier) * conjugate_vector,
            multiplier * euler * vector,
            sp.conjugate(multiplier) * euler * conjugate_vector,
        )
    )
    require(
        sp.expand(
            scaled_delta
            - (multiplier_real**2 + multiplier_imag**2) ** 2 * complex_delta
        )
        == 0,
        "Euler determinant gauge law",
    )

    normal_symbols = sp.symbols("b0:4", real=True)
    segre_determinant = sp.det(
        sp.Matrix(
            [
                [1, 1, 1, 1],
                list(support_symbols),
                list(normal_symbols),
                [
                    support_symbols[index] * normal_symbols[index]
                    for index in range(4)
                ],
            ]
        )
    )
    numerator_constant, numerator_linear = sp.symbols(
        "numerator_constant numerator_linear"
    )
    denominator_constant, denominator_linear = sp.symbols(
        "denominator_constant denominator_linear"
    )
    mobius_normal = {
        normal_symbols[index]: -(
            numerator_constant + numerator_linear * support_symbols[index]
        )
        / (
            denominator_constant
            + denominator_linear * support_symbols[index]
        )
        for index in range(4)
    }
    cleared_mobius = sp.together(segre_determinant.subs(mobius_normal))
    require(
        sp.factor(sp.fraction(cleared_mobius)[0]) == 0,
        "Mobius normals lie on the Segre coplanarity conic",
    )

    exceptional_value, common_value = sp.symbols(
        "exceptional_value common_value"
    )
    degenerate_normal = {
        normal_symbols[0]: exceptional_value,
        normal_symbols[1]: common_value,
        normal_symbols[2]: common_value,
        normal_symbols[3]: common_value,
    }
    require(
        sp.expand(segre_determinant.subs(degenerate_normal)) == 0,
        "denominator-degenerate normal",
    )

    parameter = sp.symbols("parameter", real=True)
    a0, a1, a2, a3 = support_symbols
    tangent = sp.Matrix(
        [
            (a2 - a1) + (a3 - a1) * parameter,
            -(a2 - a0) - (a3 - a0) * parameter,
            a1 - a0,
            (a1 - a0) * parameter,
        ]
    )
    require(sp.expand(sum(tangent)) == 0, "tangent pencil mean")
    require(
        sp.expand(sum(support_symbols[index] * tangent[index] for index in range(4)))
        == 0,
        "Euler tangent pencil mean",
    )
    tangent_complex = tangent + sp.I * euler * tangent
    require(
        sp.expand(
            sp.det(
                sp.Matrix.hstack(
                    tangent_complex,
                    sp.conjugate(tangent_complex),
                    euler * tangent_complex,
                    euler * sp.conjugate(tangent_complex),
                )
            )
        )
        == 0,
        "tangent pencil has zero Euler determinant",
    )
    return 7


def main() -> None:
    require(bool(sp.isprime(PRIME)), "certificate modulus is not prime")
    require(
        PRIME > MAX_MOMENT_DEGREE,
        "certificate modulus does not clear the factorial degrees",
    )
    symbolic_controls = symbolic_geometry_controls()

    # Independent Q-to-F_p agreement controls.
    constructor_controls = (
        (0, 1, 2, 3),
        (0, 2, 5, 7),
        (3, 7, 11, 17),
    )
    for support in constructor_controls:
        rational_first, rational_second = rational_tangent_remainders(support)
        modular_first, modular_second = tangent_remainders(support)
        require(
            reduce_rational_poly(rational_first) == modular_first,
            f"independent first remainder constructor {support}",
        )
        require(
            reduce_rational_poly(rational_second) == modular_second,
            f"independent second remainder constructor {support}",
        )

    # Hostile controls: the Euclidean engine must expose a planted factor, and
    # the degree guard must see a planted projective-root-at-infinity defect.
    planted_factor = [1, 1, 1]
    planted_first = multiply(planted_factor, [2, 3, 5, 7, 11, 13])
    planted_second = multiply(planted_factor, [17, 19, 23, 29, 31, 37])
    planted_gcd_degree = polynomial_gcd_degree(planted_first, planted_second)
    require(planted_gcd_degree == 2, "planted common-factor hostile")
    planted_degree_six = [1, 2, 3, 4, 5, 6, 7]
    planted_degree_seven = [1, 2, 3, 4, 5, 6, 7, 8]
    require(
        not has_full_septic_degrees(
            planted_degree_six,
            planted_degree_seven,
        ),
        "degree-drop hostile passed the census degree predicate",
    )
    require(
        has_full_septic_degrees(
            planted_degree_seven,
            planted_degree_seven,
        ),
        "full-degree positive control failed the census degree predicate",
    )

    support_count = 0
    minimum_first_degree = 100
    minimum_second_degree = 100
    maximum_gcd_degree = 0
    first_leading_coefficient = None
    last_leading_coefficient = None
    census_digest = sha256()
    census_digest.update(b"THM-2860/euler-tangent-cubic-census/v1")
    for support in combinations(range(MAX_EXPONENT + 1), 4):
        first, second = tangent_remainders(support)
        first_degree = len(first) - 1
        second_degree = len(second) - 1
        gcd_degree = polynomial_gcd_degree(first, second)
        require(
            has_full_septic_degrees(first, second),
            f"projective septic degree {support}",
        )
        require(gcd_degree == 0, f"tangent cubic common factor {support}")
        for exponent in support:
            census_digest.update(exponent.to_bytes(2, "big"))
        for label, polynomial in ((b"I1", first), (b"I2", second)):
            census_digest.update(label)
            census_digest.update(len(polynomial).to_bytes(2, "big"))
            for coefficient in polynomial:
                census_digest.update(coefficient.to_bytes(4, "big"))
        census_digest.update(gcd_degree.to_bytes(2, "big"))
        minimum_first_degree = min(minimum_first_degree, first_degree)
        minimum_second_degree = min(minimum_second_degree, second_degree)
        maximum_gcd_degree = max(maximum_gcd_degree, gcd_degree)
        if support_count == 0:
            first_leading_coefficient = (first[-1], second[-1])
        last_leading_coefficient = (first[-1], second[-1])
        support_count += 1

    expected_count = comb(MAX_EXPONENT + 1, 4)
    require(support_count == expected_count == 31_465, "support universe count")

    print("THM-2860 FOUR-SLOT EULER-TANGENT CUBIC CENSUS")
    print("status=FINITE-EXACT+SYMBOLIC-GEOMETRY")
    print(f"symbolic_geometry_controls={symbolic_controls}")
    print("delta=det(c,cbar,A*c,A*cbar)=-4*det(u,v,A*u,A*v)")
    print("gauge=delta(lambda*c)=abs(lambda)^4*delta(c)")
    print("delta_zero=plane_contains_nonzero_T_with_A*T_in_plane")
    print("dual_zero_locus=det(rows(1,a,b,a*b))=0")
    print("generic_zero_locus=b_j=-(alpha+beta*a_j)/(gamma+delta*a_j)")
    print("degenerate_zero_locus=three_equal_b_j=coordinate_three-slot_face")
    print(
        "tangent_pencil=("
        "a2-a1+(a3-a1)z,"
        "-(a2-a0)-(a3-a0)z,"
        "a1-a0,"
        "(a1-a0)z)"
    )
    print("remainders=I1,I2 are homogeneous degree7 in tangent coordinates")
    print("leading_terms=I1(T_infinity),I2(T_infinity);both nonzero in census")
    print(f"prime={PRIME} good_above_max_tensor_degree={MAX_MOMENT_DEGREE}")
    print(f"independent_constructor_controls={len(constructor_controls)}")
    print(f"planted_common_factor_degree={planted_gcd_degree}")
    print("planted_degree_drop_rejected=PASS")
    print(f"supports={support_count} expected={expected_count}")
    print(
        "degree_minima="
        f"({minimum_first_degree},{minimum_second_degree}) "
        f"maximum_gcd_degree={maximum_gcd_degree}"
    )
    print(f"first_cell_leading_pair={first_leading_coefficient}")
    print(f"last_cell_leading_pair={last_leading_coefficient}")
    print(f"census_sha256={census_digest.hexdigest()}")
    print("cubic_euler_tangent_cells=0")
    print("quartic_harmonic_zero_excluded=NO")
    print("shared_euler_transverse_line_excluded=NO")
    print("scalar_multiplier_access_obtained=NO")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
