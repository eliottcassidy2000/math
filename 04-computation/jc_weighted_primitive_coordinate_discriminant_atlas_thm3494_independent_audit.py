#!/usr/bin/env python3
"""Independent FLINT audit for THM-3494.

This companion does not import Sympy or the primary companion.  It rebuilds
the n=3 and n=4 inverse equations in ``fmpq_mpoly``, takes FLINT resultants
and discriminants, and checks the branch/index square identities.  It also
reconstructs the n=3 z eliminant, fires the flat-coordinate hostile, and
records the sharp constant-unit and tournament/XOR corrections to the
square-class interpretation.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import isqrt

from flint import fmpq, fmpq_mpoly_ctx


EXPECTED_SEMANTIC_SHA256 = "270c5b6217007b864fd3b99bb548bc55a5d26ff109c30126c1a9553c820a112f"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


R = fmpq_mpoly_ctx.get(["w", "P", "Q", "C", "X", "Y", "Z"])
w, P, Q, C, X, Y, Z = R.gens()


def digest(poly: object, variable_indices: tuple[int, ...]) -> str:
    payload = tuple(
        (
            tuple(monomial[index] for index in variable_indices),
            int(coefficient.p),
            int(coefficient.q),
        )
        for monomial, coefficient in poly.terms()
    )
    return sha256(repr(payload).encode("ascii")).hexdigest()


def term_count(poly: object) -> int:
    return sum(1 for _ in poly.terms())


def exact_quotient(numerator: object, denominator: object, label: str) -> object:
    quotient, remainder = divmod(numerator, denominator)
    require(remainder == 0, (label, "nonzero remainder", remainder))
    return quotient


def audit_coordinate(
    label: str,
    degree: int,
    inverse_equation: object,
    branch_discriminant: object,
    branch_factor: object,
    relation: object,
    coordinate_index: int,
    c_exponent: int,
    expected_core_digest: str,
) -> tuple[object, ...]:
    eliminant = inverse_equation.resultant(relation, 0)
    require(eliminant.degrees()[coordinate_index] == degree, (label, eliminant.degrees()))
    coordinate_discriminant = eliminant.discriminant(coordinate_index)
    require(coordinate_discriminant != 0, (label, "zero discriminant"))
    ratio = exact_quotient(coordinate_discriminant, branch_discriminant, label)
    index = ratio.sqrt()
    require(index * index == ratio, (label, "nonsquare ratio"))
    core = exact_quotient(index, C**c_exponent, label + " C-core")
    require(branch_factor.gcd(core) == 1, (label, "branch/index gcd"))
    core_digest = digest(core, (1, 2))
    require(core_digest == expected_core_digest, (label, core_digest, expected_core_digest))
    return (
        label,
        degree,
        term_count(eliminant),
        c_exponent,
        core_digest,
        str(ratio.factor()),
    )


def degree_row(degree: int) -> tuple[object, ...]:
    if degree == 3:
        seed = -3 * w**2 + 2 * w
        inverse_equation = -w**3 + w**2 - P * w + Q
        expected_branch_scalar = fmpq(-1)
        expected_branch_digest = "37562bfd2b6d06cedcfabcb4b4904db42f503f410ffa86a6534c6b24e6d7697e"
        expected_x = "ce9f8dadb0528fed1a40d9aa8425acc51d23e5b2c59ee0c61919e3bb0d49c454"
        expected_y = "1c9df71a6e89eb7f8040074a60bceceafc7e298505f69789f2c5707774785acb"
    elif degree == 4:
        seed = -w**3 - fmpq(3, 2) * w**2 + fmpq(3, 2) * w
        inverse_equation = (
            -fmpq(1, 4) * w**4
            - fmpq(1, 2) * w**3
            + fmpq(3, 4) * w**2
            - P * w
            + Q
        )
        expected_branch_scalar = fmpq(-1, 16)
        expected_branch_digest = "62c7064c3132534306fa01216b99211a1a1bd5651226758277fa12edab971978"
        expected_x = "e4cf75a510e1bc6ca40b2c692ec80845fa20385918cfc8a25509d314ccbc6b36"
        expected_y = "a0ce0c6d577e9982d2f88ffb7ba1a861d1628e31e0b428b92e731a00488f19c1"
    else:
        raise RuntimeError(("unsupported independent row", degree))

    require(inverse_equation.derivative(0) == seed - P, (degree, "derivative"))
    branch_discriminant = inverse_equation.discriminant(0)
    branch_scalar, factors = branch_discriminant.factor()
    require(branch_scalar == expected_branch_scalar, (degree, branch_scalar))
    require(len(factors) == 1 and factors[0][1] == 1, (degree, factors))
    branch_factor = factors[0][0]
    branch_digest = digest(branch_factor, (1, 2))
    require(branch_digest == expected_branch_digest, (degree, branch_digest))

    gamma = P - seed
    c_exponent = degree * (degree - 1) // 2
    x_row = audit_coordinate(
        "x",
        degree,
        inverse_equation,
        branch_discriminant,
        branch_factor,
        X * gamma - C,
        4,
        c_exponent,
        expected_x,
    )
    y_row = audit_coordinate(
        "y",
        degree,
        inverse_equation,
        branch_discriminant,
        branch_factor,
        C * Y - w + gamma,
        5,
        c_exponent,
        expected_y,
    )

    flat = inverse_equation.resultant(X - P, 0)
    require(flat == (X - P) ** degree, (degree, "flat resultant", flat))
    require(flat.discriminant(4) == 0, (degree, "flat hostile missed"))

    return (
        degree,
        branch_scalar,
        term_count(branch_factor),
        branch_digest,
        x_row,
        y_row,
        "flat=(X-P)^n;disc=0",
    )


def cubic_z_row() -> tuple[object, ...]:
    seed = -3 * w**2 + 2 * w
    inverse_equation = -w**3 + w**2 - P * w + Q
    gamma = P - seed
    a = -fmpq(3, 2)
    z_numerator = gamma * (gamma * (gamma - 1 + a) - a * w)
    branch_discriminant = inverse_equation.discriminant(0)
    _, factors = branch_discriminant.factor()
    row = audit_coordinate(
        "z",
        3,
        inverse_equation,
        branch_discriminant,
        factors[0][0],
        C**2 * Z - z_numerator,
        6,
        6,
        "b45429234347865c0b4c372bccccaea692692c57ae1409290c8eca3205be1a00",
    )
    require(row[2] == 44, ("z eliminant terms", row[2]))
    return row


def is_rational_square(value: Fraction) -> bool:
    require(value >= 0, ("negative square test", value))
    return (
        isqrt(value.numerator) ** 2 == value.numerator
        and isqrt(value.denominator) ** 2 == value.denominator
    )


def xor_hostile() -> tuple[object, ...]:
    # Exact vertex volumes v_1=1, v_2=2, v_3=6 give edge ratios 2,3,6.
    # All are nonsquares, so the naive nonsquare indicator is not additive:
    # 1 XOR 1 != 1.  The full square class is still a vertex coboundary, and
    # every discriminant ratio g_ij^2 is an actual square.
    volumes = (Fraction(1), Fraction(2), Fraction(6))
    g12 = volumes[1] / volumes[0]
    g23 = volumes[2] / volumes[1]
    g13 = volumes[2] / volumes[0]
    naive_bits = tuple(int(not is_rational_square(g)) for g in (g12, g23, g13))
    require(g12 * g23 == g13, ("coboundary law", g12, g23, g13))
    require(naive_bits == (1, 1, 1), naive_bits)
    require((naive_bits[0] ^ naive_bits[1]) != naive_bits[2], naive_bits)
    require(all(is_rational_square(g * g) for g in (g12, g23, g13)), "disc ratios")
    valuation_two_bits = (1, 0, 1)
    require(valuation_two_bits[0] ^ valuation_two_bits[1] == valuation_two_bits[2], valuation_two_bits)
    return (
        "volumes=(1,2,6)",
        "g=(2,3,6)",
        "naive_nonsquare_bits=(1,1,1)_not_XOR",
        "v2_character_bits=(1,0,1)_is_XOR",
        "discriminant_ratios=(4,9,36)_square_zero",
    )


def main() -> None:
    rows = tuple(degree_row(degree) for degree in (3, 4))
    z_row = cubic_z_row()
    hostile = xor_hostile()
    unit_hostile = (
        "n3:D=-B3_over_Q",
        "square_class=[-B3]_not_[B3]",
        "odd_divisor_carrier=B3",
    )
    semantic_surface = (
        "engine=python-flint-0.9.0;ring=Q[w,P,Q,C,X,Y,Z]",
        "independent_rows=n3,n4;x,y;extra_z_at_n3",
        "exact_FLINT_resultants_discriminants_factorization_gcd_square_roots",
        rows,
        z_row,
        hostile,
        unit_hostile,
    )
    semantic_sha256 = sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    require(
        semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
        (semantic_sha256, EXPECTED_SEMANTIC_SHA256),
    )

    print("THM-3494 independent weighted discriminant atlas audit")
    print("engine=python-flint-0.9.0;ring=Q[w,P,Q,C,X,Y,Z]")
    for row in rows:
        print(f"degree_row={row}")
    print(f"cubic_z_row={z_row}")
    print(f"unit_hostile={unit_hostile}")
    print(f"xor_hostile={hostile}")
    print("verdict=PASS_with_unit_and_XOR_scope_repairs;all_degree_proof_audited_separately")
    print(f"semantic_sha256={semantic_sha256}")


if __name__ == "__main__":
    main()
