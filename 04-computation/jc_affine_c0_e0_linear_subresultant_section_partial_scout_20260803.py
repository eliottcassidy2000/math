#!/usr/bin/env python3
"""Finite-exact scout for the THM-3289 rational control's critical section.

For each THM-3212 cubic accessory field, specialize the THM-3289 localized
gradient pair to C=1+x and E'=k=1.  The script computes a chosen standard
subresultant sequence in y, keeps the penultimate linear row only up to units,
and removes its exact x-content.  It then works over the squarefree saturated
resultant H and decides whether the linear coefficient is a unit there.

This is a FINITE-EXACT PARTIAL artifact for two fixed rational controls.  It is
not a theorem about all affine parameters, constructs no inverse cover, and
proves neither JC(2) nor any case of the planar Jacobian conjecture.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py":
        "03d121e57dd2edaece7cd8693d792349f03757c6e781eb5d9d0c897fcc889448",
    "05-knowledge/results/jc_centered_heptic_source_morse_obstruction_thm3212.out":
        "729e0c7b9fa51fa5c4ac5f18f50dc4413c8a8bb7bf5f0ebf2a7709650304bc85",
    "04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.py":
        "f63ff06e3f5ed30f3f6bc5be99756c347d6af5f8e9b220ce8336abff2cd2ca31",
    "05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.out":
        "1aef4341650cdfaf1043a8699e3a1725a0100af6d9848d99dfa924b6f054dba1",
    "04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289_independent_audit.py":
        "b2fa8c96854549ccb9e515485214c119b685b31456fb7c53c5e2bd83f7933831",
    "05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289_independent_audit.out":
        "48d50289c98c9dd17e099497d21cd9648cac27a097b339a1f75a1e13d8fd8837",
}
EXPECTED_DIGESTS = {
    "4111": (
        "22879e1cb10e47687ac8a0fa72473bbb3b8e6f8b56879f3f5672d3b45919e2b9",
        "c68939482975e3864fff98e8f573f78b7de72f9913ff47a30af2d778d683bf3d",
        "ae98818930edc795b9a7a73a957ddf2f2e70dee40695bb3ea489c9a8c16780d4",
    ),
    "3211": (
        "d0e81357fd33bf829d1ef2ff4be9508408c97e037765d91d8013103ffb68061a",
        "fafed43f4178e78f4f85e3cf13e856541cedfec21d2f1db4b341610b1b87a3a4",
        "621f57adfa31900c88fa44d0443d9baf81d199a201a4b9963d99914d134b1d21",
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def literal_cubic_sylvester(first, second, variable):
    """Return the standard 6-by-6 Sylvester determinant of two cubics."""

    first_poly = sp.Poly(first, variable)
    second_poly = sp.Poly(second, variable)
    require(first_poly.degree() == 3, "first Sylvester input is cubic")
    require(second_poly.degree() == 3, "second Sylvester input is cubic")
    first_row = first_poly.all_coeffs()
    second_row = second_poly.all_coeffs()
    zero = sp.S.Zero
    matrix = sp.Matrix(
        [
            first_row + [zero, zero],
            [zero] + first_row + [zero],
            [zero, zero] + first_row,
            second_row + [zero, zero],
            [zero] + second_row + [zero],
            [zero, zero] + second_row,
        ]
    )
    return sp.expand(matrix.det(method="domain-ge"))


def canonical_poly_text(poly) -> str:
    return "|".join(
        f"{degree}:{coefficient}"
        for (degree,), coefficient in sorted(poly.to_dict().items())
    )


def poly_digest(poly) -> str:
    return sha256(canonical_poly_text(poly).encode("ascii")).hexdigest()


def pair_digest(first, second) -> str:
    payload = canonical_poly_text(first) + "\n" + canonical_poly_text(second)
    return sha256(payload.encode("ascii")).hexdigest()


def specialize_sparse(poly, base_values, polynomial_ring):
    result = polynomial_ring.zero
    for monomial, coefficient in poly.terms():
        term = polynomial_ring.convert(coefficient)
        for base, exponent in zip(base_values, monomial):
            term *= base**exponent
        result += term
    return result


def derive_generic_packet():
    y, A, C, V, derivative_V, d, k = sp.symbols("y A C V Vp d k")
    L = y**2 + y + C * V
    R1 = sp.expand(2 * L * (2 * y + 1) + V * A)
    R2 = sp.expand(
        V**3 * k + V**2 * y + L * (-derivative_V * y + 2 * V**2 * d)
    )

    sequence = sp.subresultants(R1, R2, y)
    degree_profile = tuple(sp.Poly(row, y).degree() for row in sequence)
    require(degree_profile == (3, 3, 2, 1, 0), "standard subresultant degrees")

    resultant = sp.Poly(sequence[-1], y).nth(0)
    require(
        sp.expand(resultant - literal_cubic_sylvester(R1, R2, y)) == 0,
        "standard subresultant resultant equals literal Sylvester determinant",
    )
    quotient = sp.cancel(resultant / V**3)
    require(sp.denom(quotient) == 1, "resultant quotient is polynomial")
    K_poly = sp.Poly(quotient, A, C, V, derivative_V, d, k)
    require(len(K_poly.terms()) == 40, "THM-3289 quotient has forty terms")

    linear = sp.Poly(sequence[-2], y)
    raw_a = linear.nth(1)
    raw_b = linear.nth(0)
    require(sp.cancel(raw_a / (2 * V)).is_polynomial(), "generic linear a has 2V content")
    require(sp.cancel(raw_b / (2 * V)).is_polynomial(), "generic linear b has 2V content")

    # A chosen standard subresultant row, not a canonical raw PRS iterate.
    # The identity supplies its invariant ideal meaning; later normalizations
    # are used only up to nonzero field units, as required by MISTAKE-360.
    u1 = -2 * derivative_V * (4 * V**2 * d + derivative_V)
    u0 = 4 * V**2 * (4 * V**2 * d**2 + derivative_V * d + derivative_V)
    w1 = -8 * (4 * V**2 * d + derivative_V)
    w0 = -4 * (4 * V**2 * d - 4 * V**2 + derivative_V)
    bezout_row = sp.expand((u1 * y + u0) * R1 + (w1 * y + w0) * R2)
    require(sp.expand(bezout_row - linear.as_expr()) == 0, "linear Bezout row")
    require(sp.Poly(R1, y).LC() == 4, "R1 has constant projective leading coefficient")

    print(
        "generic_subresultants=(3,3,2,1,0);"
        "S0=Res_y(R1,R2)=V^3*K40;"
        "S1=U*R1+W*R2;generic_coefficient_content_contains=2V;"
        "raw_rows_used_only_up_to_units=1;MISTAKE360_respected=1"
    )
    print(
        "bezout_multipliers="
        "U=(-2Vp(4V^2d+Vp))*y+4V^2(4V^2d^2+Vp*d+Vp);"
        "W=(-8(4V^2d+Vp))*y-4(4V^2d-4V^2+Vp)"
    )
    print("projective_hostile=R1_y_lead_4=>no_common_y_infinity")
    return (
        K_poly,
        sp.Poly(raw_a, A, C, V, derivative_V, d, k),
        sp.Poly(raw_b, A, C, V, derivative_V, d, k),
    )


def build_case(name: str, generic_packet) -> tuple[str, str, str]:
    K_poly, raw_a_poly, raw_b_poly = generic_packet
    u, x = sp.symbols("u x")
    if name == "4111":
        accessory = sp.Poly(
            100 * u**3 + 244 * u**2 + 237 * u + 44, u, domain=QQ
        )
        exponent_a, exponent_b = 4, 1
    else:
        accessory = sp.Poly(
            75 * u**3 - 89 * u**2 - 31 * u + 61, u, domain=QQ
        )
        exponent_a, exponent_b = 3, 2

    field = QQ.alg_field_from_poly(accessory, alias="u")
    alpha = field.ext
    x_ring = field.poly_ring(x)
    X = x_ring.gens[0]
    if name == "4111":
        accessory_v = (8 * alpha**2 + 9 * alpha + 8) / 7
        shift = 5 * (alpha + 1) / 7
        A0 = 80 * accessory_v**2 * (alpha + 1) / 343
        extras = (9, 0)
        boundary_label = "S^3*T^8*x^9"
    else:
        accessory_v = (24 * alpha**2 - 16 * alpha - 16) / 21
        shift = (5 * alpha - 4) / 7
        A0 = 9 * accessory_v**2 * (5 * alpha - 4) / 343
        extras = (6, 3)
        boundary_label = "S^3*T^8*x^6*(x-1)^3"

    gamma = -7 * A0
    quadratic = X**2 - alpha * X + accessory_v
    D = X**exponent_a * (X - 1) ** exponent_b * quadratic
    T = X * (X - 1) * quadratic
    E = (
        exponent_a * (X - 1) * quadratic
        + exponent_b * X * quadratic
        + X * (X - 1) * (2 * X - alpha)
    ) / 7
    S = X + shift
    V = 4 * S * D * T**2 / gamma**2
    A = 2 * S * E * T / gamma
    derivative_V = V.diff(X)
    C = X + 1
    boundary = S**3 * T**8 * X**extras[0] * (X - 1) ** extras[1]
    require(boundary.degree() == 44, f"{name} boundary degree")
    require(
        2 * V * A.diff(X) - A * derivative_V == 2 * V,
        f"{name} response identity",
    )

    base_values = (A, C, V, derivative_V, field.one, field.one)
    control_K = specialize_sparse(K_poly, base_values, x_ring)
    H = control_K.exquo(boundary)
    require((control_K.degree(), H.degree()) == (97, 53), f"{name} degrees")

    raw_a = specialize_sparse(raw_a_poly, base_values, x_ring)
    raw_b = specialize_sparse(raw_b_poly, base_values, x_ring)
    primitive_a = raw_a.exquo(boundary)
    primitive_b = raw_b.exquo(boundary)
    require(
        (raw_a.degree(), raw_b.degree(), primitive_a.degree(), primitive_b.degree())
        == (80, 88, 36, 44),
        f"{name} linear row degrees",
    )

    # Normalize the pair by one common field unit.  Equality of this chosen
    # representative with another PRS iterate is neither asserted nor used.
    leading_unit = primitive_a.LC
    a = primitive_a.monic()
    b = primitive_b * (field.one / leading_unit)
    coefficient_gcd = a.gcd(b)
    require(coefficient_gcd.degree() == 0, f"{name} primitive coefficient pair")
    require(
        raw_a == boundary * leading_unit * a
        and raw_b == boundary * leading_unit * b,
        f"{name} exact boundary content",
    )

    gcd_a_H = a.gcd(H)
    gcd_b_H = b.gcd(H)
    triple_gcd = coefficient_gcd.gcd(H)
    repeated_gcd = H.gcd(H.diff(X))
    boundary_gcd = H.gcd(boundary)
    V_gcd = H.gcd(V)
    require(
        all(
            item.degree() == 0
            for item in (
                gcd_a_H,
                gcd_b_H,
                triple_gcd,
                repeated_gcd,
                boundary_gcd,
                V_gcd,
            )
        ),
        f"{name} hostile gcd gates",
    )

    inverse_coefficient, _, inverse_gcd = a.gcdex(H)
    require(inverse_gcd.degree() == 0, f"{name} a inverse exists modulo H")
    inverse_coefficient *= field.one / inverse_gcd.LC
    inverse_coefficient %= H
    require((a * inverse_coefficient) % H == x_ring.one, f"{name} a inverse")
    y_section = (-b * inverse_coefficient) % H
    require((a * y_section + b) % H == x_ring.zero, f"{name} linear section")

    section_L = y_section**2 + y_section + C * V
    reduced_R1 = (2 * section_L * (2 * y_section + 1) + V * A) % H
    reduced_R2 = (
        V**3
        + V**2 * y_section
        + section_L * (-derivative_V * y_section + 2 * V**2)
    ) % H
    require(
        reduced_R1 == x_ring.zero and reduced_R2 == x_ring.zero,
        f"{name} direct gradient-row reductions",
    )

    row_hash = pair_digest(a, b)
    section_hash = poly_digest(y_section)
    H_hash = poly_digest(H.monic())
    expected_row, expected_section, expected_H = EXPECTED_DIGESTS[name]
    if expected_row != "TBD":
        require(row_hash == expected_row, f"{name} row digest drift")
    if expected_section != "TBD":
        require(section_hash == expected_section, f"{name} section digest drift")
    require(H_hash == expected_H, f"{name} inherited H digest drift")

    print(
        f"case={name};control=(C,k)=(1+x,1);"
        "raw_S1=unit*boundary*(a*y+b);"
        f"boundary={boundary_label};boundary_degree={boundary.degree()};"
        "raw_degrees=(80,88);primitive_degrees=(36,44);"
        "primitive_coefficient_gcd_degree=0"
    )
    print(
        f"case={name};H_degree={H.degree()};"
        "gcd_degrees=(aH,bH,abH,H_Hprime,H_boundary,H_V)=(0,0,0,0,0,0);"
        f"section_degree={y_section.degree()};"
        "R1_section_mod_H=0;R2_section_mod_H=0"
    )
    print(
        f"case={name};hostiles="
        "a_zero_on_H:none|projective_y_infinity:none|"
        "boundary_intersection:none|repeated_H:none;"
        "selector_split_required=0;global_reduced_H_section=1;"
        f"row_digest={row_hash};section_digest={section_hash};H_digest={H_hash}"
    )
    return row_hash, section_hash, H_hash


def source_audit() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assert_nodes == 0 and float_literals == 0, "source AST gate")
    print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")


def main() -> None:
    print("finite-exact partial linear-subresultant critical-section scout")
    for dependency, expected_hash in DEPENDENCIES.items():
        require(
            lf_hash(ROOT / dependency) == expected_hash,
            f"dependency drift: {dependency}",
        )
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    generic_packet = derive_generic_packet()
    case_digests = tuple(
        build_case(name, generic_packet) for name in ("4111", "3211")
    )
    print(f"case_digests={case_digests}")
    print(
        "ideal_consequence=over_each_reduced_A_H=K_i[x]/(H),"
        "boundary_and_a_are_units_and_(R1,R2)=(a*y+b);"
        "unique_y=-b/a_is_one_global_section_not_a_selector_split"
    )
    print(
        "scope=FINITE-EXACT_PARTIAL_two_rational_controls_only;"
        "subresultant_representatives_only_up_to_units;"
        "no_inverse_cover_not_JC2"
    )
    source_audit()


if __name__ == "__main__":
    main()
