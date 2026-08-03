#!/usr/bin/env python3
"""Exact critical-cover/cofactor probe for the coupled affine JC clutch.

This is a non-canonical research artifact.  It studies the two fixed controls

    C=x+1,  E'=1

in the two THM-3212 cubic accessory fields.  It never calls
``sympy.resultant`` for a load-bearing calculation: resultants are literal
Sylvester determinants, and the inverse section is independently recovered
from a cubic-to-quadratic pseudo-remainder identity.

The sidecar produced here is an *elimination* cofactor for the gradient pair,
not the primitive-element/Keller cofactor of THM-3064.  Keeping that type
distinction visible is the point of the probe.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


EXPECTED_H_DIGESTS = {
    "4111": "ae98818930edc795b9a7a73a957ddf2f2e70dee40695bb3ea489c9a8c16780d4",
    "3211": "621f57adfa31900c88fa44d0443d9baf81d199a201a4b9963d99914d134b1d21",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def literal_sylvester_resultant(first, second, variable):
    """Standard Res(first,second), from the literal Sylvester determinant."""

    first_poly = sp.Poly(first, variable)
    second_poly = sp.Poly(second, variable)
    first_degree = first_poly.degree()
    second_degree = second_poly.degree()
    first_row = first_poly.all_coeffs()
    second_row = second_poly.all_coeffs()
    zero = sp.S.Zero
    rows = []
    for shift in range(second_degree):
        rows.append(
            [zero] * shift
            + first_row
            + [zero] * (second_degree - 1 - shift)
        )
    for shift in range(first_degree):
        rows.append(
            [zero] * shift
            + second_row
            + [zero] * (first_degree - 1 - shift)
        )
    return sp.expand(sp.Matrix(rows).det(method="domain-ge"))


def audit_sympy_resultant_sign() -> None:
    """Expose the installed unequal-odd-degree swap bug, then avoid it."""

    x, a, c = sp.symbols("x a c")
    low = x - a
    high = x**3 + c
    expected_low_high = a**3 + c
    determinant_low_high = literal_sylvester_resultant(low, high, x)
    determinant_high_low = literal_sylvester_resultant(high, low, x)
    sympy_low_high = sp.resultant(low, high, x)
    sympy_high_low = sp.resultant(high, low, x)
    require(
        sp.expand(determinant_low_high - expected_low_high) == 0,
        "literal Sylvester orientation",
    )
    require(
        sp.expand(determinant_high_low + expected_low_high) == 0,
        "literal Sylvester antisymmetry",
    )
    low_wrong = sp.expand(sympy_low_high - expected_low_high) != 0
    high_correct = sp.expand(sympy_high_low + expected_low_high) == 0
    antisymmetry_broken = sp.expand(sympy_low_high + sympy_high_low) != 0
    print(
        f"sympy_resultant_hazard=(version={sp.__version__},"
        f"low_odd_first_wrong={int(low_wrong)},"
        f"high_odd_first_correct={int(high_correct)},"
        f"swap_antisymmetry_broken={int(antisymmetry_broken)});"
        "load_bearing_path=literal_sylvester_only"
    )


def derive_universal_sidecar():
    """Derive the resultant, linear subresultant, and Bezout sidecar."""

    y, A, C, V, derivative_V, d, k = sp.symbols("y A C V Vp d k")
    L = y**2 + y + C * V
    first = 2 * L * (2 * y + 1) + V * A
    second = V**3 * k + V**2 * y + L * (-derivative_V * y + 2 * V**2 * d)

    resultant = literal_sylvester_resultant(first, second, y)
    critical_factor = sp.cancel(resultant / V**3)
    require(
        sp.expand(resultant - V**3 * critical_factor) == 0,
        "universal V^3 factor",
    )
    critical_poly = sp.Poly(
        critical_factor, A, C, V, derivative_V, d, k
    )
    require(len(critical_poly.terms()) == 40, "universal 40-term factor")

    # Since first has constant leading y-coefficient 4, this row operation
    # cancels the cubic term of second without a denominator.
    quadratic = sp.expand(4 * second + derivative_V * first)
    quadratic_poly = sp.Poly(quadratic, y)
    require(quadratic_poly.degree() == 2, "cubic-to-quadratic cancellation")
    a2, b2, c2 = quadratic_poly.all_coeffs()
    require(
        sp.expand(a2 - (2 * derivative_V + 8 * V**2 * d)) == 0,
        "quadratic leading coefficient",
    )

    first_poly = sp.Poly(first, y)
    leading, p2, q2, r2 = first_poly.all_coeffs()
    require(leading == 4 and p2 == 6, "normalized first cubic")

    # Fraction-free first pseudo-remainder.  If first=second=0 then
    # ell_1*y+ell_0=0.  This is the selected critical-root coordinate.
    pivot = sp.expand(a2 * p2 - 4 * b2)
    ell1 = sp.expand(a2 * (a2 * q2 - 4 * c2) - pivot * b2)
    ell0 = sp.expand(a2**2 * r2 - pivot * c2)
    quotient = sp.expand(4 * a2 * y + pivot)
    cofactor_first = sp.expand(a2**2 - derivative_V * quotient)
    cofactor_second = sp.expand(-4 * quotient)
    linear = sp.expand(ell1 * y + ell0)
    require(
        sp.expand(linear - cofactor_first * first - cofactor_second * second)
        == 0,
        "linear Bezout sidecar",
    )
    require(
        sp.expand(cofactor_first - derivative_V * cofactor_second / 4 - a2**2)
        == 0,
        "cofactor-pair unit identity",
    )

    # Independent norm identity.  It recomputes the same resultant from the
    # quadratic and its linear pseudo-remainder, with no second determinant.
    norm_numerator = sp.expand(
        a2 * ell0**2 - b2 * ell0 * ell1 + c2 * ell1**2
    )
    require(
        sp.expand(norm_numerator - 16 * a2**2 * resultant) == 0,
        "quadratic-linear norm equals cubic resultant",
    )

    variables = (A, C, V, derivative_V, d, k)
    print(
        "universal_sidecar=resultant=V^3*K40;"
        "h=4R2+Vp*R1=a*y^2+b*y+c;"
        "ell=l1*y+l0=U*R1+W*R2;"
        "U-(Vp/4)*W=a^2;"
        "a*l0^2-b*l0*l1+c*l1^2=16*a^2*Res(R1,R2)"
    )
    return {
        "variables": variables,
        "critical": critical_poly,
    }


def specialize_multivariate(poly, values, polynomial_ring):
    result = polynomial_ring.zero
    for monomial, coefficient in poly.terms():
        term = polynomial_ring.convert(coefficient)
        for value, exponent in zip(values, monomial):
            term *= value**exponent
        result += term
    return result


def canonical_poly_text(poly) -> str:
    return "|".join(
        f"{degree}:{coefficient}"
        for (degree,), coefficient in sorted(poly.to_dict().items())
    )


def poly_digest(poly) -> str:
    return sha256(canonical_poly_text(poly.monic()).encode("ascii")).hexdigest()


def projective_pair_digest(left, right) -> str:
    """Hash a polynomial pair modulo one *common* nonzero scalar.

    Normalizing the entries separately would erase their relative scale and
    therefore would not determine the rational section ``-right/left``.
    """

    require(left != 0, "the section denominator is zero")
    common_scale = left.LC
    normalized_left = left.quo_ground(common_scale)
    normalized_right = right.quo_ground(common_scale)
    return sha256(
        (
            canonical_poly_text(normalized_left)
            + "\n"
            + canonical_poly_text(normalized_right)
        ).encode("ascii")
    ).hexdigest()


def build_accessory_case(name: str, universal) -> tuple[str, str]:
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
    else:
        accessory_v = (24 * alpha**2 - 16 * alpha - 16) / 21
        shift = (5 * alpha - 4) / 7
        A0 = 9 * accessory_v**2 * (5 * alpha - 4) / 343
        extras = (6, 3)

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
    owner = S * T
    require(V.degree() == 16 and A.degree() == 8, f"{name} response degrees")
    require(
        2 * V * A.diff(X) - A * derivative_V == 2 * V,
        f"{name} response identity",
    )

    C = X + 1
    values = (A, C, V, derivative_V, field.one, field.one)
    K = specialize_multivariate(universal["critical"], values, x_ring)
    # Rebuild the pseudo-remainder after specialization.  This is exactly the
    # universal formula already checked above, but it avoids expanding large
    # multivariate intermediates before reducing to one variable.
    a2 = 2 * derivative_V + 8 * V**2
    b2 = 2 * derivative_V + 12 * V**2
    c2 = V * (4 * V**2 + 8 * C * V**2 + derivative_V * (2 * C + A))
    q2 = 2 + 4 * C * V
    r2 = V * (2 * C + A)
    pivot = 6 * a2 - 4 * b2
    ell1 = a2 * (a2 * q2 - 4 * c2) - pivot * b2
    ell0 = a2**2 * r2 - pivot * c2
    boundary = S**3 * T**8 * X ** extras[0] * (X - 1) ** extras[1]
    H = K.exquo(boundary)
    require((K.degree(), boundary.degree(), H.degree()) == (97, 44, 53), f"{name} degrees")
    require(H.gcd(owner).degree() == 0, f"{name} boundary disjoint")
    require(H.gcd(H.diff(X)).degree() == 0, f"{name} squarefree H")
    require(H.gcd(a2).degree() == 0, f"{name} quadratic pivot unit")
    require(H.gcd(ell1).degree() == 0, f"{name} linear pivot unit")
    require(poly_digest(H) == EXPECTED_H_DIGESTS[name], f"{name} H digest")

    h_digest = poly_digest(H)
    section_digest = projective_pair_digest(ell1, ell0)
    print(
        f"case={name};control=(C,Eprime)=(1+x,1);"
        "critical_scheme=Spec(K[x]/H)_graph_y=-l0/l1;"
        f"degrees=(H,a,l1,l0,Y)=({H.degree()},{a2.degree()},"
        f"{ell1.degree()},{ell0.degree()},fractional);"
        "gcds=(H,Hprime)=(1);(H,owner)=(1);(H,a)=(1);(H,l1)=(1);"
        f"H_digest={h_digest};section_projective_pair_digest={section_digest}"
    )
    return h_digest, section_digest


def hostile_checks() -> None:
    x, y = sp.symbols("x y")

    # If the resultant is not squarefree, the linear pivot can vanish on its
    # reduced support: three coincident y-branches sit over x=0.
    first = 4 * y**3
    quadratic = y**2 + x
    resultant = literal_sylvester_resultant(first, quadratic, y)
    a2, b2, c2 = sp.Poly(quadratic, y).all_coeffs()
    _, p2, q2, r2 = sp.Poly(first, y).all_coeffs()
    pivot = a2 * p2 - 4 * b2
    ell1 = sp.expand(a2 * (a2 * q2 - 4 * c2) - pivot * b2)
    require(sp.factor(resultant - 16 * x**3) == 0, "multiple hostile resultant")
    require(sp.factor(ell1 + 4 * x) == 0, "multiple hostile linear pivot")

    # A squarefree resultant need not make the quadratic leading coefficient
    # a unit.  Root recovery survives, but the one-chart cofactor unit identity
    # cannot be inferred from squarefreeness alone.
    first2 = 4 * (y**3 - 1)
    quadratic2 = x * y**2 + y - 1
    resultant2 = literal_sylvester_resultant(first2, quadratic2, y)
    require(
        sp.factor(resultant2 - 16 * x * (x**2 + 3)) == 0,
        "pivot hostile resultant",
    )
    require(sp.gcd(sp.Poly(resultant2, x), sp.Poly(x, x)).degree() == 1, "pivot hostile contact")
    print(
        "hostiles=non_squarefree:Res(4y^3,y^2+x)=16x^3_and_l1=-4x;"
        "squarefree_but_a_nonunit:Res(4(y^3-1),x*y^2+y-1)=16x(x^2+3),a=x;"
        "separation=root_inverse_needs_linear_pivot;unimodular_cofactor_chart_needs_a_sidecar"
    )


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
    print("critical inverse-cover and elimination-cofactor probe")
    audit_sympy_resultant_sign()
    universal = derive_universal_sidecar()
    digests = tuple(build_accessory_case(name, universal) for name in ("4111", "3211"))
    hostile_checks()
    source_audit()
    print(f"case_digests={digests}")
    print(
        "consequence=fixed_squarefree_degree53_controls_have_exact_reduced_critical_graphs_and_unimodular_elimination_cofactor_pairs"
    )
    print(
        "scope=FINITE-EXACT_two_controls;not_universal_affine_CE;not_Keller_primitive_cofactor;not_polynomial_inverse_cover;not_JC2"
    )


if __name__ == "__main__":
    main()
