#!/usr/bin/env python3
"""Exact referee for THM-2328's full B--W ratio-bank closure.

The two roots of the quadratic ratio factor form one smooth-triple-fibre
orbit.  The six roots of the sextic factor form one ordinary-node orbit.
Arithmetic at the exceptional fibres is performed directly in SymPy's
algebraic-number representation; this avoids expensive reconstruction of
minimal polynomials for large degree-six expressions.  Every check remains
active under optimized Python.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive_integer_vector(element) -> tuple[int, ...]:
    """Clear denominators and content in an ANP power-basis vector."""
    coefficients = [sp.Rational(value) for value in element.to_list()]
    denominator = 1
    for value in coefficients:
        denominator = sp.ilcm(denominator, int(value.q))
    integers = [int(value * denominator) for value in coefficients]
    content = 0
    for value in integers:
        content = sp.igcd(content, abs(value))
    require(content != 0, "requested nonzero field vector vanished")
    integers = [value // content for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def spectral_curve(
    u: sp.Symbol,
    y: sp.Symbol,
    bvar: sp.Expr,
    wvar: sp.Expr,
) -> sp.Expr:
    """THM-2297's G_0 at (B,C,D,W)=(bvar,0,0,wvar)."""
    return sp.expand(
        -5878656 * wvar * y
        - 26040609 * u**3
        + 49601160 * bvar * u**2
        + 1607445 * u**2 * y**2
        - 20995200 * bvar**2 * u
        - 2857680 * bvar * u * y**2
        - 138915 * u * y**4
        + 777600 * bvar**2 * y**2
        + 78120 * bvar * y**4
        + 1127 * y**6
    )


def polynomial_power_mod(
    base: sp.Poly,
    exponent: int,
    modulus: sp.Poly,
) -> sp.Poly:
    """Binary polynomial powering over a finite field."""
    result = sp.Poly(1, *modulus.gens, domain=modulus.domain)
    power = base.rem(modulus)
    while exponent:
        if exponent & 1:
            result = (result * power).rem(modulus)
        power = (power * power).rem(modulus)
        exponent //= 2
    return result


def require_global_irreducibility_certificate(
    curve: sp.Expr,
    ratio_polynomial: sp.Expr,
    u: sp.Symbol,
    y: sp.Symbol,
    t: sp.Symbol,
    label: str,
) -> None:
    """Independently rule out every polynomial root u=a*y^2+b*y+c."""
    aa, bb, cc = sp.symbols(f"aa_{label} bb_{label} cc_{label}")
    substituted = sp.Poly(
        sp.expand(curve.subs(u, aa * y**2 + bb * y + cc)),
        y,
    )
    equations = [
        substituted.coeff_monomial(y**degree) for degree in range(7)
    ]
    basis = sp.groebner(
        equations + [ratio_polynomial],
        aa,
        bb,
        cc,
        t,
        order="grevlex",
        domain=sp.QQ,
    )
    require(
        len(basis.polys) == 1 and basis.polys[0].as_expr() == 1,
        f"{label} orbit gained a polynomial u-root",
    )


def cubic_coefficients(field: sp.AlgebraicField, branch_value):
    """Return the specialized cubic coefficients as native field elements."""
    q = field.unit
    k = field.convert
    yy = branch_value
    a3 = k(-26040609)
    a2 = k(49601160) * q + k(1607445) * yy**2
    a1 = (
        k(-20995200) * q**2
        - k(2857680) * q * yy**2
        - k(138915) * yy**4
    )
    a0 = (
        k(-5878656) * q**3 * yy
        + k(777600) * q**2 * yy**2
        + k(78120) * q * yy**4
        + k(1127) * yy**6
    )
    return a3, a2, a1, a0


def curve_y_derivative(field: sp.AlgebraicField, root, branch_value):
    """Evaluate partial_y G_0 using native number-field arithmetic."""
    q = field.unit
    k = field.convert
    rr = root
    yy = branch_value
    return (
        k(-5878656) * q**3
        + k(3214890) * rr**2 * yy
        - k(5715360) * q * rr * yy
        - k(555660) * rr * yy**3
        + k(1555200) * q**2 * yy
        + k(312480) * q * yy**3
        + k(6762) * yy**5
    )


def tangent_data(field: sp.AlgebraicField, root, branch_value):
    """Return the U^2 coefficient and tangent-cone discriminant."""
    q = field.unit
    k = field.convert
    rr = root
    yy = branch_value
    tangent_u2 = k(-3 * 26040609) * rr + (
        k(49601160) * q + k(1607445) * yy**2
    )
    tangent_uy = (
        k(6429780) * rr * yy
        - k(5715360) * q * yy
        - k(555660) * yy**3
    )
    tangent_y2 = (
        k(1607445) * rr**2
        - k(2857680) * q * rr
        - k(833490) * rr * yy**2
        + k(777600) * q**2
        + k(468720) * q * yy**2
        + k(16905) * yy**4
    )
    tangent_discriminant = (
        tangent_uy**2 - k(4) * tangent_u2 * tangent_y2
    )
    return tangent_u2, tangent_discriminant


def orbit_branch_data(
    ratio_polynomial: sp.Expr,
    expected_degree: int,
    expected_local_type: str,
    u: sp.Symbol,
    y: sp.Symbol,
    t: sp.Symbol,
):
    """Check one full Galois orbit and return its normalization genus."""
    alpha = sp.CRootOf(ratio_polynomial, 0)
    field = sp.QQ.algebraic_field(alpha)
    q = field.unit
    k = field.convert
    zero = field.zero

    require(q != zero, "ratio-field generator unexpectedly vanished")
    require(
        (q**3) ** 2 / q**5 == q,
        "representative (B,W)=(alpha,alpha^3) has the wrong ratio",
    )

    curve = spectral_curve(u, y, alpha, alpha**3)
    discriminant = sp.Poly(
        sp.discriminant(curve, u),
        y,
        extension=alpha,
    )
    repeated = sp.gcd(discriminant, discriminant.diff()).monic()
    require(repeated.degree() == 1, "repeated branch factor is not linear")
    branch_value = -repeated.rep.nth(0)
    quotient, remainder = discriminant.div(repeated**2)
    require(remainder.is_zero, "double branch factor does not divide")
    residual = quotient.monic()
    require(
        quotient.rep.LC() == k(-153384762202971019112448),
        "residual branch leading coefficient changed",
    )
    require(
        residual.degree() == 10
        and sp.gcd(residual, residual.diff()).degree() == 0
        and sp.gcd(residual, repeated).degree() == 0,
        "degree-ten branch residual is not squarefree and coprime",
    )

    a3, a2, a1, a0 = cubic_coefficients(field, branch_value)
    cubic_i = a2**2 - k(3) * a3 * a1
    cubic_j = (
        k(2) * a2**3
        - k(9) * a3 * a2 * a1
        + k(27) * a3**2 * a0
    )

    if expected_local_type == "smooth_triple":
        require(
            cubic_i == zero and cubic_j == zero,
            "exceptional cubic is not a triple-root fibre",
        )
        repeated_root = -a2 / (k(3) * a3)
        require(
            a2 == -k(3) * a3 * repeated_root
            and a1 == k(3) * a3 * repeated_root**2
            and a0 == -a3 * repeated_root**3,
            "triple-root factorization changed",
        )
        derivative_y = curve_y_derivative(
            field,
            repeated_root,
            branch_value,
        )
        require(
            derivative_y != zero,
            "triple-root fibre is singular instead of smooth",
        )
        local_certificate = {
            "derivative_y_vector": primitive_integer_vector(derivative_y),
        }
        total_ramification = residual.degree() + 2
    elif expected_local_type == "ordinary_node":
        require(
            cubic_i != zero and cubic_j != zero,
            "node orbit degenerated to a triple-root fibre",
        )
        repeated_root = (
            k(9) * a3 * a0 - a2 * a1
        ) / (k(2) * cubic_i)
        simple_root = -a2 / a3 - k(2) * repeated_root
        require(
            a2 == -a3 * (simple_root + k(2) * repeated_root)
            and a1
            == a3
            * (
                k(2) * repeated_root * simple_root
                + repeated_root**2
            )
            and a0 == -a3 * repeated_root**2 * simple_root,
            "double-root fibre factorization changed",
        )
        curve_value = (
            a3 * repeated_root**3
            + a2 * repeated_root**2
            + a1 * repeated_root
            + a0
        )
        curve_u = (
            k(3) * a3 * repeated_root**2
            + k(2) * a2 * repeated_root
            + a1
        )
        curve_y = curve_y_derivative(
            field,
            repeated_root,
            branch_value,
        )
        require(
            curve_value == zero and curve_u == zero and curve_y == zero,
            "exceptional repeated point is not a curve singularity",
        )
        tangent_u2, tangent_discriminant = tangent_data(
            field,
            repeated_root,
            branch_value,
        )
        require(
            tangent_u2 != zero and tangent_discriminant != zero,
            "exceptional point is not an ordinary nonvertical node",
        )
        local_certificate = {
            "tangent_u2_vector": primitive_integer_vector(tangent_u2),
            "tangent_discriminant_vector": primitive_integer_vector(
                tangent_discriminant
            ),
        }
        total_ramification = residual.degree()
    else:
        raise RuntimeError("unknown expected local type")

    genus = (total_ramification - 4) // 2
    require(
        expected_degree == sp.degree(ratio_polynomial, t),
        "ratio orbit degree changed",
    )
    return {
        "degree": expected_degree,
        "local_type": expected_local_type,
        "residual_degree": residual.degree(),
        "total_ramification": total_ramification,
        "genus": genus,
        "local_certificate": local_certificate,
    }


def main() -> None:
    u, y, v, z, t = sp.symbols("u y v z t")

    quadratic_ratio = (
        5313800000
        + 4508659468656 * t
        - 136738899331083 * t**2
    )
    sextic_ratio = (
        5511577600000000000000000000
        + 4983290602536960000000000000000 * t
        - 6564822237254419568640000000000 * t**2
        - 3094052863483309848285092659200000 * t**3
        - 81862566455344350924421142159812608 * t**4
        - 744088924275617882256518828471658624 * t**5
        - 2973811237322720333634598763466407943 * t**6
    )

    quadratic_mod_19 = sp.Poly(quadratic_ratio, t, modulus=19).monic()
    require(
        quadratic_mod_19 == sp.Poly(t**2 + 4 * t + 5, t, modulus=19)
        and all(
            quadratic_mod_19.eval(residue) != 0
            for residue in range(19)
        ),
        "quadratic ratio factor is not irreducible modulo nineteen",
    )

    sextic_mod_13 = sp.Poly(sextic_ratio, t, modulus=13).monic()
    expected_sextic_mod_13 = sp.Poly(
        t**6 + t**5 - t**4 - t**3 - 2 * t**2 - t + 4,
        t,
        modulus=13,
    )
    finite_t = sp.Poly(t, t, modulus=13)
    frobenius_2 = polynomial_power_mod(finite_t, 13**2, sextic_mod_13)
    frobenius_3 = polynomial_power_mod(finite_t, 13**3, sextic_mod_13)
    frobenius_6 = polynomial_power_mod(finite_t, 13**6, sextic_mod_13)
    require(
        sextic_mod_13 == expected_sextic_mod_13
        and sp.gcd(
            sextic_mod_13,
            frobenius_2 - finite_t,
        ).degree()
        == 0
        and sp.gcd(
            sextic_mod_13,
            frobenius_3 - finite_t,
        ).degree()
        == 0
        and (frobenius_6 - finite_t).rem(sextic_mod_13).is_zero,
        "sextic ratio factor is not irreducible modulo thirteen",
    )

    infinity_polynomial = (
        1127 - 138915 * v + 1607445 * v**2 - 26040609 * v**3
    )
    infinity_discriminant = sp.discriminant(infinity_polynomial, v)
    require(
        infinity_discriminant == -153384762202971019112448,
        "weighted infinity cubic is not separable",
    )

    # Verify the elementary absolute-irreducibility coefficient mechanism.
    aa, bb, cc = sp.symbols("aa bb cc")
    universal_curve = spectral_curve(u, y, t, t**3)
    root_substitution = sp.Poly(
        sp.expand(universal_curve.subs(u, aa * y**2 + bb * y + cc)),
        y,
    )
    require(
        sp.expand(
            root_substitution.coeff_monomial(y**6)
            - infinity_polynomial.subs(v, aa)
        )
        == 0,
        "polynomial-root y^6 coefficient changed",
    )
    require(
        sp.expand(
            root_substitution.coeff_monomial(y**5)
            - bb * sp.diff(infinity_polynomial, v).subs(v, aa)
        )
        == 0,
        "polynomial-root y^5 coefficient changed",
    )
    require(
        sp.expand(
            root_substitution.coeff_monomial(y).subs(bb, 0)
            + 5878656 * t**3
        )
        == 0,
        "polynomial-root y coefficient changed",
    )
    require_global_irreducibility_certificate(
        universal_curve,
        quadratic_ratio,
        u,
        y,
        t,
        "quadratic",
    )
    require_global_irreducibility_certificate(
        universal_curve,
        sextic_ratio,
        u,
        y,
        t,
        "sextic",
    )

    quadratic = orbit_branch_data(
        quadratic_ratio,
        2,
        "smooth_triple",
        u,
        y,
        t,
    )
    sextic = orbit_branch_data(
        sextic_ratio,
        6,
        "ordinary_node",
        u,
        y,
        t,
    )
    require(
        quadratic["genus"] == 4 and sextic["genus"] == 3,
        "normalization genera changed",
    )
    require(
        quadratic["degree"] + sextic["degree"] == 8,
        "full B--W ratio count changed",
    )

    infinity_chart = sp.expand(
        z**6
        * spectral_curve(u, y, sp.Integer(1), sp.Integer(1)).subs(
            {u: v / z**2, y: 1 / z}
        )
    )
    require(
        sp.expand(infinity_chart.subs(z, 0) - infinity_polynomial) == 0,
        "weighted infinity chart changed",
    )

    hostile_curve = spectral_curve(
        u,
        y,
        sp.Integer(1),
        sp.Integer(1),
    )
    hostile_discriminant = sp.Poly(
        sp.discriminant(hostile_curve, u),
        y,
    )
    require(
        sp.gcd(
            hostile_discriminant,
            hostile_discriminant.diff(),
        ).degree()
        == 0,
        "B=W=1 hostile control unexpectedly has repeated branch",
    )

    print("ratio=B_W:W^2/B^5")
    print("representative=(B,W)=(alpha,alpha^3)")
    print("quadratic_ratio_orbit=degree_2_irreducible_mod_19")
    print("quadratic_branch=double_linear_times_squarefree_degree_10")
    print("quadratic_exceptional_fibre=smooth_total_ramification_e3")
    print(
        "quadratic_Gy_power_basis_vector="
        f"{quadratic['local_certificate']['derivative_y_vector']}"
    )
    print("quadratic_total_ramification=12")
    print("quadratic_normalization_genus=4")
    print("sextic_ratio_orbit=degree_6_irreducible_mod_13_Rabin")
    print("sextic_branch=double_linear_times_squarefree_degree_10")
    print("sextic_exceptional_fibre=ordinary_nonvertical_node")
    print(
        "sextic_tangent_u2_power_basis_vector="
        f"{sextic['local_certificate']['tangent_u2_vector']}"
    )
    print(
        "sextic_tangent_discriminant_power_basis_vector="
        f"{sextic['local_certificate']['tangent_discriminant_vector']}"
    )
    print("sextic_total_ramification=10")
    print("sextic_normalization_genus=3")
    print("absolute_irreducibility=PASS_coefficient_proof")
    print("absolute_irreducibility_control=PASS_two_global_Groebner_bases_[1]")
    print(f"infinity_cubic={infinity_polynomial}")
    print(f"infinity_discriminant={infinity_discriminant}")
    print("infinity=3_distinct_unramified_points")
    print("closed_BW_ratios=8")
    print("remaining_exactly_two_sparse_ratios=0")
    print("hostile_control_B=W=1_squarefree=PASS")
    print("riemann_hurwitz_and_deck_contradiction=MATHEMATICAL_PROOF_REQUIRED")
    print("status=THM2328_FULL_BW_BANK_EXACT_REFEREE")


if __name__ == "__main__":
    main()
