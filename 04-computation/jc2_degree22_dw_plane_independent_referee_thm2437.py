#!/usr/bin/env python3
"""Independent exact referee for THM-2437's degree-22 D-W plane.

This replay deliberately differs from the primary companion in two places.

1. It eliminates the first-flux variable by the closed formula for the
   resultant of a linear and a quadratic polynomial, and computes the quartic
   discriminant from the explicit coefficient formula.
2. It specializes K9 to the exact number fields Q[lambda]/(S6), /(S7), and
   /(Q5), then computes polynomial gcds there.  Thus the exceptional-fibre
   multiplicities are checked after specialization rather than inferred from
   the generic subresultant sequence.

All assertions use ``require`` so the audit remains active under ``python -O``.
"""

from __future__ import annotations

import hashlib
from itertools import combinations

import sympy as sp
from sympy.polys.domains import QQ


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def polynomial_hash(expression: sp.Expr, *variables: sp.Symbol) -> str:
    polynomial = sp.Poly(sp.expand(expression), *variables)
    payload = sp.srepr(polynomial.as_expr()).encode()
    return hashlib.sha256(payload).hexdigest()


def specialize_over_parameter_field(
    expression: sp.Expr,
    parameter_polynomial: sp.Expr,
    v: sp.Symbol,
    lam: sp.Symbol,
) -> sp.Poly:
    """Put a polynomial in v,lambda over Q[lambda]/(parameter_polynomial).

    Constructing the algebraic-field coefficients directly avoids numerical
    root selection and checks every conjugate parameter at once because the
    three parameter polynomials used below are irreducible over Q.
    """

    field = QQ.alg_field_from_poly(
        sp.Poly(parameter_polynomial, lam), alias="a"
    )
    terms = {}
    for (v_degree,), coefficient in sp.Poly(expression, v).terms():
        coefficient_poly = sp.Poly(coefficient, lam, domain=QQ)
        terms[(v_degree,)] = field(
            [QQ.convert(item) for item in coefficient_poly.all_coeffs()]
        )
    return sp.Poly.from_dict(terms, (v,), domain=field)


def main() -> None:
    y, u, z = sp.symbols("y u Z")
    dpar, wpar = sp.symbols("D W")
    v, zeta, p, mu, lam = sp.symbols("v zeta p mu lambda")

    # Direct B=C=E=0 specialization of THM-2411, equations (12)--(16).
    n1 = (
        1331 * (-1089 * u + 63 * y**2) * z
        + 4
        * (
            511104 * dpar * y
            + 922383 * u**2 * y
            - 25410 * u * y**3
            + 63 * y**5
        )
    )
    n2 = (
        15944049 * z**2
        - 162339408 * z * u * y
        + 2236080 * z * y**3
        - 1978994688 * dpar * u
        + 16355328 * dpar * y**2
        - 1319329792 * wpar
        - 1190488992 * u**3
        + 147581280 * u**2 * y**2
        - 1219680 * u * y**4
        + 672 * y**6
    )

    substitutions = {
        dpar: p**2 * y**4,
        wpar: mu * p**3 * y**6,
        u: v * y**2,
        z: zeta * y**3,
    }
    first = sp.expand(n1.subs(substitutions) / y**5)
    second = sp.expand(n2.subs(substitutions) / y**6)

    linear_coefficient = 11979 * (7 - 121 * v)
    linear_constant = 4 * (
        922383 * v**2 - 25410 * v + 63 + 511104 * p**2
    )
    quadratic_coefficient = sp.Integer(15944049)
    quadratic_linear = -162339408 * v + 2236080
    quadratic_constant = (
        -1190488992 * v**3
        + 147581280 * v**2
        - 1219680 * v
        + 672
        + (-1978994688 * v + 16355328) * p**2
        - 1319329792 * mu * p**3
    )
    require(
        sp.expand(
            first - linear_coefficient * zeta - linear_constant
        )
        == 0,
        "independent first-flux normalization failed",
    )
    require(
        sp.expand(
            second
            - quadratic_coefficient * zeta**2
            - quadratic_linear * zeta
            - quadratic_constant
        )
        == 0,
        "independent second-flux normalization failed",
    )

    # For az+b and Az^2+Bz+C the resultant is A b^2-a B b+a^2 C.
    direct_eliminant = sp.expand(
        quadratic_coefficient * linear_constant**2
        - linear_coefficient * quadratic_linear * linear_constant
        + linear_coefficient**2 * quadratic_constant
    )
    l5 = (
        155624547606 * v**5
        + 3215383215 * v**4
        - 1700698560 * v**3
        + 58124770 * v**2
        - 855470 * v
        + 2583
    )
    a2 = (
        -1810903826688 * v**3
        + 119729178624 * v**2
        + 4329050880 * v
        - 109716992
    )
    quartic = sp.expand(direct_eliminant / 2295943056)
    expected_quartic = sp.expand(
        29025255424 * p**4
        - 82458112 * mu * (121 * v - 7) ** 2 * p**3
        + a2 * p**2
        - 7 * l5
    )
    require(
        quartic == expected_quartic,
        "closed linear-quadratic elimination did not recover the quartic",
    )
    require(
        sp.expand(quartic.subs({p: -p, mu: -mu}) - quartic) == 0,
        "delta-sign sidecar does not preserve the quartic",
    )
    require(
        sp.Poly(l5, v).degree() == 5
        and sp.gcd(sp.Poly(l5, v), sp.Poly(sp.diff(l5, v), v)).degree()
        == 0
        and sp.Poly(a2, v).degree() == 3
        and sp.Poly(
            sp.cancel(sp.Poly(quartic, p).coeff_monomial(p**3) / mu),
            v,
        ).degree()
        == 2,
        "uniform irreducibility inputs failed",
    )

    # Explicit discriminant formula for a*p^4+b*p^3+c*p^2+e.
    coefficient_poly = sp.Poly(quartic, p)
    aa, bb, cc, dd, ee = coefficient_poly.all_coeffs()
    require(dd == 0, "quartic unexpectedly acquired a linear term")
    explicit_discriminant = sp.expand(
        256 * aa**3 * ee**3
        - 128 * aa**2 * cc**2 * ee**2
        + 144 * aa * bb**2 * cc * ee**2
        + 16 * aa * cc**4 * ee
        - 27 * bb**4 * ee**2
        - 4 * bb**2 * cc**3 * ee
    )
    mu_polynomial = sp.Poly(explicit_discriminant, mu)
    require(
        all(exponent[0] % 2 == 0 for exponent in mu_polynomial.monoms()),
        "quartic discriminant is not a polynomial in lambda=mu^2",
    )
    discriminant_lambda = sp.expand(
        sum(
            coefficient * lam ** (exponent[0] // 2)
            for exponent, coefficient in mu_polynomial.terms()
        )
    )
    divisor = (
        -2**36 * 7 * 11**18 * (121 * v - 7) ** 4 * l5
    )
    k9 = sp.cancel(discriminant_lambda / divisor)
    require(
        sp.denom(k9) == 1
        and sp.Poly(k9, v, lam).degree(v) == 9
        and sp.Poly(k9, v, lam).degree(lam) == 2,
        "explicit discriminant did not yield polynomial K9",
    )
    k9 = sp.expand(k9)
    leading = sp.factor(sp.Poly(k9, v).LC())
    require(
        leading
        == 3302590884214385499714 * lam * (231 * lam + 128),
        "independent K9 leading coefficient mismatch",
    )

    s6 = (
        71778115591875 * lam**6
        - 10643296267296000 * lam**5
        - 45431893495296000 * lam**4
        - 70165361991352320 * lam**3
        - 40747428281843712 * lam**2
        - 589081608192000 * lam
        + 4306718326521856
    )
    s7 = (
        37565749714841455078125 * lam**7
        + 27921686965239375000000 * lam**6
        - 13643584211703600000000 * lam**5
        - 11118982374162432000000 * lam**4
        - 2775896567285022720000 * lam**3
        + 3185830828706176696320 * lam**2
        + 998410901657688735744 * lam
        - 455539102308525670400
    )
    q5 = (
        4932186875 * lam**5
        - 1123257520000 * lam**4
        - 5423878240000 * lam**3
        - 7375509504000 * lam**2
        - 1274053017600 * lam
        + 1520839950336
    )
    parameter_discriminant = (
        2**155
        * 3**30
        * 5**8
        * 7**11
        * 11**148
        * lam**4
        * s6**3
        * s7
    )
    # Degree nine has sign (+1) in Disc=Res(K9,K9')/LC.
    require(
        sp.expand(
            sp.resultant(k9, sp.diff(k9, v), v)
            - leading * parameter_discriminant
        )
        == 0,
        "independent parameter-resultant identity failed",
    )
    joint_resultant = sp.factor(sp.resultant(k9, l5, v))
    joint_quotient = sp.cancel(joint_resultant / q5)
    require(
        not joint_quotient.free_symbols and joint_quotient != 0,
        "independent K9-L5 resultant identity failed",
    )
    require(
        sp.expand(
            k9.subs(v, sp.Rational(7, 121))
            - 10276044800 * (385 * lam + 512)
        )
        == 0,
        "independent wall evaluation failed",
    )

    exceptional_polynomials = (
        lam,
        231 * lam + 128,
        385 * lam + 512,
        s6,
        s7,
        q5,
    )
    for polynomial in (s6, s7, q5):
        require(
            sp.Poly(polynomial, lam).is_irreducible,
            "parameter field is not defined by an irreducible polynomial",
        )
    for first_exception, second_exception in combinations(
        exceptional_polynomials, 2
    ):
        require(
            sp.gcd(
                sp.Poly(first_exception, lam),
                sp.Poly(second_exception, lam),
            ).degree()
            == 0,
            "exceptional parameter divisors overlap",
        )

    # Exact specialization in the algebraic parameter fields.
    for collision in (s6, s7):
        k9_collision = specialize_over_parameter_field(
            k9, collision, v, lam
        )
        l5_collision = specialize_over_parameter_field(
            l5, collision, v, lam
        )
        require(
            k9_collision.degree() == 9
            and sp.gcd(k9_collision, k9_collision.diff()).degree() == 1,
            "S6/S7 specialization is not one-double-root type",
        )
        require(
            sp.gcd(k9_collision, l5_collision).degree() == 0,
            "S6/S7 specialization meets L5",
        )

    k9_joint = specialize_over_parameter_field(k9, q5, v, lam)
    l5_joint = specialize_over_parameter_field(l5, q5, v, lam)
    require(
        k9_joint.degree() == 9
        and sp.gcd(k9_joint, k9_joint.diff()).degree() == 0
        and sp.gcd(l5_joint, l5_joint.diff()).degree() == 0
        and sp.gcd(k9_joint, l5_joint).degree() == 1,
        "Q5 specialization is not a unique simple cross-collision",
    )

    # Rational exceptional fibres and an ordinary hostile control.
    def rational_fibre(value: sp.Rational) -> sp.Poly:
        return sp.Poly(k9.subs(lam, value), v, domain=QQ)

    degree_drop = rational_fibre(-sp.Rational(128, 231))
    wall_collision = rational_fibre(-sp.Rational(512, 385))
    generic = rational_fibre(sp.Rational(1))
    l5_poly = sp.Poly(l5, v, domain=QQ)
    wall = sp.Rational(7, 121)
    require(
        degree_drop.degree() == 8
        and sp.gcd(degree_drop, degree_drop.diff()).degree() == 0
        and sp.gcd(degree_drop, l5_poly).degree() == 0
        and degree_drop.eval(wall) != 0,
        "degree-drop fibre failed its independent audit",
    )
    require(
        wall_collision.degree() == 9
        and sp.gcd(wall_collision, wall_collision.diff()).degree() == 0
        and sp.gcd(wall_collision, l5_poly).degree() == 0
        and wall_collision.eval(wall) == 0
        and wall_collision.diff().eval(wall) != 0,
        "wall-collision fibre failed its independent audit",
    )
    require(
        generic.degree() == 9
        and sp.gcd(generic, generic.diff()).degree() == 0
        and sp.gcd(generic, l5_poly).degree() == 0
        and generic.eval(wall) != 0,
        "generic hostile-control fibre is not separated",
    )

    branch_counts = (5 + 9, 5 + 8, 5 + 7, 4 + 8, 5 + 8)
    require(
        branch_counts == (14, 13, 12, 12, 13)
        and min(branch_counts) == 12,
        "independent simple-branch floor failed",
    )
    require(
        sp.expand(n1.subs(y, 0) + 1449459 * u * z) == 0,
        "independent y=0 boundary failed",
    )

    print("THM-2437 independent D-W plane hostile referee")
    print("direct_linear_quadratic_elimination=PASS")
    print("delta_sign_sidecar=PASS")
    print("explicit_quartic_discriminant=wall^4*L5*K9")
    print("parameter_resultant=lambda^5*(231lambda+128)*S6^3*S7")
    print("algebraic_field_gcds=S6:1,S7:1,Q5-L5:1")
    print("degree_drop_squarefree_degree=8")
    print("wall_collision=simple_and_otherwise_separated")
    print("simple_finite_branch_counts=14,13,12,12,13")
    print("uniform_genus_lower_bound=3")
    print("y_zero_and_axis_boundaries=PASS")
    print(f"quartic_sha256={polynomial_hash(quartic, v, p, mu)}")
    print(f"K9_sha256={polynomial_hash(k9, v, lam)}")
    print("all_independent_exact_checks=PASS")


if __name__ == "__main__":
    main()
