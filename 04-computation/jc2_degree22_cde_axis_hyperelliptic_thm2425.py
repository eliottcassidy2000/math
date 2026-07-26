#!/usr/bin/env python3
"""Exact referee for the degree-22 C, D, and E one-sparse axes.

Starting from the hostile-audited normalized fluxes of THM-2411, this
companion independently performs the weighted substitutions

    v=u/y^2,  zeta=Z/y^3,  p=X/y^wt(X),  X in {C,D,E},

eliminates zeta, and verifies that each remaining quadratic in p completes
to a squarefree hyperelliptic curve.  The rational-map and constant-field
arguments are the mathematical part of THM-2425.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def polynomial_hash(expression: sp.Expr, variable: sp.Symbol) -> str:
    polynomial = sp.Poly(sp.expand(expression), variable)
    return hashlib.sha256(sp.srepr(polynomial.as_expr()).encode()).hexdigest()


def main() -> None:
    y, u, z = sp.symbols("y u Z")
    b, cpar, dpar, epar, wpar = sp.symbols("B C D E W")
    v, zeta, p = sp.symbols("v zeta p")

    pole_a = 616 * b - 1089 * u + 63 * y**2
    pole_k = (
        -745360 * b * u * y
        + 6160 * b * y**3
        + 2342560 * cpar * u
        - 58080 * cpar * y**2
        + 511104 * dpar * y
        - 3748096 * epar
        + 922383 * u**2 * y
        - 25410 * u * y**3
        + 63 * y**5
    )
    n1 = 1331 * pole_a * z + 4 * pole_k
    n2 = (
        15944049 * z**2
        + 65591680 * b * z * y
        - 206145280 * cpar * z
        - 162339408 * z * u * y
        + 2236080 * z * y**3
        + 1443016960 * b * u**2
        - 71554560 * b * u * y**2
        + 98560 * b * y**4
        + 449771520 * cpar * u * y
        - 1239040 * cpar * y**3
        - 1978994688 * dpar * u
        + 16355328 * dpar * y**2
        - 239878144 * epar * y
        - 1319329792 * wpar
        - 1190488992 * u**3
        + 147581280 * u**2 * y**2
        - 1219680 * u * y**4
        + 672 * y**6
    )

    base_k = 922383 * v**2 - 25410 * v + 63
    base_g = (
        15944049 * zeta**2
        - 162339408 * zeta * v
        + 2236080 * zeta
        - 1190488992 * v**3
        + 147581280 * v**2
        - 1219680 * v
        + 672
    )
    scaled_a = 9 * (7 - 121 * v)

    r_c = (
        (-5487587353600 * v**2 + 634927462400 * v - 12368716800)
        * p**2
        + (
            -4938828618240 * v**3
            + 537420744960 * v**2
            - 12593602560 * v
            + 146361600
        )
        * p
        - 9804346499178 * v**5
        - 202569142545 * v**4
        + 107144009280 * v**3
        - 3661860510 * v**2
        + 53894610 * v
        - 162729
    )
    h_c = (
        363123944414 * v**5
        - 33654344317 * v**4
        - 170069856 * v**3
        + 7408346 * v**2
        + 72842 * v
        - 6741
    )

    r_d = (
        29025255424 * p**2
        + (
            -1810903826688 * v**3
            + 119729178624 * v**2
            + 4329050880 * v
            - 109716992
        )
        * p
        - 1089371833242 * v**5
        - 22507682505 * v**4
        + 11904889920 * v**3
        - 406873390 * v**2
        + 5988290 * v
        - 18081
    )
    h_d = (
        1929229929 * v**4
        + 42517464 * v**3
        - 790614 * v**2
        - 203280 * v
        + 2485
    )

    r_e = (
        14048223625216 * p**2
        + (4938828618240 * v**2 - 571434716160 * v + 3935500800)
        * p
        - 9804346499178 * v**5
        - 202569142545 * v**4
        + 107144009280 * v**3
        - 3661860510 * v**2
        + 53894610 * v
        - 162729
    )
    h_e = 37202781 * v**3 + 6720219 * v**2 - 134673 * v + 497

    axis_data = (
        (
            "C",
            cpar,
            3,
            base_k + (2342560 * v - 58080) * p,
            base_g
            + (-206145280 * zeta + 449771520 * v - 1239040) * p,
            r_c,
            255104784,
            h_c,
            -40479436800,
            5,
            3483463546133440530410221219484583413585848637590432634557735895040,
        ),
        (
            "D",
            dpar,
            4,
            base_k + 511104 * p,
            base_g + (-1978994688 * v + 16355328) * p,
            r_d,
            2295943056,
            h_d,
            116101021696,
            4,
            -140489835872244482219052102409329102028800,
        ),
        (
            "E",
            epar,
            5,
            base_k - 3748096 * p,
            base_g - 239878144 * p,
            r_e,
            255104784,
            h_e,
            1011472101015552,
            3,
            268774001885305889280000,
        ),
    )

    zero_parameters = {b: 0, cpar: 0, dpar: 0, epar: 0, wpar: 0}
    summaries: list[tuple[str, int, int, str]] = []

    for (
        name,
        parameter,
        weight,
        scaled_k,
        scaled_g,
        resultant_polynomial,
        resultant_scalar,
        hyperelliptic_polynomial,
        discriminant_scalar,
        degree,
        hyperelliptic_discriminant,
    ) in axis_data:
        substitutions = dict(zero_parameters)
        substitutions[parameter] = p * y**weight
        substitutions.update({u: v * y**2, z: zeta * y**3})

        require(
            sp.factor(n1.subs(substitutions) / y**5)
            == sp.factor(1331 * scaled_a * zeta + 4 * scaled_k),
            f"{name}-axis first weighted flux mismatch",
        )
        require(
            sp.factor(n2.subs(substitutions) / y**6 - scaled_g) == 0,
            f"{name}-axis second weighted flux mismatch",
        )

        f1 = sp.expand(1331 * scaled_a * zeta + 4 * scaled_k)
        exact_resultant = sp.factor(sp.resultant(f1, scaled_g, zeta))
        require(
            sp.factor(
                exact_resultant - resultant_scalar * resultant_polynomial
            )
            == 0,
            f"{name}-axis resultant mismatch",
        )

        coefficients = sp.Poly(resultant_polynomial, p).all_coeffs()
        require(
            len(coefficients) == 3,
            f"{name}-axis resultant is not quadratic in p",
        )
        require(
            sp.gcd_list(coefficients) == 1,
            f"{name}-axis has a vertical identity fibre",
        )

        quadratic_discriminant = sp.factor(
            sp.discriminant(resultant_polynomial, p)
        )
        require(
            quadratic_discriminant
            == discriminant_scalar
            * (121 * v - 7) ** 2
            * hyperelliptic_polynomial,
            f"{name}-axis completed-square factorization mismatch",
        )
        require(
            sp.degree(hyperelliptic_polynomial, v) == degree,
            f"{name}-axis residual degree mismatch",
        )
        require(
            sp.gcd(
                sp.Poly(hyperelliptic_polynomial, v),
                sp.Poly(sp.diff(hyperelliptic_polynomial, v), v),
            ).degree()
            == 0,
            f"{name}-axis hyperelliptic polynomial is not squarefree",
        )
        require(
            sp.discriminant(hyperelliptic_polynomial, v)
            == hyperelliptic_discriminant,
            f"{name}-axis hyperelliptic discriminant mismatch",
        )
        require(
            sp.factor(
                hyperelliptic_polynomial.subs(v, sp.Rational(7, 121))
            )
            != 0,
            f"{name}-axis residual meets the excluded first-flux wall",
        )

        genus = (degree - 1) // 2
        summaries.append(
            (
                name,
                degree,
                genus,
                polynomial_hash(hyperelliptic_polynomial, v),
            )
        )

    # Exact y=0 hostile boundary checks.
    c_zero_z = sp.Rational(640, 99) * cpar
    require(
        sp.factor(
            n1.subs(
                {
                    y: 0,
                    b: 0,
                    dpar: 0,
                    epar: 0,
                    wpar: 0,
                    z: c_zero_z,
                }
            )
        )
        == 0,
        "C-axis y=0 first-flux reconstruction failed",
    )
    require(
        sp.factor(
            n2.subs(
                {
                    y: 0,
                    b: 0,
                    dpar: 0,
                    epar: 0,
                    wpar: 0,
                    z: c_zero_z,
                }
            )
            + sp.Rational(468512, 9)
            * (12800 * cpar**2 + 22869 * u**3)
        )
        == 0,
        "C-axis y=0 constant-field equation failed",
    )

    require(
        sp.factor(
            n1.subs(
                {
                    y: 0,
                    b: 0,
                    cpar: 0,
                    epar: 0,
                    wpar: 0,
                }
            )
            + 1331 * 1089 * u * z
        )
        == 0,
        "D-axis y=0 first flux should force Z=0",
    )

    e_zero_z = -sp.Rational(1024, 99) * epar / u
    require(
        sp.factor(
            n1.subs(
                {
                    y: 0,
                    b: 0,
                    cpar: 0,
                    dpar: 0,
                    wpar: 0,
                    z: e_zero_z,
                }
            )
        )
        == 0,
        "E-axis y=0 first-flux reconstruction failed",
    )
    require(
        sp.factor(
            n2.subs(
                {
                    y: 0,
                    b: 0,
                    cpar: 0,
                    dpar: 0,
                    wpar: 0,
                    z: e_zero_z,
                }
            )
            + sp.Rational(468512, 9) / u**2
            * (-32768 * epar**2 + 22869 * u**5)
        )
        == 0,
        "E-axis y=0 constant-field equation failed",
    )

    print("THM-2425 degree-22 C/D/E-axis exact referee")
    print("weighted_flux_reductions=PASS")
    for name, degree, genus, digest in summaries:
        print(
            f"{name}_axis=residual_degree_{degree},"
            f"genus_{genus},squarefree=PASS,sha256={digest}"
        )
    print("quadratic_coefficient_gcds=1")
    print("excluded_wall_controls=PASS")
    print("y_zero_hostile_controls=PASS")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
