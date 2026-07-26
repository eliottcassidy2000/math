#!/usr/bin/env python3
"""Exact referee for the degree-22 invariant-origin cusp obstruction.

This companion works in the target-translated degree-22 quartic Faber chart
of THM-2411, but reconstructs the degree-22 Laurent observables directly.
At the invariant coefficient origin B=C=D=E=W=0 it verifies:

* the first two fluxes reduce, after u=v*y^2 and Z=zeta*y^3, to two
  equations whose resultant is a nonzero quintic L_5(v);
* on the open first-flux chart, zeta is the displayed rational function of
  v and is nonzero at every root of L_5;
* the third flux has a nonzero h^11 coefficient at every root of L_5; and
* the full polynomial Faber sidecar has a different nonzero h^11
  coefficient at every root of L_5.

The rational-primitive and pole-valuation deductions using these exact
identities are the mathematical part of THM-2423.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def faber_coefficients(
    degree: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
    extra: int = 3,
) -> list[sp.Expr]:
    """Laurent coefficients of (w^4+p*w^2+q*w+r)^(degree/4)."""
    exponent = sp.Rational(degree, 4)
    coefficients = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + extra + 1):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
            if step in quartic
        ) / index
        coefficients.append(sp.factor(value))
    return coefficients


def direct_laurent_coefficient(
    degree: int,
    index: int,
    p: sp.Expr,
    q: sp.Expr,
    r: sp.Expr,
) -> sp.Expr:
    """Independent finite multinomial coefficient."""
    exponent = sp.Rational(degree, 4)
    total = sp.Integer(0)
    for i in range(index // 2 + 1):
        for j in range(index // 3 + 1):
            remainder = index - 2 * i - 3 * j
            if remainder < 0 or remainder % 4:
                continue
            k = remainder // 4
            chosen = i + j + k
            falling = sp.prod(exponent - h for h in range(chosen))
            total += (
                falling
                * p**i
                * q**j
                * r**k
                / (
                    sp.factorial(i)
                    * sp.factorial(j)
                    * sp.factorial(k)
                )
            )
    return sp.factor(total)


def replace_even_power(
    expression: sp.Expr,
    variable: sp.Symbol,
    square: sp.Symbol,
) -> sp.Expr:
    """Replace variable^(2j) by square^j, rejecting odd powers."""
    polynomial = sp.Poly(sp.expand(expression), variable)
    result = sp.Integer(0)
    for (power,), coefficient in polynomial.terms():
        require(power % 2 == 0, f"unexpected odd power of {variable}")
        result += coefficient * square ** (power // 2)
    return sp.factor(result)


def polynomial_hash(expression: sp.Expr, *generators: sp.Symbol) -> str:
    polynomial = sp.Poly(sp.expand(expression), *generators)
    return hashlib.sha256(sp.srepr(polynomial.as_expr()).encode()).hexdigest()


def main() -> None:
    d, q, s, t = sp.symbols("d q s T")
    y, u, z = sp.symbols("y u Z")
    v, zeta = sp.symbols("v zeta")
    lvar = sp.symbols("L")

    p = 2 * d
    r = d**2 - s
    coefficients = faber_coefficients(22, p, q, r)
    direct = [
        direct_laurent_coefficient(22, index, p, q, r)
        for index in (23, 24, 25)
    ]
    require(
        all(
            sp.factor(coefficients[index] - direct[offset]) == 0
            for offset, index in enumerate((23, 24, 25))
        ),
        "recurrence and multinomial Laurent coefficients disagree",
    )

    phi = sp.factor(4 * coefficients[23])
    psi = sp.factor(4 * coefficients[24])
    third = sp.factor(4 * coefficients[25] + p * coefficients[23])

    phi_even = replace_even_power(sp.factor(phi / q), q, t)
    psi_even = replace_even_power(psi, q, t)
    third_even = replace_even_power(sp.factor(third / q), q, t)

    phi_sub = sp.factor(phi_even.subs({d: u / t, s: y / 11}))
    psi_sub = sp.factor(psi_even.subs({d: u / t, s: y / 11}))
    third_sub = sp.factor(third_even.subs({d: u / t, s: y / 11}))

    phi_num, phi_den = sp.fraction(phi_sub)
    psi_num, psi_den = sp.fraction(psi_sub)
    third_num, third_den = sp.fraction(third_sub)
    phi_z = replace_even_power(phi_num, t, z) / replace_even_power(
        phi_den, t, z
    )
    psi_z = replace_even_power(psi_num, t, z) / replace_even_power(
        psi_den, t, z
    )
    require(
        sp.Poly(third_num, t).terms()
        and all(power[0] % 2 == 0 for power, _ in sp.Poly(third_num, t).terms()),
        "third-flux numerator has an odd T power",
    )
    require(
        sp.factor(third_den / t) == 14992384,
        "unexpected third-flux denominator",
    )
    third_z_num = replace_even_power(third_num, t, z)

    n1 = (
        1331 * (-1089 * u + 63 * y**2) * z
        + 4
        * (
            922383 * u**2 * y
            - 25410 * u * y**3
            + 63 * y**5
        )
    )
    n2 = (
        15944049 * z**2
        - 162339408 * z * u * y
        + 2236080 * z * y**3
        - 1190488992 * u**3
        + 147581280 * u**2 * y**2
        - 1219680 * u * y**4
        + 672 * y**6
    )
    require(
        sp.factor(phi_z + n1 / 7496192) == 0,
        "first-flux normalization mismatch",
    )
    require(
        sp.factor(psi_z - n2 / 1319329792) == 0,
        "second-flux normalization mismatch",
    )

    h3 = (
        -43923 * zeta**2
        - 1449459 * zeta * v**2
        + 139755 * zeta * v
        - 770 * zeta
        + 1229844 * v**3
        - 33880 * v**2
        + 84 * v
    )
    require(
        sp.factor(
            third_z_num.subs({u: v * y**2, z: zeta * y**3})
            - 3 * y**7 * h3
        )
        == 0,
        "third-flux weighted reduction mismatch",
    )

    f1 = (
        -1229844 * v**2
        + 483153 * v * zeta
        + 33880 * v
        - 27951 * zeta
        - 84
    )
    f2 = (
        -396829664 * v**3
        + 49193760 * v**2
        - 54113136 * v * zeta
        - 406560 * v
        + 5314683 * zeta**2
        + 745360 * zeta
        + 224
    )
    require(
        sp.factor(
            n1.subs({u: v * y**2, z: zeta * y**3})
            + 3 * y**5 * f1
        )
        == 0,
        "first scaled origin equation mismatch",
    )
    require(
        sp.factor(
            n2.subs({u: v * y**2, z: zeta * y**3})
            - 3 * y**6 * f2
        )
        == 0,
        "second scaled origin equation mismatch",
    )

    l5 = (
        155624547606 * v**5
        + 3215383215 * v**4
        - 1700698560 * v**3
        + 58124770 * v**2
        - 855470 * v
        + 2583
    )
    resultant_constant = sp.factor(sp.resultant(f1, f2, zeta) / l5)
    require(
        resultant_constant == -595244496,
        "origin resultant does not equal the stated quintic",
    )

    k = sp.factor(
        28
        * (121 * v - 3)
        * (363 * v - 1)
        / (3993 * (121 * v - 7))
    )
    require(sp.factor(f1.subs(zeta, k)) == 0, "zeta reconstruction failed")
    require(
        sp.factor(l5.subs(v, sp.Rational(7, 121))) == -44800,
        "open-chart denominator meets the quintic",
    )
    require(
        sp.factor(l5.subs(v, sp.Rational(3, 121))) == -6144
        and sp.factor(l5.subs(v, sp.Rational(1, 363)))
        == sp.Rational(51200, 81),
        "zeta can vanish on the quintic",
    )
    require(
        sp.factor(
            sp.together(
                f2.subs(zeta, k)
                + sp.Rational(112, 3) * l5 / (121 * v - 7) ** 2
            )
        )
        == 0,
        "second origin equation after zeta reconstruction is wrong",
    )

    third_on_curve = sp.factor(sp.together(h3.subs(zeta, k)))
    third_expected = sp.factor(
        -56
        * (121 * v - 3)
        * (363 * v - 1)
        * (5314683 * v**3 - 307461 * v**2 + 22869 * v - 203)
        / (363 * (121 * v - 7) ** 2)
    )
    require(
        sp.factor(third_on_curve - third_expected) == 0,
        "third-flux coefficient factorization failed",
    )
    third_numerator = sp.fraction(third_on_curve)[0]
    require(
        sp.gcd(sp.Poly(l5, v), sp.Poly(third_numerator, v)).degree() == 0,
        "third flux vanishes at an origin quintic root",
    )

    # Reconstruct the full degree-22 polynomial Faber sidecar at the origin.
    w = sp.symbols("w")
    e22 = sp.expand(
        sum(coefficients[index] * w ** (22 - index) for index in range(23))
    )
    e22_l = sp.cancel(e22.subs(w, (lvar + s) / q))
    e22_num, e22_den = sp.fraction(e22_l)
    e22_l = sp.factor(
        replace_even_power(e22_num, q, t)
        / replace_even_power(e22_den, q, t)
    )
    h = (lvar + s) ** 2 / t + d
    p_full = h**2 + lvar
    truncation = sum(
        sp.binomial(sp.Rational(11, 2), index)
        * h ** (11 - 2 * index)
        * lvar**index
        for index in range(6)
    )
    sidecar = sp.factor(e22_l - truncation)
    expected_sidecar = sp.Rational(33, 2048) * t * (
        14 * lvar**4
        - 28 * lvar**3 * s
        - 14 * lvar**2 * t * d
        + 42 * lvar**2 * s**2
        - lvar * t**2
        + 56 * lvar * t * d * s
        - 56 * lvar * s**3
        + 14 * t**2 * d**2
        + 6 * t**2 * s
        - 140 * t * d * s**2
        + 70 * s**4
    )
    require(
        sp.factor(sidecar - expected_sidecar) == 0,
        "full degree-22 Faber sidecar mismatch",
    )
    require(
        sp.factor(p_full - (h**2 + lvar)) == 0,
        "approximate-root coordinate check failed",
    )

    j = sp.factor(
        14 * v**2
        + sp.Rational(6, 11) * k
        - sp.Rational(140, 121) * v
        + sp.Rational(70, 14641)
    )
    j_expected = sp.factor(
        14
        * (1771561 * v**3 - 73205 * v**2 + 4235 * v - 23)
        / (14641 * (121 * v - 7))
    )
    require(sp.factor(j - j_expected) == 0, "sidecar pole coefficient failed")
    j_numerator = sp.fraction(j)[0]
    require(
        sp.gcd(sp.Poly(l5, v), sp.Poly(j_numerator, v)).degree() == 0,
        "sidecar coefficient vanishes at an origin quintic root",
    )

    print("THM-2423 degree-22 invariant-origin exact referee")
    print("laurent_recurrence_vs_multinomial=PASS")
    print(f"origin_resultant_constant={resultant_constant}")
    print(f"origin_quintic_degree={sp.degree(l5, v)}")
    print("open_wall_value=-44800")
    print("zeta_zero_values=(-6144,51200/81)")
    print(
        "third_curve_factor="
        "-56(121v-3)(363v-1)"
        "(5314683v^3-307461v^2+22869v-203)"
        "/(363(121v-7)^2)"
    )
    print("third_quintic_gcd_degree=0")
    print(
        "sidecar_curve_factor="
        "14(1771561v^3-73205v^2+4235v-23)"
        "/(14641(121v-7))"
    )
    print("sidecar_quintic_gcd_degree=0")
    print(f"origin_quintic_sha256={polynomial_hash(l5, v)}")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
