#!/usr/bin/env python3
"""Exact algebra in progress for the degree-22 D-W coefficient plane."""

from __future__ import annotations

import sympy as sp


def main() -> None:
    v, zeta, p, mu, lam = sp.symbols("v zeta p mu lambda")

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
    f1 = 11979 * (7 - 121 * v) * zeta + 4 * (
        base_k + 511104 * p**2
    )
    f2 = (
        base_g
        + (-1978994688 * v + 16355328) * p**2
        - 1319329792 * mu * p**3
    )

    resultant = sp.factor(sp.resultant(f1, f2, zeta))
    content, primitive = sp.Poly(resultant, v, p, mu).primitive()
    print(f"resultant_content={content}")
    print(f"R={sp.factor(primitive.as_expr())}")
    print(f"R_degree_p={sp.degree(primitive.as_expr(), p)}")

    disc = sp.factor(sp.discriminant(primitive.as_expr(), p))
    print(f"Disc_p_factor={disc}")
    disc_primitive = sp.factor(
        disc / (121 * v - 7) ** 4
    )
    print(f"Disc_without_wall={disc_primitive}")

    # Separate the known W-axis quintic L5 and display the remaining factor.
    l5 = (
        155624547606 * v**5
        + 3215383215 * v**4
        - 1700698560 * v**3
        + 58124770 * v**2
        - 855470 * v
        + 2583
    )
    quotient = sp.factor(disc_primitive / l5)
    q_content, q_primitive = sp.Poly(quotient, v, mu).primitive()
    print(f"disc_residual_content={q_content}")
    print(f"K9_mu={q_primitive.as_expr()}")
    require_even = sp.Poly(q_primitive.as_expr(), mu)
    assert all(power[0] % 2 == 0 for power, _ in require_even.terms())
    k9_lam = sp.expand(q_primitive.as_expr().subs(mu**2, lam))
    print(f"K9_lambda={k9_lam}")
    print(f"K9_degree_v={sp.degree(k9_lam, v)}")
    print(f"K9_degree_lambda={sp.degree(k9_lam, lam)}")
    print(f"K9_lc_v={sp.factor(sp.Poly(k9_lam, v).LC())}")
    print(
        "K9_at_wall="
        f"{sp.factor(k9_lam.subs(v, sp.Rational(7, 121)))}"
    )
    print(
        "L5_at_wall="
        f"{sp.factor(l5.subs(v, sp.Rational(7, 121)))}"
    )
    disc_k9 = sp.factor(sp.discriminant(k9_lam, v))
    res_k9_l5 = sp.factor(sp.resultant(k9_lam, l5, v))
    print(f"Disc_v_K9={disc_k9}")
    print(f"Res_K9_L5={res_k9_l5}")

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
    t5 = (
        4932186875 * lam**5
        - 1123257520000 * lam**4
        - 5423878240000 * lam**3
        - 7375509504000 * lam**2
        - 1274053017600 * lam
        + 1520839950336
    )
    for name, exceptional in (("S6", s6), ("S7", s7), ("T5", t5)):
        print(f"{name}_factor={sp.factor(exceptional)}")
        print(
            f"{name}_squarefree="
            f"{sp.gcd(sp.Poly(exceptional, lam), sp.Poly(sp.diff(exceptional, lam), lam)).degree() == 0}"
        )
    print(f"K9_lambda0_factor={sp.factor(k9_lam.subs(lam, 0))}")
    drop = sp.factor(k9_lam.subs(lam, -sp.Rational(128, 231)))
    print(f"K9_drop_factor={drop}")
    print(
        "K9_drop_degree_squarefree="
        f"{sp.degree(drop, v)},"
        f"{sp.gcd(sp.Poly(drop, v), sp.Poly(sp.diff(drop, v), v)).degree()}"
    )

    def vanishes_mod_parameter(
        expression: sp.Expr, parameter_polynomial: sp.Expr
    ) -> bool:
        coefficients = sp.Poly(expression, v).all_coeffs()
        modulus = sp.Poly(parameter_polynomial, lam)
        return all(
            sp.Poly(coefficient, lam).rem(modulus).is_zero
            for coefficient in coefficients
        )

    k9_subresultants = sp.subresultants(k9_lam, sp.diff(k9_lam, v), v)
    print(
        "K9_subresultant_degrees="
        + ",".join(str(sp.degree(item, v)) for item in k9_subresultants)
    )
    for name, exceptional in (("S6", s6), ("S7", s7)):
        zero_degrees = [
            sp.degree(item, v)
            for item in k9_subresultants
            if vanishes_mod_parameter(item, exceptional)
        ]
        print(f"{name}_zero_subresultant_degrees={zero_degrees}")

    joint_subresultants = sp.subresultants(k9_lam, l5, v)
    print(
        "K9_L5_subresultant_degrees="
        + ",".join(str(sp.degree(item, v)) for item in joint_subresultants)
    )
    joint_zero_degrees = [
        sp.degree(item, v)
        for item in joint_subresultants
        if vanishes_mod_parameter(item, t5)
    ]
    print(f"T5_joint_zero_subresultant_degrees={joint_zero_degrees}")


if __name__ == "__main__":
    main()
