#!/usr/bin/env python3
"""Exact referee for THM-2341's deep-wall local genus split.

The script verifies the normalized blow-up discriminant, identifies the
unique common P,Q ratio and its simple Q-zero, checks that y=0 is etale,
and reconstructs the repeated-root gcds at all three THM-2338 ratios.
The Hensel/Eisenstein and Riemann--Hurwitz deductions are mathematical
steps in the theorem.  Every executable check remains active under
optimized Python.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def reduce_c_square(expression: sp.Expr, cvar: sp.Symbol, ratio: sp.Expr) -> sp.Expr:
    """Reduce a polynomial in C modulo C^2-ratio."""
    polynomial = sp.Poly(sp.expand(expression), cvar)
    answer = sp.Integer(0)
    for (power,), coefficient in polynomial.terms():
        answer += (
            coefficient
            * ratio ** (power // 2)
            * (cvar if power % 2 else 1)
        )
    return sp.factor(answer)


def main() -> None:
    y, z = sp.symbols("y z")
    bvar, cvar = sp.symbols("B C")

    wall_p = 35 * y**2 * (54 * bvar + 7 * y**2)
    wall_q = 7 * y**3 * (
        1620 * bvar * y + 26244 * cvar + 77 * y**3
    )
    p_scale = sp.Rational(16, 964467)
    q_scale = sp.Rational(64, 703096443)

    # Blow up the collapsed y=0 fibre by v=yz.
    normalized_cubic = sp.expand(
        z**3
        + p_scale * (wall_p / y**2) * z
        + q_scale * (wall_q / y**3)
    )
    residual = (
        7889 * y**6
        + 211680 * bvar * y**4
        + 1047816 * cvar * y**3
        + 1814400 * bvar**2 * y**2
        + 22044960 * bvar * cvar * y
        + 2916000 * bvar**3
        + 178564176 * cvar**2
    )
    require(
        sp.expand(
            sp.discriminant(normalized_cubic, z)
            + sp.Rational(4096, 96873331012983) * residual
        )
        == 0,
        "the normalized blow-up discriminant changed",
    )

    tvar = sp.symbols("t")
    t_zero = -sp.Rational(361, 30618)
    t_plus = (
        -sp.Rational(250, 13041)
        + sp.Rational(500, 117369) * sp.sqrt(3)
    )
    t_minus = (
        -sp.Rational(250, 13041)
        - sp.Rational(500, 117369) * sp.sqrt(3)
    )
    ratio_data = (
        (t_zero, -sp.Rational(486, 19)),
        (
            t_plus,
            -sp.Rational(81, 5) - sp.Rational(27, 5) * sp.sqrt(3),
        ),
        (
            t_minus,
            -sp.Rational(81, 5) + sp.Rational(27, 5) * sp.sqrt(3),
        ),
    )

    # The y=0 chart is etale at all three ratios.
    residual_at_zero = 2916000 + 178564176 * tvar
    quadratic_ratio = 62500 + 7654500 * tvar + 199644669 * tvar**2
    require(
        sp.factor(residual_at_zero.subs(tvar, t_zero)) == 810648,
        "the t0 blow-up fibre changed",
    )
    zero_quadratic_resultant = sp.resultant(
        residual_at_zero,
        quadratic_ratio,
        tvar,
    )
    require(
        zero_quadratic_resultant == -295233008801472000000,
        "a quadratic-ratio blow-up fibre became singular at y=0",
    )

    # Away from y=0, P and Q have a common root precisely on the linear
    # discriminant factor, hence precisely at t0 on B=1.
    p_reduced = 54 * bvar + 7 * y**2
    q_reduced = 1620 * bvar * y + 26244 * cvar + 77 * y**3
    common_resultant = sp.factor(sp.resultant(p_reduced, q_reduced, y))
    require(
        sp.expand(
            common_resultant
            - 7715736 * (361 * bvar**3 + 30618 * cvar**2)
        )
        == 0,
        "the nonzero common-root resultant changed",
    )
    require(
        sp.simplify(361 + 30618 * t_plus) != 0
        and sp.simplify(361 + 30618 * t_minus) != 0,
        "a quadratic ratio acquired a common P,Q root",
    )

    # At t0 the common root is explicit and Q has a simple zero.  The
    # depressed cubic is therefore Eisenstein in the local parameter.
    common_root = -sp.Rational(486, 19) * cvar
    normalized_wall_p = wall_p.subs(bvar, 1)
    normalized_wall_q = wall_q.subs(bvar, 1)
    require(
        reduce_c_square(
            normalized_wall_p.subs(y, common_root),
            cvar,
            t_zero,
        )
        == 0
        and reduce_c_square(
            normalized_wall_q.subs(y, common_root),
            cvar,
            t_zero,
        )
        == 0,
        "the t0 repeated point stopped being a triple cubic root",
    )
    p_prime_at_common = reduce_c_square(
        sp.diff(normalized_wall_p, y).subs(y, common_root),
        cvar,
        t_zero,
    )
    q_prime_at_common = reduce_c_square(
        sp.diff(normalized_wall_q, y).subs(y, common_root),
        cvar,
        t_zero,
    )
    require(
        p_prime_at_common == sp.Rational(1837080, 19) * cvar
        and q_prime_at_common == -sp.Rational(4251528, 19) * cvar,
        "the t0 local simple-zero data changed",
    )

    # All three residual discriminants have exactly one double root.
    # At t0 this is the common P,Q root; at t+/- it is a split quadratic
    # node because P,Q are not both zero there.
    gcd_degrees: list[int] = []
    for ratio, root_scale in ratio_data:
        cvalue = sp.sqrt(ratio)
        specialized = sp.Poly(
            residual.subs({bvar: 1, cvar: cvalue}),
            y,
            extension=True,
        )
        common = sp.gcd(specialized, specialized.diff()).monic()
        require(
            common.degree() == 1
            and sp.simplify(
                common.as_expr() - (y - root_scale * cvalue)
            )
            == 0,
            f"repeated-root gcd changed at ratio {ratio}",
        )
        gcd_degrees.append(common.degree())

    # Infinity remains the same separable three-point fibre as in
    # THM-2332.
    infinity_parameter = sp.symbols("infinity_parameter")
    infinity_cubic = (
        1127
        - 138915 * infinity_parameter
        + 1607445 * infinity_parameter**2
        - 26040609 * infinity_parameter**3
    )
    infinity_discriminant = sp.discriminant(
        infinity_cubic,
        infinity_parameter,
    )
    require(
        infinity_discriminant == -153384762202971019112448,
        "the infinity fibre lost separability",
    )

    # Four simple transpositions contribute four.  The t0 Eisenstein
    # point contributes two more; the two split nodes contribute zero.
    ramification_totals = (4 + 2, 4, 4)
    genera = tuple((total - 4) // 2 for total in ramification_totals)
    require(
        ramification_totals == (6, 4, 4)
        and genera == (1, 0, 0),
        "Riemann--Hurwitz genus ledger changed",
    )

    print("THM-2341 exact deep-wall local-genus referee")
    print("blow-up discriminant=(-4096/96873331012983)*G6")
    print(
        "common-root resultant="
        "7715736*(361B^3+30618C^2)"
    )
    print(f"y=0: G6(0,t0)=810648, quadratic resultant={zero_quadratic_resultant}")
    print(
        "t0 common root=-486C/19, "
        f"Q'={q_prime_at_common}"
    )
    print(f"multiple-root gcd degrees={gcd_degrees}")
    print(f"ramification totals={ramification_totals}")
    print(f"normalization genera={genera}")
    print("VERDICT: t0 is genus one; t+ and t- are genus zero")


if __name__ == "__main__":
    main()
