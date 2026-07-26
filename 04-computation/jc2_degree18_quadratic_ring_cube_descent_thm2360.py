#!/usr/bin/env python3
"""Exact companion for THM-2360's conditional Laurent-UFD cube descent."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def reciprocal(polynomial: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Expr:
    return sp.cancel(variable**degree * polynomial.subs(variable, 1 / variable))


def main() -> None:
    x, s, x0 = sp.symbols("x s x0")
    a0, a1, a2, a3 = sp.symbols("a0 a1 a2 a3")
    b0, b1, b2, b3, b4, b5 = sp.symbols("b0 b1 b2 b3 b4 b5")
    c0, c1, c2, c3, c4 = sp.symbols("c0 c1 c2 c3 c4")

    p = a3 * x**3 + a2 * x**2 + a1 * x + a0
    q = b5 * x**5 + b4 * x**4 + b3 * x**3 + b2 * x**2 + b1 * x + b0
    square = c4 * x**4 + c3 * x**3 + c2 * x**2 + c1 * x + c0
    h = x**2 - 1

    x_of_s = (s + 1 / s) / 2
    t_of_s = (s - 1 / s) / 2
    require(
        sp.cancel(t_of_s**2 - (x_of_s**2 - 1)) == 0,
        "Laurent parametrization stopped satisfying t^2=x^2-1",
    )

    p_s = sp.cancel(p.subs(x, x_of_s))
    q_s = sp.cancel(q.subs(x, x_of_s))
    square_s = sp.cancel(square.subs(x, x_of_s))
    a_plus = sp.cancel(7 * q_s + square_s * t_of_s)
    a_minus = sp.cancel(7 * q_s - square_s * t_of_s)
    b_polynomial = sp.Poly(sp.cancel(s**5 * a_plus), s)
    b_expression = b_polynomial.as_expr()
    b_dual = sp.Poly(reciprocal(b_expression, s, 10), s).as_expr()

    require(
        b_polynomial.degree() == 10
        and b_polynomial.LC() == (7 * b5 + c4) / 32
        and b_polynomial.TC() == (7 * b5 - c4) / 32,
        "degree-ten endpoint formulas changed",
    )
    require(
        sp.cancel(b_dual - s**5 * a_minus) == 0,
        "degree-ten reciprocal stopped representing conjugation",
    )
    require(
        sp.cancel(
            b_expression * b_dual
            - s**10 * (49 * q_s**2 - (x_of_s**2 - 1) * square_s**2)
        )
        == 0,
        "Laurent norm identity changed",
    )

    endpoint_product = sp.factor(b_polynomial.LC() * b_polynomial.TC())
    require(
        sp.expand(endpoint_product - (49 * b5**2 - c4**2) / 1024) == 0,
        "endpoint product changed",
    )
    endpoint_under_leading_relation = sp.factor(
        ((49 * b5**2 - c4**2) / 1024).subs(
            c4**2,
            4 * a3**3 + 49 * b5**2,
        )
    )
    require(
        endpoint_under_leading_relation == -a3**3 / 256,
        "leading norm relation stopped forcing nonzero endpoints",
    )

    x0_pair = sp.Poly(sp.cancel(s * (x_of_s - x0)), s)
    p_six = sp.Poly(sp.cancel(s**3 * p_s), s)
    require(
        x0_pair.as_expr() == (s**2 - 2 * x0 * s + 1) / 2
        and x0_pair.degree() == 2,
        "distinguished reciprocal pair changed",
    )
    require(
        p_six.degree() == 6
        and p_six.LC() == a3 / 8
        and p_six.TC() == a3 / 8
        and sp.cancel(reciprocal(p_six.as_expr(), s, 6) - p_six.as_expr())
        == 0,
        "cubic reciprocal sextic changed",
    )
    require(
        sp.cancel(
            s**10 * (-4 * (x_of_s - x0) * p_s**3)
            + 4 * x0_pair.as_expr() * p_six.as_expr() ** 3
        )
        == 0,
        "reciprocal-pair cube norm changed",
    )

    # Exact allocation controls.  A selected root from each reciprocal pair
    # goes into B; its reciprocal goes into B^vee.
    def pair(root: sp.Rational) -> sp.Expr:
        return (s - root) * (s - 1 / root)

    distinguished = sp.Rational(2)
    ordinary_roots = (sp.Rational(3), sp.Rational(4), sp.Rational(5))
    template_b = (s - distinguished) * sp.prod(
        (s - root) ** 3 for root in ordinary_roots
    )
    template_norm = pair(distinguished) * sp.prod(
        pair(root) ** 3 for root in ordinary_roots
    )
    template_ratio = sp.factor(
        template_b * reciprocal(template_b, s, 10) / template_norm
    )
    require(
        template_ratio.is_Number and template_ratio != 0,
        "generic reciprocal-pair allocation stopped being a constant norm",
    )

    # The distinguished linear factor may also be a root of the cubic.
    overlap_roots = (distinguished, sp.Rational(3), sp.Rational(4))
    overlap_b = (s - distinguished) * sp.prod(
        (s - root) ** 3 for root in overlap_roots
    )
    overlap_norm = pair(distinguished) * sp.prod(
        pair(root) ** 3 for root in overlap_roots
    )
    overlap_ratio = sp.factor(
        overlap_b * reciprocal(overlap_b, s, 10) / overlap_norm
    )
    overlap_factorization = dict(sp.factor_list(overlap_b, s)[1])
    require(
        sp.degree(overlap_b, s) == 10
        and overlap_factorization[s - distinguished] == 4
        and overlap_ratio.is_Number
        and overlap_ratio != 0,
        "linear/cube overlap control changed",
    )

    # Repeated roots of p are absorbed with their multiplicities into U_3.
    repeated_b = (
        (s - distinguished)
        * (s - 3) ** 6
        * (s - 4) ** 3
    )
    repeated_norm = pair(distinguished) * pair(sp.Rational(3)) ** 6 * pair(
        sp.Rational(4)
    ) ** 3
    repeated_ratio = sp.factor(
        repeated_b * reciprocal(repeated_b, s, 10) / repeated_norm
    )
    require(
        sp.degree(repeated_b, s) == 10
        and repeated_ratio.is_Number
        and repeated_ratio != 0,
        "repeated-root allocation control changed",
    )

    # Hostile control: changing to a nonmaximal generator preserves the
    # function field while introducing a common p,q factor and an index square.
    beta = sp.Integer(2)
    p0 = x + 1
    q0 = x**2 + 1
    p_nonmaximal = sp.expand((x - beta) ** 2 * p0)
    q_nonmaximal = sp.expand((x - beta) ** 3 * q0)
    discriminant0 = sp.expand(-4 * p0**3 - 27 * q0**2)
    discriminant_nonmaximal = sp.expand(
        -4 * p_nonmaximal**3 - 27 * q_nonmaximal**2
    )
    hostile_gcd = sp.gcd(
        sp.Poly(p_nonmaximal, x),
        sp.Poly(q_nonmaximal, x),
    )
    require(
        sp.expand(hostile_gcd.as_expr() - (x - beta) ** 2) == 0
        and sp.expand(
            discriminant_nonmaximal - (x - beta) ** 6 * discriminant0
        )
        == 0
        and q_nonmaximal.subs(x, 0) != 0,
        "nonmaximal-generator hostile control changed",
    )

    print("THM-2360 conditional Laurent-UFD cube descent exact companion")
    print("quadratic ring: C[x,t]/(t^2-x^2+1)=C[s,s^-1]")
    print("conjugation: s -> s^-1")
    print("B=s^5*(7q+S*t): degree 10 with nonzero constant")
    print("endpoint product under leading relation: -a3^3/256")
    print("norm: B*B^vee=-4*X0*P6^3")
    print("X0=(s^2-2*x0*s+1)/2; P6=s^3*p(x) is reciprocal")
    print(f"generic allocation norm ratio: {template_ratio}")
    print(f"linear/cube overlap multiplicity: 4; norm ratio: {overlap_ratio}")
    print(f"repeated-p-root norm ratio: {repeated_ratio}")
    print("hostile nonmaximal gcd: (x-2)^2")
    print("hostile discriminant multiplier: (x-2)^6")
    print("VERDICT: cube descent is exact conditional on Res(p,q)!=0 and q(x0)!=0")


if __name__ == "__main__":
    main()
