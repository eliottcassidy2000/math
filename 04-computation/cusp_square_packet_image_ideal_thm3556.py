#!/usr/bin/env python3
"""Exact image, normalization, and conductor audit for the THM-3556 packet.

Work in integral target coordinates ``(L,T,W,R)=(L,T,2U_*,2S)``.  The
script verifies the two generators of the full image ideal, the resultant
certificate used to prove primeness, a smooth finite birational
normalization obtained by adjoining the marked root, the exact singular
locus on ``W != 0`` and on the ``W=0`` boundary, and the genus-two
unordered-pair conductor model.  Every truth-bearing gate remains active
under ``python -O``.

The commutative-algebra implications (regular sequence, unmixedness,
localization, and normalization) are stated in THM-3556.  This companion
checks their polynomial premises over QQ exactly.
"""

from __future__ import annotations

import gc
from itertools import combinations

import sympy as sp


def require(condition: bool, label: str) -> None:
    """Keep truth-bearing checks active under optimized execution."""
    if not bool(condition):
        raise ArithmeticError(f"FAILED: {label}")


def power_mod(base: sp.Poly, exponent: int, modulus: sp.Poly) -> sp.Poly:
    """Binary polynomial powering in the finite field of ``base``."""
    result = sp.Poly(1, *base.gens, modulus=base.get_modulus())
    factor = base.rem(modulus)
    while exponent:
        if exponent & 1:
            result = (result * factor).rem(modulus)
        factor = (factor * factor).rem(modulus)
        exponent >>= 1
    return result


def top_total_form(expression: sp.Expr,
                   variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Return the highest-total-degree homogeneous part."""
    polynomial = sp.Poly(expression, *variables)
    degree = polynomial.total_degree()
    return sp.expand(sum(
        coefficient * sp.prod(variable**exponent
                              for variable, exponent in zip(variables, monomial))
        for monomial, coefficient in polynomial.terms()
        if sum(monomial) == degree
    ))


def main() -> None:
    v, y = sp.symbols("v y")
    L0, T0, W0, R0 = sp.symbols("L T W R")
    variables = (L0, T0, W0, R0)

    W = sp.expand(2 + 2 * y - y**2 - 3 * v * y**2 + 9 * v * y)
    T = sp.expand(y**2 - 3 * v * W)
    R = sp.expand(2 * y**3 - 9 * v * W * y)
    L = sp.expand(v**2 * (4 * v * W - y**2))
    packet = (L, T, W, R)
    substitution = dict(zip(variables, packet))

    cubic = sp.expand(R0**2 - 4 * T0**3 - 27 * L0 * W0**2)
    quartic = sp.expand(
        27 * L0 * R0**2 - 243 * L0 * R0 * T0 - 81 * L0 * R0 * W0
        + 729 * L0 * R0 - 243 * L0 * T0 * W0**2
        + 486 * L0 * T0 * W0 - 243 * L0 * W0**3
        + 1215 * L0 * W0**2 - 1458 * L0 * W0
        + 3 * R0**2 * T0 + 7 * R0**2 * W0 - 18 * R0**2
        - 15 * R0 * T0**2 - 12 * R0 * T0 * W0 + 48 * R0 * T0
        + 3 * R0 * W0**2 - 16 * R0 * W0 + 36 * R0
        + 21 * T0**2 * W0**2 - 138 * T0**2 * W0 + 192 * T0**2
        + 6 * T0 * W0**3 - 36 * T0 * W0**2 + 48 * T0 * W0
        + W0**4 - 6 * W0**3 + 12 * W0**2 - 8 * W0
    )
    require(sp.expand(cubic.subs(substitution)) == 0,
            "cubic packet relation")
    require(sp.expand(quartic.subs(substitution)) == 0,
            "quartic packet relation")
    require(sp.Poly(quartic, *variables).total_degree() == 4,
            "quartic total target degree")
    require(len(sp.Poly(quartic, *variables).terms()) == 28,
            "quartic term count")
    cubic_top = top_total_form(cubic, variables)
    quartic_top = top_total_form(quartic, variables)
    require(sp.Poly(sp.gcd(cubic_top, quartic_top), *variables).total_degree() == 0,
            "coprime top forms for the filtered complete intersection")

    # Eliminate the marked root from its two inverse equations.
    Y = sp.symbols("Y")
    marked_cubic = Y**3 - 3 * T0 * Y + R0
    marked_quadratic = sp.expand(
        (W0 + 2 * T0) * Y**2
        - (2 * W0 + 6 * T0 + R0) * Y
        + W0**2 - 2 * W0 + 3 * R0
    )
    resultant = sp.expand(sp.resultant(marked_cubic, marked_quadratic, Y))
    coefficient_L = sp.diff(quartic, L0)
    require(sp.diff(coefficient_L, L0) == 0,
            "quartic is linear in L")
    require(sp.expand(27 * W0**2 * quartic
                      + coefficient_L * cubic - 27 * resultant) == 0,
            "resultant/localization identity")
    require(sp.Poly(resultant, R0).LC() == 1
            and sp.Poly(resultant, R0).degree() == 4,
            "resultant is monic quartic in R")

    specialized = sp.expand(resultant.subs({T0: -3, W0: -3}))
    expected_specialized = (
        R0**4 + 63 * R0**3 + 891 * R0**2 + 4320 * R0 + 216000
    )
    require(specialized == expected_specialized,
            "irreducible univariate specialization")
    modulus_seven = sp.Poly(specialized, R0, modulus=7)
    expected_mod_seven = sp.Poly(R0**4 + 2 * R0**2 + R0 + 1,
                                 R0, modulus=7)
    require(modulus_seven == expected_mod_seven,
            "specialization modulo seven")
    coordinate_mod_seven = sp.Poly(R0, R0, modulus=7)
    frobenius_7 = power_mod(coordinate_mod_seven, 7, modulus_seven)
    frobenius_49 = power_mod(coordinate_mod_seven, 49, modulus_seven)
    frobenius_2401 = power_mod(coordinate_mod_seven, 2401, modulus_seven)
    require(sp.gcd(modulus_seven, frobenius_7-coordinate_mod_seven).degree() == 0,
            "no linear factor modulo seven")
    require(sp.gcd(modulus_seven, frobenius_49-coordinate_mod_seven).degree() == 0,
            "no linear or quadratic factor modulo seven")
    require(frobenius_2401 == coordinate_mod_seven,
            "quartic Frobenius closure modulo seven")
    require(bool(modulus_seven.is_irreducible),
            "independent finite-field irreducibility gate")
    smooth_rational_point = {T0: 0, W0: 2, R0: 0}
    require(resultant.subs(smooth_rational_point) == 0
            and sp.diff(resultant, R0).subs(smooth_rational_point) == 64,
            "rational smooth point for absolute irreducibility")

    # Premises excluding a hidden W-supported minimal component.
    boundary_test = {L0: 0, T0: 1, W0: 0, R0: 2}
    require(cubic.subs(boundary_test) == 0
            and quartic.subs(boundary_test) == 270,
            "W-boundary noncontainment witness")
    jacobian_TW = sp.expand(
        sp.diff(T, v) * sp.diff(W, y) - sp.diff(T, y) * sp.diff(W, v)
    )
    require(jacobian_TW.subs({v: 0, y: 0}) == -12,
            "image transcendence degree two")

    # The linear subresultant gives the generic marked-root inverse.
    subresultants = sp.subresultants(marked_cubic, marked_quadratic, Y)
    linear_subresultant = sp.Poly(subresultants[-2], Y)
    alpha = sp.expand(linear_subresultant.coeff_monomial(Y))
    beta = sp.expand(linear_subresultant.coeff_monomial(1))
    expected_alpha = sp.expand(
        R0**2 + 6 * R0 * T0 + R0 * W0 - 12 * T0**3
        - 12 * T0**2 * W0 + 36 * T0**2 - 5 * T0 * W0**2
        + 28 * T0 * W0 - W0**3 + 6 * W0**2
    )
    expected_beta = sp.expand(
        -3 * R0**2 + 4 * R0 * T0**2 + 4 * R0 * T0 * W0
        - 18 * R0 * T0 - 4 * R0 * W0 - 6 * T0 * W0**2
        + 12 * T0 * W0 - 2 * W0**3 + 4 * W0**2
    )
    require(alpha == expected_alpha and beta == expected_beta,
            "linear subresultant coefficients")
    alpha_pullback = sp.expand(alpha.subs(substitution))
    require(alpha_pullback != 0
            and sp.expand((beta + y * alpha).subs(substitution)) == 0,
            "generic inverse y=-beta/alpha")
    require(sp.expand(marked_cubic.subs(substitution).subs(Y, y)) == 0
            and sp.expand(marked_quadratic.subs(substitution).subs(Y, y)) == 0,
            "direct marked-root substitution")

    # Hodge-normal bridge.  Defining each coefficient by the determinant
    # with two coordinate rows fixes the orientation without epsilon prose.
    # Both gradient rows annihilate both packet tangent rows because F(Z)=G(Z)=0.
    # The tangent and Hodge-normal bivectors therefore span the same line over
    # QQ(v,y).  One nonzero oriented component fixes their common multiplier;
    # all other five components then follow by linear algebra in the domain.
    gradient_F = [sp.diff(cubic, entry) for entry in variables]
    gradient_G = [sp.diff(quartic, entry) for entry in variables]
    first, second = 0, 2  # the (L,W) component
    row_first = [1, 0, 0, 0]
    row_second = [0, 0, 1, 0]
    normal_minor = sp.det(sp.Matrix(
        [gradient_F, gradient_G, row_first, row_second]
    ))
    source_minor = sp.expand(
        sp.diff(packet[first], v) * sp.diff(packet[second], y)
        - sp.diff(packet[first], y) * sp.diff(packet[second], v)
    )
    require(source_minor != 0
            and sp.expand(normal_minor.subs(substitution)
                          - 9 * alpha_pullback * source_minor) == 0,
            "oriented Hodge-normal multiplier on a nonzero component")
    del normal_minor, source_minor, gradient_F, gradient_G
    gc.collect()

    # Adjoin the integral marked root.  Put p=y^2-T=3vW.  Lexicographic
    # elimination of v gives the exact presentation of B[y].
    ell, p, w, marked_y = sp.symbols("ell p w marked_y")
    graph_ideal = [
        p - 3 * v * w,
        w - (2 + 2 * marked_y - marked_y**2)
        - 3 * v * marked_y * (3 - marked_y),
        ell - v**2 * (4 * v * w - marked_y**2),
    ]
    graph_groebner = sp.groebner(
        graph_ideal, v, ell, p, w, marked_y, order="lex", domain=sp.QQ
    )
    normalization_relations = [
        polynomial.as_expr() for polynomial in graph_groebner.polys
        if v not in polynomial.as_expr().free_symbols
    ]
    require(len(normalization_relations) == 4,
            "four-generator marked-root elimination ideal")
    del graph_groebner
    graph_substitution = {
        ell: L,
        p: 3 * v * W,
        w: W,
        marked_y: y,
    }
    require(all(sp.expand(relation.subs(graph_substitution)) == 0
                for relation in normalization_relations),
            "normalization presentation pulls back to zero")
    relation_jacobian = sp.Matrix([
        [sp.diff(relation, entry) for entry in (ell, p, w, marked_y)]
        for relation in normalization_relations
    ])
    rank_drop_minors = [
        sp.expand(relation_jacobian.extract(rows, columns).det())
        for rows in combinations(range(len(normalization_relations)), 2)
        for columns in combinations(range(4), 2)
    ]
    del relation_jacobian
    gc.collect()
    smoothness_groebner = sp.groebner(
        normalization_relations + rank_drop_minors,
        ell, p, w, marked_y, order="grevlex", domain=sp.QQ
    )
    require(len(smoothness_groebner.polys) == 1
            and smoothness_groebner.polys[0].as_expr() == 1,
            "marked-root finite model is smooth")

    # Exact singular locus of the projected resultant hypersurface.  The
    # two square-containment directions certify equality of radicals.
    resultant_derivatives = [sp.diff(resultant, entry)
                             for entry in (T0, W0, R0)]
    singular_groebner = sp.groebner(
        [resultant] + resultant_derivatives,
        R0, W0, T0, order="grevlex", domain=sp.QQ
    )
    require(singular_groebner.reduce(alpha**2)[1] == 0
            and singular_groebner.reduce(beta**2)[1] == 0,
            "singular ideal forces the linear subresultant to vanish")
    conductor_groebner = sp.groebner(
        [resultant, alpha], R0, W0, T0,
        order="grevlex", domain=sp.QQ
    )
    require(conductor_groebner.reduce(beta**2)[1] == 0
            and all(conductor_groebner.reduce(derivative**2)[1] == 0
                    for derivative in resultant_derivatives),
            "resultant-plus-alpha forces singularity set-theoretically")

    # Unordered pairs of common roots: s is the sum and p_pair the product.
    s, p_pair = sp.symbols("s p_pair")
    T_pair = (s**2 - p_pair) / 3
    R_pair = s * p_pair
    pair_remainder = sp.Poly(sp.rem(
        marked_quadratic.subs({T0: T_pair, R0: R_pair}),
        Y**2 - s * Y + p_pair, Y
    ), Y)
    equation_1 = sp.expand(
        -3 * W0 * s + 6 * W0 + 5 * p_pair * s - 6 * p_pair
        - 2 * s**3 + 6 * s**2
    )
    equation_2 = sp.expand(
        3 * W0**2 - 3 * W0 * p_pair - 6 * W0 + 2 * p_pair**2
        - 2 * p_pair * s**2 + 9 * p_pair * s
    )
    require(sp.expand(3 * pair_remainder.coeff_monomial(Y)
                      + equation_1) == 0
            and sp.expand(3 * pair_remainder.coeff_monomial(1)
                          - equation_2) == 0,
            "two-common-root remainder equations")
    pair_curve = sp.expand(
        -16 * W0**2 * s**2 + 36 * W0**2 * s - 24 * W0**2
        + 12 * W0 * s**4 - 79 * W0 * s**3 + 206 * W0 * s**2
        - 228 * W0 * s + 72 * W0 + 4 * s**6 - 42 * s**5
        + 126 * s**4 - 108 * s**3
    )
    require(sp.expand(sp.resultant(equation_1, equation_2, p_pair)
                      + 3 * pair_curve) == 0,
            "unordered-pair elimination")
    sextic = sp.expand(
        16 * s**6 - 168 * s**5 + 601 * s**4 - 1000 * s**3
        + 1048 * s**2 - 672 * s + 144
    )
    require(sp.factor(sp.discriminant(pair_curve, W0))
            == (5 * s - 6)**2 * sextic,
            "pair-curve discriminant")
    sextic_discriminant = sp.discriminant(sextic, s)
    require(sextic_discriminant == -26132324317743978381312,
            "squarefree genus-two sextic")
    require(sp.factor(pair_curve.subs(s, 3)) == -3 * W0 * (20 * W0 - 27),
            "known collision on pair curve")
    p_formula = sp.cancel(
        (3 * W0 * s - 6 * W0 + 2 * s**3 - 6 * s**2) / (5 * s - 6)
    )
    delta = sp.cancel(
        -3 * (4 * W0 * s - 8 * W0 + s**3 - 6 * s**2) / (5 * s - 6)
    )
    require(sp.cancel(s**2 - 4 * p_formula - delta) == 0,
            "ordered-sheet quadratic extension")
    known_pair = {s: 3, W0: sp.Rational(27, 20),
                  p_pair: sp.Rational(9, 20)}
    require(equation_1.subs(known_pair) == 0
            and equation_2.subs(known_pair) == 0
            and delta.subs({s: 3, W0: sp.Rational(27, 20)})
            == sp.Rational(36, 5),
            "known quadratic double fibre")
    exceptional_s = sp.Rational(6, 5)
    require(sp.expand(equation_1.subs(s, exceptional_s)
                      - 12 * (25 * W0 + 54) / 125) == 0,
            "exceptional pair-chart denominator")

    # The only W!=0 point where Q is the zero polynomial gives a triple
    # normalization fibre.
    triple_projected = {T0: sp.Rational(1, 2), W0: -1, R0: -1}
    zero_quadratic_groebner = sp.groebner(
        [W0 + 2*T0, 2*W0 + 6*T0 + R0, W0**2 - 2*W0 + 3*R0],
        R0, T0, W0, order="lex", domain=sp.QQ
    )
    expected_zero_quadratic_basis = [
        R0-W0, (2*T0+W0)/2, W0*(W0+1)
    ]
    require(len(zero_quadratic_groebner.polys) == 3
            and all(sp.expand(polynomial.as_expr()-expected) == 0
                    for polynomial, expected in zip(
                        zero_quadratic_groebner.polys,
                        expected_zero_quadratic_basis)),
            "unique W-nonzero zero-quadratic point")
    require(sp.Poly(marked_quadratic.subs(triple_projected), Y).is_zero
            and sp.discriminant(marked_cubic.subs(triple_projected), Y) != 0,
            "three distinct marked roots at the triple point")
    require(cubic.subs({L0: sp.Rational(1, 54), **triple_projected}) == 0,
            "triple point lift to the image complete intersection")

    # Boundary divisor and exact nonproper charts.  Normalizing the cusp by
    # T=t^2, R=2t^3 gives all W=0 points set-theoretically.
    t, lam = sp.symbols("t lambda")
    boundary_factor = sp.expand(
        6 * t**3 * (2 * t + 3)
        * (9 * L0 * (t - 3)**2 + t**4 - 4 * t**3 + 8 * t + 4)
    )
    require(sp.expand(quartic.subs({T0: t**2, W0: 0, R0: 2*t**3})
                      - boundary_factor) == 0,
            "complete W=0 boundary factorization")
    complete_intersection_jacobian = sp.Matrix([
        [sp.diff(cubic, entry) for entry in variables],
        [sp.diff(quartic, entry) for entry in variables],
    ])
    line_zero_jacobian = complete_intersection_jacobian.subs(
        {L0: lam, T0: 0, W0: 0, R0: 0}
    )
    require(all(line_zero_jacobian.extract((0, 1), columns).det() == 0
                for columns in combinations(range(4), 2)),
            "the first missing line is entirely singular")
    line_one_jacobian = complete_intersection_jacobian.subs(
        {L0: lam, T0: sp.Rational(9, 4), W0: 0,
         R0: -sp.Rational(27, 4)}
    )
    line_one_minors = [
        sp.factor(line_one_jacobian.extract((0, 1), columns).det())
        for columns in combinations(range(4), 2)
    ]
    nonzero_line_one_minors = [entry for entry in line_one_minors if entry != 0]
    line_one_factor = 2916 * lam + 169
    require(nonzero_line_one_minors
            and all(sp.cancel(entry / line_one_factor).free_symbols.isdisjoint({lam})
                    and sp.cancel(entry / line_one_factor) != 0
                    for entry in nonzero_line_one_minors),
            "second missing line has one singular point")
    graph_L = sp.cancel(
        -(t**4 - 4 * t**3 + 8 * t + 4) / (9 * (t - 3)**2)
    )
    graph_jacobian = complete_intersection_jacobian.subs(
        {L0: graph_L, T0: t**2, W0: 0, R0: 2 * t**3}
    )
    graph_minor_numerators = [
        sp.cancel(graph_jacobian.extract((0, 1), columns).det())
        .as_numer_denom()[0]
        for columns in combinations(range(4), 2)
    ]
    require(sp.factor(sp.gcd_list(graph_minor_numerators))
            == 8 * t**3 * (2 * t + 3),
            "no further W=0 singular points")
    require(graph_L.subs(t, -sp.Rational(3, 2))
            == -sp.Rational(169, 2916),
            "graph/second-line singular intersection")

    # For either y0=0 or y0=3, solve L=lambda after v=t^-1.  The implicit
    # equations have a simple formal root and send t=0 to the entire missing
    # target line, proving that the birational parameterization is nonfinite.
    u = sp.symbols("u")
    for y_branch, root_u, derivative_u in (
        (t * u, -sp.Rational(2, 9), 9),
        (3 + t * u, -sp.Rational(1, 9), -9),
    ):
        chart_substitution = {v: 1/t, y: y_branch}
        packet_W = sp.expand(W.subs(chart_substitution))
        packet_T = sp.expand(T.subs(chart_substitution))
        packet_R = sp.expand(R.subs(chart_substitution))
        packet_L = sp.expand(L.subs(chart_substitution))
        chart_W = sp.expand((t * y_branch**2 + t**3 * lam) / 4)
        chart_T = sp.expand((y_branch**2 - 3 * t**2 * lam) / 4)
        chart_R = sp.expand(-(y_branch**3 + 9*t**2*lam*y_branch) / 4)
        implicit_equation = sp.expand(packet_W - chart_W)
        require(implicit_equation.subs({t: 0, u: root_u}) == 0
                and sp.diff(implicit_equation, u).subs({t: 0, u: root_u})
                == derivative_u,
                "simple formal root in missing-line chart")
        require(sp.factor(packet_L - lam - 4*implicit_equation/t**3) == 0
                and sp.factor(packet_T-chart_T+3*implicit_equation/t) == 0
                and sp.factor(packet_R-chart_R
                              + 9*y_branch*implicit_equation/t) == 0,
                "exact missing-line chart identities")

    print("THM-3556 PACKET IMAGE-IDEAL AND NORMALIZATION AUDIT")
    print("coordinates=(L,T,W=2U,R=2S)")
    print("image_ideal=(degree_3_F,degree_4_G); G_terms=28")
    print("top_forms_coprime=True; complete_intersection_Hilbert_function=all_degrees")
    print("resultant_identity=27W^2G+A_F=27Res_Y(C,Q)")
    print("resultant_specialization_mod7=R^4+2R^2+R+1; irreducible=True")
    print("absolute_irreducibility_smooth_point=(T,W,R)=(0,2,0); dH/dR=64")
    print("W_local_prime=True; W_supported_minimal_prime=False")
    print("marked_root_inverse=y=-beta/alpha; alpha_nonzero=True")
    print("normalization=B[y]; elimination_generators=4; smooth=True")
    print("resultant_singular_radical=radical(H,alpha)=radical(H,dH)")
    print("hodge_normal_multiplier=N_ij(Z)/(M_ij)=9alpha for all six pairs")
    print("pair_curve_discriminant=(5s-6)^2D6; genus=2; disc(D6)="
          f"{sextic_discriminant}")
    print("known_pair=(s,W,p,delta)=(3,27/20,9/20,36/5)")
    print("triple_fibre=(L,T,W,R)=(1/54,1/2,-1,-1)")
    print("W0_components=vertical_line_0,vertical_line_1,graph_curve")
    print("W0_singular=vertical_line_0 plus (L,T,W,R)="
          "(-169/2916,9/4,0,-27/4)")
    print("parameterization_birational=True; finite=False")
    print("VERDICT=exact full image ideal, smooth normalization, and conductor audit")


if __name__ == "__main__":
    main()
