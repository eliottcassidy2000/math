#!/usr/bin/env python3
"""Exact THM-2767 certificate for split degrees 10, 14, and 18.

This checks the finite algebra used by the componentwise slope-four closure;
normalization, DVR, and projective arguments are proved in the theorem.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def coefficients(degree: int, p: sp.Expr, q: sp.Expr, r: sp.Expr, extra: int = 3):
    exponent = sp.Rational(degree, 4)
    cs = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + extra + 1):
        cs.append(sp.factor(sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * cs[index - step]
            for step in range(2, min(4, index) + 1)
        ) / index))
    return cs


def observables(degree: int, d: sp.Symbol, q: sp.Symbol, s: sp.Symbol):
    cs = coefficients(degree, 2*d, q, d**2-s)
    phi = sp.factor(4*cs[degree+1])
    psi = sp.factor(4*cs[degree+2])
    response = sp.factor(4*cs[degree+3] + 2*d*cs[degree+1])
    return phi, psi, response, cs


def local_order(expression: sp.Expr, q: sp.Symbol, s: sp.Symbol) -> int:
    return min(sum(monomial) for monomial, _ in sp.Poly(sp.expand(expression), q, s).terms())


def local_face(expression: sp.Expr, q: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), q, s)
    order = local_order(expression, q, s)
    return sp.factor(sum(
        coefficient*q**monomial[0]*s**monomial[1]
        for monomial, coefficient in polynomial.terms()
        if sum(monomial) == order
    ))


def faber_at_original_section(
    degree: int,
    cs: list[sp.Expr],
    w: sp.Symbol,
    beta: sp.Symbol,
    d: sp.Symbol,
    q: sp.Symbol,
    s: sp.Symbol,
    C: sp.Symbol,
    E: sp.Symbol,
) -> sp.Expr:
    faber = sum(cs[index]*w**(degree-index) for index in range(degree+1))
    return sp.factor(faber.subs({w: beta, d: C-beta**2, s: q*beta-E}))


def unit_ideal(expressions: list[sp.Expr], *variables: sp.Symbol) -> bool:
    return sp.groebner(expressions, *variables).contains(sp.Integer(1))


def proportional(left: sp.Expr, right: sp.Expr, *variables: sp.Symbol) -> bool:
    left_poly = sp.Poly(sp.expand(left), *variables)
    right_poly = sp.Poly(sp.expand(right), *variables)
    if left_poly.is_zero or right_poly.is_zero:
        return left_poly.is_zero and right_poly.is_zero
    pivot = right_poly.terms()[0][0]
    right_coefficient = right_poly.coeff_monomial(pivot)
    left_coefficient = left_poly.coeff_monomial(pivot)
    if left_coefficient == 0:
        return False
    ratio = sp.cancel(left_coefficient/right_coefficient)
    return sp.expand(left-ratio*right) == 0


def main() -> None:
    h, d, q, s, w = sp.symbols("h d q s w")
    omega, beta, C, E = sp.symbols("omega beta C E")

    expected_bases = {
        10: (d-3*s**2, s**3-sp.Rational(1, 64)),
        14: (
            d-160*s**5+15*s**2,
            s**6-s**3/sp.Integer(10)-sp.Rational(1, 6400),
        ),
        18: (
            d-sp.Rational(20480, 7)*s**8+896*s**5-sp.Rational(65, 7)*s**2,
            s**10-sp.Rational(5, 16)*s**7+sp.Rational(5, 1024)*s**4
            -s/sp.Integer(32768),
        ),
    }
    expected_prefix = {
        10: (
            q**3,
            5*omega**2,
            -sp.Rational(5, 32)*(32*omega**3+q),
        ),
        14: (
            q**4,
            -sp.Rational(7, 64)*(80*omega**3+q),
            sp.Rational(35, 64)*omega*(16*omega**3+q),
        ),
        18: (
            omega*q**5,
            sp.Rational(63, 128)*(32*omega**3+q),
            -sp.Rational(63, 128)*omega*(32*omega**3+3*q),
        ),
    }

    prefix_reduced = {}
    boundary_counts = {}
    for degree in (10, 14, 18):
        phi, psi, response, cs = observables(degree, d, q, s)
        k = (degree+2)//4
        require(
            sp.gcd(sp.Poly(phi, d, q, s), sp.Poly(psi, d, q, s)).total_degree() == 0,
            f"degree {degree} top flux gcd changed",
        )
        require(phi.subs(q, 0) == 0, f"degree {degree} q=0 Phi changed")
        require(response.subs(q, 0) == 0, f"degree {degree} q=0 response changed")
        psi_q0 = sp.Poly(sp.factor(psi.subs(q, 0)), d, s)
        require(
            len(psi_q0.terms()) == 1 and psi_q0.monoms()[0] == (0, k),
            f"degree {degree} q=0 boundary is not just P_infty",
        )

        phi_q1, psi_q1, response_q1 = (
            sp.factor(expression.subs(q, 1)) for expression in (phi, psi, response)
        )
        basis = sp.groebner([phi_q1, psi_q1], d, s, order="lex")
        expected = sp.groebner(list(expected_bases[degree]), d, s, order="lex")
        require(basis == expected, f"degree {degree} q=1 Groebner basis changed")
        eliminant = expected_bases[degree][1]
        boundary_counts[degree] = sp.degree(eliminant, s)

        jacobian = sp.det(sp.Matrix([
            [sp.diff(phi_q1, d), sp.diff(phi_q1, s)],
            [sp.diff(psi_q1, d), sp.diff(psi_q1, s)],
        ]))
        require(unit_ideal([phi_q1, psi_q1, response_q1], d, s),
                f"degree {degree} response vanished at a q-chart end")
        require(unit_ideal([phi_q1, psi_q1, jacobian], d, s),
                f"degree {degree} q-chart boundary lost transversality")

        phi_prefix = sp.factor(phi.subs({d: -omega**2, s: q*omega}))
        psi_prefix = sp.factor(psi.subs({d: -omega**2, s: q*omega}))
        prefix_gcd = sp.factor(sp.gcd(
            sp.Poly(phi_prefix, q, omega), sp.Poly(psi_prefix, q, omega)
        ).as_expr())
        phi_red = sp.factor(sp.cancel(phi_prefix/prefix_gcd))
        psi_red = sp.factor(sp.cancel(psi_prefix/prefix_gcd))
        expected_gcd, expected_phi_red, expected_psi_red = expected_prefix[degree]
        require(sp.expand(prefix_gcd-expected_gcd) == 0,
                f"degree {degree} prefix gcd formula changed")
        require(sp.expand(phi_red-expected_phi_red) == 0,
                f"degree {degree} reduced prefix Phi formula changed")
        require(sp.expand(psi_red-expected_psi_red) == 0,
                f"degree {degree} reduced prefix Psi formula changed")
        prefix_reduced[degree] = (prefix_gcd, phi_red, psi_red)

        if degree in (10, 14):
            require(
                unit_ideal([phi_q1, psi_q1, d+omega**2, s-omega], d, s, omega),
                f"degree {degree} acquired an exact-prefix-compatible q-end",
            )
        else:
            prefix_basis = sp.groebner(
                [phi_q1, psi_q1, d+omega**2, s-omega], d, s, omega,
                order="lex",
            )
            require(
                list(prefix_basis)
                == list(sp.groebner([d, s, omega], d, s, omega, order="lex")),
                "degree 18 exact-prefix exception is not exactly P_q",
            )
            section = faber_at_original_section(degree, cs, w, beta, d, q, s, C, E)
            section_poly = sp.Poly(sp.expand(section), beta, C, E, q)
            pure_q6 = section_poly.coeff_monomial(q**6)
            require(
                pure_q6 == -sp.Rational(21, 1024),
                "degree 18 pure q^6 section coefficient changed",
            )
            require(
                all(monomial == (0, 0, 0, 6) for monomial, coefficient in section_poly.terms()
                    if monomial[3] == 6 and coefficient != 0),
                "degree 18 acquired a second q-degree-six section term",
            )
            for monomial, coefficient in section_poly.terms():
                require(
                    monomial[0]+2*monomial[1]+4*monomial[2]+3*monomial[3]
                    == 18,
                    "degree 18 section lost weighted homogeneity",
                )
                if monomial != (0, 0, 0, 6):
                    negative_score = monomial[0]+3*monomial[3]
                    require(
                        negative_score < 18 or monomial[0] > 0,
                        "degree 18 section acquired a second exact -18alpha term",
                    )
            for lower in range(1, degree):
                if lower % 4 == 0:
                    continue
                _, _, _, lower_cs = observables(lower, d, q, s)
                lower_section = faber_at_original_section(
                    lower, lower_cs, w, beta, d, q, s, C, E
                )
                lower_poly = sp.Poly(sp.expand(lower_section), beta, C, E, q)
                for monomial, coefficient in lower_poly.terms():
                    require(
                        monomial[0]+2*monomial[1]+4*monomial[2]+3*monomial[3]
                        == lower,
                        f"E_{lower}(0) lost weighted homogeneity",
                    )
                    require(
                        monomial[0]+3*monomial[3] <= lower < 18,
                        f"E_{lower}(0) reached the degree-18 pole weight",
                    )

        # Uniform P_infty data.  Lower even rows align only at slope four;
        # odd rows have a nonzero constant Phi and the displayed R/Phi ratio.
        phi_face = local_face(phi.subs(d, 1), q, s)
        psi_face = local_face(psi.subs(d, 1), q, s)
        phi_model = sp.expand(sum(
            sp.binomial(k, exponent)
            * (-1)**((exponent-1)//2)
            * q**exponent*(-s)**(k-exponent)
            for exponent in range(1, k+1, 2)
        ))
        psi_model = sp.expand(sum(
            sp.binomial(k, exponent)
            * (-1)**(exponent//2)
            * q**exponent*(-s)**(k-exponent)
            for exponent in range(0, k+1, 2)
        ))
        require(local_order(phi.subs(d, 1), q, s) == k,
                f"degree {degree} Phi P_infty order changed")
        require(local_order(psi.subs(d, 1), q, s) == k,
                f"degree {degree} Psi P_infty order changed")
        require(sp.gcd(sp.Poly(phi_face, q, s), sp.Poly(psi_face, q, s)).total_degree() == 0,
                f"degree {degree} P_infty faces acquired a common direction")
        require(proportional(phi_face, phi_model, q, s),
                f"degree {degree} Phi face stopped being the odd part")
        require(proportional(psi_face, psi_model, q, s),
                f"degree {degree} Psi face stopped being the even part")
        require(sp.Poly(psi_face, q, s).coeff_monomial(s**k) != 0,
                f"degree {degree} Psi face lost pure s^k")
        if k % 2 == 0:
            require(sp.Poly(psi_face, q, s).coeff_monomial(q**k) != 0,
                    f"degree {degree} even-k pure q^k left Psi")
            require(sp.Poly(phi_face, q, s).coeff_monomial(q**k) == 0,
                    f"degree {degree} even-k pure q^k entered Phi")
        else:
            require(sp.Poly(phi_face, q, s).coeff_monomial(q**k) != 0,
                    f"degree {degree} odd-k pure q^k left Phi")
            require(sp.Poly(psi_face, q, s).coeff_monomial(q**k) == 0,
                    f"degree {degree} odd-k pure q^k entered Psi")
        difference = sp.factor(response+d*phi/2)
        require(local_order(difference.subs(d, 1), q, s) >= k+1,
                f"degree {degree} response lost its extra P_infty order")

        for lower in range(1, degree):
            if lower % 4 == 0:
                continue
            lower_phi, lower_psi, lower_response, _ = observables(lower, d, q, s)
            if lower % 2 == 0:
                lower_k = (lower+2)//4
                require(degree-lower+4*lower_k == degree+2,
                        f"degree {degree}, row {lower} lost slope-four alignment")
                require(local_order(lower_phi.subs(d, 1), q, s) == lower_k,
                        f"Phi_{lower} local order changed")
                require(local_order(lower_psi.subs(d, 1), q, s) == lower_k,
                        f"Psi_{lower} local order changed")
                require(local_order(lower_response.subs(d, 1), q, s) == lower_k,
                        f"R_{lower} local order changed")
                lower_difference = sp.factor(lower_response+d*lower_phi/2)
                if lower_difference != 0:
                    require(
                        local_order(lower_difference.subs(d, 1), q, s)
                        >= lower_k+1,
                        f"R_{lower}+d Phi_{lower}/2 lost its extra order",
                    )
                require(lower_phi.subs(q, 0) == 0,
                        f"even Phi_{lower}(q=0) changed")
                require(lower_response.subs(q, 0) == 0,
                        f"even R_{lower}(q=0) changed")
                lower_psi_q0 = sp.Poly(sp.factor(lower_psi.subs(q, 0)), d, s)
                require(
                    len(lower_psi_q0.terms()) == 1
                    and lower_psi_q0.monoms()[0] == (0, lower_k),
                    f"Psi_{lower}(q=0) stopped being a pure s-power",
                )
            else:
                phi0 = sp.factor(lower_phi.subs({d: 1, q: 0, s: 0}))
                response0 = sp.factor(lower_response.subs({d: 1, q: 0, s: 0}))
                require(phi0 != 0, f"odd Phi_{lower} lost its constant term")
                require(lower_psi.subs({q: 0, s: 0}) == 0,
                        f"odd Psi_{lower} acquired a constant term")
                require(
                    sp.factor(response0/phi0-sp.Rational(lower+1, 2*(lower+3))) == 0,
                    f"odd row {lower} response/Phi ratio changed",
                )

        # On q=0, the second flux is a nonzero triangular polynomial in s;
        # the odd first-flux columns have pairwise distinct leading d-degrees.
        odd_d_degrees = []
        for lower in range(1, degree, 2):
            if lower % 4 == 0:
                continue
            lower_phi, lower_psi, _, _ = observables(lower, d, q, s)
            require(lower_psi.subs(q, 0) == 0, f"odd Psi_{lower}(q=0) changed")
            phi_vertical = sp.Poly(sp.expand(lower_phi.subs(q, 0)), d)
            expected_degree = (lower+1)//2
            require(phi_vertical.degree() == expected_degree,
                    f"odd Phi_{lower} vertical degree changed")
            require(phi_vertical.LC() != 0, f"odd Phi_{lower} leading coefficient vanished")
            odd_d_degrees.append(expected_degree)
        require(len(odd_d_degrees) == len(set(odd_d_degrees)),
                f"degree {degree} odd vertical triangle collided")

    print("THM-2767 split degrees 10/14/18 exact closure certificate")
    print("status=FINITE_EXACT_SUPPORT_FOR_PROPOSED_COMPONENTWISE_CLOSURE")
    print("q1_index_cover_counts=" + ",".join(
        f"M{degree}:{boundary_counts[degree]}" for degree in (10, 14, 18)
    ))
    for degree in (10, 14, 18):
        gcd, phi_red, psi_red = prefix_reduced[degree]
        print(f"M{degree}_prefix_gcd={gcd}")
        print(f"M{degree}_prefix_reduced_Phi={phi_red}")
        print(f"M{degree}_prefix_reduced_Psi={psi_red}")
    print("M18_only_prefix_exception=P_q:[d:q:s]=[0:1:0]")
    print("M18_E18_z0_pure_q6=-21/1024")
    print("M18_E18_z0_support=all_weight18;pure_q6_unique;lower_Ej_weights_less_than18")
    print("prefix_display_formulas=hard_checked_exact")
    print("top_faces=exact_odd_even_parts;pure_qk_in_Psi_iff_k_even")
    print("P_infty_subslope_response_ratio_top=-1/2")
    print("lower_even_response_rows=strictly_later_below_slope4")
    print("odd_row_response_over_Phi=(j+1)/(2(j+3))")
    print("P_infty_no_active_odd_or_lambda=coprime_top_faces_force_unique_flux")
    print("vertical_second_flux=nonzero_triangular_polynomial_in_s")
    print("vertical_even_responses=zero_exactly_at_q0")
    print("vertical_first_flux_odd_d_degrees=1,2,...,(M/2)")
    print("finite_degree_basis=direct_exact_Faber_recurrence_not_all_degree_inheritance")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
