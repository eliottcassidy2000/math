#!/usr/bin/env python3
"""Exact hostile control for THM-2760's all-degree exact-prefix flux gcd.

The quantified proof is symbolic: after exact-prefix specialization it writes
the residual squarefreeness polynomial as a Schur--Szego composition of two
Jacobi transforms.  This companion independently checks the coefficient
identities, differential reconstruction, gcd, squarefreeness, negative-root
counts, and the nonlocalized omega boundary through k=18.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def falling(value: sp.Expr, length: int) -> sp.Expr:
    if length == 0:
        return sp.Integer(1)
    value = sp.sympify(value)
    return sp.prod(value - index for index in range(length))


def rising(value: sp.Expr, length: int) -> sp.Expr:
    if length == 0:
        return sp.Integer(1)
    value = sp.sympify(value)
    return sp.prod(value + index for index in range(length))


def faber_prefix_coefficients(k: int, r: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Return r^-k times the two exact-prefix coefficients by recurrence."""

    exponent = sp.Rational(2 * k - 1, 2)
    top_index = 4 * k
    coefficients = [sp.Integer(1)]
    quartic = {2: sp.Integer(-2), 3: r, 4: 1 - r}
    for index in range(1, top_index + 1):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
        ) / index
        coefficients.append(sp.factor(value))

    p_poly = sp.Poly(sp.expand(coefficients[4 * k - 1]), r)
    q_poly = sp.Poly(sp.expand(coefficients[4 * k]), r)
    divisor = sp.Poly(r**k, r)
    p_quotient, p_remainder = sp.div(p_poly, divisor)
    q_quotient, q_remainder = sp.div(q_poly, divisor)
    require(p_remainder.is_zero, f"k={k}: first flux lost r^k divisibility")
    require(q_remainder.is_zero, f"k={k}: second flux lost r^k divisibility")
    return sp.factor(p_quotient.as_expr()), sp.factor(q_quotient.as_expr())


def explicit_prefix_coefficients(k: int, r: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Return the independent binomial formulas for P_k and Q_k."""

    exponent = sp.Rational(2 * k - 1, 2)
    p_value = sp.Integer(0)
    q_value = sp.Integer(0)
    for ell in range((k - 1) // 3 + 1):
        common = sp.binomial(exponent, k + ell)
        p_coefficient = (
            common
            * (-1) ** (k - 1 - 3 * ell)
            * sp.Integer(2) ** (k - 1 - 3 * ell)
            * sp.binomial(k - 1 - ell, 2 * ell)
        )
        p_value += p_coefficient * r**ell
    for ell in range(k // 3 + 1):
        common = sp.binomial(exponent, k + ell)
        q_coefficient = (
            common
            * (-1) ** (k - 3 * ell)
            * sp.Integer(2) ** (k - 3 * ell - 1)
            * (
                2 * sp.binomial(k - 1 - ell, 2 * ell - 1)
                + sp.binomial(k - 1 - ell, 2 * ell)
            )
        )
        q_value += q_coefficient * r**ell
    return sp.factor(p_value), sp.factor(q_value)


def mod_three_parameters(k: int) -> tuple[int, sp.Rational, sp.Rational]:
    n, residue = divmod(k, 3)
    if residue == 0:
        return n, sp.Rational(3 * n - 1, 3), sp.Rational(3 * n - 2, 3)
    if residue == 1:
        return n, sp.Rational(3 * n + 1, 3), sp.Rational(3 * n - 1, 3)
    return n, sp.Rational(3 * n + 2, 3), sp.Rational(3 * n + 1, 3)


def schur_factors(k: int, r: sp.Symbol) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    n, a_parameter, b_parameter = mod_three_parameters(k)
    first = sp.Integer(0)
    second = sp.Integer(0)
    composition = sp.Integer(0)
    for ell in range(n + 1):
        first_coefficient = falling(a_parameter, ell) / falling(k - 1, ell)
        second_coefficient = falling(b_parameter, ell) / rising(k + 1, ell)
        binomial = sp.binomial(n, ell)
        first += binomial * first_coefficient * r**ell
        second += binomial * second_coefficient * r**ell
        composition += (
            binomial
            * first_coefficient
            * second_coefficient
            * (sp.Rational(27, 32) * r) ** ell
        )
    return sp.factor(first), sp.factor(second), sp.factor(composition)


def negative_simple(polynomial: sp.Expr, variable: sp.Symbol, degree: int) -> bool:
    if degree == 0:
        return polynomial != 0
    poly = sp.Poly(polynomial, variable)
    return (
        poly.degree() == degree
        and sp.gcd(poly, poly.diff()).degree() == 0
        and poly.count_roots(-sp.oo, 0) == degree
    )


def main() -> None:
    r, omega, q = sp.symbols("r omega q")
    checked_k = tuple(range(1, 19))
    representative: dict[str, sp.Expr] = {}

    for k in checked_k:
        recurrence_p, recurrence_q = faber_prefix_coefficients(k, r)
        formula_p, formula_q = explicit_prefix_coefficients(k, r)
        require(
            sp.expand(recurrence_p - formula_p) == 0,
            f"k={k}: first explicit coefficient formula failed",
        )
        require(
            sp.expand(recurrence_q - formula_q) == 0,
            f"k={k}: second explicit coefficient formula failed",
        )

        theta_p = sp.expand(r * sp.diff(formula_p, r))
        theta_q = sp.expand(r * sp.diff(formula_q, r))
        require(
            sp.expand(k * formula_q - 3 * theta_q + k * formula_p + theta_p)
            == 0,
            f"k={k}: Euler differential relation failed",
        )

        squarefree_carrier = sp.factor(formula_p - 3 * formula_q)
        carrier_constant = squarefree_carrier.subs(r, 0)
        require(carrier_constant != 0, f"k={k}: carrier lost its constant term")
        normalized_carrier = sp.factor(squarefree_carrier / carrier_constant)
        reconstructed_p = sp.factor(
            (k * squarefree_carrier - 3 * r * sp.diff(squarefree_carrier, r))
            / (4 * k)
        )
        reconstructed_q = sp.factor(
            -(k * squarefree_carrier + r * sp.diff(squarefree_carrier, r))
            / (4 * k)
        )
        require(
            sp.expand(reconstructed_p - formula_p) == 0,
            f"k={k}: P reconstruction from S failed",
        )
        require(
            sp.expand(reconstructed_q - formula_q) == 0,
            f"k={k}: Q reconstruction from S failed",
        )

        n, _, _ = mod_three_parameters(k)
        factorial_carrier = sum(
            sp.Rational(
                k
                * sp.factorial(k)
                * sp.factorial(k - 1 - ell),
                32**ell
                * sp.factorial(ell)
                * sp.factorial(k - 3 * ell)
                * sp.factorial(k + ell),
            )
            * r**ell
            for ell in range(n + 1)
        )
        require(
            sp.expand(normalized_carrier - factorial_carrier) == 0,
            f"k={k}: normalized carrier factorial formula failed",
        )

        first_factor, second_factor, schur_composition = schur_factors(k, r)
        require(
            sp.expand(normalized_carrier - schur_composition) == 0,
            f"k={k}: Schur--Szego factorization failed",
        )
        require(
            negative_simple(first_factor, r, n),
            f"k={k}: first Jacobi transform lost simple negative roots",
        )
        require(
            negative_simple(second_factor, r, n),
            f"k={k}: second Jacobi transform lost simple negative roots",
        )
        require(
            negative_simple(normalized_carrier, r, n),
            f"k={k}: carrier lost simple negative roots",
        )
        require(
            sp.gcd(sp.Poly(formula_p, r), sp.Poly(formula_q, r)).degree() == 0,
            f"k={k}: residual exact-prefix fluxes acquired a common root",
        )

        localized_first_flux = sp.expand(r ** (k - 1) * formula_p)
        localized_second_flux = sp.expand(r**k * formula_q)
        localized_gcd = sp.gcd(
            sp.Poly(localized_first_flux, r),
            sp.Poly(localized_second_flux, r),
        ).monic()
        require(
            sp.expand(localized_gcd.as_expr() - r ** (k - 1)) == 0,
            f"k={k}: localized flux gcd is not r^(k-1)",
        )

        # Hostile boundary: before omega is inverted and q=r*omega^3 is
        # normalized, an extra omega occurs exactly in the k=2 mod 3 lane.
        dimensional_first = sp.cancel(
            omega ** (k - 1) * q ** (k - 1) * formula_p.subs(r, q / omega**3)
        )
        dimensional_second = sp.cancel(
            omega**k * q**k * formula_q.subs(r, q / omega**3)
        )
        first_denominator = sp.cancel(dimensional_first).as_numer_denom()[1]
        second_denominator = sp.cancel(dimensional_second).as_numer_denom()[1]
        require(
            not first_denominator.has(omega, q)
            and not second_denominator.has(omega, q),
            f"k={k}: dimensional flux unexpectedly retained omega denominators",
        )
        dimensional_gcd = sp.gcd(
            sp.Poly(sp.expand(dimensional_first), omega, q),
            sp.Poly(sp.expand(dimensional_second), omega, q),
        ).monic()
        expected_omega = 1 if k % 3 == 2 else 0
        require(
            sp.expand(
                dimensional_gcd.as_expr()
                - omega**expected_omega * q ** (k - 1)
            )
            == 0,
            f"k={k}: nonlocalized omega hostile changed",
        )

        if k == 6:
            representative = {
                "P6": formula_p,
                "Q6": formula_q,
                "S6_normalized": normalized_carrier,
                "resultant6": sp.factor(sp.resultant(formula_p, formula_q, r)),
            }

        # The projective omega=0 chart cannot be reached by r=q/omega^3.
        # There q is nonzero and the quartic is 1+q*z^3, so a coefficient is
        # nonzero exactly when its index is divisible by three.  The two
        # fluxes vanish together only for k=2 mod 3; on precisely that lane
        # the third response survives.  The affine q=0 corner is excluded.
        exponent = sp.Rational(2 * k - 1, 2)

        def root_zero_coefficient(index: int) -> sp.Expr:
            if index % 3:
                return sp.Integer(0)
            power = index // 3
            return sp.binomial(exponent, power) * q**power

        root_zero_phi = 4 * root_zero_coefficient(4 * k - 1)
        root_zero_psi = 4 * root_zero_coefficient(4 * k)
        root_zero_response = 4 * root_zero_coefficient(4 * k + 1)
        common_root_zero = root_zero_phi == 0 and root_zero_psi == 0
        require(
            common_root_zero == (k % 3 == 2),
            f"k={k}: omega=0 common-flux residue classification changed",
        )
        if common_root_zero:
            require(
                root_zero_response != 0,
                f"k={k}: omega=0 third response vanished on the residual lane",
            )

    print("EXACT-PREFIX EVEN FABER FLUX GCD -- THM-2760")
    print("family=m=4k-2,alpha=k-1/2")
    print("specialization=d=-omega^2,q=r*omega^3,s=r*omega^4")
    print("normalized_quartic=(1-x^2)^2+r*x^3*(1-x)")
    print("checked_k=1..18")
    print("independent_paths=Faber_recurrence,binomial_coefficients")
    print("differential_relation=(k-3theta)Q=-(k+theta)P")
    print("carrier=S=P-3Q;P=(k-3theta)S/(4k);Q=-(k+theta)S/(4k)")
    print("carrier_formula=kk!(k-1-l)!/[32^l*l!*(k-3l)!*(k+l)!]")
    print("mod3_schur_factors=two_Jacobi_transforms")
    print("factor_and_carrier_roots=simple_negative_exact_Sturm_count")
    print("localized_gcd(Phi_m/q,Psi_m)=r^(k-1)")
    print("nonlocalized_extra_omega=iff_k_mod_3_equals_2")
    print("omega_zero_projective_chart=q_nonzero")
    print("omega_zero_common_flux=iff_k_mod_3_equals_2,iff_m_mod_12_equals_6")
    print("omega_zero_response=nonzero_on_projective_common_flux_lane")
    for label in ("P6", "Q6", "S6_normalized", "resultant6"):
        print(f"{label}={sp.sstr(representative[label])}")
    print("PASS")


if __name__ == "__main__":
    main()
