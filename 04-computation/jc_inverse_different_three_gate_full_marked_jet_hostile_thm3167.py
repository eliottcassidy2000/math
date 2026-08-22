#!/usr/bin/env python3
"""Exact companion for the inverse-different three-gate theorem (THM-3167).

The script verifies four logically separate points:

1. target shears preserve the inverse-spectral numerator q;
2. projective incidence q in K*f_T does not impose the constant-field gate;
3. a determinant-one fixed-plus-cubic hostile is produced by an actual
   supplied companion b, not by an abstract cofactor packet; and
4. that complete marked pair agrees with the diagonal pair to arbitrary
   prescribed finite deformation-jet depth.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def numerator_remainder(expression: sp.Expr, modulus: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    numerator, _denominator = sp.cancel(expression).as_numer_denom()
    return sp.expand(sp.rem(sp.Poly(numerator, variable), sp.Poly(modulus, variable)).as_expr())


def parameter_valuation(expression: sp.Expr, parameter: sp.Symbol) -> int:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_poly = sp.Poly(numerator, parameter)
    denominator_poly = sp.Poly(denominator, parameter)
    require(numerator_poly.as_expr() != 0, "valuation requested for zero numerator")
    require(denominator_poly.as_expr() != 0, "valuation requested for zero denominator")
    numerator_order = min(monomial[0] for monomial, _ in numerator_poly.terms())
    denominator_order = min(monomial[0] for monomial, _ in denominator_poly.terms())
    return numerator_order - denominator_order


def coefficientwise_parameter_order(
    expression: sp.Expr,
    polynomial_variable: sp.Symbol,
    parameter: sp.Symbol,
) -> int:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    require(not sp.Poly(denominator, polynomial_variable).degree(),
            "expected denominator independent of the primitive variable")
    orders = []
    for coefficient in sp.Poly(numerator, polynomial_variable).all_coeffs():
        if coefficient != 0:
            orders.append(parameter_valuation(coefficient / denominator, parameter))
    require(bool(orders), "coefficient order requested for zero polynomial")
    return min(orders)


def main() -> None:
    T, u, t, r, tau, v = sp.symbols("T u t r tau v", nonzero=True)

    # Target-shear covariance is a formal chain-rule cancellation.  The
    # symbols stand for coefficientwise target derivatives at fixed T.
    f_u, f_v, b_u, b_v, h_prime = sp.symbols("f_u f_v b_u b_v h_prime")
    q_before = f_v * b_u - f_u * b_v
    q_after = f_v * (b_u - h_prime * b_v) - (f_u - h_prime * f_v) * b_v
    require(sp.expand(q_after - q_before) == 0, "target shear changed q")

    # Exact fixed-plus-cubic supplied pair.
    g = T**3 - t
    d = u**3 - t
    f = sp.expand((T - u) * g)
    derivative = sp.diff(f, T)
    fixed_idempotent = g / d
    companion = sp.factor(
        r**3 * t * fixed_idempotent
        - (3 * u / r) * T**2 * (1 - fixed_idempotent)
    )
    q_raw = sp.cancel(
        sp.diff(f, t) * sp.diff(companion, u)
        - sp.diff(f, u) * sp.diff(companion, t)
    )
    q_reduced = sp.factor(
        r**-1 * derivative + (r**3 - r**-1) * g
    )
    require(numerator_remainder(q_raw - q_reduced, f, T) == 0,
            "supplied companion did not produce the claimed reduced q")

    # Fixed and cubic evaluations give the complete physical packet.
    require(sp.factor(q_reduced.subs(T, u) - r**3 * derivative.subs(T, u)) == 0,
            "fixed inverse-Jacobian scale failed")
    require(numerator_remainder(q_reduced - r**-1 * derivative, g, T) == 0,
            "cubic inverse-Jacobian scale failed")
    fixed_jacobian = r**-3
    cubic_jacobian = r
    require(sp.factor(fixed_jacobian * cubic_jacobian**3 - 1) == 0,
            "total physical-Jacobian product was not one")
    pointed_ratio = sp.factor(cubic_jacobian / fixed_jacobian)
    require(pointed_ratio == r**4, "pointed ratio was not r^4")

    # The pointed shifted resultant is the exact diagonal-line equation.
    q0 = sp.factor(q_reduced.subs(T, u))
    resultant_numerator = sp.factor(sp.resultant(
        g,
        q0 * (T - u) * sp.diff(g, T) - d * q_reduced,
        T,
    ))
    resultant_denominator = sp.factor(d**3 * sp.resultant(g, q_reduced, T))
    shifted_resultant = sp.factor(resultant_numerator / resultant_denominator)
    require(sp.factor(shifted_resultant - (r**4 - 1)**3) == 0,
            "shifted resultant was not (r^4-1)^3")

    # Arbitrarily deep hostiles in the full marked pair.  Target derivatives
    # preserve tau-order because tau is constant for partial_u,partial_t.
    jet_checks = 0
    derivative_jet_checks = 0
    reference_companion = sp.factor(companion.subs(r, 1))
    for n in range(1, 17):
        r_n = 1 + tau**n
        companion_difference = sp.cancel(companion.subs(r, r_n) - reference_companion)
        q_difference = sp.cancel(q_reduced.subs(r, r_n) - derivative)
        require(coefficientwise_parameter_order(companion_difference, T, tau) >= n,
                f"companion jet failed at N={n}")
        require(coefficientwise_parameter_order(q_difference, T, tau) >= n,
                f"q jet failed at N={n}")

        for scale in (r_n**3, r_n**-1, r_n**-3, r_n):
            require(parameter_valuation(scale - 1, tau) == n,
                    f"branch scale jet failed at N={n}")
        require(sp.factor(r_n**-3 * r_n**3 - 1) == 0,
                f"determinant-one identity failed at N={n}")
        defect = sp.factor((r_n**4 - 1)**3)
        require(parameter_valuation(defect, tau) == 3 * n,
                f"pointed defect order failed at N={n}")
        jet_checks += 1

        if n <= 8:
            for du in range(4):
                for dt in range(4 - du):
                    differentiated = sp.diff(companion_difference, u, du, t, dt)
                    require(coefficientwise_parameter_order(differentiated, T, tau) >= n,
                            f"target differential jet failed at N={n}, du={du}, dt={dt}")
                    derivative_jet_checks += 1

    # Projective incidence is not the constant-field condition.  The honest
    # polynomial map (xi,eta) -> (u=xi^4,v=xi*eta) has generic degree four.
    f4 = T**4 - u
    b4 = v * T**3 / u
    q4 = sp.factor(
        sp.diff(f4, v) * sp.diff(b4, u)
        - sp.diff(f4, u) * sp.diff(b4, v)
    )
    require(sp.factor(q4 - sp.diff(f4, T) / (4 * u)) == 0,
            "degree-four constant-field hostile failed")
    xi, eta = sp.symbols("xi eta", nonzero=True)
    P = xi**4
    Q = xi * eta
    polynomial_jacobian = sp.factor(
        sp.diff(P, xi) * sp.diff(Q, eta)
        - sp.diff(P, eta) * sp.diff(Q, xi)
    )
    require(polynomial_jacobian == 4 * xi**4,
            "direct polynomial Jacobian failed")

    # Positive punctured control with the same inverse quartic.
    b_positive = 4 * v * T**3
    q_positive = sp.factor(
        sp.diff(f4, v) * sp.diff(b_positive, u)
        - sp.diff(f4, u) * sp.diff(b_positive, v)
    )
    require(sp.factor(q_positive - sp.diff(f4, T)) == 0,
            "punctured unit-Jacobian control failed")

    print("THM-3167 inverse-different three-gate and full marked-pair jet hostile")
    print("target_shear_covariance=q_preserved_exactly")
    print("supplied_companion_q=r^-1*f_T+(r^3-r^-1)*(T^3-t) mod f")
    print("physical_packet=(r^-3,r,r,r); total_product=1; pointed_ratio=r^4")
    print("shifted_resultant=(r^4-1)^3")
    print(f"full_marked_pair_jet_checks={jet_checks} N_range=1..16")
    print(f"target_differential_jet_checks={derivative_jet_checks} order_total<=3 N_range=1..8")
    print("constant_field_hostile=f=T^4-u;b=v*T^3/u;q=f_T/(4u);Jac=4u")
    print("punctured_positive_control=f=T^4-u;b=4v*T^3;q=f_T;Jac=1")
    print("scope=local_rational_split_cover_not_connected_polynomial_A2")


if __name__ == "__main__":
    main()
