#!/usr/bin/env python3
"""Independent branch/regular-representation audit for THM-3167.

Unlike the primary companion, this script does not use a quartic resultant
or reduce the coefficientwise q formula first.  It evaluates the supplied
companion on the fixed factor and on the cubic regular representation,
computes both forward Jacobians directly on the source components, and
reconstructs the reduced q from those branch inverse Jacobians.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def tau_valuation(expression: sp.Expr, tau: sp.Symbol) -> int:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    require(numerator != 0, "valuation requested for zero")
    numerator_terms = sp.Poly(numerator, tau).terms()
    denominator_terms = sp.Poly(denominator, tau).terms()
    return min(m[0] for m, _ in numerator_terms) - min(m[0] for m, _ in denominator_terms)


def matrix_polynomial(poly: sp.Expr, variable: sp.Symbol, matrix: sp.Matrix) -> sp.Matrix:
    numerator, denominator = sp.cancel(poly).as_numer_denom()
    require(sp.Poly(denominator, variable).degree() == 0,
            "matrix evaluation needs a scalar denominator")
    result = sp.zeros(matrix.rows)
    for (exponent,), coefficient in sp.Poly(numerator, variable).terms():
        result += coefficient * matrix**exponent
    return sp.simplify(result / denominator)


def main() -> None:
    T, u, t, r, tau = sp.symbols("T u t r tau", nonzero=True)

    # The regular representation of y with y^3=t in basis 1,y,y^2.
    Y = sp.Matrix([[0, 0, t], [1, 0, 0], [0, 1, 0]])
    identity = sp.eye(3)
    require(Y**3 == t * identity, "cubic regular representation failed")

    g = T**3 - t
    d = u**3 - t
    f = sp.expand((T - u) * g)
    derivative = sp.diff(f, T)
    e = g / d
    companion = sp.factor(r**3 * t * e - (3 * u / r) * T**2 * (1 - e))

    # CRT evaluations of b use no coefficientwise target differentiation.
    require(sp.factor(companion.subs(T, u) - r**3 * t) == 0,
            "fixed companion value failed")
    cubic_b = matrix_polynomial(companion + (3 * u / r) * T**2, T, Y)
    require(all(sp.factor(entry) == 0 for entry in cubic_b),
            "cubic companion value failed")

    # Direct forward maps on the two source components.
    X, Z = sp.symbols("X Z", nonzero=True)
    fixed_u = X
    fixed_t = Z / r**3
    fixed_jacobian = sp.factor(
        sp.diff(fixed_u, X) * sp.diff(fixed_t, Z)
        - sp.diff(fixed_u, Z) * sp.diff(fixed_t, X)
    )
    cubic_u = -r * Z / (3 * X**2)
    cubic_t = X**3
    cubic_jacobian = sp.factor(
        sp.diff(cubic_u, X) * sp.diff(cubic_t, Z)
        - sp.diff(cubic_u, Z) * sp.diff(cubic_t, X)
    )
    require(fixed_jacobian == r**-3, "fixed forward Jacobian failed")
    require(cubic_jacobian == r, "cubic forward Jacobian failed")
    require(sp.factor(fixed_jacobian * cubic_jacobian**3 - 1) == 0,
            "direct total Jacobian product failed")

    # Reconstruct q=f_T/J from its two factors and verify by regular
    # representation.  This is the CRT path independent of the primary
    # coefficientwise q calculation.
    q_crt = sp.factor(r**-1 * derivative + (r**3 - r**-1) * g)
    require(sp.factor(q_crt.subs(T, u) - derivative.subs(T, u) / fixed_jacobian) == 0,
            "fixed CRT q failed")
    cubic_q = matrix_polynomial(q_crt - derivative / cubic_jacobian, T, Y)
    require(all(sp.factor(entry) == 0 for entry in cubic_q),
            "cubic CRT q failed")

    rho = sp.factor(cubic_jacobian / fixed_jacobian)
    norm_matrix = (rho - 1) * identity
    shifted_norm = sp.factor(norm_matrix.det())
    require(sp.factor(shifted_norm - (r**4 - 1)**3) == 0,
            "regular-representation shifted norm failed")

    # A longer finite-jet range checks the universal pro-jet mechanism using
    # only the four branch scales and the two CRT companion values.
    finite_jet_checks = 0
    for n in range(1, 33):
        r_n = 1 + tau**n
        scales = (r_n**3, r_n**-1, r_n**-3, r_n)
        for scale in scales:
            require(sp.expand(sp.series(scale - 1, tau, 0, n).removeO()) == 0,
                    f"truncated branch scale failed at N={n}")
            require(tau_valuation(scale - 1, tau) == n,
                    f"first branch-scale defect order failed at N={n}")
        require(tau_valuation((r_n**4 - 1)**3, tau) == 3 * n,
                f"shifted norm order failed at N={n}")
        require(sp.factor(r_n**-3 * r_n**3 - 1) == 0,
                f"total determinant failed at N={n}")
        finite_jet_checks += 1

    # Independent direct-map check of the constant-field hostile.
    xi, eta = sp.symbols("xi eta", nonzero=True)
    U = xi**4
    V = xi * eta
    jacobian = sp.factor(
        sp.diff(U, xi) * sp.diff(V, eta)
        - sp.diff(U, eta) * sp.diff(V, xi)
    )
    require(jacobian == 4 * U, "degree-four polynomial hostile Jacobian failed")
    inverse_xi_u = sp.Rational(1, 4) / xi**3
    inverse_eta_v = 1 / xi
    inverse_xi_v = 0
    inverse_eta_u = -eta / (4 * xi**4)
    inverse_determinant = sp.factor(
        inverse_xi_u * inverse_eta_v - inverse_xi_v * inverse_eta_u
    )
    require(sp.factor(inverse_determinant - 1 / (4 * U)) == 0,
            "degree-four hostile inverse determinant failed")

    print("THM-3167 independent branch and regular-representation audit")
    print("companion_CRT_values=(r^3*t,-3u*y^2/r)")
    print("direct_forward_Jacobians=(r^-3,r); total_product=1")
    print("CRT_q_branch_scales=(r^3,r^-1) relative_to_f_T")
    print("regular_representation_shifted_norm=(r^4-1)^3")
    print(f"finite_jet_checks={finite_jet_checks} N_range=1..32")
    print("polynomial_constant_field_hostile_direct_Jac=4u;inverse_det=1/(4u)")
    print("independent_audit=PASS")


if __name__ == "__main__":
    main()
