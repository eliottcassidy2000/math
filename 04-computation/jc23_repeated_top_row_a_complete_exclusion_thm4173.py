#!/usr/bin/env python3
"""Exact bridge audit for THM-4173.

The already-canonical THM-4147 resultant is checked at an exact hostile
row-A control.  The new content is the symbolic differential identity that
turns every Keller critical intersection into a reduced point, eliminating
the old discriminant and projected-fibre-separation hypotheses.

This audit intentionally does not compute Disc(Q_19), Disc(R_19), or
Q_19(-1/6).
"""

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuation(poly, variable):
    terms = sp.Poly(poly, variable).terms()
    need(bool(terms), "zero polynomial valuation")
    return min(monomial[0] for monomial, _ in terms)


def normalized_source(X, T, Delta, Theta, Phi, eta):
    P = T + X**2 * T**2
    Y = X * T * P
    K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
    G = sp.expand(
        -X**2 * T / 2
        - 3 * P
        + sp.Rational(8, 3) * P**2
        - sp.Rational(1376, 135) * P**3
        + K * Y**2
        + Phi * P**2 * Y
        + Delta * P**4
        + Theta * P * Y**2
        + eta * P**3 * Y
        - eta * Y**3
    )
    return G, P, Y, K


def reduce_in_x(expr, modulus, X):
    return sp.factor(sp.rem(sp.Poly(sp.cancel(expr), X), sp.Poly(modulus, X)).as_expr())


def main():
    X, T = sp.symbols("X T")
    Delta, Theta, Phi, eta = sp.symbols("Delta Theta Phi eta")
    G, P, Y, K = normalized_source(X, T, Delta, Theta, Phi, eta)
    C = Delta + Theta
    D_A = 4 * Theta * K**2 - 27 * eta**2

    f = sp.cancel(sp.diff(G, X) / T)
    h = sp.diff(G, T)
    need(sp.denom(f) == 1, "G_X/T polynomial")
    f = sp.expand(f)
    h = sp.expand(h)
    need(sp.Poly(f, X).degree() == 7 and sp.Poly(h, X).degree() == 8,
         "row-A X degrees")
    need(sp.factor(sp.Poly(f, X).LC() - 8 * C * T**7) == 0,
         "f leading coefficient")
    need(sp.factor(sp.Poly(h, X).LC() - 8 * C * T**7) == 0,
         "h leading coefficient")

    hessian = sp.Matrix([[sp.diff(G, X, 2), sp.diff(G, X, T)],
                         [sp.diff(G, T, X), sp.diff(G, T, 2)]]).det()
    critical_jacobian = sp.Matrix([[sp.diff(f, X), sp.diff(f, T)],
                                   [sp.diff(h, X), sp.diff(h, T)]]).det()
    bridge = sp.factor(T * critical_jacobian - hessian - f * sp.diff(G, X, T))
    need(bridge == 0, "Keller-Morse bridge identity")

    need(sp.factor(f.subs(T, 0) + X) == 0, "T=0 f row")
    need(sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0, "T=0 h row")
    need(reduce_in_x(hessian.subs(T, 0) - 6, X**2 + 6, X) == 0,
         "T=0 Hessian sign")

    t_universal = -sp.Rational(1, 6)
    need(reduce_in_x(f.subs(T, t_universal), X**2 - 6, X) == 0,
         "universal f pair")
    need(reduce_in_x(h.subs(T, t_universal), X**2 - 6, X) == 0,
         "universal h pair")
    need(reduce_in_x(hessian.subs(T, t_universal) + 6, X**2 - 6, X) == 0,
         "universal Hessian sign")

    # Exact row-A control from THM-4157.  This recomputes the complete
    # X-resultant but deliberately omits every discriminant and collision gate.
    control = {
        Delta: sp.Rational(1),
        Theta: sp.Rational(19, 11),
        Phi: sp.Rational(11, 7),
        eta: sp.Rational(23, 13),
    }
    c0 = sp.factor(C.subs(control))
    da0 = sp.factor(D_A.subs(control))
    need(c0 != 0 and da0 != 0, "control is off the inner-resultant wall")
    f0, h0 = sp.expand(f.subs(control)), sp.expand(h.subs(control))
    resultant = sp.factor(sp.resultant(f0, h0, X))
    need(valuation(resultant, T) == 42, "resultant T artifact")
    residual = sp.cancel(resultant / (T**42 * (6 * T + 1) ** 2))
    need(sp.denom(residual) == 1, "Q19 polynomial")
    q19 = sp.Poly(residual, T)
    need(q19.degree() == 19, "Q19 exact degree")
    need(sp.factor(q19.eval(0) + 12288 * c0**6) == 0, "Q19 constant")
    need(sp.factor(q19.LC() + 1458 * c0 * control[eta]**4 * da0**2) == 0,
         "Q19 leading coefficient")

    residual_length = 2 + q19.degree()
    total_length = residual_length + 2
    need((residual_length, total_length) == (21, 23), "critical length ledger")
    finite_ceiling = 2 * 19 - total_length - 1 + 3
    full_ceiling = 2 * (25 - total_length)
    need(finite_ceiling == 17 and finite_ceiling < 18, "finite response budget")
    need(full_ceiling == 4 and full_ceiling < 18, "full response budget")

    print("symbolic_bridge identity=0 degrees=(7,8) leading=8*C*T^7")
    print("universal_points T0=2 hessian=6 Tm1over6=2 hessian=-6")
    print(f"row_a_control C={c0} D_A={da0} resultant=T^42*(6T+1)^2*Q19 degree={q19.degree()}")
    print("length residual=21 restored=23")
    print("budgets finite=17<18 full=4<18")
    print("forbidden_gates discQ19=unused discR19=unused Q19(-1/6)=unused")
    print(f"checks {CHECKS}")
    print("THM4173_AUDIT_PASS")


if __name__ == "__main__":
    main()
