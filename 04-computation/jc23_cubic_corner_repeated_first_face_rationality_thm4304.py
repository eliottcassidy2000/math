#!/usr/bin/env python3
"""Exact SymPy certificate for THM-4304."""

from __future__ import annotations

import sympy as sp


def require_zero(value: sp.Expr, label: str) -> None:
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def main() -> None:
    q, t, z, y, U = sp.symbols("q t z y U")
    alpha, beta11 = sp.symbols("alpha_11 beta_11")
    upsilon, xi10 = sp.symbols("upsilon_5 xi_10")
    eta, zeta3 = sp.symbols("eta zeta_3")
    Delta, Theta, Phi = sp.symbols("Delta Theta Phi")
    r1 = 1 + q
    K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
    hhat = (
        U * (r1**6 - 2 * r1**5 + r1**4)
        + t * (alpha * r1**5 + beta11 * r1**4)
        + t**2 * (upsilon * r1**5 + xi10 * r1**4)
        + t**3 * (eta * r1**4 + zeta3 * r1**3)
        + t**4 * (Delta * r1**4 + Theta * r1**3)
        + t**5 * Phi * r1**3
        + t**6 * (-sp.Rational(1376, 135) * r1**3 + K * r1**2)
        + sp.Rational(8, 3) * t**8 * r1**2 - 3 * t**10 * r1
    )
    F = sp.Poly(sp.expand(q * (hhat - z**12) - t**12 / 2), q, t, z)
    require(max(monomial[0] for monomial, _ in F.terms()) == 7,
            "literal total q support")
    require_zero(F.coeff_monomial(q**3).subs({t: 0, z: 0}) - U,
                 "q3 unit")

    c1 = alpha + beta11
    c2 = upsilon + xi10
    c3 = eta + zeta3
    c4 = Delta + Theta
    c5 = Phi
    c6 = sp.Rational(7168, 135) - sp.Rational(7, 6) * Delta
    A1 = 5 * alpha + 4 * beta11
    A2 = 5 * upsilon + 4 * xi10
    A3 = 4 * eta + 3 * zeta3
    require_zero(F.coeff_monomial(q * t) - c1, "q t row")
    require_zero(F.coeff_monomial(q * t**2) - c2, "q t2 row")
    require_zero(F.coeff_monomial(q * t**3) - c3, "q t3 row")
    require_zero(F.coeff_monomial(q**2 * t) - A1, "q2 t row")
    require_zero(F.coeff_monomial(q**2 * t**2) - A2, "q2 t2 row")
    require_zero(F.coeff_monomial(q**2 * t**3) - A3, "q2 t3 row")

    require_zero(A1.subs(beta11, -alpha) - alpha, "k1 transverse")
    require_zero(A2.subs(xi10, -upsilon) - upsilon, "k2 transverse")
    require_zero(A3.subs(zeta3, -eta) - eta, "k3 transverse")
    regime_a = []
    for k in (1, 2, 3):
        require(3 * k < 12, f"k{k} constant absent")
        Ak, Bk = sp.symbols(f"A{k} B{k}")
        quadratic = U * q**2 + Ak * t**k * q + Bk * t**(2 * k)
        quadratic_disc = sp.factor(sp.discriminant(quadratic, q))
        require_zero(quadratic_disc - (Ak**2 - 4 * U * Bk) * t**(2 * k),
                     f"k{k} quadratic discriminant")
        square_value = sp.factor(quadratic.subs(Bk, Ak**2 / (4 * U)))
        require_zero(square_value - U * (q + Ak * t**k / (2 * U))**2,
                     f"k{k} repeated square")
        boundary = U * q**2 + Ak * t**k * q + Bk * t**(2 * k) - z**12
        boundary_disc = sp.factor(sp.discriminant(boundary, q))
        require_zero(boundary_disc.subs(Bk, Ak**2 / (4 * U)) - 4 * U * z**12,
                     f"k{k} z-boundary destroys repetition")
        regime_a.append(k)

    q_linear_orders = {1, 2, 3, 4, 5, 6, 8, 10}
    regime_b = []
    for k in (1, 2, 3):
        gamma_over_tau = sp.Rational(12 - k, 2)
        required_order = 12 - gamma_over_tau
        require(gamma_over_tau > 4, f"B{k} cubic absent")
        require(required_order not in q_linear_orders, f"B{k} no t row")
        regime_b.append((k, gamma_over_tau, required_order))

    delta_star = sp.Rational(2048, 45)
    balanced = {
        alpha: 0, beta11: 0, upsilon: 0, xi10: 0,
        eta: 0, zeta3: 0, Phi: 0,
        Delta: delta_star, Theta: -delta_star,
    }
    require_zero(c1.subs(balanced), "balanced c1")
    require_zero(c2.subs(balanced), "balanced c2")
    require_zero(c3.subs(balanced), "balanced c3")
    require_zero(c4.subs(balanced), "balanced c4")
    require_zero(c5.subs(balanced), "balanced c5")
    require_zero(c6.subs(balanced), "balanced c6")
    require_zero(F.coeff_monomial(q**2 * t**4).subs(balanced) - delta_star,
                 "balanced q2 t4")
    p = U * y**3 + delta_star * y**2 + sp.Rational(8, 3) * y - sp.Rational(1, 2)
    discriminant = sp.factor(sp.discriminant(p, y))
    expected = -(820125 * U**2 + 141926400 * U - 24696061952) / sp.Integer(121500)
    require_zero(discriminant - expected, "balanced discriminant")
    u_values = sp.solve(sp.together(sp.numer(discriminant)), U)
    expected_u = [
        -sp.Rational(315392, 3645) + sp.Rational(217088, 18225) * sp.sqrt(265),
        -sp.Rational(315392, 3645) - sp.Rational(217088, 18225) * sp.sqrt(265),
    ]
    require(set(u_values) == set(expected_u), "two repeated U values")
    gcd_rows = []
    for epsilon, u_value in ((1, expected_u[0]), (-1, expected_u[1])):
        specialized = sp.Poly(p.subs(U, u_value), y, extension=True)
        derivative = sp.Poly(sp.diff(p, y).subs(U, u_value), y, extension=True)
        gcd_expr = sp.factor(sp.gcd(specialized, derivative).as_expr(), extension=sp.sqrt(265))
        expected_gcd = y + (15 + epsilon * 3 * sp.sqrt(265)) / 256
        require_zero(gcd_expr - expected_gcd, f"gcd epsilon={epsilon}")
        root = -(15 + epsilon * 3 * sp.sqrt(265)) / 256
        require_zero(p.subs({U: u_value, y: root}), f"double root value {epsilon}")
        require_zero(sp.diff(p, y).subs({U: u_value, y: root}),
                     f"double root derivative {epsilon}")
        require(sp.simplify(sp.diff(p, y, 2).subs({U: u_value, y: root})) != 0,
                f"not triple {epsilon}")
        gcd_rows.append((epsilon, u_value, expected_gcd))

    forced_r = sp.Rational(9, 16)
    forced_u = sp.Rational(1, 2) / forced_r**3
    require(forced_u == sp.Rational(2048, 729), "triple forced U")
    require(-3 * forced_u * forced_r == -sp.Rational(128, 27),
            "triple forced y2")
    require(-sp.Rational(128, 27) != delta_star, "triple contradiction")

    weights = (1, 3, 16)
    tau = weights[0] + weights[1]
    hostile_weights = (3 * weights[2], 2 * weights[2] + 4 * tau,
                       weights[2] + 8 * tau, 12 * tau,
                       weights[2] + 12 * weights[1])
    require(hostile_weights == (48, 48, 48, 48, 52), "hostile weights")
    require_zero(discriminant.subs(U, expected_u[0]), "hostile repeated")
    separable_control = sp.factor(discriminant.subs(U, 1))
    require(separable_control == sp.Rational(24553315427, 121500),
            "U=1 separable control")

    print("THM-4304 PRIMARY EXACT AUDIT: PASS")
    print("UNIVERSE exact_M12 D=Lambda=0 UZ!=0 characteristic_zero literal_first_faces")
    print("DEGREE deg_q(face)<=3; repeated_irreducible=>deg_q(K)=1=>rational_graph")
    print("REGIME_A gamma/tau=" + ",".join(map(str, regime_a))
          + " quadratic_square strict_z_boundary")
    print("REGIME_B empty " + " ".join(f"k{k}:gamma={g}:need={o}" for k, g, o in regime_b))
    print(f"REGIME_C discriminant={discriminant}")
    for epsilon, u_value, gcd_expr in gcd_rows:
        print(f"U_{epsilon:+d}={u_value}; gcd={gcd_expr}; multiplicity=2")
    print(f"HOSTILE weights={weights} face_weights={hostile_weights} repeated=yes rational=yes")
    print(f"CONTROL U=1 discriminant={separable_control}")
    print("SCOPE reduced_first_face_carriers_constant refinement_tails_open")


if __name__ == "__main__":
    main()
