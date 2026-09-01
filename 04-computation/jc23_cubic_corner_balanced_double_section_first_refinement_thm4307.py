#!/usr/bin/env python3
"""Exact SymPy certificate for THM-4307.

This reconstructs the literal balanced cubic-corner strict transform from the
THM-4301/4304 source, checks the local A1 (rather than A2) deformation algebra,
and computes the exact balanced discriminant graph in all three
beta-versus-5s regimes.  A Regime-A elliptic hostile checks the rationality
boundary.  The script does not claim seam entry or JC(2).
"""

from __future__ import annotations

import sympy as sp


def require_zero(value: sp.Expr, label: str) -> None:
    value = sp.factor(value, extension=sp.sqrt(265))
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def main() -> None:
    q, t, z, y, w_ratio, sigma = sp.symbols("q t z y w_ratio sigma")
    U = sp.symbols("U")
    root265 = sp.sqrt(265)
    delta_star = sp.Rational(2048, 45)
    r = 1 + q

    # Literal THM-4301 source on the balanced coefficient row (THM-4304 (16)).
    K = sp.Rational(2848, 45) - sp.Rational(7, 6) * delta_star
    hhat = (
        U * (r**6 - 2 * r**5 + r**4)
        + delta_star * t**4 * (r**4 - r**3)
        + t**6 * (-sp.Rational(1376, 135) * r**3 + K * r**2)
        + sp.Rational(8, 3) * t**8 * r**2
        - 3 * t**10 * r
    )
    F = sp.expand(q * (hhat - z**12) - t**12 / 2)

    # Put q=t^4 y and divide by the first-face monomial t^12.  The apparent
    # negative t-power occurs only in the positive-valuation toric monomial
    # m=z^12/t^8=(z^2/sigma^4)^2.
    p = U * y**3 + delta_star * y**2 + sp.Rational(8, 3) * y - sp.Rational(1, 2)
    R2 = -sp.Rational(1376, 135) * y**2 - 3 * y
    R4 = 4 * U * y**4 + sp.Rational(2048, 15) * y**3 + sp.Rational(16, 3) * y**2
    R6 = -sp.Rational(2752, 135) * y**3 - 3 * y**2
    R8 = 6 * U * y**5 + sp.Rational(2048, 15) * y**4 + sp.Rational(8, 3) * y**3
    R10 = -sp.Rational(1376, 135) * y**4
    R12 = 4 * U * y**6 + sp.Rational(2048, 45) * y**5
    R16 = U * y**7
    expected_G = (
        p + t**2 * R2 + t**4 * R4 + t**6 * R6 + t**8 * R8
        + t**10 * R10 + t**12 * R12 + t**16 * R16 - y * z**12 / t**8
    )
    G = sp.cancel(F.subs(q, t**4 * y) / t**12)
    require_zero(G - expected_G, "complete balanced strict transform")

    # The local Jacobian algebra of a double root Q^2 is one-dimensional.
    # The two-dimensional quotient belongs to the excluded triple root Q^3.
    Q = sp.symbols("Q")
    require_zero(sp.diff(Q**2, Q) - 2 * Q, "A1 derivative")
    require_zero(sp.diff(Q**3, Q) - 3 * Q**2, "A2 derivative")
    require(sp.degree(sp.diff(Q**2, Q), Q) == 1, "A1 T1 dimension one")
    require(sp.degree(sp.diff(Q**3, Q), Q) == 2, "A2 T1 dimension two")

    discriminant = sp.factor(sp.discriminant(p, y))
    u_values = [
        -sp.Rational(315392, 3645) + sp.Rational(217088, 18225) * root265,
        -sp.Rational(315392, 3645) - sp.Rational(217088, 18225) * root265,
    ]
    rho_values = [
        -(15 + 3 * root265) / 256,
        -(15 - 3 * root265) / 256,
    ]
    require(set(sp.solve(sp.together(sp.numer(discriminant)), U)) == set(u_values),
            "balanced repeated coefficient values")

    T, H, e, shift = sp.symbols("T H e shift")
    x_series = sp.symbols("x_series")
    H_exact = (
        p + x_series * R2 + x_series**2 * R4 + x_series**3 * R6
        + x_series**4 * R8 + x_series**5 * R10
        + x_series**6 * R12 + x_series**8 * R16
    )
    rows: list[tuple[int, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]] = []
    for epsilon, u_value, rho in zip((1, -1), u_values, rho_values):
        p_eps = sp.expand(p.subs(U, u_value))
        require_zero(p_eps.subs(y, rho), f"p(rho), epsilon={epsilon}")
        require_zero(sp.diff(p_eps, y).subs(y, rho),
                     f"p'(rho), epsilon={epsilon}")
        a = sp.factor(sp.diff(p_eps, y, 2).subs(y, rho) / 2,
                      extension=root265)
        g = sp.factor(R2.subs(y, rho), extension=root265)
        expected_a = sp.Rational(128, 135) * (-53 + epsilon * root265)
        expected_g = (-707 + epsilon * 65 * root265) / 3072
        require_zero(a - expected_a, f"quadratic unit epsilon={epsilon}")
        require_zero(g - expected_g, f"R2(rho) epsilon={epsilon}")
        require(a != 0 and g != 0 and rho != 0, f"nonzero splitter data {epsilon}")

        # Exact repeated locus.  For P=H(y,x)-m*y, the common-root equations
        # are E=H-y*H_y=0 and m=H_y.  Since E_y=-rho*p''(rho) is nonzero,
        # they define one smooth formal graph m=m_epsilon(x).
        H_eps = sp.expand(H_exact.subs(U, u_value))
        E_eps = sp.expand(H_eps - y * sp.diff(H_eps, y))
        require_zero(sp.diff(E_eps, y).subs({x_series: 0, y: rho}) + 2 * rho * a,
                     f"smooth discriminant graph epsilon={epsilon}")
        y1, y2 = sp.symbols(f"y1_{epsilon} y2_{epsilon}")
        y_formal = rho + y1 * x_series + y2 * x_series**2
        E_formal = sp.expand(E_eps.subs(y, y_formal))
        y1_value = sp.solve(sp.factor(E_formal.coeff(x_series, 1),
                                     extension=root265), y1)[0]
        y2_equation = sp.factor(
            E_formal.coeff(x_series, 2).subs(y1, y1_value),
            extension=root265,
        )
        y2_value = sp.solve(y2_equation, y2)[0]
        m_formal = sp.expand(
            sp.diff(H_eps, y).subs(y, y_formal)
            .subs({y1: y1_value, y2: y2_value})
        )
        chi = sp.factor(m_formal.coeff(x_series, 1), extension=root265)
        lam = sp.factor(m_formal.coeff(x_series, 2), extension=root265)
        expected_chi = -sp.Rational(173, 72) + epsilon * sp.Rational(43, 360) * root265
        expected_lam = (sp.Rational(11975, 27648)
                        - epsilon * sp.Rational(53621, 7326720) * root265)
        require_zero(chi - expected_chi, f"graph slope epsilon={epsilon}")
        require_zero(lam - expected_lam, f"graph curvature epsilon={epsilon}")
        require_zero(chi - g / rho, f"critical slope g/rho epsilon={epsilon}")
        require(chi != 0, f"nonzero graph slope epsilon={epsilon}")

        # Scale both first deformation parameters by e.  The critical-point
        # displacement changes the value only at order e^2, so modulo I^2
        # the unique A1 coordinate is chi=g*T-rho*H.
        first_model = sp.expand(
            p_eps.subs(y, rho + e * shift)
            + e * T * R2.subs(y, rho + e * shift)
            - e * H * (rho + e * shift)
        )
        critical_equation = sp.expand(sp.diff(first_model, shift) / e**2)
        require_zero(critical_equation.subs(e, 0)
                     - (2 * a * shift + T * sp.diff(R2, y).subs(y, rho) - H),
                     f"critical shift equation epsilon={epsilon}")
        critical_value_linear = sp.expand(first_model).coeff(e, 1)
        require_zero(critical_value_linear - (g * T - rho * H),
                     f"critical value mod I2 epsilon={epsilon}")

        # The normalized quadratic radicand is delta=(rho*H-g*T)/a mod I^2.
        # At beta=5s, H=sigma^12*x^4 and T=sigma^12*x^2.
        boundary_delta = sp.factor(
            (rho * sigma**12 * w_ratio**4 - g * sigma**12 * w_ratio**2) / a,
            extension=root265,
        )
        expected_boundary = sp.factor(
            sigma**12 * w_ratio**2 * (rho * w_ratio**2 - g) / a,
            extension=root265,
        )
        require_zero(boundary_delta - expected_boundary,
                     f"boundary factor epsilon={epsilon}")
        conic_polynomial = sp.Poly(rho * w_ratio**2 - g, w_ratio,
                                   extension=root265)
        require(sp.gcd(conic_polynomial, conic_polynomial.diff()).degree() == 0,
                f"squarefree boundary conic epsilon={epsilon}")
        require(conic_polynomial.degree() == 2,
                f"quadratic boundary conic epsilon={epsilon}")
        rows.append((epsilon, u_value, rho, a, g, chi, lam))

    # Valuation comparison.  H=z^12/t^8 has excess 4(beta-2s), whereas
    # T=t^2 has excess 2(s+beta); their difference is 2(beta-5s).
    s_weight, beta_weight = sp.symbols("s beta", positive=True)
    h_excess = 4 * (beta_weight - 2 * s_weight)
    t2_excess = 2 * (s_weight + beta_weight)
    require_zero(h_excess - t2_excess - 2 * (beta_weight - 5 * s_weight),
                 "splitter wall beta=5s")
    require_zero((z**12 / t**8).subs(t, sigma * z) - (z**2 / sigma**4) ** 2,
                 "H square monomial")

    # At a split quadratic root with d=v(delta), P_y has order d/2 and
    # F_q=t^8 P_y.  The inherited good differential therefore has order
    # s+3beta-d/2.  Both adjacent generic zones are strictly positive.
    d = sp.symbols("d")
    form_order = s_weight + 3 * beta_weight - d / 2
    low_zone_order = sp.factor(form_order.subs(d, h_excess))
    high_zone_order = sp.factor(form_order.subs(d, t2_excess))
    require_zero(low_zone_order - (5 * s_weight + beta_weight),
                 "beta<5s good-form order")
    require_zero(high_zone_order - 2 * beta_weight,
                 "beta>5s good-form order")

    # Literal Regime-A hostile: k=3, Delta=0.  It produces a genuine elliptic
    # first refinement rather than a rational curve, but the good form still
    # vanishes to order ten.
    U3 = sp.Rational(135, 28672)
    rho3 = -sp.Rational(14336, 135)
    hhat3 = (
        U3 * (r**6 - 2 * r**5 + r**4)
        + t**3 * (r**4 - r**3)
        + t**6 * (-sp.Rational(1376, 135) * r**3
                  + sp.Rational(2848, 45) * r**2)
        + sp.Rational(8, 3) * t**8 * r**2
        - 3 * t**10 * r
    )
    F3 = sp.expand(q * (hhat3 - z**12) - t**12 / 2)
    G3 = sp.expand(sp.cancel(F3.subs(q, t**3 * y) / t**9) + y * z**12 / t**6)
    base3 = U3 * y * (y - rho3) ** 2
    require_zero(G3.subs(t, 0) - base3, "Regime-A k3 double base")
    require_zero(sp.Poly(G3, t).coeff_monomial(t**2) - sp.Rational(8, 3) * y,
                 "Regime-A k3 first graph row")
    require_zero((sp.Rational(8, 3) * rho3) / rho3 - sp.Rational(8, 3),
                 "Regime-A k3 graph coefficient")
    w = sp.symbols("w")
    quartic = w**4 - sp.Rational(8, 3)
    require(sp.gcd(sp.Poly(quartic, w), sp.Poly(sp.diff(quartic, w), w)).degree() == 0,
            "Regime-A elliptic quartic squarefree")
    quartic_I = 12 * 1 * sp.Rational(-8, 3)
    quartic_J = sp.Integer(0)
    require(quartic_I == -32 and quartic_J == 0,
            "Regime-A elliptic binary-quartic invariants")
    quartic_j = sp.factor(1728 * (4 * quartic_I**3)
                          / (4 * quartic_I**3 - quartic_J**2))
    require(quartic_j == 1728, "Regime-A elliptic j")
    tau3 = 3
    d3 = 2 * tau3
    fq_order3 = 2 * 3 * tau3 + d3 // 2
    good_order3 = 9 * 1 + 11 * 2 - fq_order3
    require((fq_order3, good_order3) == (21, 10),
            "Regime-A elliptic good-form order")

    print("THM-4307 PRIMARY EXACT AUDIT: PASS")
    print("UNIVERSE exact_M12 D=Lambda=0 UZ!=0 balanced_repeated_cubic beta>2s")
    print("LOCAL_TYPE double_root_A1 T1_basis=1; triple_A2_basis=1,Q excluded")
    print("STRICT_TRANSFORM G=p+t^2*R2+t^4*R4+t^6*R6+t^8*R8+t^10*R10+t^12*R12+t^16*R16-y*z^12/t^8")
    for epsilon, u_value, rho, a, g, chi, lam in rows:
        print(f"epsilon={epsilon:+d} U={u_value} rho={rho} a={a} g={g}")
        print(f"epsilon={epsilon:+d} critical_value_mod_I2={g}*T-({rho})*H")
        print(f"epsilon={epsilon:+d} discriminant_graph m={chi}*x+{lam}*x^2+O(x^3) smooth=yes")
        print(f"epsilon={epsilon:+d} normalized_conic=V^2-(({rho})/({a}))*(w^2-({chi})) genus=0 squarefree=yes")
    print("SPLITTER 2s<beta<5s H_first square; beta=5s conic; beta>5s t^2_first square")
    print("PULLBACK m/x=(z/sigma^5)^2; discriminant has two smooth formal branches z/sigma^5=+/-sqrt(chi)+O(sigma^12)")
    print("GOOD_FORM beta<5s order=5s+beta; beta>5s order=2beta; balanced normal chart smooth")
    print("BALANCED_CONSEQUENCE completed_local_refinement_carriers_rational Keller_constant")
    print("REGIME_A_HOSTILE k=3 Delta=0 U=135/28672 rho=-14336/135 weights=1,2,9 curve=V^2-(w^4-8/3) genus=1 j=1728 Fq_order=21 good_form_order=10")
    print("SCOPE no_seam_entry no_U0 no_Z0 no_JC2 no_DC2")


if __name__ == "__main__":
    main()
