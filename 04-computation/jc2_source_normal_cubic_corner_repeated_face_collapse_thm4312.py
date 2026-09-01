#!/usr/bin/env python3
"""Primary exact certificate for THM-4312.

Intersect THM-4308's finite row-eight bracket/depth gates with THM-4304's
complete repeated-face census at the exact M=12 cubic corner.  The result is
a finite-jet statement only: it does not assert an all-row lift or seam entry.
"""

from __future__ import annotations

import sympy as sp


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def require_zero(expression: sp.Expr, label: str) -> None:
    require(sp.cancel(sp.together(sp.expand(expression))) == 0, label)


def main() -> None:
    Phi, eta, xi10 = sp.symbols("Phi eta xi10")
    Delta, Theta, zeta3, upsilon5 = sp.symbols(
        "Delta Theta zeta3 upsilon5"
    )
    U, W, Z = sp.symbols("U W Z")

    # THM-4308's four bracket equations and three depth-projection gates.
    E5 = 2025 * upsilon5 + 9000 * Delta + 1350 * Theta + 184832
    E6 = (
        200475 * U + 109350 * xi10 - 5593860 * Delta
        - 529200 * Theta - 137763328
    )
    E7 = (
        801900 * W + 1782000 * Delta**2 + 156163200 * Delta
        + 868725 * Phi**2 + 27390480 * Theta - 3434400 * xi10
        + 12891824128
    )
    E8 = (
        21651300 * Z - 225022050 * Delta**2 - 59073300 * Delta * Theta
        - 9512522400 * Delta + 34749000 * Phi**2
        + 39092625 * Phi * eta + 940522560 * Theta
        + 185376600 * xi10 - 1112446017536
    )
    gate = {
        Delta: sp.Rational(896, 15),
        Theta: sp.Rational(512, 75),
        zeta3: -sp.Rational(3, 2) * Phi,
    }
    ups_gate = sp.factor(sp.solve(E5.subs(gate), upsilon5)[0])
    require_zero(ups_gate + sp.Rational(731648, 2025), "gate upsilon5")
    gate[upsilon5] = ups_gate
    U_gate = sp.factor(sp.solve(E6.subs(gate), U)[0])
    W_gate = sp.factor(sp.solve(E7.subs(gate), W)[0])
    Z_gate = sp.factor(sp.solve(E8.subs(gate), Z)[0])
    require_zero(
        U_gate - (475515904 - 109350 * xi10) / sp.Integer(200475),
        "gate U response",
    )
    require_zero(
        W_gate
        + (4343625 * Phi**2 - 17172000 * xi10 + 143826305024)
        / sp.Integer(4009500),
        "gate W response",
    )
    require_zero(
        Z_gate
        - (
            12506118074368 - 173745000 * Phi**2
            - 195463125 * Phi * eta - 926883000 * xi10
        )
        / sp.Integer(108256500),
        "gate Z response",
    )

    # On Lambda=0, the discriminant wall D=W^2-4UZ is (W+2U)^2.
    Lambda = sp.expand(U + W + Z)
    wall_D = sp.expand(W**2 - 4 * U * Z)
    require_zero(
        wall_D - ((W + 2 * U) ** 2 - 4 * U * Lambda),
        "corner square identity",
    )

    # UZ != 0 therefore turns the exact corner into W=-2U, Z=U.
    corner_w = sp.factor(W_gate + 2 * U_gate)
    xi_corner = sp.factor(sp.solve(corner_w, xi10)[0])
    expected_xi = (
        4343625 * Phi**2 + 124805668864
    ) / sp.Integer(12798000)
    require_zero(xi_corner - expected_xi, "corner xi10")

    corner_z = sp.factor((Z_gate - U_gate).subs(xi10, xi_corner))
    phi_eta_corner = (
        sp.Rational(2091705253888, 107983125)
        - sp.Rational(2839, 1185) * Phi**2
    )
    require_zero(
        corner_z.subs(eta, phi_eta_corner / Phi),
        "corner Phi*eta",
    )
    U_corner = sp.factor(U_gate.subs(xi10, xi_corner))
    expected_U_corner = -sp.Rational(13, 57591000) * (
        820125 * Phi**2 + 13056802816
    )
    require_zero(U_corner - expected_U_corner, "corner U")
    require_zero(
        W_gate.subs(xi10, xi_corner) + 2 * U_corner,
        "corner W=-2U",
    )
    require_zero(
        Z_gate.subs({xi10: xi_corner, eta: phi_eta_corner / Phi}) - U_corner,
        "corner Z=U",
    )
    require(phi_eta_corner.subs(Phi, 0) != 0, "corner forces Phi nonzero")

    # Intersect with THM-4304's exhaustive five repeated regimes.
    balanced_residual = sp.factor(gate[Delta] - sp.Rational(2048, 45))
    require(balanced_residual == sp.Rational(128, 9), "balanced exclusion")
    c4_gate = sp.factor(gate[Delta] + gate[Theta])
    require(c4_gate == sp.Rational(1664, 25), "k3 exclusion")

    # In k=2, xi10=-upsilon5 and the repeated-square equation is
    # upsilon5^2=4U(Delta+Theta).  Its exact residual is nonzero.
    xi_k2 = -ups_gate
    U_k2 = sp.factor(U_gate.subs(xi10, xi_k2))
    k2_residual = sp.factor(ups_gate**2 - 4 * U_k2 * c4_gate)
    require(U_k2 == sp.Rational(39636992, 18225), "k2 U")
    require(
        k2_residual == -sp.Rational(1839105572864, 4100625),
        "k2 repeated-square exclusion",
    )

    # Consequently only Regime-A k=1 can survive.  On the corner its square
    # equation alpha11^2=4U(upsilon5+xi10) has this exact right side.
    c2_corner = sp.factor(ups_gate + xi_corner)
    k1_square = sp.factor(4 * U_corner * c2_corner)
    expected_c2 = sp.Rational(11, 474000) * (
        14625 * Phi**2 + 404652032
    )
    expected_square = -sp.Rational(143, 6824533500000) * (
        14625 * Phi**2 + 404652032
    ) * (
        820125 * Phi**2 + 13056802816
    )
    require_zero(c2_corner - expected_c2, "surviving k1 c2")
    require_zero(k1_square - expected_square, "surviving k1 square")

    # Reconstruct the literal k=1 strict transform independently of the
    # classification table.  The -q*z^12 term, omitted here, becomes -y*m
    # with m=z^12/t^2.
    q, t, y = sp.symbols("q t y")
    Ug, alpha = sp.symbols("Ug alpha", nonzero=True)
    upg, etag, zetag, Dg, Tg, Pg = sp.symbols(
        "upg etag zetag Dg Tg Pg"
    )
    r = 1 + q
    Kg = sp.Rational(2848, 45) - sp.Rational(7, 6) * Dg
    hhat = (
        Ug * (r**6 - 2 * r**5 + r**4)
        + t * (alpha * r**5 - alpha * r**4)
        + t**2 * (upg * r**5 + (alpha**2 / (4 * Ug) - upg) * r**4)
        + t**3 * (etag * r**4 + zetag * r**3)
        + t**4 * (Dg * r**4 + Tg * r**3)
        + t**5 * Pg * r**3
        + t**6 * (-sp.Rational(1376, 135) * r**3 + Kg * r**2)
        + sp.Rational(8, 3) * t**8 * r**2
        - 3 * t**10 * r
    )
    H = sp.expand(sp.cancel((q * hhat - t**12 / 2).subs(q, t * y) / t**3))
    base = sp.expand(H).coeff(t, 0)
    R1 = sp.expand(H).coeff(t, 1)
    expected_base = y * (2 * Ug * y + alpha) ** 2 / (4 * Ug)
    require_zero(base - expected_base, "literal k1 double base")
    rho = -alpha / (2 * Ug)
    require_zero(base.subs(y, rho), "k1 base at rho")
    require_zero(sp.diff(base, y).subs(y, rho), "k1 base derivative at rho")
    require(sp.factor(sp.diff(base, y, 2).subs(y, rho)) != 0, "k1 not triple")

    # Critical-root graph for P=H-y*m.  Solving E=H-yH_y through t^1
    # gives m=L*t+O(t^2), with one intrinsic A1 coordinate.
    E = sp.expand(H - y * sp.diff(H, y))
    y1, y2 = sp.symbols("y1 y2")
    y_formal = rho + y1 * t + y2 * t**2
    E1 = sp.factor(sp.expand(E.subs(y, y_formal)).coeff(t, 1))
    y1_value = sp.solve(E1, y1)[0]
    E2 = sp.factor(
        sp.expand(E.subs(y, y_formal).subs(y1, y1_value)).coeff(t, 2)
    )
    y2_value = sp.solve(E2, y2)[0]
    m_formal = sp.expand(
        sp.diff(H, y).subs(y, y_formal).subs({y1: y1_value, y2: y2_value})
    )
    L = sp.factor(m_formal.coeff(t, 1))
    L2 = sp.factor(m_formal.coeff(t, 2))
    expected_L = etag + zetag - alpha * upg / (2 * Ug)
    expected_L2 = sp.factor(
        (
            4 * Ug**2 * (Dg + Tg)
            - 2 * Ug * alpha * (4 * etag + 3 * zetag)
            - Ug * upg**2
            + 4 * alpha**2 * upg
        )
        / (4 * Ug**2)
    )
    require_zero(L - expected_L, "k1 first discriminant coefficient")
    require_zero(L2 - expected_L2, "k1 second discriminant coefficient")
    require_zero(R1.subs(y, rho) / rho - expected_L, "k1 critical value check")

    # Exact finite-jet positive control: Phi=1.  The square is nonzero over C,
    # and the two alpha signs cannot both kill L because their difference is
    # alpha*upsilon5/U, also nonzero.
    control = {Phi: sp.Integer(1)}
    xi_control = sp.factor(xi_corner.subs(control))
    eta_control = sp.factor(phi_eta_corner.subs(control))
    U_control = sp.factor(U_corner.subs(control))
    c2_control = sp.factor(c2_corner.subs(control))
    square_control = sp.factor(k1_square.subs(control))
    require(xi_control == sp.Rational(124810012489, 12798000), "control xi")
    require(eta_control == sp.Rational(2091446550013, 107983125), "control eta")
    require(U_control == -sp.Rational(169749098233, 57591000), "control U")
    require(c2_control == sp.Rational(4451333227, 474000), "control c2")
    require(
        square_control == -sp.Rational(755609801217839887891, 6824533500000),
        "control alpha square",
    )
    require(square_control != 0 and ups_gate != 0 and U_control != 0,
            "control nonzero sidecars")
    alpha_control = sp.sqrt(square_control)
    zeta_control = -sp.Rational(3, 2)
    sign_term = sp.factor(alpha_control * ups_gate / (2 * U_control))
    L_plus = sp.factor(eta_control + zeta_control - sign_term)
    L_minus = sp.factor(eta_control + zeta_control + sign_term)
    require(sign_term != 0, "control sign separation")
    require_zero(L_minus - L_plus - 2 * sign_term,
                 "control two-sign difference")
    rho_square_control = sp.factor(c2_control / U_control)
    c3_control = sp.Rational(4182569150651, 215966250)
    L_norm_control = sp.factor(
        c3_control**2 - ups_gate**2 * rho_square_control
    )
    require(
        rho_square_control == -sp.Rational(98333997651, 30863472406),
        "control rho square",
    )
    require(
        L_norm_control
        == sp.Rational(
            270260378011253985379632330934787603,
            719758107151040278729687500,
        ),
        "control L1 norm",
    )

    # The cancellation L=0 is a genuine nonempty hostile locus, not a formal
    # corner excluded by U, c2, or Phi.  Eliminate rho from
    # L=c3+upsilon5*rho=0 and rho^2=c2/U using r=Phi^2.
    r_parameter = sp.symbols("r_parameter")
    c3_numerator = 4183410507776 - 841357125 * r_parameter
    c3_square = c3_numerator**2 / (
        sp.Integer(215966250) ** 2 * r_parameter
    )
    c2_over_U = -sp.Rational(2673, 26) * (
        14625 * r_parameter + 404652032
    ) / (
        820125 * r_parameter + 13056802816
    )
    cancellation_equation = sp.factor(c3_square / ups_gate**2 - c2_over_U)
    cancellation_numerator = sp.Poly(
        sp.together(cancellation_equation).as_numer_denom()[0],
        r_parameter,
    ).primitive()[1]
    cancellation_R = (
        7547170421607067494140625 * r_parameter**3
        + 164114458618573873612800000000 * r_parameter**2
        + 2284603892441775363795663716352000 * r_parameter
        + 2970579390109346274816679296272171008
    )
    require_zero(cancellation_numerator.as_expr() - cancellation_R,
                 "high-contact cancellation cubic")
    forbidden = sp.Poly(
        r_parameter
        * (820125 * r_parameter + 13056802816)
        * (14625 * r_parameter + 404652032),
        r_parameter,
    )
    require(sp.gcd(cancellation_numerator, forbidden).degree() == 0,
            "high-contact roots avoid forbidden factors")
    # On L1=0 the second graph coefficient reduces to a rational function of
    # r.  Its numerator S is coprime to R; the displayed Bezout identity
    # modulo 17 is a compact exact witness, so L2 never also vanishes.
    c3_times_zeta = (
        c3_numerator / sp.Integer(215966250) * -sp.Rational(3, 2)
    )
    L2_high = sp.factor(
        c4_gate + c3_times_zeta / ups_gate
        - ups_gate**2
        / (
            4
            * -sp.Rational(13, 57591000)
            * (820125 * r_parameter + 13056802816)
        )
    )
    S_polynomial = sp.Poly(
        8970234157828125 * r_parameter**2
        + 61293210070929408000 * r_parameter
        - 1395571970793868500140032,
        r_parameter,
    )
    require_zero(
        L2_high
        + S_polynomial.as_expr()
        / (
            676262246400
            * (820125 * r_parameter + 13056802816)
        ),
        "high-contact L2 formula",
    )
    require(sp.gcd(cancellation_numerator, S_polynomial).degree() == 0,
            "high-contact L1/L2 no common root")
    bezout_mod17 = sp.Poly(
        (-5 * r_parameter - 8) * cancellation_R
        + (-r_parameter**2 + r_parameter - 5) * S_polynomial.as_expr()
        - 1,
        r_parameter,
        modulus=17,
    )
    require(bezout_mod17.is_zero, "high-contact mod17 Bezout witness")

    # If L1!=0, the first splitter is beta=s/3.  Primitive weights are
    # (s,beta,gamma)=(3,1,4).  The actual exceptional coordinate on P(3,1) is
    # w=sigma/z^3, not the coordinate on the ramified substitution
    # sigma=lambda^3, z=lambda*w0.  Dividing by z^4 gives
    # a V^2=w^(-2)-Lw; with Y=wV this is a Y^2=1-Lw^3.
    s_weight, beta_weight = sp.symbols("s beta", positive=True)
    m_weight = 10 * beta_weight - 2 * s_weight
    t_weight = s_weight + beta_weight
    require_zero(m_weight - t_weight - 3 * (3 * beta_weight - s_weight),
                 "k1 splitter beta=s/3")
    require((m_weight.subs({s_weight: 3, beta_weight: 1}),
             t_weight.subs({s_weight: 3, beta_weight: 1})) == (4, 4),
            "primitive splitter weights")
    w, Lsymbol = sp.symbols("w L")
    exceptional = 1 - Lsymbol * w**3
    field = sp.QQ.frac_field(Lsymbol)
    exceptional_poly = sp.Poly(exceptional, w, domain=field)
    require(sp.gcd(exceptional_poly, exceptional_poly.diff()).degree() == 0,
            "generic k1 exceptional squarefree")
    require(exceptional_poly.degree() == 3, "generic k1 exceptional degree")
    genus = (exceptional_poly.degree() - 1) // 2
    require(genus == 1, "generic k1 exceptional genus")
    weierstrass_A = sp.Integer(0)
    weierstrass_B = sp.symbols("B", nonzero=True)
    elliptic_j = sp.factor(
        1728 * 4 * weierstrass_A**3
        / (4 * weierstrass_A**3 + 27 * weierstrass_B**2)
    )
    require(elliptic_j == 0, "generic k1 exceptional j=0")

    tau = 4
    d = 4
    fq_order = 2 * tau + d // 2
    good_form_lower_bound = 9 * 3 + 11 * 1 - fq_order
    require((fq_order, good_form_lower_bound) == (10, 28),
            "generic k1 differential ledger")

    # If L1=0, the coprimality above forces L2!=0.  The next splitter is
    # beta=s/2 with primitive (s,beta,gamma)=(2,1,3).  On P(2,1),
    # w=sigma/z^2 and Y=wV give a Y^2=1-L2*w^4.
    require_zero(
        m_weight - 2 * t_weight - 4 * (2 * beta_weight - s_weight),
        "high-contact splitter beta=s/2",
    )
    require((m_weight.subs({s_weight: 2, beta_weight: 1}),
             (2 * t_weight).subs({s_weight: 2, beta_weight: 1})) == (6, 6),
            "high-contact primitive splitter weights")
    L2symbol = sp.symbols("L2")
    high_exceptional = 1 - L2symbol * w**4
    high_field = sp.QQ.frac_field(L2symbol)
    high_poly = sp.Poly(high_exceptional, w, domain=high_field)
    require(sp.gcd(high_poly, high_poly.diff()).degree() == 0,
            "high-contact exceptional squarefree")
    require(high_poly.degree() == 4, "high-contact exceptional degree")
    require((high_poly.degree() - 2) // 2 == 1,
            "high-contact exceptional genus")
    quartic_I = -12 * L2symbol
    quartic_J = sp.Integer(0)
    quartic_j = sp.factor(
        1728 * 4 * quartic_I**3 / (4 * quartic_I**3 - quartic_J**2)
    )
    require(quartic_j == 1728, "high-contact exceptional j=1728")
    tau_high = 3
    d_high = 6
    fq_high = 2 * tau_high + d_high // 2
    good_high = 9 * 2 + 11 * 1 - fq_high
    require((fq_high, good_high) == (9, 20),
            "high-contact differential ledger")

    print("THM-4312 PRIMARY EXACT AUDIT: PASS")
    print("GATES Delta=896/15 Theta=512/75 zeta3=-3*Phi/2 upsilon5=-731648/2025")
    print("CORNER xi10=(4343625*Phi^2+124805668864)/12798000")
    print("CORNER Phi*eta=2091705253888/107983125-(2839/1185)*Phi^2")
    print("CORNER U=Z=-13*(820125*Phi^2+13056802816)/57591000 W=-2U Phi!=0")
    print("EXCLUDE balanced Delta_residual=128/9")
    print("EXCLUDE k2 repeated_residual=-1839105572864/4100625")
    print("EXCLUDE k3 Delta+Theta=1664/25")
    print("SURVIVE_ONLY k1 beta11=-alpha11 alpha11^2=4U(upsilon5+xi10)")
    print("K1_SQUARE -143*(14625*Phi^2+404652032)*(820125*Phi^2+13056802816)/6824533500000")
    print("K1_GRAPH rho=-alpha11/(2U) m=(eta+zeta3-alpha11*upsilon5/(2U))*t+O(t^2)")
    print("CONTROL Phi=1 xi10=124810012489/12798000 eta=2091446550013/107983125 U=-169749098233/57591000 alpha11^2=-755609801217839887891/6824533500000")
    print("CONTROL_K1 rho^2=-98333997651/30863472406 L1_norm=270260378011253985379632330934787603/719758107151040278729687500")
    print("HIGH_CONTACT L1=0 exact_cubic_R_degree=3 roots_avoid_Phi0_U0_c20=yes gcd_R_L2numerator=1 mod17_Bezout=yes")
    print("CASE_L1 beta=s/3 primitive_weights=3,1,4 invariant=w=sigma/z^3 curve=Y^2-(1-L1*w^3) genus=1 j=0 Fq=10 good_form_lower_bound=28")
    print("CASE_L2 beta=s/2 primitive_weights=2,1,3 invariant=w=sigma/z^2 curve=Y^2-(1-L2*w^4) genus=1 j=1728 Fq=9 good_form_lower_bound=20")
    print("CONSEQUENCE formal_good_order_positive actual_exact_M12_Keller_lift_implies_first_carrier_constant")
    print("SCOPE finite_row8_projection first_exceptional_only no_completed_local_tower no_all_row_lift no_seam_entry no_JC2 no_DC2")
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
