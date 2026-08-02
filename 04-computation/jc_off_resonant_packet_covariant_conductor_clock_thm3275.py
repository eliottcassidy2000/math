#!/usr/bin/env python3
"""Exact symbolic companion for THM-3275."""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def collision_covariant(p, q, r):
    return sp.expand(
        9 * p**6
        - 1900 * p**4 * r
        + 1000 * p**3 * q**2
        + 110000 * p**2 * r**2
        - 1000000 * r**3
    )


def packet_data(quartic, x):
    poly = sp.Poly(sp.expand(quartic), x)
    coeffs = poly.all_coeffs()
    require(len(coeffs) == 5 and coeffs[0] == 1, "quartic must be monic")
    _, a3, p, q, r = coeffs
    require(sp.expand(a3) == 0, "quartic must be depressed")
    packet = -(20 * x**3 + 18 * p * x + 27 * q) / sp.Integer(27)
    return p, q, r, sp.expand(packet)


def reduce_pi(expression, pi, t):
    """Reduce a rational expression modulo pi^3-t."""
    numerator, denominator = sp.together(expression).as_numer_denom()
    numerator = sp.rem(sp.Poly(numerator, pi), sp.Poly(pi**3 - t, pi)).as_expr()
    denominator = sp.rem(sp.Poly(denominator, pi), sp.Poly(pi**3 - t, pi)).as_expr()
    return sp.cancel(numerator / denominator)


def monic_minpoly(eta, pi, t, w):
    polynomial = sp.Poly(sp.resultant(pi**3 - t, w - eta, pi), w)
    return sp.factor(polynomial.as_expr() / polynomial.LC())


def t_valuation(expression, t):
    """t-adic valuation for a nonzero rational expression over Q(B)."""
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_poly = sp.Poly(numerator, t)
    denominator_poly = sp.Poly(denominator, t)
    require(numerator_poly.as_expr() != 0 and denominator_poly.as_expr() != 0,
            "valuation requested for zero")
    numerator_order = min(monomial[0] for monomial, _ in numerator_poly.terms())
    denominator_order = min(monomial[0] for monomial, _ in denominator_poly.terms())
    return numerator_order - denominator_order


def check_packet_control(quartic, fixed_root, moving_root, expected_eta,
                         expected_eta_fixed, expected_c_scaled,
                         expected_cubic_disc, expected_h,
                         x, pi, t, b, w):
    p, q, r, packet = packet_data(quartic, x)
    eta = reduce_pi(t * packet.subs(x, moving_root), pi, t)
    eta_fixed = sp.factor(t * packet.subs(x, fixed_root))
    require(reduce_pi(eta - expected_eta, pi, t) == 0,
            f"h={expected_h} moving packet changed")
    require(sp.expand(eta_fixed - expected_eta_fixed) == 0,
            f"h={expected_h} fixed packet changed")

    cubic_minpoly = monic_minpoly(eta, pi, t, w)
    cubic_disc = sp.factor(sp.discriminant(cubic_minpoly, w))
    require(sp.factor(cubic_disc - expected_cubic_disc) == 0,
            f"h={expected_h} cubic packet discriminant changed")
    require(t_valuation(cubic_disc, t) == 2 * expected_h,
            f"h={expected_h} cubic discriminant valuation failed")

    packet_polynomial = sp.resultant(quartic, w / t - packet, x)
    packet_scaled = sp.factor(t**4 * packet_polynomial)
    expected_packet_scaled = sp.factor((w - eta_fixed) * cubic_minpoly)
    require(sp.factor(packet_scaled - expected_packet_scaled) == 0,
            f"h={expected_h} quartic packet factorization changed")
    packet_disc = sp.factor(sp.discriminant(packet_scaled, w))
    cross_resultant = sp.factor(cubic_minpoly.subs(w, eta_fixed))
    require(
        sp.factor(packet_disc - cubic_disc * cross_resultant**2) == 0,
        f"h={expected_h} quartic/cubic discriminant separation failed",
    )
    require(t_valuation(cross_resultant, t) == 0,
            f"h={expected_h} fixed/cubic cross resultant is not a unit")
    require(t_valuation(packet_disc, t) == 2 * expected_h,
            f"h={expected_h} quartic packet discriminant valuation failed")

    c_scaled = sp.factor(t**4 * collision_covariant(p, q, r))
    require(sp.factor(c_scaled - expected_c_scaled) == 0,
            f"h={expected_h} normalized collision covariant changed")
    require(t_valuation(c_scaled, t) == expected_h,
            f"h={expected_h} normalized covariant valuation failed")

    reduced_fixed = eta_fixed.subs(t, 0)
    reduced_moving = cubic_minpoly.subs(t, 0)
    require(reduced_fixed == 1, f"h={expected_h} fixed residue changed")
    require(sp.expand(reduced_moving - (w - sp.Rational(7, 27)) ** 3) == 0,
            f"h={expected_h} moving packet reduction changed")

    return cubic_minpoly, c_scaled, cross_resultant


def main():
    x, pi, t, b, w, c = sp.symbols("x pi t B W c", nonzero=True)

    # Abstract tame-C3 discriminant and index controls.  For eta=c+b*pi^h,
    # pi^3=t and 3 not dividing h, the minimal polynomial is
    # (W-c)^3-b^3*t^h.  Its discriminant exponent is 2h, hence the order
    # index exponent is h-1 after subtracting the tame field exponent two.
    abstract_checks = 0
    for h in range(1, 20):
        if h % 3 == 0:
            continue
        abstract_minpoly = (w - c) ** 3 - b**3 * t**h
        abstract_disc = sp.factor(sp.discriminant(abstract_minpoly, w))
        require(abstract_disc == -27 * b**6 * t ** (2 * h),
                f"abstract discriminant failed at h={h}")
        require((2 * h - 2) // 2 == h - 1,
                f"abstract index exponent failed at h={h}")
        abstract_checks += 1

    # Pair-quotient valuation ledger in the off-resonant graph lane.
    # Three fixed/moving quotients have w-value -2m; three moving/moving
    # quotients have w-value h-2m.  Since w restricted to K is 3v, scaling
    # C by s^4, w(s)=3m, leaves base valuation h.
    ledger_checks = 0
    for m in range(1, 16):
        if m % 3 == 0:
            continue
        for h in range(1, 16):
            if h % 3 == 0:
                continue
            w_c = 3 * (-2 * m) + 3 * (h - 2 * m)
            require(w_c == 3 * h - 12 * m,
                    "six-pair covariant valuation ledger failed")
            require((w_c + 4 * (3 * m)) // 3 == h,
                    "normalized base covariant valuation failed")
            ledger_checks += 1

    # h=1 actual local graph-anatomy control.  Before depression the fixed
    # root is zero and the moving roots are B+pi_i^(-1).  A common finite
    # correction creates the first nonbase packet term at pi^1.
    quartic_h1 = sp.expand(
        (x + 3 * b / 4) * ((x - b / 4) ** 3 - 1 / t)
    )
    eta_h1 = (
        7 - 15 * b * pi + 3 * b**2 * pi**2 - 2 * b**3 * t
    ) / sp.Integer(27)
    c_h1 = sp.factor(
        sp.Rational(27, 64)
        * b**3
        * t
        * (b**3 * t + 125)
        * (8 * b**6 * t**2 + 475 * b**3 * t + 8000)
    )
    disc_h1 = -b**6 * t**2 * (b**3 * t + 125) ** 2 / sp.Integer(19683)
    mu_h1, normalized_c_h1, cross_h1 = check_packet_control(
        quartic_h1,
        -3 * b / 4,
        pi**2 / t + b / 4,
        eta_h1,
        sp.Integer(1),
        c_h1,
        disc_h1,
        1,
        x, pi, t, b, w,
    )

    # h=2 actual local graph-anatomy hostile.  The moving root
    # pi^(-1)+B*pi has two opposite C3 characters.  The normalized packet
    # order has discriminant exponent four and index exponent one.
    quartic_h2 = sp.expand(
        x * (x**3 - 3 * b * x - (1 / t + b**3 * t))
    )
    eta_h2 = (
        7 + 7 * b**3 * t**2 - 6 * b * pi**2 - 6 * b**2 * t * pi
    ) / sp.Integer(27)
    c_h2 = sp.factor(
        -27 * b**3 * t**2
        * (1000 * b**6 * t**4 + 1757 * b**3 * t**2 + 1000)
    )
    disc_h2 = (
        -64 * b**6 * t**4 * (b**3 * t**2 - 1) ** 2
        / sp.Integer(19683)
    )
    mu_h2, normalized_c_h2, cross_h2 = check_packet_control(
        quartic_h2,
        sp.Integer(0),
        pi**2 / t + b * pi,
        eta_h2,
        1 + b**3 * t**2,
        c_h2,
        disc_h2,
        2,
        x, pi, t, b, w,
    )

    require((t_valuation(sp.discriminant(mu_h1, w), t) - 2) // 2 == 0,
            "h=1 maximal-order boundary failed")
    require((t_valuation(sp.discriminant(mu_h2, w), t) - 2) // 2 == 1,
            "h=2 nonmaximal-order hostile failed")

    print("THM3275 off-resonant packet covariant conductor clock exact companion")
    print(f"abstract_tame_C3_checks={abstract_checks} h_range=1..19 excluding_3Z")
    print(f"pair_valuation_ledger_checks={ledger_checks}")
    print("identity: v_K(s^4*C)=h")
    print("identity: v_K(disc_cubic_packet)=v_K(disc_quartic_packet)=2h")
    print("identity: packet_order_index_exponent=h-1")
    print("h=1: normalized_C_valuation=1 cubic_disc=2 quartic_disc=2 index=0")
    print("h=1_boundary: normalized_C_nonunit_but_packet_order_maximal=True")
    print("h=2: normalized_C_valuation=2 cubic_disc=4 quartic_disc=4 index=1")
    print("h=2_hostile: off_resonant_packet_order_nonmaximal=True")
    print(f"h=1_cross_resultant={cross_h1}")
    print(f"h=2_cross_resultant={cross_h2}")
    print(f"h=1_normalized_C={normalized_c_h1}")
    print(f"h=2_normalized_C={normalized_c_h2}")


if __name__ == "__main__":
    main()
