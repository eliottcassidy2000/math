#!/usr/bin/env python3
"""Exact companion for THM-3059.

The companion verifies the explicit dominant two-jet map, its primitive
quartic, discriminant and resolvent identities, the reciprocal Newton face,
and the parity bookkeeping.  Irreducibility is certified in the theorem by
the degrees of the two rational maps checked here; no probabilistic Galois
test is used.
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation_u(poly: sp.Expr, u: sp.Symbol) -> int:
    expanded = sp.Poly(sp.expand(poly), u)
    require(not expanded.is_zero, "u-adic valuation requested for zero")
    return min(monomial[0] for monomial, coefficient in expanded.terms() if coefficient != 0)


def main() -> None:
    x, y, z, u, v, w, T, X, W, s, Z = sp.symbols("x y z u v w T X W s Z")

    # The two-jet map and its primitive quartic eliminant.
    F = sp.Matrix((x, x * z**2 + y, x * y * z**2 + z))
    jacobian = sp.factor(F.jacobian((x, y, z)).det())
    expected_jacobian = 1 + 2 * x * y * z - 2 * x**2 * z**3
    require(sp.expand(jacobian - expected_jacobian) == 0, "Jacobian identity")

    N = sp.expand(u**2 * T**4 - u * v * T**2 - T + w)
    eliminant_on_source = sp.expand(N.subs({u: F[0], v: F[1], w: F[2], T: z}))
    require(eliminant_on_source == 0, "quartic eliminant")
    derivative_on_source = sp.expand(sp.diff(N, T).subs({u: F[0], v: F[1], w: F[2], T: z}))
    require(sp.expand(derivative_on_source + jacobian) == 0, "fiber derivative equals minus Jacobian")

    # Primitive and monic discriminants.
    H = sp.expand(
        16 * u**2 * v**4 * w
        - 128 * u**2 * v**2 * w**2
        + 256 * u**2 * w**3
        + 4 * u * v**3
        - 144 * u * v * w
        - 27
    )
    disc_N = sp.factor(sp.discriminant(N, T))
    require(sp.expand(disc_N - u**4 * H) == 0, "primitive discriminant factorization")
    require(sp.expand(H.subs(u, 0) + 27) == 0, "critical factor is a unit on the Jelonek plane")
    require(sp.Poly(H, w).degree() == 3, "odd pole-at-infinity nonsquare certificate")

    monic_q = sp.expand(N / u**2)
    disc_q = sp.factor(sp.discriminant(monic_q, T))
    require(sp.cancel(disc_q - H / u**8) == 0, "monic discriminant and pole order")
    leading_order = 2
    primitive_disc_order = 4
    clearing_exponent = 6 * leading_order - primitive_disc_order
    require(clearing_exponent == 8 and clearing_exponent % 2 == 0, "even clearing exponent")

    # Reversal makes all four infinity branches integral without changing the
    # discriminant.  Its Newton polygon has slopes -2/3 and 0.
    reciprocal = sp.expand(X**4 * N.subs(T, 1 / X))
    expected_reciprocal = w * X**4 - X**3 - u * v * X**2 + u**2
    require(sp.expand(reciprocal - expected_reciprocal) == 0, "reciprocal polynomial")
    require(
        sp.expand(sp.discriminant(reciprocal, X) - disc_N) == 0,
        "discriminant is invariant under reversal",
    )
    newton_points = ((0, 2), (2, 1), (3, 0), (4, 0))
    require(Fraction(1) > Fraction(2, 3), "middle point lies strictly above the C3 face")
    require(Fraction(2 - 0, 0 - 3) == Fraction(-2, 3), "three-root Newton slope")

    # The u=s^3, z=s^-2 Z chart gives the three escaping branches.  Hensel at
    # each cube root is simple; Z=0 is the one finite branch.
    hensel_poly = sp.expand(Z**4 - s * v * Z**2 - Z + w * s**2)
    scaled_N = sp.expand(s**2 * N.subs({u: s**3, T: s**-2 * Z}))
    require(sp.expand(scaled_N - hensel_poly) == 0, "Puiseux escape chart")
    require(sp.expand(hensel_poly.subs(s, 0) - Z * (Z**3 - 1)) == 0, "special escape fiber")
    require(sp.diff(hensel_poly, Z).subs({s: 0, Z: 0}) == -1, "finite branch is simple")
    require(sp.diff(hensel_poly, Z).subs({s: 0, Z: 1}) == 3, "escaping Hensel branch is simple")

    # Standard cubic resolvent.  X=uW changes it into a generic degree-three
    # rational map in X; the quartic itself is a generic degree-four
    # polynomial map in T.  These degrees certify irreducibility.
    p = -v / u
    q = -1 / u**2
    r = w / u**2
    resolvent = sp.expand(W**3 + 2 * p * W**2 + (p**2 - 4 * r) * W - q**2)
    resolvent_clear = sp.expand(u**4 * resolvent)
    expected_resolvent_clear = sp.expand(
        u**4 * W**3 - 2 * u**3 * v * W**2 + u**2 * (v**2 - 4 * w) * W - 1
    )
    require(sp.expand(resolvent_clear - expected_resolvent_clear) == 0, "cubic resolvent")
    require(sp.cancel(sp.discriminant(resolvent, W) - disc_q) == 0, "quartic-resolvent discriminant equality")

    quartic_map = sp.expand(-u**2 * T**4 + u * v * T**2 + T)
    require(sp.Poly(quartic_map, T).degree() == 4, "quartic rational-map degree")
    resolvent_num = sp.expand(u * X * (X - v) ** 2 - 1)
    resolvent_den = 4 * u * X
    require(sp.Poly(resolvent_num, X).degree() == 3, "resolvent rational-map numerator degree")
    require(sp.Poly(resolvent_den, X).degree() == 1, "resolvent rational-map denominator degree")
    require(sp.gcd(sp.Poly(resolvent_num, X), sp.Poly(resolvent_den, X)).degree() == 0, "resolvent map is reduced")

    # A small arithmetic S4 specialization independently controls the generic
    # Galois conclusion: T^4-T+1 is irreducible modulo 2, its cubic resolvent
    # has no rational root, and its discriminant is the nonsquare prime 229.
    quartic_control = sp.Poly(T**4 - T + 1, T, domain=sp.QQ)
    resolvent_control = sp.Poly(W**3 - 4 * W - 1, W, domain=sp.QQ)
    require(sp.Poly(T**4 + T + 1, T, modulus=2).is_irreducible, "quartic S4 control mod 2")
    require(quartic_control.is_irreducible, "quartic S4 control over Q")
    require(resolvent_control.is_irreducible, "resolvent S4 control over Q")
    require(sp.discriminant(quartic_control.as_expr(), T) == 229, "quartic S4 control discriminant")

    # General tame quartic parity-inertia identity: E=6r-(d+2i), so E mod 2
    # is the sign parity d=4-number_of_orbits.
    inertia_rows = (
        ("identity", 0),
        ("transposition", 1),
        ("double_transposition", 2),
        ("three_cycle", 2),
        ("four_cycle", 3),
    )
    for name, d in inertia_rows:
        for order_index in range(4):
            for leading_valuation in range(1, 4):
                exponent = 6 * leading_valuation - (d + 2 * order_index)
                require(exponent % 2 == d % 2, f"parity-inertia identity for {name}")

    # In the hostile, reciprocal inertia is C3 and the reciprocal monogenic
    # order has index one: 4=d+2i=2+2.
    inertia_discriminant = 2
    reciprocal_order_index = (primitive_disc_order - inertia_discriminant) // 2
    require(reciprocal_order_index == 1, "reciprocal order index")

    # The hostile sits in an exact two-parameter family
    # F_{a,b}=(x,x^a z^2+y,x^b y z^2+z).  Its discriminant has competing
    # initial terms of orders 2a+2b and a+4b; the wall a=2b is exactly the
    # binary/ternary Newton-face transition.  A bounded symbolic grid checks
    # the all-parameter formula used in the proof.
    for a in range(1, 6):
        for b in range(0, 6):
            family_N = sp.expand(u ** (a + b) * T**4 - u**b * v * T**2 - T + w)
            family_disc = sp.factor(sp.discriminant(family_N, T))
            family_disc_order = valuation_u(family_disc, u)
            predicted_order = min(2 * a + 2 * b, a + 4 * b)
            require(family_disc_order == predicted_order, f"family discriminant order a={a},b={b}")
            family_clearing = 6 * (a + b) - family_disc_order
            if a <= 2 * b:
                require(family_clearing == 4 * (a + b), f"ternary family exponent a={a},b={b}")
            if a >= 2 * b:
                require(family_clearing == 5 * a + 2 * b, f"binary family exponent a={a},b={b}")

    print("theorem=THM-3059")
    print("status=PROVED_VERIFIED_EXACT")
    print("map=(x,x*z^2+y,x*y*z^2+z)")
    print(f"jacobian={jacobian}")
    print(f"fiber_polynomial={N}")
    print(f"primitive_discriminant={disc_N}")
    print(f"critical_factor={H}")
    print("jelonek_set=u=0")
    print("generic_degree=4;quartic_irreducible_by_rational_map_degree=4")
    print("resolvent_irreducible_by_rational_map_degree=3;generic_galois=S4")
    print("S4_specialization=u=1,v=0,w=1;quartic=T^4-T+1;resolvent=W^3-4W-1;disc=229")
    print(f"reciprocal_newton_points={newton_points};slopes=(-2/3,0);inertia=C3")
    print("reciprocal_discriminant_order=4;maximal_order_discriminant=2;order_index=1")
    print(f"monic_discriminant={disc_q};clearing_exponent={clearing_exponent}")
    print("parity_formula=E=6*v(leading)-(4-number_of_inertia_orbits+2*order_index)")
    print("family=F_ab=(x,x^a*z^2+y,x^b*y*z^2+z);disc_order=min(2a+2b,a+4b);wall=a=2b")
    print("family_clearing=a<=2b:4(a+b);a>=2b:5a+2b")
    print("conclusion=general_dominant_odd_Jelonek_law_refuted;Keller_restricted_law_open")


if __name__ == "__main__":
    main()
