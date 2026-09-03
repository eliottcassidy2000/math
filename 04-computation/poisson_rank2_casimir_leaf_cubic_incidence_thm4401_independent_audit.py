#!/usr/bin/env python3
"""Clean-room exact audit of the THM-4401 Casimir-leaf calculation.

The calculation starts from the three core polynomials in Long's arXiv
source and does not import the THM-4401 primary certificate.  It reconstructs
the nonzero Casimir leaves, their marked-cubic incidence model, the special
R=0 Chinese-remainder/Kummer boundary, and the divisor lost by the Laurent
chart.
"""

from __future__ import annotations

import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def exact(value: sp.Expr) -> sp.Expr:
    """Normalize an exact rational expression."""

    return sp.cancel(sp.together(sp.expand(value)))


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def main() -> None:
    x, y, beta = sp.symbols("x y beta")
    rho, t, u, U = sp.symbols("rho t u U")
    sigma, tau = sp.symbols("S T")

    # Long's three-dimensional core, copied directly from the source formulas.
    xy = x * y
    R = 2 * x - 3 * x**2 * y - x**3 * beta
    S = y + 3 * x * (1 + xy) ** 2 * beta + 3 * x * y**2 * (4 + 3 * xy)
    T = -sp.Rational(1, 2) * (
        (1 + xy) ** 3 * beta + y**2 * (1 + xy) * (4 + 3 * xy)
    )
    w = 1 + xy
    alpha = 2 - 3 * xy - x**2 * beta

    check(exact(R - x * alpha) == 0, "Casimir factorization")
    ambient_jacobian = sp.factor(
        sp.det(sp.Matrix([[sp.diff(f, v) for v in (x, y, beta)] for f in (R, T, S)]))
    )
    check(ambient_jacobian == 1, "ambient core Keller identity")

    # On R=rho != 0, x is a unit.  The following substitutions are the
    # independently reconstructed Laurent coordinates t=1/x and u=t*w.
    leaf_substitution = {
        x: 1 / t,
        y: u - t,
        beta: 5 * t**2 - 3 * t * u - rho * t**3,
    }
    check(exact(R.subs(leaf_substitution) - rho) == 0, "nonzero leaf equation")
    check(exact(w.subs(leaf_substitution) - u / t) == 0, "w coordinate")
    check(exact(alpha.subs(leaf_substitution) - rho * t) == 0, "alpha coordinate")

    # Verify the inverse reconstruction on the localization at x.  Here rho
    # is replaced by the actual polynomial R.
    t_back = 1 / x
    u_back = t_back * w
    x_back = 1 / t_back
    y_back = u_back - t_back
    beta_back = 5 * t_back**2 - 3 * t_back * u_back - R * t_back**3
    check(exact(x_back - x) == 0, "inverse x")
    check(exact(y_back - y) == 0, "inverse y")
    check(exact(beta_back - beta) == 0, "inverse beta")

    S_leaf = 2 * t + 4 * u - 3 * rho * u**2
    T_leaf = (
        rho * u**3 / 2
        - u**2 / 2
        - u * t / 2
    )
    check(exact(S.subs(leaf_substitution) - S_leaf) == 0, "reduced S")
    check(exact(T.subs(leaf_substitution) - T_leaf) == 0, "reduced T")

    leaf_jacobian = sp.factor(
        sp.det(
            sp.Matrix(
                [
                    [sp.diff(S_leaf, t), sp.diff(S_leaf, u)],
                    [sp.diff(T_leaf, t), sp.diff(T_leaf, u)],
                ]
            )
        )
    )
    check(leaf_jacobian == -t, "oriented leaf Jacobian")

    # The source coordinate u is a marked root of this cubic.  Conversely,
    # its derivative reconstructs t, and the incidence equation reconstructs
    # T.  Thus the full (t,u)-plane is exactly the marked-root incidence
    # surface; t != 0 is its simple-marked-root locus.
    P = rho * U**3 - 2 * U**2 + sigma * U + 4 * tau
    P_leaf = sp.expand(P.subs({sigma: S_leaf, tau: T_leaf}))
    check(exact(P_leaf.subs(U, u)) == 0, "marked root incidence")
    P_derivative = sp.diff(P, U)
    check(
        exact(P_derivative.subs({U: u, sigma: S_leaf}) - 2 * t) == 0,
        "marked root derivative",
    )
    t_incidence = (3 * rho * U**2 - 4 * U + sigma) / 2
    S_from_incidence = 2 * t_incidence + 4 * U - 3 * rho * U**2
    T_from_incidence = rho * U**3 / 2 - U**2 / 2 - U * t_incidence / 2
    check(exact(S_from_incidence - sigma) == 0, "incidence inverse S")
    check(exact(T_from_incidence - tau + P / 4) == 0, "incidence inverse T modulo P")
    check(sp.Poly(P, U).LC() == rho, "cubic leading coefficient")

    quotient, remainder = sp.div(
        sp.Poly(P_leaf, U, domain=sp.EX),
        sp.Poly(U - u, U, domain=sp.EX),
    )
    expected_quotient = (
        rho * U**2
        + (rho * u - 2) * U
        + 2 * (t + u - rho * u**2)
    )
    check(remainder.as_expr() == 0, "cubic marked-root division remainder")
    check(exact(quotient.as_expr() - expected_quotient) == 0, "unmarked quadratic")

    cubic_discriminant = sp.factor(sp.discriminant(P, U))
    expected_discriminant = 4 * (
        sigma**2
        - rho * sigma**3
        + 32 * tau
        - 36 * rho * sigma * tau
        - 108 * rho**2 * tau**2
    )
    check(exact(cubic_discriminant - expected_discriminant) == 0, "cubic discriminant")
    unmarked_discriminant = (
        9 * rho**2 * u**2
        - 8 * rho * t
        - 12 * rho * u
        + 4
    )
    check(
        exact(
            sp.discriminant(expected_quotient, U)
            - unmarked_discriminant
        )
        == 0,
        "unmarked-root discriminant",
    )
    pulled_discriminant = exact(
        cubic_discriminant.subs({sigma: S_leaf, tau: T_leaf})
    )
    check(
        exact(pulled_discriminant - 4 * t**2 * unmarked_discriminant) == 0,
        "pulled cubic discriminant",
    )
    check(
        exact(unmarked_discriminant.subs(t, 0) - (3 * rho * u - 2) ** 2) == 0,
        "ramification-line discriminant",
    )

    # The unique triple-root target for rho != 0 is also the cusp point of
    # the discriminant.  Every other geometric cubic has at least one simple
    # root, while this one has none.
    triple_root = sp.Rational(2, 3) / rho
    triple_S = sp.Rational(4, 3) / rho
    triple_T = -sp.Rational(2, 27) / rho**2
    check(
        exact(
            P.subs({sigma: triple_S, tau: triple_T})
            - rho * (U - triple_root) ** 3
        )
        == 0,
        "unique triple-root target",
    )
    check(
        exact(S_leaf.subs({t: 0, u: triple_root}) - triple_S) == 0
        and exact(T_leaf.subs({t: 0, u: triple_root}) - triple_T) == 0,
        "triple root lies on deleted line",
    )

    # The special Casimir fibre R=0 is not obtained by simply restoring
    # t=0: it splits by CRT into x=0 and alpha=0, and the Laurent chart sees
    # only the second component.
    bezout = alpha / 2 + x * (3 * y + x * beta) / 2
    check(exact(bezout - 1) == 0, "CRT Bezout identity")

    S_plane = sp.expand(S.subs(x, 0))
    T_plane = sp.expand(T.subs(x, 0))
    check(S_plane == y, "R=0 plane component S")
    check(T_plane == -beta / 2 - 2 * y**2, "R=0 plane component T")
    plane_inverse = {y: sigma, beta: -2 * tau - 4 * sigma**2}
    check(exact(S_plane.subs(plane_inverse) - sigma) == 0, "plane inverse S")
    check(exact(T_plane.subs(plane_inverse) - tau) == 0, "plane inverse T")
    plane_jacobian = sp.det(
        sp.Matrix(
            [
                [sp.diff(S_plane, y), sp.diff(S_plane, beta)],
                [sp.diff(T_plane, y), sp.diff(T_plane, beta)],
            ]
        )
    )
    check(plane_jacobian == -sp.Rational(1, 2), "plane component automorphism")

    S_zero = sp.expand(S_leaf.subs(rho, 0))
    T_zero = sp.expand(T_leaf.subs(rho, 0))
    zero_jacobian = sp.factor(
        sp.det(
            sp.Matrix(
                [
                    [sp.diff(S_zero, t), sp.diff(S_zero, u)],
                    [sp.diff(T_zero, t), sp.diff(T_zero, u)],
                ]
            )
        )
    )
    check(zero_jacobian == -t, "R=0 Laurent-component Jacobian")
    P_zero = sp.expand(P.subs(rho, 0))
    check(
        exact(P_zero.subs({U: u, sigma: S_zero, tau: T_zero})) == 0,
        "quadratic marked root",
    )
    check(
        exact(sp.diff(P_zero, U).subs({U: u, sigma: S_zero}) - 2 * t) == 0,
        "quadratic root derivative",
    )

    quadratic_discriminant = sigma**2 + 32 * tau
    check(
        exact(sp.discriminant(P_zero, U) - quadratic_discriminant) == 0,
        "standard quadratic discriminant",
    )
    check(
        exact(quadratic_discriminant.subs({sigma: S_zero, tau: T_zero}) - 4 * t**2)
        == 0,
        "pulled quadratic discriminant",
    )
    check(
        exact(expected_discriminant.subs(rho, 0) - 4 * quadratic_discriminant) == 0,
        "cubic-formula specialization scaling",
    )

    # Source and target automorphisms expose the exact punctured Kummer form.
    b = sp.symbols("b")
    u_from_b = (b - 2 * t) / 4
    check(exact(S_zero.subs(u, u_from_b) - b) == 0, "Kummer target b")
    check(
        exact(
            quadratic_discriminant.subs({sigma: S_zero, tau: T_zero})
            - 4 * t**2
        )
        == 0,
        "Kummer square",
    )
    source_change_jacobian = sp.det(
        sp.Matrix([[1, 0], [sp.diff(S_zero, t), sp.diff(S_zero, u)]])
    )
    target_change_jacobian = sp.det(
        sp.Matrix([[1, 0], [2 * sigma, 32]])
    )
    normal_form_jacobian = sp.det(sp.Matrix([[0, 1], [8 * t, 0]]))
    check(source_change_jacobian == 4, "Kummer source automorphism")
    check(target_change_jacobian == 32, "Kummer target automorphism")
    check(normal_form_jacobian == -8 * t, "Kummer normal-form Jacobian")
    check(
        exact(normal_form_jacobian - target_change_jacobian * zero_jacobian / source_change_jacobian)
        == 0,
        "Kummer Jacobian chain rule",
    )

    collision_points = (
        {x: 0, y: 0, beta: -sp.Rational(1, 4)},
        {x: 1, y: -sp.Rational(3, 2), beta: sp.Rational(13, 2)},
        {x: -1, y: sp.Rational(3, 2), beta: sp.Rational(13, 2)},
    )
    for index, point in enumerate(collision_points):
        check(exact(R.subs(point)) == 0, f"collision point {index} R")
        check(exact(S.subs(point)) == 0, f"collision point {index} S")
        check(exact(T.subs(point) - sp.Rational(1, 8)) == 0, f"collision point {index} T")

    print("THM-4401 INDEPENDENT CASIMIR-LEAF AUDIT")
    print("source: Long arXiv:2608.23777 core formulas; no THM-4401 primary import")
    print("base field: characteristic zero; rho is inverted in the nonzero-leaf block")
    print()
    print("NONZERO CASIMIR LEAF R=rho")
    print("coordinate ring = K[t,t^-1,u], with t=1/x and u=t(1+xy)")
    print("x=1/t, y=u-t, beta=5*t^2-3*t*u-rho*t^3")
    print(f"S = {S_leaf}")
    print(f"T = {T_leaf}")
    print("Jac_(t,u)(S,T) = -t (reversing S,T changes the sign)")
    print("P_rho(U) = rho*U^3-2*U^2+S*U+4*T")
    print("P_rho(u)=0 and P_rho'(u)=2*t")
    print("the Laurent leaf is exactly the simple-marked-root locus t!=0")
    print("for rho!=0 the full incidence completion is finite flat of degree 3")
    print(f"Disc_U(P_rho) = {expected_discriminant}")
    print(f"pullback Disc = 4*t^2*({unmarked_discriminant})")
    print("the second factor is the discriminant of the two unmarked roots")
    print(
        "unique triple-root target: "
        f"(S,T)=({triple_S},{triple_T}), root={triple_root}"
    )
    print("geometric image of the punctured leaf = A^2 minus that triple-root point")
    print()
    print("R=0 BOUNDARY")
    print("R=x*alpha and 1=alpha/2+x(3*y+x*beta)/2, so CRT gives two components")
    print("x=0 component: (S,T)=(y,-beta/2-2*y^2), a polynomial automorphism")
    print("alpha=0 component: K[t,t^-1,u] with the rho=0 formulas above")
    print("P_0(U)=-2*U^2+S*U+4*T")
    print("standard quadratic discriminant Delta_2=S^2+32*T, pullback Delta_2=4*t^2")
    print("the rho=0 specialization of the cubic discriminant is 4*Delta_2=16*t^2")
    print("after b=S, Delta_2=S^2+32*T, the Laurent map is (t,b)->(b,4*t^2)")
    print("target (S,T)=(0,1/8) has one plane point and the Kummer pair t=+/-1")
    print()
    print("PUNCTURE / COMPLETION FIREWALL")
    print("adjoining t=0 gives the unique polynomial extension on this (t,u)-plane")
    print("its Jacobian is -t, so the added line is ramified, not Keller")
    print("on that line P_rho(u)=P_rho'(u)=0: it is the marked-root ramification divisor")
    print("unmarked-root collisions may occur with t!=0 and do not spoil etaleness")
    print("scope: this does not exclude another affine modification or another cubic filling")
    print()
    print(f"CHECKS {CHECKS}")
    print("PASS")


if __name__ == "__main__":
    main()
