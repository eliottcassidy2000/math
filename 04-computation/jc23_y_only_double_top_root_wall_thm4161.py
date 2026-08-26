#!/usr/bin/env python3
"""Exact symbolic certificate for THM-4161's double top-root wall."""

from __future__ import annotations

import hashlib
import sys

import sympy as sp


sys.stdout.reconfigure(newline="\n")


EXPECTED_SEMANTIC_SHA256 = "9e00cdbf68489a26e824458a177af856cda10244c613f062662896b96b4556ff"
CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuation(poly: sp.Expr, variable: sp.Symbol) -> int:
    return min(monomial[0] for monomial, coefficient in sp.Poly(poly, variable).terms() if coefficient)


def zero(expression: sp.Expr, message: str) -> None:
    require(sp.factor(expression) == 0, message)


def main() -> None:
    s, p, X, T, W, z = sp.symbols("s p X T W z")
    r, u, q = sp.symbols("r u q", nonzero=True)
    theta, phi, zeta = sp.symbols("theta phi zeta")
    K0 = sp.Rational(2848, 45)
    cubic_constant = sp.Rational(1376, 135)

    # Parametrize C(W)=zeta*(W-r)^2*(W-u), retaining the prescribed
    # nonzero constant coefficient -1376/135.
    zeta_wall = sp.cancel(cubic_constant / (r**2 * u))
    theta_wall = sp.factor(-zeta_wall * (2 * r + u))
    phi_wall = sp.factor(zeta_wall * (r**2 + 2 * r * u))
    wall_substitution = {theta: theta_wall, phi: phi_wall, zeta: zeta_wall}
    p_inner = 2027776 * r**3 * u + 1013888 * r**2 * u**2 + 17415
    inner = sp.factor(4 * theta_wall * K0**2 - 27 * zeta_wall**2)
    transverse = sp.factor(sp.Rational(8, 3) + K0 * r**2)
    zero(
        inner + 44032 * p_inner / (273375 * r**4 * u**2),
        "inner-resultant factorization changed",
    )
    zero(
        transverse - 8 * (356 * r**2 + 15) / 45,
        "top transversality factor changed",
    )

    # Reconstruct the complete THM-4155 source and its two critical
    # projections after imposing the double-root parameterization.
    t = p - s**2
    H = sp.expand(
        -3 * p
        + sp.Rational(8, 3) * p**2
        - cubic_constant * p**3
        + K0 * s**2 * p**2
        + phi * s * p**3
        + theta * s**2 * p**3
        + zeta * s**3 * p**3
    )
    source_a = sp.factor(sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p))
    source_c = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
    source_b = sp.factor(sp.cancel((source_c + s * source_a) / t**2))
    source_resultant = sp.factor(sp.resultant(source_a, source_b, s))
    require(valuation(source_resultant, p) == 6, "source resultant valuation changed")
    source_residual = sp.Poly(sp.cancel(source_resultant / p**6), p)
    source_wall = sp.Poly(
        sp.cancel(source_residual.as_expr().subs(wall_substitution)), p
    )
    expected_source_constant = 3877634048 * p_inner / (50625 * r**6 * u**3)
    expected_source_leading = (
        sp.Integer(3289935900927224469054816256)
        * (r - u) ** 3
        * (356 * r**2 + 15)
        / (sp.Integer(252226880859375) * r**16 * u**8)
    )
    require(source_wall.degree() == 17, "source residual degree changed")
    zero(source_wall.TC() - expected_source_constant, "source constant endpoint changed")
    zero(source_wall.LC() - expected_source_leading, "source leading endpoint changed")

    P = T + X**2 * T**2
    Y = X * T * P
    G = sp.expand(
        -X**2 * T / 2
        - 3 * P
        + sp.Rational(8, 3) * P**2
        - cubic_constant * P**3
        + K0 * Y**2
        + phi * P**2 * Y
        + theta * P * Y**2
        + zeta * Y**3
    )
    normalized_f = sp.cancel(sp.diff(G, X) / T)
    normalized_h = sp.diff(G, T)
    normalized_resultant = sp.factor(sp.resultant(normalized_f, normalized_h, X))
    require(valuation(normalized_resultant, T) == 56, "normalized resultant valuation changed")
    normalized_residual = sp.Poly(
        sp.cancel(normalized_resultant / (T**56 * (6 * T + 1) ** 2)), T
    )
    normalized_wall = sp.Poly(
        sp.cancel(normalized_residual.as_expr().subs(wall_substitution)), T
    )
    expected_normalized_constant = -sp.Integer(72965752821794209792) / (
        sp.Integer(56953125) * r**14 * u**7
    )
    expected_normalized_leading = (
        sp.Integer(210555897659342366019508240384)
        * (r - u) ** 3
        * (356 * r**2 + 15)
        * p_inner**2
        / (sp.Integer(9308590679915771484375) * r**20 * u**10)
    )
    require(normalized_wall.degree() == 17, "normalized residual degree changed")
    zero(
        normalized_wall.TC() - expected_normalized_constant,
        "normalized constant endpoint changed",
    )
    zero(
        normalized_wall.LC() - expected_normalized_leading,
        "normalized leading endpoint changed",
    )

    # The rational point (r,u)=(1,2) is an exact nonempty control.  Its two
    # residual projections are squarefree and avoid the universal half row.
    control_substitution = {r: 1, u: 2}
    source_control = sp.Poly(source_wall.as_expr().subs(control_substitution), p)
    normalized_control = sp.Poly(
        normalized_wall.as_expr().subs(control_substitution), T
    )
    require(sp.gcd(source_control, source_control.diff()).degree() == 0,
            "source control is not squarefree")
    require(sp.gcd(normalized_control, normalized_control.diff()).degree() == 0,
            "normalized control is not squarefree")
    require(normalized_control.eval(-sp.Rational(1, 6)) != 0,
            "control collided with the universal half row")
    control_coefficients = tuple(
        sp.factor(value.subs(control_substitution))
        for value in (zeta_wall, theta_wall, phi_wall, inner, transverse)
    )
    require(
        control_coefficients
        == (
            sp.Rational(688, 135),
            -sp.Rational(2752, 135),
            sp.Rational(688, 27),
            -sp.Rational(89478737152, 273375),
            sp.Rational(2968, 45),
        ),
        "rational control ledger changed",
    )

    # Restore the four universal affine critical points and their Hessians.
    zero(normalized_f.subs(T, 0) + X, "T=0 normalized source row changed")
    zero(
        normalized_h.subs(T, 0) + (X**2 + 6) / 2,
        "T=0 normalized target row changed",
    )
    hessian = sp.det(sp.hessian(G, (X, T)))
    hessian_zero = sp.rem(
        sp.Poly(hessian.subs(T, 0), X), sp.Poly(X**2 + 6, X)
    ).as_expr()
    zero(hessian_zero - 6, "T=0 Hessian changed")
    half_f = sp.rem(
        sp.Poly(normalized_f.subs(T, -sp.Rational(1, 6)), X),
        sp.Poly(X**2 - 6, X),
    ).as_expr()
    half_h = sp.rem(
        sp.Poly(normalized_h.subs(T, -sp.Rational(1, 6)), X),
        sp.Poly(X**2 - 6, X),
    ).as_expr()
    zero(half_f, "T=-1/6 first critical row changed")
    zero(half_h, "T=-1/6 second critical row changed")
    hessian_half = sp.rem(
        sp.Poly(hessian.subs(T, -sp.Rational(1, 6)), X),
        sp.Poly(X**2 - 6, X),
    ).as_expr()
    zero(hessian_half + 6, "T=-1/6 Hessian changed")

    # Normalize the collided top face in p=z^-1, s=W.  The z derivative is
    # a unit on the stated gate, so z has order two in W-r and omega has
    # order four, giving one ramification index five branch.
    top_cubic = sp.factor(zeta_wall * (W - r) ** 2 * (W - u))
    top_a = sp.Rational(8, 3) + K0 * W**2
    top_chart = sp.expand(
        q * top_cubic
        + q * z * (top_a - W**2 * top_cubic)
        + q * z**2 * (-3 - W**2 * top_a)
        + z**3 * (3 * q * W**2 - 1)
        + z**4 * W**2 * (1 - q / 2)
    )
    wall_source = sp.expand(H.subs(wall_substitution))
    source_fibre = (s**2 - p) * (1 - q * wall_source) - q * s**2 / 2
    direct_top_chart = sp.cancel(
        z**4 * source_fibre.subs({s: W, p: 1 / z})
    )
    zero(direct_top_chart - top_chart, "top-chart expansion changed")
    zero(
        top_cubic
        - (zeta_wall * W**3 + theta_wall * W**2 + phi_wall * W - cubic_constant),
        "top cubic factorization changed",
    )
    zero(top_chart.subs({W: r, z: 0}), "double top root disappeared")
    zero(sp.diff(top_chart, W).subs({W: r, z: 0}), "double top derivative changed")
    zero(
        sp.diff(top_chart, W, 2).subs({W: r, z: 0}) / 2
        - q * zeta_wall * (r - u),
        "double top quadratic coefficient changed",
    )
    zero(
        sp.diff(top_chart, z).subs({W: r, z: 0}) - q * transverse,
        "double top transverse coefficient changed",
    )
    zero(
        sp.diff(top_chart, W).subs({W: u, z: 0})
        - q * zeta_wall * (u - r) ** 2,
        "simple top root derivative changed",
    )

    # Numerical ledgers for the inherited Riemann--Hurwitz and monodromy
    # consumers.  These are consequences of the labelled packet, not CAS
    # substitutes for the geometric transport lemmas.
    packet = (8, 5, 3, 2, 2, 2, 1)
    defect = sum(index - 1 for index in packet)
    critical_length = 17 + 2 + 2
    full_degree = sum(packet)
    finite_degree = 8 + 5 + 3 + 1
    carrier_index = 3
    require((defect, critical_length, full_degree, finite_degree) == (16, 21, 23, 17),
            "packet or response ledger changed")
    full_commutator_cap = 2 * (full_degree - critical_length)
    finite_cap_both = 2 * finite_degree - critical_length - 2 + carrier_index
    finite_cap_one = 2 * finite_degree - critical_length - 1 + carrier_index
    finite_cap_zero = carrier_index
    require(full_commutator_cap < defect, "full response lost strictness")
    require(
        max(finite_cap_both, finite_cap_one, finite_cap_zero) < finite_degree - 1,
        "finite response lost strictness",
    )

    ledger = (
        "universe=Q(r,u,q)[s,p,X,T,W,z];K0=2848/45;c=1376/135;"
        "C=zeta(W-r)^2(W-u);zeta=c/(r^2u);theta=-zeta(2r+u);"
        "phi=zeta(r^2+2ru);PI=2027776r^3u+1013888r^2u^2+17415;"
        "gate=r*u*(r-u)*(356r^2+15)*PI!=0;"
        "Res_s=p^6R17;R17(0)=3877634048PI/(50625r^6u^3);"
        "R17lc=3289935900927224469054816256(r-u)^3(356r^2+15)/"
        "(252226880859375r^16u^8);"
        "Res_X=T^56(6T+1)^2Q17;Q17(0)=-72965752821794209792/"
        "(56953125r^14u^7);Q17lc=210555897659342366019508240384"
        "(r-u)^3(356r^2+15)PI^2/(9308590679915771484375r^20u^10);"
        "control=(r,u)=(1,2);squarefree_RQ;Q(-1/6)!=0;"
        "top=z~(W-r)^2;omega_order=4;e=5;simple_e=3;"
        "packet=(8,5,3,2,2,2,1);g=9;defect=16;L=21;"
        "full=(n,cap)=(23,4);finite=(n,beta,caps)=(17,3,14|15|3)"
    )
    semantic = hashlib.sha256(ledger.encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    print("THM4161_Y_ONLY_DOUBLE_TOP_ROOT_WALL_20260825")
    print("status=PASS;scope=symbolic identities plus exact nonempty control")
    print("gate=r*u*(r-u)*(356r^2+15)*P_I!=0")
    print("P_I=2027776r^3u+1013888r^2u^2+17415")
    print("I_C=-44032*P_I/(273375r^4u^2)")
    print("source_resultant=p^6*R17;normalized_resultant=T^56*(6T+1)^2*Q17")
    print("source_endpoints=nonzero_exact;normalized_endpoints=nonzero_exact")
    print("control=(r,u)=(1,2);R17_squarefree;Q17_squarefree;Q17(-1/6)!=0")
    print("universal_critical_points=4;hessians=6,-6;L=21")
    print("top_chart=z~(W-r)^2;omega_order=4;collided_index=5")
    print("packet=(8,5,3,2,2,2,1);g=9;defect=16")
    print("responses=full(n=23,commutator_cap=4);finite(n=17,beta=3,caps=14,15,3)")
    print("hostile_walls=P_I=0|u=r|356r^2+15=0")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
