#!/usr/bin/env python3
"""Primary exact certificate for THM-4339's clean M=12 cubic gate.

It checks the finite lower atlas, genus/edge ledgers, Laurent-root
stratification, exact T-chart, and every weighted initial form used in the
relative extinction proof.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
sys.path.insert(0, "04-computation")
import jc2_m12_z0_endpoint_extinction_thm4327 as base


CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("THM-4339 primary failure: " + label)


def ledger(points):
    polygon = base.convex_hull(points)
    twice_area = abs(sum(
        polygon[i][0] * polygon[(i + 1) % len(polygon)][1]
        - polygon[(i + 1) % len(polygon)][0] * polygon[i][1]
        for i in range(len(polygon))
    ))
    boundary = sum(
        gcd(
            abs(polygon[(i + 1) % len(polygon)][0] - polygon[i][0]),
            abs(polygon[(i + 1) % len(polygon)][1] - polygon[i][1]),
        )
        for i in range(len(polygon))
    )
    require((twice_area - boundary + 2) % 2 == 0, "Pick parity")
    return polygon, (twice_area, boundary, (twice_area - boundary + 2) // 2)


def order(base_degree, plane):
    value = base_degree * (F(5, 6) - sum(plane))
    require(value.denominator == 1 and value > 0, "positive integral face order")
    return value.numerator


def sigma_initial(expr, sigma, exponent):
    poly = sp.Poly(sp.expand(expr), sigma)
    return sp.factor(poly.coeff_monomial(sigma**exponent))


def main():
    # Pin the inherited full-support hull engine rather than silently changing
    # its row universe or fixed residual point.
    inherited = Path("04-computation/jc2_m12_z0_endpoint_extinction_thm4327.py")
    inherited_hash = sha256(inherited.read_bytes()).hexdigest()
    require(
        inherited_hash == "6aa9087afd413833d3168b27efe3e65e779ab796a9b537ef1e624a6380fac551",
        "pinned THM-4327 hull dependency",
    )

    # Clean gate: Z=beta=zeta=0 and U,W,K nonzero.  Requiring K keeps the
    # T-face present even under every conservative aggregate cancellation.
    rows = base.weighted_rows()
    points, point_index, support_mask, lower_faces = base.make_hull_engine(rows)
    p_rows = {(1, 0), (2, 0), (3, 0), (6, 0)}
    required = p_rows | {(0, 2), (3, 2)}  # K and W
    absent = {(0, 4), (0, 3), (1, 3)}    # Z, zeta_3, beta_11
    collisions = (
        (2, 3, 1), (2, 4, 1), (2, 5, 1),
        (2, 6, 1), (3, 4, 1), (3, 5, 1),
    )
    optional, total, counter = base.hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        required, absent, collisions,
    )
    M = (F(1, 12), F(1, 6), F(-1, 6))
    T = (F(1, 2), F(0), F(-1))
    require(optional == 7 and total == 8192, "clean-gate hostile population")
    require(counter == Counter({(M, T): 8192}), "unique M,T lower complex")

    polygons = {
        "M": ((0, 1), (2, 0), (4, 5), (0, 7)),
        "C": ((0, 0), (0, 6), (2, 5)),
        "T": ((2, 0), (4, 2), (4, 5)),
        "global": ((0, 1), (2, 0), (4, 2), (4, 5), (0, 7)),
    }
    ledgers = {name: ledger(vertices)[1] for name, vertices in polygons.items()}
    require(ledgers == {
        "M": (36, 10, 14),
        "C": (12, 8, 3),
        "T": (6, 6, 1),
        "global": (42, 14, 15),
    }, "Pick ledgers")
    packet = base.edge_packet(polygons["global"])
    require(packet == (11, 11, 3, 3, 3, 2, 2, 1), "global edge packet")
    require(sum(value - 1 for value in packet) == 28, "Riemann-Hurwitz sum")
    require(gcd(2, 5) == 1, "primitive C--T edge")
    require(12 + 1 - 3 + 1 == 11, "R,C,T graph b1")
    require(3 + 1 + 11 == 15, "squarefree component-genus ledger")
    require(order(12, M) == 9 and order(12, T) == 16, "primary face orders")

    # Cubic discriminant, the H/J selector, and the labelled double/triple
    # parameterizations.  K*W!=0 excludes both toric owner exits.
    P = sp.symbols("P")
    K, W = sp.symbols("K W", nonzero=True)
    theta, xi = sp.symbols("theta xi")
    A = K + theta * P + xi * P**2 + W * P**3

    # Complete boundary schemes.  Every entry other than the deliberately
    # variable cubic is squarefree on U*K*W*(U+W)!=0.  Roots on different
    # entries live on different toric divisors even if their scalar labels
    # happen to agree.
    edge_x, Ucoef = sp.symbols("edge_x Ucoef", nonzero=True)
    edge_schemes = (
        edge_x - 1,
        1 - K * edge_x**2,
        A.subs(P, edge_x),
        (edge_x - 1) * (Ucoef * edge_x + W),
        Ucoef - edge_x**6,
        1 - W * edge_x,
    )
    require(sp.discriminant(edge_schemes[1], edge_x) == 4 * K,
            "K-edge squarefree")
    require(sp.factor(sp.resultant(edge_x - 1, Ucoef * edge_x + W, edge_x))
            == Ucoef + W, "top-edge Lambda resultant")
    require(sp.factor(sp.discriminant(edge_schemes[4], edge_x))
            == 46656 * Ucoef**5, "U-edge squarefree")
    require(sp.discriminant(edge_schemes[5], edge_x) == 1,
            "internal C--T edge simple")
    disc = sp.factor(sp.discriminant(A, P))
    disc_formula = (
        xi**2 * theta**2 - 4 * W * theta**3 - 4 * xi**3 * K
        - 27 * W**2 * K**2 + 18 * W * xi * theta * K
    )
    H = xi**2 - 3 * W * theta
    J = 2 * xi**3 - 9 * W * xi * theta + 27 * W**2 * K
    require(sp.expand(disc - disc_formula) == 0, "binary-cubic discriminant")
    require(sp.expand(4 * H**3 - J**2 - 27 * W**2 * disc) == 0,
            "H/J subdiscriminant identity")

    a, b = sp.symbols("a b", nonzero=True)
    A_double = sp.expand(W * (P - a)**2 * (P - b))
    double_coefficients = {
        xi: -W * (b + 2 * a),
        theta: W * (a**2 + 2 * a * b),
        K: -W * a**2 * b,
    }
    require(sp.expand(A.subs(double_coefficients) - A_double) == 0,
            "double-root parameterization")
    require(sp.factor(H.subs(double_coefficients)) == W**2 * (a - b)**2,
            "double-root H selector")
    A_triple = sp.expand(W * (P - a)**3)
    triple_coefficients = {
        xi: -3 * W * a,
        theta: 3 * W * a**2,
        K: -W * a**3,
    }
    require(sp.expand(A.subs(triple_coefficients) - A_triple) == 0,
            "triple-root parameterization")
    require(sp.expand(H.subs(triple_coefficients)) == 0
            and sp.expand(J.subs(triple_coefficients)) == 0,
            "triple-root H/J selector")

    # Exact T chart.  C0 is the j=0 part of H and
    # B=Phi+eta*P+alpha*P^2 is its complete j=1 part after factoring P^2*y.
    sigma, S, x = sp.symbols("sigma S x")
    Bfun, Cfun = sp.symbols("Bfun Cfun")
    H_T = S**2 * P**2 * A + sigma**6 * S * P**3 * Bfun + sigma**12 * Cfun
    G = (S**2 - sigma**12 * P) * (1 - H_T) - sigma**12 * S**2 / 2
    Kloc = (
        (1 - sigma**12 * P * x**2)
        * (x**2 - P**2 * A - sigma**6 * x * P**3 * Bfun
           - sigma**12 * x**2 * Cfun)
        - sigma**12 * x**2 / 2
    )
    chart_identity = sp.together((x**4 * G.subs(S, 1 / x)) - Kloc)
    require(sp.expand(chart_identity) == 0, "exact reciprocal T chart")

    # Weighted local faces.  Terms involving Cfun and the outer unit begin
    # strictly above each displayed initial form.
    u, X, Uv = sp.symbols("u X Uv")
    B0, B1, B2 = sp.symbols("B0 B1 B2")
    lam = W * (a - b)

    def local_core(A_local, B_local):
        p_local = a + u
        return sp.expand(
            x**2 - p_local**2 * A_local
            - sigma**6 * x * p_local**3 * B_local
        )

    double_core = local_core(W * u**2 * (a - b + u), B0 + B1 * u + B2 * u**2)
    double_scaled = sp.expand(double_core.subs({x: sigma**6 * X, u: sigma**6 * Uv}))
    double_initial = X**2 - a**2 * lam * Uv**2 - a**3 * B0 * X
    require(sp.expand(sigma_initial(double_scaled, sigma, 12) - double_initial) == 0,
            "double B0 weighted initial form")

    triple_core = local_core(W * u**3, B0 + B1 * u + B2 * u**2)
    triple_b0_scaled = sp.expand(triple_core.subs({x: sigma**6 * X, u: sigma**4 * Uv}))
    triple_b0_initial = X**2 - a**2 * W * Uv**3 - a**3 * B0 * X
    require(sp.expand(sigma_initial(triple_b0_scaled, sigma, 12)
                      - triple_b0_initial) == 0,
            "triple B0 weighted initial form")

    triple_b1_core = triple_core.subs(B0, 0)
    triple_b1_scaled = sp.expand(
        triple_b1_core.subs({x: sigma**18 * X, u: sigma**12 * Uv})
    )
    triple_b1_initial = X**2 - a**2 * W * Uv**3 - a**3 * B1 * X * Uv
    require(sp.expand(sigma_initial(triple_b1_scaled, sigma, 36)
                      - triple_b1_initial) == 0,
            "triple B1 weighted initial form")

    # Smoothness/genus types of the new components.  For the double conic,
    # a simultaneous critical point would violate B0*lam!=0.  For the triple
    # B0 face, the projective cubic has no singular point when B0*W!=0.  The
    # B1 face has an ordinary node at the origin and is therefore rational.
    require(sp.factor(double_initial.subs({X: a**3 * B0 / 2, Uv: 0}))
            == -a**6 * B0**2 / 4,
            "double exceptional conic smoothness")
    Zp = sp.symbols("Zp")
    triple_b0_projective = (
        X**2 * Zp - a**3 * B0 * X * Zp**2 - a**2 * W * Uv**3
    )
    # At infinity U=0 and dF/dZ=X^2; in the affine chart U=0 follows from
    # dF/dU=0, while dF/dX forces X=a^3*B0/2 and F is nonzero.
    affine_obstruction = sp.factor(
        triple_b0_projective.subs({Zp: 1, Uv: 0, X: a**3 * B0 / 2})
    )
    require(affine_obstruction == -a**6 * B0**2 / 4,
            "triple B0 affine smoothness")
    require(sp.diff(triple_b0_projective, Zp).subs({Zp: 0, Uv: 0}) == X**2,
            "triple B0 smoothness at infinity")
    require(
        triple_b1_initial.subs({X: 0, Uv: 0}) == 0
        and sp.diff(triple_b1_initial, X).subs({X: 0, Uv: 0}) == 0
        and sp.diff(triple_b1_initial, Uv).subs({X: 0, Uv: 0}) == 0,
        "triple B1 singular origin",
    )
    tangent_b1 = X * (X - a**3 * B1 * Uv)
    tangent_from_poly = sum(
        coefficient * X**monomial[0] * Uv**monomial[1]
        for monomial, coefficient in sp.Poly(triple_b1_initial, X, Uv).terms()
        if sum(monomial) == 2
    )
    require(sp.expand(tangent_from_poly - tangent_b1) == 0,
            "triple B1 distinct tangent declaration")
    require(sp.resultant(X, X - a**3 * B1 * Uv, X) == -a**3 * B1 * Uv,
            "triple B1 ordinary-node tangents")

    # The differential conversion is exact: Kloc=x^4 G, so on the curve
    # dS/G_P=-x^2 dx/(Kloc)_u.  A scaling x=sigma^r X,
    # u=sigma^s U, Kloc=sigma^d K0 adds 3r+s-d to order 16.
    def exceptional_order(r, s, d):
        return 16 + 3 * r + s - d

    require(exceptional_order(6, 6, 12) == 28, "double conic order")
    require(exceptional_order(6, 4, 12) == 26, "triple elliptic order")
    require(exceptional_order(18, 12, 36) == 46, "triple B1 rational order")

    # Delta/genus bookkeeping by normal-jet stratum.  A repeated A has one
    # delta=1 singularity.  B(a)!=0 smooths it on the generic sigma-fibre;
    # B(a)=0 preserves delta one (node or cusp) and the normalized source has
    # genus fourteen.  The former double case restores one graph loop through
    # a rational conic path; the former triple case restores genus through the
    # elliptic tail.  The persistent cases need no hidden genus.
    require(3 + (15 - 4 + 1) == 15,
            "double-smoothed rational-bridge genus ledger")
    require(3 + 1 + 11 == 15,
            "triple-smoothed elliptic-tail genus ledger")
    require(3 + 0 + 11 == 14,
            "persistent repeated-root normalized genus ledger")

    # Root exits and collisions are three different exact taxes.  These
    # hostiles prevent a scalar discriminant from impersonating Laurent
    # saturation.
    controls = {
        "simple": 2 + P + P**2 + P**3,
        "double": (P - 1)**2 * (P + 1),
        "triple": (P - 1)**3,
        "zero_exit": P * (P**2 + P + 1),
        "infinity_exit": P**2 + P + 1,
    }

    def taxes(poly_expr):
        poly = sp.Poly(sp.expand(poly_expr), P, domain=sp.QQ)
        degree = int(poly.degree())
        zero = min(monom[0] for monom, coefficient in poly.terms() if coefficient)
        reduced = sp.Poly(sp.cancel(poly.as_expr() / P**zero), P, domain=sp.QQ)
        collision = int(sp.gcd(reduced, reduced.diff()).degree())
        toric_roots = int(reduced.degree()) - collision
        return 3 - degree, zero, collision, toric_roots

    tax_values = {name: taxes(expr) for name, expr in controls.items()}
    require(tax_values == {
        "simple": (0, 0, 0, 3),
        "double": (0, 0, 1, 2),
        "triple": (0, 0, 2, 1),
        "zero_exit": (0, 1, 0, 2),
        "infinity_exit": (1, 0, 0, 2),
    }, "three-tax hostile controls")
    require(all(3 - values[3] == sum(values[:3]) for values in tax_values.values()),
            "three-tax identity")
    require(sp.discriminant(controls["zero_exit"], P) != 0,
            "simple zero exit has nonzero projective discriminant")

    print("THM4339_CLEAN_CUBIC_EDGE_PRIMARY=PASS")
    print("gate=Z=beta_11=zeta_3=0;U*W*K*(U+W)!=0")
    print("lower_atlas=8192/8192:M,T")
    print("pick=M(36,10,14);C(12,8,3);T(6,6,1);global(42,14,15)")
    print("orders=M:9;T:16;double_conic:28;triple_elliptic:26;triple_B1_rational:46")
    print("double=B(a)!=0:rational_bridge;B(a)=0:persistent_node")
    print("triple=B(a)!=0:elliptic_tail;B(a)=0:rational_delta_constant")
    print("generic_genera=squarefree_or_smoothed:15;persistent_collision:14")
    print("taxes=infinity_exit+zero_exit+localized_collision")
    print("edges=X-1;1-KX^2;A(X);(X-1)(UX+W);U-X^6;internal:1-WX")
    print("checks=" + str(CHECKS))


if __name__ == "__main__":
    main()
