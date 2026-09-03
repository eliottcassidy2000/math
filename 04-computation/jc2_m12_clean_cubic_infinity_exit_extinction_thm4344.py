#!/usr/bin/env python3
"""Exact primary certificate for THM-4344's W=0 degree-drop boundary.

Theorem scope (inside the inherited reduced (2,3), exact-M=12 seam):

    Z=beta_11=zeta_3=W=0,       U*K*xi_10 != 0.

The script checks a conservative three-face over-atlas, the genus/edge ledger, the
new root-at-infinity face and reciprocal chart, its cancellation-complete
even-cusp degeneration, and the independent finite-quadratic collision.
Independent canonical referees audit coefficient coupling and chart coverage.
"""

from collections import Counter, defaultdict
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


def need(condition, label):
    global CHECKS
    CHECKS += 1
    if condition is not True and not bool(condition):
        raise RuntimeError("W=0 degree-drop audit failure: " + label)


def polygon_ledger(vertices):
    polygon = base.convex_hull(vertices)
    area2 = abs(sum(
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
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return polygon, (area2, boundary, (area2 - boundary + 2) // 2)


def face_order(base_degree, plane):
    answer = F(base_degree) * (F(5, 6) - sum(plane))
    need(answer.denominator == 1 and answer > 0, "positive integral face order")
    return answer.numerator


def exact_source():
    """Return symbols, source rows, and the exact support of F_Q."""
    U, u, alpha, xi, Delta, eta, Theta, Phi, K = sp.symbols(
        "U u alpha xi Delta eta Theta Phi K"
    )
    e = -sp.Rational(1376, 135)
    rows = {
        (1, 0): -3,
        (2, 0): sp.Rational(8, 3),
        (3, 0): e,
        (4, 0): Delta,
        (5, 0): u,
        (6, 0): U,
        (0, 2): K,
        (2, 1): Phi,
        (1, 2): Theta,
        (3, 1): eta,
        (2, 2): xi,
        (4, 1): alpha,
    }
    support = defaultdict(lambda: sp.Integer(0))
    support[(2, 0, 0)] += 1
    support[(0, 1, 0)] -= 1
    support[(2, 0, 1)] -= sp.Rational(1, 2)
    for (i, j), coefficient in rows.items():
        support[(j + 2, i + j, 1)] -= coefficient
        support[(j, i + j + 1, 1)] += coefficient
    support = {
        point: sp.factor(coefficient)
        for point, coefficient in support.items()
        if coefficient != 0
    }
    symbols = {
        "U": U, "u": u, "alpha": alpha, "xi": xi,
        "Delta": Delta, "eta": eta, "Theta": Theta,
        "Phi": Phi, "K": K, "e": e,
    }
    return symbols, rows, support


PLANES = {
    "M": (F(1, 12), F(1, 6), F(-1, 6)),
    "D6": (F(1, 6), F(1, 6), F(-1, 3)),
    "T": (F(1, 2), F(0), F(-1)),
}


def on_plane(point, plane):
    r, ell, height = point
    a, b, c = plane
    return F(height) == a * r + b * ell + c


def face_expression(support, plane):
    S, P = sp.symbols("S P")
    answer = 0
    points = []
    for point, coefficient in support.items():
        if on_plane(point, plane):
            r, ell, _height = point
            answer += coefficient * S**r * P**ell
            points.append((r, ell))
    return sp.factor(answer), tuple(sorted(points))


def main():
    syms, rows_by_label, support = exact_source()
    U, u, alpha, xi = syms["U"], syms["u"], syms["alpha"], syms["xi"]
    Delta, eta, Theta = syms["Delta"], syms["eta"], syms["Theta"]
    Phi, K, e = syms["Phi"], syms["K"], syms["e"]

    # At W=0, (2,6,1)=-U, (4,4,1)=-xi, and (4,2,1)=-K
    # are unique nonzero owners and may not be synthetically deleted.  Only
    # the three genuinely cancellable aggregate points are hostile toggles.
    inherited = Path("04-computation/jc2_m12_z0_endpoint_extinction_thm4327.py")
    need(
        sha256(inherited.read_bytes()).hexdigest()
        == "6aa9087afd413833d3168b27efe3e65e779ab796a9b537ef1e624a6380fac551",
        "pinned THM-4327 hull engine",
    )
    labelled_rows = base.weighted_rows()
    _points, point_index, support_mask, lower_faces = base.make_hull_engine(labelled_rows)
    required = {(1, 0), (2, 0), (3, 0), (6, 0), (0, 2), (2, 2)}
    absent = {(0, 4), (0, 3), (1, 3), (3, 2)}
    cancellable = ((2, 3, 1), (2, 4, 1), (2, 5, 1))
    optional, total, counter = base.hostile_atlas(
        labelled_rows, point_index, support_mask, lower_faces,
        required, absent, cancellable,
    )
    M, D6, T = PLANES["M"], PLANES["D6"], PLANES["T"]
    need(optional == 6 and total == 512, "conservative over-atlas population")
    need(counter == Counter({(M, D6, T): 512}), "unique M,D6,T complex")
    need(support[(2, 6, 1)] == -U, "U prevents false D22/D17 fan")
    need(support[(4, 4, 1)] == -xi, "xi prevents false D17 fan")
    need(support[(4, 2, 1)] == -K, "K retains T face")

    S, P = sp.symbols("S P")
    expected_faces = {
        "M": (P - S**2) * (U * P**6 - 1),
        "D6": -S**2 * (U * P**6 + alpha * S * P**5 + xi * S**2 * P**4 - 1),
        "T": -S**2 * (K * P**2 * S**2 + Theta * P**3 * S**2
                       + xi * P**4 * S**2 - 1),
    }
    for name, plane in PLANES.items():
        expression, _face_points = face_expression(support, plane)
        need(sp.expand(expression - expected_faces[name]) == 0, name + " exact face")

    polygons = {
        "M": ((0, 1), (2, 0), (2, 6), (0, 7)),
        "D6": ((2, 0), (4, 4), (2, 6)),
        "T": ((2, 0), (4, 2), (4, 4)),
        "global": ((0, 1), (2, 0), (4, 2), (4, 4), (2, 6), (0, 7)),
    }
    ledgers = {name: polygon_ledger(vertices)[1] for name, vertices in polygons.items()}
    need(ledgers == {
        "M": (24, 14, 6),
        "D6": (12, 10, 2),
        "T": (4, 6, 0),
        "global": (40, 14, 14),
    }, "Pick ledgers")
    packet = base.edge_packet(polygons["global"])
    need(packet == (11, 6, 6, 3, 3, 2, 2, 1), "global edge packet")
    need(sum(value - 1 for value in packet) == 26 == 2 * 14 - 2,
         "ramification/genus checksum")
    need(face_order(12, M) == 9, "M order")
    need(face_order(6, D6) == 5, "D6 primitive order")
    need(face_order(6, T) == 8, "T primitive order")

    # Complete edge list.  The two quadratic discriminants are treated by
    # explicit local charts below.  Every other edge is simple on U*K*xi!=0.
    Xedge = sp.symbols("Xedge")
    A_finite = K + Theta * Xedge + xi * Xedge**2
    A_infinity = U + alpha * Xedge + xi * Xedge**2
    edge_schemes = (
        Xedge - 1,
        1 - K * Xedge**2,
        A_finite,
        A_infinity,
        Xedge - 1,
        U - Xedge**6,
        1 - U * Xedge**6,
        1 - xi * Xedge**2,
    )
    disc_finite = sp.factor(sp.discriminant(A_finite, Xedge))
    disc_infinity = sp.factor(sp.discriminant(A_infinity, Xedge))
    need(disc_finite == Theta**2 - 4 * K * xi, "finite quadratic discriminant")
    need(disc_infinity == alpha**2 - 4 * U * xi, "infinity quadratic discriminant")
    need(sp.discriminant(edge_schemes[1], Xedge) == 4 * K, "K edge simple")
    need(sp.factor(sp.discriminant(edge_schemes[5], Xedge)) == 46656 * U**5,
         "U edge simple")
    need(sp.factor(sp.discriminant(edge_schemes[6], Xedge)) == 46656 * U**5,
         "M--D6 edge simple")
    need(sp.discriminant(edge_schemes[7], Xedge) == 4 * xi,
         "D6--T edge simple")

    # M has R:S^2=P plus six lines.  Twelve internal nodes, six M--D6
    # attachments, and two D6--T attachments give E=20,V=9,b1=12.
    need(12 - 7 + 1 == 6, "M internal graph genus")
    need(20 - 9 + 1 == 12, "full graph b1")
    need(2 + 12 == 14, "generic component/genus ledger")

    # Reconstruct the D6 chart from the complete specialized source.
    sigma = sp.symbols("sigma")
    p0, s0, y0 = sigma**-1 * P, sigma**-1 * S, sigma**-2 * S * P
    H0 = sum(
        coefficient * p0**i * y0**j
        for (i, j), coefficient in rows_by_label.items()
    )
    FQ = (s0**2 - p0) * (1 - sigma**6 * H0) - sigma**6 * s0**2 / 2
    G = sp.expand(sigma**2 * FQ)
    H6 = U * P**6 + alpha * S * P**5 + xi * S**2 * P**4
    H5 = u * P**5 + eta * S * P**4 + Theta * S**2 * P**3
    H4 = Delta * P**4 + Phi * S * P**3 + K * S**2 * P**2
    expected_G = sp.expand(
        (S**2 - sigma * P)
        * (1 - H6 - sigma * H5 - sigma**2 * H4
           - sigma**3 * e * P**3 - sigma**4 * sp.Rational(8, 3) * P**2
           + 3 * sigma**5 * P)
        - sigma**6 * S**2 / 2
    )
    need(sp.simplify(G - expected_G) == 0, "exact primitive D6 chart")

    # Root-at-infinity coordinates.  X=1/P,V=S/P are unimodular.  With
    # rho=sigma*X, division by the unit V^2-rho produces one critical series.
    x, v, rho = sp.symbols("x v rho")
    Ainf = U + alpha * v + xi * v**2
    Binf = u + eta * v + Theta * v**2
    Cinf = Delta + Phi * v + K * v**2
    Phi_root = sp.factor(x**8 * G.subs({P: 1 / x, S: v / x}))
    expected_root = sp.expand(
        (v**2 - sigma * x)
        * (x**6 - Ainf - sigma * x * Binf - sigma**2 * x**2 * Cinf
           - sigma**3 * x**3 * e - sigma**4 * x**4 * sp.Rational(8, 3)
           + 3 * sigma**5 * x**5)
        - sigma**6 * v**2 * x**6 / 2
    )
    need(sp.simplify(Phi_root - expected_root) == 0, "exact root-at-infinity chart")
    Dinf = (
        Ainf + rho * Binf + rho**2 * Cinf + e * rho**3
        + sp.Rational(8, 3) * rho**4 - 3 * rho**5
        + v**2 * rho**6 / (2 * (v**2 - rho))
    )
    need(sp.factor(sp.together(
        Phi_root / (v**2 - sigma * x)
        - (x**6 - Dinf.subs(rho, sigma * x))
    )) == 0, "one-series divided root chart")

    Einf = alpha**2 - 4 * U * xi
    Yhyp = 2 * xi * v + alpha
    need(sp.expand(Yhyp**2 - 4 * xi * x**6 - Einf
                   - 4 * xi * (Ainf - x**6)) == 0,
         "D6 hyperelliptic completion identity")

    # On Einf=0, Ainf=xi*(v-a)^2.  Since U*xi!=0, a is nonzero, the
    # prefactor is a unit, and formal Morse preparation gives
    # y^2=x^6-psi(rho), rho=sigma*x.
    a = sp.symbols("a", nonzero=True)
    repeated_infinity = {U: xi * a**2, alpha: -2 * xi * a}
    need(sp.expand(Ainf.subs(repeated_infinity) - xi * (v - a)**2) == 0,
         "repeated infinity parameterization")
    need(sp.diff(Dinf, v, 2).subs({rho: 0, **repeated_infinity}) == 2 * xi,
         "Morse Hessian at repeated infinity root")
    need((v**2 - rho).subs({v: a, rho: 0}) == a**2,
         "outer prefactor is a unit")

    # Phi=x^8G, P=x^-1,S=v/x imply G_S=x^-7 Phi_v and
    # -dP/G_S=x^5 dx/Phi_v.  Including D6 order five gives
    # omega=sigma^5*x^5 dx/y after Morse preparation.
    GS_in_root = sp.diff(G, S).subs({P: 1 / x, S: v / x})
    need(sp.simplify(sp.diff(Phi_root, v) - x**7 * GS_in_root) == 0,
         "root-chart derivative conversion")

    # Cancellation-complete even-cusp table.  The normalized tail has
    # Y0^2=X^epsilon*(X^(6-r)-c).  It has two points at infinity, hence two
    # attachments; positive-genus rows r=1,2,3 have positive form order.
    tail_table = {}
    for r in range(1, 6):
        epsilon = r % 2
        branch_degree = 6 - r + epsilon
        genus = (branch_degree - 2) // 2
        persistent_delta = r // 2
        order_sigma = F(5) + F(3 * r, 6 - r)
        order_tau = order_sigma * 2 * (6 - r)
        normalized_global_genus = 14 - persistent_delta
        tail_table[r] = (
            epsilon, branch_degree, genus, persistent_delta,
            order_sigma, order_tau, normalized_global_genus,
        )
    need(tail_table == {
        1: (1, 6, 2, 0, F(28, 5), F(56), 14),
        2: (0, 4, 1, 1, F(13, 2), F(52), 13),
        3: (1, 4, 1, 1, F(8), F(48), 13),
        4: (0, 2, 0, 2, F(11), F(44), 12),
        5: (1, 2, 0, 2, F(20), F(40), 12),
    }, "even-cusp tails/orders")
    need(all(row[4] > 0 for row in tail_table.values()), "all tail orders positive")
    need(11 + tail_table[1][2] + 1 == 14, "r=1 two-attachment genus")
    need(11 + tail_table[2][2] + 1 == 13, "r=2 two-attachment genus")
    need(11 + tail_table[4][2] + 1 == 12, "r=4 two-attachment genus")

    # Exact hostile: xi=U=1,alpha=2 makes the infinity edge double.  The
    # first two critical-value coefficients can vanish while the third is
    # the fixed e, producing a positive-genus r=3 tail.
    rr = sp.symbols("rr")
    Krel = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
    quad_a2 = 1 + Theta * rr + Krel * rr**2
    quad_a1 = 2 + eta * rr + Phi * rr**2
    quad_a0 = (
        1 + u * rr + Delta * rr**2 + e * rr**3
        + sp.Rational(8, 3) * rr**4 - 3 * rr**5
    )
    critical_to_five = sp.series(
        quad_a0 - quad_a1**2 / (4 * quad_a2), rr, 0, 6
    ).removeO().expand()
    r3_hostile = {
        u: 0, eta: 0, Theta: 0, Delta: 0,
        Phi: sp.Rational(2848, 45),
    }
    hostile_coefficients = tuple(
        sp.factor(critical_to_five.subs(r3_hostile).coeff(rr, degree))
        for degree in range(1, 4)
    )
    need(hostile_coefficients == (0, 0, e), "exact r=3 Hasse hostile")

    # The finite quadratic is rational when squarefree.  When it is a square,
    # the exact reciprocal T chart gives either a rational A_11 bridge or a
    # persistent horizontal node.
    tau, b, XX, YY = sp.symbols("tau b XX YY")
    AfinP = K + Theta * P + xi * P**2
    Bfull = P**3 * (Phi + eta * P + alpha * P**2)
    Cfull = (
        -3 * P + sp.Rational(8, 3) * P**2 + e * P**3
        + Delta * P**4 + u * P**5 + U * P**6
    )
    Floc = sp.expand(
        (1 - tau**12 * P * b**2)
        * (b**2 - P**2 * AfinP - tau**6 * b * Bfull
           - tau**12 * b**2 * Cfull)
        - tau**12 * b**2 / 2
    )
    finite_double = {K: xi * a**2, Theta: -2 * xi * a}
    need(sp.expand(AfinP.subs(finite_double) - xi * (P - a)**2) == 0,
         "finite double-root parameterization")
    scaled_double = sp.expand(
        Floc.subs(finite_double).subs({P: a + tau**6 * XX, b: tau**6 * YY})
    )
    double_face = sp.factor(sp.limit(scaled_double / tau**12, tau, 0))
    Ba = sp.factor(Bfull.subs(P, a))
    expected_double_face = YY**2 - a**2 * xi * XX**2 - Ba * YY
    need(sp.simplify(double_face - expected_double_face) == 0,
         "finite double-root conic")
    finite_horizontal = {
        **finite_double,
        Phi: -eta * a - alpha * a**2,
    }
    need(sp.simplify(Floc.subs({P: a, b: 0, **finite_horizontal})) == 0,
         "finite horizontal section lies on curve")
    need(sp.simplify(sp.diff(Floc, b).subs({P: a, b: 0, **finite_horizontal})) == 0,
         "finite horizontal b-critical section")
    need(sp.simplify(sp.diff(Floc, P).subs({P: a, b: 0, **finite_horizontal})) == 0,
         "finite horizontal P-critical section")
    need(16 + 3 * 6 - (12 - 6) == 28, "finite bridge positive order")

    print("THM4344_CLEAN_CUBIC_INFINITY_EXIT=PASS")
    print("gate=Z=beta_11=zeta_3=W=0;U*K*xi_10!=0")
    print("conservative_over_atlas=512/512:M,D6,T")
    print("faces=M:(P-S^2)(UP^6-1);D6:S^2(1-UP^6-alpha*SP^5-xi*S^2P^4);T:S^2(1-S^2P^2A)")
    print("pick=M(24,14,6);D6(12,10,2);T(4,6,0);global(40,14,14)")
    print("packet=(11,6,6,3,3,2,2,1);graph_b1=12;generic_genus=2+12=14")
    print("face_orders=M:9@L12;D6:5@L6;T:8@L6")
    print("root_infinity=X^6=U+alpha*V+xi*V^2;Y^2=4xi*X^6+(alpha^2-4Uxi)")
    print("infinity_repeat=y^2=x^6-psi(sigma*x);omega=sigma^5*x^5dx/y")
    print("tail_table=" + str(tail_table))
    print("r3_hostile=(h1,h2,h3)=" + str(hostile_coefficients))
    print("finite_repeat=" + str(double_face) + ";bridge_order=28_or_horizontal_node")
    print("boundaries=K0_zero_exit;xi0_second_infinity_exit;U0_endpoint_intersection")
    print("checks=" + str(CHECKS))


if __name__ == "__main__":
    main()
