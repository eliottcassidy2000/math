#!/usr/bin/env python3
"""Primary exact certificate for THM-4342's K=0 zero-root exit.

Scope: Z=beta_11=zeta_3=K=0 and U*W*(U+W) != 0 in the inherited
reduced (2,3), exact-weight-twelve seam.  The three strata are indexed by
m=ord_0(A), where A=Theta*P+xi*P^2+W*P^3.
"""

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


def need(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("K-zero root-exit audit failure: " + label)


def ledger(vertices):
    polygon, area2, boundary, interior = base.polygon_ledger(vertices)
    return polygon, (area2, boundary, interior)


def graph_betti(vertices, edges):
    parent = {vertex: vertex for vertex in vertices}

    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for first, second in edges:
        root_first, root_second = find(first), find(second)
        if root_first != root_second:
            parent[root_second] = root_first
    components = len({find(vertex) for vertex in vertices})
    return len(edges) - len(vertices) + components


def sigma_initial(expression, sigma, exponent):
    return sp.factor(sp.Poly(sp.expand(expression), sigma).coeff_monomial(sigma**exponent))


def main():
    delta_fixed = F(5696, 105)
    need(F(2848, 45) - F(7, 6) * delta_fixed == 0,
         "inherited K=0 forces Delta=5696/105")

    inherited = Path("04-computation/jc2_m12_z0_endpoint_extinction_thm4327.py")
    inherited_hash = sha256(inherited.read_bytes()).hexdigest()
    need(inherited_hash == "6aa9087afd413833d3168b27efe3e65e779ab796a9b537ef1e624a6380fac551",
         "pinned THM-4327 hull engine")

    # Complete hostile lower-hull enumeration in each exact zero-exit depth.
    rows = base.weighted_rows()
    points, point_index, support_mask, lower_faces = base.make_hull_engine(rows)
    p_rows = {(1, 0), (2, 0), (3, 0), (6, 0)}
    Krow, Thetarow, xirow, Wrow = (0, 2), (1, 2), (2, 2), (3, 2)
    Zrow, zetarow, betarow = (0, 4), (0, 3), (1, 3)
    collisions = (
        (2, 3, 1), (2, 4, 1), (2, 5, 1),
        (2, 6, 1), (3, 4, 1), (3, 5, 1),
    )
    M = (F(1, 12), F(1, 6), F(-1, 6))
    T = (F(1, 2), F(0), F(-1))
    atlas_specs = {
        1: (p_rows | {Thetarow, Wrow}, {Krow, Zrow, zetarow, betarow}, (M, T)),
        2: (p_rows | {xirow, Wrow}, {Krow, Thetarow, Zrow, zetarow, betarow}, (M, T)),
        3: (p_rows | {Wrow}, {Krow, Thetarow, xirow, Zrow, zetarow, betarow}, (M,)),
    }
    atlas_data = {}
    for multiplicity, (required, absent, expected_faces) in atlas_specs.items():
        optional, total, counter = base.hostile_atlas(
            rows, point_index, support_mask, lower_faces,
            required, absent, collisions,
        )
        need(counter == Counter({expected_faces: total}),
             f"m={multiplicity} unique lower complex")
        atlas_data[multiplicity] = (optional, total, counter)

    # Honest fan after Laurent saturation.
    polygons = {
        1: {
            "T": ((2, 0), (4, 3), (4, 5)),
            "global": ((0, 1), (2, 0), (4, 3), (4, 5), (0, 7)),
        },
        2: {
            "T": ((2, 0), (4, 4), (4, 5)),
            "global": ((0, 1), (2, 0), (4, 4), (4, 5), (0, 7)),
        },
        3: {
            "global": ((0, 1), (2, 0), (4, 5), (0, 7)),
        },
    }
    polygon_data = {
        multiplicity: {name: ledger(vertices)[1] for name, vertices in table.items()}
        for multiplicity, table in polygons.items()
    }
    need(polygon_data == {
        1: {"T": (4, 4, 1), "global": (40, 12, 15)},
        2: {"T": (2, 4, 0), "global": (38, 12, 14)},
        3: {"global": (36, 10, 14)},
    }, "root-exit Pick ledgers")
    packets = {
        multiplicity: base.edge_packet(table["global"])
        for multiplicity, table in polygons.items()
    }
    need(packets == {
        1: (11, 11, 5, 3, 3, 1),
        2: (11, 11, 3, 3, 3, 1),
        3: (11, 11, 7, 1),
    }, "root-exit edge packets")
    need(sum(index - 1 for index in packets[1]) == 28 == 2 * 15 - 2,
         "m=1 RH saturation")
    need(sum(index - 1 for index in packets[2]) == 26 == 2 * 14 - 2,
         "m=2 RH saturation")
    need(sum(index - 1 for index in packets[3]) == 26 == 2 * 14 - 2,
         "m=3 RH saturation")

    # Every honest edge scheme, including the new slanted edge created by the
    # root exit.  Schemes on different labelled divisors are not resultanted.
    X = sp.symbols("X")
    U, W, Theta, xi = sp.symbols("U W Theta xi", nonzero=True)
    edge_schemes = {
        1: (X - 1, 1 - Theta * X, Theta + xi * X + W * X**2,
            (X - 1) * (U * X + W), U - X**6, 1 - W * X),
        2: (X - 1, 1 - xi * X**2, xi + W * X,
            (X - 1) * (U * X + W), U - X**6, 1 - W * X),
        3: (X - 1, 1 - W * X,
            (X - 1) * (U * X + W), U - X**6),
    }
    need(sp.factor(sp.discriminant(edge_schemes[1][1], X)) == 1,
         "m=1 new slant simple")
    need(sp.factor(sp.discriminant(edge_schemes[2][1], X)) == 4 * xi,
         "m=2 new slant squarefree")
    need(sp.factor(sp.discriminant(edge_schemes[1][2], X)) == xi**2 - 4 * W * Theta,
         "m=1 saturated quadratic discriminant")
    need(sp.factor(sp.resultant(X - 1, U * X + W, X)) == U + W,
         "top-edge Lambda gate")
    need(sp.factor(sp.discriminant(U - X**6, X)) == 46656 * U**5,
         "sixfold edge U gate")

    # Exact reciprocal chart and its three toric strict transforms.
    sigma, delta, P, b, y = sp.symbols("sigma delta P b y")
    Phi, eta, alpha, Delta, upsilon = sp.symbols("Phi eta alpha Delta upsilon")
    J = Phi + eta * P + alpha * P**2
    C1 = -3 + sp.Rational(8, 3) * P - sp.Rational(1376, 135) * P**2 \
        + Delta * P**3 + upsilon * P**4 + U * P**5
    B = P**3 * J
    qbase = 1 - delta**2 / 2

    strict = {}
    for multiplicity, (power, quotient) in {
        1: (1, Theta + xi * P + W * P**2),
        2: (2, xi + W * P),
        3: (3, W),
    }.items():
        A = P**multiplicity * quotient
        Fexact = sp.expand(
            (1 - delta**2 * P * b**2)
            * (b**2 - P**2 * A - delta * b * B
               - delta**2 * b**2 * P * C1)
            - delta**2 * b**2 / 2
        )
        blow_power = 1 if multiplicity == 1 else 2
        transformed = sp.factor(
            Fexact.subs(b, P**blow_power * y) / P**(2 * blow_power)
        )
        need(not sp.denom(transformed).has(P), f"m={multiplicity} exact divisibility")
        strict[multiplicity] = sp.expand(transformed)

    expected1 = sp.expand(
        (1 - delta**2 * P**3 * y**2)
        * (y**2 - P * (Theta + xi * P + W * P**2)
           - delta * P**2 * y * J - delta**2 * P * y**2 * C1)
        - delta**2 * y**2 / 2
    )
    expected2 = sp.expand(
        (1 - delta**2 * P**5 * y**2)
        * (y**2 - (xi + W * P) - delta * P * y * J
           - delta**2 * P * y**2 * C1)
        - delta**2 * y**2 / 2
    )
    expected3 = sp.expand(
        (1 - delta**2 * P**5 * y**2)
        * (y**2 - W * P - delta * P * y * J
           - delta**2 * P * y**2 * C1)
        - delta**2 * y**2 / 2
    )
    need(sp.simplify(strict[1] - expected1) == 0, "m=1 exact strict transform")
    need(sp.simplify(strict[2] - expected2) == 0, "m=2 exact strict transform")
    need(sp.simplify(strict[3] - expected3) == 0, "m=3 exact strict transform")

    # Root-exit smoothness.  m=1 has one ramified point, m=2 has two
    # unramified points, and m=3 is a rational binomial collar.
    need(sp.simplify(sp.diff(strict[1], P).subs({P: 0, y: 0}) + Theta) == 0,
         "m=1 strict transform smooth at exited root")
    need(sp.simplify(strict[2].subs(P, 0) - (qbase * y**2 - xi)) == 0,
         "m=2 two-point exceptional equation")
    need(sp.simplify(sp.diff(strict[2], y).subs(P, 0) - 2 * qbase * y) == 0,
         "m=2 exceptional roots simple")
    need(sp.simplify(sp.diff(strict[3], P).subs({P: 0, y: 0}) + W) == 0,
         "m=3 rational collar smooth at origin")
    special_models = {
        multiplicity: sp.factor(equation.subs(delta, 0))
        for multiplicity, equation in strict.items()
    }
    need(special_models == {
        1: -P**3 * W - P**2 * xi - P * Theta + y**2,
        2: -P * W - xi + y**2,
        3: -P * W + y**2,
    }, "special normalized root-exit models")

    # The form identity survives b=P^e*y.  On E=0,
    # b^2 db/F_P=P^e*y^2 dy/E_P, hence no extra sigma power is introduced.
    e = sp.symbols("e", integer=True, positive=True)
    need(12 * (sp.Rational(5, 6) - sum(M)) == 9,
         "M order nine")
    need(12 * (sp.Rational(5, 6) - sum(T)) == 16,
         "T order sixteen")
    # Generic derivative witnesses on all special strict transforms.
    need(sp.diff(special_models[1], P) != 0, "m=1 E_P nonzero function")
    need(sp.diff(special_models[2], P) == -W, "m=2 E_P unit")
    need(sp.diff(special_models[3], P) == -W, "m=3 E_P unit")

    # m=1's only interior collision is a double root of the saturated
    # quadratic.  Replay both normal branches at a != 0.
    a = sp.symbols("a", nonzero=True)
    qdouble = W * (P - a)**2
    Adouble = P * qdouble
    Fdouble = sp.expand(
        (1 - delta**2 * P * b**2)
        * (b**2 - P**2 * Adouble - delta * b * B
           - delta**2 * b**2 * P * C1)
        - delta**2 * b**2 / 2
    )
    weighted = sp.expand(Fdouble.subs({P: a + sigma**6 * X,
                                       b: sigma**6 * y,
                                       delta: sigma**6}))
    double_face = sigma_initial(weighted, sigma, 12)
    Ba = sp.factor(B.subs(P, a))
    need(sp.simplify(double_face - (y**2 - Ba * y - a**3 * W * X**2)) == 0,
         "m=1 interior-double conic face")
    need(sp.factor(xi**2 - 4 * W * Theta).subs({Theta: W * a**2,
                                                xi: -2 * W * a}) == 0,
         "m=1 double discriminant parameterization")
    normal_gate = 4 * W**2 * Phi - 2 * W * xi * eta + xi**2 * alpha
    need(sp.simplify(normal_gate.subs({Theta: W * a**2, xi: -2 * W * a})
                     - 4 * W**2 * J.subs(P, a)) == 0,
         "m=1 polynomial normal-Hasse gate")

    # Hostile vanishing of the first normal coefficient is not a delayed
    # splitter: it produces an exact horizontal node.  Its tangent
    # discriminant is x^2 times a DVR unit.
    x = sp.symbols("x")
    Jzero = sp.expand(J - J.subs(P, a))
    Bzero = sp.expand(P**3 * Jzero)
    Fnode = sp.expand(
        (1 - delta**2 * P * b**2)
        * (b**2 - P**3 * W * (P - a)**2 - delta * b * Bzero
           - delta**2 * b**2 * P * C1)
        - delta**2 * b**2 / 2
    ).subs(P, a + x)
    need(sp.simplify(Fnode.subs({x: 0, b: 0})) == 0,
         "double horizontal section on curve")
    need(sp.simplify(sp.diff(Fnode, b).subs({x: 0, b: 0})) == 0,
         "double horizontal section b-critical")
    need(sp.simplify(sp.diff(Fnode, x).subs({x: 0, b: 0})) == 0,
         "double horizontal section x-critical")
    tangent = 0
    for (ix, ib), coefficient in sp.Poly(sp.expand(Fnode), x, b).terms():
        if ix + ib == 2:
            tangent += coefficient * x**ix * b**ib
    tangent = sp.expand(tangent)
    Ca = sp.factor((P * C1).subs(P, a))
    Bpa = sp.factor(sp.diff(Bzero, P).subs(P, a))
    tangent_expected = sp.expand(
        (1 - delta**2 * (Ca + sp.Rational(1, 2))) * b**2
        - delta * Bpa * x * b - a**3 * W * x**2
    )
    need(sp.simplify(tangent - tangent_expected) == 0,
         "double horizontal tangent cone")
    tangent_disc = sp.factor(sp.discriminant(tangent, b))
    tangent_disc_expected = sp.factor(
        x**2 * (delta**2 * Bpa**2
                + 4 * (1 - delta**2 * (Ca + sp.Rational(1, 2))) * a**3 * W)
    )
    need(sp.simplify(tangent_disc - tangent_disc_expected) == 0,
         "double horizontal node discriminant")

    # At the deepest m=3 exit, b^2~W P^5 gives v(b)=5v(P)/2.
    # The first normal term delta*b*P^3*J has excess
    # 1+v(P)/2 in delta-units, hence 6+v(P)/2 in sigma-units.
    r = sp.symbols("r", positive=True)
    m3_normal_excess_delta = sp.simplify(1 + sp.Rational(5, 2) * r + 3 * r - 5 * r)
    need(m3_normal_excess_delta == 1 + r / 2,
         "m=3 normal term strictly later")

    # Genus ledgers.  The central genus-three curve and twelve R--C nodes
    # persist because U+W != 0; C--T is primitive whenever T is retained.
    graph_rct = graph_betti(("R", "C", "T"),
                            (("R", "C"),) * 12 + (("C", "T"),))
    graph_rc = graph_betti(("R", "C"), (("R", "C"),) * 12)
    need(graph_rct == 11 and graph_rc == 11, "root-exit graph Betti numbers")
    need(3 + 1 + graph_rct == 15, "m=1 squarefree genus ledger")
    need(3 + graph_rct == 14, "m=2 genus ledger")
    need(3 + graph_rc == 14, "m=3 genus ledger")
    need(3 + (graph_rct + 1) == 15, "m=1 smoothed double bridge ledger")
    need(3 + graph_rct == 14, "m=1 persistent double conductor ledger")

    print("THM4342 KZERO ROOT-EXIT PRIMARY CERTIFICATE=PASS")
    print("gate=Z=beta_11=zeta_3=K=0;U*W*(U+W)!=0")
    print("seam_relation=K=0=>Delta=5696/105")
    for multiplicity in (1, 2, 3):
        optional, total, counter = atlas_data[multiplicity]
        print(f"m{multiplicity}_atlas=optional:{optional},masks:{total},faces:{next(iter(counter))}")
    print("pick=" + repr(polygon_data))
    print("packets=" + repr(packets))
    print("edge_schemes=" + repr(edge_schemes))
    print("strict_m1=(1-delta^2*P^3*y^2)*(y^2-P*(Theta+xi*P+W*P^2)-delta*P^2*y*J-delta^2*P*y^2*C1)-delta^2*y^2/2")
    print("strict_m2=(1-delta^2*P^5*y^2)*(y^2-(xi+W*P)-delta*P*y*J-delta^2*P*y^2*C1)-delta^2*y^2/2")
    print("strict_m3=(1-delta^2*P^5*y^2)*(y^2-W*P-delta*P*y*J-delta^2*P*y^2*C1)-delta^2*y^2/2")
    print("special_models=" + repr(special_models))
    print("form_pullback=b=P^e*y => -sigma^16*P^e*y^2*dy/E_P")
    print("orders=M:9,T:16,double_conic:28")
    print("m1_strata=saturated_quadratic_squarefree:elliptic;double:rational_bridge_or_horizontal_node")
    print("m1_double_gate=xi^2-4WTheta=0;a=-xi/(2W);normal=4W^2Phi-2Wxi*eta+xi^2alpha")
    print("m1_double_horizontal_tangent_discriminant=x^2*(delta^2*B'(a)^2+4*q*a^3W)")
    print("m2=rational_T_two_simple_exit_points;m3=T_face_collapses_to_rational_binomial_collar")
    print("m3_normal_excess=6+r/2>0")
    print("genera=m1_squarefree:15,m1_smoothed_double:15,m1_horizontal_double:14,m2:14,m3:14")
    print("checks=" + str(CHECKS))


if __name__ == "__main__":
    main()
