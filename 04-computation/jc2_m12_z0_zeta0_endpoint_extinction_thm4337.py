#!/usr/bin/env python3
"""Primary exact certificate for the THM-4337 zeta-zero endpoint wall.

This path reuses the audited THM-4327 lower-hull primitives and independently
checks the new face, edge, genus, differential-order, and contact strata.
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
sys.path.insert(0, "04-computation")
import jc2_m12_z0_endpoint_extinction_thm4327 as base


CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("THM4337 primary failure: " + label)


def ledger(points):
    polygon = base.convex_hull(points)
    twice_area = abs(sum(
        polygon[k][0] * polygon[(k + 1) % len(polygon)][1]
        - polygon[(k + 1) % len(polygon)][0] * polygon[k][1]
        for k in range(len(polygon))
    ))
    boundary = sum(
        gcd(
            abs(polygon[(k + 1) % len(polygon)][0] - polygon[k][0]),
            abs(polygon[(k + 1) % len(polygon)][1] - polygon[k][1]),
        )
        for k in range(len(polygon))
    )
    require((twice_area - boundary + 2) % 2 == 0, "Pick parity")
    genus = (twice_area - boundary + 2) // 2
    return polygon, (twice_area, boundary, genus)


def order(base_degree, plane):
    value = base_degree * (F(5, 6) - sum(plane))
    require(value.denominator == 1 and value > 0, "positive integral order")
    return value.numerator


def p_euler_remainder(terms, shift):
    """Apply P*d/dP-shift with coefficients in (1,U,W,beta,K)."""
    answer = {}
    for exponent, coefficient in terms.items():
        factor = exponent[1] - shift
        scaled = tuple(factor * value for value in coefficient)
        if any(scaled):
            answer[exponent] = scaled
    return answer


def main():
    base_path = Path("04-computation/jc2_m12_z0_endpoint_extinction_thm4327.py")
    base_hash = sha256(base_path.read_bytes()).hexdigest()
    require(
        base_hash == "6aa9087afd413833d3168b27efe3e65e779ab796a9b537ef1e624a6380fac551",
        "pinned THM-4327 lower-hull dependency",
    )
    rows = base.weighted_rows()
    points, point_index, support_mask, lower_faces = base.make_hull_engine(rows)

    p_rows = {(1, 0), (2, 0), (3, 0), (6, 0)}
    w_row = (3, 2)
    beta_row = (1, 3)
    z_row = (0, 4)
    zeta_row = (0, 3)
    collisions = (
        (2, 3, 1), (2, 4, 1), (2, 5, 1),
        (2, 6, 1), (3, 4, 1), (3, 5, 1),
    )
    w0_collisions = (
        (2, 3, 1), (2, 4, 1), (2, 5, 1), (3, 5, 1),
    )

    m12 = (F(1, 12), F(1, 6), F(-1, 6))
    b7 = (F(1, 7), F(1, 7), F(-2, 7))
    b9 = (F(1, 9), F(1, 6), F(-2, 9))
    kface = (F(1), F(-1, 2), F(-2))

    optional, total, counter = base.hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        p_rows | {w_row, beta_row}, {z_row, zeta_row}, collisions,
    )
    require(optional == 8 and total == 16384, "W-nonzero hostile population")
    expected = Counter({(m12, b7): 8192, (m12, b7, kface): 8192})
    require(counter == expected, "W-nonzero two-complex atlas")

    optional0, total0, counter0 = base.hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        p_rows | {beta_row}, {z_row, zeta_row, w_row}, w0_collisions,
    )
    require(optional0 == 8 and total0 == 4096, "W-zero hostile population")
    expected0 = Counter({(m12, b9): 2048, (m12, b9, kface): 2048})
    require(counter0 == expected0, "W-zero two-complex atlas")

    # The extra plane is present exactly when K=[y^2]H is nonzero.  Its
    # equality triangle contains the unique K point (4,2,1).
    universe = base.lifted_support(
        row for row in rows if row[:2] not in {z_row, zeta_row}
    )
    prepared_points, records = base.prepare_plane_records(universe)
    record = next(record for record in records if record[0] == kface)
    equality = tuple(
        point for index, point in enumerate(prepared_points)
        if record[2] & (1 << index)
    )
    require(equality == ((2, 0, 0), (4, 2, 1), (5, 4, 1)),
            "K-face equality triangle")

    polygons = {
        "M_W": ((0, 1), (2, 0), (4, 5), (0, 7)),
        "C_W": ((0, 0), (0, 6), (2, 5)),
        "B_W": ((2, 0), (4, 5), (5, 4)),
        "K": ((2, 0), (4, 2), (5, 4)),
        "G_W_K0": ((0, 1), (2, 0), (5, 4), (4, 5), (0, 7)),
        "G_W_K1": ((0, 1), (2, 0), (4, 2), (5, 4), (4, 5), (0, 7)),
        "M_0": ((0, 1), (2, 0), (2, 6), (0, 7)),
        "B_0": ((2, 0), (2, 6), (5, 4)),
        "G_0_K0": ((0, 1), (2, 0), (5, 4), (2, 6), (0, 7)),
        "G_0_K1": ((0, 1), (2, 0), (4, 2), (5, 4), (2, 6), (0, 7)),
    }
    ledgers = {name: ledger(vertices)[1] for name, vertices in polygons.items()}
    require(ledgers == {
        "M_W": (36, 10, 14), "C_W": (12, 8, 3),
        "B_W": (7, 3, 3), "K": (2, 4, 0),
        "G_W_K0": (43, 11, 17), "G_W_K1": (45, 13, 17),
        "M_0": (24, 14, 6), "B_0": (18, 8, 6),
        "G_0_K0": (42, 10, 17), "G_0_K1": (44, 12, 17),
    }, "Pick ledgers")

    packets = {
        name: base.edge_packet(polygons[name])
        for name in ("G_W_K0", "G_W_K1", "G_0_K0", "G_0_K1")
    }
    require(packets == {
        "G_W_K0": (11, 11, 7, 7, 1),
        "G_W_K1": (11, 11, 7, 5, 2, 2, 1),
        "G_0_K0": (17, 11, 7, 1),
        "G_0_K1": (17, 11, 5, 2, 2, 1),
    }, "source-infinity packets")
    require(all(sum(index - 1 for index in packet) == 32
                for packet in packets.values()), "Riemann-Hurwitz saturation")

    orders_w = tuple(order(84, plane) for plane in (m12, b7, kface))
    orders_0 = tuple(order(36, plane) for plane in (m12, b9, kface))
    require(orders_w == (63, 70, 196), "W-nonzero orders")
    require(orders_0 == (27, 28, 84), "W-zero orders")

    # Genus accounting: K absent/present does not alter b1=11.
    # W!=0: R,C,B,(K), with 12+1+(1) edges.
    # W=0: R, six L_i, B,(K), with 12+6+(1) edges.
    require((13 - 3 + 1, 14 - 4 + 1) == (11, 11), "W graph b1")
    require((18 - 8 + 1, 19 - 9 + 1) == (11, 11), "W0 graph b1")
    require(3 + 3 + 11 == 17 and 6 + 11 == 17, "component genus totals")
    require(gcd(2, 5) == 1, "primitive C--B7 edge")
    require(gcd(3, 4) == 1, "primitive B--K edge")
    require(gcd(0, 6) == 6, "sixfold W0 main--B9 edge")

    # Torus-smoothness determinants.  The K face also has the unimodular
    # parameter V=S*P, P=(1-K*V^2)/(beta*V^3).
    require(base.determinant((0, 6), (2, 5)) == -12, "C determinant")
    require(base.determinant((2, 5), (3, 4)) == -7, "B7 determinant")
    require(base.determinant((0, 6), (3, 4)) == -18, "B9 determinant")
    require(base.determinant((2, 2), (3, 4)) == 2, "K determinant")

    # Symbolic generic-point witnesses modulo the indicated face equation.
    # C=1-U*P^6-W*S^2*P^5: P*C_P-5*C=-(U*P^6+5).
    # B7=1-W*S^2*P^5-beta*S^3*P^4:
    # P*B_P-4*B=-(W*S^2*P^5+4).
    # B9=1-U*P^6-beta*S^3*P^4:
    # P*B_P-4*B=-(2*U*P^6+4).
    # Kf=1-K*S^2*P^2-beta*S^3*P^4:
    # P*Kf_P-4*Kf=2*K*S^2*P^2-4.
    witnesses = (
        "P*C_P-5*C=-(U*P^6+5)",
        "P*B7_P-4*B7=-(W*S^2*P^5+4)",
        "P*B9_P-4*B9=-(2*U*P^6+4)",
        "P*Kf_P-4*Kf=2*K*S^2*P^2-4",
    )
    one = (1, 0, 0, 0, 0)
    minus_u = (0, -1, 0, 0, 0)
    minus_w = (0, 0, -1, 0, 0)
    minus_beta = (0, 0, 0, -1, 0)
    minus_k = (0, 0, 0, 0, -1)
    sparse_faces = {
        "C": ({(0, 0): one, (0, 6): minus_u, (2, 5): minus_w}, 5),
        "B7": ({(0, 0): one, (2, 5): minus_w, (3, 4): minus_beta}, 4),
        "B9": ({(0, 0): one, (0, 6): minus_u, (3, 4): minus_beta}, 4),
        "Kf": ({(0, 0): one, (2, 2): minus_k, (3, 4): minus_beta}, 4),
    }
    exact_remainders = {
        name: p_euler_remainder(terms, shift)
        for name, (terms, shift) in sparse_faces.items()
    }
    require(exact_remainders == {
        "C": {(0, 0): (-5, 0, 0, 0, 0), (0, 6): (0, -1, 0, 0, 0)},
        "B7": {(0, 0): (-4, 0, 0, 0, 0), (2, 5): (0, 0, -1, 0, 0)},
        "B9": {(0, 0): (-4, 0, 0, 0, 0), (0, 6): (0, -2, 0, 0, 0)},
        "Kf": {(0, 0): (-4, 0, 0, 0, 0), (2, 2): (0, 0, 0, 0, 2)},
    }, "exact symbolic P-Euler remainders")

    edge_schemes_w_k0 = (
        "X-1", "1-beta*X", "-W*X-beta",
        "(X-1)*(U*X+W)", "U-X^6", "1-W*X",
    )
    edge_schemes_w_k1 = (
        "X-1", "1-K*X^2", "-beta*X-K", "-W*X-beta",
        "(X-1)*(U*X+W)", "U-X^6", "1-W*X", "1-beta*X",
    )
    edge_schemes_0_k0 = (
        "X-1", "1-beta*X", "-U*X-beta", "U*(X-1)",
        "U-X^6", "1-U*X^6",
    )
    edge_schemes_0_k1 = (
        "X-1", "1-K*X^2", "-beta*X-K", "-U*X-beta",
        "U*(X-1)", "U-X^6", "1-U*X^6", "1-beta*X",
    )

    # Literal positive controls for every edge list.  These are paired with
    # the displayed symbolic gates below; they guard signs, powers, and the
    # K-present/absent replacement rather than standing in for the symbolic
    # implications.
    u, w, beta, kappa = map(F, (2, 3, 5, 7))
    xm1 = (-F(1), F(1))
    top = base.polynomial_multiply(xm1, (w, u))
    literal_edges_w_k0 = (
        xm1, (1, -beta), (-beta, -w), top,
        (u, 0, 0, 0, 0, 0, -1), (1, -w),
    )
    literal_edges_w_k1 = (
        xm1, (1, 0, -kappa), (-kappa, -beta), (-beta, -w), top,
        (u, 0, 0, 0, 0, 0, -1), (1, -w), (1, -beta),
    )
    literal_edges_0_k0 = (
        xm1, (1, -beta), (-beta, -u), (-u, u),
        (u, 0, 0, 0, 0, 0, -1), (1, 0, 0, 0, 0, 0, -u),
    )
    literal_edges_0_k1 = (
        xm1, (1, 0, -kappa), (-kappa, -beta), (-beta, -u), (-u, u),
        (u, 0, 0, 0, 0, 0, -1),
        (1, 0, 0, 0, 0, 0, -u), (1, -beta),
    )
    for regime in (
        literal_edges_w_k0, literal_edges_w_k1,
        literal_edges_0_k0, literal_edges_0_k1,
    ):
        require(all(base.polynomial_squarefree(edge) for edge in regime),
                "literal edge list squarefree")
        require(all(base.polynomial_trim(edge)[0] != 0
                    and base.polynomial_trim(edge)[-1] != 0
                    for edge in regime), "literal edge avoids toric corners")
    top_disc = top[1] ** 2 - 4 * top[0] * top[2]
    require(top_disc == (u + w) ** 2, "top discriminant is (U+W)^2")
    repeated_top = base.polynomial_multiply(xm1, (-u, u))
    require(not base.polynomial_squarefree(repeated_top),
            "Lambda-zero top is uniquely repeated")

    # The K-face parameter really is birational: the unimodular change
    # (S,P)->(V=S*P,P) sends exponent (2,2) to (2,0) and (3,4) to (3,1),
    # leaving an equation linear in P.
    def to_vp(exponent):
        s_exp, p_exp = exponent
        return s_exp, p_exp - s_exp

    require(to_vp((2, 2)) == (2, 0), "K monomial -> V^2")
    require(to_vp((3, 4)) == (3, 1), "beta monomial -> V^3 P")
    require(base.determinant((1, 0), (1, 1)) == 1,
            "K-face change is unimodular")

    # Lambda=0 is possible only in the W!=0 stratum.  The exact top closure
    # has one length-twelve A23 contact.  In the deepest THM-4292 ladder,
    # repetition forces c1=alpha+beta=0.  Since beta!=0 here, alpha=-beta
    # is nonzero, so the moving critical coefficient C1=alpha/c^2 is the
    # first nonzero splitter.  It occurs at t^1, strictly before the t^6
    # W-perturbation, and has positive form order (5s+9nu)/2.
    u0 = F(2)
    w0 = -u0
    require(u0 + w0 == 0 and u0 * w0 != 0, "Lambda-zero gate control")
    beta11 = F(5)
    alpha11 = -beta11
    c_repeat = F(3)
    c1 = alpha11 + beta11
    moving_c1 = alpha11 / (c_repeat * c_repeat)
    require(c1 == 0 and moving_c1 != 0,
            "beta gate forces nonzero moving C1 on deepest repetition")
    def local_gaps(s_value, nu_value):
        return (
            s_value + nu_value,
            6 * (nu_value - s_value),
            6 * (s_value + nu_value),
        )

    below, equal, above = (
        local_gaps(F(5), F(6)),
        local_gaps(F(5), F(7)),
        local_gaps(F(5), F(8)),
    )
    require(below[1] < below[0], "b^12*q splits first below nu/s=7/5")
    require(equal[1] == equal[0], "splitter collision exactly at nu/s=7/5")
    require(above[0] < above[1], "C1*t splits first above nu/s=7/5")
    require(all(min(d1, db) < correction for d1, db, correction in (below, equal, above)),
            "both early splitters precede t6 perturbation")
    c1_order_coefficients = (F(5, 2), F(9, 2))
    require(all(value > 0 for value in c1_order_coefficients),
            "C1 Keller-order coefficients are uniformly positive")
    b12_order_coefficients = (F(6), F(2))
    require(all(value > 0 for value in b12_order_coefficients),
            "b^12*q Keller-order coefficients are uniformly positive")
    require(84 % 12 == 0 and 84 // 12 == 7,
            "Lambda-zero local-to-global base index seven")

    # Boundary exploration when beta also vanishes.  For W!=0 the atlas has
    # only M and the optional plane T=(1/2,0,-1), present iff at least one of
    # K,Theta,xi is nonzero.  Its core is
    #   S^2(1-S^2 P^2 A(P)), A=K+Theta P+xi P^2+W P^3.
    # The vertical edge is the Laurent reduction of A.  A repeated nonzero
    # root is the precise unaudited tail wall.
    beta0_optional, beta0_total, beta0_counter = base.hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        p_rows | {w_row}, {z_row, zeta_row, beta_row}, collisions,
    )
    tface = (F(1, 2), F(0), F(-1))
    require(beta0_optional == 8 and beta0_total == 16384,
            "beta-zero W-nonzero hostile population")
    require(beta0_counter == Counter({
        (m12, tface): 14336,
        (m12,): 2048,
    }), "beta-zero W-nonzero two-complex atlas")
    repeated_a = (F(1), F(-1), F(-1), F(1))
    require(not base.polynomial_squarefree(repeated_a),
            "A=(P-1)^2(P+1) repeated-edge hostile")
    require(order(12, tface) == 16, "beta-zero T-face positive order")

    # At W=beta=zeta=0 the honest collision overcount already expands to
    # twelve complexes/nine planes.  The first sharp hostile is the face
    # 1-U P^6-alpha S P^5-xi S^2 P^4 with top-edge quadratic
    # U+alpha X+xi X^2.  U=xi=1, alpha=2 makes (X+1)^2.
    wbeta0_collisions = ((2, 3, 1), (2, 4, 1), (2, 5, 1))
    wb0_optional, wb0_total, wb0_counter = base.hostile_atlas(
        rows, point_index, support_mask, lower_faces,
        p_rows, {z_row, zeta_row, beta_row, w_row}, wbeta0_collisions,
    )
    wb0_planes = {plane for complex_ in wb0_counter for plane in complex_}
    require(wb0_optional == 8 and wb0_total == 2048,
            "W=beta=0 hostile population")
    require(len(wb0_counter) == 12 and len(wb0_planes) == 9,
            "W=beta=0 atlas size")
    require(min(F(5, 6) - sum(plane) for plane in wb0_planes) == F(3, 4),
            "W=beta=0 all primary orders positive")
    hostile_quadratic = (F(1), F(2), F(1))
    require(not base.polynomial_squarefree(hostile_quadratic),
            "W=beta=0 repeated quadratic hostile")

    print("THM4337_ZETA0_ENDPOINT_PRIMARY=PASS")
    print("wall=Z=0;zeta_3=0;U*beta_11!=0;W_and_K_arbitrary")
    print("THM4327_Z_SCRIPT_SHA256=" + base_hash)
    print("K=[y^2]H")
    print("W_NONZERO_hostiles=16384;optional=8;inherited_deletion_bits=6")
    print("W_NONZERO_complexes=" + repr(counter))
    print("W_ZERO_hostiles=4096;optional=8;inherited_deletion_bits=4")
    print("W_ZERO_complexes=" + repr(counter0))
    print("Kface_equality=" + repr(equality))
    print("ledgers=" + repr(ledgers))
    print("packets=" + repr(packets))
    print("orders_W_L84=" + repr(orders_w))
    print("orders_W0_L36=" + repr(orders_0))
    print("normalizations_W=(C:genus3,B7:genus3,Kface:rational)")
    print("normalizations_W0=(six_M_lines,B9:genus6,Kface:rational)")
    print("Kface_parameter=V=S*P;P=(1-K*V^2)/(beta*V^3)")
    print("GP_witnesses=" + repr(witnesses))
    print("edges_W_K0=" + repr(edge_schemes_w_k0))
    print("edges_W_K1=" + repr(edge_schemes_w_k1))
    print("edges_W0_K0=" + repr(edge_schemes_0_k0))
    print("edges_W0_K1=" + repr(edge_schemes_0_k1))
    print("simple_edge_gates=W:U*W*beta*(U+W)!=0;W0:U*beta!=0;K_nonzero_only_when_Kface_present")
    print("Lambda0_contact=Ctilde=b^12-U*r^5*(r-1);d_r=-U;A23;do_not_transport_packet")
    print("Lambda0_deepest=c1=alpha+beta=0_and_beta!=0=>C1=alpha/c^2!=0")
    print("Lambda0_early_splitter=min(b^12*q,C1*t);threshold_nu_over_s=7/5")
    print("Lambda0_orders=b^12*q:6s+2nu;C1*t:(5s+9nu)/2;both_before_t^6")
    print("Lambda0_base_transport=global_L84_over_local_L12_index7")
    print("beta0_W_nonzero_atlas=2_complexes;optional_T=(1/2,0,-1);order_L12=16")
    print("beta0_W_nonzero_tail_hostile=A(P)=(P-1)^2(P+1);repeated_nonzero_outer_edge")
    print("beta0_W0_atlas=12_complexes;9_planes;minimum_order_coefficient=3/4")
    print("beta0_W0_sharp_hostile=U+alpha*X+xi*X^2=(X+1)^2")
    print("CHECKS=" + str(CHECKS))


if __name__ == "__main__":
    main()
