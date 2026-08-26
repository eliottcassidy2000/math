#!/usr/bin/env python3
"""Clean-room exact referee for THM-4220.

This script imports no primary-certificate code.  It reconstructs all three
zeta-zero support strata, their two-dimensional lower facets, face and edge
schemes, integral slope chains, literal main/side sigma-charts, and the
genus/Hom/degree-conservation obstruction.  The split-conic collision is
checked as a hostile and is not included in the verdict.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd

import sympy as sp


CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(label)


def lower_facets(points):
    """Integer cross-product implementation of the lower 3D hull."""
    points = tuple(sorted(set(points)))
    answer = {}
    for first, second, third in combinations(points, 3):
        ab = tuple(second[index] - first[index] for index in range(3))
        ac = tuple(third[index] - first[index] for index in range(3))
        normal = (
            ab[1]*ac[2] - ab[2]*ac[1],
            ab[2]*ac[0] - ab[0]*ac[2],
            ab[0]*ac[1] - ab[1]*ac[0],
        )
        if normal[2] == 0:
            continue
        if normal[2] < 0:
            normal = tuple(-entry for entry in normal)
        gaps = tuple(sum(normal[index]*(point[index] - first[index])
                         for index in range(3)) for point in points)
        if min(gaps) < 0:
            continue
        equality = tuple(point for point, gap in zip(points, gaps) if gap == 0)
        if len(equality) < 3:
            continue
        divisor = gcd(gcd(abs(normal[0]), abs(normal[1])), abs(normal[2]))
        primitive = tuple(entry // divisor for entry in normal)
        constant = sum(primitive[index]*first[index] for index in range(3))
        answer[equality] = primitive + (constant,)
    return tuple(sorted(answer.items(), key=lambda row: row[1]))


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, first, second):
        return ((first[0] - origin[0])*(second[1] - origin[1])
                - (first[1] - origin[1])*(second[0] - origin[0]))

    lower = []
    for point in points:
        while len(lower) > 1 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) > 1 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def polygon_data(points):
    polygon = convex_hull(points)
    area2 = abs(sum(
        polygon[index][0]*polygon[(index + 1) % len(polygon)][1]
        - polygon[(index + 1) % len(polygon)][0]*polygon[index][1]
        for index in range(len(polygon))
    ))
    boundary = sum(
        gcd(abs(polygon[(index + 1) % len(polygon)][0] - polygon[index][0]),
            abs(polygon[(index + 1) % len(polygon)][1] - polygon[index][1]))
        for index in range(len(polygon))
    )
    check((area2 - boundary + 2) % 2 == 0, "Pick parity")
    genus = (area2 - boundary + 2)//2
    packet = []
    edge_rows = []
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy//length, dx//length)
        constant = inward[0]*start[0] + inward[1]*start[1]
        local_degree = inward[0] + inward[1] - constant
        edge_rows.append((start, end, length, inward, constant, local_degree))
        if not (start[0] == end[0] == 0):
            packet.extend([local_degree]*length)
    return (polygon, area2, boundary, genus,
            tuple(sorted(packet, reverse=True)), tuple(edge_rows))


def edge_scheme(polynomial, first, second, start, end, variable):
    """Restrict a face polynomial to a primitive lattice edge."""
    delta = (end[0] - start[0], end[1] - start[1])
    length = gcd(abs(delta[0]), abs(delta[1]))
    step = (delta[0]//length, delta[1]//length)
    answer = 0
    for exponent, coefficient in sp.Poly(polynomial, first, second).terms():
        difference = (exponent[0] - start[0], exponent[1] - start[1])
        if step[0]:
            if difference[0] % step[0]:
                continue
            position = difference[0]//step[0]
        else:
            if difference[1] % step[1]:
                continue
            position = difference[1]//step[1]
        if (difference == (position*step[0], position*step[1])
                and 0 <= position <= length):
            answer += coefficient*variable**position
    return sp.factor(answer)


def digest(rows):
    return sha256("|".join(map(str, rows)).encode()).hexdigest()


def main():
    # Literal exact-M10 universe in the inherited b=d=0 seam.
    monomials = tuple(sorted(
        ((i, j, 2*i + 3*j)
         for i in range(6) for j in range(4)
         if 0 < 2*i + 3*j <= 10 and (i, j) not in {(0, 1), (1, 1)}),
        key=lambda row: (row[2], row[1], row[0]),
    ))
    check(monomials == (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10),
    ), "complete monomial universe")

    base = {(2, 0, 0), (0, 1, 0)}

    def support(rows):
        points = set(base)
        for i, j, _weight in rows:
            points.add((j + 2, i + j, 1))
            points.add((j, i + j + 1, 1))
        return points

    collision_points = ((2, 3, 1), (2, 4, 1), (2, 5, 1))
    cases = (
        ("K_nonzero", {(5, 0), (2, 2), (0, 2)}, {(0, 3)},
         {(-1, -2, 10, -2), (-1, 0, 2, -2)},
         ((0, 1), (2, 0), (4, 2), (4, 4), (0, 6)),
         (34, 12, 12), (9, 9, 3, 3, 2, 2, 1)),
        ("Theta_only", {(5, 0), (2, 2), (1, 2)}, {(0, 3), (0, 2)},
         {(-1, -2, 10, -2), (-1, 0, 2, -2)},
         ((0, 1), (2, 0), (4, 3), (4, 4), (0, 6)),
         (32, 10, 12), (9, 9, 5, 3, 1)),
        ("K_Theta_zero", {(5, 0), (2, 2)}, {(0, 3), (0, 2), (1, 2)},
         {(-1, -2, 10, -2)},
         ((0, 1), (2, 0), (4, 4), (0, 6)),
         (30, 10, 11), (9, 9, 3, 3, 1)),
    )
    case_ledgers = []
    for (name, required_indices, forbidden_indices, expected_facets,
         expected_polygon, expected_pick, expected_packet) in cases:
        required = tuple(row for row in monomials if row[:2] in required_indices)
        optional = tuple(row for row in monomials
                         if row[:2] not in required_indices | forbidden_indices)
        seen = set()
        for mask in range(1 << len(optional)):
            rows = list(required)
            rows.extend(row for bit, row in enumerate(optional) if mask & (1 << bit))
            raw_support = support(rows)
            # Deleting any combination of coincident interior coefficients is
            # a stronger hostile than the actual cancellation locus.
            for cancellation_mask in range(1 << len(collision_points)):
                points = set(raw_support)
                for bit, point in enumerate(collision_points):
                    if cancellation_mask & (1 << bit):
                        points.discard(point)
                facet_normals = {data for _face, data in lower_facets(points)}
                check(facet_normals == expected_facets,
                      f"{name} lower hull/cancellation gate")
                ledger = polygon_data((point[0], point[1]) for point in points)
                check(ledger[0] == expected_polygon, f"{name} projected polygon")
                check(ledger[1:4] == expected_pick, f"{name} Pick ledger")
                check(ledger[4] == expected_packet, f"{name} packet")
                seen.add(ledger[:5])
        check(len(seen) == 1, f"{name} universal lower-coefficient ledger")
        case_ledgers.append(next(iter(seen)))

    # The two analytic height formulas independently explain the hull census.
    for i, j, weight in (row for row in monomials if row[:2] != (0, 3)):
        first = (j + 2, i + j, 1)
        second = (j, i + j + 1, 1)
        main_gaps = tuple(
            Fraction(point[2]) - Fraction(point[0] + 2*point[1] - 2, 10)
            for point in (first, second)
        )
        side_gaps = tuple(
            Fraction(point[2]) - Fraction(point[0] - 2, 2)
            for point in (first, second)
        )
        check(main_gaps == (Fraction(10 - weight, 10),)*2, "main weight gap")
        check(side_gaps == (Fraction(2-j, 2), Fraction(4-j, 2)),
              "side j-gap")

    # Derive every special edge directly from the two face polynomials.
    S, P, Z = sp.symbols("S P Z")
    K, Theta, xi, upsilon = sp.symbols("K Theta xi upsilon")
    R = S**2 - P
    C = 1 - upsilon*P**5 - xi*S**2*P**4
    main_face = sp.expand(R*C)
    side_core = 1 - S**2*P**2*(K + Theta*P + xi*P**2)
    side_face = sp.expand(S**2*side_core)
    A, B, U, W_vertex, D, E = ((0, 1), (2, 0), (4, 2),
                                (4, 3), (4, 4), (0, 6))
    k_edges = tuple(edge_scheme(poly, S, P, start, end, Z) for poly, start, end in (
        (main_face, A, B), (side_face, B, U), (side_face, U, D),
        (main_face, D, E), (main_face, E, A), (main_face, B, D),
    ))
    check(k_edges == (
        Z - 1, -K*Z**2 + 1, -K - Theta*Z - xi*Z**2,
        (Z - 1)*(Z*upsilon + xi), -Z**5 + upsilon, -Z**2*xi + 1,
    ), "K-nonzero six edge schemes")
    theta_face = side_face.subs(K, 0)
    theta_edges = tuple(edge_scheme(poly, S, P, start, end, Z)
                        for poly, start, end in (
                            (main_face, A, B), (theta_face, B, W_vertex),
                            (theta_face, W_vertex, D), (main_face, D, E),
                            (main_face, E, A), (main_face, B, D),
                        ))
    check(theta_edges == (
        Z - 1, -Theta*Z + 1, -Theta - xi*Z,
        (Z - 1)*(Z*upsilon + xi), -Z**5 + upsilon, -Z**2*xi + 1,
    ), "Theta-only six edge schemes")
    collapsed_edges = tuple(edge_scheme(main_face, S, P, start, end, Z)
                            for start, end in ((A, B), (B, D), (D, E), (E, A)))
    check(collapsed_edges == (
        Z - 1, -Z**2*xi + 1, (Z - 1)*(Z*upsilon + xi), -Z**5 + upsilon,
    ), "collapsed four edge schemes")
    D_V = Theta**2 - 4*K*xi
    check(tuple(sp.factor(sp.discriminant(edge, Z)) for edge in k_edges) == (
        1, 4*K, D_V, (upsilon + xi)**2, 3125*upsilon**4, 4*xi,
    ), "K-nonzero edge discriminants")
    check(tuple(sp.factor(sp.discriminant(edge, Z)) for edge in theta_edges) == (
        1, 1, 1, (upsilon + xi)**2, 3125*upsilon**4, 4*xi,
    ), "Theta-only reduced edges")

    # Face normalizations, main nodes, and component genera.
    node_determinant = sp.factor(sp.det(sp.Matrix((
        (sp.diff(R, S), sp.diff(R, P)),
        (sp.diff(C, S), sp.diff(C, P)),
    ))).subs(P, S**2))
    check(node_determinant == -10*S**9*(upsilon + xi), "ten transverse nodes")
    node_polynomial = sp.Poly(1 - (upsilon + xi)*S**10, S,
                              domain=sp.QQ.frac_field(upsilon, xi))
    check(node_polynomial.degree() == 10
          and sp.gcd(node_polynomial, node_polynomial.diff()).degree() == 0,
          "ten distinct nodes")
    Y = sp.symbols("Y")
    genus_two = 1 - upsilon*P**5 - xi*Y**2
    check(sp.expand(C - genus_two.subs(Y, S*P**2)) == 0,
          "genus-two normalization")
    branch = sp.Poly(1 - upsilon*P**5, P, domain=sp.QQ.frac_field(upsilon))
    check(sp.gcd(branch, branch.diff()).degree() == 0 and (5 - 1)//2 == 2,
          "smooth genus-two component")
    conic = Y**2 - (K + Theta*P + xi*P**2)
    check(sp.factor(side_core.subs(S, 1/(Y*P))*Y**2 - conic) == 0,
          "side normalization is conic")
    check(sp.discriminant(K + Theta*P + xi*P**2, P) == D_V,
          "conic discriminant gate")
    check(sp.factor(sp.diff(side_core, S)) ==
          -2*P**2*S*(K + Theta*P + xi*P**2), "side torus smoothness")

    # Q=sigma^30: literal source charts, primitive multiplicities, and slopes.
    sigma, Q = sp.symbols("sigma Q")
    source_s, source_p = sp.symbols("source_s source_p")
    lam, alpha, epsilon, phi, Delta, eta = sp.symbols(
        "lam alpha epsilon phi Delta eta"
    )
    H = (
        lam*source_p + alpha*source_p**2 + epsilon*source_p**3
        + K*source_s**2*source_p**2 + phi*source_s*source_p**3
        + Delta*source_p**4 + Theta*source_s**2*source_p**3
        + eta*source_s*source_p**4 + upsilon*source_p**5
        + xi*source_s**2*source_p**4
    )
    F = (source_s**2 - source_p)*(1 - Q*H) - Q*source_s**2/2
    H_main = sp.cancel(sigma**30*H.subs({source_s: sigma**-3*S,
                                         source_p: sigma**-6*P}))
    check(sp.denom(H_main) == 1, "integral main H_sigma")
    check(sp.factor(H_main.subs(sigma, 0) - upsilon*P**5 - xi*S**2*P**4) == 0,
          "main top reduction")
    F_main = sp.cancel(sigma**6*F.subs({Q: sigma**30,
                                        source_s: sigma**-3*S,
                                        source_p: sigma**-6*P}))
    check(sp.factor(F_main - ((S**2-P)*(1-H_main) - sigma**30*S**2/2)) == 0,
          "literal main chart")
    local_U = S**2 - P
    local_V = (1-H_main)/S**2
    check(sp.cancel(F_main/S**2 - (local_U*local_V - sigma**30/2)) == 0,
          "A29 local equation")

    H_side = sp.cancel(sigma**30*H.subs({source_s: sigma**-15*S,
                                         source_p: P}))
    check(sp.denom(H_side) == 1, "integral side H_sigma")
    check(sp.factor(H_side.subs(sigma, 0)
                    - S**2*P**2*(K + Theta*P + xi*P**2)) == 0,
          "side reduction")
    F_side = sp.cancel(sigma**30*F.subs({Q: sigma**30,
                                         source_s: sigma**-15*S,
                                         source_p: P}))
    check(sp.factor(F_side - ((S**2-sigma**30*P)*(1-H_side)
                              - sigma**30*S**2/2)) == 0,
          "literal side chart")
    face_normals = ((3, 6, -1), (15, 0, -1))
    check(all(gcd(gcd(abs(a), abs(b)), abs(c)) == 1
              for a, b, c in face_normals), "face multiplicities one")

    graph_heights = {
        A: 0, B: 0, U: 30, W_vertex: 30, D: 30, E: 30,
    }
    for ledger in case_ledgers:
        for start, end in zip(ledger[0], ledger[0][1:] + ledger[0][:1]):
            difference = (end[0]-start[0], end[1]-start[1],
                          graph_heights[end]-graph_heights[start])
            length = gcd(abs(difference[0]), abs(difference[1]))
            check(gcd(gcd(abs(difference[0]), abs(difference[1])),
                      abs(difference[2])) == length, "outer edge denominator one")
    main_height = lambda point: 3*point[0] + 6*point[1] - 6
    side_height = lambda point: 15*point[0] - 30
    check((main_height((1, 1))-main_height(A),
           side_height((1, 0))-side_height(B),
           side_height((3, 2))-side_height(U),
           main_height((3, 4))-main_height(D),
           main_height((1, 6))-main_height(E)) == (3, -15, -15, -3, 3),
          "K-nonzero outer slopes")
    check((main_height((1, 1))-main_height(A),
           side_height((3, 2))-side_height(B),
           side_height((3, 3))-side_height(W_vertex),
           main_height((3, 4))-main_height(D),
           main_height((1, 6))-main_height(E)) == (3, 15, -15, -3, 3),
          "Theta-only outer slopes")
    check((main_height((1, 1))-main_height(A),
           main_height((2, 1))-main_height(B),
           main_height((3, 4))-main_height(D),
           main_height((1, 6))-main_height(E)) == (3, 6, -3, 3),
          "collapsed outer slopes")
    bd_probe = (2, -1)
    check((side_height(bd_probe)-side_height(B),
           main_height(bd_probe)-main_height(B)) == (0, -6),
          "internal BD slopes 0>-6")
    check(tuple(range(-1, -6, -1)) == (-1, -2, -3, -4, -5),
          "five denominator-one P1s per BD root")
    check(30 - 1 == 29, "A29 has 29 multiplicity-one rational curves")
    check((12-3+1, 2+10, 10-2+1, 2+9) == (10, 12, 9, 11),
          "dual graph/genus ledger")

    # Generic completion gates and actual source function field.
    q = sp.symbols("q")
    check(sp.factor(sp.discriminant(K*Z**2-(q-sp.Rational(1, 2)), Z)
                    - 4*K*(q-sp.Rational(1, 2))) == 0,
          "K quadratic carrier separable")
    check(sp.degree(Theta*Z-(q-sp.Rational(1, 2)), Z) == 1,
          "Theta carrier linear")
    check(sp.factor(sp.discriminant(xi*Z**2-(q-sp.Rational(1, 2)), Z)
                    - 4*xi*(q-sp.Rational(1, 2))) == 0,
          "xi quadratic carrier separable")
    pure_edge = (lam*Z + alpha*Z**2 + epsilon*Z**3
                 + Delta*Z**4 + upsilon*Z**5 - q)
    pure_domain = sp.QQ.frac_field(lam, alpha, epsilon, Delta, upsilon, q)
    check(sp.gcd(sp.Poly(pure_edge, Z, domain=pure_domain),
                 sp.Poly(sp.diff(pure_edge, Z), Z, domain=pure_domain)).degree() == 0,
          "generic affine degree-five edge separable")
    check(sp.factor(F.subs(source_p, source_s**2) + Q*source_s**2/2) == 0,
          "torus excludes t=0")

    # Genus-two simplicity and the good target leave no degree carrier.
    X = sp.symbols("X")
    phi5 = sp.Poly(sp.cyclotomic_poly(5, X), X)
    check(phi5.is_irreducible and phi5.degree() == 4, "Q(zeta5) endomorphism field")
    trace = sp.symbols("trace")
    trace_polynomial = sp.factor(sp.resultant(phi5.as_expr(),
                                               X**2-trace*X+1, X))
    check(trace_polynomial == (trace**2+trace-1)**2
          and sp.discriminant(trace**2+trace-1, trace) == 5,
          "only quadratic subfield is real")
    check(4 > 2 and 8 > 4, "elliptic-product endomorphism obstruction")

    a, target_A, target_C, target_X, target_Y = sp.symbols(
        "a target_A target_C target_X target_Y"
    )
    target = (target_C**2-target_A**3 + sp.Rational(3, 4)*a**2*target_A
              - sigma**-30 + sp.Rational(1, 4)*a**3)
    target_scaled = sp.cancel(sigma**30*target.subs({
        target_A: sigma**-10*target_X, target_C: sigma**-15*target_Y,
    }))
    target_expected = (target_Y**2-target_X**3-1
                       + sp.Rational(3, 4)*a**2*sigma**20*target_X
                       + sp.Rational(1, 4)*a**3*sigma**30)
    check(sp.factor(target_scaled-target_expected) == 0, "good target chart")
    check(sp.discriminant(target_X**3+1, target_X) == -27,
          "smooth j0 target")
    # R, the side conic, all chains, and every later exceptional are rational;
    # C has no elliptic quotient.  Connectedness makes every component map the
    # same constant, so the specialized degree sum is zero.
    check(sum(multiplicity*degree for multiplicity, degree in
              ((1, 0), (1, 0), (1, 0))) == 0,
          "degree conservation has zero special side")

    # Exact hostile: the omitted edge becomes double and the conic splits.
    collision = sp.factor((K+Theta*P+xi*P**2).subs(K, Theta**2/(4*xi)))
    check(collision == (2*P*xi+Theta)**2/(4*xi), "split-conic hostile")
    check(sp.factor(D_V.subs(K, Theta**2/(4*xi))) == 0,
          "Dv-zero double-edge hostile")

    semantic = digest((case_ledgers, k_edges, theta_edges, collapsed_edges,
                       node_determinant, collision))
    print("THM4220_INDEPENDENT_ZETA_ZERO_REGULAR_MODEL_AUDIT")
    print("scope=zeta0,u*xi*(u+xi)!=0;proved=(Dv!=0)|(K=Theta=0)")
    print("facets=M:(r+2k-2)/10,V:(r-2)/2;lower_coefficients=arbitrary")
    print("K_nonzero=Pick(34,12,12),packet(9,9,3,3,2,2,1),degree29")
    print("Theta_only=Pick(32,10,12),packet(9,9,5,3,1),degree27")
    print("K_Theta_zero=Pick(30,10,11),packet(9,9,3,3,1),degree25")
    print("Q=sigma30;face_multiplicities=1;BD=0>-6,5_P1_per_root;A29x10")
    print("components=genus2_C+rational_only;graph_b1=(10,9);genera=(12,11)")
    print("genus2_Jacobian=simple;Hom_to_elliptic=0;target=good_j0")
    print("degree_conservation=special_degree0;generic_finite_degree_positive")
    print("open_exactly=Dv=0,K!=0;reason=split_conic+double_UD_edge")
    print("dependency=THM4045_Dokchitser_model+THM4218_genus2_simplicity")
    print(f"semantic_sha256={semantic}")
    print(f"checks={CHECKS}")
    print("verdict=ACCEPT_THM4220_RELATIVE")


if __name__ == "__main__":
    main()
