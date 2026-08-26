#!/usr/bin/env python3
"""Exact certificate for THM-4232's exact-M=11 U=0 planar-Jacobian wall.

The certificate checks the complete inherited monomial universe, both
first-surviving pure-p rows, support/collision hostiles,
Newton/edge/regular-chart ledgers, and the unchanged primitive-CM component
from THM-4222.  Milne's primitive-CM simplicity theorem and proper-flat degree
conservation remain cited proof steps rather than computations performed here.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def monomials_through(weight):
    rows = []
    for i in range(weight // 2 + 1):
        for j in range(weight // 3 + 1):
            w = 2 * i + 3 * j
            if 0 < w <= weight and (i, j) not in {(0, 1), (1, 1)}:
                rows.append((i, j, w))
    return tuple(sorted(rows, key=lambda row: (row[2], row[1], row[0])))


def expanded_support(rows):
    # The fixed residual term -Q*s^2/2 is the third base point.  It is
    # strictly above every lower face here, but retaining it is essential for
    # a literally complete valued-support audit (the dense THM-4222 helper
    # omitted this harmless point).
    support = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j, _weight in rows:
        support.add((j + 2, i + j, 1))
        support.add((j, i + j + 1, 1))
    return support


def candidate_plane_records(points):
    points = tuple(sorted(points))
    records = []
    planes = set()
    for first, second, third in combinations(points, 3):
        det = ((second[0] - first[0]) * (third[1] - first[1])
               - (second[1] - first[1]) * (third[0] - first[0]))
        if det == 0:
            continue
        slope_s = Fraction(
            (second[2] - first[2]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[2] - first[2]), det)
        slope_p = Fraction(
            (second[0] - first[0]) * (third[2] - first[2])
            - (second[2] - first[2]) * (third[0] - first[0]), det)
        constant = (Fraction(first[2]) - slope_s * first[0]
                    - slope_p * first[1])
        planes.add((slope_s, slope_p, constant))
    for plane in planes:
        slope_s, slope_p, constant = plane
        below = 0
        equal = 0
        for index, (r, ell, height) in enumerate(points):
            gap = Fraction(height) - slope_s * r - slope_p * ell - constant
            if gap < 0:
                below |= 1 << index
            elif gap == 0:
                equal |= 1 << index
        records.append((plane, below, equal))
    return points, tuple(records)


def projected_rank_two(points, bits):
    chosen = [point for index, point in enumerate(points) if bits & (1 << index)]
    for first, second, third in combinations(chosen, 3):
        det = ((second[0] - first[0]) * (third[1] - first[1])
               - (second[1] - first[1]) * (third[0] - first[0]))
        if det:
            return True
    return False


def support_bits(points, support):
    index = {point: position for position, point in enumerate(points)}
    return sum(1 << index[point] for point in support)


def lower_planes(points, records, support):
    bits = support_bits(points, support)
    answer = set()
    for plane, below, equal in records:
        if below & bits:
            continue
        equality = equal & bits
        if projected_rank_two(points, equality):
            answer.add(plane)
    return answer


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, first, second):
        return ((first[0] - origin[0]) * (second[1] - origin[1])
                - (first[1] - origin[1]) * (second[0] - origin[0]))

    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(points):
    vertices = convex_hull(points)
    area2 = abs(sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))))
    boundary = sum(
        gcd(abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]))
        for index in range(len(vertices)))
    need((area2 - boundary + 2) % 2 == 0, "Pick parity")
    interior = (area2 - boundary + 2) // 2
    return vertices, area2, boundary, interior


def edge_packet(vertices):
    packet = []
    rows = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        rows.append((start, end, length, inward, constant, index))
        if start[0] == end[0] == 0:
            continue
        packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True)), tuple(rows)


def edge_polynomial(polynomial, first, second, start, end, variable):
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = gcd(abs(dx), abs(dy))
    ux, uy = dx // length, dy // length
    answer = 0
    for (i, j), coefficient in sp.Poly(polynomial, first, second).terms():
        vx, vy = i - start[0], j - start[1]
        if vx * dy - vy * dx:
            continue
        if ux:
            if vx % ux:
                continue
            position = vx // ux
        else:
            if vy % uy:
                continue
            position = vy // uy
        if 0 <= position <= length:
            answer += coefficient * variable**position
    return sp.factor(answer)


def triangle_interiors(vertices):
    def area2(first, second, point):
        return abs((second[0] - first[0]) * (point[1] - first[1])
                   - (second[1] - first[1]) * (point[0] - first[0]))

    total = area2(vertices[0], vertices[1], vertices[2])
    answer = []
    for i in range(max(point[0] for point in vertices) + 1):
        for j in range(max(point[1] for point in vertices) + 1):
            pieces = [area2(vertices[k], vertices[(k + 1) % 3], (i, j))
                      for k in range(3)]
            if sum(pieces) == total and all(piece > 0 for piece in pieces):
                answer.append((i, j))
    return tuple(answer)


def determinant_one(sequence):
    return all(abs(first[0] * second[1] - first[1] * second[0]) == 1
               for first, second in zip(sequence, sequence[1:]))


def unit_step_sequence(first, second):
    step = 1 if second >= first else -1
    return tuple((value, 1) for value in range(first, second + step, step))


def check_case(case, monomials, universal_points, plane_records, collisions):
    name = case["name"]
    deleted = case["deleted"]
    required = case["required"]
    expected_planes = case["planes"]
    retained = tuple(row for row in monomials if row[:2] not in deleted)
    required_rows = tuple(row for row in retained if row[:2] in required)
    optional_rows = tuple(row for row in retained if row[:2] not in required)

    hostile_count = 0
    for mask in range(1 << len(optional_rows)):
        chosen = list(required_rows)
        chosen.extend(row for bit, row in enumerate(optional_rows)
                      if mask & (1 << bit))
        support = expanded_support(chosen)
        for collision_mask in range(1 << len(collisions)):
            hostile_support = support - {
                point for bit, point in enumerate(collisions)
                if collision_mask & (1 << bit)
            }
            need(lower_planes(universal_points, plane_records, hostile_support)
                 == expected_planes,
                 f"{name}: lower hull changed under support/collision deletion")
            hostile_count += 1
    need(hostile_count == case["hostile_count"], f"{name}: hostile count")

    # Analytic gap identities independently explain the exhaustive result.
    main = case["main_plane"]
    tail = case["tail_plane"]
    replacement = case["replacement_plane"]
    for i, j, weight in retained:
        first = (j + 2, i + j, 1)
        second = (j, i + j + 1, 1)
        main_gaps = tuple(Fraction(point[2]) - main[0] * point[0]
                          - main[1] * point[1] - main[2]
                          for point in (first, second))
        tail_gaps = tuple(Fraction(point[2]) - tail[0] * point[0]
                          - tail[1] * point[1] - tail[2]
                          for point in (first, second))
        replacement_gaps = tuple(
            Fraction(point[2]) - replacement[0] * point[0]
            - replacement[1] * point[1] - replacement[2]
            for point in (first, second))
        need(main_gaps == (Fraction(11 - weight, 11),) * 2,
             f"{name}: main exact-weight gap")
        need(tail_gaps == (Fraction(3 - j, 3), Fraction(5 - j, 3)),
             f"{name}: tail exact gap")
        expected_replacement_gaps = (
            (Fraction(7 - i, 4), Fraction(4 - i, 4))
            if name == "Delta4"
            else (Fraction(8 - i + j, 3), Fraction(3 - i + j, 3))
        )
        need(replacement_gaps == expected_replacement_gaps
             and all(gap >= 0 for gap in replacement_gaps),
             f"{name}: replacement exact gap/nonnegativity")
    residual = (2, 0, 1)
    residual_gaps = tuple(
        Fraction(residual[2]) - plane[0] * residual[0]
        - plane[1] * residual[1] - plane[2]
        for plane in (main, tail, replacement))
    need(residual_gaps == case["residual_gaps"],
         f"{name}: fixed -Q*s^2/2 residual lies strictly above hull")

    # Polygon, face, and packet data.
    main_points = ((0, 1), (2, 0), (5, 4), (1, 6))
    tail_points = ((2, 0), (5, 3), (5, 4))
    replacement_points = case["replacement_points"]
    need(polygon_ledger(main_points)[1:] == (33, 5, 15),
         f"{name}: main Pick")
    need(polygon_ledger(tail_points)[1:] == (3, 5, 0),
         f"{name}: tail Pick")
    need(polygon_ledger(replacement_points)[1:] == case["replacement_ledger"],
         f"{name}: replacement Pick")
    polygon, area2, boundary, genus = polygon_ledger(
        main_points + tail_points + replacement_points)
    need((polygon, area2, boundary, genus) == case["global_ledger"],
         f"{name}: global Pick")
    packet, edge_rows = edge_packet(polygon)
    need(packet == (10, 10, 5, 4, 2, 2, 2, 1),
         f"{name}: packet")
    need(sum(packet) == 36 and sum(value - 1 for value in packet) == 28,
         f"{name}: packet degree/defect")

    # Face equations and every outer/internal edge scheme.
    S, P, X = sp.symbols("S P X")
    A, B, Z = sp.symbols("A B Z")
    c = case["coefficient"]
    k = case["pure_degree"]
    rational = S**2 - P
    cyclotomic = 1 - A*S*P**5 - B*S**3*P**4
    main_face = sp.expand(rational * cyclotomic)
    tail_core = 1 - Z*S**3*P**3 - B*S**3*P**4
    tail_face = sp.expand(S**2 * tail_core)
    replacement_core = -1 + c*P**k + A*S*P**5
    replacement_face = sp.expand(P * replacement_core)
    point_A, point_B, point_C, point_D, point_E, point_F = polygon
    schemes = (
        edge_polynomial(main_face, S, P, point_A, point_B, X),
        edge_polynomial(tail_face, S, P, point_B, point_C, X),
        edge_polynomial(tail_face, S, P, point_C, point_D, X),
        edge_polynomial(main_face, S, P, point_D, point_E, X),
        edge_polynomial(replacement_face, S, P, point_E, point_F, X),
        edge_polynomial(replacement_face, S, P, point_F, point_A, X),
        edge_polynomial(main_face, S, P, point_A, point_E, X),
        edge_polynomial(main_face, S, P, point_B, point_D, X),
    )
    expected_schemes = (
        X - 1,
        1 - Z*X**3,
        -Z - B*X,
        (X - 1)*(A*X + B),
        A + c*X,
        c - X**k,
        A*X - 1,
        1 - B*X,
    )
    need(tuple(sp.factor(actual - expected)
               for actual, expected in zip(schemes, expected_schemes))
         == (0,) * 8, f"{name}: eight edge schemes")
    discriminants = tuple(sp.factor(sp.discriminant(scheme, X))
                          for scheme in schemes)
    expected_pure_disc = -256*c**3 if k == 4 else -27*c**2
    need(discriminants == (1, -27*Z**2, 1, (A + B)**2,
                           1, expected_pure_disc, 1, 1),
         f"{name}: edge discriminants")

    # Face smoothness and unchanged primitive-CM main component.
    node_determinant = sp.factor(sp.det(sp.Matrix((
        (sp.diff(rational, S), sp.diff(rational, P)),
        (sp.diff(cyclotomic, S), sp.diff(cyclotomic, P)),
    ))).subs(P, S**2))
    need(sp.factor(node_determinant + 11*(A + B)*S**10) == 0,
         f"{name}: transverse main nodes")
    node_polynomial = sp.Poly(1 - (A + B)*S**11, S,
                              domain=sp.QQ.frac_field(A, B))
    need(node_polynomial.degree() == 11
         and sp.gcd(node_polynomial, node_polynomial.diff()).degree() == 0,
         f"{name}: eleven distinct nodes")
    need(sp.det(sp.Matrix(((A, 3*B), (5*A, 4*B)))) == -11*A*B,
         f"{name}: C torus smooth")
    need(sp.diff(replacement_core, S) == A*P**5,
         f"{name}: replacement rational/smooth")
    T0 = sp.symbols("T0")
    need(sp.diff(1 - T0**3*(Z + B*P), P) == -B*T0**3,
         f"{name}: tail rational/smooth")

    x = sp.symbols("x")
    x_expression = A*S*P**5
    need(sp.factor(A**3*P**11*(1 - x_expression) - B*x_expression**3
                   - A**3*P**11*cyclotomic) == 0,
         f"{name}: cyclic degree-eleven identity")
    interiors = triangle_interiors(((0, 0), (1, 5), (3, 4)))
    need(interiors == ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4)),
         f"{name}: CM interiors")
    cm_type = {(6*i + j) % 11 for i, j in interiors}
    need(cm_type == {4, 5, 8, 9, 10}, f"{name}: CM type")
    stabilizer = tuple(unit for unit in range(1, 11)
                       if {(unit*value) % 11 for value in cm_type} == cm_type)
    need(stabilizer == (1,), f"{name}: primitive CM stabilizer")
    need(sp.Poly(sp.cyclotomic_poly(11, x), x).degree() == 10,
         f"{name}: CM field degree 2g")
    branch_residues = (3, 10, 9)
    need(sum(branch_residues) % 11 == 0
         and all(gcd(11, value) == 1 for value in branch_residues),
         f"{name}: three primitive cyclic branch residues")
    need(2*5 - 2 == 11*(-2) + 3*(11 - 1),
         f"{name}: Riemann-Hurwitz genus five")
    chevalley_weil = set()
    for character in range(1, 11):
        dimension = -1 + sum(
            Fraction((-character*value) % 11, 11)
            for value in branch_residues)
        need(dimension in (0, 1), f"{name}: Chevalley-Weil multiplicity")
        if dimension == 1:
            chevalley_weil.add(character)
    need(chevalley_weil == cm_type,
         f"{name}: Chevalley-Weil/Newton CM agreement")
    complement = {(-value) % 11 for value in cm_type}
    need(cm_type.isdisjoint(complement)
         and cm_type | complement == set(range(1, 11)),
         f"{name}: complete nontrivial H1 character spectrum")

    # Integral lower graphs, primitive multiplicities, and toric chains.
    N = 132
    heights = case["height_functions"]
    normals = case["normals"]
    need(all(gcd(gcd(abs(a), abs(b)), abs(cn)) == 1
             for a, b, cn in normals), f"{name}: primitive normals")
    current_support = expanded_support(retained)
    need(len(current_support) == case["support_count"],
         f"{name}: complete active valued-support count")
    for r, ell, height in current_support:
        for height_function in heights:
            need(N*height - height_function((r, ell)) >= 0,
                 f"{name}: integral height gap")

    graph_vertices = {
        (0, 1): (0, 1, 0), (2, 0): (2, 0, 0),
        (5, 3): (5, 3, N), (5, 4): (5, 4, N),
        (1, 6): (1, 6, N), point_F: (point_F[0], point_F[1], N),
    }
    for start, end, length, _normal, _constant, _index in edge_rows:
        difference = tuple(graph_vertices[end][axis] - graph_vertices[start][axis]
                           for axis in range(3))
        need(gcd(gcd(abs(difference[0]), abs(difference[1])),
                 abs(difference[2])) == length,
             f"{name}: outer edge denominator")
    for start, end in ((point_A, point_E), (point_B, point_D)):
        difference = tuple(graph_vertices[end][axis] - graph_vertices[start][axis]
                           for axis in range(3))
        need(gcd(gcd(abs(difference[0]), abs(difference[1])),
                 abs(difference[2])) == 1,
             f"{name}: internal edge primitive")

    hm, ht, hv = heights
    outer_slopes = (
        hm((1, 1)) - hm(point_A),
        ht((1, 0)) - ht(point_B),
        ht((4, 3)) - ht(point_C),
        hm((4, 4)) - hm(point_D),
        hv((1, 5)) - hv(point_E),
        hv((1, point_F[1])) - hv(point_F),
    )
    need(outer_slopes == case["outer_slopes"], f"{name}: outer slopes")
    need(all(determinant_one(((slope, 1), (slope - 1, 1)))
             for slope in outer_slopes), f"{name}: outer determinant-one")
    ae_slopes = (hm((0, 0)) - hm(point_A), hv((0, 0)) - hv(point_A))
    bd_slopes = (hm((1, -1)) - hm(point_B), ht((1, -1)) - ht(point_B))
    ae_sequence = unit_step_sequence(*ae_slopes)
    bd_sequence = unit_step_sequence(*bd_slopes)
    need(ae_slopes == case["ae_slopes"] and determinant_one(ae_sequence)
         and len(ae_sequence) - 2 == case["ae_intermediate"],
         f"{name}: AE chain")
    need(bd_slopes == (-36, -44) and determinant_one(bd_sequence)
         and len(bd_sequence) - 2 == 7, f"{name}: BD chain")

    # Exact main, tail, and replacement charts.
    source_s, source_p, Q, sigma = sp.symbols("source_s source_p Q sigma")
    Delta, Phi, Theta, eta, Xi = sp.symbols("Delta Phi Theta eta Xi")
    epsilon = -sp.Rational(1376, 135)
    K = sp.Rational(2848, 45) - sp.Rational(7, 6)*Delta
    H = (
        -3*source_p + sp.Rational(8, 3)*source_p**2
        + epsilon*source_p**3 + K*source_s**2*source_p**2
        + Phi*source_s*source_p**3 + Delta*source_p**4
        + Theta*source_s**2*source_p**3 + eta*source_s*source_p**4
        + Z*source_s**3*source_p**3 + Xi*source_s**2*source_p**4
        + A*source_s*source_p**5 + B*source_s**3*source_p**4)
    if name == "Delta0":
        H = sp.expand(H.subs(Delta, 0))
    F_source = sp.expand((source_s**2 - source_p)*(1 - Q*H)
                         - Q*source_s**2/2)
    valued_coefficients = {
        (i, j, h): coefficient
        for (i, j, h), coefficient
        in sp.Poly(F_source, source_s, source_p, Q).terms()
    }
    need(valued_coefficients[(2, 0, 1)] == -sp.Rational(1, 2),
         f"{name}: fixed residual coefficient ledger")
    expected_collisions = case["collision_coefficients"]
    need(all(sp.factor(valued_coefficients[point] - coefficient) == 0
             for point, coefficient in expected_collisions.items()),
         f"{name}: complete aggregate collision coefficient ledger")
    need(set(valued_coefficients) == current_support,
         f"{name}: exact symbolic valued support equals hull universe")

    H_main = sp.cancel(sigma**N * H.subs({
        source_s: sigma**-12*S, source_p: sigma**-24*P}))
    need(sp.denom(H_main) == 1, f"{name}: main H integral")
    need(sp.factor(H_main.subs(sigma, 0)
                   - A*S*P**5 - B*S**3*P**4) == 0,
         f"{name}: main initial form")
    F_main = sp.cancel(sigma**24 * F_source.subs({
        Q: sigma**N, source_s: sigma**-12*S, source_p: sigma**-24*P}))
    expected_main = ((S**2 - P)*(1 - H_main) - sigma**N*S**2/2)
    need(sp.factor(F_main - expected_main) == 0, f"{name}: main chart")
    local_u = S**2 - P
    local_v = (1 - H_main)/S**2
    need(sp.cancel(F_main/S**2 - (local_u*local_v - sigma**N/2)) == 0,
         f"{name}: A131 local equation")
    need(N - 1 == 131, f"{name}: A131 resolution length")

    H_tail = sp.cancel(sigma**N * H.subs({
        source_s: sigma**-44*S, source_p: P}))
    need(sp.denom(H_tail) == 1, f"{name}: tail H integral")
    F_tail = sp.cancel(sigma**88 * F_source.subs({
        Q: sigma**N, source_s: sigma**-44*S, source_p: P}))
    expected_tail = ((S**2 - sigma**88*P)*(1 - H_tail)
                     - sigma**N*S**2/2)
    need(sp.factor(F_tail - expected_tail) == 0, f"{name}: tail chart")
    need(sp.factor(F_tail.subs(sigma, 0) - tail_face) == 0,
         f"{name}: tail special face")

    s_power, p_power, chart_factor = case["replacement_scaling"]
    H_replacement = sp.cancel(sigma**N * H.subs({
        source_s: sigma**s_power*S, source_p: sigma**p_power*P}))
    need(sp.denom(H_replacement) == 1, f"{name}: replacement H integral")
    need(sp.factor(H_replacement.subs(sigma, 0)
                   - (c*P**k + A*S*P**5)) == 0,
         f"{name}: replacement initial form")
    F_replacement = sp.cancel(sigma**chart_factor * F_source.subs({
        Q: sigma**N, source_s: sigma**s_power*S,
        source_p: sigma**p_power*P}))
    expected_replacement = (
        (sigma**(2*s_power - p_power)*S**2 - P)
        * (1 - H_replacement)
        - sigma**(N + 2*s_power - p_power)*S**2/2)
    need(sp.factor(F_replacement - expected_replacement) == 0,
         f"{name}: replacement chart")
    need(sp.factor(F_replacement.subs(sigma, 0) - replacement_face) == 0,
         f"{name}: replacement special face")

    # Generic moving edges and target good reduction.
    q = sp.symbols("q")
    carrier = Z*X**3 + K*X**2 - (q - sp.Rational(1, 2))
    if name == "Delta0":
        carrier = sp.expand(carrier.subs(Delta, 0))
    need(sp.Poly(carrier, X).degree() == 3, f"{name}: cubic carrier degree")
    carrier_field = sp.QQ.frac_field(Z, q) if name == "Delta0" \
        else sp.QQ.frac_field(Delta, Z, q)
    need(sp.gcd(sp.Poly(carrier, X, domain=carrier_field),
                sp.Poly(sp.diff(carrier, X), X, domain=carrier_field)).degree() == 0,
         f"{name}: cubic carrier separable")
    vertical_generic = (-3*X + sp.Rational(8, 3)*X**2
                        + epsilon*X**3)
    if name == "Delta4":
        vertical_generic += Delta*X**4
        vertical_field = sp.QQ.frac_field(Delta, q)
    else:
        vertical_field = sp.QQ.frac_field(q)
    vertical_generic -= q
    need(sp.Poly(vertical_generic, X).degree() == k,
         f"{name}: affine endpoint degree")
    need(sp.gcd(sp.Poly(vertical_generic, X, domain=vertical_field),
                sp.Poly(sp.diff(vertical_generic, X), X,
                        domain=vertical_field)).degree() == 0,
         f"{name}: affine endpoint separable")
    need(sp.factor(F_source.subs(source_p, source_s**2)
                   + Q*source_s**2/2) == 0,
         f"{name}: source inverse excludes t=0")
    need((sum(packet), sum(value - 1 for value in packet),
          sum(packet) - 3*2) == (36, 28, 30),
         f"{name}: full/defect/finite response ledger")

    target_a, target_A, target_C, target_X, target_Y = sp.symbols(
        "target_a target_A target_C target_X target_Y")
    target_equation = (target_C**2 - target_A**3
                       + sp.Rational(3, 4)*target_a**2*target_A
                       - sigma**-N + sp.Rational(1, 4)*target_a**3)
    scaled_target = sp.cancel(sigma**N * target_equation.subs({
        target_A: sigma**-44*target_X,
        target_C: sigma**-66*target_Y}))
    expected_target = (target_Y**2 - target_X**3 - 1
                       + sp.Rational(3, 4)*target_a**2*sigma**88*target_X
                       + sp.Rational(1, 4)*target_a**3*sigma**N)
    need(sp.factor(scaled_target - expected_target) == 0,
         f"{name}: target scaling")
    need(sp.discriminant(target_X**3 + 1, target_X) == -27,
         f"{name}: smooth j0 target")

    # Complete component/genus/degree inventory.
    core_vertices = 4
    core_edges = 11 + 1 + 1
    graph_rank = core_edges - core_vertices + 1
    need(graph_rank == 10 and 5 + graph_rank == genus,
         f"{name}: complete special genus")
    need((15, 0, 0) == (polygon_ledger(main_points)[3],
                         polygon_ledger(tail_points)[3],
                         polygon_ledger(replacement_points)[3]),
         f"{name}: face genus inventory")
    # Milne Prop. 3.13 + the exact primitive type above makes J(C) simple;
    # dim J(C)=5 then gives Hom(J(C),E0)=0.  Every other component is rational.
    component_degrees_after_cited_gate = (0, 0, 0, 0)
    need(sum(component_degrees_after_cited_gate) == 0,
         f"{name}: specialized degree zero")

    return {
        "name": name,
        "hostiles": hostile_count,
        "polygon": polygon,
        "pick": (area2, boundary, genus),
        "packet": packet,
        "outer_slopes": outer_slopes,
        "chains": (case["ae_intermediate"], 7),
        "pure_degree": k,
    }


def main():
    monomials = monomials_through(11)
    expected_monomials = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11),
        (1, 3, 11))
    need(monomials == expected_monomials, "complete exact-M11 universe")
    need(tuple(row for row in monomials if row[2] == 11)
         == ((4, 1, 11), (1, 3, 11)), "exactly two M11 top monomials")

    universal_support = expanded_support(monomials)
    universal_points, plane_records = candidate_plane_records(universal_support)
    collisions = ((2, 3, 1), (2, 4, 1), (2, 5, 1),
                  (3, 4, 1), (3, 5, 1))
    main_plane = (Fraction(1, 11), Fraction(2, 11), Fraction(-2, 11))
    tail_plane = (Fraction(1, 3), Fraction(0), Fraction(-2, 3))

    Delta = sp.symbols("Delta")
    epsilon = -sp.Rational(1376, 135)
    hm = lambda point: 12*point[0] + 24*point[1] - 24
    ht = lambda point: 44*point[0] - 88
    hv4 = lambda point: -33*point[0] + 33*point[1] - 33
    hv3 = lambda point: -88*point[0] + 44*point[1] - 44
    cases = (
        {
            "name": "Delta4",
            "deleted": {(5, 0)},
            "required": {(4, 0), (0, 3), (4, 1), (1, 3)},
            "hostile_count": 8192,
            "main_plane": main_plane,
            "tail_plane": tail_plane,
            "replacement_plane": (Fraction(-1, 4), Fraction(1, 4),
                                  Fraction(-1, 4)),
            "planes": {main_plane, tail_plane,
                       (Fraction(-1, 4), Fraction(1, 4), Fraction(-1, 4))},
            "replacement_points": ((0, 1), (1, 6), (0, 5)),
            "replacement_ledger": (4, 6, 0),
            "global_ledger": (((0, 1), (2, 0), (5, 3), (5, 4),
                               (1, 6), (0, 5)), 40, 12, 15),
            "coefficient": Delta,
            "pure_degree": 4,
            "height_functions": (hm, ht, hv4),
            "normals": ((12, 24, -1), (44, 0, -1), (-33, 33, -1)),
            "outer_slopes": (12, -44, -44, -12, -33, -33),
            "ae_slopes": (-24, -33),
            "ae_intermediate": 8,
            "replacement_scaling": (33, -33, 33),
            "support_count": 23,
            "residual_gaps": (Fraction(1), Fraction(1), Fraction(7, 4)),
            "collision_coefficients": {
                (2, 3, 1): (sp.Rational(2848, 45)
                            - sp.Rational(7, 6)*Delta - epsilon),
                (2, 4, 1): sp.symbols("Theta") - Delta,
                (2, 5, 1): sp.symbols("Xi"),
                (3, 4, 1): sp.symbols("Z") - sp.symbols("eta"),
                (3, 5, 1): sp.symbols("B") - sp.symbols("A"),
            },
        },
        {
            "name": "Delta0",
            "deleted": {(5, 0), (4, 0)},
            "required": {(3, 0), (0, 3), (4, 1), (1, 3)},
            "hostile_count": 4096,
            "main_plane": main_plane,
            "tail_plane": tail_plane,
            "replacement_plane": (Fraction(-2, 3), Fraction(1, 3),
                                  Fraction(-1, 3)),
            "planes": {main_plane, tail_plane,
                       (Fraction(-2, 3), Fraction(1, 3), Fraction(-1, 3))},
            "replacement_points": ((0, 1), (1, 6), (0, 4)),
            "replacement_ledger": (3, 5, 0),
            "global_ledger": (((0, 1), (2, 0), (5, 3), (5, 4),
                               (1, 6), (0, 4)), 39, 11, 15),
            "coefficient": epsilon,
            "pure_degree": 3,
            "height_functions": (hm, ht, hv3),
            "normals": ((12, 24, -1), (44, 0, -1), (-88, 44, -1)),
            "outer_slopes": (12, -44, -44, -12, -44, -88),
            "ae_slopes": (-24, -44),
            "ae_intermediate": 19,
            "replacement_scaling": (88, -44, 44),
            "support_count": 22,
            "residual_gaps": (Fraction(1), Fraction(1), Fraction(8, 3)),
            "collision_coefficients": {
                (2, 3, 1): sp.Rational(2848, 45) - epsilon,
                (2, 4, 1): sp.symbols("Theta"),
                (2, 5, 1): sp.symbols("Xi"),
                (3, 4, 1): sp.symbols("Z") - sp.symbols("eta"),
                (3, 5, 1): sp.symbols("B") - sp.symbols("A"),
            },
        },
    )

    results = [check_case(case, monomials, universal_points,
                          plane_records, collisions) for case in cases]
    need(sum(row["hostiles"] for row in results) == 12288,
         "combined hostile census")

    print("THM4232_WEIGHT11_U0_PRIMITIVE_CM_EXACT_CERTIFICATE")
    print("scope=exact_M11,b=d=0,reduced_(2,3),U=0,A*B*Z*(A+B)!=0,all_lower_coefficients")
    print("partition=Delta!=0|Delta=0;next_pure_owner=(p4|fixed_-1376/135*p3)")
    print("monomials=13;universal_valued_support=24;active_support=(Delta4:23,Delta0:22);fixed_residual_included")
    print("top=(p4y,py3);support_collision_hostiles=(256*32)+(128*32)=8192+4096=12288")
    print("faces_Delta4=M11,T3,V4;Pick=((33,5,15),(3,5,0),(4,6,0));global=(40,12,15)")
    print("faces_Delta0=M11,T3,V3;Pick=((33,5,15),(3,5,0),(3,5,0));global=(39,11,15)")
    print("packet_both=(10,10,5,4,2,2,2,1);full_degree=36;finite_degree=30;defect=28")
    print("edge_schemes_both=8_reduced;gate=A*B*Z*(A+B)*c_k!=0;c4=Delta;c3=-1376/135")
    print("base_change=Q:sigma^132;face_multiplicities=1;main_nodes=11*A131")
    print("chains_Delta4=AE:8,BD:7;chains_Delta0=AE:19,BD:7;all_inserted_components=rational")
    print("positive_genus=C11:g5;dual_graph_rank=10;special_genus=15")
    print("C11=P^11-(B/A^3)x^3/(1-x);CM_type=(4,5,8,9,10);stabilizer=(1)")
    print("CITED_INPUT=Milne_CM_Prop_3.13;primitive_CM_implies_simple;Hom(JC11,E0)=0")
    print("target_special=Y^2-X^3-1;degree_conservation=special_sum_0")
    print(f"checks={CHECKS}")
    print("verdict=EXACT_M11_U0_WALL_EXCLUDED_RELATIVE")


if __name__ == "__main__":
    main()
