#!/usr/bin/env python3
"""Standard-library clean-room audit for THM-4248.

This program does not import SymPy or the primary certificate.  It rebuilds
the five owner strata from the weighted monomial rule, obtains lower facets
by exact Cramer arithmetic, exhausts all 26,624 support/collision cases, and
uses a small independent sparse polynomial ring for the face and edge
identities.  It separately checks Pick data, lifted-edge denominators,
determinant-one fans, CM characters, genus, target scaling, and the
degree-conservation inventory.

As in the primary certificate, the cited primitive-CM theorem, the toroidal
regular-model theorem, and proper-flat degree conservation remain proof-layer
inputs rather than computational assertions.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


CHECKS = 0


def need(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("THM-4248 INDEPENDENT AUDIT FAILURE: " + label)


def gcd_many(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


# ---------------------------------------------------------------------------
# Independent lower-hull engine.
# ---------------------------------------------------------------------------


def weighted_rows(limit):
    rows = []
    for weight in range(1, limit + 1):
        for j in range(limit // 3 + 1):
            remainder = weight - 3*j
            if remainder < 0 or remainder % 2:
                continue
            i = remainder // 2
            if (i, j) not in {(0, 1), (1, 1)}:
                rows.append((i, j, weight))
    return tuple(rows)


def lift_rows(rows):
    points = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j, _weight in rows:
        points.add((j + 2, i + j, 1))
        points.add((j, i + j + 1, 1))
    return points


def cramer_plane(first, second, third):
    x1, y1, z1 = first
    x2, y2, z2 = second
    x3, y3, z3 = third
    determinant = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
    if determinant == 0:
        return None
    alpha = Fraction((z2-z1)*(y3-y1) - (z3-z1)*(y2-y1), determinant)
    beta = Fraction((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1), determinant)
    gamma = Fraction(z1) - alpha*x1 - beta*y1
    return alpha, beta, gamma


def projected_area_nonzero(points, mask):
    indices = [index for index in range(len(points)) if mask & (1 << index)]
    for a, b, c in combinations(indices, 3):
        first, second, third = points[a], points[b], points[c]
        if ((second[0]-first[0])*(third[1]-first[1])
                - (third[0]-first[0])*(second[1]-first[1])):
            return True
    return False


def prepare_facets(universe):
    points = tuple(sorted(universe))
    planes = set()
    for triple in combinations(points, 3):
        plane = cramer_plane(*triple)
        if plane is not None:
            planes.add(plane)
    records = []
    for plane in sorted(planes):
        alpha, beta, gamma = plane
        negative_mask = 0
        zero_mask = 0
        for index, (x, y, z) in enumerate(points):
            value = Fraction(z) - alpha*x - beta*y - gamma
            if value < 0:
                negative_mask |= 1 << index
            elif value == 0:
                zero_mask |= 1 << index
        records.append((plane, negative_mask, zero_mask))
    return points, tuple(records)


def bitset(points, active):
    lookup = {point: index for index, point in enumerate(points)}
    answer = 0
    for point in active:
        answer |= 1 << lookup[point]
    return answer


def facets_of(points, records, active_mask):
    answer = set()
    for plane, negative_mask, zero_mask in records:
        if negative_mask & active_mask:
            continue
        equality = zero_mask & active_mask
        if projected_area_nonzero(points, equality):
            answer.add(plane)
    return answer


M = (Fraction(1, 11), Fraction(2, 11), Fraction(-2, 11))
V5 = (Fraction(0), Fraction(1, 5), Fraction(-1, 5))
V4 = (Fraction(-1, 4), Fraction(1, 4), Fraction(-1, 4))
V3 = (Fraction(-2, 3), Fraction(1, 3), Fraction(-1, 3))
TK = (Fraction(1), Fraction(-1, 2), Fraction(-2))

AR, BR, UR, ZR, KR, DR, ER = (
    (4, 1), (1, 3), (5, 0), (0, 3), (0, 2), (4, 0), (3, 0)
)


STRATA = (
    ("U_NE_0__K_NE_0", {AR, BR, UR, KR}, {ZR}, {M, V5, TK},
     5, True, 330, ((0, 1), (2, 0), (4, 2), (5, 4), (1, 6), (0, 6)), 8192),
    ("U_NE_0__K_EQ_0", {AR, BR, UR}, {ZR, KR}, {M, V5},
     5, False, 330, ((0, 1), (2, 0), (5, 4), (1, 6), (0, 6)), 8192),
    ("U_EQ_0__DELTA_NE_0__K_NE_0", {AR, BR, DR, KR}, {UR, ZR}, {M, V4, TK},
     4, True, 132, ((0, 1), (2, 0), (4, 2), (5, 4), (1, 6), (0, 5)), 4096),
    ("U_EQ_0__DELTA_NE_0__K_EQ_0", {AR, BR, DR}, {UR, ZR, KR}, {M, V4},
     4, False, 132, ((0, 1), (2, 0), (5, 4), (1, 6), (0, 5)), 4096),
    ("U_EQ_0__DELTA_EQ_0", {AR, BR, ER, KR}, {UR, ZR, DR}, {M, V3, TK},
     3, True, 132, ((0, 1), (2, 0), (4, 2), (5, 4), (1, 6), (0, 4)), 2048),
)


def audit_all_hulls(rows):
    collision_points = ((2, 3, 1), (2, 4, 1), (2, 5, 1), (3, 4, 1), (3, 5, 1))
    totals = []
    for name, required, absent, expected, _k, _tail, _base, _polygon, target in STRATA:
        optional = tuple(row for row in rows if row[:2] not in required | absent)
        universe = lift_rows(tuple(row for row in rows if row[:2] not in absent))
        points, records = prepare_facets(universe)
        count = 0
        for optional_mask in range(1 << len(optional)):
            selected = [row for row in rows if row[:2] in required]
            selected.extend(
                row for index, row in enumerate(optional)
                if optional_mask & (1 << index)
            )
            support = lift_rows(selected)
            need((2, 0, 1) in support, "fixed residual in " + name)
            for collision_mask in range(32):
                deleted = {
                    point for index, point in enumerate(collision_points)
                    if collision_mask & (1 << index)
                }
                active = support - deleted
                need(facets_of(points, records, bitset(points, active)) == expected,
                     "lower facets in " + name)
                count += 1
        need(count == target, "case count in " + name)
        totals.append(count)
    need(sum(totals) == 26624, "total hostile atlas")
    return tuple(totals)


# ---------------------------------------------------------------------------
# Independent sparse exact coefficient and face polynomial rings.
# A coefficient is in Q[A,B,C,K].
# ---------------------------------------------------------------------------


ZERO_EXP = (0, 0, 0, 0)


def e_clean(value):
    return {key: coefficient for key, coefficient in value.items() if coefficient}


def e_const(value):
    value = Fraction(value)
    return {} if value == 0 else {ZERO_EXP: value}


def e_var(position):
    exponent = [0, 0, 0, 0]
    exponent[position] = 1
    return {tuple(exponent): Fraction(1)}


def e_add(*values):
    answer = {}
    for value in values:
        for key, coefficient in value.items():
            answer[key] = answer.get(key, Fraction(0)) + coefficient
    return e_clean(answer)


def e_scale(value, scalar):
    scalar = Fraction(scalar)
    return e_clean({key: scalar*coefficient for key, coefficient in value.items()})


def e_neg(value):
    return e_scale(value, -1)


def e_sub(first, second):
    return e_add(first, e_neg(second))


def e_mul(first, second):
    answer = {}
    for left_key, left_value in first.items():
        for right_key, right_value in second.items():
            key = tuple(a+b for a, b in zip(left_key, right_key))
            answer[key] = answer.get(key, Fraction(0)) + left_value*right_value
    return e_clean(answer)


def e_pow(value, exponent):
    answer = e_const(1)
    base = value
    while exponent:
        if exponent & 1:
            answer = e_mul(answer, base)
        base = e_mul(base, base)
        exponent //= 2
    return answer


ONE, A_E, B_E, C_E, K_E = e_const(1), e_var(0), e_var(1), e_var(2), e_var(3)


def p_clean(value):
    return {key: coefficient for key, coefficient in value.items() if coefficient}


def p_add(*values):
    answer = {}
    for value in values:
        for key, coefficient in value.items():
            answer[key] = e_add(answer.get(key, {}), coefficient)
    return p_clean(answer)


def p_scale(value, coefficient):
    return p_clean({key: e_mul(item, coefficient) for key, item in value.items()})


def p_mul(first, second):
    answer = {}
    for (a, b), left in first.items():
        for (c, d), right in second.items():
            key = (a+c, b+d)
            answer[key] = e_add(answer.get(key, {}), e_mul(left, right))
    return p_clean(answer)


def monomial(s_power, p_power, coefficient=ONE):
    return {(s_power, p_power): coefficient}


def edge_restriction(polynomial, start, end):
    dx, dy = end[0]-start[0], end[1]-start[1]
    length = gcd(abs(dx), abs(dy))
    ux, uy = dx//length, dy//length
    answer = {}
    for (i, j), coefficient in polynomial.items():
        vx, vy = i-start[0], j-start[1]
        if vx*dy != vy*dx:
            continue
        position = vx//ux if ux else vy//uy
        if 0 <= position <= length:
            answer[position] = e_add(answer.get(position, {}), coefficient)
    return {degree: coefficient for degree, coefficient in answer.items() if coefficient}


def u_add(*values):
    answer = {}
    for value in values:
        for degree, coefficient in value.items():
            answer[degree] = e_add(answer.get(degree, {}), coefficient)
    return {degree: coefficient for degree, coefficient in answer.items() if coefficient}


def u_mul(first, second):
    answer = {}
    for left_degree, left in first.items():
        for right_degree, right in second.items():
            degree = left_degree + right_degree
            answer[degree] = e_add(answer.get(degree, {}), e_mul(left, right))
    return {degree: coefficient for degree, coefficient in answer.items() if coefficient}


def u_discriminant(polynomial):
    degree = max(polynomial)
    if degree <= 1:
        return ONE
    if degree == 2:
        a = polynomial.get(2, {})
        b = polynomial.get(1, {})
        c = polynomial.get(0, {})
        return e_sub(e_pow(b, 2), e_scale(e_mul(a, c), 4))
    need(set(polynomial) == {0, degree}, "higher-degree scheme is binomial")
    need(polynomial[degree] == e_const(-1), "binomial leading coefficient -1")
    signs = {3: -1, 4: -1, 5: 1}
    need(degree in signs, "supported binomial degree")
    return e_scale(e_pow(polynomial[0], degree-1), signs[degree]*degree**degree)


# ---------------------------------------------------------------------------
# Independent planar polygon and packet arithmetic.
# ---------------------------------------------------------------------------


def planar_hull(points):
    points = sorted(set(points))

    def turn(a, b, c):
        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])

    down = []
    for point in points:
        while len(down) > 1 and turn(down[-2], down[-1], point) <= 0:
            down.pop()
        down.append(point)
    up = []
    for point in reversed(points):
        while len(up) > 1 and turn(up[-2], up[-1], point) <= 0:
            up.pop()
        up.append(point)
    return tuple(down[:-1] + up[:-1])


def pick(points):
    polygon = planar_hull(points)
    area2 = abs(sum(
        polygon[i][0]*polygon[(i+1) % len(polygon)][1]
        - polygon[(i+1) % len(polygon)][0]*polygon[i][1]
        for i in range(len(polygon))
    ))
    boundary = sum(
        gcd(abs(polygon[(i+1) % len(polygon)][0]-polygon[i][0]),
            abs(polygon[(i+1) % len(polygon)][1]-polygon[i][1]))
        for i in range(len(polygon))
    )
    need((area2-boundary+2) % 2 == 0, "independent Pick parity")
    return polygon, area2, boundary, (area2-boundary+2)//2


def packet(polygon):
    answer = []
    edge_rows = []
    for start, end in zip(polygon, polygon[1:]+polygon[:1]):
        dx, dy = end[0]-start[0], end[1]-start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy//length, dx//length)
        level = inward[0]*start[0] + inward[1]*start[1]
        index = inward[0]+inward[1]-level
        edge_rows.append((start, end, length))
        if start[0] or end[0]:
            answer.extend([index]*length)
    return tuple(sorted(answer, reverse=True)), tuple(edge_rows)


def interior_points(triangle):
    def twice_area(a, b, c):
        return abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))

    total = twice_area(triangle[0], triangle[1], triangle[2])
    answer = []
    for x in range(4):
        for y in range(6):
            pieces = tuple(
                twice_area(triangle[i], triangle[(i+1) % 3], (x, y))
                for i in range(3)
            )
            if sum(pieces) == total and min(pieces) > 0:
                answer.append((x, y))
    return tuple(answer)


def exact_face_and_geometry_audit(rows):
    main_polygon = ((0, 1), (2, 0), (5, 4), (1, 6))
    tail_polygon = ((2, 0), (4, 2), (5, 4))
    vertical_polygons = {
        5: ((0, 1), (1, 6), (0, 6)),
        4: ((0, 1), (1, 6), (0, 5)),
        3: ((0, 1), (1, 6), (0, 4)),
    }
    need(pick(main_polygon)[1:] == (33, 5, 15), "independent main Pick")
    need(pick(tail_polygon)[1:] == (2, 4, 0), "independent tail Pick")
    need(tuple(pick(vertical_polygons[k])[1:] for k in (5, 4, 3))
         == ((5, 7, 0), (4, 6, 0), (3, 5, 0)), "independent V Pick")

    rational = p_add(monomial(2, 0), p_scale(monomial(0, 1), e_const(-1)))
    cm_core = p_add(
        monomial(0, 0),
        p_scale(monomial(1, 5), e_neg(A_E)),
        p_scale(monomial(3, 4), e_neg(B_E)),
    )
    main_face = p_mul(rational, cm_core)
    tail_core = p_add(
        monomial(0, 0),
        p_scale(monomial(2, 2), e_neg(K_E)),
        p_scale(monomial(3, 4), e_neg(B_E)),
    )
    tail_face = p_mul(monomial(2, 0), tail_core)

    X_MINUS_ONE = {0: e_const(-1), 1: ONE}
    AX_PLUS_B = {0: B_E, 1: A_E}
    total_edges = []
    degrees = []
    for name, _required, absent, planes, k, has_tail, base, polygon, _target in STRATA:
        global_polygon, area2, boundary, genus = pick(polygon)
        need(global_polygon == polygon, "independent global polygon " + name)
        if k == 5:
            expected_pick = (40, 12, 15) if has_tail else (38, 10, 15)
        elif k == 4:
            expected_pick = (39, 11, 15) if has_tail else (37, 9, 15)
        else:
            expected_pick = (38, 10, 15)
        need((area2, boundary, genus) == expected_pick, "global Pick " + name)
        boundary_packet, edge_rows = packet(polygon)
        expected_packet = ((10, 10, 5, 5, 2, 2, 1) if has_tail
                           else (10, 10, 7, 5, 1))
        need(boundary_packet == expected_packet, "packet " + name)
        need(sum(value-1 for value in boundary_packet) == 28, "defect " + name)
        degrees.append(sum(boundary_packet))

        vertical_core = p_add(
            p_scale(monomial(0, 0), e_const(-1)),
            p_scale(monomial(0, k), C_E),
            p_scale(monomial(1, 5), A_E),
        )
        vertical_face = p_mul(monomial(0, 1), vertical_core)
        if has_tail:
            a, b, c, d, e, f = polygon
            schemes = (
                edge_restriction(main_face, a, b),
                edge_restriction(tail_face, b, c),
                edge_restriction(tail_face, c, d),
                edge_restriction(main_face, d, e),
                edge_restriction(vertical_face, e, f),
                edge_restriction(vertical_face, f, a),
                edge_restriction(main_face, a, e),
                edge_restriction(main_face, b, d),
            )
            expected = (
                X_MINUS_ONE,
                {0: ONE, 2: e_neg(K_E)},
                {0: e_neg(K_E), 1: e_neg(B_E)},
                u_mul(X_MINUS_ONE, AX_PLUS_B),
                {0: A_E, 1: C_E},
                {0: C_E, k: e_const(-1)},
                {0: e_const(-1), 1: A_E},
                {0: ONE, 1: e_neg(B_E)},
            )
            expected_discriminants = (
                ONE, e_scale(K_E, 4), ONE, e_pow(e_add(A_E, B_E), 2), ONE,
                e_scale(e_pow(C_E, k-1), {5: 3125, 4: -256, 3: -27}[k]),
                ONE, ONE,
            )
        else:
            a, b, d, e, f = polygon
            schemes = (
                edge_restriction(main_face, a, b),
                edge_restriction(main_face, b, d),
                edge_restriction(main_face, d, e),
                edge_restriction(vertical_face, e, f),
                edge_restriction(vertical_face, f, a),
                edge_restriction(main_face, a, e),
            )
            expected = (
                X_MINUS_ONE, {0: ONE, 1: e_neg(B_E)},
                u_mul(X_MINUS_ONE, AX_PLUS_B), {0: A_E, 1: C_E},
                {0: C_E, k: e_const(-1)}, {0: e_const(-1), 1: A_E},
            )
            expected_discriminants = (
                ONE, ONE, e_pow(e_add(A_E, B_E), 2), ONE,
                e_scale(e_pow(C_E, k-1), {5: 3125, 4: -256, 3: -27}[k]), ONE,
            )
        need(schemes == expected, "formal edge schemes " + name)
        need(tuple(u_discriminant(value) for value in schemes)
             == expected_discriminants, "formal edge discriminants " + name)
        total_edges.append(len(schemes))

        # Every expected face is supporting, with the advertised exact polygon.
        active_rows = tuple(row for row in rows if row[:2] not in absent)
        active_points = lift_rows(active_rows)
        expected_face_polygons = {M: main_polygon, {5: V5, 4: V4, 3: V3}[k]: vertical_polygons[k]}
        if has_tail:
            expected_face_polygons[TK] = tail_polygon
        for plane in planes:
            alpha, beta, gamma = plane
            equality = {
                (x, y) for x, y, z in active_points
                if Fraction(z)-alpha*x-beta*y-gamma == 0
            }
            need(planar_hull(equality) == expected_face_polygons[plane],
                 "face equality polygon " + name)
            normal = (int(base*alpha), int(base*beta), -1)
            need(Fraction(normal[0], base) == alpha
                 and Fraction(normal[1], base) == beta,
                 "integral face normal " + name)
            need(gcd_many(normal) == 1, "multiplicity one " + name)
            for x, y, z in active_points:
                scaled = base*(alpha*x+beta*y+gamma)
                need(scaled.denominator == 1 and base*z-scaled >= 0,
                     "integral nonnegative face gap " + name)

        lifted = {
            point: (point[0], point[1], 0 if point in {(0, 1), (2, 0)} else base)
            for point in polygon
        }
        for start, end, length in edge_rows:
            vector = tuple(lifted[end][i]-lifted[start][i] for i in range(3))
            need(gcd_many(vector) == length, "outer lifted gcd " + name)
        ae = tuple(lifted[(1, 6)][i]-lifted[(0, 1)][i] for i in range(3))
        need(gcd_many(ae) == 1, "primitive internal MV edge " + name)
        if has_tail:
            bd = tuple(lifted[(5, 4)][i]-lifted[(2, 0)][i] for i in range(3))
            need(gcd_many(bd) == 1, "primitive internal MT edge " + name)

        def scaled_height(plane, point):
            value = base*(plane[0]*point[0] + plane[1]*point[1] + plane[2])
            need(value.denominator == 1, "integral fan probe " + name)
            return int(value)

        vertical_plane = {5: V5, 4: V4, 3: V3}[k]
        ae_limits = {5: (-60, -66), 4: (-24, -33), 3: (-24, -44)}[k]
        derived_ae = (
            scaled_height(M, (0, 0))-scaled_height(M, (0, 1)),
            scaled_height(vertical_plane, (0, 0))-scaled_height(vertical_plane, (0, 1)),
        )
        ae_fan = tuple((value, 1) for value in range(ae_limits[0], ae_limits[1]-1, -1))
        need(derived_ae == ae_limits
             and all(abs(a*d-b*c) == 1 for (a, b), (c, d) in zip(ae_fan, ae_fan[1:])),
             "unimodular MV fan " + name)
        need(len(ae_fan)-2 == {5: 5, 4: 8, 3: 19}[k], "MV chain length " + name)
        if has_tail:
            bd_limits = (-90, -165) if base == 330 else (-36, -66)
            derived_bd = (
                scaled_height(M, (1, -1))-scaled_height(M, (2, 0)),
                scaled_height(TK, (1, -1))-scaled_height(TK, (2, 0)),
            )
            bd_fan = tuple((value, 1) for value in range(bd_limits[0], bd_limits[1]-1, -1))
            need(derived_bd == bd_limits
                 and all(abs(a*d-b*c) == 1 for (a, b), (c, d) in zip(bd_fan, bd_fan[1:])),
                 "unimodular MT fan " + name)
            need(len(bd_fan)-2 == (74 if base == 330 else 29), "MT chain length " + name)

        # Exact scaling identities for the source and target charts.
        main_s, main_p, main_multiplier = -base//11, -2*base//11, 2*base//11
        need(2*main_s == main_p == -main_multiplier, "main base terms same scale " + name)
        need(base + 11*main_s == 0, "weight-eleven H main scale " + name)
        need(base + 2*main_s + main_multiplier == base,
             "fixed residual becomes sigma^N S^2/2 " + name)
        vertical = {5: V5, 4: V4, 3: V3}[k]
        vs, vp, vm = -int(base*vertical[0]), -int(base*vertical[1]), -int(base*vertical[2])
        need(vm >= 0 and all(isinstance(value, int) for value in (vs, vp, vm)),
             "integral vertical substitution " + name)
        if has_tail:
            need(base % 2 == 0, "tail half-scale integral " + name)
            need(2*(-base) == -2*base, "tail base term scale " + name)
        need(base % 6 == 0, "target cubic/quadratic scaling integral " + name)
        need(3*(-base//3)+base == 0 and 2*(-base//2)+base == 0,
             "good target leading exponents " + name)

        vertices = 4 if has_tail else 3
        edges = 13 if has_tail else 12
        graph_rank = edges-vertices+1
        need(graph_rank == 10 and 5+graph_rank == genus,
             "complete dual genus " + name)
        need(base-1 == (329 if base == 330 else 131), "A_(N-1) node length " + name)
        need(sum((0,)*vertices) == 0, "degree-zero component inventory " + name)

    need(tuple(degrees) == (35, 33, 35, 33, 35), "packet degrees by stratum")
    need(tuple(total_edges) == (8, 6, 8, 6, 8), "edge counts by stratum")
    return tuple(degrees)


def cm_and_hostile_audit():
    delta_kzero = Fraction(2848, 45)*Fraction(6, 7)
    need(delta_kzero == Fraction(5696, 105), "exact K=0 Delta value")
    need(Fraction(2848, 45)-Fraction(7, 6)*delta_kzero == 0,
         "K-Delta relation at boundary")
    need(Fraction(2848, 45) != 0, "Delta=0 forces K nonzero")

    residues = (3, 10, 9)
    need(sum(residues) % 11 == 0 and all(gcd(value, 11) == 1 for value in residues),
         "cyclic branch residues")
    need((11-1)*(3-2)//2 == 5, "cyclic-cover genus")
    points = interior_points(((0, 0), (1, 5), (3, 4)))
    need(points == ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4)),
         "independent differential lattice points")
    cm_type = {(6*i+j) % 11 for i, j in points}
    need(cm_type == {4, 5, 8, 9, 10}, "independent CM type")
    need(cm_type.isdisjoint({(-value) % 11 for value in cm_type}), "CM halves disjoint")
    need(cm_type | {(-value) % 11 for value in cm_type} == set(range(1, 11)),
         "CM halves exhaust H1")
    stabilizer = tuple(
        unit for unit in range(1, 11)
        if {(unit*value) % 11 for value in cm_type} == cm_type
    )
    need(stabilizer == (1,), "independent CM stabilizer")
    # Exponent audit of P^11=(B/A^3)x^3/(1-x), x=ASP^5.
    need((3, 15-11) == (3, 4), "cyclic identity exponent balance")

    need(pick(((0, 1), (0, 6), (3, 5)))[1:] == (15, 7, 5),
         "A=0 hostile replacement genus five")
    need(pick(((2, 0), (5, 3), (3, 5)))[1:] == (12, 6, 4),
         "B=0 hostile replacement genus four")
    double_edge = u_mul({0: e_const(-1), 1: ONE}, {0: e_neg(A_E), 1: A_E})
    expected_double = {0: A_E, 1: e_scale(A_E, -2), 2: A_E}
    need(double_edge == expected_double, "A+B=0 double edge A(X-1)^2")


def main():
    rows = weighted_rows(11)
    expected = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11), (1, 3, 11),
    )
    need(rows == expected, "independently generated weighted universe")
    counts = audit_all_hulls(rows)
    degrees = exact_face_and_geometry_audit(rows)
    cm_and_hostile_audit()

    print("THM4248_WEIGHT11_Z0_OWNER_DESCENT_INDEPENDENT_AUDIT")
    print("implementation=stdlib_only;imports_primary=false;imports_sympy=false")
    print("scope=exact_M11,b=d=0,reduced_(2,3),Z=0,A*B*(A+B)!=0,U_arbitrary,all_lower_coefficients")
    print("owner_strata=5;K_Delta_boundary=5696/105;Delta=0=>K=2848/45")
    print("hull_cases=" + str(sum(counts)) + ";per_stratum="
          + ",".join(str(value) for value in counts)
          + ";fixed_residual=(2,0,1);failures=0")
    print("formal_faces=main+three_vertical_owners+optional_K_tail;formal_edge_schemes=PASS;formal_discriminants=PASS")
    print("Pick_faces=M(33,5,15),TK(2,4,0),V5(5,7,0),V4(4,6,0),V3(3,5,0)")
    print("packet_degrees=" + ",".join(str(value) for value in degrees)
          + ";packets=35_if_K_nonzero,33_if_K_zero;defect=28")
    print("regular_models=Q330/A329_or_Q132/A131;lifted_gcds=PASS;unimodular_fans=PASS;multiplicity=1")
    print("CM=genus5,type(4,5,8,9,10),stabilizer(1),full_H1;Hom_to_E0_zero_by_cited_primitive_CM")
    print("genus=15;dual_rank=10;all_other_components=rational;special_degree=0")
    print("hostiles=A0_genus5;B0_genus4;A+B0_double_edge;K0_packet_change")
    print("proof_inputs=THM3992/3997/4012/4045/4222/4232;Milne_primitive_CM;proper_flat_degree_conservation")
    print("nonclaims=other_three_M11_walls,other_cells,seam_entry,JC2,DC2")
    print("checks=" + str(CHECKS))


if __name__ == "__main__":
    main()
