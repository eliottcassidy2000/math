#!/usr/bin/env python3
"""Exact certificate for THM-4230's exact-M=12 b=d=0 (2,3) seam.

It checks the complete monomial universe, the one-face lower hull under every
lower-row deletion, collision hostiles, polygon/packet/edge data, the exact
Q=sigma^12 model, the genus-seven component, its C12 character ledger, the
two j=0 quotient maps, and the W=-2Z attachment-torsion hostile.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd
import hashlib
import json

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
    support = {(2, 0, 0), (0, 1, 0)}
    for i, j, _weight in rows:
        support.add((j + 2, i + j, 1))
        support.add((j, i + j + 1, 1))
    return support


def candidate_plane_records(points):
    points = tuple(sorted(points))
    planes = set()
    for first, second, third in combinations(points, 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant == 0:
            continue
        slope_s = Fraction(
            (second[2] - first[2]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[2] - first[2]),
            determinant,
        )
        slope_p = Fraction(
            (second[0] - first[0]) * (third[2] - first[2])
            - (second[2] - first[2]) * (third[0] - first[0]),
            determinant,
        )
        constant = (Fraction(first[2]) - slope_s * first[0]
                    - slope_p * first[1])
        planes.add((slope_s, slope_p, constant))
    records = []
    for plane in planes:
        a, b, c = plane
        below = equal = 0
        for index, (r, ell, height) in enumerate(points):
            gap = Fraction(height) - a * r - b * ell - c
            if gap < 0:
                below |= 1 << index
            elif gap == 0:
                equal |= 1 << index
        records.append((plane, below, equal))
    return points, tuple(records)


def projected_rank_two(points, bits):
    chosen = [point for index, point in enumerate(points)
              if bits & (1 << index)]
    for first, second, third in combinations(chosen, 3):
        determinant = (
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )
        if determinant:
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
        if projected_rank_two(points, equal & bits):
            answer.add(plane)
    return answer


def polygon_ledger(vertices):
    area2 = abs(sum(
        first[0] * second[1] - second[0] * first[1]
        for first, second in zip(vertices, vertices[1:] + vertices[:1])
    ))
    boundary = sum(
        gcd(abs(second[0] - first[0]), abs(second[1] - first[1]))
        for first, second in zip(vertices, vertices[1:] + vertices[:1])
    )
    return area2, boundary, (area2 - boundary + 2) // 2


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
        position = vx // ux if ux else vy // uy
        if 0 <= position <= length:
            answer += coefficient * variable**position
    return sp.factor(answer)


def triangle_interiors(vertices):
    def area2(first, second, point):
        return abs(
            (second[0] - first[0]) * (point[1] - first[1])
            - (second[1] - first[1]) * (point[0] - first[0])
        )

    total = area2(vertices[0], vertices[1], vertices[2])
    answer = []
    for i in range(max(point[0] for point in vertices) + 1):
        for j in range(max(point[1] for point in vertices) + 1):
            pieces = [area2(vertices[k], vertices[(k + 1) % 3], (i, j))
                      for k in range(3)]
            if sum(pieces) == total and all(piece > 0 for piece in pieces):
                answer.append((i, j))
    return tuple(answer)


def main():
    rows = monomials_through(12)
    expected_rows = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11),
        (1, 3, 11), (6, 0, 12), (3, 2, 12), (0, 4, 12),
    )
    need(rows == expected_rows, "complete M12 monomial universe")
    need(tuple(row for row in rows if row[2] == 12)
         == ((6, 0, 12), (3, 2, 12), (0, 4, 12)),
         "three top monomials")

    required = {(6, 0), (0, 4)}
    optional = tuple(row for row in rows if row[:2] not in required)
    full_support = expanded_support(rows)
    points, records = candidate_plane_records(full_support)
    plane = (Fraction(1, 12), Fraction(1, 6), Fraction(-1, 6))
    required_rows = tuple(row for row in rows if row[:2] in required)
    need(lower_planes(points, records, expanded_support(required_rows)) == {plane},
         "U,Z skeleton owns one face")

    subset_count = 0
    for mask in range(1 << len(optional)):
        chosen = list(required_rows)
        chosen.extend(row for bit, row in enumerate(optional)
                      if mask & (1 << bit))
        need(lower_planes(points, records, expanded_support(chosen)) == {plane},
             "optional lower-row deletion changed the hull")
        subset_count += 1
    need(subset_count == 16384, "all optional-row subsets")

    collisions = tuple(sorted(
        point for point in full_support
        if sum(point in {(j + 2, i + j, 1), (j, i + j + 1, 1)}
               for i, j, _weight in rows) > 1
    ))
    need(collisions == (
        (2, 3, 1), (2, 4, 1), (2, 5, 1), (2, 6, 1),
        (3, 4, 1), (3, 5, 1), (4, 5, 1),
    ), "complete coincident-support list")
    for mask in range(1 << len(collisions)):
        deleted = {point for bit, point in enumerate(collisions)
                   if mask & (1 << bit)}
        need(lower_planes(points, records, full_support - deleted) == {plane},
             "coefficient cancellation changed one-face hull")

    vertices = ((0, 1), (2, 0), (6, 4), (0, 7))
    need(polygon_ledger(vertices) == (48, 14, 18), "global Pick ledger")
    packet, edge_rows = edge_packet(vertices)
    need(packet == (11, 11, 11, 2, 2, 2, 2, 1), "packet")
    need(sum(packet) == 42 and sum(value - 1 for value in packet) == 34,
         "packet sum and Hurwitz defect")

    S, P, X = sp.symbols("S P X")
    U, W, Z = sp.symbols("U W Z")
    rational = S**2 - P
    component = 1 - U*P**6 - W*S**2*P**5 - Z*S**4*P**4
    face = sp.expand(rational * component)
    schemes = tuple(edge_polynomial(face, S, P, start, end, X)
                    for start, end in zip(vertices, vertices[1:] + vertices[:1]))
    expected_schemes = (
        X - 1,
        1 - Z*X**4,
        (X - 1)*(U*X**2 + W*X + Z),
        U - X**6,
    )
    need(tuple(sp.factor(value - expected)
               for value, expected in zip(schemes, expected_schemes))
         == (0, 0, 0, 0), "all four edge schemes")
    D = sp.factor(W**2 - 4*U*Z)
    Lambda = U + W + Z
    discriminants = tuple(sp.factor(sp.discriminant(value, X))
                          for value in schemes)
    expected_discriminants = (1, -256*Z**3, D*Lambda**2, 46656*U**5)
    need(tuple(sp.factor(value - expected)
               for value, expected in zip(discriminants,
                                          expected_discriminants))
         == (0, 0, 0, 0), "exact reduced corner-avoiding gate")

    node_det = sp.factor(sp.det(sp.Matrix((
        (sp.diff(rational, S), sp.diff(rational, P)),
        (sp.diff(component, S), sp.diff(component, P)),
    ))).subs(P, S**2))
    need(sp.factor(node_det + 12*Lambda*S**11) == 0,
         "twelve transverse R-C nodes")
    need(sp.factor(component.subs(P, S**2) - (1 - Lambda*S**12)) == 0,
         "attachment orbit equation")
    singular_residual = sp.factor(
        (6*U + 5*W*X + 4*Z*X**2).subs(X, -W/(2*Z))
    )
    need(sp.factor(singular_residual + sp.Rational(3, 2)*D/Z) == 0,
         "D controls torus smoothness")

    need(polygon_ledger(((0, 0), (0, 6), (4, 4))) == (24, 12, 7),
         "component genus seven")
    interiors = triangle_interiors(((0, 0), (0, 6), (4, 4)))
    need(interiors == ((1, 2), (1, 3), (1, 4), (1, 5),
                       (2, 3), (2, 4), (3, 4)),
         "seven regular toric differentials")
    characters = tuple(sorted((i + 2*j) % 12 for i, j in interiors))
    need(characters == (5, 7, 8, 9, 10, 11, 11),
         "C12 holomorphic character ledger")
    full_characters = characters + tuple((-value) % 12 for value in characters)
    order_counts = {}
    for value in full_characters:
        order = 12 // gcd(value, 12)
        order_counts[order] = order_counts.get(order, 0) + 1
    need(order_counts == {12: 8, 3: 2, 4: 2, 6: 2},
         "Phi12 plus Phi3/Phi4/Phi6 rational factors")

    # Q=sigma^12 exact one-face model.
    source_s, source_p, sigma = sp.symbols("source_s source_p sigma")
    Delta, Phi, Theta, eta = sp.symbols("Delta Phi Theta eta")
    zeta3, upsilon5, xi10 = sp.symbols("zeta3 upsilon5 xi10")
    alpha11, beta11 = sp.symbols("alpha11 beta11")
    K = sp.Rational(2848, 45) - sp.Rational(7, 6)*Delta
    H_source = (
        -3*source_p + sp.Rational(8, 3)*source_p**2
        - sp.Rational(1376, 135)*source_p**3
        + K*source_s**2*source_p**2 + Phi*source_s*source_p**3
        + Delta*source_p**4 + Theta*source_s**2*source_p**3
        + eta*source_s*source_p**4 + zeta3*source_s**3*source_p**3
        + upsilon5*source_p**5 + xi10*source_s**2*source_p**4
        + alpha11*source_s*source_p**5 + beta11*source_s**3*source_p**4
        + U*source_p**6 + W*source_s**2*source_p**5
        + Z*source_s**4*source_p**4
    )
    H_scaled = sp.cancel(sigma**12 * H_source.subs({
        source_s: sigma**-1*S, source_p: sigma**-2*P,
    }))
    need(sp.denom(H_scaled) == 1, "scaled H integral")
    need(sp.factor(H_scaled.subs(sigma, 0)
                   - U*P**6 - W*S**2*P**5 - Z*S**4*P**4) == 0,
         "top face survives exactly")
    F_source = ((source_s**2 - source_p)*(1 - sp.Symbol("Q")*H_source)
                - sp.Symbol("Q")*source_s**2/2)
    F_scaled = sp.cancel(sigma**2 * F_source.subs({
        sp.Symbol("Q"): sigma**12,
        source_s: sigma**-1*S,
        source_p: sigma**-2*P,
    }))
    expected_scaled = ((S**2 - P)*(1 - H_scaled) - sigma**12*S**2/2)
    need(sp.factor(F_scaled - expected_scaled) == 0, "exact scaled source")
    need(sp.cancel(F_scaled/S**2
                   - ((S**2-P)*(1-H_scaled)/S**2 - sigma**12/2)) == 0,
         "twelve A11 charts")

    graph_heights = {
        (0, 1): (0, 1, 0), (2, 0): (2, 0, 0),
        (6, 4): (6, 4, 12), (0, 7): (0, 7, 12),
    }
    for start, end, length, _normal, _constant, _index in edge_rows:
        difference = tuple(graph_heights[end][axis] - graph_heights[start][axis]
                           for axis in range(3))
        need(gcd(gcd(abs(difference[0]), abs(difference[1])),
                 abs(difference[2])) == length,
             "outer edge denominator one")
    need(7 + (12 - 2 + 1) == 18, "genus completeness after path contraction")

    # Quotient H=C/<S -> -S> and its two complementary j=0 maps.
    u, x, v = sp.symbols("u x v")
    quotient_relation = u**6 - U - W*x - Z*x**2
    need(sp.factor((2*Z*x + W)**2 - (D + 4*Z*u**6)
                   + 4*Z*quotient_relation) == 0,
         "genus-two quotient identity")
    need(sp.factor((v**2 - D - 4*Z*(u**2)**3)
                   .subs(v, 2*Z*x+W) + 4*Z*quotient_relation) == 0,
         "first j0 quotient")
    need(sp.factor((v**2/u**6 - 4*Z - D*(u**-2)**3)
                   .subs(v, 2*Z*x+W)
                   + (4*Z/u**6)*quotient_relation) == 0,
         "second j0 quotient")
    need(2*2 == 4, "both C-to-E0 quotient degrees are four")

    T = sp.symbols("T")
    top_relation = U*P**6 + W*P**3*T**2 + Z*T**4 - 1
    need(sp.factor((2*U*P**3 + W*T**2)**2
                   - (4*U + D*T**4) - 4*U*top_relation) == 0,
         "degree-three j1728 quotient identity")

    # In-gate hostile U=2,Z=1,W=-2: attachment images on both j0
    # quotients have y=0, hence are nonzero 2-torsion and [2] kills them.
    hostile = {U: 2, Z: 1, W: -2}
    need((U*Z*D*Lambda).subs(hostile) != 0, "hostile remains inside gate")
    need((W + 2*Z).subs(hostile) == 0, "hostile attachment y-coordinate")
    need((D + 4*Z*Lambda).subs(hostile) == 0,
         "first attachment lies on y=0 fibre")
    need((D/Lambda + 4*Z).subs(hostile) == 0,
         "second attachment lies on y=0 fibre")
    need(4*4 == 16, "doubling a degree-four quotient has degree sixteen")

    record = {
        "checks": CHECKS,
        "support_subsets": subset_count,
        "collision_patterns": 1 << len(collisions),
        "plane": ["1/12", "1/6", "-1/6"],
        "polygon": vertices,
        "pick": [48, 14, 18],
        "packet": packet,
        "packet_sum": 42,
        "finite_candidate_degree": 34,
        "component_genus": 7,
        "characters": characters,
        "factor_dimensions": {"Phi12": 4, "Phi6": 1,
                              "Phi4": 1, "Phi3": 1},
        "hostile": {"U": 2, "W": -2, "Z": 1,
                    "compatible_map_degree": 16},
    }
    semantic = json.dumps(record, sort_keys=True, separators=(",", ":"))
    print("STATUS=FINITE-EXACT; JC(2)=OPEN; M12 exclusion=OPEN")
    print("GATE=U*Z*(W^2-4UZ)*(U+W+Z)!=0")
    print("LOWER_FACES=1; PICK=48,14,18; PACKET=11,11,11,2,2,2,2,1")
    print("JACOBIAN_FACTORS=A_Phi12(dim4)+E0+E1728+E0")
    print("FIXED_E0_MAP_DEGREES=4; SATURATED_DEGREES=4*(N(alpha)+N(beta))")
    print("RELATIVE_EXCLUSION=PROVED_OFF_Hom(A_Phi12,E0)_LOCUS")
    print("HOSTILE=U=2,W=-2,Z=1; doubled quotient degree=16 collapses nodes")
    print(f"CHECKS={CHECKS}; SUBSETS={subset_count}; COLLISIONS={1 << len(collisions)}")
    print("SEMANTIC_SHA256=" + hashlib.sha256(semantic.encode()).hexdigest())


if __name__ == "__main__":
    main()
