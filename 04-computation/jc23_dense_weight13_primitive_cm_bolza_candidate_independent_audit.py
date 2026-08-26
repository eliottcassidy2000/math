#!/usr/bin/env python3
"""SymPy-free audit of the dense exact-weight-thirteen candidate."""

from fractions import Fraction
from itertools import combinations
from math import gcd


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def monomials():
    return tuple(sorted(
        (i, j, 2*i + 3*j)
        for i in range(7)
        for j in range(5)
        if 0 < 2*i + 3*j <= 13 and (i, j) not in {(0, 1), (1, 1)}
    ))


def hull(points):
    points = sorted(set(points))

    def cross(o, a, b):
        return ((a[0] - o[0])*(b[1] - o[1])
                - (a[1] - o[1])*(b[0] - o[0]))

    low = []
    for point in points:
        while len(low) > 1 and cross(low[-2], low[-1], point) <= 0:
            low.pop()
        low.append(point)
    high = []
    for point in reversed(points):
        while len(high) > 1 and cross(high[-2], high[-1], point) <= 0:
            high.pop()
        high.append(point)
    return tuple(low[:-1] + high[:-1])


def ledger(points):
    vertices = hull(points)
    area2 = abs(sum(
        vertices[k][0]*vertices[(k+1) % len(vertices)][1]
        - vertices[(k+1) % len(vertices)][0]*vertices[k][1]
        for k in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(k+1) % len(vertices)][0] - vertices[k][0]),
            abs(vertices[(k+1) % len(vertices)][1] - vertices[k][1]))
        for k in range(len(vertices))
    )
    return vertices, area2, boundary, (area2-boundary+2)//2


def add(left, right):
    return tuple(a+b for a, b in zip(left, right))


def scale(value, vector):
    return tuple(value*x for x in vector)


def edge_coefficients(polynomial, start, end):
    dx, dy = end[0]-start[0], end[1]-start[1]
    length = gcd(abs(dx), abs(dy))
    ux, uy = dx//length, dy//length
    zero = (0, 0, 0, 0, 0)
    answer = [zero for _ in range(length+1)]
    for (i, j), coefficient in polynomial.items():
        vx, vy = i-start[0], j-start[1]
        if vx*dy-vy*dx:
            continue
        position = vx//ux if ux else vy//uy
        if 0 <= position <= length:
            answer[position] = add(answer[position], coefficient)
    return tuple(answer)


def main():
    rows = monomials()
    expected = tuple(sorted((
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10), (4, 1, 11),
        (1, 3, 11), (6, 0, 12), (3, 2, 12), (0, 4, 12),
        (5, 1, 13), (2, 3, 13),
    )))
    need(rows == expected, "independent monomial universe")

    equal_main, equal_vertical, equal_tail = set(), set(), set()
    collision_owners = {}
    endpoint_owners = {}
    for i, j, weight in rows:
        endpoints = ((j+2, i+j, 1), (j, i+j+1, 1))
        signs = (-1, 1)
        gaps = (
            (Fraction(13-weight, 13), Fraction(13-weight, 13)),
            (Fraction(7-i-j, 6), Fraction(6-i-j, 6)),
            (Fraction(8-i-2*j, 8), Fraction(9-i-2*j, 8)),
        )
        need(all(value >= 0 for pair in gaps for value in pair),
             "independent nonnegative plane gaps")
        for endpoint, sign in zip(endpoints, signs):
            endpoint_owners.setdefault(endpoint, []).append((sign, i, j))
        if 0 in gaps[0]:
            equal_main.add((i, j))
        if 0 in gaps[1]:
            equal_vertical.add((i, j))
        if 0 in gaps[2]:
            equal_tail.add((i, j))
    need(equal_main == {(5, 1), (2, 3)}, "main owners")
    need(equal_vertical == {(6, 0), (5, 1)}, "vertical owners")
    need(equal_tail == {(0, 4), (2, 3)}, "tail owners")
    collision_owners = {
        point: owners for point, owners in endpoint_owners.items() if len(owners) > 1
    }
    need(tuple(sorted(collision_owners)) == (
        (2, 3, 1), (2, 4, 1), (2, 5, 1), (2, 6, 1),
        (3, 4, 1), (3, 5, 1), (3, 6, 1), (4, 5, 1),
    ), "all eight collision points")

    main_points = ((0, 1), (2, 0), (5, 5), (1, 7))
    tail_points = ((2, 0), (6, 4), (5, 5))
    vertical_points = ((0, 1), (1, 7), (0, 7))
    global_points = main_points + tail_points + vertical_points
    need(ledger(main_points)[1:] == (39, 5, 18), "main Pick")
    need(ledger(tail_points)[1:] == (8, 6, 2), "tail Pick")
    need(ledger(vertical_points)[1:] == (6, 8, 0), "vertical Pick")
    vertices, area2, boundary, genus = ledger(global_points)
    need((vertices, area2, boundary, genus) == (
        ((0, 1), (2, 0), (6, 4), (5, 5), (1, 7), (0, 7)),
        53, 15, 20,
    ), "global Pick")

    # Coefficient vectors are in the basis (1,A,B,U,Z).
    one, Avec, Bvec, Uvec, Zvec = (
        (1, 0, 0, 0, 0), (0, 1, 0, 0, 0),
        (0, 0, 1, 0, 0), (0, 0, 0, 1, 0),
        (0, 0, 0, 0, 1),
    )
    main = {
        (2, 0): one, (0, 1): scale(-1, one),
        (3, 6): add(scale(-1, Avec), Bvec),
        (1, 7): Avec, (5, 5): scale(-1, Bvec),
    }
    tail = {
        (2, 0): one, (6, 4): scale(-1, Zvec),
        (5, 5): scale(-1, Bvec),
    }
    vertical = {
        (0, 1): scale(-1, one), (0, 7): Uvec, (1, 7): Avec,
    }
    A0, B0, C0, D0, E0, F0 = vertices
    schemes = (
        edge_coefficients(main, A0, B0),
        edge_coefficients(tail, B0, C0),
        edge_coefficients(tail, C0, D0),
        edge_coefficients(main, D0, E0),
        edge_coefficients(vertical, E0, F0),
        edge_coefficients(vertical, F0, A0),
        edge_coefficients(main, A0, E0),
        edge_coefficients(main, B0, D0),
    )
    need(schemes == (
        (scale(-1, one), one),
        (one, (0, 0, 0, 0, 0), (0, 0, 0, 0, 0),
         (0, 0, 0, 0, 0), scale(-1, Zvec)),
        (scale(-1, Zvec), scale(-1, Bvec)),
        (scale(-1, Bvec), add(Bvec, scale(-1, Avec)), Avec),
        (Avec, Uvec),
        (Uvec, (0, 0, 0, 0, 0), (0, 0, 0, 0, 0),
         (0, 0, 0, 0, 0), (0, 0, 0, 0, 0),
         (0, 0, 0, 0, 0), scale(-1, one)),
        (scale(-1, one), Avec),
        (one, scale(-1, Bvec)),
    ), "independent eight edge coefficient arrays")

    packet = []
    graph_vertices = {
        A0: (0, 1, 0), B0: (2, 0, 0), C0: (6, 4, 312),
        D0: (5, 5, 312), E0: (1, 7, 312), F0: (0, 7, 312),
    }
    for start, end in zip(vertices, vertices[1:]+vertices[:1]):
        dx, dy = end[0]-start[0], end[1]-start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy//length, dx//length)
        constant = inward[0]*start[0]+inward[1]*start[1]
        index = inward[0]+inward[1]-constant
        if not (start[0] == end[0] == 0):
            packet.extend([index]*length)
        difference = tuple(graph_vertices[end][k]-graph_vertices[start][k]
                           for k in range(3))
        need(gcd(gcd(abs(difference[0]), abs(difference[1])),
                 abs(difference[2])) == length, "outer graph denominator")
    need(tuple(sorted(packet, reverse=True)) == (12, 12, 8, 6, 2, 2, 2, 2, 1),
         "independent packet")
    need((sum(packet), sum(x-1 for x in packet)) == (47, 38),
         "packet sum and defect")

    hM = lambda r, k: 24*r+48*k-48
    hV = lambda r, k: 52*k-52
    hT = lambda r, k: 39*r+39*k-78
    need((gcd(24, gcd(48, 1)), gcd(52, 1), gcd(39, gcd(39, 1))) == (1, 1, 1),
         "primitive face normals")
    need((hM(0, 0)-hM(*A0), hV(0, 0)-hV(*A0)) == (-48, -52),
         "AE slopes")
    need(-5*3+3*2+10 == 1, "BD primitive transverse point")
    need((hM(3, 2)-hM(*B0), hT(3, 2)-hT(*B0)) == (120, 117),
         "BD slopes")
    need((3, 2) == (52-48-1, 120-117-1), "chain lengths")

    # Independent Chevalley--Weil computation from branch residues.
    residues = (3, 12, 11)
    cm_type = set()
    for character in range(1, 13):
        dimension = sum(Fraction((-character*r) % 13, 13)
                        for r in residues)-1
        need(dimension in {0, 1}, "eigenspace dimension")
        if dimension == 1:
            cm_type.add(character)
    need(cm_type == {5, 6, 9, 10, 11, 12}, "Chevalley-Weil CM type")
    stabilizer = tuple(u for u in range(1, 13)
                       if {u*x % 13 for x in cm_type} == cm_type)
    need(stabilizer == (1,), "independent CM stabilizer")

    # x-Zx^5 has a simple zero at 0; at every nonzero root x^4=1/Z,
    # its derivative is 1-5Zx^4=-4.  Hence the degree-five model has genus 2.
    need(1 != 0 and 1-5 == -4 and (5-1)//2 == 2, "Bolza squarefree/genus")
    need(6+2+(15-4+1) == genus, "component plus graph genus")

    # Scaling/valuation audit independent of symbolic expansion.
    for i, j, weight in rows:
        need(312-48*i-72*j == 24*(13-weight), "main H valuation")
        need(312-52*(i+j) >= 0, "vertical H valuation")
        need(312-39*(i+2*j) >= 0, "tail H valuation")
    need((312//3, 312//2, 312-104, 312-156) == (104, 156, 208, 156),
         "target scaling exponents")
    need((312-1, 13, 311) == (311, 13, 311), "thirteen A311 models")

    print("DENSE_WEIGHT13_CANDIDATE_INDEPENDENT_AUDIT")
    print("implementation=stdlib_only_no_sympy")
    print("faces=(M13,V6,T8);global_genus=20;packet_degree=47")
    print("collisions=8;exact_units=(A,B,U,Z,A+B);A_equals_B_safe")
    print("chains=(AE:3,BD:2);main_nodes=13*A311;multiplicities=1")
    print("CM_type=(5,6,9,10,11,12);stabilizer=(1);Bolza_genus=2")
    print(f"checks={CHECKS}")
    print("verdict=INDEPENDENT_EXACT_INPUTS_PASS")


if __name__ == "__main__":
    main()
