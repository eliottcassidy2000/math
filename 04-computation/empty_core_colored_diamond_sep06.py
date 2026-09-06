#!/usr/bin/env python3
"""Exact controls for the proper-sublattice colored-diamond theorem.

The analytic proof is separate. This census does not prove unbounded scope.
Run: python3 -B 04-computation/empty_core_colored_diamond_sep06.py
"""
from collections import Counter
from functools import cmp_to_key
from math import gcd


def require(value, label):
    if not value:
        raise RuntimeError(label)


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def positive(v):
    return v[1] > 0 or v[1] == 0 and v[0] > 0


def direction(v):
    return v if positive(v) else (-v[0], -v[1])


def normalized_diamond(r, D):
    # u=(1,0), v=(r,D), 0<=r<D, gcd(r,D)=1.
    # Exact horizontal section: (r+1)y/D-1 <= x <= (r-1)y/D+1.
    points = [(1, 0)]
    for y in range(1, D + 1):
        lo = ((r + 1) * y - D + D - 1) // D
        hi = ((r - 1) * y + D) // D
        for x in range(lo, hi + 1):
            if gcd(x, y) == 1:
                points.append((x, y))
    return points


def brute_diamond(r, D):
    # Independent bounding-box scan and determinant membership test.
    u, v = (1, 0), (r, D)
    radius = max(1, r)
    return [(x, y) for x in range(-radius, radius + 1) for y in range(D + 1)
            if positive((x, y)) and gcd(x, y) == 1
            and abs(det((x, y), v)) + abs(det(u, (x, y))) <= D]


def angle_compare(u, v):
    hu, hv = 0 if positive(u) else 1, 0 if positive(v) else 1
    if hu != hv:
        return -1 if hu < hv else 1
    value = det(u, v)
    return -1 if value > 0 else 1 if value < 0 else 0


def audit_fan(points, invisible):
    signed = points + [(-x, -y) for x, y in points]
    fan = sorted(signed, key=cmp_to_key(angle_compare))
    for u, v in zip(fan, fan[1:] + fan[:1]):
        require(det(u, v) == 1, "consecutive primitive fan is unimodular")
        require(not (invisible(u) and invisible(v)), "no adjacent invisible rays")
    indices = [i for i, v in enumerate(fan) if not invisible(v)]
    live_directions = {direction(fan[i]) for i in indices}
    forced_segments = 0
    for offset, i in enumerate(indices):
        j = indices[(offset + 1) % len(indices)]
        gap = (j - i) % len(fan)
        require(gap in (1, 2), "at most one invisible separator")
        u, v = fan[i], fan[j]
        k = det(u, v)
        require(k >= 1, "consecutive live rays span less than pi")
        if gap == 1:
            require(k == 1, "direct live adjacency")
        else:
            w = fan[(i + 1) % len(fan)]
            require(v == (k * w[0] - u[0], k * w[1] - u[1]), "separator relation")
            forced = {(u[0] - t * w[0], u[1] - t * w[1]) for t in range(k + 1)}
            require(all(gcd(*p) == 1 and not invisible(p) for p in forced),
                    "segment points primitive and live")
            require(len({direction(p) for p in forced}) == k + 1,
                    "segment direction count")
            require(all(direction(p) in live_directions for p in forced),
                    "forced segment lies in original diamond")
            forced_segments += 1
    return forced_segments


def convex_hull(points):
    points = sorted(set(points))
    def turn(a, b, c):
        return det((b[0] - a[0], b[1] - a[1]), (c[0] - a[0], c[1] - a[1]))
    halves = []
    for ordered in (points, list(reversed(points))):
        half = []
        for p in ordered:
            while len(half) >= 2 and turn(half[-2], half[-1], p) <= 0:
                half.pop()
            half.append(p)
        halves.append(half[:-1])
    return halves[0] + halves[1]


def polygon_primitive_points(vertices, live):
    hull = convex_hull(vertices)
    edges = list(zip(hull, hull[1:] + hull[:1]))
    result = []
    for x in range(min(p[0] for p in hull), max(p[0] for p in hull) + 1):
        for y in range(min(p[1] for p in hull), max(p[1] for p in hull) + 1):
            if not positive((x, y)) or gcd(x, y) != 1 or not live((x, y)):
                continue
            if all(det((b[0] - a[0], b[1] - a[1]), (x - a[0], y - a[1])) >= 0
                   for a, b in edges):
                result.append((x, y))
    return sorted(result)


def circuit_type(points):
    require(len(points) == 3, "circuit takes exactly three directions")
    a, b, c = points
    values = [abs(det(b, c)), abs(det(c, a)), abs(det(a, b))]
    content = gcd(gcd(values[0], values[1]), values[2])
    return tuple(sorted(v // content for v in values))


def main():
    print("FINITE-EXACT colored-diamond controls; analytic theorem is separate")
    cases, segments = 0, 0
    minima, first = {}, {}
    types = Counter()
    for D in range(1, 301):
        minima[D] = 10**9
        for r in range(D):
            if gcd(r, D) != 1:
                continue
            points = normalized_diamond(r, D)
            if D <= 30:
                require(set(points) == set(brute_diamond(r, D)), "independent membership scan")
            for color in range(3):
                if (r + color * D) % 3 == 0:
                    continue
                invisible = lambda p, color=color: (p[0] + color * p[1]) % 3 == 0
                live = [p for p in points if not invisible(p)]
                cases += 1
                if D <= 40:
                    segments += audit_fan(points, invisible)
                if len(live) < minima[D]:
                    minima[D] = len(live)
                    first[D] = (r, color, live)
                require(D <= 2 or len(live) >= 4, "target colored-diamond consequence")
                if len(live) == 3:
                    kind = circuit_type(live)
                    require(kind in ((1, 1, 1), (1, 1, 2)), "three-direction circuit classification")
                    types[kind] += 1
    print("normalized universe 1<=D<=300, 0<=r<D, gcd(r,D)=1, f=x+c*y mod3, live endpoints")
    print("colored_diamonds", cases)
    print("independent membership replay D<=30; full primitive-fan audit D<=40")
    print("audited_invisible_separator_segments", segments)
    print("minima_D_1_to_12", [(D, minima[D]) for D in range(1, 13)])
    print("sharp_D3_example", first[3])
    print("exactly_three_diamond_circuit_counts", sorted(types.items()))

    composite_cases = 0
    for modulus in (2, 4, 5, 6, 7):
        for D in range(1, 26):
            for r in range(D):
                if gcd(r, D) != 1:
                    continue
                points = normalized_diamond(r, D)
                for a in range(1, modulus):
                    for b in range(modulus):
                        if gcd(gcd(a, b), modulus) != 1 or (a * r + b * D) % modulus == 0:
                            continue
                        invisible = lambda p, a=a, b=b, modulus=modulus: (a*p[0]+b*p[1])%modulus == 0
                        live = [p for p in points if not invisible(p)]
                        require(D <= 2 or len(live) >= 4, "other-index target")
                        audit_fan(points, invisible)
                        composite_cases += 1
    print("other_cyclic_indices 2,4,5,6,7; determinants<=25; cases", composite_cases)

    for expected, base in [((1, 1, 1), [(1, 0), (1, 1), (2, 1)]),
                           ((1, 1, 2), [(1, 0), (1, 1), (1, 2)])]:
        signed = base + [(-x, -y) for x, y in base]
        actual = polygon_primitive_points(signed, lambda p: p[0] % 3 != 0)
        require(set(actual) == set(base) and circuit_type(actual) == expected, "sharp circuit realization")
        print("sharp_symmetric_hull", actual, "primitive_circuit", expected)

    hostile = polygon_primitive_points([(0, 0), (1, 0), (1, 3)],
                                      lambda p: (p[0] + p[1]) % 3 != 0)
    require(hostile == [(1, 0), (1, 1), (1, 3)], "central symmetry hostile")
    require(circuit_type(hostile) == (1, 2, 3), "forbidden circuit without symmetry")
    print("hostile_without_central_symmetry", hostile, "primitive_circuit", circuit_type(hostile))
    print("PASS: complete diamond consequence, constructive segment proof, and scope controls")


if __name__ == "__main__":
    main()
