#!/usr/bin/env python3
"""All-parity norm-five proof head, with exact physical/network separation.

Complete parameter head: I=(x,y,2(x+y)), x<y, xy<84; and
II=sort(2x,y,x+y), xy<84. Retain distinct, primitive, ternary-unit speeds.
The full raw dictionary and literal six-sheet network are independent paths.
"""

from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import gcd
from pathlib import Path


SPEC = spec_from_file_location("boundary_referee", Path(__file__).with_name(
    "lrc14_sum_ray_boundary_discrepancy_overnight_hexagon_sep05.py"))
BOUNDARY = module_from_spec(SPEC)
SPEC.loader.exec_module(BOUNDARY)
BASE = BOUNDARY.BASE
CHECKS = 0
NETWORK = Q(46, 665)
PHYSICAL = Q(51, 770)


def need(value, payload):
    global CHECKS
    CHECKS += 1
    if not value:
        raise AssertionError(payload)


def family(kind, x, y):
    if kind == "I":
        w, d = (x, y, 2 * (x + y)), (2, 2, -1)
        ending = 2
    else:
        ordered = sorted(zip((2 * x, y, x + y), (1, 2, -2)))
        w = tuple(item[0] for item in ordered)
        d = tuple(item[1] for item in ordered)
        ending = w.index(x + y)
    return w, d, ending


def continuous_physical(w, d, B):
    q = Q(3, 7 * max(w))
    lines = [(q, Q(0))]
    for i in range(3):
        j, k = [s for s in range(3) if s != i]
        lines.append((Q(3 * (w[j] + w[k]), 14 * w[j] * w[k]),
                      Q(abs(d[i]), w[j] * w[k])))
    points = {Q(0), B}
    for (a, b), (c, e) in combinations(lines, 2):
        if b != e:
            t = (a - c) / (b - e)
            if 0 < t < B:
                points.add(t)
    points = sorted(points)
    f = lambda t: min(a - b * t for a, b in lines)
    area = sum(((right - left) * (f(left) + f(right)) / 2
                for left, right in zip(points, points[1:])), Q())
    return Q(4, 3) * area


def audit(kind, x, y, literal=True):
    w, d, ending = family(kind, x, y)
    h = x * y
    need(len(set(w)) == 3 and gcd(gcd(*w[:2]), w[2]) == 1 and all(t % 3 for t in w),
         ("typed family", kind, x, y, w))
    need(sum(a * b for a, b in zip(w, d)) == 0, ("norm-five relation", w, d))
    B = min(Q(3 * (sum(w) - w[i]), 14 * abs(d[i])) for i in range(3))
    K = (B.numerator - 1) // B.denominator
    expected = {tuple(sign * k * t for t in d)
                for k in range(1, K + 1) if k % 3 for sign in (-1, 1)}
    raw = BASE.carriers(w)
    need(raw == expected, ("complete live ray", w, raw ^ expected))
    projections, physical = BASE.projection_data(w, raw)
    if literal:
        need(BASE.literal_projection_data(w) == (projections, physical),
             ("literal network and overlap", w))
    q = Q(3, 7 * max(w))
    j, k = [s for s in range(3) if s != ending]
    slope = Q(abs(d[ending]), w[j] * w[k])
    cutoff = Q(3 * (w[j] + w[k]), 14 * abs(d[ending]))
    need(cutoff == B and slope == Q(1, h), ("zero-boundary projection", w, ending))
    plateau = B - q / slope
    need(0 <= plateau <= B, ("plateau", w))
    continuum = Q(4, 3) * (q * B - slope * (B - plateau)**2 / 2)
    need(continuum <= Q(3, 49), ("continuous selected ceiling", w, continuum))
    exact = continuum + 2 * slope * (BOUNDARY.primitive(B) - BOUNDARY.primitive(plateau))
    need(exact == projections[ending], ("zero-boundary quadrature", w))
    need(projections[ending] <= Q(3, 49) + Q(2, 3 * h), ("network product tail", w))
    physical_continuum = continuous_physical(w, d, B)
    need(physical_continuum == Q(3, 56), ("independent trapezoid/coarea identity", w, physical_continuum))
    maximum_slope = max(Q(abs(d[i]), w[(i + 1) % 3] * w[(i + 2) % 3]) for i in range(3))
    need(maximum_slope == Q(1, h), ("physical final slope", w))
    need(abs(physical - Q(3, 56)) <= Q(2, 3 * h), ("physical product error", w))
    return w, projections, physical


def main():
    need(Q(3, 49) + Q(2, 3 * 84) == Q(61, 882) < NETWORK < Q(6, 77),
         "network analytic product cutoff")
    need(Q(3, 56) + Q(2, 3 * 53) == Q(589, 8904) < PHYSICAL < NETWORK,
         "physical analytic product cutoff")
    rows = []
    counts = {"I": 0, "II": 0}
    digest = sha256()
    for x in range(1, 84):
        for y in range(1, 84 // x + 1):
            if x * y >= 84 or gcd(x, y) != 1:
                continue
            for kind in ("I", "II"):
                if kind == "I" and x >= y:
                    continue
                w, _, _ = family(kind, x, y)
                if len(set(w)) != 3 or not all(t % 3 for t in w):
                    continue
                w, projections, physical = audit(kind, x, y)
                need(min(projections) <= NETWORK and physical <= PHYSICAL,
                     ("complete finite head", w, projections, physical))
                record = (kind, x, y, w, projections, physical)
                rows.append(record)
                counts[kind] += 1
                digest.update((repr(record) + "\n").encode())
    need(counts == {"I": 41, "II": 82}, ("family universe", counts))
    need(len({row[3] for row in rows}) == 123, "distinct complete head triples")
    eq_network = [row[3] for row in rows if min(row[4]) == NETWORK]
    eq_physical = [row[3] for row in rows if row[5] == PHYSICAL]
    need(eq_network == [(2, 19, 20)], ("network equality", eq_network))
    need(eq_physical == [(1, 11, 20)], ("physical equality", eq_physical))
    # Uncolored geometric roofs genuinely allow a nonzero affine defect.
    w, d, C = (2, 19, 20), (1, 2, -2), (-8, 4, -3)
    need(sum(a * b for a, b in zip(w, C)) == 0, "uncolored integer kernel hostile")
    need(all(14 * abs(C[i]) < 3 * (sum(w) - w[i]) for i in range(3)),
         "uncolored strict roofs hostile")
    cross = (d[1] * C[2] - d[2] * C[1], d[2] * C[0] - d[0] * C[2],
             d[0] * C[1] - d[1] * C[0])
    need(cross == w and C[2] % 3 == 0, "owner gate deletes defect one")
    print("COMPLETE PRODUCT HEAD h<84", counts, "DISTINCT TRIPLES", len(rows))
    print("SEMANTIC SHA256", digest.hexdigest())
    print("NETWORK SHARP", NETWORK, "EQUALITY", eq_network,
          "TAIL AT84", Q(61, 882), "MARGIN", NETWORK - Q(61, 882))
    print("PHYSICAL SHARP", PHYSICAL, "EQUALITY", eq_physical,
          "TAIL AT53", Q(589, 8904), "MARGIN", PHYSICAL - Q(589, 8904))
    print("UNCOLORED HOSTILE", w, "relation", d, "carrier", C, "cross", cross)
    for kind, x, y in (("II", 1, 19), ("II", 10, 1), ("I", 1, 100),
                       ("I", 5, 17), ("II", 100, 1), ("II", 1, 100)):
        w, projections, physical = audit(kind, x, y)
        print("CONTROL", kind, (x, y), w, "E", projections, "PHYSICAL", physical)
    print("CHECKS", CHECKS, "INDEPENDENT RAW/LITERAL CHECKS", BASE.CHECKS)
    print("PROVED all-parity signed(1,2,2): min E<=46/665; actual mass<=51/770; distinct equality triples")
    print("OPEN other low circuits and actual arbitrary-body entry/synchronization")


if __name__ == "__main__":
    main()
