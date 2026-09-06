#!/usr/bin/env python3
"""Independent direct-box census and abstract colored-convex hostile probe.

The whole H=79 universe is checked using direct boxes, independent of the C++
congruence enumeration and integer-denominator accumulator. The literal-sheet
engine is explicitly imported as a third representation for named controls.
The abstract convex probe is a bounded experiment, not a classification proof.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import gcd
from pathlib import Path
import importlib.util

CHECKS = 0


def need(ok, detail):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(detail)


def primitive(v):
    g = gcd(*map(abs, v))
    return tuple(x//g*(1 if v[0] > 0 else -1) for x in v)


def det(u, v):
    return u[0]*v[1]-u[1]*v[0]


def carriers(w):
    bounds = [(3*(sum(w)-v)-1)//14 for v in w]
    points = set()
    for x in range(-bounds[0], bounds[0]+1):
        for y in range(-bounds[1], bounds[1]+1):
            if not x % 3 or not y % 3 or (w[0]*x+w[1]*y) % w[2]:
                continue
            z = -(w[0]*x+w[1]*y)//w[2]
            if z % 3 and abs(z) <= bounds[2]:
                points.add((x, y, z))
    return points


def projections(w, live):
    sums = [Q(0)]*3
    for p in live:
        for i in range(3):
            j, k = [t for t in range(3) if t != i]
            sums[i] += min(Q(3, 7*w[2]), Q(3*(w[j]+w[k])-14*abs(p[i]), 14*w[j]*w[k]))
    return tuple(sums)


def circuit_type(vs):
    vals = [abs(det(vs[(i+1) % 3], vs[(i+2) % 3])) for i in range(3)]
    common = gcd(*vals)
    return tuple(sorted(v//common for v in vals))


def has_additive_circuit(live):
    return any(primitive(u) != primitive(v) and tuple(x+y for x, y in zip(u, v)) in live
               for u in live for v in live)


def hull(points):
    def half(ps):
        h = []
        for v in ps:
            while len(h) > 1 and det(tuple(y-x for x, y in zip(h[-2], h[-1])),
                                    tuple(y-x for x, y in zip(h[-1], v))) <= 0:
                h.pop()
            h.append(v)
        return h[:-1]
    points = sorted(points)
    return half(points)+half(list(reversed(points)))


def inside(p, polygon):
    return all(det(tuple(y-x for x, y in zip(u, v)), tuple(y-x for x, y in zip(u, p))) >= 0
               for u, v in zip(polygon, polygon[1:]+polygon[:1]))


def main():
    classes, types = Counter(), Counter()
    leaders = [(Q(0), None)]*4
    misses, hostiles = [], []
    for w in combinations([v for v in range(1, 80, 2) if v % 3], 3):
        if gcd(*w) != 1:
            continue
        live = carriers(w)
        need(live == {tuple(-x for x in p) for p in live}, ('central symmetry', w))
        rays = sorted({primitive(p) for p in live})
        classes[len(rays)] += 1
        if len(rays) < 2:
            continue
        values = projections(w, live)
        need(max(values) < Q(6, 77), ('multi-ray every projection', w, values))
        stats = [max(values), min(values), sum(values)/3, Q(len(live), w[2])]
        for i, value in enumerate(stats):
            if leaders[i][0] < value:
                leaders[i] = value, w
        if 11*len(live) > 2*w[2]:
            misses.append(w)
        if len(rays) == 3:
            types[circuit_type(rays)] += 1
        if len(rays) >= 3 and not has_additive_circuit(live):
            hostiles.append(w)
    need(classes == {0: 36, 1: 1268, 2: 698, 3: 796, 4: 108, 5: 4}, ('C++ class match', classes))
    need(types == {(1, 1, 1): 505, (1, 1, 2): 291}, ('C++ circuit match', types))
    need(misses == [(19, 23, 29)], ('count-gate exceptions', misses))
    need(len(hostiles) == 266 and min(hostiles, key=lambda w: (w[2], w[0], w[1])) == (41, 47, 49), 'minimal additive hostile')
    expected = [(Q(18, 301), (5, 37, 43)), (Q(12, 343), (35, 47, 49)),
                (Q(1088, 24087), (11, 31, 37)), (Q(6, 29), (19, 23, 29))]
    need(leaders == expected, ('all C++ leader comparisons', leaders))
    print('INDEPENDENT H79 DIRECT BOX:', sum(classes.values()), 'rows; direction counts', sorted(classes.items()))
    print('MULTI LEADERS maxE,minE,meanE,N/c:', leaders)
    print('THREE-RAY CIRCUITS', sorted(types.items()), 'ADDITIVE HOSTILES', len(hostiles))

    # Independent literal interval path: separate source, explicit dependency.
    path = Path(__file__).with_name('lrc14_one_ray_overnight_hexagon_sep05.py')
    spec = importlib.util.spec_from_file_location('literal_sheet_dependency', path)
    literal = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(literal)
    controls = ((1, 5, 11), (1, 5, 7), (17, 23, 25), (19, 23, 29),
                (41, 47, 49), (5, 37, 43), (35, 47, 49), (11, 31, 37), (5, 49, 251))
    for w in controls:
        live = carriers(w)
        values = projections(w, live)
        need(values == literal.literal_projection_data(w)[0], ('literal control', w))
        print('LITERAL CONTROL', w, 'N', len(live), 'E', tuple(map(str, values)))
    hostile = carriers((41, 47, 49))
    base = {(11, 5, -14), (14, -7, -5), (17, -19, 4)}
    need(hostile == base | {tuple(-x for x in v) for v in base}, 'exact hostile full support')
    need(not has_additive_circuit(hostile), 'no independent additive circuit')
    need(tuple(x+y for x, y in zip((11, 5, -14), (17, -19, 4))) == tuple(2*x for x in (14, -7, -5)), 'retained AP circuit')

    # Every live point in the symmetric hull has a primitive representative
    # within the same coordinate box, so this finite primitive check is complete.
    B = 5
    vs = [(x, y) for x in range(1, B+1) if x % 3
          for y in range(-B, B+1) if gcd(x, abs(y)) == 1]
    abstract = Counter()
    for triple in combinations(vs, 3):
        polygon = hull(triple+tuple((-x, -y) for x, y in triple))
        if any(p not in triple and inside(p, polygon) for p in vs):
            continue
        abstract[circuit_type(triple)] += 1
    need(abstract == {(1, 1, 1): 10, (1, 1, 2): 19}, ('abstract box5', abstract))
    print('FINITE ABSTRACT BOX5: live x mod3!=0; symmetric closed hull; three primitive directions:', sorted(abstract.items()))
    print('OPEN: arbitrary height, universal colored circuit classification, LRC14')
    print('EXPLICIT CHECKS', CHECKS)


if __name__ == '__main__':
    main()
