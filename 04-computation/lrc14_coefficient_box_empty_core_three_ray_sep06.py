#!/usr/bin/env python3
"""Complete exact width-slope coefficient box, maximum coefficient<=18.

Independent rational polygon clipping handles both actual-zero coefficients
and full support. Every signed coefficient permutation is retained. The frozen
full-support producer's cube-edge compiler provides a declared second path for
full-support patterns; a closed rectangle formula checks support-two patterns.
No enumeration over speed heights is used in this artifact.
"""

import importlib.util
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement, permutations, product
from math import gcd
from pathlib import Path

CHECKS = 0
R = Q(3, 14)
TARGET = Q(15, 98)
EXCLUDED = {(0, 1, 1), (1, 1, 2)}


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def keep_halfplane(polygon, A, B, C):
    """Sutherland-Hodgman clipping of an exact rational convex polygon."""
    if not polygon:
        return []
    result = []
    for u, v in zip(polygon, polygon[1:]+polygon[:1]):
        fu, fv = A*u[0]+B*u[1]-C, A*v[0]+B*v[1]-C
        if fu <= 0:
            result.append(u)
        if fu*fv < 0:
            t = fu/(fu-fv)
            result.append((u[0]+t*(v[0]-u[0]), u[1]+t*(v[1]-u[1])))
    return list(dict.fromkeys(result))


def plane_polygon(v, delta, lo, hi):
    need(v[0] != 0, ('nonzero elimination pivot', v))
    polygon = [(lo, lo), (hi, lo), (hi, hi), (lo, hi)]
    lower, upper = sorted((v[0]*lo, v[0]*hi))
    polygon = keep_halfplane(polygon, -v[1], -v[2], upper-delta)
    polygon = keep_halfplane(polygon, v[1], v[2], delta-lower)
    vertices = tuple(((delta-v[1]*y-v[2]*z)/v[0], y, z) for y, z in polygon)
    need(bool(vertices), ('nonempty closed plane slice', v, delta, lo, hi))
    for e in vertices:
        need(all(lo <= x <= hi for x in e) and sum(v[i]*e[i] for i in range(3)) == delta,
             ('exact clipped vertex', v, delta, e))
    return vertices


def sectors(pattern):
    return tuple(sorted({v for a in permutations(pattern) for k in range(3)
                         if a[k] != 0
                         for v in [tuple(-a[i] if i == k else a[i] for i in range(3))]
                         if any(x > 0 for x in v) and any(x < 0 for x in v)}))


def compile_by_clipping(pattern):
    unit = all(x % 3 for x in pattern)
    bound = (3*sum(pattern)-1)//14
    defects = tuple(d for d in range(-bound, bound+1) if (d % 3 == 0) == unit)
    rho = Q(2, 3) if unit else Q(1, 3)
    intercept = (Q(4, 3) if unit else Q(1))*len(defects)
    best = Q(-1)
    witness = None
    for original in sectors(pattern):
        pivot = next(i for i in range(3) if original[i])
        v = original[pivot:]+original[:pivot]
        # A second nonzero coordinate exists for both permitted support sizes.
        i = next(i for i in range(1, 3) if v[i])
        j, k = (i+1) % 3, (i+2) % 3
        error_polygons = [plane_polygon(v, Q(d), -R, R) for d in defects]
        speed_polygon = plane_polygon(v, Q(0), Q(0), Q(1))
        for w in speed_polygon:
            value = Q(0)
            for polygon in error_polygons:
                scalars = [(w[j]*e[k]-w[k]*e[j])/v[i] for e in polygon]
                value += rho*(max(scalars)-min(scalars))
            if value > best:
                best, witness = value, (v, w)
    return defects, best, intercept, witness


def support_two_formula(pattern, defects):
    _, p, q = pattern
    return (2*R*p*len(defects)+sum(min(2*p*R, R*(p+q)-abs(d)) for d in defects))/(3*p*q)


def main():
    frozen_path = Path(__file__).with_name('lrc14_empty_core_certificate_sep06.py')
    spec = importlib.util.spec_from_file_location('frozen_cube_edge_compiler', frozen_path)
    other = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(other)
    patterns = tuple(p for p in combinations_with_replacement(range(19), 3)
                     if sum(x != 0 for x in p) >= 2 and sum(p) % 2 == 0
                     and gcd(gcd(p[0], p[1]), p[2]) == 1
                     and sum(x % 3 == 0 for x in p) <= 1 and p not in EXCLUDED)
    second_universe = {tuple(sorted(p)) for p in product(range(19), repeat=3)
                       if sum(x != 0 for x in p) >= 2 and sum(p) % 2 == 0
                       and gcd(gcd(p[0], p[1]), p[2]) == 1
                       and sum(x % 3 == 0 for x in p) <= 1 and tuple(sorted(p)) not in EXCLUDED}
    need(set(patterns) == second_universe, 'independent coefficient universe construction')
    need(len(patterns) == 308, ('complete box cardinality', len(patterns)))
    need(Counter(sum(x != 0 for x in p) for p in patterns) == {2: 15, 3: 293}, 'support counts')
    # The forbidden equal-speed coefficient row is geometrically impossible
    # in the distinct positive-speed domain. The norm-four row is a real
    # hostile to this count-slope target and has its own proved network gate.
    norm4 = compile_by_clipping((1, 1, 2))
    need(norm4[1] > TARGET, ('norm-four slope hostile', norm4))
    equalities = []
    maximum = Q(0)
    digest = sha256()
    for pattern in patterns:
        defects, slope, intercept, witness = compile_by_clipping(pattern)
        if pattern[0]:
            alternate = other.compile_pattern(pattern)
            need((defects, slope, intercept) == (alternate['deltas'], alternate['slope'], alternate['intercept']),
                 ('independent cube-edge compiler agreement', pattern, slope))
        else:
            need(slope == support_two_formula(pattern, defects),
                 ('independent rectangle formula agreement', pattern, slope))
        need(slope <= TARGET, ('uniform finite coefficient-box slope', pattern, slope))
        maximum = max(maximum, slope)
        if slope == TARGET:
            equalities.append(pattern)
        print('pattern', pattern, 'defects', defects, 'slope', slope, 'intercept', intercept,
              'maximizer', witness)
        digest.update(repr((pattern, defects, slope, intercept)).encode())
    print('PASS COMPLETE COEFFICIENT BOX max<=18: A_pattern<=15/98')
    print('UNIVERSE 308 patterns:293 full-support and15 support-two; primitive,even l1,at mostone coordinate0mod3')
    print('EXCLUSIONS equal-speed(0,1,1) and inherited norm-four(1,1,2)')
    print('MAXIMUM', maximum, 'EQUALITY PATTERNS', equalities)
    print('HOSTILE norm-four slope', norm4[1])
    print('COMPLETE SEMANTIC SHA256', digest.hexdigest())
    print('CHECKS', CHECKS, 'INDEPENDENT CUBE-EDGE CHECKS', other.CHECKS)
    print('SCOPE finite coefficient-box theorem; no speed-height census; universal extension requires separate zonotope lemma')


if __name__ == '__main__':
    main()
