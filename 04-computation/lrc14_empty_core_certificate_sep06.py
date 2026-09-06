#!/usr/bin/env python3
"""Exact proof compiler for full-support l1<=20 network certificates.

The infinite tail is reduced analytically to coefficient-dependent finite
heads. This script checks the rational polytope maxima and every head row;
it is not an extrapolation of a height census. No repository math is imported
by the producer. The literal-sheet cross-check is explicitly an independent
implementation inherited from the one-ray verifier.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement, permutations, product
from math import gcd

CHECKS = 0
R = Q(3, 14)
TARGET = Q(6, 77)


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def cube_slice(v, delta, low, high):
    """All vertices of v.x=delta intersect [low,high]^3, exact rationals."""
    vertices = set()
    for k in range(3):
        i, j = [t for t in range(3) if t != k]
        for xi, xj in product((low, high), repeat=2):
            point = [Q(0)]*3
            point[i], point[j] = xi, xj
            point[k] = (delta-v[i]*xi-v[j]*xj)/v[k]
            if low <= point[k] <= high:
                vertices.add(tuple(point))
    return tuple(sorted(vertices))


def cross(u, v):
    return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0])


def width(v, w, vertices):
    # All differences w cross (e-e') are scalar multiples of primitive v.
    values = [cross(w, e)[0]/v[0] for e in vertices]
    return max(values)-min(values)


def coefficient_patterns(max_norm=20):
    return tuple(a for norm in range(4, max_norm+1, 2)
                 for a in combinations_with_replacement(range(1, norm-1), 3)
                 if sum(a) == norm and gcd(gcd(a[0], a[1]), a[2]) == 1
                 and sum(x % 3 == 0 for x in a) <= 1)


def sectors(pattern):
    return tuple(sorted(set(tuple(-a[i] if i == k else a[i]
                                  for i in range(3))
                            for a in permutations(pattern) for k in range(3))))


def compile_pattern(pattern):
    unit = all(x % 3 for x in pattern)
    bound = (3*sum(pattern)-1)//14
    deltas = tuple(d for d in range(-bound, bound+1)
                   if (d % 3 == 0) == unit)
    need(bool(deltas), ('nonempty defect list', pattern))
    rho = Q(2 if unit else 1, 3)
    remainder = Q(4, 3) if unit else Q(1)
    best = Q(0)
    witness = None
    # Coordinate permutations preserve the maximum; three isolated signs
    # suffice for the full speed cube, which already allows all orderings.
    for k in range(3):
        v = tuple(-pattern[i] if i == k else pattern[i] for i in range(3))
        speed_vertices = cube_slice(v, Q(0), Q(0), Q(1))
        error_polygons = [cube_slice(v, Q(d), -R, R) for d in deltas]
        need(all(error_polygons), ('nonempty error slice', pattern, v))
        for w in speed_vertices:
            value = rho*sum(width(v, w, p) for p in error_polygons)
            if value > best:
                best, witness = value, (v, w)
    intercept = remainder*len(deltas)
    cutoff = None if best >= Q(2, 11) else intercept/(Q(2, 11)-best)
    return dict(pattern=pattern, unit=unit, deltas=deltas, slope=best,
                intercept=intercept, cutoff=cutoff, witness=witness)


def raw_carriers(w):
    a, b, c = w
    bounds = tuple((3*(sum(w)-wi)-1)//14 for wi in w)
    d = gcd(b, c)
    modulus = c//d
    inverse = pow(b//d, -1, modulus) if modulus > 1 else 0
    answer = set()
    for x in range(-bounds[0], bounds[0]+1):
        if x % 3 == 0 or a*x % d:
            continue
        residue = (-a*x//d*inverse) % modulus if modulus > 1 else 0
        y0 = -bounds[1]+(residue+bounds[1]) % modulus
        for y in range(y0, bounds[1]+1, modulus):
            if y % 3 == 0:
                continue
            z = -(a*x+b*y)//c
            if z % 3 and abs(z) <= bounds[2]:
                answer.add((x, y, z))
    return answer


def projection_data(w, live):
    cap = Q(3, 7*w[2])
    sums = [Q(0)]*3
    physical = Q(0)
    for C in live:
        terms = []
        for i in range(3):
            j, k = [t for t in range(3) if t != i]
            term = min(cap, Q(3*(w[j]+w[k])-14*abs(C[i]), 14*w[j]*w[k]))
            need(term > 0, ('positive roof', w, C, i))
            sums[i] += term
            terms.append(term)
        physical += min(terms)
    return tuple(sums), physical


def primitive(C):
    divisor = gcd(gcd(abs(C[0]), abs(C[1])), abs(C[2]))
    v = tuple(x//divisor for x in C)
    return v if v[0] > 0 else tuple(-x for x in v)


def head_rows(record):
    # Tail is strict already at equality c=cutoff; head is c<cutoff.
    cutoff = record['cutoff']
    limit = (cutoff.numerator-1)//cutoff.denominator
    eligible = tuple(x for x in range(1, limit+1, 2) if x % 3)
    rows = set()
    for v in sectors(record['pattern']):
        for i, a in enumerate(eligible):
            for b in eligible[i+1:]:
                numerator = -(v[0]*a+v[1]*b)
                if numerator % v[2]:
                    continue
                c = numerator//v[2]
                if b < c <= limit and c % 2 and c % 3 and gcd(gcd(a, b), c) == 1:
                    rows.add((a, b, c))
    return rows


def relation_layer_audit(w, v, live):
    unit = all(x % 3 for x in v)
    layers = Counter()
    for C in live:
        vector = cross(v, C)
        delta = Q(vector[0], w[0])
        need(delta.denominator == 1 and vector == tuple(delta*t for t in w),
             ('integral affine defect', w, v, C, delta))
        delta = int(delta)
        need(abs(delta) < R*sum(map(abs, v)) and (delta % 3 == 0) == unit,
             ('complete allowed defect', w, v, C, delta))
        layers[delta] += 1
    for delta, count in layers.items():
        vertices = cube_slice(v, Q(delta), -R, R)
        length = width(v, w, vertices)
        rho, remainder = (Q(2, 3), Q(4, 3)) if unit else (Q(1, 3), Q(1))
        need(count < rho*length+remainder,
             ('strict residue interval bound', w, v, delta, count, length))


def main():
    records = [compile_pattern(p) for p in coefficient_patterns()]
    need(len(records) == 73, ('complete patterns count', len(records)))
    exceptional = [r for r in records if r['cutoff'] is None]
    need([r['pattern'] for r in exceptional] == [(1, 1, 2)], 'sole norm-four exception')
    records = [r for r in records if r['cutoff'] is not None]
    need(max(r['cutoff'] for r in records) == Q(4312, 31), 'largest finite cutoff')

    proof_rows = set()
    memberships = {}
    for record in records:
        rows = head_rows(record)
        record['head_count'] = len(rows)
        proof_rows.update(rows)
        for w in rows:
            memberships.setdefault(w, []).append(record)
        print('sector', record['pattern'], 'rho', '2/3' if record['unit'] else '1/3',
              'defects', record['deltas'], 'slope', record['slope'],
              'intercept', record['intercept'], 'cutoff', record['cutoff'],
              'head', len(rows), 'maximizer', record['witness'])

    # The external implementation uses literal six-sheet contact graphs and
    # is independent of all polytope, relation-box, and affine-layer code here.
    from lrc14_one_ray_overnight_hexagon_sep05 import literal_projection_data
    digest = sha256()
    kinds = Counter()
    equalities = set()
    max_non_norm4 = (Q(0), None)
    max_multi = (Q(0), None)
    for w in sorted(proof_rows, key=lambda t:(t[2], t[0], t[1])):
        live = raw_carriers(w)
        projected, physical = projection_data(w, live)
        need((projected, physical) == literal_projection_data(w),
             ('independent literal-sheet agreement', w))
        need(min(projected) <= TARGET, ('finite head network target', w, projected))
        if min(projected) == TARGET:
            equalities.add(w)
        directions = {primitive(C) for C in live}
        norm4 = any(sum(map(abs, v)) == 4 for v in directions)
        kinds['empty' if not directions else 'norm4' if norm4 else
              'other_one' if len(directions) == 1 else 'multi'] += 1
        if not norm4:
            need(max(projected) < TARGET, ('all projections outside norm4', w, projected))
            if max(projected) > max_non_norm4[0]:
                max_non_norm4 = max(projected), w
        if len(directions) >= 2 and max(projected) > max_multi[0]:
            max_multi = max(projected), w
        for record in memberships[w]:
            candidates = [v for v in sectors(record['pattern'])
                          if sum(v[i]*w[i] for i in range(3)) == 0]
            need(bool(candidates), ('honest coefficient membership', w, record['pattern']))
            for v in candidates:
                relation_layer_audit(w, v, live)
            need(len(live) < record['slope']*w[2]+record['intercept'],
                 ('compiled affine count envelope', w, record['pattern'], len(live)))
        digest.update(repr((w, tuple(sorted(live)), projected, physical)).encode())

    need(equalities == {(1, 5, 11)}, ('unique target equality in complete heads', equalities))
    # Independent positive and hostile controls, including unbounded directions.
    controls = ((1, 5, 7), (1, 5, 11), (19, 23, 29), (41, 47, 49),
                (1, 401, 601), (1, 1201, 1801))
    for w in controls:
        live = raw_carriers(w)
        projected, physical = projection_data(w, live)
        need((projected, physical) == literal_projection_data(w), ('wide control', w))
        for v in ((1, -3, 2),):
            if sum(v[i]*w[i] for i in range(3)) == 0:
                relation_layer_audit(w, v, live)
                need(len(live) < Q(2*w[2], 21)+2, ('norm-six triangle count', w))
        print('control', w, 'N', len(live), 'directions', len({primitive(C) for C in live}),
              'projections', projected, 'physical', physical)
    need(projection_data((1, 5, 7), raw_carriers((1, 5, 7)))[0][1] > TARGET,
         'norm-four exception is necessary')
    print('patterns', len(records)+1, 'certified_count_patterns', len(records))
    print('head_memberships', sum(r['head_count'] for r in records))
    print('unique_head_rows', len(proof_rows), 'support_classes', dict(sorted(kinds.items())))
    print('head_equalities', tuple(sorted(equalities)))
    print('head_max_non_norm4', max_non_norm4, 'head_max_multi', max_multi)
    print('complete_head_digest', digest.hexdigest())
    print('explicit_checks', CHECKS)
    print('PASS: all full-support primitive relation shells through l1=20 satisfy the network target')


if __name__ == '__main__':
    main()
