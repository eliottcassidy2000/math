#!/usr/bin/env python3
"""Exact controls for the quadratic quotient of three-dimensional collision periods.

Standard-library only. No inherited producer, large census, or actual Keller map.
The universal statements are proved in the matching note.
"""
from fractions import Fraction as F
from itertools import product
import hashlib
import json

GATES = 0


def check(test, label):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(label)


def rref(rows, width=None):
    a = [[F(x) for x in row] for row in rows]
    n = len(a[0]) if a else width or 0
    pivots = []
    k = 0
    for j in range(n):
        p = next((i for i in range(k, len(a)) if a[i][j]), None)
        if p is None:
            continue
        a[k], a[p] = a[p], a[k]
        d = a[k][j]
        a[k] = [x / d for x in a[k]]
        for i in range(len(a)):
            if i != k and a[i][j]:
                d = a[i][j]
                a[i] = [x - d*y for x, y in zip(a[i], a[k])]
        pivots.append(j)
        k += 1
        if k == len(a):
            break
    return a, pivots


def rank(rows):
    return len(rref(rows)[1])


def nullspace(rows, width):
    a, pivots = rref(rows, width)
    out = []
    for j in range(width):
        if j not in pivots:
            v = [F(0)] * width
            v[j] = F(1)
            for i, p in enumerate(pivots):
                v[p] = -a[i][j]
            out.append(v)
    return out


def matvec(a, x):
    return [sum(u*v for u, v in zip(row, x)) for row in a]


def transpose(a):
    return list(map(list, zip(*a)))


def cross(a, b):
    x, y, z = a
    u, v, w = b
    return (y*w-z*v, z*u-x*w, x*v-y*u)


def cross_matrix(t):
    x, y, z = t
    return ((0, -z, y), (z, 0, -x), (-y, x, 0))


def veronese(t):
    x, y, z = t
    return (x*x, y*y, z*z, x*y, x*z, y*z)


def collision_data(ts):
    m = len(ts)
    check(rank(transpose(ts)) == 3, 'declared tangent universe spans dimension three')
    rels = nullspace(transpose(ts), m)
    per = []
    for ell in rels:
        for j in range(3):
            per.append([ell[i]*cross_matrix(ts[i])[j][k]
                        for i in range(m) for k in range(3)])
    common = []
    for i, t in enumerate(ts):
        for j in range(3):
            common.append([int(j == k) for k in range(3)] +
                          [t[j] if k == i else 0 for k in range(m)])
    check(rank(common) == m+3, 'common motion and independent tangent gauges')
    check(all(not any(matvec(per, col)) for col in transpose(common)),
          'all collision motions have zero full scalar periods')
    rq, rp = rank([veronese(t) for t in ts]), rank(per)
    miss = 3*m-rp-rank(common)
    check(miss == 6-rq, 'exact quadratic quotient dimension')
    return dict(m=m, relation_dim=len(rels), period_rank=rp,
                quadratic_rank=rq, undetected=miss), per, common


def test_hostile(ts, ns):
    data, per, common = collision_data(ts)
    flat = [x for n in ns for x in n]
    check(not any(matvec(per, flat)), 'hostile keeps every scalar relation and form slot')
    check(rank([row+[v] for row, v in zip(common, flat)]) > rank(common),
          'hostile has no common motion even with tangent reparametrization')
    return data


def trim(p):
    p = list(map(F, p))
    while p and not p[-1]:
        p.pop()
    return p


def poly_eval(p, x):
    out = F(0)
    for a in reversed(p):
        out = out*x+a
    return out


def poly_div(p, q):
    p, q = trim(p), trim(q)
    quot = [F(0)] * max(0, len(p)-len(q)+1)
    while len(p) >= len(q):
        j, a = len(p)-len(q), p[-1]/q[-1]
        quot[j] = a
        for k, b in enumerate(q):
            p[j+k] -= a*b
        p = trim(p)
    return trim(quot), p


def main():
    e1, e2, e3 = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    complete = [e1, e2, e3, (1, 1, 0), (1, 0, 1), (0, 1, 1)]
    universes = {
        'three_independent': [e1, e2, e3],
        'four_general': [e1, e2, e3, (1, 1, 1)],
        'five_coordinate_controls': complete[:5],
        'six_complete': complete,
        'seven_complete': complete+[(1, 1, 1)],
        'eight_with_repeated_directions': complete+[(2, 0, 0), (0, -3, 0)],
        'five_conic': [(1, a, a*a) for a in range(-2, 3)],
        'nine_conic': [(1, a, a*a) for a in range(-4, 5)],
    }
    dims = {name: collision_data(ts)[0] for name, ts in universes.items()}
    check([dims[name]['undetected'] for name in universes] == [3, 2, 1, 0, 0, 0, 1, 1],
          'named universe dimensions')
    check(rank([veronese(t) for t in complete]) == 6, 'sharp six-direction positive control')
    check(cross(e1, (0, 0, -1)) == (0, 1, 0), 'three-branch hostile first wedge')
    test_hostile([e1, e2, e3], [(0, 0, -1), (0, 0, 1), (0, 0, 0)])

    conic_controls = []
    for parameters in ((-1, 0, 1), (-2, -1, 0, 1, 2), tuple(range(-4, 5))):
        ts = [(1, a, a*a) for a in parameters]
        ns = [(0, 1, 2*a) for a in parameters]
        data = test_hostile(ts, ns)
        for a, t, n in zip(parameters, ts, ns):
            check(cross(t, n) == (a*a, -2*a, 1), 'literal conic wedge identity')
            check(sum(u*v for u, v in zip((0, -2, 0, 0, 2, 0), veronese(t))) == 0,
                  'missing quadratic evaluated at actual tangent')
        w = (2, -3, 5)
        gauged = [tuple(n[j]+w[j]+(i-2)*t[j] for j in range(3))
                  for i, (t, n) in enumerate(zip(ts, ns))]
        test_hostile(ts, gauged)
        g = ((1, 2, 0), (0, 1, 1), (1, 0, 1))
        test_hostile([matvec(g, t) for t in ts], [matvec(g, n) for n in ns])
        test_hostile([tuple((i+1)*x for x in t) for i, t in enumerate(ts)], ns)
        conic_controls.append(data)

    # The cofactor identity is symbolic in the unspecified actual r' values.
    c, e = (3, 3, 3), (-9, 4, 9)
    cofactors = (c[1]*e[2]-c[2]*e[1],
                 c[2]*e[0]-c[0]*e[2],
                 c[0]*e[1]-c[1]*e[0])
    check(cofactors == (15, -54, 39), 'determinant equals three times retained period')
    # These are synthetic jet controls, not evaluations of the degree175 r.
    graph_cases = 0
    for jets in ((1, 0, 0), (0, 0, 1), (1, 2, 3)):
        ts = [(3, e[i], jets[i]) for i in range(3)]
        data, per, common = collision_data(ts)
        check(data['undetected'] == 3 and not per,
              'independent suspension tangents lose all scalar relation tests')
        vertical = [[int(j == 2 and i == k) for k in range(3)]
                    for i in range(3) for j in range(3)]
        check(rank([row+extra for row, extra in zip(common, vertical)])-rank(common) == 2,
              'vertical graph slice has exactly two unpaid differences')
        for hs in product((-1, 0, 1), repeat=3):
            flat = [hs[i] if j == 2 else 0 for i in range(3) for j in range(3)]
            compatible = rank([row+[v] for row, v in zip(common, flat)]) == rank(common)
            check(compatible == (hs[0] == hs[1] == hs[2]),
                  'all retained cube controls match equal-value graph criterion')
            graph_cases += 1

    L = (0, -1, 0, 1)
    graph_polys = {
        'zero': (), 'constant': (1,), 'x_hostile': (0, 1),
        'L_positive': L, 'old_surface_normal_is_different': (0, -9, 4),
        'one_plus_L': (1, -1, 0, 1),
        'xL_positive': (0, 0, -1, 0, 1),
        'x_squared_hostile': (0, 0, 1),
    }
    graph_values = {}
    for name, p in graph_polys.items():
        values = [poly_eval(p, i) for i in (-1, 0, 1)]
        shifted = list(p) or [F(0)]
        shifted[0] -= values[1]
        quotient, remainder = poly_div(shifted, L)
        check((not remainder) == (values[0] == values[1] == values[2]),
              'monic divisor exactly detects fixed-section graph collision')
        graph_values[name] = [str(v) for v in values]
    check(graph_values['x_hostile'] == ['-1', '0', '1'], 'actual source graph hostile')
    check(graph_values['L_positive'] == ['0', '0', '0'], 'nonconstant exact collision control')

    manifest = dict(named_tangent_universe=dims, conic_controls=conic_controls,
                    suspension_synthetic_jet_controls=3, retained_graph_cube_controls=graph_cases,
                    graph_polynomial_controls=graph_values,
                    scope='universal first-order quotient proved analytically; exact suspension graph iff; no Keller claim')
    encoded = json.dumps(manifest, sort_keys=True, separators=(',', ':'))
    print('PROVED ANALYTICALLY: scalar collision periods lose exactly the tangent quadrics in dimension three')
    print(encoded)
    print('EXPLICIT_GATES', GATES)
    print('SEMANTIC_SHA256', hashlib.sha256(encoded.encode()).hexdigest())


if __name__ == '__main__':
    main()
