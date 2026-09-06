#!/usr/bin/env python3
"""Small exact controls for the variable-radius coarea/count lemma.

No speed census: complete F3 fiber table, bounded interval words, named
threshold/parity hostiles, and four samples of each asymptotically sharp family.
All counts retain the complete raw integer support and all defect layers.
"""

from fractions import Fraction as Q
from itertools import product
from math import gcd

CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def cross(x, y):
    return (x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0])


def multiplicities(w, v):
    out = []
    for d in range(3):
        live = [C for C in product(range(3), repeat=3)
                if all(C) and sum(C[i]*w[i] for i in range(3)) % 3 == 0
                and all((cross(v, C)[i]-d*w[i]) % 3 == 0 for i in range(3))]
        out.append(len(live))
    need(out[1] == out[2], ('opposite defect symmetry', w, v, out))
    return tuple(out)


def clip(poly, A, B, C):
    if not poly:
        return []
    out = []
    for x, y in zip(poly, poly[1:]+poly[:1]):
        fx, fy = A*x[0]+B*x[1]-C, A*y[0]+B*y[1]-C
        if fx <= 0:
            out.append(x)
        if fx*fy < 0:
            t = fx/(fx-fy)
            out.append((x[0]+t*(y[0]-x[0]), x[1]+t*(y[1]-x[1])))
    return list(dict.fromkeys(out))


def section_width(w, v, d, r):
    pivot = next(i for i in range(3) if v[i])
    w, v = w[pivot:]+w[:pivot], v[pivot:]+v[:pivot]
    poly = [(-r, -r), (r, -r), (r, r), (-r, r)]
    lo, hi = sorted((-v[0]*r, v[0]*r))
    poly = clip(poly, -v[1], -v[2], hi-d)
    poly = clip(poly, v[1], v[2], d-lo)
    points = [((d-v[1]*y-v[2]*z)/v[0], y, z) for y, z in poly]
    scalars = [cross(w, e)[0]/v[0] for e in points]
    return max(scalars)-min(scalars)


def carriers(w, r):
    bounds = [(r*(sum(w)-x)).__ceil__()-1 for x in w]
    result = set()
    for x in range(-bounds[0], bounds[0]+1):
        if x % 3 == 0:
            continue
        for y in range(-bounds[1], bounds[1]+1):
            if y % 3 == 0 or (w[0]*x+w[1]*y) % w[2]:
                continue
            z = -(w[0]*x+w[1]*y)//w[2]
            if z % 3 and abs(z) <= bounds[2]:
                result.add((x, y, z))
    return result


def audit(w, v, r):
    need(gcd(gcd(w[0], w[1]), w[2]) == 1 and all(x > 0 for x in w), ('speed type', w))
    need(gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2])) == 1
         and sum(w[i]*v[i] for i in range(3)) == 0, ('relation type', w, v))
    S, M = sum(map(abs, v)), max(map(abs, v))
    mult = multiplicities(tuple(x % 3 for x in w), tuple(x % 3 for x in v))
    h = sum(mult)
    kappa = Q(mult[1]+abs(mult[0]-mult[1]), 3)
    bound = (r*S).__ceil__()-1
    F = Q(0)
    B = Q(0)
    D = []
    for d in range(-bound, bound+1):
        m = mult[d % 3]
        if not m:
            continue
        D.append(d)
        F += Q(m, 3)*section_width(w, v, Q(d), r)
        B += Q(1) if m == 1 else Q(4, 3)
    live = carriers(w, r)
    N = len(live)
    f0 = section_width(w, v, Q(0), r)
    integral = 4*r*r*sum(w)
    need(abs(N-F) <= B, ('complete affine integer-count discrepancy', w, v, r, N, F, B))
    need(abs(F-Q(h, 9)*integral) <= kappa*f0,
         ('coarea quadrature discrepancy', w, v, r, F, integral, f0, mult))
    need(abs(N-Q(h, 9)*integral) <= kappa*f0+B,
         ('complete two-sided count', w, v, r))
    i = next(i for i in range(3) if abs(v[i]) == M)
    need(f0 <= 2*r*(sum(w)-w[i])/M, ('central section peak', w, v, r, f0))
    need((not D) == (h == 0 or (mult[0] == 0 and r*S <= 1)), ('exact empty defect list', w, v, r, D))
    print('CONTROL', w, 'v', v, 'r', str(r), 'owner multiplicities', mult,
          'N', N, 'N/c', str(Q(N, max(w))), 'F', str(F), 'B', str(B), 'I', str(integral))
    return live


def main():
    table = set()
    for w in product(range(3), repeat=3):
        if not any(w):
            continue
        for v in product(range(3), repeat=3):
            if not any(v) or sum(w[i]*v[i] for i in range(3)) % 3:
                continue
            mult = multiplicities(w, v)
            zeros = sum(x == 0 for x in w)
            expected = ((2, 0, 0) if all(v) else (0, 1, 1)) if zeros == 0 else (
                       ((2, 1, 1) if all(v) else (0, 2, 2)) if zeros == 1 else (0, 0, 0))
            need(mult == expected, ('complete finite-field table', w, v, mult, expected))
            table.add((zeros, all(v), mult))
    # Sharp interval-word errors, including open endpoints.
    for mask in ((0,), (0, 1), (0, 2)):
        m = len(mask)
        beta = Q(1) if m == 1 else Q(4, 3)
        for A in range(-8, 9):
            for B in range(A+1, 10):
                a, b = Q(A, 3), Q(B, 3)
                count = sum(a < n < b and n % 3 in mask for n in range(-5, 6))
                need(-beta <= count-Q(m, 3)*(b-a) < beta, ('open interval residue error', mask, a, b))
    controls = [((1, 5, 11), (1, 2, -1), Q(3, 14)),
                ((2, 5, 7), (1, 1, -1), Q(3, 14)),
                ((1, 2, 5), (2, -1, 0), Q(1, 3)),
                ((1, 2, 5), (2, -1, 0), Q(3, 8)),
                ((2, 3, 5), (1, 1, -1), Q(3, 14)),
                ((2, 3, 5), (3, -2, 0), Q(3, 14)),
                ((3, 6, 7), (2, -1, 0), Q(1, 2))]
    for w, v, r in controls:
        live = audit(w, v, r)
        if w == (1, 2, 5):
            expected = set() if r == Q(1, 3) else {(-1, -2, 1), (1, 2, -1)}
            need(live == expected, ('sharp empty-list threshold witness', r, live))
    for n in (5, 11, 17, 23):
        audit((n*n, n*(n+2), (n+2)**2), (n+2, -n, 0), Q(3, 14))
    for n in (1, 4, 7, 10):
        audit((n*n+1, n*n+n+1, n*n+2*n+2), (n+1, -2*n-1, n), Q(3, 14))
    print('FINITE-FIELD TABLE', sorted(table))
    print('CHECKS', CHECKS)
    print('PASS: variable-radius coarea/count controls; no speed-height census')


if __name__ == '__main__':
    main()
