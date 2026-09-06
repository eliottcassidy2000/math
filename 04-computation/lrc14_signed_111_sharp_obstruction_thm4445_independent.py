#!/usr/bin/env python3
"""Clean-room referee for the signed (1,1,1) all-height candidate.

No candidate script or repository computation module is imported.  Carrier
rows are solved by a modular one-coordinate compiler, while the proposed ray
values and continuum profiles are rebuilt separately from their definitions.
"""

from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd

R = Q(3, 14)
TARGET = Q(6, 77)
UPPER = Q(6, 55)
LOWER = Q(31, 392)
H = 223
CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def typed(w):
    a, b, c = w
    return (0 < a < b < c and a+b == c and gcd(a, b) == 1
            and all(x % 3 for x in w))


def all_rows(height):
    for c in range(1, height+1):
        if c % 3 == 0:
            continue
        for a in range(1, (c+1)//2):
            w = (a, c-a, c)
            if typed(w):
                yield w


def bounds(w):
    return tuple((3*(sum(w)-w[i])-1)//14 for i in range(3))


def modular_raw_carriers(w):
    """Solve b*y=-a*x mod c; its y-roof interval is shorter than c."""
    a, b, c = w
    bx, by, bz = bounds(w)
    need(gcd(b, c) == 1, ('invertible row', w))
    inv_b = pow(b, -1, c)
    live = set()
    for x in range(-bx, bx+1):
        if x % 3 == 0:
            continue
        residue = (-a*x*inv_b) % c
        y = -by + (residue+by) % c
        if y > by or y % 3 == 0:
            continue
        z = -(a*x+b*y)//c
        need(a*x+b*y+c*z == 0, ('kernel reconstruction', w, x, y, z))
        if z % 3 and abs(z) <= bz:
            live.add((x, y, z))
    return live


def raw_values(w, live):
    c = w[2]
    E = [Q(0)]*3
    physical = Q(0)
    for C in live:
        terms = []
        for i in range(3):
            j, k = [h for h in range(3) if h != i]
            term = min(Q(3, 7*c),
                       Q(3*(w[j]+w[k])-14*abs(C[i]),
                         14*w[j]*w[k]))
            need(term > 0, ('strict raw roof', w, C, i, term))
            E[i] += term
            terms.append(term)
        physical += min(terms)
    return tuple(E), physical


def ray_values(w):
    a, b, c = w
    t = Q(a, c)
    addresses = tuple(k for k in range(1, (3*c-1)//14+1) if k % 3)
    E = [Q(0)]*3
    physical = Q(0)
    for k in addresses:
        s = Q(k, c)
        terms = (
            min(2*R, (R*(2-t)-s)/(1-t)),
            min(2*R, (R*(1+t)-s)/t),
            min(2*R, (R-s)/(t*(1-t))),
        )
        need(min(terms) > 0, ('positive ray sample', w, k, terms))
        for i in range(3):
            E[i] += Q(2, c)*terms[i]
        physical += Q(2, c)*min(terms)
    ray = {(sgn*k, sgn*k, -sgn*k)
           for k in addresses for sgn in (-1, 1)}
    return tuple(E), physical, ray, addresses


def continuum(w):
    """Generic capped-line integration, with crossings generated afresh."""
    t = Q(w[0], w[2])
    lines = (
        (R*(2-t)/(1-t), Q(1)/(1-t)),
        (R*(1+t)/t, Q(1)/t),
        (R/(t*(1-t)), Q(1)/(t*(1-t))),
    )
    points = {Q(0), R}
    for alpha, beta in lines:
        cap = (alpha-2*R)/beta
        if 0 < cap < R:
            points.add(cap)
    for i in range(3):
        for j in range(i):
            ai, bi = lines[i]
            aj, bj = lines[j]
            if bi != bj:
                x = (ai-aj)/(bi-bj)
                if 0 < x < R:
                    points.add(x)
    points = sorted(points)
    integrals = [Q(0)]*3
    physical = Q(0)
    owners = []

    def value(i, x):
        alpha, beta = lines[i]
        return min(2*R, alpha-beta*x)

    for lo, hi in zip(points, points[1:]):
        mid = (lo+hi)/2
        owner = min(range(3), key=lambda i: value(i, mid))
        owners.append((lo, hi, owner))
        physical += (hi-lo)*(value(owner, lo)+value(owner, hi))/2
        for i in range(3):
            integrals[i] += (hi-lo)*(value(i, lo)+value(i, hi))/2
        need(value(1, mid) >= value(0, mid),
             ('b profile dominates a', w, lo, hi))
    return tuple(integrals), physical, tuple(owners)


def update_max(table, key, value, row):
    if key not in table or value > table[key][0]:
        table[key] = (value, [row])
    elif value == table[key][0]:
        table[key][1].append(row)


def update_min(table, key, value, row):
    if key not in table or value < table[key][0]:
        table[key] = (value, [row])
    elif value == table[key][0]:
        table[key][1].append(row)


def first_upper_cutoff(bulk, value):
    c = 1
    while bulk+Q(4, 7*c) >= value:
        c += 1
    return c


def first_lower_cutoff(value):
    c = 1
    while Q(9, 98)-Q(4, 7*c) <= value:
        c += 1
    return c


# Sign classification: four full-support signs modulo global reversal.
signs = {s if s[0] == 1 else tuple(-x for x in s)
         for s in product((-1, 1), repeat=3)}
need(len(signs) == 4, ('sign quotient', signs))
viable = set()
for w in combinations(range(1, 51), 3):
    for u in signs:
        if sum(u[i]*w[i] for i in range(3)) == 0:
            viable.add(u)
need(viable == {(1, 1, -1)}, ('sorted positive sign classification', viable))

# The exact count formula is constant on (m,m+1] and its error extrema occur
# at interval endpoints.  This rational sweep is an orthogonal hostile probe.
for denominator in range(1, 20):
    for numerator in range(1, 400):
        T = Q(numerator, denominator)
        M = (T.numerator+T.denominator-1)//T.denominator-1
        exact = M-M//3
        need(Q(2, 3)*T-Q(2, 3) <= exact < Q(2, 3)*T+Q(2, 3),
             ('two-sided address discrepancy', T, M, exact))

maxima = {}
minima = {}
at_target = {'N': [], 'P': []}
parity_counts = defaultdict(int)
rows = list(all_rows(H))
for w in rows:
    raw_carriers = modular_raw_carriers(w)
    raw_E, raw_P = raw_values(w, raw_carriers)
    ray_E, ray_P, ray_carriers, addresses = ray_values(w)
    need((raw_carriers, raw_E, raw_P) == (ray_carriers, ray_E, ray_P),
         ('clean-room raw/ray mismatch', w, raw_carriers, ray_carriers,
          raw_E, ray_E, raw_P, ray_P))
    need(all(14*k != 3*w[2] for k in addresses),
         ('strict endpoint collision', w, addresses))
    integrals, physical_integral, owners = continuum(w)
    t = Q(w[0], w[2])
    expected = (
        R*R*(3+t)/2,
        None,
        2*R*R*(1-t+t*t),
    )
    need(integrals[0] == expected[0] and integrals[2] == expected[2],
         ('continuum projection formulas', w, integrals, expected))
    need(physical_integral == Q(3, 2)*R*R,
         ('continuum physical formula', w, physical_integral))
    A = Q(4, 3)*integrals[0]
    C = Q(4, 3)*integrals[2]
    P = Q(4, 3)*physical_integral
    err = Q(4, 7*w[2])
    need(A == Q(9, 98)+Q(3, 98)*t and C == Q(12, 98)*(1-t+t*t)
         and P == Q(9, 98), ('bulk formulas', w, A, C, P))
    need(min(A, C)-err <= min(raw_E) < min(A, C)+err,
         ('network two-sided quadrature', w, raw_E, A, C, err))
    need(P-err <= raw_P < P+err,
         ('physical two-sided quadrature', w, raw_P, P, err))
    mask = sum((x % 2) << i for i, x in enumerate(w))
    label = {3: 'c_even', 5: 'b_even', 6: 'a_even'}.get(mask)
    need(label is not None, ('unique even coordinate', w, mask))
    parity_counts[label] += 1
    for kind, value in (('N', min(raw_E)), ('P', raw_P)):
        update_max(maxima, (kind, 'all'), value, w)
        update_max(maxima, (kind, label), value, w)
        if w != (1, 4, 5):
            update_min(minima, kind, value, w)
        if value <= TARGET:
            at_target[kind].append((w, value))

expected_maxima = {
    ('N', 'all'): (UPPER, [(1, 10, 11)]),
    ('N', 'c_even'): (Q(2946, 28861), [(7, 31, 38)]),
    ('N', 'b_even'): (UPPER, [(1, 10, 11)]),
    ('N', 'a_even'): (Q(223, 2156), [(4, 7, 11)]),
    ('P', 'all'): (UPPER, [(1, 10, 11)]),
    ('P', 'c_even'): (Q(222, 2275), [(1, 25, 26)]),
    ('P', 'b_even'): (UPPER, [(1, 10, 11)]),
    ('P', 'a_even'): (Q(102, 1001), [(2, 11, 13)]),
}
expected_minima = {
    'N': (LOWER, [(1, 7, 8)]),
    'P': (LOWER, [(1, 7, 8)]),
}
expected_at_target = {
    'N': [((1, 4, 5), Q(1, 28))],
    'P': [((1, 4, 5), Q(1, 28))],
}
need(maxima == expected_maxima, ('sharp maxima', maxima))
need(minima == expected_minima, ('sharp nonexception minima', minima))
need(at_target == expected_at_target, ('target threshold locus', at_target))
need(dict(parity_counts) == {'b_even': 639, 'a_even': 642, 'c_even': 620},
     ('parity counts', dict(parity_counts)))

network_cutoffs = {label: first_upper_cutoff(Q(39, 392), maxima[('N', label)][0])
                   for label in ('c_even', 'b_even', 'a_even')}
physical_cutoffs = {label: first_upper_cutoff(Q(9, 98), maxima[('P', label)][0])
                    for label in ('c_even', 'b_even', 'a_even')}
need(network_cutoffs == {'c_even': 221, 'b_even': 60, 'a_even': 145},
     ('network cutoff table', network_cutoffs))
need(physical_cutoffs == {'c_even': 100, 'b_even': 34, 'a_even': 57},
     ('physical cutoff table', physical_cutoffs))
need(first_lower_cutoff(TARGET) == 42,
     ('target lower cutoff', first_lower_cutoff(TARGET)))
need(first_lower_cutoff(LOWER) == 45,
     ('sharp lower cutoff', first_lower_cutoff(LOWER)))
need(Q(9, 98)-Q(4, 7*43)-TARGET == Q(29, 46354),
     ('c43 exact margin', Q(9, 98)-Q(4, 7*43)-TARGET))
need(Q(9, 98)-Q(4, 7*45)-LOWER == Q(1, 17640),
     ('c45 exact margin', Q(9, 98)-Q(4, 7*45)-LOWER))

# Independently type-check the three claimed cofinal parity families.
cofinal = defaultdict(int)
for n in range(1, 1000):
    if n % 6 == 1:
        w = (n, 3*n+1, 4*n+1)
        need(typed(w), ('b-even cofinal row', n, w))
        cofinal['b_even'] += 1
    if n % 6 == 4:
        w = (n, 3*n+1, 4*n+1)
        need(typed(w), ('a-even cofinal row', n, w))
        cofinal['a_even'] += 1
    if n % 6 == 5:
        w = (n, 3*n+2, 4*n+2)
        need(typed(w), ('c-even cofinal row', n, w))
        cofinal['c_even'] += 1

print('CLEANROOM_SIGNED_111_REFEREE')
print('sign_classes', len(signs), 'viable', sorted(viable))
print('raw_modular_rows_c_le_223', len(rows), 'parity_counts', dict(parity_counts))
print('sharp_maxima', maxima)
print('sharp_nonexception_minima', minima)
print('rows_at_or_below_6_77', at_target)
print('network_cutoffs', network_cutoffs)
print('physical_cutoffs', physical_cutoffs)
print('lower_cutoffs_arithmetic_target_nonexception', (42, 45),
      'first_admissible_target_height', 43)
print('exact_margins', (Q(29, 46354), Q(1, 17640)))
print('cofinal_family_rows_checked', dict(cofinal))
print('explicit_checks', CHECKS)
print('PASS')
