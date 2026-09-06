#!/usr/bin/env python3
"""Finite-exact obligations for the all-parity signed (1,1,2) ray theorem.

The proof outside this file supplies the norm-four covector argument and the
deleted-third quadrature.  This program independently reconstructs the raw
carrier head, normalized continuum profiles, exact family leaders, and an
H611 diagnostic extension.  It uses only Python's standard library.
"""

from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import gcd

R = Q(3, 14)
HEAD = 68
EXTENSION = 611
TARGET = Q(6, 77)
GLOBAL_SHARP = Q(11, 140)
CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def family(w):
    a, b, c = w
    out = []
    if 2*a + b == c:
        out.append((1, (2, 1, -1), 0))
    if a + c == 2*b:
        out.append((2, (1, -2, 1), 1))
    if a + 2*b == c:
        out.append((3, (1, 2, -1), 1))
    return out


def typed(w):
    return (w[0] < w[1] < w[2] and gcd(gcd(w[0], w[1]), w[2]) == 1
            and all(x % 3 for x in w))


def rows_to(height):
    rows = []
    for c in range(1, height + 1):
        if c % 3 == 0:
            continue
        for a in range(1, c):
            if a % 3 == 0:
                continue
            candidates = []
            if c - 2*a > a:
                candidates.append(c - 2*a)
            if (a + c) % 2 == 0:
                candidates.append((a + c)//2)
            if (c - a) % 2 == 0:
                candidates.append((c - a)//2)
            for b in set(candidates):
                w = (a, b, c)
                if not typed(w):
                    continue
                ff = family(w)
                need(len(ff) == 1, ('unique family', w, ff))
                rows.append((w, *ff[0]))
    return rows


def roofs(w):
    return tuple((3*(sum(w)-w[i])-1)//14 for i in range(3))


def ray_row(w, u):
    c = w[2]
    bounds = roofs(w)
    K = min(bounds[i]//abs(u[i]) for i in range(3))
    addresses = tuple(k for k in range(1, K+1) if k % 3)
    projections = [Q(0)]*3
    physical = Q(0)
    for k in addresses:
        terms = []
        for i in range(3):
            j, ell = [z for z in range(3) if z != i]
            term = min(Q(3, 7*c),
                       Q(3*(w[j]+w[ell])-14*abs(u[i])*k,
                         14*w[j]*w[ell]))
            need(term > 0, ('positive ray term', w, u, k, i, term))
            terms.append(term)
            projections[i] += 2*term
        physical += 2*min(terms)
    carriers = {(s*k*u[0], s*k*u[1], s*k*u[2])
                for k in addresses for s in (-1, 1)}
    return tuple(projections), physical, carriers, addresses


def raw_row(w):
    a, b, c = w
    bounds = roofs(w)
    live = set()
    d = gcd(b, c)
    modulus = c//d
    inv = pow(b//d, -1, modulus) if modulus > 1 else 0
    for x in range(-bounds[0], bounds[0]+1):
        if x % 3 == 0 or a*x % d:
            continue
        residue = (-a*x//d*inv) % modulus if modulus > 1 else 0
        low = -bounds[1] + (residue+bounds[1]) % modulus
        for y in range(low, bounds[1]+1, modulus):
            if y % 3 == 0:
                continue
            z = -(a*x+b*y)//c
            if z % 3 and abs(z) <= bounds[2]:
                live.add((x, y, z))
    projections = [Q(0)]*3
    physical = Q(0)
    for C in live:
        terms = []
        for i in range(3):
            j, ell = [z for z in range(3) if z != i]
            term = min(Q(3, 7*c),
                       Q(3*(w[j]+w[ell])-14*abs(C[i]),
                         14*w[j]*w[ell]))
            need(term > 0, ('positive raw term', w, C, i, term))
            projections[i] += term
            terms.append(term)
        physical += min(terms)
    return tuple(projections), physical, live


def continuum(w, u):
    """Integrate exact capped-affine profiles on the positive ray."""
    W = tuple(Q(v, w[2]) for v in w)
    B = min(R*(W[j]+W[k])/abs(u[i])
            for i in range(3)
            for j, k in [tuple(z for z in range(3) if z != i)])
    lines = []
    points = {Q(0), B}
    for i in range(3):
        j, k = [z for z in range(3) if z != i]
        alpha = R*(W[j]+W[k])/(W[j]*W[k])
        beta = Q(abs(u[i]), 1)/(W[j]*W[k])
        lines.append((alpha, beta))
        cap = (alpha-2*R)/beta
        if 0 < cap < B:
            points.add(cap)
    for i in range(3):
        for j in range(i):
            ai, bi = lines[i]
            aj, bj = lines[j]
            if bi != bj:
                x = (ai-aj)/(bi-bj)
                if 0 < x < B:
                    points.add(x)
    points = sorted(points)
    integrals = [Q(0)]*3
    envelope = Q(0)
    for lo, hi in zip(points, points[1:]):
        mid = (lo+hi)/2

        def value(i, x):
            alpha, beta = lines[i]
            return min(2*R, alpha-beta*x)

        owner = min(range(3), key=lambda i: value(i, mid))
        envelope += (hi-lo)*(value(owner, lo)+value(owner, hi))/2
        for i in range(3):
            integrals[i] += (hi-lo)*(value(i, lo)+value(i, hi))/2
    return B, tuple(integrals), envelope


def expected_integrals(fam, t):
    if fam == 1:
        normalized = (Q(1), 2*(1-t), 2-4*t+4*t*t)
    elif fam == 2:
        normalized = ((1+2*t)/(1+t), Q(1),
                      (1+2*t+t*t-t*t*t)/(1+t))
    else:
        normalized = (1+t, Q(1), 1+t*t)
    return tuple(R*R*x for x in normalized)


def update(leaders, key, value, payload):
    if key not in leaders or value > leaders[key][0]:
        leaders[key] = (value, [payload])
    elif value == leaders[key][0]:
        leaders[key][1].append(payload)


# Twelve coefficient sign placements modulo global reversal, only three cones.
directions = set()
for mag in set(permutations((1, 1, 2))):
    for signs in product((-1, 1), repeat=3):
        u = tuple(mag[i]*signs[i] for i in range(3))
        if u[0] < 0:
            u = tuple(-x for x in u)
        directions.add(u)
need(len(directions) == 12, ('normalized directions', directions))
viable = set()
for w in combinations([x for x in range(1, 101) if x % 3], 3):
    if gcd(gcd(w[0], w[1]), w[2]) != 1:
        continue
    for u in directions:
        if sum(u[i]*w[i] for i in range(3)) == 0:
            viable.add(u)
expected_viable = {(2, 1, -1), (1, -2, 1), (1, 2, -1)}
need(viable == expected_viable, ('viable sign cones', viable))

rows = rows_to(EXTENSION)
head = [row for row in rows if row[0][2] <= HEAD]
counts = defaultdict(int)
head_counts = defaultdict(int)
leaders = {}
head_leaders = {}
target_gt = []
target_eq = []

for w, fam, u, selected in rows:
    counts[fam] += 1
    projections, physical, carriers, addresses = ray_row(w, u)
    need(sum(u[i]*w[i] for i in range(3)) == 0, ('relation', w, u))
    need(physical == projections[selected] == min(projections),
         ('pointwise selector consequence', w, projections, physical, selected))
    bulk = Q(3, 49)
    error = Q(4, 7*w[2])
    need(bulk-error <= physical < bulk+error,
         ('two-sided deleted-third band', w, physical, bulk, error))
    B, integrals, envelope = continuum(w, u)
    t = Q(w[0], w[2])
    need(integrals == expected_integrals(fam, t),
         ('continuum formula', w, integrals, expected_integrals(fam, t)))
    need(envelope == integrals[selected] == R*R,
         ('continuum selector', w, selected, envelope, integrals))
    expected_B = ((1-t)*R if fam == 1 else (1+t)*R/2)
    need(B == expected_B, ('support endpoint', w, B, expected_B))
    update(leaders, fam, physical, (w, projections, addresses))
    if w[2] <= HEAD:
        head_counts[fam] += 1
        update(head_leaders, fam, physical, (w, projections, addresses))
        raw_projections, raw_physical, raw_carriers = raw_row(w)
        need((projections, physical, carriers) ==
             (raw_projections, raw_physical, raw_carriers),
             ('raw versus ray', w, (projections, physical, carriers),
              (raw_projections, raw_physical, raw_carriers)))
    if physical > TARGET:
        target_gt.append((w, physical))
    elif physical == TARGET:
        target_eq.append((w, physical))

expected_counts = {1: 9482, 2: 14214, 3: 4742}
need(dict(counts) == expected_counts, ('H611 family counts', dict(counts)))
need(target_gt == [((2, 11, 20), Q(11, 140))], ('target failures', target_gt))
need(target_eq == [((1, 5, 11), Q(6, 77))], ('target equalities', target_eq))
need(leaders[1][0] == Q(58, 833) and
     leaders[1][1] == [((5, 7, 17),
                        (Q(58, 833), Q(12, 119), Q(346, 4165)), (1, 2))],
     ('family 1 sharp head', leaders[1]))
need(leaders[2][0] == GLOBAL_SHARP and
     leaders[2][1] == [((2, 11, 20),
                        (Q(131, 1540), Q(11, 140), Q(3, 35)), (1, 2))],
     ('family 2 sharp head', leaders[2]))
need(leaders[3][0] == TARGET and
     leaders[3][1] == [((1, 5, 11),
                        (TARGET, TARGET, TARGET), (1,))],
     ('family 3 sharp head', leaders[3]))
need(head_leaders == leaders, ('H68 contains every H611 leader',
                               head_leaders, leaders))

# The deleted-third bound 3/49+4/(7c) beats each frozen leader beyond
# c=68,32,34 respectively.  HEAD=68 therefore proves all three sharp bounds.
cutoffs = {1: (68, Q(58, 833)), 2: (Q(560, 17), GLOBAL_SHARP),
           3: (Q(308, 9), TARGET)}
for fam, (threshold, leader) in cutoffs.items():
    need(Q(3, 49) + Q(4, 7)/(threshold if isinstance(threshold, int)
                              else threshold) == leader,
         ('tail crossing', fam, threshold, leader))

print('normalized_sign_directions', len(directions),
      'viable_sorted_cones', sorted(viable))
print('H611_rows', len(rows), 'counts_F1_F2_F3', dict(counts))
print('H68_rows', len(head), 'counts_F1_F2_F3', dict(head_counts))
for fam in (1, 2, 3):
    print('family', fam, 'sharp', leaders[fam])
print('global_sharp', GLOBAL_SHARP, 'unique_at', (2, 11, 20))
print('target_6_over_77_failures', target_gt)
print('target_6_over_77_equalities', target_eq)
print('selected_and_physical_continuum_integral', R*R)
print('deleted_third_tail', '3/49 + 4/(7c)',
      'crossing_heights_F1_F2_F3', (68, Q(560, 17), Q(308, 9)))
print('explicit_checks', CHECKS)
print('PASS')
