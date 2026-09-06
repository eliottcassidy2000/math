#!/usr/bin/env python3
"""Exact finite heads for the 3--6 direction network reduction.

The main enumerator solves the b-coordinate congruence. A separately coded
c-coordinate sweep solves for the a-coordinate and checks the complete support
before any ray-count filter. Literal six-sheet interval graphs, supplied by
the earlier one-ray verifier, independently check all selected projection rows.
Only Python's standard library is required. All gates survive python -O.
"""

import importlib.util
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path

CHECKS = 0
TARGET = Q(6, 77)
CUTOFFS = {3: 99, 4: 159, 5: 233, 6: 319}


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def primitive(C):
    g = gcd(gcd(abs(C[0]), abs(C[1])), abs(C[2]))
    v = tuple(x // g for x in C)
    return v if v[0] > 0 else tuple(-x for x in v)


def carriers_by_first_row(w):
    a, b, c = w
    bx, by, bz = [(3 * (sum(w) - v) - 1) // 14 for v in w]
    d = gcd(b, c)
    modulus = c // d
    inverse = pow(b // d, -1, modulus) if modulus > 1 else 0
    out = set()
    for x in range(-bx, bx + 1):
        if x % 3 == 0 or a * x % d:
            continue
        residue = -(a * x // d) * inverse % modulus
        first = -by + (residue + by) % modulus
        for y in range(first, by + 1, modulus):
            z = -(a * x + b * y) // c
            if y % 3 and z % 3 and abs(z) <= bz:
                out.add((x, y, z))
    return out


def carriers_by_last_row(w):
    """Independent cyclic row sweep; no use of first-row enumerator."""
    a, b, c = w
    d = gcd(a, b)
    period = b // d
    inv = pow(a // d, -1, period) if period > 1 else 0
    out = set()
    zbound = (3 * (a + b) - 1) // 14
    xbound = (3 * (b + c) - 1) // 14
    for z in range(1, zbound + 1):
        if z % 3 == 0 or c * z % d:
            continue
        residue = inv * (-c * z // d) % period
        low = (-xbound - residue + period - 1) // period
        high = (xbound - residue) // period
        for t in range(low, high + 1):
            x = residue + period * t
            y = (-c * z - a * x) // b
            if x % 3 and y % 3 and 14 * abs(y) < 3 * (a + c):
                out.add((x, y, z))
                out.add((-x, -y, -z))
    return out


def projections(w, live):
    denominator = 14 * w[0] * w[1] * w[2]
    cap = 6 * w[0] * w[1]
    totals = [0, 0, 0]
    for C in live:
        for i in range(3):
            margin = 3 * (sum(w) - w[i]) - 14 * abs(C[i])
            need(margin > 0, ('positive roof', w, C, i))
            totals[i] += min(cap, w[i] * margin)
    return tuple(Q(total, denominator) for total in totals)


def tail_gate(r, c):
    return (c >= 33 and 6 * c > 44 * r
            and 49 * (6 * c - 44 * r) ** 2 > 11616 * r * r * c
            and 162 * c > 539 * (12 * r - 8))


def audit(w, live, directions):
    c = w[2]
    r = len(directions)
    maxima = []
    reciprocal = Q(0)
    for v in directions:
        M = max(map(abs, v))
        maxima.append(M)
        reciprocal += Q(1, M)
        need(sum(map(abs, v)) >= 16 and M >= 7,
             ('short-relation obstruction', w, v))
        B = min(Q(3 * (sum(w) - w[i]), 14 * abs(v[i])) for i in range(3))
        K = (B.numerator - 1) // B.denominator
        actual = {C for C in live if primitive(C) == v}
        predicted = {tuple(sign * k * x for x in v)
                     for k in range(1, K + 1) if k % 3
                     for sign in (-1, 1)}
        need(actual == predicted, ('all raw multipliers', w, v))
        need(B < Q(3 * c, 7 * M), ('strict coordinate roof', w, v))
        need(len(actual) < Q(4, 3) * (B + 1), ('mod-three block count', w, v))
    for i, j in combinations(range(r), 2):
        u, v = directions[i], directions[j]
        cross = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2],
                 u[0]*v[1]-u[1]*v[0])
        k = Q(cross[2], c)
        need(k and k.denominator == 1 and k.numerator % 3 == 0,
             ('owner determinant', w, u, v, cross))
        need(cross == tuple(k * x for x in w), ('oriented determinant', w, u, v))
        need(2 * maxima[i] * maxima[j] >= 3 * c, ('pairwise product', w, u, v))
    P = Q(3 * c, 2)
    # The max(r/sqrt(P), 1/7+7(r-1)/P) bound is checked without radicals.
    if P >= 49:
        linear_branch = Q(1, 7) + Q(7 * (r - 1), 1) / P
        need(reciprocal <= linear_branch or reciprocal * reciprocal * P <= r * r,
             ('reciprocal extremal lemma', w, maxima, reciprocal))
    values = projections(w, live)
    need(max(values) < Q(12, 49) * reciprocal + Q(4 * r, 7 * c),
         ('complete-count projection envelope', w, values))
    if tail_gate(r, c):
        need(max(values) < TARGET, ('analytic tail control', w, values))
    return values


def main():
    literal_path = Path(__file__).with_name('lrc14_one_ray_overnight_hexagon_sep05.py')
    spec = importlib.util.spec_from_file_location('literal_sheets', literal_path)
    literal = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(literal)
    for r, c in CUTOFFS.items():
        need(tail_gate(r, c) and not tail_gate(r, c-1), ('exact cutoff', r, c))
    # Positive coefficients after the shift r=s+2 prove the general cutoff.
    for r in range(2, 101):
        c = r * (7*r + 13)
        need(49*(6*c-44*r)**2-11616*r*r*c
             == 4*r*r*(1281*r*r-2766*r+14161), ('first polynomial identity', r))
        need(162*c-539*(12*r-8)
             == 2*(567*r*r-2181*r+2156), ('second polynomial identity', r))
        need(tail_gate(r, c), ('general sufficient cutoff', r))
    counts = Counter()
    universe = Counter()
    leaders = {r: (Q(0), None, None) for r in CUTOFFS}
    raw_leaders = {r: (0, None) for r in CUTOFFS}
    dense = []
    controls = {}
    digest = sha256()
    totals = 0
    eligible = [v for v in range(1, max(CUTOFFS.values()), 2) if v % 3]
    for w in combinations(eligible, 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        totals += 1
        for r, c in CUTOFFS.items():
            if w[2] < c:
                universe[r] += 1
        live = carriers_by_first_row(w)
        need(live == carriers_by_last_row(w), ('complete two-row enumeration', w))
        if w[2] < 99:
            need(live == literal.carriers(w), ('independent low integer box', w))
        directions = tuple(sorted({primitive(C) for C in live}))
        r = len(directions)
        if r in CUTOFFS and w[2] >= CUTOFFS[r] and r not in controls:
            controls[r] = w
        if r not in CUTOFFS or w[2] >= CUTOFFS[r]:
            continue
        counts[r] += 1
        values = audit(w, live, directions)
        need(values == literal.literal_projection_data(w)[0],
             ('independent literal six-sheet graph', w, values))
        need(max(values) < TARGET, ('complete finite proof head', w, r, values))
        if max(values) > leaders[r][0]:
            leaders[r] = max(values), w, values
        if len(live) > raw_leaders[r][0]:
            raw_leaders[r] = len(live), w
        if 11 * len(live) > 2 * w[2]:
            dense.append((r, w, len(live), values))
        digest.update(repr((w, tuple(sorted(live)), directions, values)).encode())
    need(dict(counts) == {3: 1791, 4: 5004, 5: 13656, 6: 22971},
         ('complete head cardinalities', counts))
    print('PROVED: complete supports with 3,4,5,6 primitive directions have every E_i<6/77')
    print('GENERAL TAIL: r>=2 and c>=7r^2+13r imply every E_i<6/77')
    print('EXACT ENVELOPE TAIL CUTOFFS', CUTOFFS)
    print('COMPLETE ENUMERATION', totals, 'primitive sorted distinct positive odd ternary-unit triples c<319')
    for r in CUTOFFS:
        print('HEAD', r, 'rays; universe', universe[r], 'selected', counts[r],
              'every-projection leader', leaders[r], 'maximum raw count', raw_leaders[r])
    print('COUNT-GATE MISSES WITHIN DECLARED PROOF HEADS', dense)
    for w in list(controls.values()) + [(1, 5, 319), (107, 127, 149), (41, 47, 49)]:
        live = carriers_by_first_row(w)
        directions = tuple(sorted({primitive(C) for C in live}))
        values = audit(w, live, directions)
        need(values == literal.literal_projection_data(w)[0], ('literal tail/hostile control', w))
        if w == (41, 47, 49):
            need(directions == ((11, 5, -14), (14, -7, -5), (17, -19, 4)),
                 ('arithmetic progression hostile', w, directions))
            need(not any(tuple(s*u[i]+t*v[i] for i in range(3)) in live
                         for u, v in combinations(directions, 2)
                         for s in (-1, 1) for t in (-1, 1)),
                 ('no additive three-ray circuit', w, directions))
        print('CONTROL', w, 'rays', len(directions), 'raw count', len(live),
              'directions', directions, 'E', tuple(map(str, values)))
    print('SEMANTIC SHA256', digest.hexdigest())
    print('CHECKS', CHECKS, 'LITERAL CHECKS', literal.CHECKS)
    print('OPEN: seven or more directions below their quadratic height cutoff; entry; synchronization; LRC14')


if __name__ == '__main__':
    main()
