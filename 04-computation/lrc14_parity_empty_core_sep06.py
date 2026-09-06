#!/usr/bin/env python3
"""Bounded exact H79 parity probe; no all-height mixed-parity claim.

The primary path builds all literal sheet intervals and contact capacities.
An independent integer-relation path checks every projection and physical
mass. The universe includes all odd rows as inherited positive controls.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd

TARGET = Q(6, 77)
CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def intersect(A, B):
    answer = []
    i = j = 0
    while i < len(A) and j < len(B):
        left, right = max(A[i][0], B[j][0]), min(A[i][1], B[j][1])
        if left < right:
            answer.append((left, right))
        ar, br = A[i][1], B[j][1]
        i += ar <= br
        j += br <= ar
    return answer


def native(w):
    denominator = 42*w[0]*w[1]*w[2]
    sheet = {}
    for i, speed in enumerate(w):
        unit = denominator//(42*speed)
        for s in range(3):
            intervals = []
            for k in range(speed):
                center = 14*((3*k-speed*s) % (3*speed))*unit
                left, right = center-3*unit, center+3*unit
                if left < 0:
                    intervals.extend(((0, right), (denominator+left, denominator)))
                elif right > denominator:
                    intervals.extend(((left, denominator), (0, right-denominator)))
                else:
                    intervals.append((left, right))
            sheet[i, s] = sorted(intervals)
    projections, masses = [], []
    for omitted in range(3):
        i, j = [s for s in range(3) if s != omitted]
        capacity = physical = 0
        for pi in permutations(range(3)):
            pair = intersect(sheet[i, pi[i]], sheet[j, pi[j]])
            third = sheet[omitted, pi[omitted]]
            x = y = 0
            while x < len(pair) and y < len(third):
                left, right = max(pair[x][0], third[y][0]), min(pair[x][1], third[y][1])
                if left < right:
                    capacity += min(pair[x][1]-pair[x][0], third[y][1]-third[y][0])
                    physical += right-left
                xr, yr = pair[x][1], third[y][1]
                x += xr <= yr
                y += yr <= xr
        projections.append(Q(capacity, denominator))
        masses.append(Q(physical, denominator))
    need(len(set(masses)) == 1, ('three native pairings', w))
    return tuple(projections), masses[0]


def raw(w):
    a, b, c = w
    bounds = [(3*(sum(w)-wi)-1)//14 for wi in w]
    d = gcd(b, c)
    modulus = c//d
    inverse = pow(b//d, -1, modulus) if modulus > 1 else 0
    live = set()
    for x in range(-bounds[0], bounds[0]+1):
        if x % 3 == 0 or a*x % d:
            continue
        residue = (-a*x//d*inverse) % modulus if modulus > 1 else 0
        low = -bounds[1]+(residue+bounds[1]) % modulus
        for y in range(low, bounds[1]+1, modulus):
            if y % 3 == 0:
                continue
            z = -(a*x+b*y)//c
            if z % 3 and abs(z) <= bounds[2]:
                live.add((x, y, z))
    E = [Q(0)]*3
    mass = Q(0)
    for C in live:
        terms = []
        for i in range(3):
            j, k = [s for s in range(3) if s != i]
            value = min(Q(3, 7*c), Q(3*(w[j]+w[k])-14*abs(C[i]), 14*w[j]*w[k]))
            need(value > 0, ('raw strict roof', w, C, i))
            terms.append(value)
            E[i] += value
        mass += min(terms)
    return tuple(E), mass, live


def direction(C):
    divisor = gcd(gcd(abs(C[0]), abs(C[1])), abs(C[2]))
    v = tuple(x//divisor for x in C)
    return v if v[0] > 0 else tuple(-x for x in v)


def main():
    eligible = [s for s in range(1, 80) if s % 3]
    rows = sorted((w for w in combinations(eligible, 3)
                   if gcd(gcd(w[0], w[1]), w[2]) == 1),
                  key=lambda w:(w[2], w[0], w[1]))
    counts, failures, mechanisms = Counter(), Counter(), Counter()
    leaders = {k: {'physical': (Q(-1), []), 'network': (Q(-1), [])} for k in range(3)}
    first_failure = None
    other_failure = []
    norm5 = (Q(0), None)
    frozen = {}
    digest = sha256()
    for w in rows:
        even = sum(s % 2 == 0 for s in w)
        need(even < 3, ('primitive parity universe', w))
        counts[even] += 1
        E, physical = native(w)
        raw_E, raw_mass, C = raw(w)
        need((E, physical) == (raw_E, raw_mass), ('native versus raw', w))
        need(physical <= min(E), ('actual projection consumer', w))
        for name, value in (('physical', physical), ('network', min(E))):
            old, leaders_at = leaders[even][name]
            if value > old:
                leaders[even][name] = value, [w]
            elif value == old:
                leaders_at.append(w)
        directions = {direction(v) for v in C}
        if len(directions) == 1 and sum(map(abs, next(iter(directions)))) == 5 and min(E) > norm5[0]:
            norm5 = min(E), (w, E, physical, tuple(sorted(C)))
        if min(E) > TARGET:
            failures[even] += 1
            need(even != 0, ('inherited odd positive control', w))
            key = tuple(sorted({sum(map(abs, v)) for v in directions}))
            mechanisms[key] += 1
            record = w, E, physical, tuple(sorted(C))
            if first_failure is None:
                first_failure = record
            if key != (3,):
                other_failure.append(record)
        if w in ((1, 5, 11), (2, 5, 7), (1, 10, 11), (2, 11, 20), (2, 19, 20)):
            frozen[w] = E, physical, tuple(sorted(C))
        digest.update(repr((w, E, physical, tuple(sorted(C)))).encode())
    need(counts == {0:2910, 1:9044, 2:8694}, ('complete parity universe counts', counts))
    need(first_failure[0] == (2, 5, 7), ('minimal mixed hostile', first_failure))
    need(mechanisms == {(3,):243, (4,):1}, ('all failed row mechanisms', mechanisms))
    need(other_failure[0][0] == (2, 11, 20) and len(other_failure) == 1,
         ('sole nonadditive failure in finite head', other_failure))
    need(leaders[1]['network'] == (Q(6, 55), [(1, 10, 11)]), 'one-even finite leader')
    need(leaders[2]['network'] == (Q(11, 140), [(2, 11, 20)]), 'two-even finite leader')
    print('universe primitive sorted distinct positive 3-units, max<=79; no parity prefilter')
    print('row_counts_by_even_count', dict(sorted(counts.items())))
    print('network_failures_by_even_count', dict(sorted(failures.items())))
    print('failure_live_direction_norm_sets', dict(sorted(mechanisms.items())))
    print('first_failure', first_failure)
    print('sole_non_norm3_failure', other_failure[0])
    print('finite_norm5_network_leader', norm5)
    for parity in range(3):
        print('parity_even_count', parity, 'finite_leaders', leaders[parity])
    for w, values in sorted(frozen.items()):
        print('control', w, values)
    print('complete_semantic_digest', digest.hexdigest())
    print('explicit_checks', CHECKS)
    print('PASS: mixed-parity 6/77 is refuted; H79 parity maxima are FINITE-EXACT only')


if __name__ == '__main__':
    main()
