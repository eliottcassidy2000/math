#!/usr/bin/env python3
"""Exact controls for the all-height one-ray projection closure.

The complete proof head is c<43; an independent sheet-interval path checks
every raw projection in that head. No repository mathematics is imported.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd


CHECKS = 0
TARGET = Q(6, 77)
SHARP = Q(12, 161)


def need(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)


def direction(v):
    g = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
    result = tuple(x//g for x in v)
    return result if next(x for x in result if x) > 0 else tuple(-x for x in result)


def carriers(w):
    a, b, c = w
    bounds = [(3*(sum(w)-v)-1)//14 for v in w]
    result = set()
    for x in range(-bounds[0], bounds[0]+1):
        if x % 3 == 0:
            continue
        for y in range(-bounds[1], bounds[1]+1):
            if y % 3 == 0 or (a*x+b*y) % c:
                continue
            z = -(a*x+b*y)//c
            if z % 3 and abs(z) <= bounds[2]:
                result.add((x, y, z))
    return result


def projection_data(w, live):
    q = Q(3, 7*w[2])
    totals = [Q(0)]*3
    physical = Q(0)
    for C in live:
        caps = []
        for i in range(3):
            j, k = [t for t in range(3) if t != i]
            cap = min(q, Q(3*(w[j]+w[k])-14*abs(C[i]), 14*w[j]*w[k]))
            need(cap > 0, ("strict raw roof", w, C, i))
            totals[i] += cap
            caps.append(cap)
        physical += min(caps)
    return tuple(totals), physical


def sheets(w):
    denominator = 42*w[0]*w[1]*w[2]
    result = {}
    for i, speed in enumerate(w):
        multiplier = denominator//(42*speed)
        for s in range(3):
            pieces = []
            for owner in range(speed):
                center = 14*((3*owner-speed*s) % (3*speed))*multiplier
                left, right = center-3*multiplier, center+3*multiplier
                if left < 0:
                    pieces.extend(((0, right), (denominator+left, denominator)))
                elif right > denominator:
                    pieces.extend(((left, denominator), (0, right-denominator)))
                else:
                    pieces.append((left, right))
            result[i, s] = sorted(pieces)
    return denominator, result


def intersections(first, second):
    i = j = 0
    result = []
    while i < len(first) and j < len(second):
        left, right = max(first[i][0], second[j][0]), min(first[i][1], second[j][1])
        if left < right:
            result.append((left, right))
        ir, jr = first[i][1], second[j][1]
        i += ir <= jr
        j += jr <= ir
    return result


def literal_projection_data(w):
    denominator, ds = sheets(w)
    totals = []
    masses = []
    for omitted in range(3):
        i, j = [v for v in range(3) if v != omitted]
        total = mass = 0
        for assignment in permutations(range(3)):
            pair = intersections(ds[i, assignment[i]], ds[j, assignment[j]])
            third = ds[omitted, assignment[omitted]]
            x = y = 0
            while x < len(pair) and y < len(third):
                l, r = max(pair[x][0], third[y][0]), min(pair[x][1], third[y][1])
                if l < r:
                    total += min(pair[x][1]-pair[x][0], third[y][1]-third[y][0])
                    mass += r-l
                xr, yr = pair[x][1], third[y][1]
                x += xr <= yr
                y += yr <= xr
        totals.append(Q(total, denominator))
        masses.append(Q(mass, denominator))
    need(len(set(masses)) == 1, ("three literal groupings", w))
    return tuple(totals), masses[0]


def audit_ray(w, live, projection):
    directions = {direction(C) for C in live}
    if len(directions) != 1:
        return None
    v = next(iter(directions))
    need(sum(abs(t) for t in v) % 2 == 0, ("odd-speed parity", w, v))
    need(all(t % 3 for t in v), ("primitive direction owner residue", w, v))
    B = min(Q(3*(sum(w)-w[i]), 14*abs(v[i])) for i in range(3))
    K = (B.numerator-1)//B.denominator
    expected = {tuple(sign*k*t for t in v) for k in range(1, K+1)
                if k % 3 for sign in (-1, 1)}
    need(live == expected, ("complete raw ray with strict endpoint", w, live ^ expected))
    count = 2*(K-K//3)
    need(count == len(live), ("ray count", w))
    need(count < Q(4, 3)*(B+1), ("residue count envelope", w, B, count))
    if sum(abs(t) for t in v) != 4:
        need(max(abs(t) for t in v) >= 4, ("max coefficient", w, v))
        need(B < Q(3*w[2], 28), ("ray cutoff", w, v, B))
        need(max(projection) <= Q(3*count, 7*w[2]), ("count pays every projection", w))
        need(max(projection) < Q(3, 49)+Q(4, 7*w[2]), ("all-height envelope", w))
        if w[2] >= 43:
            need(max(projection) < SHARP, ("sharp infinite tail", w))
    return v, B, count


def main():
    need(Q(3, 49)+Q(4, 7*43) < SHARP < TARGET, "tail thresholds")
    # The finite coefficient classification behind max |v_i|>=4.
    small = []
    for v in __import__('itertools').product((-2, -1, 1, 2), repeat=3):
        if gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2])) == 1 and sum(v) % 2 == 0:
            small.append(tuple(sorted(map(abs, v))))
    need(set(small) == {(1, 1, 2)}, "complete small direction universe")

    kinds = Counter()
    leader = (Q(0), None, None)
    equalities = []
    first_multi = None
    digest = sha256()
    eligible = [a for a in range(1, 43, 2) if a % 3]
    rows = sorted((w for w in combinations(eligible, 3)
                   if gcd(gcd(w[0], w[1]), w[2]) == 1), key=lambda w:(w[2],w[0],w[1]))
    for w in rows:
        live = carriers(w)
        projections, physical = projection_data(w, live)
        literal = literal_projection_data(w)
        need((projections, physical) == literal, ("raw/literal projection equality", w))
        ray = audit_ray(w, live, projections)
        if not live:
            kinds['empty'] += 1
            need(projections == (0,0,0), ("empty positive control", w))
        elif ray is None:
            kinds['multiple'] += 1
            if first_multi is None:
                v = min({direction(C) for C in live}, key=lambda C:(sum(map(abs,C)),C))
                B = min(Q(3*(sum(w)-w[i]),14*abs(v[i])) for i in range(3))
                first_multi = w, tuple(sorted(live)), v, B, Q(4,3)*(B+1), projections
        elif sum(map(abs,ray[0])) == 4:
            kinds['norm4'] += 1
            need(min(projections) <= TARGET, ("inherited norm4 control", w))
        else:
            kinds['non_norm4_ray'] += 1
            need(max(projections) <= SHARP, ("finite sharp head", w, projections))
            if max(projections) > leader[0]:
                leader = max(projections), w, projections
            if max(projections) == SHARP:
                equalities.append((w, tuple(i+1 for i,v in enumerate(projections) if v == SHARP)))
        digest.update(repr((w,tuple(sorted(live)),projections,physical)).encode())

    need(kinds['non_norm4_ray'] == 236, ("non-norm4 proof head", kinds))
    need(equalities == [((1,19,23),(2,3))], ("sharp equality universe", equalities))
    need(first_multi[0] == (17,23,25), ("smallest multi-ray height", first_multi))
    need(len(first_multi[1]) > first_multi[4], "hostile destroys one-ray counting implication")
    print('PROVED one-ray non-norm4 envelope: every E_i<3/49+4/(7c)')
    print('PROVED sharp every-projection cap:12/161; equality only(1,19,23),coordinates2,3')
    print('PROOF HEAD primitive odd ternary-unit triples c<43:', len(rows), dict(kinds))
    print('HEAD LEADER', leader)
    print('SMALLEST MULTI-RAY HOSTILE height then lex:', first_multi)

    controls = ((1,5,7),(1,5,11),(1,19,23),(17,23,25),(19,23,29),
                (1,1201,1205),(1,599,607),(5,1001,1003))
    for w in controls:
        live = carriers(w)
        result = projection_data(w, live)
        need(result == literal_projection_data(w), ("wide raw/literal", w))
        ray = audit_ray(w, live, result[0])
        print('CONTROL', w, 'N', len(live), 'ray', ray, 'E', tuple(map(str,result[0])), 'mass', str(result[1]))
    live = carriers((1,1201,1205))
    need(len(live) > 90, 'unbounded-count family control beyond ninety-carrier gate')
    print('SEMANTIC SHA256', digest.hexdigest())
    print('CHECKS', CHECKS)
    print('OPEN multi-ray target; entry; synchronization; LRC14')


if __name__ == '__main__':
    main()
