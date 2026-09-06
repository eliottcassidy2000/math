#!/usr/bin/env python3
"""Independent centered-division audit of coprime-maximum actual decoder entry.

Standard library only; no mathematical producer imports and no height census.
Incoming overnight12 gcd_semigroup Corollary4 strictly subsumes this closure.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd, prod
import json

Q = 91**6
GATES = 0


def check(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def centered_division(x, u, K, bound):
    if not (1 <= u < K <= bound and gcd(u, K) == 1):
        raise ValueError('coprime anchor and bounded maximum required')
    if 2*abs(x) > bound*K:
        raise ValueError('outside the proved sufficient interval')
    a = (x*pow(u, -1, K)) % K
    if 2*a > K:
        a -= K
    b = (x-u*a)//K
    return a, b


def inert_pair(p, q):
    if not (1 <= p < q and gcd(p, q) == 1 and p+q <= 356):
        return False
    n, prime = p+q, 2
    while prime*prime <= n:
        if n % prime == 0:
            exponent = 0
            while n % prime == 0:
                n //= prime
                exponent += 1
            if prime % 3 != 2 or exponent > 2:
                return False
        prime += 1
    return n == 1 or n % 3 == 2


def norm(v, x):
    residue = (v*x.numerator) % x.denominator
    return F(min(residue, x.denominator-residue), x.denominator)


def ceiling(x):
    return -((-x.numerator)//x.denominator)


def grid_witness(V, u, t, g, p, q, y):
    K = max(V)
    check(len(V) == len(set(V)) == 11 and min(V) > 1,
          'nonunit eleven-shape witness control')
    check(u in V and u < K and gcd(u, K) == 1,
          'actual coprime-to-maximum anchor')
    check(gcd(t, g) == 1 and inert_pair(p, q), 'scales and pair domain')
    check(min(norm(v, y) for v in V) >= F(1, 12), 'literal core supplier')
    n = tuple(t*v for v in V)+(g*p, g*q)
    check(len(set(n)) == 13 and gcd(*n) == 1 and sum(n) <= Q*Q,
          'physical control type and finite box')
    order = p+q
    z = F((pow(p, -1, order)*(order//2)) % order, order)
    check(min(norm(p, z), norm(q, z)) >= F(1, 3), 'explicit pair phase')
    if 11*t >= 21*q:
        left, right = z-F(11, 42*q), z+F(11, 42*q)
        check(right-left >= F(1, t), 'closed pair arc meets t-grid')
        k = ceiling(t*left-g*y)
        j = (pow(g, -1, t)*k) % t
        x = (y+j)/t
        address = (g*y+k)/t
        check(left <= address <= right and (g*x-address).denominator == 1,
              'coherent pair-clock lift')
        branch = 'large core scale'
    else:
        delta = gcd(t, p)
        check(t <= 677 < Q and 2*g*p > delta*Q*K,
              'exact surviving small-scale inequality')
        check(g > 42*K, 'forced core-grid density')
        left, right = y-F(1, 84*K), y+F(1, 84*K)
        k = ceiling(g*left-t*z)
        if (t*z+k)/g == left:
            k += 1
        j = (pow(t, -1, g)*k) % g
        x = (z+j)/g
        address = (t*z+k)/g
        check(left < address < right and (t*x-address).denominator == 1,
              'coherent core-clock lift')
        branch = 'forced large pair scale'
    clearance = min(norm(v, x) for v in n)
    check(clearance >= F(1, 14), 'direct full-row clearance')
    return dict(core=list(V), anchor=u, t=t, g=g, p=p, q=q,
                branch=branch, phase=str(x), clearance=str(clearance))


def main():
    division_rows = 0
    for bound in range(2, 13):
        for K in range(2, bound+1):
            for u in range(1, K):
                if gcd(u, K) != 1:
                    continue
                radius = bound*K//2
                for x in range(-radius, radius+1):
                    a, b = centered_division(x, u, K, bound)
                    check(x == a*u+b*K and 2*abs(a) <= K
                          and abs(b) <= bound, 'complete centered integer interval')
                    division_rows += 1
    # Named controls retain inclusive endpoints, large exact integers and signs.
    for K in (11, 13, Q):
        u = K-1
        for x in (-(Q*K//2), -1, 0, 1, Q*K//2):
            a, b = centered_division(x, u, K, Q)
            check(a*u+b*K == x and max(abs(a), abs(b)) <= Q,
                  'large exact compiler endpoints')
    check(all(2*a+4*b != 1 for a in range(-4, 5) for b in range(-4, 5)),
          'dropping coprimality loses an integer residue')
    # The radius is sufficient, not sharp: this point just outside still works.
    check(2*7 > 4*3 and 7 == 2*2+1*3, 'unclaimed sharpness control')

    # Independent literal controls for the stronger incoming signed-box lemma.
    box_rows = 0
    multiplier_rows = 0
    for bound in range(2, 10):
        for K in range(2, bound+1):
            for u in range(1, K):
                if gcd(u, K) != 1:
                    continue
                values = {u*a+K*b for a in range(-bound, bound+1)
                          for b in range(-bound, bound+1)}
                radius = bound*(u+K)-(u-1)*(K-1)
                check(all(x in values for x in range(-radius, radius+1))
                      and radius+1 not in values and -radius-1 not in values,
                      'incoming exact centered radius by full coefficient box')
                check(radius >= bound*K and 2*radius > bound*(u+K),
                      'incoming radius dominates our sufficient half-radius')
                box_rows += 1
                if bound > 6:
                    continue
                for scale in (1, 2, 5):
                    physical_values = {scale*v for v in values}
                    for outside in range(1, 21):
                        delta = gcd(scale, outside)
                        minimal = scale//delta
                        x = outside//delta
                        exists = any(c*outside in physical_values
                                     for c in range(1, bound+1))
                        check(exists == (minimal <= bound and x in values),
                              'incoming minimal outside multiplier iff')
                        multiplier_rows += 1
    outside_box = {a+6*b for a in range(-2, 3) for b in range(-2, 3)}
    check(3 not in outside_box and 6 in outside_box,
          'incoming minimal-multiplier lemma needs reduced height bound')
    check(all(c not in {3*a+6*b for a in range(-2, 3) for b in range(-2, 3)}
              for c in (1, 2)), 'outside coefficient gate cannot be dropped')

    crossing_rows = 0
    for bound in (8, 13, 21):
        for u, K in ((2, 3), (2, 5), (3, 5), (5, 8), (7, 11)):
            if K > bound:
                continue
            for t in (1, 2, 3):
                for p in (1, 2, 3, 5):
                    for g in range(1, 32):
                        if gcd(t, g) != 1:
                            continue
                        delta = gcd(t, p)
                        x = g*p//delta
                        if 2*x > bound*K:
                            continue
                        a, b = centered_division(x, u, K, bound)
                        c = t//delta
                        check(c*g*p-a*t*u-b*t*K == 0
                              and max(c, abs(a), abs(b)) <= bound
                              and c*g*p != 0, 'literal bounded crossing relation')
                        crossing_rows += 1

    pairs = [(p, q) for q in range(2, 356)
             for p in range(1, min(q, 357-q)) if inert_pair(p, q)]
    check(len(pairs) == 5855 and max(p for p, q in pairs) == 177
          and max(q for p, q in pairs) == 355, 'complete inherited pair atlas')
    for p, q in pairs:
        check(Q > 84*p and (21*q-1)//11 <= 677 < Q,
              'all-atlas crossing-to-gluing constants')
    endpoint_cutoff = Q//(42*177)
    check(endpoint_cutoff == 76388115 and 677*endpoint_cutoff < Q
          and 42*177*endpoint_cutoff <= Q,
          'strictly stronger incoming endpoint-gcd cutoff')

    # A genuine actual-entry control outside the literal-unit theorem.
    V = tuple(range(2, 12))+(13,)
    g = 2**45
    n = V+(g, 3*g)
    parent = list(range(13))

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i

    for i, j in combinations(range(13), 2):
        d = gcd(n[i], n[j])
        p, q = sorted((n[i]//d, n[j]//d))
        if inert_pair(p, q):
            parent[find(i)] = find(j)
    components = sorted((sorted(i for i in range(13) if find(i) == r)
                         for r in {find(i) for i in range(13)}), key=len)
    check(components == [[11, 12], list(range(11))], 'actual decoder graph is 2+11')
    check(sum(n) <= Q*Q and g > 2*Q*max(V),
          'actual finite box and no-crossing dominance')
    # Any crossing row has one or two core terms, of size at most2QK<g;
    # its nonzero pair partial sum is a nonzero multiple of g. Thus no
    # support<=3, height<=Q crossing relation exists. Connected internal
    # decoder rows span all internal relations, proving W=V_dec analytically.
    actual = grid_witness(V, 2, 1, g, 1, 3, F(1, 12))

    # Canonical primitive nonunit hostile to replacing the anchor by gcd(V)=1.
    primes = (37, 43, 61, 67, 73, 79, 97, 103, 127)
    P = 15*prod(primes)
    hostile = tuple(2*P//r for r in primes)+(P//3, P//5)
    K = max(hostile)
    check(gcd(*hostile) == 1 and min(hostile) > 1 and K > Q,
          'primitive nonunit maximum need not be bounded')
    check(all(gcd(u, K) > 1 for u in hostile if u != K),
          'canonical hostile lies outside coprime-maximum hypothesis')
    check(max(max(a, b)//gcd(a, b) for a, b in combinations(hostile, 2)) == 127,
          'hostile retains every inherited internal pair-height bound')
    check((gcd(*hostile[:7]), gcd(*hostile[:8]), gcd(*hostile[:9]),
           gcd(*hostile[:9], hostile[9])) == (392430, 3810, 30, 5),
          'normalization hostile already violates inherited large-subset caps')

    witnesses = []
    for core, anchor, y in ((V, 2, F(1, 12)),
                            (tuple(range(3, 24, 2)), 3, F(1, 2))):
        for p, q in ((1, 3), (1, 4), (1, 355), (177, 178)):
            t = (21*q+10)//11
            witnesses.append(grid_witness(core, anchor, t, 1, p, q, y))
            for t in (1, 2, 17, 677):
                if 11*t >= 21*q:
                    continue
                delta = gcd(t, p)
                g = delta*Q*max(core)//(2*p)+1
                while gcd(t, g) != 1:
                    g += 1
                witnesses.append(grid_witness(core, anchor, t, g, p, q, y))

    manifest = dict(Q=Q, division_rows=division_rows, crossing_rows=crossing_rows,
                    incoming_exact_box_rows=box_rows,
                    incoming_minimal_multiplier_rows=multiplier_rows,
                    incoming_endpoint_cutoff=endpoint_cutoff,
                    atlas_rows=len(pairs), actual_entry=actual,
                    generic_grid_witnesses=witnesses, hostile_maximum=K,
                    components=components)
    raw = json.dumps(manifest, sort_keys=True, separators=(',', ':')).encode()
    print('STATUS: independent centered-division audit; coprime specialization already')
    print('        subsumed by incoming overnight12 gcd_semigroup Corollary4; LRC14 OPEN')
    print('DOMAIN: full W=V_dec, actual 11+2 decoder split and finite box;')
    print('        inert pair atlas, primitive scales, and a label coprime to maxV')
    print('centered_division_rows', division_rows, 'cleared_crossing_rows', crossing_rows)
    print('incoming_exact_box_rows', box_rows, 'incoming_minimal_multiplier_rows',
          multiplier_rows, 'incoming_endpoint_cutoff', endpoint_cutoff)
    print('inherited_pair_atlas', len(pairs), 'generic_grid_witnesses', len(witnesses))
    print('actual_nonunit_entry', json.dumps(actual, sort_keys=True))
    print('nonunit_hostile_maximum', K, 'coprime_maximum_anchors', 0)
    print('normalization_hostile_already_excluded_by_large_subset_caps', True)
    print('generic_grid_branches', sorted({r['branch'] for r in witnesses}))
    print('gates', GATES)
    print('semantic_sha256', sha256(raw).hexdigest())


if __name__ == '__main__':
    main()
