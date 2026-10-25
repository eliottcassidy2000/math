#!/usr/bin/env python3
"""Exact, standard-library controls. All gates remain active under python -O."""
import hashlib
import itertools as it
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def det(rows):
    a = [list(r) for r in rows]
    sign, prev = 1, 1
    for k in range(len(a)-1):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        v = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                num = v*a[i][j]-a[i][k]*a[k][j]
                if num % prev:
                    raise RuntimeError('nonexact Bareiss division')
                a[i][j] = num//prev
            a[i][k] = 0
        prev = v
    return sign*a[-1][-1]


def inverse(rows, modulus=None):
    n = len(rows)
    a = [[Fraction(x) for x in r] + [Fraction(i == j) for j in range(n)]
         for i, r in enumerate(rows)]
    if modulus:
        a = [[int(x) % modulus for x in r] for r in a]
    for k in range(n):
        i = next(i for i in range(k, n) if a[i][k])
        a[k], a[i] = a[i], a[k]
        v = pow(a[k][k], -1, modulus) if modulus else 1/a[k][k]
        a[k] = [x*v % modulus if modulus else x*v for x in a[k]]
        for i in range(n):
            if i != k:
                v = a[i][k]
                a[i] = [(x-v*y) % modulus if modulus else x-v*y
                        for x, y in zip(a[i], a[k])]
    return [r[n:] for r in a]


def product(a, b):
    return [[sum(x*y for x, y in zip(r, c)) for c in zip(*b)] for r in a]


def bits(v):
    return tuple((v >> i) & 1 for i in range(6))


def chars(labels):
    return frozenset(sum((((a & v).bit_count() % 2) << i)
                         for i, v in enumerate(labels[1:])) for a in range(1, 8))


def bases(code):
    return [tuple(bits(x) for x in t) for t in it.combinations(sorted(code), 3)
            if t[0] ^ t[1] ^ t[2]]


def rows(labels, masks):
    return [[(m & v).bit_count() % 2 for v in labels[1:]] for m in masks]


def matrix_repair(c):
    """Independent integral lift via the mod-two elimination word."""
    a = [[x % 2 for x in r] for r in c]
    n = len(c)
    ops = []
    for k in range(n):
        i = next(i for i in range(k, n) if a[i][k])
        if i != k:
            a[k], a[i] = a[i], a[k]
            ops.append(('swap', k, i))
        for i in range(n):
            if i != k and a[i][k]:
                a[i] = [x ^ y for x, y in zip(a[i], a[k])]
                ops.append(('add', i, k))
    lift = [[int(i == j) for j in range(n)] for i in range(n)]
    for op, i, j in reversed(ops):
        if op == 'swap':
            lift[i], lift[j] = lift[j], lift[i]
        else:
            lift[i] = [x-y for x, y in zip(lift[i], lift[j])]
    need(all((lift[i][j]-c[i][j]) % 2 == 0 for i in range(n) for j in range(n)),
         'marked-character lift parity')
    need(abs(det(lift)) == 1, 'marked-character lift unimodular')
    return lift


def monomial(vals, r):
    out = Fraction(1)
    for v, exponent in zip(vals, r):
        out *= v**exponent
    return out


def main():
    U = (1, 2, 4, 5, 6, 3, 7)
    D = (1, 2, 4, 6, 3, 7, 5)
    U4 = (1, 2, 4, 5, 7, 3, 6)
    D4 = (1, 2, 4, 7, 3, 6, 5)
    c0 = rows(U, (1, 2, 4)) + rows(D, (1, 2, 4))
    c1 = rows(U, (1, 2, 4)) + rows(D, (1, 2, 5))
    c4 = rows(U4, (1, 2, 4)) + rows(D4, (1, 2, 4))
    need([det(c) for c in (c0, c1, c4)] == [3, 1, -3], 'norm pair determinants')
    inv = inverse(c1)
    need(all(v.denominator == 1 for row in inv for v in row), 'integer inverse')
    expected = [[0,1,0,-1,0,0],[1,2,-2,-1,-2,3],[0,-1,1,0,1,-1],
                [-1,-1,1,1,1,-1],[1,1,-2,0,-1,2],[0,0,1,0,0,-1]]
    need(inv == expected, 'printed inverse')
    lifts = [matrix_repair(c) for c in (c0, c1, c4)]
    for c in (c0, c1, c4):
        b = inverse(c, 2)
        bc = product(b, c)
        h = [[(bc[i][j]-int(i == j))//2 for j in range(6)] for i in range(6)]
        need(all(bc[i][j] == int(i == j)+2*h[i][j] for i in range(6) for j in range(6)),
             'exact root decoder exponents')
        # Every basis monomial maps to a different parity, with a unit coefficient.
        images = set()
        for word in it.product(range(2), repeat=6):
            exponent = [sum(word[i]*c[i][j] for i in range(6)) for j in range(6)]
            images.add(tuple(v % 2 for v in exponent))
        need(len(images) == 64, 'universal algebra basis permutation')
        for signs in it.product((-1, 1), repeat=6):
            target = [int(monomial(signs, row)) for row in c]
            full = 1
            pair = 1
            for s, t in zip(signs, target):
                full *= 1+s
                pair *= 1+t
            need(full == pair, 'common-source quadratic character count')
        for signs in it.product((-1, 1), repeat=6):
            u = [Fraction(s*(i+2), i+1) for i, s in enumerate(signs)]
            rad = [v*v for v in u]
            w = [monomial(u, row) for row in c]
            recovered = [monomial(w, row)/monomial(rad, hrow) for row, hrow in zip(b, h)]
            need(recovered == u, 'rational signed-root recovery')
    u = [Fraction(i+2, i+1) for i in range(6)]
    w0 = [monomial(u, row) for row in c0]
    w1 = [monomial(u, row) for row in c1]
    need(w1[5] == w0[3]*w0[5]/(u[4]*u[5])**2, 'base correction K4')
    need([monomial(w1, row) for row in inv] == u, 'Laurent inverse without radicands')
    # Dropping the radicands really loses information: old torus map over F_7.
    kernel3 = [v for v in it.product(range(3), repeat=6)
               if all(sum(x*y for x,y in zip(row,v)) % 3 == 0 for row in c0)]
    need(len(kernel3) == 3, 'old mu3 kernel')
    need(all(all(sum(row[j] for j in range(6) if j != axis) >= 1 for row in c1)
             for axis in range(6)), 'affine axes collapse')
    code_set = sorted({chars((1,)+perm) for perm in it.permutations((2,3,4,5,6,7))},
                      key=lambda x: tuple(sorted(x)))
    need(len(code_set) == 30, 'all labelled Fano structures')
    by_code = {c: bases(c) for c in code_set}
    need(all(len(b) == 28 for b in by_code.values()), 'all unordered character bases')
    pairs = [(a,b) for a,b in it.combinations(code_set, 2) if not a.intersection(b)]
    need(len(pairs) == 120, 'all complementary Fano pairs')
    target = Counter({1:462,3:231,5:63,7:21,9:7})
    total = Counter()
    representatives = {}
    for a, b in pairs:
        distribution = Counter()
        for ba in by_code[a]:
            for bb in by_code[b]:
                c = ba+bb
                value = abs(det(c))
                distribution[value] += 1
                representatives.setdefault(value, c)
        need(distribution == target, 'Boolean lift census at fixed labelled branch zero')
        total.update(distribution)
    for c in representatives.values():
        matrix_repair(c)
    # Six/seven comparison, with a hostile showing seven is not the first case.
    from math import gcd
    eligible = {m: [a for a in range(1,m) if gcd(a,m)==1 and a in (1,m-1)]
                for m in (5,6,7)}
    need(eligible == {5:[1,4],6:[1,5],7:[1,6]}, 'marked scalar integral lifts')
    result = {'status':'FINITE-EXACT','checks':CHECKS,'Fano_codes':30,
              'complementary_pairs':120,'Boolean_matrices':sum(total.values()),
              'per_pair_abs_det':dict(sorted(target.items())),
              'all_pairs_abs_det':dict(sorted(total.items())),
              'C0':c0,'C1':c1,'C1_inverse':expected,'C4':c4,
              'same_character_unimodular_lifts':lifts,
              'old_mod3_kernel_exponents':kernel3,
              'determinant_representatives':{k:list(v) for k,v in sorted(representatives.items())},
              'scope':'All 120 labelled Fano pairs, 784 unordered Boolean bases each; all 64 sign patterns for three displayed matrices. No point-existence or whole-surface F1 claim.',
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
