#!/usr/bin/env python3
"""Projective two-jet interpolation: exact inverse, Smith, and tree audits.

Finite universe is printed below. Standard-library arithmetic only. Independent
paths: polynomial cardinal columns, rational Gauss-Jordan inverse, rational DVR
Smith elimination, and Bareiss determinants/all minors in small cases. No
optimization-disabled assertions or floating-point arithmetic.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm, prod
import json
import random
import sys
sys.stdout.reconfigure(newline='\n')
GATES = 0


def require(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def det(v, w):
    return v[0]*w[1]-v[1]*w[0]


def tangent(v):
    # Extended Euclid, including signed primitive coordinates.
    a, b = v
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r//r
        old_r, r, old_s, s, old_t, t = r, old_r-q*r, s, old_s-q*s, t, old_t-q*t
    require(abs(old_r) == 1, 'primitive direction')
    w = (-old_t*old_r, old_s*old_r)
    require(det(v, w) == 1, 'unimodular tangent')
    return w


def vp(x, p):
    if x == 0:
        return None
    x = Fraction(x)
    a, b, ans = abs(x.numerator), x.denominator, 0
    while a % p == 0:
        a //= p
        ans += 1
    while b % p == 0:
        b //= p
        ans -= 1
    return ans


def observer(vectors, shifts=None):
    n = len(vectors)
    degree = 2*n-1
    shifts = shifts or [0]*n
    matrix = []
    for v, shift in zip(vectors, shifts):
        w0 = tangent(v)
        w = tuple(w0[k]+shift*v[k] for k in (0, 1))
        a, b = v
        c, d = w
        matrix.append([a**q*b**(degree-q) for q in range(degree+1)])
        matrix.append([(q*c*a**(q-1)*b**(degree-q) if q else 0)
                       + ((degree-q)*d*a**q*b**(degree-q-1) if q < degree else 0)
                       for q in range(degree+1)])
    return matrix


def multiply(a, b):
    ans = [Fraction(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            ans[i+j] += x*y
    return ans


def cardinal(vectors):
    columns, largest, local = [], 1, []
    for i, v in enumerate(vectors):
        w = tangent(v)
        differences = [det(v, u) for j, u in enumerate(vectors) if j != i]
        d = prod(differences)
        c = 2*sum(det(w, u)*prod(det(v, z) for k, z in enumerate(vectors)
                                if k not in (i, j))
                  for j, u in enumerate(vectors) if j != i)
        q = [Fraction(1)]
        for j, u in enumerate(vectors):
            if j != i:
                line = [-u[0], u[1]]  # det((X,Y),u), ordered by X degree
                q = multiply(multiply(q, line), line)
        t, s = [v[0], -v[1]], [-w[0], w[1]]
        value = multiply(q, [Fraction(d*s[k]-c*t[k], d**3) for k in (0, 1)])
        derivative = multiply(q, [Fraction(x, d*d) for x in t])
        columns.extend((value, derivative))
        denominator = abs(d)**3//gcd(abs(d), c)
        require(lcm(*(x.denominator for x in value)) == denominator,
                'value column exact denominator')
        require(lcm(*(x.denominator for x in derivative)) == d*d,
                'derivative column exact denominator')
        largest = lcm(largest, denominator)
        local.append((d, c))
    return list(map(list, zip(*columns))), largest, local


def inverse(matrix):
    n = len(matrix)
    a = [[Fraction(x) for x in row]+[Fraction(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for k in range(n):
        i = next(i for i in range(k, n) if a[i][k])
        a[k], a[i] = a[i], a[k]
        pivot = a[k][k]
        a[k] = [x/pivot for x in a[k]]
        for i in range(n):
            if i != k:
                factor = a[i][k]
                a[i] = [x-factor*y for x, y in zip(a[i], a[k])]
    return [row[n:] for row in a]


def bareiss(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    sign, previous = 1, 1
    if n == 0:
        return 1
    for k in range(n-1):
        i = next((i for i in range(k, n) if a[i][k]), None)
        if i is None:
            return 0
        if i != k:
            a[k], a[i] = a[i], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = pivot*a[i][j]-a[i][k]*a[k][j]
                require(numerator % previous == 0, 'Bareiss exact division')
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def smith(matrix, p):
    a = [[Fraction(x) for x in row] for row in matrix]
    n, exponents = len(a), []
    for k in range(n):
        i, j = min(((i, j) for i in range(k, n) for j in range(k, n) if a[i][j]),
                   key=lambda ij: vp(a[ij[0]][ij[1]], p))
        a[k], a[i] = a[i], a[k]
        for row in a:
            row[k], row[j] = row[j], row[k]
        pivot = a[k][k]
        exponents.append(vp(pivot, p))
        for i in range(k+1, n):
            factor = a[i][k]/pivot
            require(not factor or vp(factor, p) >= 0, 'DVR integral elimination')
            for j in range(k, n):
                a[i][j] -= factor*a[k][j]
    require(exponents == sorted(exponents), 'Smith ordering')
    return tuple(exponents)


def tree_loss(vectors, p):
    n = len(vectors)
    if n == 1:
        return 0
    h = [[None if i == j else vp(det(v, u), p)
          for j, u in enumerate(vectors)] for i, v in enumerate(vectors)]
    scores = [sum(x for x in row if x is not None) for row in h]
    result = 0
    for i in range(n):
        for f in sorted(set(x for x in h[i] if x is not None)):
            cluster = [j for j in range(n) if i == j or h[i][j] >= f]
            if all(h[j][k] == f for j, k in combinations(cluster, 2)):
                require(len({scores[j] for j in cluster}) == 1, 'terminal mass constant')
                require(len(cluster) <= (p+1 if f == 0 else p), 'projective child count')
                result = max(result, 2*scores[i]+max(0, f-(len(cluster) == p)))
    return result


def apply(g, v):
    return (g[0][0]*v[0]+g[0][1]*v[1], g[1][0]*v[0]+g[1][1]*v[1])


def run():
    cases = [((1, 0),), ((1, 0), (0, 1)),
             ((1, 0), (0, 1), (1, 1)),
             ((1, 0), (1, 2), (1, 4)),
             ((0, 1), (8, 1), (16, 1), (40, 1)),
             ((0, 1), (8, 1), (24, 1), (32, 1))]
    for p in (2, 3, 5, 7):
        cases.append(tuple([(1, 0)]+[(a, 1) for a in range(p)]))
        for e in (1, 2, 4):
            # A deep cluster at infinity plus every finite residue class.
            cases.append(tuple([(1, 0), (1, p**e)]+[(a, 1) for a in range(p)]))
    pool = [(a, b) for b in range(0, 5) for a in range(-4, 5)
            if gcd(a, b) == 1 and (b > 0 or a > 0)]
    rng = random.Random(20260906)
    for n in (2, 3, 4, 5):
        for _ in range(12):
            cases.append(tuple(rng.sample(pool, n)))
    primes = (2, 3, 5, 7)
    records = []
    for vectors in cases:
        a = observer(vectors)
        expected, largest, local = cardinal(vectors)
        require(inverse(a) == expected, 'cardinal inverse equals independent Gauss Jordan')
        require(lcm(*(x.denominator for row in inverse(a) for x in row)) == largest,
                'global largest denominator')
        determinant = prod(det(v, u)**4 for v, u in combinations(vectors, 2))
        require(abs(bareiss(a)) == determinant, 'projective confluent determinant')
        spectra = []
        for p in primes:
            alpha = smith(a, p)
            require(sum(alpha) == vp(determinant, p), 'determinant mass')
            require(alpha[-1] == vp(largest, p) == tree_loss(vectors, p), 'sharp projective tree precision')
            spectra.append(alpha)
        records.append([vectors, largest, spectra])
    transformations = (((0, -1), (1, 0)), ((1, 2), (1, 3)), ((1, 0), (0, -1)))
    symmetry_rows = 0
    for vectors in cases[:18]+cases[-12:]:
        old = observer(vectors)
        for g in transformations:
            newvectors = tuple(apply(g, v) for v in vectors)
            shifted = observer(newvectors, [(-1)**i*(i+1) for i in range(len(vectors))])
            for p in primes:
                require(smith(old, p) == smith(shifted, p), 'full Smith projective and tangent covariance')
                symmetry_rows += 1
    all_minor_rows = 0
    for vectors in [((1, 0), (1, 4)), ((1, 0), (1, 2), (1, 4))]:
        a = observer(vectors)
        n = len(a)
        divisors = [1]
        for k in range(1, n+1):
            divisor = 0
            for rows in combinations(range(n), k):
                for cols in combinations(range(n), k):
                    divisor = gcd(divisor, bareiss([[a[i][j] for j in cols] for i in rows]))
                    all_minor_rows += 1
            divisors.append(divisor)
        for p in primes:
            direct = tuple(vp(divisors[k], p)-vp(divisors[k-1], p) for k in range(1, n+1))
            require(direct == smith(a, p), 'independent all-minor full Smith')
    twin_a, twin_b = cases[4], cases[5]
    alpha_a, alpha_b = smith(observer(twin_a), 2), smith(observer(twin_b), 2)
    require(alpha_a == (0, 0, 4, 7, 12, 16, 19, 26), 'inherited twin A')
    require(alpha_b == (0, 0, 4, 7, 12, 17, 18, 26), 'inherited twin B')
    # Same projective directions without primitive normalization change data.
    primitive = observer(((1, 0), (0, 1)))
    nonprimitive = [row[:] for row in primitive]
    nonprimitive[0] = [8*x for x in nonprimitive[0]]
    nonprimitive[1] = [4*x for x in nonprimitive[1]]
    require(smith(primitive, 2) == (0, 0, 0, 0), 'primitive basis control')
    require(smith(nonprimitive, 2) == (0, 0, 2, 3), 'nonprimitive scaling hostile')
    print('PROVED SCOPE: primitive integer directions, complete two Hasse jets, homogeneous degree 2n-1.')
    print('Finite universe:', len(cases), 'direction configurations; primes', primes)
    print('Seeded pool: primitive (a,b), -4<=a<=4, 0<=b<=4, canonical sign; 12 samples at each n=2,3,4,5; seed=20260906.')
    print('Structured controls: n=1, infinity, all p+1 residue classes, depth 1/2/4 infinity cluster at p=2,3,5,7.')
    print('Independent inverse/determinant/tree/Smith comparisons:', len(cases)*len(primes))
    print('Full projective/tangent Smith symmetry rows:', symmetry_rows)
    print('Exhaustive nonempty minors:', all_minor_rows)
    print('Dyadic full-partition hostile A:', alpha_a)
    print('Dyadic full-partition hostile B:', alpha_b)
    print('Nonprimitive same-direction hostile: (0,0,0,0) versus (0,0,2,3).')
    print('semantic_sha256:', sha256(json.dumps(records, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
    print('Exact gates:', GATES)
    print('PASS')


if __name__ == '__main__':
    run()
