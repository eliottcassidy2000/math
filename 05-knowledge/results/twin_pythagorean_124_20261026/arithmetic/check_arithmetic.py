#!/usr/bin/env python3
"""Exact finite controls. No third-party packages or assertion dependence."""
from collections import defaultdict
from itertools import product
from math import gcd, isqrt, prod
import json

gates = 0


def need(condition, label):
    global gates
    gates += 1
    if not condition:
        raise RuntimeError(label)


def factors(n):
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d, e))
        d += 1
    if n > 1:
        out.append((n, 1))
    return out


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def rad(n):
    return prod(p for p, _ in factors(n))


def omega(n):
    return sum(e for _, e in factors(n))


def aliquot(n):
    return prod(sum(p ** j for j in range(e + 1))
                for p, e in factors(n)) - n


def literal_aliquot(n):
    result = 0
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            q = n // d
            if d < n:
                result += d
            if q != d and q < n:
                result += q
    return result


def triple(m, n):
    return (m * m - n * n, 2 * m * n, m * m + n * n)


def radii(T):
    A, B, C = T
    h = (A + B + C) // 2
    area = A * B // 2
    need(area % h == 0 and area % (h - B) == 0, 'integral radii')
    return area // h, area // (h - B)


def abc_check(a, b, c):
    need(a > 0 and b > 0 and c > 0 and a + b == c, 'abc sum')
    need(gcd(a, b) == gcd(a, c) == gcd(b, c) == 1, 'abc primitive')


euclid = set()
by_radius = defaultdict(list)
prime_radii = []
for m in range(2, 201):
    for n in range(1, m):
        if gcd(m, n) != 1 or (m - n) % 2 != 1:
            continue
        T = triple(m, n)
        A, B, C = T
        euclid.add(T)
        need(A * A + B * B == C * C, 'Euclid null cone')
        need(A % 2 == 1 and B % 2 == 0 and gcd(A, B) == 1, 'primitive parity')
        r, rho = radii(T)
        need((r, rho) == (n * (m - n), n * (m + n)), 'radius factorization')
        need(rho - r == 2 * n * n == C - A and rho + r == B, 'radius gap')
        nn = isqrt((rho - r) // 2)
        mm = (rho + r) // (2 * nn)
        need((mm, nn) == (m, n), 'radius inverse parameters')
        need(r * rho % (n * n) == 0 and
             (r * rho // (n * n), rho + r, r * rho // (n * n) + rho - r) == T,
             'radius inverse sides')
        predicates = [prime(r) and prime(rho),
                      n == 1 and prime(m - 1) and prime(m + 1),
                      C - A == 2 and omega(A) == 2]
        need(len(set(predicates)) == 1, 'three-way twin characterization')
        if predicates[0]:
            prime_radii.append([r, rho, list(T)])
        if prime(r):
            by_radius[r].append(T)

# A second, side-based universe does not call the Euclid generator.
direct = set()
for C in range(5, 1001, 2):
    for A in range(1, C, 2):
        bb = C * C - A * A
        B = isqrt(bb)
        if B > 0 and B % 2 == 0 and B * B == bb and gcd(A, B) == 1:
            direct.add((A, B, C))
need(direct == {T for T in euclid if T[2] <= 1000}, 'direct side universe equals Euclid')
for T in sorted(direct):
    A, B, C = T
    m, n = isqrt((C + A) // 2), isqrt((C - A) // 2)
    need(triple(m, n) == T, 'direct side inverse')
    need((prime(radii(T)[0]) and prime(radii(T)[1])) ==
         (C - A == 2 and omega(A) == 2), 'side-based twin characterization')

for p in range(2, 200):
    if not prime(p):
        continue
    expected = {triple(p + 1, p)}
    if p != 2:
        expected.add(triple(p + 1, 1))
    need(set(by_radius[p]) == expected, 'complete prime-inradius fibre')
need(triple(3, 1) == (8, 6, 10) and gcd(8, gcd(6, 10)) == 2,
     'p2 odd-odd content hostile')

# Twin enumeration uses Eratosthenes, independent of primality by trial division.
flags = bytearray(b'\x01') * 10003
flags[0] = flags[1] = 0
for p in range(2, isqrt(10002) + 1):
    if flags[p]:
        for k in range(p * p, 10003, p):
            flags[k] = 0
twins = []
for p in range(3, 10001, 2):
    if not flags[p] or not flags[p + 2]:
        continue
    q, m = p + 2, p + 1
    need(prime(p) and prime(q), 'twin sieve/trial agreement')
    T = triple(m, 1)
    A, B, C = T
    need((A, B, C) == (p * q, 2 * m, m * m + 1), 'twin sides')
    need(radii(T) == (p, q), 'twin radii')
    need(aliquot(2 * p) == literal_aliquot(2 * p) == p + 3,
         'semiprime incoming, independent divisor path')
    need(aliquot(q * q) == literal_aliquot(q * q) == p + 3,
         'prime-square incoming, independent divisor path')
    need(2 * p != q * q and (2 * p == p + 3) == (p == 3), 'fixed boundary')
    if p >= 5:
        need(aliquot(6 * p) == 6 * q, 'guarded perfect-core edge')
        need(m % 6 == 0 and B % 12 == 0, 'nonexceptional centre modulus')
    abc_check(1, A, m * m)
    R0 = rad(A * m * m)
    need(R0 == A * rad(m) and R0 > m * m, 'centre radical bound')
    abc_check(A * A, B * B, C * C)
    R1 = rad(A) * rad(B) * rad(C)
    need(gcd(A, C) == gcd(m, C) == 1, 'new norm prime support disjoint')
    need(R1 == R0 * rad(C), 'square-abc exact radical cost')
    need((m * m - 1, 2 * m) == (A, B), 'Gaussian square')
    # Formal monomial support is retained before either specialization.
    semiprime_support = [(0, 0), (1, 0), (0, 1)]
    square_support = [0, 1]
    Eleft = sum(2 ** e * p ** f for e, f in semiprime_support)
    Eright = sum(q ** e for e in square_support)
    need(Eleft == Eright == p + 3, 'weighted monoid collision')
    need(len(semiprime_support) == 3 and len(square_support) == 2,
         'augmentation distinguishes incoming witnesses')
    twins.append({'p': p, 'q': q, 'triple': list(T), 'target': p + 3,
                  'centre_radical': R0, 'square_radical': R1})

# Polynomial coefficients of (t^2-1)^2+(2t)^2=(t^2+1)^2.
def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out
left = poly_mul([-1, 0, 1], [-1, 0, 1])
left[2] += 4
need(left == poly_mul([1, 0, 1], [1, 0, 1]), 'symbolic polynomial abc identity')

# Proper-divisor multicomplex evaluated literally, not via sigma product.
for N in range(2, 1001):
    F = factors(N)
    top = tuple(e for _, e in F)
    D = [es for es in product(*(range(e + 1) for _, e in F)) if es != top]
    weighted = sum(prod(p ** e for (p, _), e in zip(F, es)) for es in D)
    need(weighted == literal_aliquot(N), 'divisor polynomial evaluation')
    need(len(D) == prod(e + 1 for _, e in F) - 1, 'augmentation counts divisors')

# Target-complete fibres, using the proved global source bound 99^2.
bound = 99 * 99
sieve = [0] * (bound + 1)
for d in range(1, bound // 2 + 1):
    for k in range(2 * d, bound + 1, d):
        sieve[k] += d
fibres = {M: [] for M in range(2, 101)}
for N in range(1, bound + 1):
    need(sieve[N] == aliquot(N), 'additive sieve versus factor-product sigma')
    if sieve[N] in fibres:
        M = sieve[N]
        need(N <= (M - 1) ** 2, 'targetwise inverse bound')
        fibres[M].append(N)
need(fibres[6] == [6, 25] and fibres[8] == [10, 49], 'first twin fibres')
need(fibres[74] == [70, 142, 5329], 'complete extra-preimage hostile')
nonexhaustive = []
for row in twins:
    p, q, M = row['p'], row['q'], row['target']
    if M > 100:
        continue
    need(2 * p in fibres[M] and q * q in fibres[M], 'marked witnesses in complete fibre')
    if len(fibres[M]) > 2:
        nonexhaustive.append(p)
need(nonexhaustive[0] == 71, 'first bounded nonexhaustion')

# Hostile controls for every dropped predicate or sidecar.
need(radii((3, 4, 5)) == (1, 3), 'root is not prime-radius triangle')
need(radii((7, 24, 25)) == (3, 21), 'second prime-inradius branch')
need(radii((63, 16, 65)) == (7, 9), 'prime inradius alone fails')
need(omega(15) == 2 and triple(8, 7) == (15, 112, 113), 'semiprime leg alone fails')
need(triple(12, 1) == (143, 24, 145) and triple(9, 8) == (17, 144, 145),
     'equal norms distinct primitive lifts')
need(radii((143, 24, 145)) == (11, 13) and
     radii((17, 144, 145)) == (8, 136), 'norm forgets twin radii')
need(not prime(145) and prime(11) and prime(13), 'twins do not force prime norm')
need(prime(101) and not prime(9) and triple(10, 1) == (99, 20, 101),
     'prime norm does not force twins')
need(aliquot(18) == 21 != 30, 'perfect-core domain at p3')
need(aliquot(9) == 4 != 6 == aliquot(25), 'equal augmentation unequal arithmetic')
need(aliquot(6) == aliquot(25) == 6, 'ordinary source abc distinction')
need((6 // 6, aliquot(6) // 6, (6 + aliquot(6)) // 6) == (1, 1, 2), 'perfect source normalization')
need((25, aliquot(25), 25 + aliquot(25)) == (25, 6, 31), 'square source normalization')
need((1, 3 * 5, 4 * 4) == (1, 15, 16), 'distinct centre abc at same twin pair')

print(json.dumps({
    'status': 'PASS', 'gates': gates,
    'universes': {'euclid_m_max': 200, 'direct_hypotenuse_max': 1000,
                  'twin_lower_prime_max': 10000, 'prime_radius_max': 199,
                  'divisor_polynomial_N_max': 1000, 'aliquot_target_max': 100,
                  'proved_common_source_bound': bound},
    'euclid_triples': len(euclid), 'direct_side_triples': len(direct),
    'prime_radius_pairs_in_euclid_universe': prime_radii,
    'twin_pairs': len(twins), 'first_twin_rows': twins[:8],
    'first_nonexhaustive_twin_p_through_97': nonexhaustive[0],
    'complete_fibre_74': fibres[74],
    'normal_optimized_contract': 'same source; explicit RuntimeError gates; deterministic JSON'
}, indent=2, sort_keys=True))
