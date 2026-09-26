"""Exact controls: Bernoulli endpoints, dyadic jumps, cylinders, power sums.

All arithmetic in the consequence checks is integral or Fraction arithmetic.
No convergence conclusion is inferred from finite universes.
"""
from fractions import Fraction as Q
from math import comb, factorial, prod
from itertools import product
from pathlib import Path
import hashlib


def require(ok, context):
    if not ok:
        raise RuntimeError(context)


def psi(x):
    x = Q(x)
    return x - x.numerator // x.denominator - Q(1, 2)


def jump(z, m):
    return Q(1, m) + psi(Q(z - 1, m)) - psi(Q(z, m))


def v2(z):
    z = abs(z)
    if z == 0:
        raise ValueError('zero has infinite valuation')
    return (z & -z).bit_length() - 1


def bernoulli(N):
    b = [Q(1)]
    for n in range(1, N + 1):
        b.append(-sum(comb(n + 1, j) * b[j] for j in range(n)) / (n + 1))
    return b


def bp(n, x, b):
    return sum(comb(n, j) * b[j] * Q(x)**(n - j) for j in range(n + 1))


def step(n, multiplier=3, sign=1):
    z = multiplier * n + sign
    k = v2(z)
    return z >> k, k


def main():
    b = bernoulli(40)
    require(b[1] == -Q(1, 2), 'B1 convention')
    endpoint_count = 0
    for n in range(1, 41):
        require(bp(n, 1, b) - bp(n, 0, b) == int(n == 1), ('endpoints', n))
        require(bp(n, 1, b) == (-1)**n * b[n], ('reflection', n))
        if n % 2:
            require(bp(n, Q(1, 2), b) == 0, ('midpoint', n))
        endpoint_count += 1
    require(bp(3, Q(1, 4), b) == Q(3, 64), 'odd polynomial hostile')
    print('Bernoulli endpoint/reflection degrees 1..40:', endpoint_count, 'PASS')
    print('Shift hostile: B3=0 but B3(1/4)=3/64')

    count = 0
    for z in range(-1024, 1025):
        for m in range(2, 65):
            require(jump(z, m) == int(z % m == 0), ('jump', z, m))
            count += 1
    print('Integer divisibility jump controls:', count, 'PASS')
    for z in range(-2048, 2049):
        total = sum(jump(z, 1 << s) for s in range(1, 14))
        require(total == (13 if z == 0 else v2(z)), ('valuation', z))
    print('Dyadic valuation controls: 4097 inputs, 13 scales, zero retained PASS')

    count = 0
    for B in range(1, 33):
        for a in range(1, B + 1):
            for N in range(97):
                observed = sum(j % B == a % B for j in range(1, N + 1))
                formula = Q(N, B) + psi(-Q(a, B)) - psi(Q(N - a, B))
                require(observed == formula, ('cylinder', B, a, N))
                count += 1
    print('Exact cylinder counts:', count, 'PASS')
    for m in range(2, 10):
        for j in range(-34, 35):
            x = Q(j, 17)
            require(sum(psi((x + r) / m) for r in range(m)) == psi(x),
                    ('multiplication', m, x))
    print('Complete-child Bernoulli cancellation: 552 cases PASS')

    for q, sign, n, word in ((3, -1, 5, (1, 2)),
                             (5, 1, 13, (1, 1, 5)), (3, 1, 1, (2,))):
        state = n
        for k in word:
            state, actual = step(state, q, sign)
            require(k == actual, ('hostile itinerary', q, sign, n))
        require(state == n, ('hostile cycle', n))
        a = (n + 1) // 2
        for t in range(1, 21):
            B = 2**(sum(word) * t)
            correction = psi(-Q(a, B)) - psi(Q(a - a, B))
            require(correction == 1 - Q(a, B), ('boundary atom', q, n, t))
        print('Persistent boundary atom:', (q, sign, n, word), '20 depths PASS')

    count = 0
    for k in range(1, 41):
        for m in range(2, 41):
            literal = Q(sum(j**k for j in range(1, m)), m**k)
            expression = Q(m, k + 1) - Q(1, 2)
            for r in range(1, k // 2 + 1):
                expression += (b[2*r] / factorial(2*r)
                               * prod(range(k - 2*r + 2, k + 1))
                               / m**(2*r - 1))
            require(literal == expression, ('Faulhaber', k, m))
            count += 1
    print('Normalized power sums:', count, 'PASS')
    require(2*b[1] == -1, 'even prime 2')
    require(all(b[n-1] == 0 for n in range(4, 42, 2)), 'even composite gate')
    require(all((30//p - 1) % p == 0 for p in (2, 3, 5)), 'Giuga 30')
    require(all(560 % (p-1) == 0 for p in (3, 11, 17)), 'Carmichael 561')
    require((561//11) % 11 == 7, '561 not Giuga')
    print('Agoh-Giuga parity and separated-condition hostiles PASS')
    print('script_sha256', hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
