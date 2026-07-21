#!/usr/bin/env python3
"""Exact checks for THM-2041.

The theorem is algebraic; this script stress-checks both exact-period kernels,
the cyclic group-algebra Frobenius identity, and a non-invariant negative
control. Tournament Analysis is not natural here: the observable is equality
under a group automorphism, and orienting it would discard the tied projector.
"""

from math import gcd


PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23)


def mobius(n: int) -> int:
    value = 1
    d = 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            value = -value
            if n % d == 0:
                return 0
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        value = -value
    return value


def divisors(n: int):
    return [d for d in range(1, n + 1) if n % d == 0]


def ramanujan(q: int, n: int) -> int:
    """c_q(n)=sum_{d|gcd(q,n)} d*mu(q/d), exact over Z."""
    g = gcd(q, abs(n))
    return sum(d * mobius(q // d) for d in divisors(g))


def mul(a, b, p):
    q = len(a)
    out = [0] * q
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[(i + j) % q] = (out[(i + j) % q] + x * y) % p
    return out


def power(a, n, p):
    out = [0] * len(a)
    out[0] = 1
    base = list(a)
    while n:
        if n & 1:
            out = mul(out, base, p)
        base = mul(base, base, p)
        n //= 2
    return out


def sigma(a, p):
    q = len(a)
    out = [0] * q
    for i, x in enumerate(a):
        out[(p * i) % q] = x
    return out


def primitive_kernel(q, p):
    return [int(gcd(a, q) == 1) % p for a in range(q)]


def ramanujan_kernel(q, p):
    return [ramanujan(q, a) % p for a in range(q)]


ramanujan_checks = 0
for q in range(2, 61):
    for p in PRIMES:
        if gcd(p, q) != 1:
            continue
        for n in range(-80, 81):
            assert ramanujan(q, p * n) == ramanujan(q, n)
            ramanujan_checks += 1


frobenius_checks = 0
for q in range(2, 25):
    for p in PRIMES:
        if gcd(p, q) != 1:
            continue
        lam = [((a + 1) * (a + 3) + 2 * q + p) % p for a in range(q)]
        for theta in (primitive_kernel(q, p), ramanujan_kernel(q, p)):
            assert sigma(theta, p) == theta
            assert power(theta, p, p) == theta
            for m in range(1, 6):
                lhs = mul(theta, power(lam, p * m, p), p)[0]
                base = mul(theta, power(lam, m, p), p)[0]
                rhs = pow(base, p, p)
                assert lhs == rhs
                frobenius_checks += 1


negative_control = None
for q in range(5, 16):
    for p in PRIMES:
        if gcd(p, q) != 1:
            continue
        theta = [0] * q
        theta[1] = 1
        theta[2] = 2 % p
        if sigma(theta, p) == theta:
            continue
        lam = [((3 * a * a + a + 1) % p) for a in range(q)]
        for m in range(1, 6):
            lhs = mul(theta, power(lam, p * m, p), p)[0]
            rhs = pow(mul(theta, power(lam, m, p), p)[0], p, p)
            if lhs != rhs:
                negative_control = (q, p, m, lhs, rhs)
                break
        if negative_control:
            break
    if negative_control:
        break
assert negative_control is not None


units14 = [a for a in range(14) if gcd(a, 14) == 1]
orbits14 = {}
for p in PRIMES:
    if gcd(p, 14) == 1:
        orbits14[p] = [((p * a) % 14) for a in units14]
        assert sorted(orbits14[p]) == units14


print("THM-2041 FROBENIUS EXACT-PERIOD AUDIT")
print(f"Ramanujan dilation identities checked: {ramanujan_checks}")
print(f"Cyclic twisted-moment identities checked: {frobenius_checks}")
print(f"Non-invariant twist failure (q,p,m,lhs,rhs): {negative_control}")
print(f"q=14 primitive phase set: {units14}")
for p, image in orbits14.items():
    print(f"  multiply by p={p}: {image}")
print("PASS: primitive/Ramanujan projectors survive as whole Frobenius layers.")
print("SCOPE: preservation only; no base nonvanishing or pointwise LRC lift is proved.")
