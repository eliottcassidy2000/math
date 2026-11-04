#!/usr/bin/env python3
"""Exact, stdlib-only checks for the Euler-aware odd-order support sieve."""
from itertools import combinations
from math import gcd
from fractions import Fraction
import json

def isprime(n):
    return n >= 2 and all((n % d for d in range(2, 1 + __import__('math').isqrt(n))))

def order(p, r):
    if not (isprime(r) and p % r):
        raise AssertionError('Verification failed')
    x, d = (p % r, 1)
    while x != 1:
        x = x * p % r
        d += 1
    return d

def v2(n):
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k

def graph(S, q):
    return {(p, r) for p in S for r in S if p != r and v2(order(p, r)) <= (1 if p == q else 0)}

def closure(S, E):
    R = {p: {p} | {r for a, r in E if a == p} for p in S}
    for k in S:
        for p in S:
            if k in R[p]:
                R[p] |= R[k]
    return R

def strongly_connected(S, E):
    return all((len(v) == len(S) for v in closure(S, E).values()))

def lcm(a, b):
    return a // gcd(a, b) * b

def realizing_exponents(S, q, E):
    out = {}
    for p in S:
        n = 2 if p == q else 1
        for a, r in E:
            if a == p:
                d = order(p, r)
                n = lcm(n, r if d == 1 else d)
        if p != q and n == 1:
            n = 3
        out[p] = n - 1
    return out

def sigma_mod(p, a, r):
    n = a + 1
    power, total = (1, 0)
    bp, bs = (p % r, 1)
    while n:
        if n & 1:
            total = (total + power * bs) % r
            power = power * bp % r
        bs = bs * (1 + bp) % r
        bp = bp * bp % r
        n //= 2
    return total

def run():
    checks = 0
    P = [p for p in range(3, 30, 2) if isprime(p)]
    for size in range(2, 6):
        for S in combinations(P, size):
            for q in S:
                if q % 4 != 1:
                    continue
                E = graph(S, q)
                A = realizing_exponents(S, q, E)
                if not A[q] % 4 == 1:
                    raise AssertionError('Verification failed')
                if not all((a >= 2 and a % 2 == 0 for p, a in A.items() if p != q)):
                    raise AssertionError('Verification failed')
                actual = {(p, r) for p in S for r in S if p != r and sigma_mod(p, A[p], r) == 0}
                if not E == actual:
                    raise AssertionError((S, q, A, E, actual))
                checks += 1
    direct = 0
    for p in P:
        for r in P:
            if p == r:
                continue
            for a in range(1, 42):
                if not sigma_mod(p, a, r) == sum((pow(p, i, r) for i in range(a + 1))) % r:
                    raise AssertionError('Verification failed')
                if sigma_mod(p, a, r) == 0:
                    if a % 2 == 0:
                        if not v2(order(p, r)) == 0:
                            raise AssertionError('Verification failed')
                    if a % 4 == 1:
                        if not v2(order(p, r)) <= 1:
                            raise AssertionError('Verification failed')
                direct += 1
    S = (3, 5, 23, 31, 43)
    q = 5
    E = graph(S, q)
    C = {3, 23, 43}
    R = closure(S, E)
    if not all((R[p] == C for p in C)):
        raise AssertionError('Verification failed')
    if not not strongly_connected(S, E):
        raise AssertionError('Verification failed')
    if not all((any((a == p for a, r in E)) and any((r == p for a, r in E)) for p in S)):
        raise AssertionError('Verification failed')
    H = Fraction(1)
    for p in S:
        H *= Fraction(p, p - 1)
    if not H > 2:
        raise AssertionError('Verification failed')
    return {'realized_support_graphs': checks, 'direct_sigma_checks': direct, 'excluded_support': S, 'Euler_candidate': q, 'abundancy_ceiling': str(H), 'closed_sink': sorted(C), 'Euler_aware_edges': sorted(E), 'realizing_exponents': realizing_exponents(S, q, E), 'status': 'PASS; arithmetic relaxation only, not an OPN nonexistence proof'}
if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
