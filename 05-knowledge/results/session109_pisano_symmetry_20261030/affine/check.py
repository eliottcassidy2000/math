#!/usr/bin/env python3
"""Stdlib exact independent finite checks; checks remain active under -O."""
from math import gcd, lcm
from collections import Counter
import json

def require(condition, message):
    if not condition:
        raise AssertionError(message)
I = (1, 0, 0, 1)
Q = (1, 1, 1, 0)

def mul(A, B, n=None):
    a, b, c, d = A
    e, f, g, h = B
    C = (a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h)
    return tuple((x % n for x in C)) if n else C

def power(A, r, n=None):
    B = I
    while r:
        if r & 1:
            B = mul(B, A, n)
        A = mul(A, A, n)
        r //= 2
    return B

def mv(A, v, n):
    return ((A[0] * v[0] + A[1] * v[1]) % n, (A[2] * v[0] + A[3] * v[1]) % n)

def add(v, w, n):
    return tuple(((a + b) % n for a, b in zip(v, w)))

def smith(r):
    A = power(Q, r)
    entries = (A[0] - 1, A[1], A[2], A[3] - 1)
    d = gcd(*entries)
    D = abs(entries[0] * entries[3] - entries[1] * entries[2])
    return (d, D // d)

def fib(n):
    a, b = (0, 1)
    for _ in range(n):
        a, b = (b, a + b)
    return a

def luc(n):
    return 2 if not n else fib(n - 1) + fib(n + 1)

def order(A, n):
    B = I
    t = 0
    while True:
        B = mul(B, A, n)
        t += 1
        if B == I:
            return t
        require(t <= 6 * n * n, 'exact check failed at original line 32')

def affine_order_brute(A, b, n):
    B = I
    c = (0, 0)
    t = 0
    while True:
        c = add(mv(A, c, n), b, n)
        B = mul(A, B, n)
        t += 1
        if B == I and c == (0, 0):
            return t
        require(t <= 6 * n ** 3, 'exact check failed at original line 39')

def return_carry(A, b, n, t):
    c = (0, 0)
    for _ in range(t):
        c = add(mv(A, c, n), b, n)
    return c

def predicted(A, b, n):
    t = order(A, n)
    c = return_carry(A, b, n, t)
    return t * (n // gcd(n, *c))
checks = 0
for r in range(1, 201):
    d, e = smith(r)
    expected = (2 if r % 3 == 0 else 1) if r % 2 else luc(r // 2) if r % 4 == 2 else fib(r // 2)
    require(d == expected and d * e == luc(r) - 1 - (-1) ** r and (e % d == 0), 'exact check failed at original line 53')
for n in range(2, 16):
    for r in range(1, 13):
        A = power(Q, r, n)
        d, e = smith(r)
        image = {((x - mv(A, (x, y), n)[0]) % n, (y - mv(A, (x, y), n)[1]) % n) for x in range(n) for y in range(n)}
        require(n * n // len(image) == gcd(n, d) * gcd(n, e), 'exact check failed at original line 58')
        fix = sum((mv(A, (x, y), n) == (x, y) for x in range(n) for y in range(n)))
        require(fix == gcd(n, d) * gcd(n, e), 'exact check failed at original line 60')
        for x in range(n):
            for y in range(n):
                b = (x, y)
                actual = affine_order_brute(A, b, n)
                require(actual == predicted(A, b, n), 'exact check failed at original line 64')
                if b in image:
                    require(actual == order(A, n), 'exact check failed at original line 65')
                checks += 1
for n, m in [(3, 4), (4, 5), (7, 11), (8, 9)]:
    for r in range(1, 13):
        for b in [(0, 0), (1, 0), (1, 3), (7, 9)]:
            require(predicted(power(Q, r, n * m), b, n * m) == lcm(predicted(power(Q, r, n), b, n), predicted(power(Q, r, m), b, m)), 'exact check failed at original line 70')
examples = []
for n, r in [(11, 5), (121, 5), (4, 6), (19, 9), (5, 4)]:
    A = power(Q, r, n)
    t = order(A, n)
    hist = Counter((predicted(A, (x, y), n) for x in range(n) for y in range(n)))
    examples.append({'N': n, 'stride': r, 'linear_order': t, 'affine_order_distribution': dict(sorted(hist.items())), 'carry_b10': return_carry(A, (1, 0), n, t)})
result = {'status': 'PASS', 'smith_strides_checked': 200, 'all_carry_bruteforce_cases': checks, 'smith_table': {r: smith(r) for r in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 18]}, 'examples': examples}
print(json.dumps(result, indent=2))
