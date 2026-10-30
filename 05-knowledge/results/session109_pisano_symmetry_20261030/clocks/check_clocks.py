#!/usr/bin/env python3
"""Independent exact finite checks; stdlib only; python3 check_clocks.py."""
import json
from math import gcd, lcm
from pathlib import Path
Q = (1, 1, 1, 0)
A = (2, 1, 1, 1)

def mul(a, b, n):
    return tuple((sum((a[2 * i + k] * b[2 * k + j] for k in range(2))) % n for i in range(2) for j in range(2)))

def powers(b, n):
    ident = (1 % n, 0, 0, 1 % n)
    ans = [ident]
    x = mul(ident, b, n)
    while x != ident:
        ans.append(x)
        x = mul(x, b, n)
        if not len(ans) < 100 * n * n:
            raise AssertionError('verification failed')
    return ans

def period(word):
    m = len(word)
    return next((d for d in range(1, m + 1) if m % d == 0 and all((word[i] == word[(i + d) % m] for i in range(m)))))

def clocks(b, n):
    p = powers(b, n)
    m = len(p)
    scalar = next((i for i, x in enumerate(p[1:], 1) if x[1] == x[2] == 0 and x[0] == x[3]), m)
    tr = period([(x[0] + x[3]) % n for x in p])
    vector = period([(x[0], x[2]) for x in p])
    return (m, scalar, tr, vector)
cache = {}
tests = 0
examples = []
for n in range(2, 401):
    for name, b in [('Q', Q), ('A', A)]:
        m, s, t, v = clocks(b, n)
        cache[name, n] = (m, s, t, v)
        if not (v == m and m % s == 0 and (m // s in ((1, 2, 4) if name == 'Q' else (1, 2)))):
            raise AssertionError('verification failed')
        if not m // t in (1, 5):
            raise AssertionError('verification failed')
        u = n
        e = 0
        while u % 5 == 0:
            u //= 5
            e += 1
        non5 = clocks(b, u)[0] if u > 1 else 1
        expect = non5 if e == 0 else lcm(non5, (4 if name == 'Q' else 2) * 5 ** (e - 1))
        if not t == expect:
            raise AssertionError((name, n, m, t, expect))
        for shift, x in enumerate(powers(b, n)):
            pair_ok = (x[0] + x[3] - 2) % n == 0 and (sum(mul(b, x, n)[::3]) - sum(b[::3])) % n == 0
            if not pair_ok == (shift % t == 0):
                raise AssertionError('verification failed')
        tests += 1
for name, b in [('Q', Q), ('A', A)]:
    for n in (5, 25, 125, 625, 3125):
        m, s, t, v = clocks(b, n)
        if not (m == (4 if name == 'Q' else 2) * n and t == m // 5):
            raise AssertionError('verification failed')
        examples.append(dict(matrix=name, modulus=n, matrix_clock=m, scalar_clock=s, trace_clock=t, marked_clock=v))
    for a in range(2, 31):
        for c in range(2, 31):
            if gcd(a, c) != 1:
                continue
            ca = clocks(b, a)
            cc = clocks(b, c)
            cac = clocks(b, a * c)
            if not cac == tuple((lcm(x, y) for x, y in zip(ca, cc))):
                raise AssertionError('verification failed')
            tests += 1
pa = powers(A, 5)
scalar = (4, 0, 0, 4)
if not (len(pa) == 10 and len(powers(scalar, 5)) == 2):
    raise AssertionError('verification failed')
x = (1, 0, 0, 1)
y = x
for i in range(20):
    if not (x[0] + x[3]) % 5 == (y[0] + y[3]) % 5:
        raise AssertionError('verification failed')
    x = mul(x, A, 5)
    y = mul(y, scalar, 5)
p = powers(A, 15)
if not (len(p) == 20 and p[10] == (11, 0, 0, 11)):
    raise AssertionError('verification failed')
if not not any((x == (14, 0, 0, 14) for x in p)):
    raise AssertionError('verification failed')
sidecar_cases = 0
for n in [5, 10, 15, 25, 50, 55, 125]:
    u = n
    e = 0
    while u % 5 == 0:
        u //= 5
        e += 1
    for name, b in [('Q', Q), ('A', A)]:
        seen = set()
        for aa in range(n):
            for bb in range(n):
                x = tuple(((aa * (i in (0, 3)) + bb * b[i]) % n for i in range(4)))
                key = ((x[0] + x[3]) % n, sum(mul(b, x, n)[::3]) % n, bb % 5 ** e // 5 ** (e - 1))
                if not key not in seen:
                    raise AssertionError((n, name, aa, bb, key))
                seen.add(key)
                sidecar_cases += 1

def amul(x, y, n):
    a, b = x
    c, d = y
    return ((a * c + b * d) % n, (a * d + b * c + b * d) % n)

def atr(x, n):
    return (2 * x[0] + x[1]) % n
ramified_cases = 0
for k in range(1, 4):
    n = 5 ** k
    g = (2, 3)
    ideals = []
    S = {(a, b) for a in range(n) for b in range(n)}
    for j in range(2 * k + 1):
        ideals.append(S)
        S = {amul(x, (-1, 2), n) for x in S}
    for a in range(n):
        for b in range(n):
            x = (a, b)
            if x == (0, 0):
                continue
            val = max((j for j, I in enumerate(ideals) if x in I))
            r = 2 * k - val
            orbit = [x]
            y = amul(g, x, n)
            while y != x:
                orbit.append(y)
                y = amul(g, y, n)
            if not len(orbit) == 5 ** (r // 2):
                raise AssertionError('verification failed')
            if not period([atr(z, n) for z in orbit]) == 5 ** ((r - 1) // 2):
                raise AssertionError('verification failed')
            ramified_cases += 1
summary = {'status': 'PASS', 'sidecar_algebra_cases': sidecar_cases, 'ramified_all_element_cases': ramified_cases, 'moduli_full_sweep': '2..400 for Q and A', 'count_sweep_plus_CRT_cases': tests, 'prime_power_examples': examples, 'hostile_equal_trace': 'A versus -I mod 5: all traces agree, orders 10 versus 2', 'crt_sign_mismatch': 'A^10=11I mod15; scalar clock10, ±I clock20', 'masked_example': dict(zip(['matrix', 'scalar', 'trace', 'marked'], clocks(A, 55)))}
print(json.dumps(summary, indent=2))
Path(__file__).with_name('results.json').write_text(json.dumps(summary, indent=2) + '\n')
