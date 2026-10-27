#!/usr/bin/env python3
"""Recover prime-tournament aggregate cycle index from odd-partition Burnside data."""
from fractions import Fraction as F
from collections import Counter
from math import gcd, factorial
import json, sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 12

def add(a, b, scale=1):
    c = a.copy()
    for k, v in b.items():
        c[k] = c.get(k, 0) + scale * v
        if not c[k]:
            del c[k]
    return c

def mul(a, b):
    c = {}
    for p, v in a.items():
        for q, w in b.items():
            if sum(p) + sum(q) <= N:
                k = tuple(sorted(p + q))
                c[k] = c.get(k, 0) + v * w
    return {k: v for k, v in c.items() if v}

def scale_vars(a, d):
    return {tuple((d * i for i in p)): v for p, v in a.items() if d * sum(p) <= N}

def compose(p, a):
    cache = {d: scale_vars(a, d) for d in range(1, N + 1)}
    c = {}
    for part, v in p.items():
        term = {(): v}
        for d in part:
            term = mul(term, cache[d])
        c = add(c, term)
    return c

def odd_parts(n, minpart=1):
    if n == 0:
        yield ()
    for k in range(minpart, n + 1, 2):
        for p in odd_parts(n - k, k):
            yield ((k,) + p)

def burnside(part):
    e = sum(((k - 1) // 2 for k in part)) + sum((gcd(k, l) for i, k in enumerate(part) for l in part[i + 1:]))
    z = 1
    for k, m in Counter(part).items():
        z *= k ** m * factorial(m)
    return F(2 ** e, z)
ZT = {p: burnside(p) for n in range(N + 1) for p in odd_parts(n)}
A = add(ZT, {(): 1}, -1)
inv = {(): F(1)}
term = inv.copy()
for k in range(1, N + 1):
    term = mul(term, A)
    inv = add(inv, term, (-1) ** k)
ZS = add({(): F(1)}, inv, -1)
R = add(ZS, {(1,): 1}, -1)
P = {}
for n in range(3, N + 1):
    residual = add(R, compose(P, A), -1)
    P = add(P, {p: v for p, v in residual.items() if sum(p) == n})
if not compose(P, A) == R:
    raise RuntimeError('Verification failed: compose(P, A) == R')
if not all((v >= 0 for v in P.values())):
    raise RuntimeError('Verification failed: all((v >= 0 for v in P.values()))')
if not not any((1 in p and sum((d > 1 for d in p)) == 1 for p in P)):
    raise RuntimeError('Verification failed: not any((1 in p and sum((d > 1 for d in p)) == 1 for p in P))')
for part in P:
    for i, d in enumerate(part):
        if d > 1 and len(part) > 1 and all(gcd(d, e) == 1 for j, e in enumerate(part) if i != j):
            raise RuntimeError('Forbidden coprime-isolated nontrivial cycle in prime action: ' + str(part))
prime_counts = [sum((v for p, v in P.items() if sum(p) == n)) for n in range(N + 1)]
if not all((v.denominator == 1 for v in prime_counts)):
    raise RuntimeError('Verification failed: all((v.denominator == 1 for v in prime_counts))')
if not list(map(int, prime_counts[3:11])) == [1, 0, 3, 15, 197, 4008, 136655, 7906161][:max(0, min(N - 2, 8))]:
    raise RuntimeError('Verification failed: list(map(int, prime_counts[3:11])) == [1, 0, 3, 15, 197, 4008, 136655, 7906161][:max(0, min(N - 2, 8))]')
if N >= 6:
    from pathlib import Path
    import subprocess
    receipt = json.loads(subprocess.check_output([sys.executable, str(Path(__file__).with_name('verify_prime_cycle_index.py'))], text=True))
    explicit = {}
    for obj in receipt['prime_cycle_indices']:
        for ct, count in obj['cycles']:
            key = tuple(ct)
            explicit[key] = explicit.get(key, 0) + F(count, obj['aut_order'])
    if not explicit == {p: v for p, v in P.items() if sum(p) <= 6}:
        raise RuntimeError('Verification failed: explicit == {p: v for p, v in P.items() if sum(p) <= 6}')
# Closed scalar boundary: full-cycle coefficients obey a Dirichlet convolution.
f = {n: 2 ** ((n - 1) // 2) for n in range(1, N + 1, 2)}
g = {1: 0}
for n in range(3, N + 1, 2):
    g[n] = f[n] - sum(g[d] * f[n // d] for d in range(3, n, 2) if n % d == 0)
    if g[n] != n * P.get((n,), 0):
        raise RuntimeError(('Full-cycle coefficient mismatch', n, g[n], P.get((n,), 0)))
    if n <= 15 and N >= 6 and g[n] != receipt['specified_full_cycle_prime_counts'][str(n)]:
        raise RuntimeError(('Direct circulant enumeration mismatch', n))
report = {'N': N, 'specified_full_cycle_prime_counts': g, 'all_tournaments': [int(sum((v for p, v in ZT.items() if sum(p) == n))) for n in range(N + 1)], 'prime_tournaments': list(map(int, prime_counts)), 'prime_cycle_index_monomials': len(P), 'aggregate_prime_cycle_index': {','.join(map(str, p)): str(v) for p, v in sorted(P.items(), key=lambda kv: (sum(kv[0]), kv[0]))}, 'checks': 'PASS: plethystic identity in every coefficient; nonnegative rational coefficients; integral prime counts; coprime-isolated-cycle exclusion; full-cycle divisor recurrence and direct circulant census through15; independent labelled-universe cycle-index comparison through6; prior census through10'}
print(json.dumps(report, indent=2))
