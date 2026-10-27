#!/usr/bin/env python3
"""Exhaustive labelled-universe -> isotype census, then independent composition check."""
from itertools import permutations, combinations
from collections import Counter
from fractions import Fraction as F
import json
N = 6

def cycle_type(p):
    unseen = set(range(len(p)))
    out = []
    while unseen:
        v = next(iter(unseen))
        k = 0
        while v in unseen:
            unseen.remove(v)
            k += 1
            v = p[v]
        out.append(k)
    return tuple(sorted(out))

def census(n):
    pairs = list(combinations(range(n), 2))
    index = {p: i for i, p in enumerate(pairs)}
    transforms = []
    for p in permutations(range(n)):
        table = []
        for i, j in pairs:
            a, b = (p[i], p[j])
            table.append((index[min(a, b), max(a, b)], a > b))
        transforms.append((p, table))
    seen = set()
    reps = []
    prime = []
    strong = 0
    for code in range(1 << len(pairs)):
        if code in seen:
            continue
        autos = []
        for p, table in transforms:
            image = sum(((code >> k & 1 ^ flip) << dest for k, (dest, flip) in enumerate(table)))
            seen.add(image)
            if image == code:
                autos.append(p)
        reps.append(code)

        def arc(a, b):
            if a > b:
                return not arc(b, a)
            return bool(code >> index[a, b] & 1)
        reach = [set([v] + [w for w in range(n) if w != v and arc(v, w)]) for v in range(n)]
        for k in range(n):
            for v in range(n):
                if k in reach[v]:
                    reach[v] |= reach[k]
        strong += all((len(r) == n for r in reach))
        isprime = True
        for size in range(2, n):
            for M in combinations(range(n), size):
                if all((len({arc(v, m) for m in M}) == 1 for v in range(n) if v not in M)):
                    isprime = False
                    break
            if not isprime:
                break
        if isprime and n >= 3:
            ct = Counter((cycle_type(p) for p in autos))
            prime.append({'n': n, 'code': code, 'aut_order': len(autos), 'cycles': [(list(c), v) for c, v in sorted(ct.items())], 'vertex_orbits': sum((sum((k == 1 for k in c)) * v for c, v in ct.items())) // len(autos)})
    if not len(seen) == 1 << len(pairs):
        raise RuntimeError('Verification failed: len(seen) == 1 << len(pairs)')
    return (len(reps), strong, prime)

def mul(a, b):
    c = [F(0)] * (N + 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:N + 1 - i]):
            c[i + j] += x * y
    return c

def subst(t, d):
    return [F(0) if i % d else t[i // d] for i in range(N + 1)]
counts = [1]
strongs = [0]
primes = []
for n in range(1, N + 1):
    a, s, p = census(n)
    counts.append(a)
    strongs.append(s)
    primes += p
T = [F(1)] + [F(0)] * N
for iteration in range(N + 1):
    A = T.copy()
    A[0] = 0
    S = [F(0)] * (N + 1)
    S[1] = 1
    for p in primes:
        for ct, multiplicity in p['cycles']:
            term = [F(1)] + [F(0)] * N
            for d in ct:
                term = mul(term, subst(A, d))
            for i in range(N + 1):
                S[i] += F(multiplicity, p['aut_order']) * term[i]
    newT = [F(1)] + [F(0)] * N
    for n in range(1, N + 1):
        newT[n] = sum((S[k] * newT[n - k] for k in range(1, n + 1)))
    T = newT
if not T == counts:
    raise RuntimeError((T, counts))
if not S == strongs:
    raise RuntimeError((S, strongs))
if not counts == [1, 1, 1, 2, 4, 12, 56]:
    raise RuntimeError('Verification failed: counts == [1, 1, 1, 2, 4, 12, 56]')
if not [sum((p['n'] == n for p in primes)) for n in range(3, 7)] == [1, 0, 3, 15]:
    raise RuntimeError("Verification failed: [sum((p['n'] == n for p in primes)) for n in range(3, 7)] == [1, 0, 3, 15]")
if not sorted((p['vertex_orbits'] for p in primes if p['n'] == 5)) == [1, 5, 5]:
    raise RuntimeError("Verification failed: sorted((p['vertex_orbits'] for p in primes if p['n'] == 5)) == [1, 5, 5]")
# Independent fixed-n-cycle check: enumerate all difference-sign choices and ALL modules.
def circulant_prime_count(n):
    full = (1 << n) - 1
    total = 0
    modules = [m for m in range(1, full) if m.bit_count() >= 2]
    for choice in range(1 << ((n - 1) // 2)):
        differences = [d if (choice >> (d - 1)) & 1 else n - d for d in range(1, (n + 1) // 2)]
        rows = [sum(1 << ((v + d) % n) for d in differences) for v in range(n)]
        prime = True
        for module in modules:
            outside = full ^ module
            uniform = True
            while outside:
                low = outside & -outside
                v = low.bit_length() - 1
                if rows[v] & module not in (0, module):
                    uniform = False
                    break
                outside ^= low
            if uniform:
                prime = False
                break
        total += prime
    return total
circulant_counts = {n: circulant_prime_count(n) for n in range(3, 16, 2)}
if circulant_counts != {3: 2, 5: 4, 7: 8, 9: 12, 11: 32, 13: 64, 15: 112}:
    raise RuntimeError(('Unexpected full-cycle prime count', circulant_counts))
report = {'max_n': N, 'specified_full_cycle_prime_counts': circulant_counts, 'all_isotypes': counts, 'strong_isotypes': strongs, 'prime_counts': [sum((p['n'] == n for p in primes)) for n in range(N + 1)], 'prime_cycle_indices': primes, 'composition_counts': list(map(int, T)), 'checks': 'PASS: complete labelled universe, exact rational composition, strong counts, five-prime vertex orbit hostile'}
print(json.dumps(report, indent=2))
