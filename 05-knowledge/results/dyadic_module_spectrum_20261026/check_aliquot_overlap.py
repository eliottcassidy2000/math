#!/usr/bin/env python3
"""Exact tests of the factor-pair interpretation; no primality oracle imports."""
from collections import Counter
from math import isqrt, prod
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def factor(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            a = 0
            while n % p == 0:
                a += 1
                n //= p
            out.append((p, a))
        p = 3 if p == 2 else p + 2
    if n > 1:
        out.append((n, 1))
    return out


def sigma(n):
    return prod((p ** (a + 1) - 1) // (p - 1) for p, a in factor(n))


def pair_data(n):
    pairs = [(a, n // a) for a in range(1, isqrt(n) + 1) if n % a == 0]
    divisors = sorted({d for pair in pairs for d in pair})
    s = sum(divisors[:-1])
    return s, divisors, [list(pair) for pair in pairs if sum(pair) + 1 == s]


def run():
    limit = 1_000_000
    # First engine: the actual divisor incidence table, including n itself.
    sig = [0] * (limit + 1)
    tau = [0] * (limit + 1)
    for d in range(1, limit + 1):
        for n in range(d, limit + 1, d):
            sig[n] += d
            tau[n] += 1
    counts = Counter()
    exceptions = []
    for n in range(1, limit + 1):
        I = sig[n] - n - 1
        D = I * I - 4 * n
        cubic = I > 0 and D >= 0 and isqrt(D) ** 2 == D
        expected = tau[n] == 4 or sig[n] == 2 * n + 2
        require(cubic == expected, ('classification', n))
        if cubic:
            r = isqrt(D)
            require((I - r) % 2 == 0, ('root parity', n))
            a, b = (I - r) // 2, (I + r) // 2
            require(a > 0 and a * b == n and a + b + 1 == sig[n] - n,
                    ('decoded factor ports', n))
            counts['interior' if a > 1 else 'outer'] += 1
        if sig[n] == 2 * n + 2:
            f = factor(n)
            a, m = 0, n
            while m % 2 == 0:
                m //= 2
                a += 1
            c, x = sigma(m), sigma(m) - m
            require((2 ** (a + 1) - 1) * x == m + 2, ('tower', n))
            require(sig[n] == sigma(n), ('factored check', n))
            exceptions.append(dict(n=n, factors=f, proper_divisor_count=tau[n]-1,
                                   odd_part=m, v2=a, carry=c, odd_aliquot=x))
    # Independent engine: pair the actual divisors, including square fixed points.
    controls = sorted(set(range(1, 20_001)) | {x['n'] for x in exceptions})
    for n in controls:
        s, divisors, witnesses = pair_data(n)
        require(s == sig[n] - n and len(divisors) == tau[n], ('paired sieve', n))
        require(sigma(n) == sig[n], ('multiplicative sigma', n))
        require(bool(witnesses) == (tau[n] == 4 or sig[n] == 2*n+2), ('direct ports', n))
        require(all(a != b for a, b in witnesses), ('no false diagonal', n))
    require([x['n'] for x in exceptions] == [20, 104, 464, 650, 1952, 130304, 522752],
            'complete finite abundance-two list')
    require([x['n'] for x in exceptions if x['proper_divisor_count'] == 5] == [20],
            'quintic overlap')
    require([x['n'] for x in exceptions if x['proper_divisor_count'] == 11] == [650,1952],
            'eleven-term overlap')
    require([x['n'] for x in exceptions if x['proper_divisor_count'] <= 11] == [20,104,464,650,1952],
            'at-most-eleven overlap')
    # Separate exact algebraic checks for the three-prime shape (2,1,1).
    primes = [n for n in range(2,100) if tau[n] == 2]
    triples = 0
    for p in primes:
        for i,q in enumerate(primes):
            if q == p:
                continue
            for r in primes[i+1:]:
                if r == p:
                    continue
                n = p*p*q*r
                excess = (p*p+p+1)*(q+1)*(r+1)-2*n
                require((excess == 2) == (n == 650), ('p2qr classification control',p,q,r))
                triples += 1
    sample = {}
    for n in (4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 24, 25, 28, 35, 49, 104, 343, 650, 1001):
        s, ds, ports = pair_data(n)
        sample[str(n)] = dict(aliquot=s, proper_divisors=ds[:-1], cubic_ports=ports,
                             abundance=s-n, square_total=isqrt(n+s)**2 == n+s)
    # Check the inherited doubling law, and the new one-hit equation, independently.
    tower_cases = 0
    for m in range(1, 10_000, 2):
        hits = []
        for a in range(8):
            n = m << a
            if n > limit:
                break
            defect = sig[n] - 2*n
            require(defect == (1 << (a+1)) * (sig[m]-m) - sig[m], ('defect carry', n))
            if defect == 2:
                hits.append(a)
            tower_cases += 1
        require(len(hits) <= 1, ('one-hit', m))
    # Pin the already-proved complete fiber of 343, using its proved cutoff.
    fiber = [n for n in range(2, 342**2+1) if sig[n]-n == 343]
    three_prime = [n for n in fiber if len(factor(n)) == 3 and all(a == 1 for _, a in factor(n))]
    require(len(fiber) == 19 and three_prime == [1001], 'inherited anchor scope')
    orbits = []
    for seed in (1001, 1000001):
        chain = [seed]
        while chain[-1] > 0:
            require(len(chain) < 100, 'unexpected unresolved orbit')
            n = chain[-1]
            chain.append(sigma(n) - n)
        orbits.append(chain)
    require(orbits[0] == [1001,343,57,23,1,0], 'anchor orbit')
    require(orbits[1] == [1000001,10003,1437,483,285,195,141,51,21,11,1,0], 'decimal orbit')
    require(1001+1001000 == 1001**2 and sigma(1001000)-1001000 == 2143960,
            'square-sum versus aliquot edge')
    # Local scalar compression does not preserve named square-sum ports.
    require(sigma(27)-27 == sigma(35)-35 == 13, 'same observation')
    require(isqrt(9+27)**2 == 9+27 and isqrt(9+35)**2 != 9+35, 'port hostile')
    # Compression at 20 does not persist under native aliquot or doubling steps.
    require(pair_data(20)[2] == [[1,20]] and pair_data(40)[2] == [] and
            pair_data(18)[2] == [] and pair_data(16)[2] == [], 'class hostiles')
    return dict(status='PASS', full_integer_universe=[1,limit], classification_counts=dict(counts),
                direct_pair_controls=len(controls), tower_cases=tower_cases, p2qr_controls=triples, exceptions=exceptions,
                samples=sample, complete_343_fiber=fiber, orbits=orbits,
                scope='Finite validation supplements the elementary proofs; no all-odd exclusion or infinitude claim.')


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
