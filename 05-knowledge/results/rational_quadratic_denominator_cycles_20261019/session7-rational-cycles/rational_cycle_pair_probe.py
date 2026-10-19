#!/usr/bin/env python3
"""Complete search for rational cycles of length >=3 at a given denominator.

The theorem assumes a certified prime factorization is supplied. The CLI's
--factorization option has this same precondition; it checks the product and
basic syntax, not primality. The bounded --max-b census uses trial division.
There is no bound on the parameter numerator. All orbit arithmetic has an
explicit real-escape guard. No assertion is used as an acceptance condition.
"""
import argparse
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, prod
import json


def factor(n):
    result = []
    p = 2
    while p*p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            result.append((p, e))
        p += 1
    if n > 1:
        result.append((n, 1))
    return result


def search(b, certified_prime_factors):
    fs = sorted(certified_prime_factors)
    if (b % 4 or len({p for p, _ in fs}) != len(fs)
            or any(p < 2 or e < 1 for p, e in fs)
            or prod(p**e for p, e in fs) != b):
        raise ValueError('Expected b divisible by four and its prime factorization')
    N = 2**len(fs)+2
    pairs = [(0, k) for k in range(N)]+[(k, k-1) for k in range(1, N)]
    h = b//2
    magnitudes = [1]
    local = []
    for p, e in fs:
        powers = [(p**(e-(p == 2)))**k for k in range(N)]
        magnitudes = [a*z for a in magnitudes for z in powers]
        local.append([(powers[i], powers[j]) for i, j in pairs])
    smooth = set(magnitudes)
    found = set()
    primitive = escaped = unsigned = 0
    for parts in product(*local):
        A = prod(a for a, _ in parts)
        B = prod(v for _, v in parts)
        hv_over_u, remainder = divmod(h*B, A)
        if remainder:
            raise ArithmeticError('Universal pair divisibility lemma failed')
        unsigned += 1
        for su, sv in product((1, -1), repeat=2):
            u, v = su*A, sv*B
            w = v-h+su*sv*hv_over_u
            if abs(w) not in smooth:
                continue
            a0, a1, a2 = u+v-w, u+w-v, v+w-u
            if len({a0, a1, a2}) != 3 or any(gcd(a, b) != 1 for a in (a0, a1, a2)):
                continue
            C = b*a1-a0*a0
            if a1*a1+C != b*a2 or gcd(C, b) != 1:
                raise ArithmeticError('Reconstruction identity failed')
            primitive += 1
            escape_bound = b+isqrt(abs(C))+1
            orbit, seen, a = [], set(), a0
            for step in range(N+1):
                if a in seen:
                    if a == a0 and len(orbit) >= 3:
                        canonical = min(tuple(orbit[k:]+orbit[:k]) for k in range(len(orbit)))
                        found.add((C, canonical))
                    break
                if abs(a) > escape_bound:
                    escaped += 1
                    break
                if step == N:
                    break
                seen.add(a)
                orbit.append(a)
                a, rem = divmod(a*a+C, b)
                if rem or gcd(a, b) != 1:
                    break
    return dict(b=b, prime_factorization=fs, omega=len(fs), period_bound=N,
                local_pair_count=len(pairs), unsigned_pair_count=unsigned,
                signed_seed_count=4*unsigned, primitive_triangle_candidates=primitive,
                escape_rejections=escaped,
                cycles=[dict(C=C, c=str(Fraction(C, b*b)), numerators=list(a), period=len(a))
                        for C, a in sorted(found)])


def main():
    p = argparse.ArgumentParser()
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--max-b', type=int)
    mode.add_argument('--factorization', nargs='+', help='Certified factors, e.g. 2:2 3:20 5:20')
    p.add_argument('--output', required=True)
    p.add_argument('--compare', help='Primary exact census JSON to compare')
    args = p.parse_args()
    if args.factorization:
        fs = [tuple(map(int, item.split(':'))) for item in args.factorization]
        inputs = [(prod(q**e for q, e in fs), fs)]
    else:
        inputs = [(b, factor(b)) for b in range(4, args.max_b+1, 4) if len(factor(b)) >= 3]
    results = [search(b, fs) for b, fs in inputs]
    out = dict(method='Universal two-edge enumeration with certified factorization and real escape guard',
               results=results)
    if args.compare:
        expected = {r['b']: r['cycles'] for r in json.load(open(args.compare))['results']}
        actual = {r['b']: r['cycles'] for r in results}
        if expected != actual:
            raise ArithmeticError('Complete denominator/cycle comparison failed')
        out['comparison'] = 'PASS: all denominators and exact oriented cycles match'
    with open(args.output, 'w') as f:
        json.dump(out, f, indent=2)
        f.write('\n')
    print(json.dumps(dict(denominator_count=len(results), cycle_count=sum(len(r['cycles']) for r in results),
                          periods=sorted({c['period'] for r in results for c in r['cycles']}))))


if __name__ == '__main__':
    main()
