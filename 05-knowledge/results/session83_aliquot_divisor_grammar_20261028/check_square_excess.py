#!/usr/bin/env python3
"""Exact certificate: Δ(n)=2 and square σ(n) imply ω(n)>=356.

Standard library only. Normal/-O identical; explicit checks are not assertions.
This does not prove existence or nonexistence of such n.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt, prod
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


def trial_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def sieve_primes(limit):
    prime = bytearray(b'\1') * (limit + 1)
    prime[0:2] = b'\0\0'
    for p in range(2, isqrt(limit) + 1):
        if prime[p]:
            prime[p*p:limit+1:p] = b'\0' * ((limit-p*p)//p + 1)
    return [p for p in range(2, limit+1) if prime[p]]


def certificate():
    cap = 5441
    sieve = sieve_primes(cap)
    trial = [n for n in range(2, cap+1) if trial_prime(n)]
    require(sieve == trial, 'Independent primality universes disagree')
    allowed = [p for p in sieve if p % 8 in (1, 7)]
    require(len(allowed) == 356 and allowed[-2:] == [5431, 5441], 'Prime-list boundary changed')
    ratio = Fraction(1)
    boundaries = []
    crossing = None
    for j, p in enumerate(allowed, 1):
        ratio *= Fraction(p, p-1)
        if crossing is None and ratio > 2:
            crossing = j
        if j in (355, 356):
            boundaries.append({
                'prime_count': j,
                'last_prime': p,
                'ratio_numerator': str(ratio.numerator),
                'ratio_denominator': str(ratio.denominator),
                'numerator_minus_twice_denominator': str(ratio.numerator-2*ratio.denominator),
            })
    require(crossing == 356, 'First abundancy-envelope crossing is not 356')
    r355 = Fraction(int(boundaries[0]['ratio_numerator']), int(boundaries[0]['ratio_denominator']))
    require(Fraction(1999836, 1000000) < r355 < Fraction(1999837, 1000000), '355 decimal enclosure')
    require(Fraction(2000203, 1000000) < ratio < Fraction(2000204, 1000000), '356 decimal enclosure')
    lower = prod(allowed)
    require(10**1155 <= lower < 10**1156, 'Primorial decimal size changed')
    # Independent direct quadratic-residue control of the supplementary law.
    for p in sieve:
        if p > 2:
            residues = {x*x % p for x in range(1, p)}
            require((2 in residues) == (p % 8 in (1, 7)), f'Residue law failed at {p}')
    # Hostile: satisfying the norm shape does not imply the divisor-sum equation.
    require(7 == 2*2**2-1 and sum(d for d in range(1,8) if 7 % d == 0) != 4*2**2,
            'Necessary norm shape was accidentally treated as sufficient')
    return {
        'status': 'FINITE-EXACT CERTIFICATE FOR A PROVED NECESSARY CONDITION',
        'statement': 'If sigma(n)=2*n+2 is a square, then n is odd and omega(n)>=356.',
        'not_claimed': ['existence', 'nonexistence', 'sharpness for actual solutions'],
        'prime_universe': {'all_primes_through': cap, 'independent_sieve_and_trial': True},
        'allowed_primes': allowed,
        'allowed_prime_sha256': sha256(','.join(map(str, allowed)).encode()).hexdigest(),
        'boundary_certificates': boundaries,
        'first_envelope_crossing': crossing,
        'source_lower_bound_decimal_digits': len(str(lower)),
        'source_lower_bound_exact': str(lower),
        'quadratic_residue_controls': len(sieve)-1,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('square_excess_results.json'))
    args = parser.parse_args()
    result = certificate()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('Verified: first allowed-prime envelope crossing is 356; every candidate has >=356 distinct primes and >=1156 decimal digits.')


if __name__ == '__main__':
    main()
