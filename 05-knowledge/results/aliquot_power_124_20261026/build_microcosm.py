#!/usr/bin/env python3
"""Deterministic comparison of declared tables; no source-cutoff conjecture."""
import json
from math import gcd, isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent


def factor(n):
    result = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = 1
    return result


def odd(n):
    while n % 2 == 0:
        n //= 2
    return n


def aliquot(n):
    return sum(d for d in range(1, n) if n % d == 0)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def first_power(n):
    exponents = factor(n)
    g = 0
    for e in exponents.values():
        g = gcd(g, e)
    return __import__('math').prod(p ** (e // g) for p, e in exponents.items())


def combined_addresses(value):
    addresses = []
    for r in range(2, value + 1):
        power, k = r, 1
        while power <= value:
            v, j = power, 0
            while v <= value:
                if v == value:
                    addresses.append([r, k, j])
                v *= 2
                j += 1
            power *= r
            k += 1
    return sorted(addresses)


def main():
    inherited = json.loads((HERE / 'aliquot/results.normal.json').read_text())['fibres']
    # Completeness follows from the proved bound n<=(m-1)^2, m>1.
    independent = {m: [] for m in range(2, 33)}
    for source in range(1, 31 ** 2 + 1):
        target = aliquot(source)
        if target in independent:
            independent[target].append(source)
    rows = []
    for n in range(1, 33):
        if n >= 2:
            require(independent[n] == inherited[str(n)], 'complete fibre cross-check')
        addresses = combined_addresses(n * n) if n > 1 else []
        first_combined = addresses[0][0] if addresses else None
        rows.append({
            'n': n, 's_n': aliquot(n),
            'inverse_fibre': independent[n] if n >= 2 else 'all primes',
            'doubling_first_row_n': odd(n),
            'n_already_in_strictly_earlier_doubling_row': odd(n) < n,
            'power_first_base_n': first_power(n) if n > 1 else None,
            'square': n * n, 'square_doubling_first_row': odd(n) ** 2,
            'square_row_available_by_base_n': odd(n) ** 2 <= n,
            'cube_doubling_first_row': odd(n) ** 3,
            'square_combined_first_base': first_combined,
            'square_combined_address_count': len(addresses),
            'square_combined_addresses': addresses,
            'perfect_fixed_point': aliquot(n) == n,
            'prime_square_preimages': [x for x in independent.get(n, [])
                                      if isqrt(x) ** 2 == x and
                                      factor(isqrt(x)) == {isqrt(x): 1}],
        })
    require(rows[5]['inverse_fibre'] == [6, 25], 'six fibre')
    require(rows[5]['square_doubling_first_row'] == 9, '36 row9')
    require(rows[5]['square_combined_first_base'] == 3, '36 combined base3')
    require(rows[27]['square_doubling_first_row'] == 49, '28 perfect diagram')
    require(rows[11]['square_doubling_first_row'] == 9 and
            rows[11]['cube_doubling_first_row'] == 27, '12 exponent boundary')
    # Literal face incidence on three prime coordinates: triangle boundary.
    proper_faces = [mask for mask in range(8) if mask != 7]
    prime_weights = [2, 3, 5]
    evaluation = sum(__import__('math').prod(prime_weights[i]
                     for i in range(3) if mask >> i & 1) for mask in proper_faces)
    require(evaluation == 42 == aliquot(30), 'triangle boundary weighted at30')
    require([sum(mask.bit_count() == k for mask in proper_faces)
             for k in range(4)] == [1, 3, 3, 0], 'triangle face profile')
    print(json.dumps({'status': 'PASS', 'base_universe': [1, 32],
        'fibre_universe': {'targets': [2, 32], 'proved_source_bound': 961},
        'zero_fibre': [1], 'first_row_order': 'minimum numerical base',
        'power_and_combined_base_domain': 'base>=2; exponent>=1; doubling height>=0',
        'shape_control': {'n': 30, 'prime_weights': prime_weights,
                          'face_profile': [1, 3, 3, 0], 'weighted_sum': evaluation},
        'rows': rows}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
