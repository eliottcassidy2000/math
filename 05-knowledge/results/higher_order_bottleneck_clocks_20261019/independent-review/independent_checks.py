#!/usr/bin/env python3
"""Exact-integer-grid escape certificates, independent of Decimal power/rounding.

Only standard-library integers and Fractions are used for certification.
The true iterate lies in [lower / D, upper / D], with D = 2**192.
At each step the rational polynomial at each endpoint is evaluated exactly,
then its scaled value is rounded outwards with integer floor/ceiling division.
Since 1 + f'(y) > 0 on [-1/8, 1/8], induction preserves the enclosure.
Both endpoint iterates must be below the threshold before every step, and
the final lower endpoint must be at or above it. Ambiguity raises an error.
"""

from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


D = 1 << 192
A = D // 8
J_RIGHT = 5 * D // 32


def escape(r, power, lam):
    need(r >= 2 and r % 2 == 0, 'even degree required')
    lam = F(lam)
    need(0 <= lam <= 1 and (r == 4 or lam == 0), 'unsupported deformation')
    need(3*r <= power <= 192, 'epsilon outside theorem/grid range')
    m = max(r, 6) if lam else r
    denominator = lam.denominator * D**(m-1)
    coeff_r = lam.denominator * D**(m-r)
    coeff_6 = lam.numerator * D**(m-6) if lam else 0
    epsilon_scaled = 1 << (192-power)

    def numerator(x):
        return ((x + epsilon_scaled) * denominator
                + coeff_r * x**r + (coeff_6*x**6 if lam else 0))

    lower = upper = -A
    n = 0
    while True:
        need(-A <= lower <= upper < A, 'pre-exit enclosure outside interval')
        last_upper = upper
        old_lower = lower
        lower = numerator(lower) // denominator
        upper = -((-numerator(upper)) // denominator)
        need(old_lower < lower <= upper <= J_RIGHT, 'step/enclosure failure')
        n += 1
        if lower >= A:
            return {
                'r': r,
                'epsilon': f'2^(-{power})',
                'lambda': str(lam),
                'certified_escape_count': n,
                'previous_upper_threshold_gap_grid_units': str(A-last_upper),
                'final_lower_threshold_gap_grid_units': str(lower-A),
                'final_interval_width_grid_units': str(upper-lower),
            }
        need(upper < A, 'threshold ambiguous: increase grid precision')
        need(n < 1000000, 'resource guard: count not certified')


def main():
    a, c = F(1,8), F(5,32)
    clock = F(1,3) + (F(8,3)+F(1,80))/F(11,16)
    clock += F(1,2) + (2+F(1,64))/F(11,16)
    need(clock == F(5057,660) < 8, 'clock coefficient')
    # The following ratios bound the r -> r+2 induction for every r >= 2.
    need(2*c*c < 1, 'first derivative induction ratio')
    need(6*c*c < 1, 'second derivative induction ratio')
    need(4*a*a < 1, 'squared derivative induction ratio')
    need(2*c == F(5,16) and 2 == 2 and 4*a**0 == 4, 'r=2 bases')
    need(2*a**4+a**6 <= min(a/4,F(1,32)), 'quartic step')
    need(4*c**3+6*c**5 < F(5,16), 'quartic first derivative')
    need(12*c*c+30*c**4 < 2, 'quartic second derivative')
    need((4+6*a*a)**2*a*a < 4, 'quartic squared derivative')

    author_dir = Path(__file__).resolve().parent.parent
    author = json.loads((author_dir / 'verification.json').read_text())
    cases = [(2,10,0),(2,14,0),(4,12,0),(4,16,0),(4,20,0),
             (4,16,F(1,2)),(4,20,F(1,2)),(4,16,1),(4,20,1),(6,18,0)]
    rows = [escape(*case) for case in cases]
    need(len(rows) == len(author['rows']), 'receipt length mismatch')
    for row, reference in zip(rows, author['rows']):
        for key in ('r', 'epsilon', 'certified_escape_count'):
            need(row[key] == reference[key], f'author mismatch: {key}')
        need(F(row['lambda']) == F(reference['lambda']), 'lambda mismatch')
    extras = [escape(2,6,0), escape(4,12,1), escape(4,16,F(1,3)), escape(4,18,1)]
    source_hashes = {name: hashlib.sha256((author_dir/name).read_bytes()).hexdigest()
                     for name in ['HIGHER_ORDER_BOTTLENECKS.md',
                                  'verify_higher_bottlenecks.py', 'verification.json']}
    result = {
        'status': 'PASS',
        'method': 'Exact integer polynomial arithmetic and outward integer division',
        'grid_denominator': str(D),
        'clock_constant': str(clock),
        'author_counts_matched': len(rows),
        'source_sha256_at_review': source_hashes,
        'matched_rows': rows,
        'additional_parameter_checks': extras,
        'scope': 'Finite numerical controls; the written proof supplies uniform parameter claims.',
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
