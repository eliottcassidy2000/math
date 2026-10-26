"""Exact, bounded controls for the Busy Beaver / ternary / divisor atlas.

No imported research engines; all gates survive python -O.
The cylinder claims are proved in README.md, not extrapolated from this test.
"""
from fractions import Fraction
import itertools
import json
from math import prod


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def raw_next(x):
    r = x % 3
    return None if r == 2 else (5 * x + (18 if r == 0 else 22)) // 3


def next_z(z):
    return None if z % 3 == 0 else (5 * z + (1 if z % 3 == 1 else -1)) // 3


def prefix(z, depth):
    signs = []
    for _ in range(depth):
        if z % 3 == 0:
            return None
        signs.append(1 if z % 3 == 1 else -1)
        z = (5 * z + signs[-1]) // 3
    return z, tuple(signs)


def predecessor(y):
    if y % 5 == 2:
        return (3 * y - 1) // 5
    if y % 5 == 3:
        return (3 * y + 1) // 5
    return None


def root(y):
    while predecessor(y) is not None:
        p = predecessor(y)
        need(0 < p < y and next_z(p) == y, 'backward path')
        y = p
    return y


def compose(signs):
    a = 0
    denominator = 1
    for e in signs:
        a = 5 * a + e * denominator
        denominator *= 3
    return a


def ternary(n, width=0):
    digits = []
    while n or len(digits) < max(1, width):
        digits.append(n % 3)
        n //= 3
    return tuple(reversed(digits))


def factors(n):
    result = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def sigma(n):
    return prod(sum(p ** i for i in range(a + 1)) for p, a in factors(n).items())


def macro_control():
    expected = [0, 6, 16, 34, 64, 114, 196, 334, 564, 946, 1584,
                2646, 4416, 7366, 12284]
    x, values, clocks = 0, [], []
    while True:
        values.append(x)
        k, r = divmod(x, 3)
        clocks.append(5 * k * k + 19 * k + 15 if r == 0 else
                      5 * k * k + 25 * k + 27 if r == 1 else 6 * k + 12)
        y = raw_next(x)
        if y is None:
            break
        need(x % 2 == y % 2 == 0, 'even invariant')
        need(next_z(x // 2 + 5) == y // 2 + 5, 'affine normalization')
        x = y
    need(values == expected and sum(clocks) == 47176870, 'source orbit or macro clock')
    need(values[-1] // 3 + 4 == 4098, 'halting tape ones')
    for x in range(0, 10000, 2):
        y, w = raw_next(x), next_z(x // 2 + 5)
        need((y is None and w is None) or (y is not None and w == y // 2 + 5),
             'normalization off the selected source')
    return {'raw_orbit': values, 'normalized_orbit': [x // 2 + 5 for x in values],
            'macro_steps_before_halt': len(values) - 1, 'macro_clock_costs': clocks,
            'cited_clock_sum': sum(clocks), 'cited_halting_ones': 4098,
            'clock_scope': 'Exact sum of Michel cited transition formulas, not a new full Turing simulation.'}


def cylinder_control():
    data = []
    for depth in range(11):
        m = 3 ** depth
        algebraic = {}
        for signs in itertools.product((-1, 1), repeat=depth):
            a = compose(signs)
            residue = (-a * pow(5 ** depth, -1, m)) % m if m > 1 else 0
            need(residue not in algebraic, 'duplicate itinerary cylinder')
            algebraic[residue] = signs
            out = prefix(residue, depth)
            need(out is not None and out[1] == signs, 'algebraic guard')
        direct = {z for z in range(m) if prefix(z, depth) is not None}
        need(direct == set(algebraic) and len(direct) == 2 ** depth, 'complete cylinder count')
        for z in direct:
            targets = [prefix(z + t * m, depth)[0] for t in range(3)]
            need(len({w % 3 for w in targets}) == 3, 'next-digit bijection')
            need(targets[1] - targets[0] == targets[2] - targets[1] == 5 ** depth,
                 'exact affine carry')
        data.append({'depth': depth, 'classes': len(direct), 'modulus': m})
    a = prefix(5, 14)
    b = prefix(5 + 3 ** 14, 14)
    c = prefix(5 + 3 ** 15, 14)
    need(a[1] == b[1] == c[1], 'source-prefix controls')
    need(a[0] % 3 == c[0] % 3 == 0 and b[0] % 3 != 0, 'halting precision control')
    need(root(5) == 5 and root(5 + 3 ** 14) != 5, 'distinct integer components')
    return {'levels': data, 'same_14_step_itinerary_sources': [5, 5 + 3 ** 14, 5 + 3 ** 15],
            'their_14_step_targets': [a[0], b[0], c[0]],
            'halting_on_next_test': [True, False, True],
            'their_integer_roots': [root(5), root(5 + 3 ** 14), root(5 + 3 ** 15)]}


def forest_and_periodic_control():
    for z in range(1, 20001):
        y = next_z(z)
        if y is not None:
            need(y > z and y % 5 in (2, 3) and predecessor(y) == z, 'integer forward forest')
        p = predecessor(z)
        if p is not None:
            need(0 < p < z and next_z(p) == z, 'exact image characterization')
    counts = []
    for length in range(1, 10):
        values = set()
        for signs in itertools.product((-1, 1), repeat=length):
            z = Fraction(-compose(signs), 5 ** length - 3 ** length)
            need(abs(z) <= Fraction(1, 2), 'periodic real bound')
            initial = z
            for e in signs:
                residue = (z.numerator * pow(z.denominator, -1, 3)) % 3
                need(residue == (1 if e == 1 else 2), 'rational periodic guard')
                z = (5 * z + e) / 3
            need(z == initial, 'periodic fixed-point identity')
            values.add(initial)
        need(len(values) == 2 ** length, 'periodic coding uniqueness')
        counts.append(len(values))
    need(root(27) == 16 and next_z(27) is None, 'literal 27 is not the Collatz source')
    return {'integer_range': [1, 20000], 'periodic_word_lengths': list(range(1, 10)),
            'distinct_periodic_points': counts, 'all_in_real_interval': ['-1/2', '1/2'],
            'literal_27_control': {'root': 16, 'halts_immediately': True}}


def sloane_control():
    levels = []
    parents = None
    for k in range(1, 11):
        m, period = 3 ** k, 2 * 3 ** (k - 1)
        unit, seen, good = 1, set(), set()
        for exponent in range(period):
            need(unit not in seen, 'primitive root period too short')
            seen.add(unit)
            if 0 not in ternary(unit, k):
                good.add(exponent)
            unit = 2 * unit % m
        need(unit == 1 and len(seen) == period and len(good) == 2 ** k, 'Sloane level')
        if parents is not None:
            old_period = period // 3
            for parent in parents:
                need(sum(parent + j * old_period in good for j in range(3)) == 2,
                     'two-of-three exponent lift')
        parents = good
        levels.append({'digits': k, 'exponent_period': period, 'zeroless_classes': len(good)})
    exceptions = [n for n in range(513) if 0 not in ternary(2 ** n)]
    need(exceptions == [0, 1, 2, 3, 4, 15], 'finite Sloane control')
    need(0 not in ternary(2 ** 15) and 0 in ternary(2 ** 15, len(ternary(2 ** 15)) + 1),
         'leading-zero padding hostile')
    return {'levels': levels, 'finite_exponent_range': [0, 512], 'zero_free_exponents': exceptions,
            'ternary_2_to_15': ''.join(map(str, ternary(2 ** 15))),
            'padding_a_leading_zero_is_not_a_digit_counterexample': True}


def divisor_controls():
    records = []
    for n, k in ((6, 2), (28, 2), (120, 3), (30240, 4)):
        f = factors(n)
        local = {p: factors(sum(p ** j for j in range(a + 1))) for p, a in f.items()}
        all_primes = set(factors(k)) | set(f)
        for row in local.values():
            all_primes |= set(row)
        need(sigma(n) == k * n, 'multiperfect control')
        need(all(sum(row.get(q, 0) for row in local.values())
                 == f.get(q, 0) + factors(k).get(q, 0) for q in all_primes), 'valuation balance')
        records.append({'n': n, 'multiplier': k, 'aliquot_image': sigma(n) - n,
                        'factorization': f, 'local_sigma_prime_flows': local})
    need(sigma(120) - 120 == 240 != 120, 'triperfect is not aliquot fixed')
    need(sigma(32445) - 2 * 32445 == 6, 'old odd near-perfect hostile')
    def prime(p):
        return p >= 2 and all(p % d for d in range(2, int(p ** 0.5) + 1))
    supports = {n: [d + 1 for d in range(1, n + 1) if n % d == 0 and prime(d + 1)]
                for n in (6, 114)}
    need(supports[6] == supports[114] == [2, 3, 7], 'Bernoulli denominator nonuniqueness')
    need(Fraction(5, 3) ** 2 - Fraction(5, 3) - 1 == Fraction(1, 9), 'phi versus 5/3')
    # A tempting cubic fit to 6,16,34,64 predicts 110; the real fifth term is114.
    need(64 + (30 + (12 + 4)) == 110 != 114, 'sequence-fit hostile')
    return {'multiperfect_controls': records, 'odd_near_perfect_32445_defect': 6,
            'von_staudt_supports': supports, 'golden_polynomial_at_5_over_3': '1/9',
            'first_four_terms_cubic_prediction': 110, 'actual_fifth_term': 114}


def main():
    result = {'status': 'PASS', 'macro': macro_control(), 'cylinders': cylinder_control(),
              'integer_forest_and_periodic': forest_and_periodic_control(),
              'sloane': sloane_control(), 'divisors_and_negative_controls': divisor_controls()}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
