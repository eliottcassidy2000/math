"""Exact Euler-product decoding and level-11 local arithmetic. Stdlib only.

No modularity theorem is inferred from a finite match. Print deterministic JSON;
--write saves the companion .out. All controls survive python -O.
"""
from collections import Counter
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import argparse
import json

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(label)


def euler_product(exponents):
    bound = len(exponents) - 1
    out = [1] + [0] * bound
    for n, multiplicity in enumerate(exponents[1:], 1):
        for _ in range(abs(multiplicity)):
            if multiplicity >= 0:
                for k in range(bound, n - 1, -1):
                    out[k] -= out[k - n]
            else:
                for k in range(n, bound + 1):
                    out[k] += out[k - n]
    return out


def recover_exponents(coefficients):
    """F*c=-qF', then divisor-triangular inversion; integers throughout."""
    bound = len(coefficients) - 1
    c = [0] * (bound + 1)
    b = [0] * (bound + 1)
    for n in range(1, bound + 1):
        c[n] = -n * coefficients[n] - sum(
            coefficients[i] * c[n - i] for i in range(1, n))
        numerator = c[n] - sum(d * b[d] for d in range(1, n) if n % d == 0)
        check(numerator % n == 0, "Euler multiplicity integrality")
        b[n] = numerator // n
    return b, c


def coefficients_from_divisors(b):
    bound = len(b) - 1
    c = [sum(d * b[d] for d in range(1, n + 1) if n % d == 0)
         for n in range(bound + 1)]
    out = [1] + [0] * bound
    for n in range(1, bound + 1):
        numerator = -sum(c[k] * out[n - k] for k in range(1, n + 1))
        check(numerator % n == 0, "independent coefficient integrality")
        out[n] = numerator // n
    return out


class Field:
    """Small F_p[X]/(monic polynomial), base-p integer representation."""
    def __init__(self, p, polynomial):
        self.p, self.modulus = p, polynomial
        self.degree = len(polynomial) - 1
        self.size = p ** self.degree

    def digits(self, a):
        return [(a // self.p ** j) % self.p for j in range(self.degree)]

    def encode(self, digits):
        return sum((a % self.p) * self.p ** j for j, a in enumerate(digits))

    def add(self, a, b):
        return self.encode([x + y for x, y in zip(self.digits(a), self.digits(b))])

    def mul(self, a, b):
        if self.degree == 1:
            return a * b % self.p
        aa, bb = self.digits(a), self.digits(b)
        product = [0] * (2 * self.degree - 1)
        for i, x in enumerate(aa):
            for j, y in enumerate(bb):
                product[i + j] += x * y
        for k in range(len(product) - 1, self.degree - 1, -1):
            lead = product[k] % self.p
            for j, coefficient in enumerate(self.modulus):
                product[k - self.degree + j] -= lead * coefficient
        return self.encode(product[:self.degree])

    def power(self, a, n):
        out = 1
        while n:
            if n & 1:
                out = self.mul(out, a)
            a = self.mul(a, a)
            n //= 2
        return out


def count_curve(field):
    # y^2+y=x^3-x^2-10x-20; add the single point at infinity.
    left = Counter(field.add(field.mul(y, y), y) for y in range(field.size))
    total = 1
    for x in range(field.size):
        square = field.mul(x, x)
        right = field.add(field.mul(square, x), field.mul((-1) % field.p, square))
        right = field.add(right, field.mul((-10) % field.p, x))
        right = field.add(right, (-20) % field.p)
        total += left[right]
    return total


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def frobenius_traces(p, a, length):
    traces = [2, a]
    for _ in range(2, length + 1):
        traces.append(a * traces[-1] - p * traces[-2])
    return traces


def add_points(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    xx, yy = right
    if x == xx and y + yy == -1:
        return None
    slope = ((3 * x * x - 2 * x - 10) / (2 * y + 1)
             if left == right else (yy - y) / (xx - x))
    intercept = y - slope * x
    xxx = slope * slope + 1 - x - xx
    return xxx, -slope * xxx - intercept - 1


def native_word_series(bound):
    """Two independent counts of no-LB words by dyadic denominator cost."""
    costs = {'H': 10, 'G': 3, 'A': 7, 'B': 7, 'L': 4}
    counts = [[0, 0] for _ in range(bound + 1)]
    counts[0][0] = 1  # previous letter is not L (including the empty word).
    for cost in range(bound + 1):
        for previous_l in range(2):
            for letter, weight in costs.items():
                if cost + weight <= bound and not (previous_l and letter == 'B'):
                    counts[cost + weight][int(letter == 'L')] += counts[cost][previous_l]
    total = [sum(row) for row in counts]
    denominator = {0: 1, 3: -1, 4: -1, 7: -2, 10: -1, 11: 1}
    reciprocal = [1] + [0] * bound
    for n in range(1, bound + 1):
        reciprocal[n] = -sum(v * reciprocal[n - k] for k, v in denominator.items()
                             if 0 < k <= n)
    check(total == reciprocal, "automaton vs rational-series word census")
    exponents, _ = recover_exponents(total)
    check(euler_product(exponents) == total, "history census Euler inversion")
    return total, exponents


def run():
    bound = 256
    b = [0] + [2 + 2 * (n % 11 == 0) for n in range(1, bound + 1)]
    coefficients = euler_product(b)
    check(coefficients == coefficients_from_divisors(b), "two expansion paths")
    decoded, logarithmic = recover_exponents(coefficients)
    check(decoded == b, "lossless Euler exponent decoding")
    a = [0] + coefficients  # a[n] is the coefficient of q^n, through 257.
    check(a[1:20] == [1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1,
                      -2, 4, 4, -1, -4, -2, 4, 0], "initial series")
    # General integer exponents, including inverse factors: not eta-specific.
    mixed = [0] + [((7 * n) % 9) - 4 for n in range(1, 65)]
    check(recover_exponents(euler_product(mixed))[0] == mixed, "signed exponent control")
    altered = b.copy()
    altered[65] += 1
    changed = euler_product(altered)
    check(changed[:65] == coefficients[:65] and changed[65] != coefficients[65],
          "finite-prefix does not determine unseen factors")
    for m in range(1, bound + 1):
        for n in range(1, bound // m + 1):
            if gcd(m, n) == 1:
                check(a[m * n] == a[m] * a[n], "coprime multiplicativity")
    prime_counts = {}
    for p in range(2, 100):
        if not prime(p) or p == 11:
            continue
        brute = 1 + sum((y * y + y - x ** 3 + x * x + 10 * x + 20) % p == 0
                        for x in range(p) for y in range(p))
        check(brute == count_curve(Field(p, [0, 1])), "two point-count paths")
        check(a[p] == p + 1 - brute, "good-prime modularity finite control")
        prime_counts[str(p)] = brute
        power, previous, current = p, 1, a[p]
        while power * p <= bound:
            nxt = a[p] * current - p * previous
            check(a[power * p] == nxt, "Hecke recurrence finite control")
            previous, current, power = current, nxt, power * p
    check(a[11] == a[121] == 1, "bad-prime coefficient control")
    fields = [(2, 1, [0, 1]), (2, 2, [1, 1, 1]), (2, 3, [1, 1, 0, 1]),
              (2, 4, [1, 1, 0, 0, 1]), (3, 1, [0, 1]), (3, 2, [2, 2, 1])]
    extension_counts = {}
    for p, degree, polynomial in fields:
        field = Field(p, polynomial)
        for x in range(1, field.size):
            check(field.power(x, field.size - 1) == 1, "quotient really a field")
        count = count_curve(field)
        trace = frobenius_traces(p, a[p], degree)[degree]
        check(count == field.size + 1 - trace, "extension field direct count")
        extension_counts[str(field.size)] = count
    check(extension_counts == {'2': 5, '4': 5, '8': 5, '16': 25, '3': 5, '9': 15},
          "small field point counts")
    golden = Field(3, [2, 2, 1])
    check(golden.power(3, 8) == 1 and golden.power(3, 4) == 2,
          "phi has order8 in F9")
    check(a[9] == -2 and frobenius_traces(3, a[3], 2)[2] == -5,
          "Hecke coefficient is not the Frobenius power trace")
    generator = (F(5), F(5))
    multiples, point = [], None
    for _ in range(5):
        point = add_points(point, generator)
        multiples.append(None if point is None else [str(v) for v in point])
        if point is not None:
            x, y = point
            check(y * y + y == x ** 3 - x * x - 10 * x - 20, "rational curve point")
    check(multiples == [['5', '5'], ['16', '-61'], ['16', '60'], ['5', '-6'], None],
          "exact order5 point")
    word_counts, word_exponents = native_word_series(64)
    check(word_counts[1] == 0 != coefficients[1], "controller census differs from eta product")
    return {'status': 'PASS', 'checks': CHECKS, 'product_degree': bound + 1,
            'coefficients_a1_to_a32': a[1:33], 'euler_multiplicities_b1_to_b33': b[1:34],
            'logarithmic_coefficients_c1_to_c11': logarithmic[1:12],
            'good_prime_point_counts_below100': prime_counts,
            'extension_field_point_counts': extension_counts,
            'five_torsion_multiples': multiples,
            'local_factors_p2_p3_p11': ['1+2T+2T^2', '1+T+3T^2', '1-T'],
            'native_word_series_denominator': '1-z^3-z^4-2z^7-z^10+z^11',
            'native_word_counts_cost0_to24': word_counts[:25],
            'native_word_euler_exponents_b1_to24': word_exponents[1:25],
            'trace_vs_hecke_at3': {'a9': a[9], 'trace_F9': -5, 'points_F9': 15},
            'scope': 'Finite exact controls; formal Euler decoding proved in note; modularity and general Frobenius facts cited; no Collatz payment transfer.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, sort_keys=True) + '\n'
    if args.write:
        root = Path(__file__).resolve().parents[2]
        (root / '05-knowledge/results/eta11_lossless_coordinates_20261005.out').write_text(
            result, encoding='utf-8', newline='\n')
    print(result, end='')
