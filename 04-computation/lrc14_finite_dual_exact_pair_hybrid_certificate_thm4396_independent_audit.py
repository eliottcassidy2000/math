#!/usr/bin/env python3
"""Independent rigorous audit of the THM-4396 two-Fejer/exact-pair certificate.

No producer or repository implementation is imported.  Geometry is exact
Fraction arithmetic; transcendental quantities are enclosed by rational
intervals derived from Machin's identity and alternating Taylor bounds.
Every theorem check is explicit and remains active under python -O.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from math import factorial
import json
import sys


LAM = Fraction(1, 14)
RAW_RADIUS = Fraction(3, 14)
CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


class Interval:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = Fraction(lo)
        self.hi = self.lo if hi is None else Fraction(hi)
        if self.lo > self.hi:
            raise RuntimeError(f"invalid interval [{self.lo},{self.hi}]")

    def __repr__(self):
        return f"Interval({self.lo!r},{self.hi!r})"


def iv_add(x, y):
    return Interval(x.lo + y.lo, x.hi + y.hi)


def iv_neg(x):
    return Interval(-x.hi, -x.lo)


def iv_sub(x, y):
    return iv_add(x, iv_neg(y))


def iv_mul(x, y):
    corners = (x.lo * y.lo, x.lo * y.hi, x.hi * y.lo, x.hi * y.hi)
    return Interval(min(corners), max(corners))


def iv_scale(x, q):
    return iv_mul(x, Interval(q))


def iv_reciprocal(x):
    if x.lo <= 0 <= x.hi:
        raise RuntimeError(f"division interval meets zero: {x}")
    values = (1 / x.lo, 1 / x.hi)
    return Interval(min(values), max(values))


def iv_div(x, y):
    return iv_mul(x, iv_reciprocal(y))


def iv_sum(values):
    result = Interval(0)
    for value in values:
        result = iv_add(result, value)
    return result


def iv_square(x):
    if x.lo <= 0 <= x.hi:
        return Interval(0, max(x.lo * x.lo, x.hi * x.hi))
    values = (x.lo * x.lo, x.hi * x.hi)
    return Interval(min(values), max(values))


def arctan_alternating(x, last_index):
    """Consecutive alternating partial sums rigorously bracket atan(x)."""
    check(0 < x < 1, f"arctan input {x}")
    partial = Fraction(0)
    snapshots = []
    for k in range(last_index + 2):
        term = x ** (2 * k + 1) / (2 * k + 1)
        partial += term if k % 2 == 0 else -term
        if k in (last_index, last_index + 1):
            snapshots.append(partial)
    return Interval(min(snapshots), max(snapshots))


ATAN_FIFTH = arctan_alternating(Fraction(1, 5), 50)
ATAN_239 = arctan_alternating(Fraction(1, 239), 12)
PI_FINE = iv_sub(iv_scale(ATAN_FIFTH, 16), iv_scale(ATAN_239, 4))
PI_SCALE = 10**24
PI_INTERVAL = Interval(
    Fraction((PI_FINE.lo * PI_SCALE).numerator // (PI_FINE.lo * PI_SCALE).denominator, PI_SCALE),
    Fraction(-((-PI_FINE.hi * PI_SCALE).numerator // (-PI_FINE.hi * PI_SCALE).denominator), PI_SCALE),
)
check(PI_INTERVAL.lo > Fraction(103993, 33102), "Machin pi lower improves convergent")
check(PI_INTERVAL.hi < Fraction(104348, 33215), "Machin pi upper improves convergent")
check(PI_INTERVAL.lo <= PI_FINE.lo <= PI_FINE.hi <= PI_INTERVAL.hi,
      "outward fixed-denominator pi enclosure")
check(PI_INTERVAL.hi - PI_INTERVAL.lo <= Fraction(2, PI_SCALE), "Machin pi interval width")


def reduce_mod_two(value):
    quotient = value // 2
    result = value - 2 * quotient
    if result > 1:
        result -= 2
    check(Fraction(-1) < result <= 1, f"mod-two reduction {value}, {result}")
    return result


@lru_cache(maxsize=None)
def sin_small_pi(x):
    """Enclose sin(pi*x) for rational 0<=x<=1/2."""
    check(0 <= x <= Fraction(1, 2), f"small sine argument {x}")
    if x == 0:
        return Interval(0)
    z = iv_scale(PI_INTERVAL, x)
    check(z.lo >= 0 and z.hi * z.hi < 3, f"alternating sine range {x}")
    total = Interval(0)
    last_k = 14
    for k in range(last_k + 1):
        power = Interval(z.lo ** (2 * k + 1), z.hi ** (2 * k + 1))
        term = iv_scale(power, Fraction(1, factorial(2 * k + 1)))
        total = iv_add(total, term if k % 2 == 0 else iv_neg(term))
    next_power = 2 * (last_k + 1) + 1
    remainder = z.hi ** next_power / factorial(next_power)
    return Interval(total.lo - remainder, total.hi + remainder)


@lru_cache(maxsize=None)
def cos_small_pi(x):
    """Enclose cos(pi*x) for rational 0<=x<=1/2."""
    check(0 <= x <= Fraction(1, 2), f"small cosine argument {x}")
    if x == 0:
        return Interval(1)
    if x == Fraction(1, 2):
        return Interval(0)
    z = iv_scale(PI_INTERVAL, x)
    check(z.lo >= 0 and z.hi * z.hi < 3, f"alternating cosine range {x}")
    total = Interval(0)
    last_k = 14
    for k in range(last_k + 1):
        power = Interval(z.lo ** (2 * k), z.hi ** (2 * k))
        term = iv_scale(power, Fraction(1, factorial(2 * k)))
        total = iv_add(total, term if k % 2 == 0 else iv_neg(term))
    next_power = 2 * (last_k + 1)
    remainder = z.hi ** next_power / factorial(next_power)
    return Interval(total.lo - remainder, total.hi + remainder)


@lru_cache(maxsize=None)
def sin_pi(value):
    x = reduce_mod_two(Fraction(value))
    sign = 1
    if x < 0:
        x = -x
        sign = -1
    if x > Fraction(1, 2):
        x = 1 - x
    result = sin_small_pi(x)
    return result if sign > 0 else iv_neg(result)


@lru_cache(maxsize=None)
def cos_pi(value):
    x = abs(reduce_mod_two(Fraction(value)))
    sign = 1
    if x > Fraction(1, 2):
        x = 1 - x
        sign = -1
    result = cos_small_pi(x)
    return result if sign > 0 else iv_neg(result)


@lru_cache(maxsize=None)
def hhat(n):
    if n == 0:
        return Interval(Fraction(1, 7))
    magnitude = abs(n)
    if magnitude % 7 == 0:
        return Interval(0)
    numerator = sin_pi(Fraction(magnitude, 7))
    denominator = iv_scale(PI_INTERVAL, magnitude)
    return iv_div(numerator, denominator)


def tau(n, degree):
    if abs(n) > degree:
        return Fraction(0)
    return Fraction(degree + 1 - abs(n), degree + 1)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def sheet_character(w, n):
    residues = tuple((w[i] * n[i]) % 3 for i in range(3))
    check(sum(residues) % 3 == 0, f"resonance residues sum to zero {w}, {n}")
    if len(set(residues)) == 1:
        return 6
    if set(residues) == {0, 1, 2}:
        return -3
    raise RuntimeError(f"bad ternary character class {w}, {n}, {residues}")


def finite_spectrum(w, pair, degrees):
    i, j = pair
    k = next(index for index in range(3) if index not in pair)
    rows = []
    for ni in range(-degrees[0], degrees[0] + 1):
        for nj in range(-degrees[1], degrees[1] + 1):
            numerator = -(w[i] * ni + w[j] * nj)
            if numerator % w[k]:
                continue
            n = [0, 0, 0]
            n[i], n[j], n[k] = ni, nj, numerator // w[k]
            check(dot(w, n) == 0, f"finite resonance {w}, {n}")
            rows.append(tuple(n))
    return tuple(rows)


EXPECTED_SPECTRUM = (
    (-5, -1, 4), (-4, 6, -2), (-3, -4, 5), (-2, 3, -1),
    (-1, -7, 6), (0, 0, 0), (1, 7, -6), (2, -3, 1),
    (3, 4, -5), (4, -6, 2), (5, 1, -4),
)


def finite_M(w, pair, degrees):
    spectrum = finite_spectrum(w, pair, degrees)
    check(spectrum == EXPECTED_SPECTRUM, f"complete eleven-site spectrum {spectrum}")
    weights = {6: 0, -3: 0}
    killed = []
    terms = []
    for n in spectrum:
        character = sheet_character(w, n)
        weights[character] += 1
        factors = [hhat(n[index]) for index in range(3)]
        product_factor = iv_mul(iv_mul(factors[0], factors[1]), factors[2])
        term = iv_scale(
            product_factor,
            character * tau(n[pair[0]], degrees[0]) * tau(n[pair[1]], degrees[1]),
        )
        if term.lo == 0 and term.hi == 0:
            killed.append(n)
        terms.append(term)
    check(weights == {6: 3, -3: 8}, f"character counts {weights}")
    check(tuple(killed) == ((-1, -7, 6), (1, 7, -6)), f"exact sinc-zero sites {killed}")
    check(hhat(0).lo == hhat(0).hi == Fraction(1, 7), "zero mode survives")
    return iv_sum(terms), spectrum, weights, tuple(killed)


def sheet_intervals(speed, sheet):
    radius = LAM / speed
    rows = []
    for nearest in range(speed):
        center = Fraction(nearest, speed) - Fraction(sheet, 3)
        center -= center.numerator // center.denominator
        left, right = center - radius, center + radius
        if left < 0:
            rows.append((Fraction(0), right))
            rows.append((left + 1, Fraction(1)))
        elif right > 1:
            rows.append((left, Fraction(1)))
            rows.append((Fraction(0), right - 1))
        else:
            rows.append((left, right))
    rows.sort()
    for previous, current in zip(rows, rows[1:]):
        check(previous[1] <= current[0], f"disjoint sheet intervals {speed}, {sheet}")
    check(sum((right - left for left, right in rows), Fraction(0)) == Fraction(1, 7),
          f"sheet measure {speed}, {sheet}")
    return tuple(rows)


def intersect_interval_lists(left_rows, right_rows):
    i = j = 0
    result = []
    while i < len(left_rows) and j < len(right_rows):
        left_end, right_end = left_rows[i][1], right_rows[j][1]
        left = max(left_rows[i][0], right_rows[j][0])
        right = min(left_end, right_end)
        if left < right:
            result.append((left, right))
        if left_end <= right_end:
            i += 1
        if right_end <= left_end:
            j += 1
    return tuple(result)


def ordered_pair_pieces(a, b):
    rows = []
    for s in range(3):
        for t in range(3):
            if s == t:
                continue
            pieces = intersect_interval_lists(sheet_intervals(a, s), sheet_intervals(b, t))
            check(len(pieces) == 4, f"four pair pieces {a}, {b}, {s}, {t}")
            check(sum((right - left for left, right in pieces), Fraction(0)) == Fraction(20, 1001),
                  f"ordered pair measure {a}, {b}, {s}, {t}")
            rows.extend((s, t, left, right) for left, right in pieces)
    check(len(rows) == 24, f"complete pair piece count {len(rows)}")
    measure = sum((right - left for _, _, left, right in rows), Fraction(0))
    check(measure == Fraction(120, 1001), f"complete pair measure {measure}")
    return tuple(rows), measure


def integrate_cos(q, phase, left, right):
    if q == 0:
        return iv_scale(cos_pi(Fraction(2 * phase, 3)), right - left)
    at_right = sin_pi(2 * q * right + Fraction(2 * phase, 3))
    at_left = sin_pi(2 * q * left + Fraction(2 * phase, 3))
    numerator = iv_sub(at_right, at_left)
    denominator = iv_scale(PI_INTERVAL, 2 * q)
    return iv_div(numerator, denominator)


def exact_pair_remainder(a, b, H, K):
    pieces, pair_measure = ordered_pair_pieces(a, b)
    terms = []
    for s, t, left, right in pieces:
        for n in range(-H, H + 1):
            for m in range(-K, K + 1):
                coefficient = iv_scale(
                    iv_mul(hhat(n), hhat(m)), tau(n, H) * tau(m, K)
                )
                q = n * a + m * b
                phase = n * a * s + m * b * t
                terms.append(iv_mul(coefficient, integrate_cos(q, phase, left, right)))
    pair_product = iv_sum(terms)
    check(pair_product.lo > 0 and pair_product.hi < pair_measure,
          "pair product lies strictly inside pair domain")
    remainder = Interval(pair_measure - pair_product.hi, pair_measure - pair_product.lo)
    return remainder, pair_product, pieces, pair_measure


def epsilon(degree):
    summands = []
    for n in range(1, degree + 1):
        summands.append(iv_scale(iv_square(hhat(n)), tau(n, degree)))
    bracket = iv_sub(Interval(Fraction(1, 7) - Fraction(1, 49)), iv_scale(iv_sum(summands), 2))
    result = iv_scale(bracket, 2)
    check(result.lo > 0, f"positive Fejer L1 error H={degree}")
    return result


def literal_raw_components(w):
    interval_lists = []
    for speed in w:
        inverse = pow(speed, -1, 3)
        rows = []
        for nearest in range(speed + 1):
            left = max(Fraction(0), (Fraction(nearest) - RAW_RADIUS) / speed)
            right = min(Fraction(1), (Fraction(nearest) + RAW_RADIUS) / speed)
            if left < right:
                rows.append((left, right, nearest, (-inverse * nearest) % 3))
        interval_lists.append(rows)

    indices = [0, 0, 0]
    result = {}
    while all(indices[i] < len(interval_lists[i]) for i in range(3)):
        current = tuple(interval_lists[i][indices[i]] for i in range(3))
        left = max(row[0] for row in current)
        right = min(row[1] for row in current)
        if left < right and len({row[3] for row in current}) == 3:
            n = tuple(row[2] for row in current)
            C = cross(w, n)
            check(dot(w, C) == 0 and all(x % 3 for x in C), f"raw physical carrier {w}, {C}")
            result[C] = result.get(C, Fraction(0)) + right - left
        endpoint = min(row[1] for row in current)
        for i in range(3):
            if current[i][1] == endpoint:
                indices[i] += 1
    return result


def open_subset_of_sheet(interval, speed, sheet):
    left, right = interval
    return any(sheet_left <= left < right <= sheet_right
               for sheet_left, sheet_right in sheet_intervals(speed, sheet))


def open_disjoint_from_sheet(interval, speed, sheet):
    left, right = interval
    return all(max(left, sheet_left) >= min(right, sheet_right)
               for sheet_left, sheet_right in sheet_intervals(speed, sheet))


HOSTILE_INTERVALS = (
    ((0, 1), 2, 0, (Fraction(67, 462), Fraction(73, 462))),
    ((0, 2), 1, 0, (Fraction(1, 14), Fraction(17, 210))),
    ((1, 2), 0, 1, (Fraction(0), Fraction(11, 210))),
)


def audit_equality_obstruction():
    w = (1, 5, 11)
    sheets = (0, 1, 2)
    raw = literal_raw_components(w)
    expected = {
        (-1, -2, 1): Fraction(3, 77),
        (1, 2, -1): Fraction(3, 77),
    }
    check(raw == expected, f"equality raw carrier dictionary {raw}")
    check(sum(raw.values(), Fraction(0)) == Fraction(6, 77), "equality comb mass")
    for pair, exact_index, outside_index, interval in HOSTILE_INTERVALS:
        check(pair == tuple(i for i in range(3) if i != exact_index),
              f"hostile pair/exact partition {pair}")
        check(interval[0] < interval[1], f"positive hostile interval {pair}")
        check(open_subset_of_sheet(interval, w[exact_index], sheets[exact_index]),
              f"hostile interval inside exact third sheet {pair}")
        check(outside_index in pair, f"outside coordinate belongs to smoothed pair {pair}")
        check(open_disjoint_from_sheet(interval, w[outside_index], sheets[outside_index]),
              f"hostile interval outside one exact pair sheet {pair}")
    return raw


def decimal_endpoint(value, digits, upper):
    scale = 10 ** digits
    scaled_numerator = value.numerator * scale
    integer = (-((-scaled_numerator) // value.denominator)
               if upper else scaled_numerator // value.denominator)
    sign = "-" if integer < 0 else ""
    integer = abs(integer)
    whole, fractional = divmod(integer, scale)
    return f"{sign}{whole}.{fractional:0{digits}d}"


def decimal_interval(interval, digits=18):
    return f"[{decimal_endpoint(interval.lo, digits, False)},{decimal_endpoint(interval.hi, digits, True)}]"


def fraction_text(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def interval_json(interval):
    return (fraction_text(interval.lo), fraction_text(interval.hi))


def main():
    # Make the frozen byte stream portable across Windows and POSIX runners.
    sys.stdout.reconfigure(newline="\n")
    w = (11, 13, 17)
    pair = (0, 1)
    degrees = (5, 9)
    M, spectrum, weights, killed = finite_M(w, pair, degrees)
    remainder, pair_product, pieces, pair_measure = exact_pair_remainder(11, 13, 5, 9)
    total = iv_add(M, remainder)
    target = Fraction(6, 77)
    slack = Interval(target - total.hi, target - total.lo)
    check(slack.lo > 0, f"certified strict local bound {slack}")

    epsilon5, epsilon9 = epsilon(5), epsilon(9)
    sheet_blind = iv_add(M, iv_scale(iv_add(epsilon5, epsilon9), 3))
    check(sheet_blind.lo > target, f"sheet-blind fallback misses target {sheet_blind}")

    raw_hostile = literal_raw_components(w)
    expected_hostile = {
        (-5, -1, 4): Fraction(10, 1547),
        (5, 1, -4): Fraction(10, 1547),
    }
    check(raw_hostile == expected_hostile, f"hostile raw carrier dictionary {raw_hostile}")
    check(sum(raw_hostile.values(), Fraction(0)) == Fraction(20, 1547), "hostile exact mass")
    raw_equality = audit_equality_obstruction()

    semantic = {
        "pi": interval_json(PI_INTERVAL),
        "spectrum": spectrum,
        "weights": weights,
        "killed": killed,
        "M": interval_json(M),
        "pair_product": interval_json(pair_product),
        "remainder": interval_json(remainder),
        "total": interval_json(total),
        "slack": interval_json(slack),
        "epsilon5": interval_json(epsilon5),
        "epsilon9": interval_json(epsilon9),
        "sheet_blind": interval_json(sheet_blind),
        "pair_measure": fraction_text(pair_measure),
        "pieces": tuple((s, t, fraction_text(left), fraction_text(right))
                        for s, t, left, right in pieces),
        "raw_hostile": tuple(sorted((C, fraction_text(v)) for C, v in raw_hostile.items())),
        "raw_equality": tuple(sorted((C, fraction_text(v)) for C, v in raw_equality.items())),
        "hostile_intervals": HOSTILE_INTERVALS,
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, default=str).encode()).hexdigest()

    print("LRC14 HYBRID FOURIER/EXACT-PAIR CLEAN-ROOM AUDIT")
    print("status=PASS local_certificate_only; LRC14=OPEN")
    print("identity=two_finite_Fejer_factors_plus_one_exact_indicator")
    print("pointwise=f_k*X<=X_+; quotient_sharp=sup_[0,1](zX)=X_+")
    print("spectrum_sites=11 character_counts=-3:8,6:3 sinc_killed=2 live=9")
    print(f"pair_pieces={len(pieces)} pair_measure={fraction_text(pair_measure)}")
    print("M=" + decimal_interval(M))
    print("pair_product=" + decimal_interval(pair_product))
    print("R=" + decimal_interval(remainder))
    print("M_plus_R=" + decimal_interval(total))
    print("six_over_77_minus_bound=" + decimal_interval(slack))
    print("epsilon_5=" + decimal_interval(epsilon5))
    print("epsilon_9=" + decimal_interval(epsilon9))
    print("sheet_blind_bound=" + decimal_interval(sheet_blind))
    print("raw_controls=(11,13,17):20/1547;(1,5,11):6/77")
    print("equality_obstruction=all_3_pairs_all_finite_degrees_strict_slack")
    print("equality_nonclaim=does_not_obstruct_other_hybrid_or_primal_certificates")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
