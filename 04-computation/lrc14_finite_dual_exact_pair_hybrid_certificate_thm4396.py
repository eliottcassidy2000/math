#!/usr/bin/env python3
"""Dependency-free rigorous interval verifier for the THM-4396 hybrid bound.

The proof object has two pieces.  Two of the three danger indicators are
Fejer-smoothed, so their Fourier expansion is finite.  The third indicator is
left exact.  The lost three-way incidence is paid for by an exactly integrated
two-way primal remainder.  All transcendental comparisons use outward-rounded
fixed-point intervals whose sine enclosure is certified by a Taylor remainder.

No repository mathematical implementation is imported.  Every check remains
live under ``python -O``.
"""

from functools import lru_cache
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
from json import dumps
from math import factorial, gcd
import sys


Q = Fraction
SCALE_DIGITS = 36
SCALE = 10 ** SCALE_DIGITS
LAMBDA = Q(1, 14)
RAW_RADIUS = Q(3, 14)
TARGET = Q(6, 77)
W = (11, 13, 17)
SMOOTHED = (0, 1)
EXACT_INDEX = 2
DEGREES = (5, 9)
CHECKS = 0


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise VerificationError(message)


def guard(condition, message):
    """Optimization-safe internal invariant check without inflating CHECKS."""
    if not condition:
        raise VerificationError(message)


def floor_fraction(value):
    return value.numerator // value.denominator


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def floor_div(a, b):
    guard(b != 0, ("nonzero-floor-divisor", b))
    if b < 0:
        a, b = -a, -b
    return a // b


def ceil_div(a, b):
    guard(b != 0, ("nonzero-ceil-divisor", b))
    if b < 0:
        a, b = -a, -b
    return -((-a) // b)


class RI:
    """Closed real interval [lo/SCALE, hi/SCALE], with integer endpoints."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = int(lo)
        self.hi = int(lo if hi is None else hi)
        guard(self.lo <= self.hi, ("ordered-interval", self.lo, self.hi))

    @staticmethod
    def rational(value):
        value = Q(value)
        return RI(
            floor_fraction(value * SCALE),
            ceil_fraction(value * SCALE),
        )

    def __add__(self, other):
        other = as_ri(other)
        return RI(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return RI(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-as_ri(other))

    def __rsub__(self, other):
        return as_ri(other) - self

    def __mul__(self, other):
        other = as_ri(other)
        products = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return RI(
            floor_div(min(products), SCALE),
            ceil_div(max(products), SCALE),
        )

    __rmul__ = __mul__

    def reciprocal(self):
        guard(not (self.lo <= 0 <= self.hi), ("nonzero-denominator", self.lo, self.hi))
        if self.lo > 0:
            return RI(
                floor_div(SCALE * SCALE, self.hi),
                ceil_div(SCALE * SCALE, self.lo),
            )
        return RI(
            floor_div(SCALE * SCALE, self.hi),
            ceil_div(SCALE * SCALE, self.lo),
        )

    def __truediv__(self, other):
        return self * as_ri(other).reciprocal()

    def __rtruediv__(self, other):
        return as_ri(other) / self

    def contains_fraction(self, value):
        value = Q(value)
        return self.lo * value.denominator <= value.numerator * SCALE <= self.hi * value.denominator

    def strictly_below_fraction(self, value):
        value = Q(value)
        return self.hi * value.denominator < value.numerator * SCALE

    def strictly_above_fraction(self, value):
        value = Q(value)
        return self.lo * value.denominator > value.numerator * SCALE


def as_ri(value):
    return value if isinstance(value, RI) else RI.rational(value)


def atan_reciprocal_bounds(denominator, even_index=28):
    """Alternating-series enclosure for atan(1/denominator)."""
    require(denominator >= 2, ("atan-denominator", denominator))
    require(even_index % 2 == 0, ("atan-even-index", even_index))
    partial = Q(0)
    even_partial = None
    odd_partial = None
    for k in range(even_index + 2):
        term = Q(1, (2 * k + 1) * denominator ** (2 * k + 1))
        partial += term if k % 2 == 0 else -term
        if k == even_index:
            even_partial = partial
        if k == even_index + 1:
            odd_partial = partial
    require(odd_partial < even_partial, ("atan-alternation", denominator))
    return odd_partial, even_partial


ATAN5_LO, ATAN5_HI = atan_reciprocal_bounds(5)
ATAN239_LO, ATAN239_HI = atan_reciprocal_bounds(239)
PI_LO_Q = 16 * ATAN5_LO - 4 * ATAN239_HI
PI_HI_Q = 16 * ATAN5_HI - 4 * ATAN239_LO
PI = RI.rational(PI_LO_Q)
PI = RI(PI.lo, RI.rational(PI_HI_Q).hi)


def ri_power(base, exponent):
    guard(exponent >= 0, ("nonnegative-exponent", exponent))
    result = RI.rational(1)
    factor = base
    power = exponent
    while power:
        if power & 1:
            result = result * factor
        factor = factor * factor
        power //= 2
    return result


SIN_TERMS_LAST_INDEX = 18
SIN_REMAINDER = Q(11, 7) ** (2 * SIN_TERMS_LAST_INDEX + 3) / factorial(
    2 * SIN_TERMS_LAST_INDEX + 3
)


def sin_point_scaled(x_scaled):
    """Rigorous sine interval at the exactly represented x_scaled/SCALE."""
    guard(0 <= x_scaled and 7 * x_scaled <= 11 * SCALE, ("sine-reduced-domain", x_scaled))
    x = RI(x_scaled)
    x2 = x * x
    term = x
    total = x
    for k in range(1, SIN_TERMS_LAST_INDEX + 1):
        term = -term * x2 / Q((2 * k) * (2 * k + 1))
        total = total + term
    remainder = RI.rational(SIN_REMAINDER)
    return total + RI(-remainder.hi, remainder.hi)


@lru_cache(maxsize=None)
def sin_pi(value):
    """Rigorous interval for sin(pi*value), value a Fraction."""
    value = Q(value)
    reduced = value % 2
    sign = 1
    if reduced > 1:
        reduced -= 1
        sign = -1
    if reduced > Q(1, 2):
        reduced = 1 - reduced
    if reduced == 0:
        return RI(0)
    x = PI * reduced
    lower = sin_point_scaled(x.lo)
    upper = sin_point_scaled(x.hi)
    result = RI(lower.lo, upper.hi)
    return result if sign > 0 else -result


@lru_cache(maxsize=None)
def cos_pi(value):
    return sin_pi(Q(value) + Q(1, 2))


@lru_cache(maxsize=None)
def hhat(n):
    if n == 0:
        return RI.rational(Q(1, 7))
    n = abs(n)
    return sin_pi(Q(n, 7)) / (PI * n)


def fejer_coefficient(n, degree):
    if abs(n) > degree:
        return RI(0)
    return hhat(n) * Q(degree + 1 - abs(n), degree + 1)


@lru_cache(maxsize=None)
def integral_cos_pi(frequency, phase, left, right):
    """Integral of cos(2*pi*frequency*x + pi*phase) on [left,right]."""
    frequency = int(frequency)
    phase = Q(phase)
    left = Q(left)
    right = Q(right)
    require(left <= right, ("integral-order", left, right))
    if frequency == 0:
        return (right - left) * cos_pi(phase)
    numerator = sin_pi(2 * frequency * right + phase) - sin_pi(
        2 * frequency * left + phase
    )
    return numerator / (2 * frequency * PI)


def sheet_intervals(speed, sheet):
    radius = LAMBDA / speed
    result = []
    for owner in range(speed):
        center = Q(owner, speed) - Q(sheet, 3)
        center -= floor_fraction(center)
        left, right = center - radius, center + radius
        if left < 0:
            result.extend(((Q(0), right), (left + 1, Q(1))))
        elif right > 1:
            result.extend(((left, Q(1)), (Q(0), right - 1)))
        else:
            result.append((left, right))
    result.sort()
    return tuple(result)


def intersect_interval_lists(first, second):
    i = j = 0
    result = []
    while i < len(first) and j < len(second):
        first_right = first[i][1]
        second_right = second[j][1]
        left = max(first[i][0], second[j][0])
        right = min(first_right, second_right)
        if left < right:
            result.append((left, right))
        if first_right <= second_right:
            i += 1
        if second_right <= first_right:
            j += 1
    return tuple(result)


def ordered_distinct_pair_pieces(speed_a, speed_b):
    cache_a = tuple(sheet_intervals(speed_a, sheet) for sheet in range(3))
    cache_b = tuple(sheet_intervals(speed_b, sheet) for sheet in range(3))
    return tuple(
        (sheet_a, sheet_b, left, right)
        for sheet_a in range(3)
        for sheet_b in range(3)
        if sheet_a != sheet_b
        for left, right in intersect_interval_lists(
            cache_a[sheet_a], cache_b[sheet_b]
        )
    )


def pair_product_integral(speed_a, speed_b, degree_a, degree_b, pieces):
    """Integral over the ordered distinct-sheet pair domain of g_a*g_b."""
    total = RI(0)
    # The (-n,-m) term is the complex conjugate and has the same real
    # integral.  Retaining one representative halves the interval workload.
    for n in range(-degree_a, degree_a + 1):
        coefficient_a = fejer_coefficient(n, degree_a)
        for m in range(-degree_b, degree_b + 1):
            if n < 0 or (n == 0 and m < 0):
                continue
            coefficient_b = fejer_coefficient(m, degree_b)
            multiplicity = 1 if n == 0 and m == 0 else 2
            mode_integral = RI(0)
            for sheet_a, sheet_b, left, right in pieces:
                frequency = n * speed_a + m * speed_b
                phase = Q(2 * (n * speed_a * sheet_a + m * speed_b * sheet_b), 3)
                mode_integral += integral_cos_pi(
                    frequency, phase, left, right
                )
            total += multiplicity * coefficient_a * coefficient_b * mode_integral
    return total


def dot(first, second):
    return sum(x * y for x, y in zip(first, second))


def sheet_weight(w, n):
    residues = tuple((w[index] * n[index]) % 3 for index in range(3))
    require(sum(residues) % 3 == 0, ("resonance-residues", w, n, residues))
    equal = len(set(residues)) == 1
    distinct = len(set(residues)) == 3
    require(equal or distinct, ("ternary-dichotomy", w, n, residues))
    return 6 if equal else -3


def finite_hybrid_sum(w, smoothed, exact_index, degrees):
    first, second = smoothed
    degree_first, degree_second = degrees
    result = RI(0)
    records = []
    for n_first in range(-degree_first, degree_first + 1):
        for n_second in range(-degree_second, degree_second + 1):
            numerator = -(w[first] * n_first + w[second] * n_second)
            if numerator % w[exact_index]:
                continue
            n_exact = numerator // w[exact_index]
            n = [0, 0, 0]
            n[first] = n_first
            n[second] = n_second
            n[exact_index] = n_exact
            n = tuple(n)
            require(dot(w, n) == 0, ("exact-resonance", w, n))
            require(abs(n[first]) <= degree_first and abs(n[second]) <= degree_second,
                    ("frequency-window", n, degrees))
            weight = sheet_weight(w, n)
            multiplier_first = Q(degree_first + 1 - abs(n_first), degree_first + 1)
            multiplier_second = Q(degree_second + 1 - abs(n_second), degree_second + 1)
            term = (
                weight
                * multiplier_first
                * multiplier_second
                * hhat(n[0])
                * hhat(n[1])
                * hhat(n[2])
            )
            killed = any(coordinate and coordinate % 7 == 0 for coordinate in n)
            if killed:
                require(term.lo == term.hi == 0, ("sinc-zero", n, term.lo, term.hi))
            result += term
            records.append(
                (
                    n,
                    weight,
                    multiplier_first,
                    multiplier_second,
                    killed,
                )
            )
    return result, tuple(records)


def epsilon_l1(degree):
    """Exact interval for ||h-Fejer_degree*h||_1, h=1_{||x||<1/14}."""
    inside = RI.rational(Q(1, 49))
    for n in range(1, degree + 1):
        inside += 2 * Q(degree + 1 - n, degree + 1) * hhat(n) * hhat(n)
    return 2 * (Q(1, 7) - inside)


def strict_carrier_bound(value):
    return (value.numerator - 1) // value.denominator


def carrier_length(w, carrier):
    w1, w2, w3 = w
    c1, c2, c3 = carrier
    return max(
        Q(0),
        min(
            2 * RAW_RADIUS / w1,
            2 * RAW_RADIUS / w2,
            2 * RAW_RADIUS / w3,
            RAW_RADIUS / w1 + RAW_RADIUS / w2 - Q(abs(c3), w1 * w2),
            RAW_RADIUS / w1 + RAW_RADIUS / w3 - Q(abs(c2), w1 * w3),
            RAW_RADIUS / w2 + RAW_RADIUS / w3 - Q(abs(c1), w2 * w3),
        ),
    )


def raw_carriers(w):
    bounds = (
        strict_carrier_bound(RAW_RADIUS * (w[1] + w[2])),
        strict_carrier_bound(RAW_RADIUS * (w[0] + w[2])),
        strict_carrier_bound(RAW_RADIUS * (w[0] + w[1])),
    )
    result = {}
    for c1 in range(-bounds[0], bounds[0] + 1):
        for c2 in range(-bounds[1], bounds[1] + 1):
            numerator = -(w[0] * c1 + w[1] * c2)
            if numerator % w[2]:
                continue
            c3 = numerator // w[2]
            if abs(c3) > bounds[2]:
                continue
            carrier = (c1, c2, c3)
            if any(coordinate % 3 == 0 for coordinate in carrier):
                continue
            length = carrier_length(w, carrier)
            if length:
                result[carrier] = length
    return result


def subtract_interval_lists(first, second):
    result = []
    for left, right in first:
        cursor = left
        for cut_left, cut_right in second:
            if cut_right <= cursor or cut_left >= right:
                continue
            if cut_left > cursor:
                result.append((cursor, min(cut_left, right)))
            cursor = max(cursor, cut_right)
            if cursor >= right:
                break
        if cursor < right:
            result.append((cursor, right))
    return tuple((left, right) for left, right in result if left < right)


def equality_hostile_witnesses():
    """Positive-measure exact-sheet regions omitted by every pair quotient."""
    w = (1, 5, 11)
    assignment = (0, 1, 2)
    result = []
    for first, second in ((0, 1), (0, 2), (1, 2)):
        exact_index = 3 - first - second
        difference = subtract_interval_lists(
            sheet_intervals(w[exact_index], assignment[exact_index]),
            sheet_intervals(w[first], assignment[first]),
        )
        require(difference, ("hostile-positive-difference", first, second))
        left, right = difference[0]
        require(left < right, ("hostile-open-witness", first, second, left, right))
        result.append(((first, second), exact_index, assignment, left, right))
    return tuple(result)


def decimal_outward(interval, digits=18):
    require(0 <= digits <= SCALE_DIGITS, ("decimal-digits", digits))
    divisor = 10 ** (SCALE_DIGITS - digits)
    low = floor_div(interval.lo, divisor)
    high = ceil_div(interval.hi, divisor)

    def render(value):
        sign = "-" if value < 0 else ""
        value = abs(value)
        if digits == 0:
            return sign + str(value)
        whole, fractional = divmod(value, 10**digits)
        return f"{sign}{whole}.{fractional:0{digits}d}"

    return "[" + render(low) + "," + render(high) + "]"


def fraction_text(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    # Freeze LF output on every supported host so replay is genuinely bytewise.
    sys.stdout.reconfigure(newline="\n")
    require(gcd(gcd(*W[:2]), W[2]) == 1, "primitive fixed comb")
    require(len(set(W)) == 3 and all(value > 0 and value % 2 and value % 3 for value in W),
            "distinct positive odd ternary units")
    require(Q(103993, 33102) < PI_LO_Q < PI_HI_Q < Q(104348, 33215), "pi enclosure")
    require(sin_pi(0).contains_fraction(0), "sin zero")
    require(sin_pi(Q(1, 2)).contains_fraction(1), "sin half")
    require(cos_pi(0).contains_fraction(1), "cos zero")

    pieces = ordered_distinct_pair_pieces(W[0], W[1])
    pair_measure = sum((right - left for _, _, left, right in pieces), Q(0))
    require(len(pieces) == 24, ("pair-piece-count", len(pieces)))
    require(pair_measure == Q(120, 1001), ("pair-measure", pair_measure))

    pair_product = pair_product_integral(W[0], W[1], DEGREES[0], DEGREES[1], pieces)
    pair_remainder = pair_measure - pair_product
    require(pair_remainder.lo > 0, ("positive-pair-remainder", pair_remainder.lo))

    low_sum, records = finite_hybrid_sum(W, SMOOTHED, EXACT_INDEX, DEGREES)
    expected_vectors = (
        (-5, -1, 4),
        (-4, 6, -2),
        (-3, -4, 5),
        (-2, 3, -1),
        (-1, -7, 6),
        (0, 0, 0),
        (1, 7, -6),
        (2, -3, 1),
        (3, 4, -5),
        (4, -6, 2),
        (5, 1, -4),
    )
    require(tuple(record[0] for record in records) == expected_vectors, "finite spectrum")
    require(len(records) == 11, ("finite-spectrum-count", len(records)))
    require(sum(record[4] for record in records) == 2, "two exact sinc zeros")
    require(sum(not record[4] for record in records) == 9, "nine nonzero terms")
    require(sum(record[1] == 6 for record in records) == 3, "character plus count")
    require(sum(record[1] == -3 for record in records) == 8, "character minus count")

    hybrid_upper = low_sum + pair_remainder
    require(hybrid_upper.strictly_below_fraction(TARGET), ("hybrid-closes", hybrid_upper.hi))
    target_gap = TARGET - hybrid_upper
    require(target_gap.lo > 0, ("positive-target-gap", target_gap.lo))

    epsilon_first = epsilon_l1(DEGREES[0])
    epsilon_second = epsilon_l1(DEGREES[1])
    coarse_upper = low_sum + 3 * (epsilon_first + epsilon_second)
    require(coarse_upper.strictly_above_fraction(TARGET), ("coarse-does-not-close", coarse_upper.lo))

    exact_raw = raw_carriers(W)
    expected_raw = {
        (-5, -1, 4): Q(10, 1547),
        (5, 1, -4): Q(10, 1547),
    }
    require(exact_raw == expected_raw, ("physical-positive-control", exact_raw))
    exact_measure = sum(exact_raw.values(), Q(0))
    require(exact_measure == Q(20, 1547), ("physical-measure", exact_measure))
    require(hybrid_upper.strictly_above_fraction(exact_measure), "certificate above exact value")

    equality_raw = raw_carriers((1, 5, 11))
    require(dot((1, 2, -1), (1, 5, 11)) == 0, "minimal-hostile relation")
    require(sum(abs(value) for value in (1, 2, -1)) == 4, "minimal-hostile l1 norm")
    require(
        equality_raw
        == {
            (-1, -2, 1): Q(3, 77),
            (1, 2, -1): Q(3, 77),
        },
        ("sharp-equality-control", equality_raw),
    )
    require(sum(equality_raw.values(), Q(0)) == TARGET, "sharp equality measure")
    hostile_witnesses = equality_hostile_witnesses()
    expected_hostile_intervals = (
        ((0, 1), 2, (0, 1, 2), Q(67, 462), Q(73, 462)),
        ((0, 2), 1, (0, 1, 2), Q(1, 14), Q(17, 210)),
        ((1, 2), 0, (0, 1, 2), Q(0), Q(11, 210)),
    )
    require(hostile_witnesses == expected_hostile_intervals, "equality hostile intervals")

    semantic = {
        "fixed_comb": W,
        "smoothed_indices": SMOOTHED,
        "exact_index": EXACT_INDEX,
        "degrees": DEGREES,
        "spectrum": [
            (
                record[0],
                record[1],
                fraction_text(record[2]),
                fraction_text(record[3]),
                record[4],
            )
            for record in records
        ],
        "pair_pieces": len(pieces),
        "pair_measure": fraction_text(pair_measure),
        "low_interval_scaled": (low_sum.lo, low_sum.hi, SCALE),
        "pair_product_interval_scaled": (pair_product.lo, pair_product.hi, SCALE),
        "pair_remainder_interval_scaled": (pair_remainder.lo, pair_remainder.hi, SCALE),
        "upper_interval_scaled": (hybrid_upper.lo, hybrid_upper.hi, SCALE),
        "gap_interval_scaled": (target_gap.lo, target_gap.hi, SCALE),
        "coarse_interval_scaled": (coarse_upper.lo, coarse_upper.hi, SCALE),
        "exact_physical": fraction_text(exact_measure),
        "equality_hostile": [
            (pair, exact_index, assignment, fraction_text(left), fraction_text(right))
            for pair, exact_index, assignment, left, right in hostile_witnesses
        ],
        "scope": "fixed local three-speed comb; LRC14 OPEN",
    }
    digest = sha256(dumps(semantic, sort_keys=True).encode()).hexdigest()

    print("LRC THM4396 TWO-FACTOR FEJER / EXACT-PAIR HYBRID CERTIFICATE")
    print("status=PASS RIGOROUS_INTERVAL; scope=fixed_local_three_speed_comb; LRC14=OPEN")
    print("comb=(11,13,17); smoothed_speeds=(11,13); exact_speed=17; degrees=(5,9)")
    print("finite_resonance_sites=11; nonzero_sinc_terms=9; exact_sinc_zeros=2")
    print("character_counts=weight6:3,weight-3:8")
    print("finite_spectrum=" + str(expected_vectors))
    print("low_sum_interval=" + decimal_outward(low_sum))
    print("ordered_distinct_pair_pieces=24; pair_domain_measure=120/1001")
    print("pair_product_interval=" + decimal_outward(pair_product))
    print("exact_pair_remainder_interval=" + decimal_outward(pair_remainder))
    print("hybrid_upper_interval=" + decimal_outward(hybrid_upper))
    print("target=6/77; certified_gap_interval=" + decimal_outward(target_gap))
    print("coarse_L1_upper_interval=" + decimal_outward(coarse_upper) + "; closes_6/77=no")
    print("physical_control_exact=20/1547; used_by_certificate=no")
    print("sharp_hostile=(1,5,11); exact=6/77; strict_for=all_3_coordinate_pairs_and_all_finite_degree_pairs")
    print("hostile_open_intervals=" + str(expected_hostile_intervals))
    print("preserved=dual_character,sinc_signs_and_zeros,resonance_equation,two_frequency_windows,ordered_pair_sheet_geometry")
    print("lost=third_sheet_incidence_in_remainder,negative_tail_cancellation,carrier_owner_and_entry_data")
    print("optimization_safe_checks=yes")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
