#!/usr/bin/env python3
"""Exact referee for THM-2257's 169-image fibre sieve."""

from fractions import Fraction
from math import gcd


M = 169
TARGET = Fraction(457, 1183)
UNIFORM_LOWER = Fraction(164775, 426496)


def comb_intervals(speed, denominator):
    """Components in [0,1] of {x: ||speed*x|| < 1/denominator}."""
    radius = Fraction(1, denominator * speed)
    intervals = [(Fraction(), radius)]
    for index in range(1, speed):
        center = Fraction(index, speed)
        intervals.append((center - radius, center + radius))
    intervals.append((1 - radius, Fraction(1)))
    return intervals


def intersect_intervals(left, right):
    output = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            output.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return output


def merge(intervals):
    output = []
    for lo, hi in sorted(intervals):
        if output and lo <= output[-1][1]:
            output[-1] = (output[-1][0], max(output[-1][1], hi))
        else:
            output.append((lo, hi))
    return output


def multiply_image(intervals, multiplier=M):
    """Exact set image of an interval union under x -> multiplier*x mod 1."""
    pieces = []
    for lo, hi in intervals:
        length = multiplier * (hi - lo)
        if length >= 1:
            return [(Fraction(), Fraction(1))]
        raw_lo = multiplier * lo
        raw_hi = multiplier * hi
        floor_lo = raw_lo.numerator // raw_lo.denominator
        floor_hi = raw_hi.numerator // raw_hi.denominator
        frac_lo = raw_lo - floor_lo
        frac_hi = raw_hi - floor_hi
        if floor_lo == floor_hi:
            pieces.append((frac_lo, frac_hi))
        else:
            pieces.append((frac_lo, Fraction(1)))
            if frac_hi:
                pieces.append((Fraction(), frac_hi))
    return merge(pieces)


def exact_image_mass(a, b):
    """Measure of 169(E_a intersection D_b), for a coprime reduced pair."""
    assert gcd(a, b) == 1
    assert gcd(a * b, M) == 1
    source = intersect_intervals(
        comb_intervals(a, 7),
        comb_intervals(b, 14),
    )
    image = multiply_image(source)
    return sum((hi - lo for lo, hi in image), Fraction())


def fibre_multiplicity(r):
    """Maximum 25-block / r*(49-block) intersection in Z/169."""
    danger = {value % M for value in range(-12, 13)}
    eligibility = [value % M for value in range(-24, 25)]
    return max(
        len(
            danger
            & {(shift + r * value) % M for value in eligibility}
        )
        for shift in range(M)
    )


def overlap_bound(k_value, product):
    """THM-2080's F >= -1/8 bound, followed by the fibre cap."""
    return Fraction(M, k_value) * (
        Fraction(2, 49) - Fraction(1, 4 * product)
    )


units = [r for r in range(1, M) if gcd(r, M) == 1]
k_by_r = {r: fibre_multiplicity(r) for r in units}
classes = {}
for r, k_value in k_by_r.items():
    classes.setdefault(k_value, []).append(r)

expected_counts = {8: 42, 9: 72, 10: 26, 12: 2, 13: 8, 17: 2, 25: 4}
assert {key: len(value) for key, value in classes.items()} == expected_counts
assert classes[12] == [75, 94]
assert classes[13] == [2, 4, 42, 57, 112, 127, 165, 167]
assert classes[17] == [56, 113]
assert classes[25] == [1, 84, 85, 168]

assert Fraction(M, 42 * 10) - TARGET == Fraction(163, 10140)

thresholds = {12: 19, 13: 23, 17: 128}
threshold_gaps = {
    12: Fraction(24295, 7552272),
    13: Fraction(2287, 761852),
    17: Fraction(2879, 72077824),
}
for k_value, threshold in thresholds.items():
    assert overlap_bound(k_value, threshold) > TARGET
    assert overlap_bound(k_value, threshold - 1) <= TARGET
    assert overlap_bound(k_value, threshold) - TARGET == threshold_gaps[k_value]
assert overlap_bound(17, 128) == UNIFORM_LOWER
assert UNIFORM_LOWER - TARGET == Fraction(2879, 72077824)


def small_candidates(k_value, product_cap):
    output = []
    for a in range(1, product_cap + 1, 2):
        if gcd(a, M) != 1:
            continue
        for b in range(1, product_cap // a + 1):
            if gcd(b, M) != 1 or gcd(a, b) != 1:
                continue
            r = b * pow(a, -1, M) % M
            if k_by_r[r] == k_value:
                output.append((a, b, r, a * b))
    return output


finite_cases = {
    12: small_candidates(12, 18),
    13: small_candidates(13, 22),
    17: small_candidates(17, 127),
}
assert finite_cases == {
    12: [(9, 1, 94, 9)],
    13: [(1, 2, 2, 2), (1, 4, 4, 4), (3, 2, 57, 6)],
    17: [(1, 56, 56, 56), (1, 113, 113, 113), (3, 1, 113, 3)],
}

expected_finite_masses = {
    (9, 1): Fraction(1),
    (1, 2): Fraction(1),
    (1, 4): Fraction(1),
    (3, 2): Fraction(1),
    (1, 56): Fraction(267, 392),
    (1, 113): Fraction(555, 791),
    (3, 1): Fraction(1),
}
for pair, expected in expected_finite_masses.items():
    assert exact_image_mass(*pair) == expected
    assert expected > UNIFORM_LOWER

# The exceptional mesh-strip bounds.  The second has the necessary factor
# two: in the variable s=b*x the overlap length is 1/7-||ell*y||/2.
same_sign_strip = 2 * (Fraction(3, 14) - Fraction(1, M))
double_sign_strip = 2 * (Fraction(2, 7) - Fraction(2, M))
assert same_sign_strip == Fraction(493, 1183)
assert double_sign_strip == Fraction(648, 1183)
assert same_sign_strip - TARGET == Fraction(36, 1183)
assert double_sign_strip - TARGET == Fraction(191, 1183)

# Hostile resonant probes: the first is the minimum in the exact a,b<=300
# scout that motivated the proof, but no sharpness claim relies on that scan.
hostile_probes = {
    (159, 10): Fraction(464, 1113),
    (135, 34): Fraction(394, 945),
    (87, 82): Fraction(254, 609),
    (1, 168): Fraction(491, 1176),
}
for pair, expected in hostile_probes.items():
    assert exact_image_mass(*pair) == expected
    assert expected >= same_sign_strip

print("THM-2257 DEPTH-THREE COMMON-CORE 169-IMAGE SIEVE")
print("m", M)
print("target", TARGET)
print("uniform_lower", UNIFORM_LOWER)
print("uniform_gap", UNIFORM_LOWER - TARGET)
print()
print("fibre_K_distribution")
for k_value in sorted(classes):
    print(k_value, len(classes[k_value]))
print("K12_residues", classes[12])
print("K13_residues", classes[13])
print("K17_residues", classes[17])
print("K25_residues", classes[25])
print()
print("large_product_thresholds")
for k_value in (12, 13, 17):
    threshold = thresholds[k_value]
    print(
        k_value,
        threshold,
        overlap_bound(k_value, threshold),
        overlap_bound(k_value, threshold) - TARGET,
    )
print()
print("finite_cases")
for k_value in (12, 13, 17):
    print("K", k_value)
    for a, b, r, product in finite_cases[k_value]:
        print(a, b, r, product, exact_image_mass(a, b))
print()
print("exceptional_strip_bounds")
print("r=+-1", same_sign_strip, same_sign_strip - TARGET)
print("r=84,85", double_sign_strip, double_sign_strip - TARGET)
print()
print("hostile_resonant_probes")
for pair, expected in hostile_probes.items():
    print(pair[0], pair[1], expected, expected - same_sign_strip)
print("status=PROVED_VERIFIED_EXACT")
