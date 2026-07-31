"""Exact companion for THM-2638.

The script checks the radial-shell lift, the private-height decoder, the
equal-height source-character hostile, and the external mixed-radix selector.
All arithmetic is exact.
"""

from fractions import Fraction
from itertools import product
from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def weak_compositions(total, slots):
    if slots == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, slots - 1):
            yield (first,) + tail


def channel_rows(support, moment):
    charges = tuple(a - b for a, b in support)
    rows = []
    for occ in weak_compositions(moment, len(support)):
        if sum(q * r for q, r in zip(charges, occ)) != 0:
            continue
        hol = sum(a * r for (a, _), r in zip(support, occ))
        anti = sum(b * r for (_, b), r in zip(support, occ))
        require(hol == anti, f"unbalanced retained row: {occ}")
        multinomial = factorial(moment)
        for r in occ:
            multinomial //= factorial(r)
        wick = multinomial * factorial(hol)
        rows.append((occ, hol, multinomial, wick))
    return tuple(rows)


def grade_fibres(rows):
    fibres = {}
    for occ, height, _, _ in rows:
        fibres.setdefault(height, []).append(occ)
    return {height: tuple(occs) for height, occs in fibres.items()}


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


hostile = ((6, 0), (0, 2), (0, 18))
private_levels = 0
channel_counts = []
height_ranges = []
for j in range(1, 13):
    moment = 4 * j
    rows = channel_rows(hostile, moment)
    expected_occ = tuple((j + 2 * t, 3 * (j - t), t) for t in range(j + 1))
    expected_heights = tuple(6 * j + 12 * t for t in range(j + 1))
    got_occ = tuple(row[0] for row in rows)
    got_heights = tuple(row[1] for row in rows)
    require(got_occ == expected_occ, f"occupation parametrization failed at j={j}")
    require(got_heights == expected_heights, f"height parametrization failed at j={j}")
    fibres = grade_fibres(rows)
    require(all(len(fibre) == 1 for fibre in fibres.values()),
            f"radial height is not private at j={j}")
    require(all(wick == multinomial * factorial(height)
                for _, height, multinomial, wick in rows),
            f"Wick/Laplace factorization failed at j={j}")
    private_levels += 1
    channel_counts.append(len(rows))
    height_ranges.append((got_heights[0], got_heights[-1]))


rows4 = channel_rows(hostile, 4)
require(rows4 == (
    ((1, 3, 0), 6, 4, 4 * factorial(6)),
    ((3, 0, 1), 18, 4, 4 * factorial(18)),
), "fourth-moment row mismatch")

a = Fraction(1)
b = Fraction(1)
c = -Fraction(factorial(6), factorial(18))
shell4 = {
    height: multinomial * (a ** occ[0]) * (b ** occ[1]) * (c ** occ[2])
    for occ, height, multinomial, _ in rows4
}
require(shell4[6] != 0 and shell4[18] != 0,
        "fourth shell unexpectedly vanished")
require(sum(factorial(height) * coefficient
            for height, coefficient in shell4.items()) == 0,
        "fourth scalar cancellation failed")


equal_height = ((2, 0), (1, 1), (0, 2))
equal_rows = channel_rows(equal_height, 2)
require(equal_rows == (
    ((0, 2, 0), 2, 1, 2),
    ((1, 0, 1), 2, 2, 4),
), "equal-height hostile rows mismatch")
require(grade_fibres(equal_rows) == {2: ((0, 2, 0), (1, 0, 1))},
        "equal-height hostile did not collide")

# Finite controls for the universal identity that every source-torus
# one-parameter grading restricts to (alpha+beta)A.
source_character_checks = 0
for alpha, beta in product(range(-3, 4), repeat=2):
    monomial_grades = tuple(alpha * aa + beta * bb for aa, bb in equal_height)
    channel_grades = tuple(dot(monomial_grades, occ)
                           for occ, _, _, _ in equal_rows)
    expected = tuple((alpha + beta) * height
                     for _, height, _, _ in equal_rows)
    require(channel_grades == expected,
            f"source-torus grade formula failed at {(alpha, beta)}")
    if alpha + beta == 0:
        require(channel_grades == (0, 0),
                f"angular grade failed at {(alpha, beta)}")
    source_character_checks += 1

# Mixed-radix coefficient weights separate every bounded occupation vector.
moment = 2
base = moment + 1
weights = tuple(base ** i for i in range(len(equal_height)))
modulus = base ** len(equal_height)
all_occ = tuple(weak_compositions(moment, len(equal_height)))
codes = tuple(dot(weights, occ) % modulus for occ in all_occ)
require(len(set(codes)) == len(all_occ), "mixed-radix code collision")

# Exact cyclic root orthogonality is an exponent-count statement.
selector_checks = 0
for target in all_occ:
    target_code = dot(weights, target) % modulus
    for occ in all_occ:
        code = dot(weights, occ) % modulus
        exponent_histogram = [0] * modulus
        for u in range(modulus):
            exponent_histogram[(u * (code - target_code)) % modulus] += 1
        if occ == target:
            require(exponent_histogram[0] == modulus,
                    f"target selector lost its trivial phase: {target}")
            require(sum(exponent_histogram[1:]) == 0,
                    f"target selector gained a nontrivial phase: {target}")
        else:
            # Sum_u xi^(u d)=0 because d is nonzero modulo N.  For this
            # fixture the histogram is uniform on a nontrivial subgroup;
            # the geometric-series criterion N does not divide d is exact.
            require((code - target_code) % modulus != 0,
                    f"off-target mixed-radix collision: {target}, {occ}")
        selector_checks += 1


print("THM-2638 exact radial-height graded Wick decoder")
print(f"hostile_private_levels={private_levels}/12 channel_counts={tuple(channel_counts)}")
print(f"hostile_height_ranges={tuple(height_ranges)}")
print(f"fourth_shell_heights={tuple(shell4)} scalar_pushforward=0 shell_nonzero=1")
print(f"equal_height_channels={tuple(row[0] for row in equal_rows)} height=2 radial_private=0")
print(f"source_character_checks={source_character_checks} angular_classes=1")
print(f"mixed_radix_modulus={modulus} codes={codes} selector_checks={selector_checks}")
print("verdict=PASS: radial height restores the private hostile rows; scalar L and source characters retain the stated boundaries")
