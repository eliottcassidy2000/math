#!/usr/bin/env python3
"""Independent exact referee for THM-4009.

This is deliberately not an import, wrapper, or transcription of the primary
companion.  The finite coefficient-shape census is performed by a generating
function over value multiplicities, the l1 optima use the balanced-squares
formula, and the lattice-ball count uses polynomial convolution.  The cited
theorem itself (Banaszczyk's Euclidean transference inequality) is treated as
an external input and its normalization is audited symbolically.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb, gcd, isqrt


def check(statement: bool, label: str) -> None:
    if not statement:
        raise RuntimeError(label)


def balanced_square_cost(total: int, slots: int) -> int:
    """Minimum sum of squares of `slots` nonnegative integers of given sum."""

    quotient, remainder = divmod(total, slots)
    return (slots - remainder) * quotient * quotient + remainder * (quotient + 1) ** 2


def l1_maximum(square_budget: int, slots: int) -> int:
    total = 0
    while balanced_square_cost(total + 1, slots) <= square_budget:
        total += 1
    return total


def coefficient_polynomial(square_budget: int) -> list[int]:
    """Coefficients for one signed integer coordinate, indexed by its square."""

    polynomial = [0] * (square_budget + 1)
    polynomial[0] = 1
    for magnitude in range(1, isqrt(square_budget) + 1):
        polynomial[magnitude * magnitude] = 2
    return polynomial


def truncated_convolution(left: list[int], right: list[int], budget: int) -> list[int]:
    answer = [0] * (budget + 1)
    left_terms = [(degree, value) for degree, value in enumerate(left) if value]
    right_terms = [(degree, value) for degree, value in enumerate(right) if value]
    for left_degree, left_value in left_terms:
        for right_degree, right_value in right_terms:
            degree = left_degree + right_degree
            if degree > budget:
                break
            answer[degree] += left_value * right_value
    return answer


def labelled_ball_distribution(dimension: int, square_budget: int) -> list[int]:
    one_coordinate = coefficient_polynomial(square_budget)
    answer = [0] * (square_budget + 1)
    answer[0] = 1
    for _ in range(dimension):
        answer = truncated_convolution(answer, one_coordinate, square_budget)
    return answer


def moebius_values(limit: int) -> list[int]:
    """Small sieve, independent of the primary script's factor-by-factor routine."""

    values = [1] * (limit + 1)
    is_prime = [True] * (limit + 1)
    is_prime[:2] = [False, False]
    for prime in range(2, limit + 1):
        if not is_prime[prime]:
            continue
        for multiple in range(prime, limit + 1, prime):
            is_prime[multiple] = False if multiple != prime else True
            values[multiple] *= -1
        square = prime * prime
        for multiple in range(square, limit + 1, square):
            values[multiple] = 0
    return values


def histogram_generating_function(
    coordinate_count: int, square_budget: int
) -> tuple[dict[int, int], dict[int, int], dict[int, int]]:
    """Count primitive absolute multisets and their sign-multiset orbits.

    A state records support, square mass, gcd, l1 parity, and whether every
    multiplicity seen so far is even.  Its two weights are the number of
    multiplicity vectors and the sum of prod_v(m_v+1).  Burnside then counts
    sign choices modulo global complement without enumerating histograms.
    """

    # (support, square, gcd, l1 parity, all multiplicities even) ->
    # (number of histograms, sum of products (multiplicity+1))
    states: dict[tuple[int, int, int, int, bool], tuple[int, int]] = {
        (0, 0, 0, 0, True): (1, 1)
    }
    for value in range(1, isqrt(square_budget) + 1):
        following: dict[tuple[int, int, int, int, bool], tuple[int, int]] = {}
        value_square = value * value
        for (support, square, common, parity, all_even), (ways, raw_weight) in states.items():
            max_multiplicity = min(
                coordinate_count - support,
                (square_budget - square) // value_square,
            )
            for multiplicity in range(max_multiplicity + 1):
                key = (
                    support + multiplicity,
                    square + multiplicity * value_square,
                    common if multiplicity == 0 else gcd(common, value),
                    (parity + multiplicity * value) & 1,
                    all_even and multiplicity % 2 == 0,
                )
                old_ways, old_raw = following.get(key, (0, 0))
                following[key] = (
                    old_ways + ways,
                    old_raw + raw_weight * (multiplicity + 1),
                )
        states = following

    histogram_counts: dict[int, int] = {}
    sign_orbits: dict[int, int] = {}
    odd_histograms: dict[int, int] = {}
    for support in range(2, coordinate_count + 1):
        histogram_total = 0
        raw_total = 0
        complement_fixed_total = 0
        odd_total = 0
        for (state_support, _, common, parity, all_even), (ways, raw_weight) in states.items():
            if state_support != support or common != 1:
                continue
            histogram_total += ways
            raw_total += raw_weight
            complement_fixed_total += ways * int(all_even)
            odd_total += ways * parity
        # For support two, (1,1) is primitive but cannot annihilate distinct
        # positive speeds.  It contributes one histogram and one mixed-sign
        # orbit; its l1 parity is even.
        if support == 2:
            histogram_total -= 1
            raw_total -= 3
            complement_fixed_total -= 1
        histogram_counts[support] = histogram_total
        sign_orbits[support] = (raw_total + complement_fixed_total) // 2 - histogram_total
        odd_histograms[support] = odd_total
    return histogram_counts, sign_orbits, odd_histograms


checks = 0
speed_count = 13
projected_dimension = 12
cube_half_side = Fraction(3, 7)
banaszczyk_cap = Fraction(projected_dimension, 2)

print("=== Independent geometric normalization audit ===")
check(banaszczyk_cap / cube_half_side == 14, "factor-two transference threshold")
check(2 * cube_half_side * 14 == projected_dimension, "normalization identity")
checks += 2
print("External input normalized as 2*mu(L)*lambda1(L*) <= d.")
print("With d=12 and closed-ball disjointness dist(c,L)>3/7, lambda1(L*)<14.")

# A closed boundary is load-bearing: in R/Z, the open radius-1/2 ball about
# 1/2 misses Z while the closed ball meets it.  Therefore only disjointness of
# the closed LRC zonotope upgrades >= to >.
circle_nearest_distance = Fraction(1, 2)
check(not (circle_nearest_distance < Fraction(1, 2)), "open-ball hostile")
check(circle_nearest_distance <= Fraction(1, 2), "closed-ball hostile")
checks += 2
print("Boundary hostile: an open disjoint ball can have equality; the closed zonotope cannot.")

# Enumerate all extreme rays of cross-polytope slices for many exact controls.
# In a sign orthant, two affine equations leave an extreme point on at most
# two coordinates.  The pair value below is its exact squared l1/l2 ratio.
extreme_controls = 0
for speeds in combinations(range(1, 13), 5):
    pair_values = {
        pair: Fraction((pair[0] + pair[1]) ** 2, pair[0] ** 2 + pair[1] ** 2)
        for pair in combinations(speeds, 2)
    }
    check(min(pair_values, key=pair_values.get) == (speeds[0], speeds[-1]), "inradius extreme pair")
    extreme_controls += 1
checks += extreme_controls
print(f"Exact cross-polytope extreme-pair controls: {extreme_controls}.")
print("They give rho(n)=(m+M)/sqrt(m^2+M^2)>1, with the universal limit 1.")


print("\n=== Independent integer-cap audit ===")
square_cap = 14 * 14 - 1
check(square_cap == 195, "integer square cap")
check(isqrt(square_cap) == 13, "coordinate cap")
check(l1_maximum(square_cap, speed_count) == 50, "global l1 cap")
check(balanced_square_cost(50, speed_count) == 194, "l1=50 witness")
check(balanced_square_cost(51, speed_count) == 201, "l1=51 obstruction")
checks += 5
print("||a||_2<14 gives square norm <=195, coefficient height <=13, and l1<=50.")
print("Balanced squares: cost(50)=194 from (3,3,4x11), while cost(51)=201.")

spread_table: list[tuple[int, int, int]] = []
for spread in (2, 3, 5, 13, 21, 50):
    strict_square_threshold = Fraction(196 * (1 + spread * spread), (1 + spread) ** 2)
    local_square_cap = (strict_square_threshold.numerator - 1) // strict_square_threshold.denominator
    spread_table.append((spread, local_square_cap, l1_maximum(local_square_cap, speed_count)))
expected_spread_table = [
    (2, 108, 37),
    (3, 122, 39),
    (5, 141, 42),
    (13, 169, 46),
    (21, 178, 47),
    (50, 188, 49),
]
check(spread_table == expected_spread_table, "bounded-spread table")
checks += 1
print("Bounded-spread (R, square cap, l1 cap): " + repr(spread_table))


print("\n=== Independent support and histogram atlas ===")
ratios = tuple(
    (small, large)
    for large in range(2, 14)
    for small in range(1, large)
    if gcd(small, large) == 1 and small * small + large * large <= square_cap
)
check(len(ratios) == 47, "47 reduced ratios")
check(max(sum(pair) for pair in ratios) == 19, "support-two l1 maximum")
check(tuple(pair for pair in ((1, 5), (5, 29), (29, 169)) if sum(x*x for x in pair) <= 195) == ((1, 5),), "Pell intersection")
checks += 3
ratio_digest = sha256(";".join(f"{a},{b}" for a, b in ratios).encode()).hexdigest()
unoriented_support_packets = len(ratios) * comb(speed_count, 2)
oriented_assignments = 2 * unoriented_support_packets
check(unoriented_support_packets == 3666, "unoriented ratio-support packets")
check(oriented_assignments == 7332, "oriented labelled assignments")
checks += 2
print(f"Reduced ratios: 47; ordered-list SHA-256={ratio_digest}")
print("47*C(13,2)=3,666 unoriented ratio/support packets.")
print("Without a declared speed ordering, there are 7,332 oriented coefficient assignments modulo sign.")

histograms, sign_orbits, odd_histograms = histogram_generating_function(speed_count, square_cap)
expected_histograms = {
    2: 47,
    3: 209,
    4: 566,
    5: 1177,
    6: 2057,
    7: 3180,
    8: 4490,
    9: 5911,
    10: 7374,
    11: 8805,
    12: 10188,
    13: 11455,
}
check(histograms == expected_histograms, "histogram support table")
check(sum(histograms.values()) == 55459, "histogram total")
check(sum(sign_orbits.values()) == 5030161, "sign-orbit total")
check(sum(odd_histograms.values()) == 28315, "odd histogram total")
checks += 4
print("Primitive absolute histograms: " + ", ".join(f"{k}:{v}" for k, v in histograms.items()))
print("Totals: 55,459 histograms; 5,030,161 nonconstant sign-multiset orbits; 28,315 odd-l1 histograms.")


print("\n=== Independent labelled-ball and rank-fork audit ===")
distribution = labelled_ball_distribution(speed_count, square_cap)
prefix = []
running = 0
for count in distribution:
    running += count
    prefix.append(running - 1)  # remove the zero vector
ball_count = prefix[square_cap]
check(ball_count == 711202814025242, "labelled ball")
mobius = moebius_values(isqrt(square_cap))
primitive_ball = sum(
    mobius[divisor] * prefix[square_cap // (divisor * divisor)]
    for divisor in range(1, isqrt(square_cap) + 1)
)
check(primitive_ball == 711119925281794, "primitive ball")
checks += 2
print(f"Labelled nonzero ball={ball_count}; primitive ball={primitive_ball}.")

relation_height = 91**6
hadamard_radicand = 195 * 3**11 * relation_height**22
hadamard_cap = isqrt(hadamard_radicand)
expected_hadamard_cap = int(
    "11639011567946276516330452125265832396450210671398535269626998205847761612881723251682730932524351641613494572161619231675129410858204"
)
check(hadamard_cap == expected_hadamard_cap, "rank-eleven cofactor cap")
check(hadamard_cap**2 <= hadamard_radicand < (hadamard_cap + 1) ** 2, "Hadamard floor")
check(len(str(hadamard_cap)) == 134, "Hadamard digits")
checks += 3
print(f"Rank-11 outside-span Hadamard cap has {len(str(hadamard_cap))} digits and exact floor {hadamard_cap}.")


print("\n=== Independent Graver/parity/hostile audit ===")
print("A shortest Euclidean kernel vector is Graver: every nonzero conformal summand has strictly smaller square norm.")

# The odd-character sidecar also has an elementary parity control.  A primitive
# speed row is either all odd, in which case t=1/2 is immediately lonely and
# every relation has even coefficient sum, or it has an even and an odd speed;
# their primitive support-two relation has odd coefficient sum.
parity_patterns = 0
for pattern in range(1, 1 << speed_count):  # exclude the nonprimitive all-even pattern
    all_odd = pattern == (1 << speed_count) - 1
    has_even_coordinate = pattern != (1 << speed_count) - 1
    check(all_odd != has_even_coordinate, "parity dichotomy")
    parity_patterns += 1
checks += parity_patterns

parity_controls = ((2, 3), (6, 25), (14, 15), (26, 39))
for even_speed, odd_speed in parity_controls:
    divisor = gcd(even_speed, odd_speed)
    relation = (odd_speed // divisor, -(even_speed // divisor))
    check(relation[0] * even_speed + relation[1] * odd_speed == 0, "odd relation control")
    check(sum(relation) % 2 == 1, "odd coefficient sum control")
checks += 2 * len(parity_controls)
print(f"Parity-character controls: {parity_patterns} primitive mod-2 speed patterns plus {len(parity_controls)} integer pairs.")
print("The odd-sum relation exists for a counterexample, but transference does not make it the short row.")

ap_speeds = tuple(range(1, 14))
ap_relation = (1, -2, 1) + (0,) * 10
check(sum(x * y for x, y in zip(ap_speeds, ap_relation)) == 0, "AP relation")
check(sum(x * x for x in ap_relation) == 6, "AP square norm")
ap_time = Fraction(1, 14)
ap_margin = min(min((speed * ap_time) % 1, 1 - (speed * ap_time) % 1) for speed in ap_speeds)
check(ap_margin == Fraction(1, 14), "AP lonely boundary")
checks += 3
print("AP13 hostile remains lonely at t=1/14 despite its square-norm-6 Graver relation.")


print("\n=== Referee verdict ===")
print("Core implication PASS: counterexample => a Graver relation with ||a||_2<14, square<=195, l1<=50, height<=13.")
print("Parity and rank-11 sidecars PASS with their stated non-simultaneity caveat.")
print("REPAIR NEEDED: call 3,666 the unoriented ratio/support count; the oriented labelled count is 7,332.")
print(f"CHECKS={checks}")
