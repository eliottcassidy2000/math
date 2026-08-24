#!/usr/bin/env python3
"""Exact arithmetic audit for THM-4009.

The cited input is Banaszczyk's Euclidean transference inequality

    2 mu(Lambda) lambda_1(Lambda*) <= d.

For the 13-speed lonely-runner zonotope, the in-repo algebra proves that a
counterexample makes a translate of the closed radius-3/7 Euclidean ball in
the 12-dimensional projected space disjoint from Lambda.  Thus mu>3/7 and
lambda_1(Lambda*)<14.  This script independently audits every integer
consequence, the exact coefficient-shape atlases, the rank-eleven Hadamard
terminal, and hostile controls.  It does not prove the cited transference
theorem or the lonely-runner zonotope equivalence.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import comb, gcd, isqrt


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def mobius(value: int) -> int:
    """Return the exact Moebius function for a positive integer."""

    remaining = value
    prime_count = 0
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            exponent = 0
            while remaining % divisor == 0:
                remaining //= divisor
                exponent += 1
            if exponent > 1:
                return 0
            prime_count += 1
        divisor += 1
    if remaining > 1:
        prime_count += 1
    return -1 if prime_count % 2 else 1


def lattice_ball_count(dimension: int, square_budget: int) -> int:
    """Count nonzero labelled integer vectors with squared norm at most B."""

    counts = [0] * (square_budget + 1)
    counts[0] = 1
    coordinate_cap = isqrt(square_budget)
    for _ in range(dimension):
        next_counts = [0] * (square_budget + 1)
        for used_square, multiplicity in enumerate(counts):
            if multiplicity == 0:
                continue
            for coordinate in range(-coordinate_cap, coordinate_cap + 1):
                new_square = used_square + coordinate * coordinate
                if new_square <= square_budget:
                    next_counts[new_square] += multiplicity
        counts = next_counts
    return sum(counts) - 1


def primitive_lattice_ball_count(dimension: int, square_budget: int) -> int:
    """Count nonzero primitive vectors by exact Moebius inversion."""

    return sum(
        mobius(divisor)
        * lattice_ball_count(dimension, square_budget // (divisor * divisor))
        for divisor in range(1, isqrt(square_budget) + 1)
        if mobius(divisor) != 0
    )


def primitive_absolute_histograms(
    support: int, square_budget: int
) -> list[tuple[int, ...]]:
    """Primitive nondecreasing positive coefficient shapes of fixed support."""

    answer: list[tuple[int, ...]] = []

    def recurse(
        lower: int,
        places_left: int,
        budget_left: int,
        common_divisor: int,
        prefix: tuple[int, ...],
    ) -> None:
        if places_left == 0:
            if common_divisor == 1:
                answer.append(prefix)
            return
        upper = isqrt(budget_left // places_left)
        for coefficient in range(lower, upper + 1):
            recurse(
                coefficient,
                places_left - 1,
                budget_left - coefficient * coefficient,
                gcd(common_divisor, coefficient),
                prefix + (coefficient,),
            )

    recurse(1, support, square_budget, 0, ())
    return answer


def exact_l1_cap(dimension: int, square_budget: int) -> tuple[int, tuple[int, ...]]:
    """Maximize l1 over sorted absolute coordinates under a square budget."""

    states: dict[tuple[int, int], tuple[int, ...]] = {(0, 0): ()}
    coordinate_cap = isqrt(square_budget)
    for _ in range(dimension):
        next_states: dict[tuple[int, int], tuple[int, ...]] = {}
        for (used_square, used_sum), values in states.items():
            for coordinate in range(coordinate_cap + 1):
                new_square = used_square + coordinate * coordinate
                if new_square <= square_budget:
                    next_states.setdefault(
                        (new_square, used_sum + coordinate), values + (coordinate,)
                    )
        states = next_states
    ell_one_cap = max(total for _, total in states)
    witness = min(
        values for (_, total), values in states.items() if total == ell_one_cap
    )
    return ell_one_cap, tuple(sorted(witness))


def sign_multiset_orbit_count(histogram: tuple[int, ...]) -> int:
    """Count nonconstant sign multisets modulo global sign and equal entries."""

    multiplicities = tuple(Counter(histogram).values())
    raw = 1
    for multiplicity in multiplicities:
        raw *= multiplicity + 1
    complement_fixed = int(all(multiplicity % 2 == 0 for multiplicity in multiplicities))
    # Burnside for sign-complement, then remove the all-positive/all-negative orbit.
    return (raw + complement_fixed) // 2 - 1


checks = 0

print("=== Euclidean inball and transference constants ===")
speed_count = 13
dimension = speed_count - 1
half_side = Fraction(3, 7)
transference_product_cap = Fraction(dimension, 2)
critical_dual_length = transference_product_cap / half_side
require(dimension == 12, "projected dimension")
require(critical_dual_length == 14, "exact dual-length threshold")
require(2 * half_side * critical_dual_length == dimension, "transference normalization")
checks += 3

print("[1/14,13/14]^13 = (1/2)1 + (3/7)[-1,1]^13.")
print("For x in n^perp with ||x||_2<=1, every |x_i|<=1 and projection(x)=x;")
print("hence the projected zonotope contains the closed ball B(c,3/7).")
print("A counterexample has dist(c,Lambda)>3/7, so mu(Lambda)>3/7.")
print("Using 2 mu(Lambda) lambda_1(Lambda*)<=12 gives lambda_1(Lambda*)<14.")

spread_controls = (
    (1, 2, 3, 4),
    (1, 3, 7, 13),
    (2, 5, 11, 21),
    (7, 8, 9, 10, 11, 12, 13),
)
for speeds in spread_controls:
    pair_inradius_squares = {
        (first, second): Fraction((first + second) ** 2, first * first + second * second)
        for index, first in enumerate(speeds)
        for second in speeds[index + 1 :]
    }
    minimizing_pair = min(pair_inradius_squares, key=pair_inradius_squares.get)
    require(minimizing_pair == (min(speeds), max(speeds)), "extreme-pair inradius control")
checks += len(spread_controls)

print("More exactly, the projected cube's centred Euclidean inradius factor is")
print("  rho(n)=min_(i<j) (n_i+n_j)/sqrt(n_i^2+n_j^2)")
print("        =(n_min+n_max)/sqrt(n_min^2+n_max^2).")
print("Thus lambda_1<14 sqrt(n_min^2+n_max^2)/(n_min+n_max); the universal")
print("bound is its unbounded-spread limit.")


print("\n=== Exact integer coefficient consequences ===")
strict_square_threshold = critical_dual_length * critical_dual_length
require(strict_square_threshold == 196, "strict squared-norm threshold")
square_cap = int(strict_square_threshold) - 1
require(square_cap == 195, "integer squared-norm cap")
coefficient_cap = isqrt(square_cap)
require(coefficient_cap == 13, "coordinate-height cap")

ell_one_cap, ell_one_witness = exact_l1_cap(speed_count, square_cap)
require(ell_one_cap == 50, "exact l1 cap")
require(
    ell_one_witness == (3, 3) + (4,) * 11,
    "unique sorted arithmetic l1 extremizer",
)
checks += 5

print("Every transference relation a satisfies sum a_i^2<=195 and |a_i|<=13.")
print("The exact integer optimization gives ||a||_1<=50.")
print("The unique sorted absolute arithmetic extremizer is (3,3,4,...,4),")
print("with squared norm 194 and l1 norm 50 (realizability is not asserted).")

ratio_caps: list[tuple[int, Fraction, int, int]] = []
for ratio_cap in (2, 3, 5, 13, 21, 50):
    strict_square_cap = Fraction(
        196 * (1 + ratio_cap * ratio_cap), (1 + ratio_cap) ** 2
    )
    integer_square_cap = (strict_square_cap.numerator - 1) // strict_square_cap.denominator
    local_l1_cap, _ = exact_l1_cap(speed_count, integer_square_cap)
    ratio_caps.append(
        (ratio_cap, strict_square_cap, integer_square_cap, local_l1_cap)
    )
require(
    ratio_caps
    == [
        (2, Fraction(980, 9), 108, 37),
        (3, Fraction(245, 2), 122, 39),
        (5, Fraction(1274, 9), 141, 42),
        (13, Fraction(170, 1), 169, 46),
        (21, Fraction(21658, 121), 178, 47),
        (50, Fraction(490196, 2601), 188, 49),
    ],
    "bounded-spread cap table",
)
checks += 1

print("Bounded-spread integer caps (R=n_max/n_min):")
print("  R : square-norm cap, l1 cap")
for ratio_cap, _, local_square_cap, local_l1_cap in ratio_caps:
    print(f"  {ratio_cap}: {local_square_cap}, {local_l1_cap}")


print("\n=== Primitive relation atlases ===")
pair_ratios = [
    (first, second)
    for first in range(1, coefficient_cap + 1)
    for second in range(first + 1, coefficient_cap + 1)
    if gcd(first, second) == 1
    and first * first + second * second <= square_cap
]
require(len(pair_ratios) == 47, "support-two ratio count")
require(max(first + second for first, second in pair_ratios) == 19, "pair l1 cap")
unoriented_pair_support_packets = len(pair_ratios) * comb(speed_count, 2)
oriented_pair_assignments = 2 * unoriented_pair_support_packets
require(unoriented_pair_support_packets == 3666, "unoriented support-two packet count")
require(oriented_pair_assignments == 7332, "oriented support-two assignment count")
pell_pairs = ((1, 5), (5, 29), (29, 169))
pell_survivors = tuple(
    pair for pair in pell_pairs if pair[0] * pair[0] + pair[1] * pair[1] <= square_cap
)
require(pell_survivors == ((1, 5),), "Pell-selector intersection")
checks += 5

histogram_counts: dict[int, int] = {}
sign_orbit_counts: dict[int, int] = {}
odd_sum_histogram_counts: dict[int, int] = {}
for support in range(2, speed_count + 1):
    histograms = primitive_absolute_histograms(support, square_cap)
    if support == 2:
        # (1,1) cannot annihilate two distinct positive speeds.
        histograms.remove((1, 1))
    histogram_counts[support] = len(histograms)
    sign_orbit_counts[support] = sum(
        sign_multiset_orbit_count(histogram) for histogram in histograms
    )
    odd_sum_histogram_counts[support] = sum(
        sum(histogram) % 2 for histogram in histograms
    )

require(sum(histogram_counts.values()) == 55459, "absolute histogram atlas")
require(sum(sign_orbit_counts.values()) == 5030161, "signed-multiset orbit atlas")
require(sum(odd_sum_histogram_counts.values()) == 28315, "odd-l1 histogram atlas")
checks += 3

print(f"Support-two reduced ratios p<q with p^2+q^2<=195: {len(pair_ratios)}.")
print(f"Unoriented ratio/support packets: {unoriented_pair_support_packets}.")
print(f"Oriented labelled coefficient assignments: {oriented_pair_assignments}.")
print("Their maximum p+q is 19; the prior Pell packet leaves only (1,5).")
print("Primitive absolute coefficient histograms by support:")
print("  " + ", ".join(f"{support}:{histogram_counts[support]}" for support in histogram_counts))
print(f"Total absolute histograms: {sum(histogram_counts.values())}.")
print(f"Signed-multiset types modulo permutation/global sign: {sum(sign_orbit_counts.values())}.")
print(f"Absolute histograms with odd coefficient sum: {sum(odd_sum_histogram_counts.values())}.")


print("\n=== Labelled Euclidean-ball universe ===")
ball_count = lattice_ball_count(speed_count, square_cap)
primitive_ball_count = primitive_lattice_ball_count(speed_count, square_cap)
require(ball_count == 711202814025242, "labelled Euclidean-ball count")
require(primitive_ball_count == 711119925281794, "primitive Euclidean-ball count")
primitive_support_at_least_two = primitive_ball_count - 2 * speed_count
require(primitive_support_at_least_two == 711119925281768, "remove support-one vectors")
orientation_quotient = primitive_support_at_least_two // 2
require(orientation_quotient == 355559962640884, "global-sign quotient")

old_l1_universe = sum(
    2**support * comb(speed_count, support) * comb(356, support)
    for support in range(1, speed_count + 1)
)
require(old_l1_universe == 1978967793896659449022201064, "THM-3743 universe replay")
reduction_floor = old_l1_universe // ball_count
require(reduction_floor == 2782564628359, "universe reduction floor")
checks += 6

print(f"All nonzero labelled integer vectors with sum a_i^2<=195: {ball_count}.")
print(f"Primitive vectors: {primitive_ball_count}; primitive support>=2 modulo sign: {orientation_quotient}.")
print(f"This is more than {reduction_floor} times smaller than THM-3743's l1<=356 ambient universe.")


print("\n=== Rank-eleven Hadamard join ===")
relation_height = 91**6
mixed_cofactor_square_bound = square_cap * 3**11 * relation_height**22
mixed_cofactor_cap = isqrt(mixed_cofactor_square_bound)
require(
    mixed_cofactor_cap
    == 11639011567946276516330452125265832396450210671398535269626998205847761612881723251682730932524351641613494572161619231675129410858204,
    "mixed Hadamard cofactor cap",
)
require(
    mixed_cofactor_cap * mixed_cofactor_cap
    <= mixed_cofactor_square_bound
    < (mixed_cofactor_cap + 1) * (mixed_cofactor_cap + 1),
    "exact square-root floor",
)
require(len(str(mixed_cofactor_cap)) == 134, "Hadamard digit count")
checks += 3

print("If the short row lies outside THM-2052's rank-11 span, Hadamard gives")
print(f"  max(n_i) <= {mixed_cofactor_cap}")
print("(134 digits). If it lies inside, that span contains a row with squared norm<=195.")
print("Support at most three is automatically inside because coefficient height<=13<91^6.")


print("\n=== Half-lattice parity sidecar and safe hostile ===")
print("The inball center c=(1/2)projection(1) satisfies 2c in Lambda.")
print("Because a counterexample has c notin Lambda=Lambda**, some (not necessarily short)")
print("dual relation has <a,c>=1/2 mod 1, equivalently sum a_i is odd.")
print("Banaszczyk transference alone does not force the norm-<14 row to be this odd row.")

ap_speeds = tuple(range(1, speed_count + 1))
ap_relation = (1, -2, 1) + (0,) * (speed_count - 3)
require(sum(a * n for a, n in zip(ap_relation, ap_speeds)) == 0, "AP relation")
require(sum(a * a for a in ap_relation) == 6, "AP relation squared norm")
lonely_time = Fraction(1, 14)
ap_minimum = min(
    min((speed * lonely_time) % 1, 1 - ((speed * lonely_time) % 1))
    for speed in ap_speeds
)
require(ap_minimum == Fraction(1, 14), "AP boundary lonely witness")
checks += 3

print("The safe AP (1,...,13) has relation (1,-2,1), squared norm 6, and")
print("a valid boundary loneliness time t=1/14. Short resonance remains necessary,")
print("never sufficient; owner, phase, sign partition, and arrival are still lost.")


print("\n=== Scope verdict ===")
print("CITED + PROVED-ALGEBRA target: counterexample => a nonzero Graver relation")
print("with ||a||_2<14, sum a_i^2<=195, ||a||_1<=50, and |a_i|<=13.")
print("The citation supplies only Euclidean transference; all printed finite constants are exact.")
print("This does not prove LRC(14).")
print(f"CHECKS={checks}")
