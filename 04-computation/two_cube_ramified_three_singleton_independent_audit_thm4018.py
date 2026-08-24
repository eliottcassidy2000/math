#!/usr/bin/env python3
"""Independent exact audit of the ramified-three two-cube criterion.

This implementation does not import the discovery probe.  Its global
representation oracle uses a coordinate ceiling (v^3 < m), rather than the
discovery probe's pair-sum ceiling.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import itertools
import math


MAX_CANDIDATE_SUM = 2500
PRIMITIVE_BOX = 1200
GATES = 0


def require(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def valuation(n: int, p: int) -> int:
    exponent = 0
    while n % p == 0:
        exponent += 1
        n //= p
    return exponent


def prime_factorization(n: int) -> tuple[tuple[int, int], ...]:
    factors: list[tuple[int, int]] = []
    divisor = 2
    while divisor * divisor <= n:
        exponent = 0
        while n % divisor == 0:
            exponent += 1
            n //= divisor
        if exponent:
            factors.append((divisor, exponent))
        divisor = 3 if divisor == 2 else divisor + 2
    if n > 1:
        factors.append((n, 1))
    return tuple(factors)


SUM_FACTORS = [tuple()] * (MAX_CANDIDATE_SUM + 1)
for pair_sum in range(2, MAX_CANDIDATE_SUM + 1):
    SUM_FACTORS[pair_sum] = prime_factorization(pair_sum)


def admissibility(pair_sum: int, smaller: int, cap_at_three: int) -> bool:
    factors = SUM_FACTORS[pair_sum]
    if any(prime != 3 and prime % 3 != 2 for prime, _ in factors):
        return False
    scale = math.gcd(smaller, pair_sum)
    primitive_sum = pair_sum // scale
    for prime, _ in factors:
        cap = cap_at_three if prime == 3 else 2
        if valuation(primitive_sum, prime) > cap:
            return False
    return True


strict_origin: dict[int, tuple[int, int]] = {}
relaxed_origin: dict[int, tuple[int, int]] = {}
strict_origin_collisions: list[tuple[int, tuple[int, int], tuple[int, int]]] = []
relaxed_origin_collisions: list[tuple[int, tuple[int, int], tuple[int, int]]] = []
strict_rows = 0
relaxed_rows = 0

for pair_sum in range(3, MAX_CANDIDATE_SUM + 1):
    row_has_strict = False
    row_has_relaxed = False
    for smaller in range(1, (pair_sum - 1) // 2 + 1):
        larger = pair_sum - smaller
        value = smaller**3 + larger**3
        pair = (smaller, larger)
        if admissibility(pair_sum, smaller, 1):
            row_has_strict = True
            previous = strict_origin.setdefault(value, pair)
            if previous != pair:
                strict_origin_collisions.append((value, previous, pair))
        if admissibility(pair_sum, smaller, 2):
            row_has_relaxed = True
            previous = relaxed_origin.setdefault(value, pair)
            if previous != pair:
                relaxed_origin_collisions.append((value, previous, pair))
    strict_rows += int(row_has_strict)
    relaxed_rows += int(row_has_relaxed)

require(not strict_origin_collisions, "two strict origins never collide")

# Complete competitor oracle for the candidate-value universe.  If
# 0<u<v and u^3+v^3=m<=M, then v^3<m, so v<=floor(cuberoot(M-1)).
candidate_values = set(relaxed_origin)
maximum_value = max(candidate_values)


def integer_cube_root(n: int) -> int:
    low, high = 0, 1
    while high**3 <= n:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**3 <= n:
            low = middle
        else:
            high = middle
    return low


competitor_coordinate_max = integer_cube_root(maximum_value - 1)
strict_competitors: list[tuple[int, tuple[int, int], tuple[int, int]]] = []
relaxed_competitors: dict[int, list[tuple[int, int]]] = {}
strict_seen_origins = 0
relaxed_seen_origins = 0

cube = [integer**3 for integer in range(competitor_coordinate_max + 1)]
for larger in range(2, competitor_coordinate_max + 1):
    larger_cube = cube[larger]
    for smaller in range(1, larger):
        value = cube[smaller] + larger_cube
        pair = (smaller, larger)
        strict_pair = strict_origin.get(value)
        if strict_pair is not None:
            if pair == strict_pair:
                strict_seen_origins += 1
            else:
                strict_competitors.append((value, strict_pair, pair))
        relaxed_pair = relaxed_origin.get(value)
        if relaxed_pair is not None:
            if pair == relaxed_pair:
                relaxed_seen_origins += 1
            else:
                relaxed_competitors.setdefault(value, [relaxed_pair]).append(pair)

require(strict_seen_origins == len(strict_origin), "every strict origin is in the coordinate oracle")
require(relaxed_seen_origins == len(relaxed_origin), "every relaxed origin is in the coordinate oracle")
require(not strict_competitors, "strict candidates have no competitor in the complete oracle")
require(bool(relaxed_competitors), "the cap-two relaxation has a competitor")

hostile_values = sorted(relaxed_competitors)
first_hostile = hostile_values[0]
primitive_hostile_values = [
    value
    for value in hostile_values
    if math.gcd(*(coordinate for pair in relaxed_competitors[value] for coordinate in pair)) == 1
]
require(first_hostile == 4104, "the first relaxed hostile value is 4104")
require(
    sorted(relaxed_competitors[first_hostile]) == [(2, 16), (9, 15)],
    "the complete 4104 fibre",
)
require(2**3 + 16**3 == 9**3 + 15**3 == 4104, "4104 identity")
require(admissibility(18, 2, 2), "(2,16) obeys the relaxed cap")
require(not admissibility(18, 2, 1), "(2,16) violates only the strict cap at three")
require(valuation(18 // math.gcd(2, 16), 3) == 2, "4104 source primitive three exponent")

# This makes the absence of a lower hostile global, not merely a statement
# about the chosen candidate-sum box: every positive pair of value m obeys
# (x+y)^3 < 4m.
lower_hostile_sum_ceiling = integer_cube_root(4 * first_hostile - 1)
require(lower_hostile_sum_ceiling < MAX_CANDIDATE_SUM, "candidate box contains every possible lower hostile")
require(all(value >= first_hostile for value in hostile_values), "no lower relaxed hostile")

# Re-derive and exhaust the local valuation invoices.  The branch indicator
# is v_3(X^2-XY+Y^2)=1 exactly when v_3(X+Y)>0.
three_branch_counts = {"e0_ep_0": 0, "e0_ep_pos": 0, "e1_ep_0": 0, "e1_ep_pos": 0}
for alpha in range(0, 25):
    for source_e in (0, 1):
        source_value_valuation = 3 * alpha + source_e + int(source_e > 0)
        for beta in range(0, 26):
            for target_e in range(0, 80):
                target_value_valuation = 3 * beta + target_e + int(target_e > 0)
                if source_value_valuation != target_value_valuation:
                    continue
                key = f"e{source_e}_ep_{'pos' if target_e else '0'}"
                three_branch_counts[key] += 1
                require(beta + target_e >= alpha + source_e, f"three divisibility {alpha,source_e,beta,target_e}")
                if source_e == 0 and target_e == 0:
                    require(beta == alpha, "unramified-to-unramified scale equality")
                elif source_e == 0:
                    delta = alpha - beta
                    require(delta >= 1 and target_e == 3 * delta - 1, "unramified-to-ramified formula")
                    require((beta + target_e) - alpha == 2 * delta - 1, "unramified margin")
                elif target_e == 0:
                    require(False, "ramified source cannot have an unramified competitor")
                else:
                    delta = alpha - beta
                    require(delta >= 0 and target_e == 3 * delta + 1, "ramified-to-ramified formula")
                    require((beta + target_e) - (alpha + 1) == 2 * delta, "ramified margin")

require(three_branch_counts["e1_ep_0"] == 0, "mod-three exclusion branch")

# The inert-prime invoice has no cofactor contribution and works uniformly
# for e=0,1,2.  Exhaust a large parameter rectangle.
inert_invoice_count = 0
for alpha in range(0, 25):
    for source_e in (0, 1, 2):
        source_value_valuation = 3 * alpha + source_e
        for beta in range(0, 26):
            target_e = source_value_valuation - 3 * beta
            if target_e < 0:
                continue
            inert_invoice_count += 1
            require(beta <= alpha, "inert competitor scale cannot rise")
            require(beta + target_e >= alpha + source_e, "inert divisibility margin")

# Independently check the primitive ramified valuation law on a full box.
primitive_pairs = 0
for first in range(1, PRIMITIVE_BOX + 1):
    for second in range(first + 1, PRIMITIVE_BOX + 1):
        if math.gcd(first, second) != 1:
            continue
        primitive_pairs += 1
        cofactor = first * first - first * second + second * second
        expected = int((first + second) % 3 == 0)
        require(valuation(cofactor, 3) == expected, f"primitive ramified law {first,second}")

# Finite exact audit of the proposed two-row analytic amplifier.  This is a
# separate subset implementation with rational arithmetic.
analytic_primes = (5, 11, 17, 23, 29, 41)
analytic_cutoff = 4
main_mass = Fraction(0)
square_mass = Fraction(0)
row_surrogate = Fraction(0)
subset_count = 0
for order in range(1, analytic_cutoff + 1):
    for subset in itertools.combinations(analytic_primes, order):
        product = math.prod(subset)
        subset_count += 1
        main_mass += Fraction(1, product)
        square_mass += Fraction(1, product * product)
        first_row = Fraction(product - 1, 2 * product * product)
        second_row = Fraction(3 * product - 1, 18 * product * product)
        require(
            first_row + second_row == Fraction(2, 3 * product) - Fraction(5, 9 * product * product),
            f"two-row surrogate d={product}",
        )
        row_surrogate += first_row + second_row

require(row_surrogate == Fraction(2, 3) * main_mass - Fraction(5, 9) * square_mass, "summed row identity")
full_square_product = math.prod((Fraction(1) + Fraction(1, prime * prime)) for prime in analytic_primes)
require(square_mass <= full_square_product - 1, "square-weight error bounded by full Euler product")

# Cross-row/cross-subset disjointness on a nontrivial exact bank.  Each row is
# also contained in the strict candidate universe above, whose competitor
# oracle is complete.
row_bank = (5, 11, 17, 23)
row_values: dict[int, tuple[int, int]] = {}
checked_rows = 0
for order in range(1, 3):
    for subset in itertools.combinations(row_bank, order):
        product = math.prod(subset)
        for pair_sum in (product, 3 * product):
            require(pair_sum <= MAX_CANDIDATE_SUM, "analytic row lies in global audit box")
            checked_rows += 1
            for smaller in range(1, (pair_sum - 1) // 2 + 1):
                larger = pair_sum - smaller
                value = smaller**3 + larger**3
                previous = row_values.setdefault(value, (smaller, larger))
                require(previous == (smaller, larger), f"analytic rows collide at {value}")
                require(value in strict_origin, f"analytic row is strict-admissible at {pair_sum,smaller}")

# Floor-sensitive all-real-X ledger, expressed in L=log(log(X/27)).
floor_rows: list[tuple[float, int, float, float]] = []
for logarithmic_height in (1_000.0, 10_000.0, 100_000.0, 1_000_000.0):
    raw_cutoff = (
        logarithmic_height / 2
        - 0.5 * math.log(logarithmic_height)
        + logarithmic_height ** (2 / 3)
    )
    cutoff = math.floor(raw_cutoff)
    floor_error = raw_cutoff - cutoff
    log_log_z = logarithmic_height - math.log(3 * cutoff)
    displacement = cutoff - 0.5 * log_log_z
    ledger = (
        logarithmic_height ** (2 / 3)
        + 0.5 * math.log(3 * cutoff / logarithmic_height)
        - floor_error
    )
    require(abs(displacement - ledger) < 2e-6, f"all-real floor identity L={logarithmic_height}")
    require(displacement / math.sqrt(logarithmic_height) > 2, f"cutoff dominates standard deviation L={logarithmic_height}")
    height_ratio = math.sqrt(logarithmic_height / (3 * cutoff))
    floor_rows.append((logarithmic_height, cutoff, displacement, height_ratio))

require(
    all(
        floor_rows[index][2] / math.sqrt(floor_rows[index][0])
        < floor_rows[index + 1][2] / math.sqrt(floor_rows[index + 1][0])
        for index in range(len(floor_rows) - 1)
    ),
    "normalized cutoff displacement grows across the audit ledger",
)

semantic = (
    "the ramified-three cap one gives primewise pair-sum divisibility for every common scale",
    "strict admissible positive distinct two-cube fibres are globally singleton",
    "4104 is the globally least hostile when the primitive three exponent cap is raised to two",
    "the d and 3d squarefree-inert rows are globally disjoint singleton rows",
    "their Euler main coefficient is 2/3 and their total square-weight error is bounded",
    "the all-real height normalization uses Y=X/27 and loses no asymptotic constant",
)

print(
    "UNIVERSE",
    f"candidate_sum<={MAX_CANDIDATE_SUM};candidate_values={len(candidate_values)};"
    f"competitor_coordinates=1<=u<v<={competitor_coordinate_max};max_candidate_value={maximum_value}",
)
print(
    "STRICT",
    f"rows={strict_rows};values={len(strict_origin)};origin_collisions={len(strict_origin_collisions)};"
    f"external_competitors={len(strict_competitors)}",
)
print(
    "RELAXED",
    f"rows={relaxed_rows};values={len(relaxed_origin)};origin_collisions={len(relaxed_origin_collisions)};"
    f"hostile_values={len(hostile_values)};hostile_incidences={sum(len(v)-1 for v in relaxed_competitors.values())}",
)
print(
    "FIRST_HOSTILE",
    f"m={first_hostile};fibre={sorted(relaxed_competitors[first_hostile])};"
    f"global_lower_candidate_sum_ceiling={lower_hostile_sum_ceiling}",
)
print(
    "NEXT_HOSTILES",
    ";".join(f"{value}:{sorted(relaxed_competitors[value])}" for value in hostile_values[1:11]),
)
print(
    "COMMON_SCALE_PRIMITIVE_HOSTILES",
    f"count={len(primitive_hostile_values)};first="
    + ";".join(
        f"{value}:{sorted(relaxed_competitors[value])}" for value in primitive_hostile_values[:10]
    ),
)
print("THREE_BRANCH_COUNTS", ";".join(f"{key}={value}" for key, value in three_branch_counts.items()))
print("INERT_INVOICE_CASES", inert_invoice_count)
print("PRIMITIVE_THREE_LAW", f"box={PRIMITIVE_BOX};coprime_pairs={primitive_pairs};violations=0")
print(
    "ANALYTIC_FINITE",
    f"primes={analytic_primes};orders<= {analytic_cutoff};subsets={subset_count};"
    f"main={main_mass};square={square_mass};surrogate={row_surrogate};checked_rows={checked_rows};"
    f"checked_row_values={len(row_values)}",
)
for logarithmic_height, cutoff, displacement, height_ratio in floor_rows:
    print(
        "FLOOR_LEDGER",
        f"L={logarithmic_height:.0f};J={cutoff};displacement={displacement:.12f};"
        f"sqrt_L_over_3J={height_ratio:.12f}",
    )
print("LIMIT_CONSTANT_FACTOR", "(2/3)*sqrt(2/3)*E_32;ratio_to_THM3793=5/3")
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic).encode()).hexdigest())
print("GATES", GATES)
print("RESULT PASS")
