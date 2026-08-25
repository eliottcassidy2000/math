#!/usr/bin/env python3
"""Independent exact audit for the THM-4026 Sun 2-4-6-8 candidate.

The load-bearing calculation is integer-only and dependency-free.  It uses
the identity

    1 + 8 * C(w, 2) = (2*w - 1)^2

and bit masks of quadratic-residue conditions.  In the first terminal route,
324 candidates survive primes through 89; the height test rejects 31 and
``isqrt`` rejects the remaining 293.  The second terminal route continues
the modular sieve through 137 and leaves no candidate at all.

The script also audits the finite exceptional-prime rows used in the proof of
THM-4027 (universal modular solubility), exact small-prime local-density
factors, and finite controls for the binomial-period formula.  The all-prime
period and local-solubility statements themselves are proved in the companion
report; these bounded checks are hostile controls, not proof by sampling.

Reproduction:
    python3 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
    python3 -O 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
"""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, isqrt


TARGET = 896_315_812_331_399
DEGREES = (2, 4, 6, 8)
SIEVE_PRIMES = (
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
    131, 137,
)
EXACT_TAIL_PRIME = 89

EXPECTED_TRIPLE_COUNTS = {
    3: 1_794_052_919,
    5: 990_109_497,
    7: 515_638_619,
    11: 183_814_972,
    13: 88_390_982,
    17: 38_944_175,
    19: 18_074_135,
    23: 8_289_810,
    29: 4_244_318,
    31: 2_163_573,
    37: 1_118_833,
    41: 569_601,
    43: 294_987,
    47: 148_763,
    53: 75_256,
    59: 38_100,
    61: 19_601,
    67: 9_971,
    71: 5_051,
    73: 2_580,
    79: 1_306,
    83: 675,
    89: 324,
    97: 152,
    101: 84,
    103: 45,
    107: 21,
    109: 8,
    113: 3,
    127: 3,
    131: 1,
    137: 0,
}


def require(condition: bool, message: str) -> None:
    """Optimization-safe replacement for a load-bearing assertion."""

    if not condition:
        raise RuntimeError(message)


def largest_top_index(k: int, bound: int) -> int:
    """Largest t with C(t,k) <= bound, by exact monotone binary search."""

    low, high = k - 1, k
    while comb(high, k) <= bound:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if comb(middle, k) <= bound:
            low = middle
        else:
            high = middle
    require(comb(low, k) <= bound < comb(low + 1, k), "exact top-index bound")
    return low


def floor_log_p(k: int, p: int) -> int:
    """Return floor(log_p(k)) using integer arithmetic."""

    exponent = 0
    power = p
    while power <= k:
        exponent += 1
        power *= p
    return exponent


def prime_power_period(k: int, p: int, a: int) -> int:
    """The proved least period p^(a+floor(log_p(k)))."""

    return p ** (a + floor_log_p(k, p))


def audit_period_controls() -> int:
    """Finite hostile checks of the exact period formula."""

    rows = 0
    for k in DEGREES:
        for p in (2, 3, 5, 7, 11, 13):
            for a in (1, 2, 3):
                modulus = p**a
                predicted = prime_power_period(k, p, a)
                values = [comb(t, k) % modulus for t in range(2 * predicted)]
                require(
                    all(values[t] == values[t + predicted] for t in range(predicted)),
                    "predicted binomial period",
                )
                if predicted > 1:
                    require(
                        any(
                            values[t] != values[t + predicted // p]
                            for t in range(predicted)
                        ),
                        "proper p-divisor is not a period",
                    )
                rows += 1
    return rows


def add_sets(left: set[int], right: set[int], modulus: int) -> set[int]:
    """Cyclic sumset modulo ``modulus``."""

    return {(a + b) % modulus for a in left for b in right}


def audit_exceptional_regular_rows() -> dict[int, tuple[int, ...]]:
    """Exact finite rows not discharged by the Cauchy--Davenport bound.

    The triangular coordinate is restricted to inputs with nonzero formal
    derivative 2*w-1.  Hence every represented residue is Hensel-regular in
    that coordinate.
    """

    exceptional = (11, 13, 17, 19, 23, 29, 31, 47)
    traces: dict[int, tuple[int, ...]] = {}
    for p in exceptional:
        regular_two = {
            comb(w, 2) % p for w in range(p) if (2 * w - 1) % p != 0
        }
        running = regular_two
        sizes = [len(running)]
        for k in (4, 6, 8):
            image = {comb(t, k) % p for t in range(p)}
            running = add_sets(running, image, p)
            sizes.append(len(running))
        require(running == set(range(p)), "exceptional prime regular coverage")
        traces[p] = tuple(sizes)
    return traces


def cyclic_convolution(left: list[int], right: list[int]) -> list[int]:
    """Exact cyclic convolution of equally sized integer arrays."""

    modulus = len(left)
    result = [0] * modulus
    left_support = [(i, value) for i, value in enumerate(left) if value]
    right_support = [(i, value) for i, value in enumerate(right) if value]
    for i, a in left_support:
        for j, b in right_support:
            result[(i + j) % modulus] += a * b
    return result


def local_factor_small_prime(p: int, a: int) -> Fraction:
    """Normalized exact target density at p^a using true sequence periods."""

    modulus = p**a
    distribution = [1] + [0] * (modulus - 1)
    denominator = 1
    for k in DEGREES:
        period = prime_power_period(k, p, a)
        atom = [0] * modulus
        for t in range(period):
            atom[comb(t, k) % modulus] += 1
        distribution = cyclic_convolution(distribution, atom)
        denominator *= period
    return Fraction(modulus * distribution[TARGET % modulus], denominator)


def local_factor_large_prime(p: int) -> Fraction:
    """Normalized exact target density modulo a prime p > 8."""

    atoms = []
    for k in DEGREES:
        atom = [0] * p
        for t in range(p):
            atom[comb(t, k) % p] += 1
        atoms.append(atom)
    distribution = [1] + [0] * (p - 1)
    for atom in atoms:
        distribution = cyclic_convolution(distribution, atom)
    return Fraction(distribution[TARGET % p], p**3)


def is_critical_mod_prime(t: int, k: int, p: int) -> bool:
    """Whether the formal derivative of C(t,k) vanishes modulo p>k."""

    derivative_numerator = 0
    for omitted in range(k):
        term = 1
        for j in range(k):
            if j != omitted:
                term = term * (t - j) % p
        derivative_numerator = (derivative_numerator + term) % p
    return derivative_numerator == 0


def critical_tuple_audit(p: int) -> tuple[int, int, Fraction]:
    """Count fully critical target tuples and their first lift.

    If no fully critical tuple lifts from p to p^2, all surviving p^2
    solutions are nonsingular and the displayed p-adic factor is stable at
    every subsequent level.
    """

    critical_inputs = [
        [t for t in range(p) if is_critical_mod_prime(t, k, p)]
        for k in DEGREES
    ]
    target_tuples = []
    lifting_tuples = []
    for inputs in product(*critical_inputs):
        value = sum(comb(t, k) for t, k in zip(inputs, DEGREES))
        if value % p == TARGET % p:
            target_tuples.append(inputs)
            if value % (p * p) == TARGET % (p * p):
                lifting_tuples.append(inputs)
    first_level = local_factor_large_prime(p)
    representation_count = first_level.numerator * p**3 // first_level.denominator
    require(
        Fraction(representation_count, p**3) == first_level,
        "first-level representation count recovery",
    )
    stable_factor = Fraction(
        representation_count - len(target_tuples), p**3
    ) + Fraction(len(lifting_tuples), p**2)
    return len(target_tuples), len(lifting_tuples), stable_factor


def build_quadratic_residue_masks(
    x_values: range, primes: tuple[int, ...]
) -> dict[int, tuple[int, ...]]:
    """Build masks indexed by p and the y,z contribution modulo p."""

    discriminant_base = 8 * TARGET + 1
    byte_count = (len(x_values) + 7) // 8
    masks: dict[int, tuple[int, ...]] = {}
    for p in primes:
        quadratic_residues = {a * a % p for a in range(p)}
        rows = [bytearray(byte_count) for _ in range(p)]
        for index, x in enumerate(x_values):
            base = (discriminant_base - 8 * comb(x, 4)) % p
            for square in quadratic_residues:
                contribution = (base - square) % p
                rows[contribution][index >> 3] |= 1 << (index & 7)
        masks[p] = tuple(int.from_bytes(row, "little") for row in rows)
    return masks


def audit_counterexample() -> tuple[dict[int, int], int, int, int, int, int]:
    """Run both terminal routes over the complete canonical universe."""

    x_max = largest_top_index(4, TARGET)
    y_max = largest_top_index(6, TARGET)
    z_max = largest_top_index(8, TARGET)
    require((x_max, y_max, z_max) == (12_112, 932, 281), "canonical maxima")

    # For k>2 all legal top indices below k have value zero, so k-1 is the
    # canonical representative of that zero fibre.  The w coordinate starts
    # at 2 and is recovered from the discriminant rather than enumerated.
    x_values = range(3, x_max + 1)
    y_values = range(5, y_max + 1)
    z_values = range(7, z_max + 1)
    universe = len(x_values) * len(y_values) * len(z_values)
    require(universe == 3_090_472_000, "complete canonical triple universe")

    masks = build_quadratic_residue_masks(x_values, SIEVE_PRIMES)
    full_x_mask = (1 << len(x_values)) - 1
    triple_counts = {p: 0 for p in SIEVE_PRIMES}
    terminal_candidates = 0
    height_rejections = 0
    exact_tail_checks = 0
    exact_tail_solutions = 0
    final_survivors = 0

    for y in y_values:
        y_atom = comb(y, 6)
        for z in z_values:
            high_atom = y_atom + comb(z, 8)
            surviving = full_x_mask
            for p in SIEVE_PRIMES:
                contribution = (8 * high_atom) % p
                surviving &= masks[p][contribution]
                triple_counts[p] += surviving.bit_count()

                if p == EXACT_TAIL_PRIME:
                    tail = surviving
                    while tail:
                        lowest_bit = tail & -tail
                        index = lowest_bit.bit_length() - 1
                        tail -= lowest_bit
                        terminal_candidates += 1
                        x = x_values[index]
                        remainder = TARGET - comb(x, 4) - high_atom
                        if remainder < 1:
                            height_rejections += 1
                            continue
                        exact_tail_checks += 1
                        discriminant = 8 * remainder + 1
                        root = isqrt(discriminant)
                        if root >= 3 and root % 2 == 1 and root * root == discriminant:
                            exact_tail_solutions += 1

                if not surviving:
                    break
            final_survivors += surviving.bit_count()

    require(triple_counts == EXPECTED_TRIPLE_COUNTS, "frozen sieve fingerprint")
    require(terminal_candidates == 324, "terminal residue-sieve tail size")
    require(height_rejections == 31, "terminal negative-height count")
    require(exact_tail_checks == 293, "exact integer-square test count")
    require(exact_tail_solutions == 0, "candidate has no representation")
    require(final_survivors == 0, "pure modular terminal route")
    return (
        triple_counts,
        universe,
        terminal_candidates,
        height_rejections,
        exact_tail_checks,
        final_survivors,
    )


def positive_control() -> tuple[int, int, int, int, int]:
    """A planted representation checks the square-discriminant consequence."""

    w, x, y, z = 42_321, 511, 71, 19
    value = comb(w, 2) + comb(x, 4) + comb(y, 6) + comb(z, 8)
    remainder = value - comb(x, 4) - comb(y, 6) - comb(z, 8)
    root = isqrt(8 * remainder + 1)
    require(root == 2 * w - 1, "planted discriminant root")
    require((root + 1) // 2 == w, "planted w recovery")
    return value, w, x, y, z


period_rows = audit_period_controls()
regular_traces = audit_exceptional_regular_rows()
positive = positive_control()

small_local_levels = {2: 4, 3: 2, 5: 2, 7: 2}
small_local = {
    p: local_factor_small_prime(p, level)
    for p, level in small_local_levels.items()
}
small_local_moduli = {2: 16, 3: 9, 5: 25, 7: 49}
expected_small = {
    2: Fraction(1, 1),
    3: Fraction(68, 81),
    5: Fraction(566, 625),
    7: Fraction(310, 343),
}
require(small_local == expected_small, "small-prime local density factors")

large_local = {
    p: local_factor_large_prime(p) for p in (11, 13, 17, 19, 23, 31, 499)
}
expected_large = {
    11: Fraction(72, 121),
    13: Fraction(154, 169),
    17: Fraction(240, 289),
    19: Fraction(316, 361),
    23: Fraction(472, 529),
    31: Fraction(942, 961),
    499: Fraction(248_998, 249_001),
}
require(large_local == expected_large, "large-prime first-level density factors")

critical_rows = {
    p: critical_tuple_audit(p) for p in (11, 13, 17, 19, 23, 31, 499)
}
expected_critical_rows = {
    11: (0, 0, Fraction(72, 121)),
    13: (0, 0, Fraction(154, 169)),
    17: (0, 0, Fraction(240, 289)),
    19: (0, 0, Fraction(316, 361)),
    23: (0, 0, Fraction(472, 529)),
    31: (4, 0, Fraction(29_198, 29_791)),
    499: (2, 0, Fraction(124_250_000, 124_251_499)),
}
require(critical_rows == expected_critical_rows, "critical p-adic lift rows")

# The mod-33 count uses the true periods 33,99,99,99.  It is deliberately
# computed directly as a separate local fingerprint.
modulus_33 = 33
distribution_33: Counter[int] = Counter({0: 1})
denominator_33 = 1
for k in DEGREES:
    period = 1
    for p, a in ((3, 1), (11, 1)):
        period *= prime_power_period(k, p, a)
    atom = Counter(comb(t, k) % modulus_33 for t in range(period))
    next_distribution: Counter[int] = Counter()
    for left, left_count in distribution_33.items():
        for right, right_count in atom.items():
            next_distribution[(left + right) % modulus_33] += left_count * right_count
    distribution_33 = next_distribution
    denominator_33 *= period
probability_33 = Fraction(distribution_33[TARGET % 33], denominator_33)
require(TARGET % 33 == 20, "hostile mod-33 target class")
require(probability_33 == Fraction(16, 1089), "mod-33 target density")
minimum_classes_33 = [
    residue
    for residue, count in sorted(distribution_33.items())
    if count == min(distribution_33.values())
]
require(
    minimum_classes_33 == [20],
    "mod-33 minimum-density class is not uniquely 20",
)

(
    triple_counts,
    universe,
    terminal_candidates,
    height_rejections,
    exact_tail_checks,
    final_survivors,
) = audit_counterexample()

print("THM-4026 SUN 2-4-6-8 COUNTEREXAMPLE INDEPENDENT EXACT AUDIT")
print(f"target = {TARGET}")
print(f"canonical (x,y,z) universe = {universe}")
print("canonical maxima = x:12112 y:932 z:281")
print("quadratic-residue sieve triple counts:")
for p in SIEVE_PRIMES:
    print(f"  p={p:3d} survivors={triple_counts[p]}")
print(f"terminal candidates after p<=89 = {terminal_candidates}")
print(f"negative-height rejections after p<=89 = {height_rejections}")
print(f"exact isqrt checks after p<=89 = {exact_tail_checks}")
print("exact isqrt solutions = 0")
print(f"pure modular survivors after p<=137 = {final_survivors}")
print(f"positive control (value,w,x,y,z) = {positive}")
print(f"finite exact-period hostile rows = {period_rows}")
print("exceptional-prime regular-coverage traces:")
for p in sorted(regular_traces):
    print(f"  p={p:2d} partial_sizes={regular_traces[p]}")
print("small-prime exact local factors at audited levels:")
for p in sorted(small_local):
    print(
        f"  p={p} q={small_local_moduli[p]} "
        f"level={p}^{small_local_levels[p]} sigma={small_local[p]}"
    )
print("selected p>8 first-level local factors:")
for p in sorted(large_local):
    print(f"  p={p} sigma_1={large_local[p]}")
print("fully critical p-adic lift rows (critical,lifting,stable_sigma):")
for p in sorted(critical_rows):
    print(f"  p={p} row={critical_rows[p]}")
print(f"target mod 33 = {TARGET % 33}")
print(f"target probability mod 33 = {probability_33}")
print(f"minimum-density classes mod 33 = {minimum_classes_33}")
print("PERIOD-CONTROLS: PASS")
print("REGULAR-LOCAL-COVERAGE-CONTROLS: PASS")
print("EXACT-SQUARE-TAIL: PASS")
print("PURE-MODULAR-TERMINAL: PASS")
print("RESULT: COUNTEREXAMPLE VERIFIED")
