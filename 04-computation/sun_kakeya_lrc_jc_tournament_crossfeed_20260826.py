#!/usr/bin/env python3
"""Exact hostile probes for the Sun/Kakeya -> LRC14/planar-JC crossfeed.

This is a maintained reflection companion, not a canonical theorem
dependency.  It checks:

1. compatible p-adic solution branches for the THM-4026 target;
2. the unique projective collision of THM-4035's two AP tail moments;
3. independence of cyclotomic-CM stabilizers and finite Fourier full-spark;
4. the natural bad-prime critical fibre of the M=11 and M=13 JC faces.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

from sturmian_fibonacci_triangular_kakeya_sixty_clock_thm4035 import (  # noqa: E402
    aggregate_ap_laws,
)


N = 896_315_812_331_399
ROLES = (2, 4, 6, 8)
LOWERS = (2, 3, 5, 7)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def binomial_sum(indices: tuple[int, int, int, int]) -> int:
    return sum(comb(index, degree) for index, degree in zip(indices, ROLES))


def floor_log(prime: int, degree: int) -> int:
    exponent = 0
    power = prime
    while power <= degree:
        exponent += 1
        power *= prime
    return exponent


def exact_period(prime: int, exponent: int, degree: int) -> int:
    return prime ** (exponent + floor_log(prime, degree))


def regular_base_solution(prime: int) -> tuple[int, int, int, int]:
    """First solution mod p whose triangular derivative is a unit."""
    for indices in product(range(prime), repeat=4):
        if (2 * indices[0] - 1) % prime == 0:
            continue
        if (binomial_sum(indices) - N) % prime == 0:
            return indices
    raise RuntimeError((prime, "no regular base solution"))


def local_branch(prime: int, maximum_exponent: int = 3):
    """A compatible solution branch through p^maximum_exponent.

    For p=2 use the even triangular permutation.  For p=3,5,7 use the
    degree-(p+1) permutation from THM-4027.  For p>8 lift only w along a
    regular triangular derivative.
    """
    branch = []
    previous = None
    if prime == 2:
        for exponent in range(1, maximum_exponent + 1):
            modulus = prime**exponent
            roots = [t for t in range(modulus) if comb(2 * t, 2) % modulus == N % modulus]
            require(len(roots) == 1, (prime, exponent, roots))
            indices = (2 * roots[0], 3, 5, 7)
            if previous is not None:
                require(roots[0] % (modulus // prime) == previous, "2-adic incompatibility")
            previous = roots[0]
            require((binomial_sum(indices) - N) % modulus == 0, (prime, exponent))
            branch.append(indices)
        return tuple(branch)

    if prime in (3, 5, 7):
        role = {3: 1, 5: 2, 7: 3}[prime]
        degree = prime + 1
        for exponent in range(1, maximum_exponent + 1):
            modulus = prime**exponent
            roots = [
                q
                for q in range(modulus)
                if comb(prime * q + 1, degree) % modulus == (N - 1) % modulus
            ]
            require(len(roots) == 1, (prime, exponent, len(roots)))
            if previous is not None:
                require(roots[0] % (modulus // prime) == previous, (prime, "incompatibility"))
            previous = roots[0]
            mutable = list(LOWERS)
            mutable[role] = prime * roots[0] + 1
            indices = tuple(mutable)
            require((binomial_sum(indices) - N) % modulus == 0, (prime, exponent))
            branch.append(indices)
        return tuple(branch)

    indices = list(regular_base_solution(prime))
    branch.append(tuple(indices))
    require((binomial_sum(tuple(indices)) - N) % prime == 0, (prime, 1))
    for exponent in range(1, maximum_exponent):
        old_modulus = prime**exponent
        new_modulus = old_modulus * prime
        candidates = []
        for digit in range(prime):
            candidate = list(indices)
            candidate[0] += digit * old_modulus
            if (binomial_sum(tuple(candidate)) - N) % new_modulus == 0:
                candidates.append(candidate)
        require(len(candidates) == 1, (prime, exponent + 1, len(candidates)))
        indices = candidates[0]
        branch.append(tuple(indices))
    return tuple(branch)


def crt_pair(residue_a: int, modulus_a: int, residue_b: int, modulus_b: int):
    require(__import__("math").gcd(modulus_a, modulus_b) == 1, "CRT moduli not coprime")
    step = ((residue_b - residue_a) * pow(modulus_a, -1, modulus_b)) % modulus_b
    value = (residue_a + modulus_a * step) % (modulus_a * modulus_b)
    return value, modulus_a * modulus_b


def sun_profinite_probe():
    primes = (2, 3, 5, 7, 11, 13)
    exponent = 3
    branches = {prime: local_branch(prime, exponent) for prime in primes}
    for prime, branch in branches.items():
        for level, indices in enumerate(branch, 1):
            require((binomial_sum(indices) - N) % (prime**level) == 0, (prime, level))
            if level > 1:
                for coordinate, degree in enumerate(ROLES):
                    old_period = exact_period(prime, level - 1, degree)
                    require(
                        indices[coordinate] % old_period
                        == branch[level - 2][coordinate] % old_period,
                        (prime, level, degree, "index-period incompatibility"),
                    )

    global_indices = []
    coordinate_periods = []
    for coordinate, degree in enumerate(ROLES):
        residue = 0
        modulus = 1
        for prime in primes:
            period = exact_period(prime, exponent, degree)
            local_residue = branches[prime][-1][coordinate] % period
            residue, modulus = crt_pair(residue, modulus, local_residue, period)
        if residue < LOWERS[coordinate]:
            residue += ((LOWERS[coordinate] - residue + modulus - 1) // modulus) * modulus
        global_indices.append(residue)
        coordinate_periods.append(modulus)

    global_indices_tuple = tuple(global_indices)
    target_modulus = 1
    for prime in primes:
        target_modulus *= prime**exponent
    difference = binomial_sum(global_indices_tuple) - N
    require(difference != 0 and difference % target_modulus == 0, "CRT hostile failed")

    branch_digest = sha256(repr(branches).encode()).hexdigest()
    print(
        "sun_profinite_control=(primes:2,3,5,7,11,13,levels:3,"
        f"branch_sha256:{branch_digest})"
    )
    print(
        "sun_crt_archimedean_hostile=("
        f"target_modulus:{target_modulus},coordinate_periods:{tuple(coordinate_periods)},"
        f"indices:{global_indices_tuple},difference_quotient_bits:{abs(difference // target_modulus).bit_length()})"
    )
    return branches, global_indices_tuple


def ap_projective_probe():
    _, _, signatures, moments = aggregate_ap_laws()
    fibres = defaultdict(list)
    for phase, (first, second) in enumerate(moments):
        require(first != 0, (phase, "zero first moment"))
        fibres[second / first].append(phase)
    collisions = tuple(tuple(phases) for phases in fibres.values() if len(phases) > 1)
    require(len(fibres) == 59 and collisions == ((18, 25),), (len(fibres), collisions))
    scale_first = moments[25][0] / moments[18][0]
    scale_second = moments[25][1] / moments[18][1]
    require(scale_first == scale_second == Fraction(25, 21), "projective scale mismatch")
    require(signatures[18] != signatures[25], "full AP laws unexpectedly agree")
    print(
        "ap_two_moment_projectivization=(affine_images:60,projective_images:59,"
        "unique_collision:(18,25),scale:25/21,full_laws_equal:False)"
    )
    return collisions


def multiplicative_order(value: int, prime: int) -> int:
    current = 1
    for exponent in range(1, prime):
        current = current * value % prime
        if current == 1:
            return exponent
    raise RuntimeError((value, prime, "order"))


def root_of_order(order: int, prime: int) -> int:
    return next(value for value in range(2, prime) if multiplicative_order(value, prime) == order)


def determinant_mod(matrix, prime: int) -> int:
    work = [list(row) for row in matrix]
    determinant = 1
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column] % prime),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            determinant = -determinant
        value = work[column][column] % prime
        determinant = determinant * value % prime
        inverse = pow(value, -1, prime)
        for row in range(column + 1, len(work)):
            multiplier = work[row][column] * inverse % prime
            for index in range(column, len(work)):
                work[row][index] = (
                    work[row][index] - multiplier * work[column][index]
                ) % prime
    return determinant % prime


def cm_stabilizer(order: int, cm_type: tuple[int, ...]):
    packet = set(cm_type)
    return tuple(
        unit
        for unit in range(1, order)
        if {(unit * character) % order for character in packet} == packet
    )


def zero_fourier_minors(order: int, cm_type: tuple[int, ...], prime: int):
    root = root_of_order(order, prime)
    zeros = []
    for rows in combinations(range(order), len(cm_type)):
        matrix = [[pow(root, row * character, prime) for character in cm_type] for row in rows]
        if determinant_mod(matrix, prime) == 0:
            zeros.append(rows)
    return root, tuple(zeros)


def cm_fourier_probe():
    packets = {
        7: ((1, 2, 4), (29,)),
        11: ((4, 5, 8, 9, 10), (23, 67, 89, 199)),
        13: ((5, 6, 9, 10, 11, 12), (53, 79, 131, 157)),
    }
    expected = {
        (7, 29): 0,
        (11, 23): 22,
        (11, 67): 11,
        (11, 89): 0,
        (11, 199): 0,
        (13, 53): 26,
        (13, 79): 26,
        (13, 131): 0,
        (13, 157): 13,
    }
    records = []
    for order, (cm_type, primes) in packets.items():
        stabilizer = cm_stabilizer(order, cm_type)
        for prime in primes:
            root, zeros = zero_fourier_minors(order, cm_type, prime)
            require(len(zeros) == expected[(order, prime)], (order, prime, len(zeros)))
            records.append(
                (order, cm_type, stabilizer, prime, root, len(zeros), zeros[:1])
            )
    require(cm_stabilizer(7, packets[7][0]) == (1, 2, 4), "order-seven hostile")
    require(cm_stabilizer(11, packets[11][0]) == (1,), "M11 primitivity")
    require(cm_stabilizer(13, packets[13][0]) == (1,), "M13 primitivity")
    print(f"cm_fourier_independence={tuple(records)}")
    return tuple(records)


def face_value(m: int, coefficient_a: int, coefficient_b: int, s: int, pvar: int) -> int:
    return 1 - coefficient_a * s * pvar**m - coefficient_b * s**3 * pvar ** (m - 1)


def face_derivatives_mod(
    m: int, coefficient_a: int, coefficient_b: int, s: int, pvar: int, prime: int
):
    derivative_s = (
        -coefficient_a * pow(pvar, m, prime)
        - 3 * coefficient_b * pow(s, 2, prime) * pow(pvar, m - 1, prime)
    ) % prime
    derivative_p = (
        -m * coefficient_a * s * pow(pvar, m - 1, prime)
        - (m - 1) * coefficient_b * pow(s, 3, prime) * pow(pvar, m - 2, prime)
    ) % prime
    return derivative_s, derivative_p


def face_bad_prime_probe():
    """Use least nonnegative integer lifts of every residue-class label.

    The first residual depends on this integral-lift convention.  These rows
    are bad-reduction diagnostics, not characteristic-zero JC obstructions.
    """
    records = []
    for m, prime, expected_point, expected_residual in (
        (5, 11, (4, 7), 5),
        (6, 13, (8, 3), 4),
    ):
        residual_histogram = defaultdict(int)
        collision_lifts = 0
        collision_nonlifts = 0
        all_pairs = 0
        unit_pair_unique = True
        chosen_record = None
        for coefficient_a in range(1, prime):
            for coefficient_b in range(1, prime):
                all_pairs += 1
                singular = []
                for s in range(1, prime):
                    for pvar in range(1, prime):
                        value = face_value(m, coefficient_a, coefficient_b, s, pvar) % prime
                        ds, dp = face_derivatives_mod(
                            m, coefficient_a, coefficient_b, s, pvar, prime
                        )
                        if value == ds == dp == 0:
                            singular.append((s, pvar))
                unit_pair_unique &= len(singular) == 1
                require(len(singular) == 1, (m, prime, coefficient_a, coefficient_b, singular))
                s, pvar = singular[0]
                integer_value = face_value(m, coefficient_a, coefficient_b, s, pvar)
                require(integer_value % prime == 0, "singular representative residual")
                residual = (integer_value // prime) % prime
                residual_histogram[residual] += 1
                lift_count = 0 if residual else prime * prime
                if (coefficient_a + coefficient_b) % prime == 0:
                    if lift_count:
                        collision_lifts += 1
                    else:
                        collision_nonlifts += 1
                if (coefficient_a, coefficient_b) == (1, 1):
                    explicit_lifts = 0
                    for digit_s, digit_p in product(range(prime), repeat=2):
                        if (
                            face_value(
                                m,
                                coefficient_a,
                                coefficient_b,
                                s + prime * digit_s,
                                pvar + prime * digit_p,
                            )
                            % (prime * prime)
                            == 0
                        ):
                            explicit_lifts += 1
                    require((s, pvar) == expected_point, (m, "chosen singular point"))
                    require(residual == expected_residual and explicit_lifts == 0, (m, residual))
                    chosen_record = ((s, pvar), residual, explicit_lifts)
        require(unit_pair_unique, (m, "unit-pair uniqueness"))
        records.append(
            (
                m,
                prime,
                all_pairs,
                tuple(sorted(residual_histogram.items())),
                collision_lifts,
                collision_nonlifts,
                chosen_record,
            )
        )
    print(f"jc_face_natural_bad_prime={tuple(records)}")
    return tuple(records)


def main() -> None:
    print("SUN/KAKEYA CROSSFEED EXACT HOSTILE PROBE")
    print("JC_RESIDUAL_LIFT_CONVENTION=least_nonnegative_integer_representatives")
    semantic = (
        sun_profinite_probe(),
        ap_projective_probe(),
        cm_fourier_probe(),
        face_bad_prime_probe(),
    )
    print(f"semantic_sha256={sha256(repr(semantic).encode()).hexdigest()}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
