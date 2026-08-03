#!/usr/bin/env python3
"""Exact controls for the prime right-boundary rigidity theorem.

The accompanying proof shows that, at every prime modulus, the cells just to
the right of alpha=a/p already force a micro-staircase blocker to be a scalar
ramp.  The key interpolation object is

    A(y) = product_{l in supp(w)} (l-y*w_l).

This program checks the boundary convention, exhausts the normalized systems
for p=3,5,7, verifies a constructive witness extracted from the proof, audits
the older critical-support moment argument, and records the composite n=14
hostile.  It uses explicit truth gates rather than Python asserts.

Scope: exact controls for an elementary proof.  This is a theorem about the
finite residue/cell model, not a speed-to-residue lift and not LRC.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path


SMALL_PRIMES = (3, 5, 7)
CONVENTION_PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
MOMENT_PRIMES = (7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def floor_bin(p: int, index: int, alpha: Fraction) -> int:
    fractional = (index * alpha) % 1
    scaled = p * fractional
    return scaled.numerator // scaled.denominator


def right_boundary_pattern(p: int, a: int) -> tuple[int, ...]:
    epsilon = Fraction(1, 2 * p * p * (p - 1))
    alpha = Fraction(a, p) + epsilon
    return tuple(floor_bin(p, index, alpha) for index in range(1, p))


def blocks_candidate(
    p: int, vector: tuple[int, ...], a: int, shift: int
) -> bool:
    return any(
        (a * index + shift * vector[index - 1]) % p in (0, p - 1)
        for index in range(1, p)
    )


def boundary_blocker(p: int, vector: tuple[int, ...]) -> bool:
    return all(
        blocks_candidate(p, vector, a, shift)
        for a in range(p)
        for shift in range(p)
    )


def normalize(p: int, vector: tuple[int, ...]) -> tuple[int, ...]:
    multiplier = vector[0]
    return tuple(
        (value - multiplier * index) % p
        for index, value in enumerate(vector, start=1)
    )


def proof_witness(
    p: int, vector: tuple[int, ...]
) -> tuple[int, int, int, int]:
    """Return (a,s,y,j) exposed by the interpolation/sum proof.

    The input must be normalized and nonzero.  First find y outside
    X={l/w_l} for which A(y) differs from L=product(S).  If every such value
    equalled L, summing A over F_p would give both zero and -|X|L.  For this y,
    some j in S fails the two-incidence cover, yielding the boundary miss.
    """
    require(vector[0] == 0, (p, "not normalized", vector))
    support = tuple(
        index for index in range(1, p) if vector[index - 1] % p != 0
    )
    require(support, (p, "zero vector has no proof witness"))
    require(len(support) <= p - 2, (p, support))

    q = {index: pow(vector[index - 1], -1, p) for index in support}
    x_set = {index * q[index] % p for index in support}
    support_product = 1
    for index in support:
        support_product = support_product * index % p

    chosen_y = None
    for y in range(p):
        if y in x_set:
            continue
        moving_product = 1
        for index in support:
            moving_product = moving_product * (
                index - y * vector[index - 1]
            ) % p
        if moving_product != support_product:
            chosen_y = y
            break
    require(chosen_y is not None, (p, "sum contradiction failed", vector))

    chosen_j = None
    for j in support:
        incidences = {
            index * q[index] % p for index in support
        } | {
            (index - j) * q[index] % p for index in support
        }
        if chosen_y not in incidences:
            chosen_j = j
            break
    require(chosen_j is not None, (p, "interpolation witness failed", vector))

    inverse_j = pow(chosen_j, -1, p)
    a = -inverse_j % p
    shift = chosen_y * inverse_j % p
    require(
        not blocks_candidate(p, vector, a, shift),
        (p, "translated witness blocks", vector, a, shift, chosen_y, chosen_j),
    )
    return a, shift, chosen_y, chosen_j


def audit_boundary_convention() -> int:
    cases = 0
    for p in CONVENTION_PRIMES:
        for a in range(p):
            exact = right_boundary_pattern(p, a)
            expected = tuple(a * index % p for index in range(1, p))
            require(exact == expected, (p, a, exact, expected))
            cases += p - 1
    return cases


def exhaustive_small_primes() -> tuple[tuple[int, int, int, int], ...]:
    rows = []
    for p in SMALL_PRIMES:
        scanned = 0
        blockers = []
        witness_checksum = 0
        for tail in product(range(p), repeat=p - 2):
            vector = (0,) + tail
            scanned += 1
            if boundary_blocker(p, vector):
                blockers.append(vector)
            if any(vector):
                a, shift, y, j = proof_witness(p, vector)
                witness_checksum = (
                    witness_checksum + a + 3 * shift + 5 * y + 7 * j
                ) % 1_000_003
        require(blockers == [(0,) * (p - 1)], (p, blockers[:3]))
        rows.append((p, scanned, len(blockers), witness_checksum))
    return tuple(rows)


def deterministic_dense_controls() -> tuple[tuple[int, int, int, int, int], ...]:
    rows = []
    for p in (11, 13, 17, 19):
        vectors = []
        for degree in range(1, 6):
            for offset in range(1, min(p, 8)):
                vector = [0]
                for index in range(2, p):
                    vector.append((pow(index, degree, p) + offset * index + 1) % p)
                if any(vector):
                    vectors.append(tuple(vector))
        for mask in range(1, 1 << min(p - 2, 10)):
            vector = [0] * (p - 1)
            for bit in range(min(p - 2, 10)):
                if (mask >> bit) & 1:
                    index = bit + 2
                    vector[index - 1] = (3 * index + 1) % p or 1
            vectors.append(tuple(vector))

        checksum = 0
        for vector in vectors:
            require(not boundary_blocker(p, vector), (p, "dense hostile survived"))
            a, shift, y, j = proof_witness(p, vector)
            checksum = (checksum + a + 3 * shift + 5 * y + 7 * j) % 1_000_003
        rows.append((p, len(vectors), checksum, min(map(sum_support, vectors)),
                     max(map(sum_support, vectors))))
    return tuple(rows)


def sum_support(vector: tuple[int, ...]) -> int:
    return sum(value != 0 for value in vector)


def audit_scalar_gauge() -> int:
    cases = 0
    for p in SMALL_PRIMES:
        for vector in product(range(p), repeat=p - 1):
            normalized = normalize(p, vector)
            for a in range(p):
                for shift in range(p):
                    reindexed_a = (a - shift * vector[0]) % p
                    require(
                        blocks_candidate(p, normalized, a, shift)
                        == blocks_candidate(p, vector, reindexed_a, shift),
                        (p, vector, a, shift),
                    )
                    cases += 1
    return cases


def audit_critical_moment_argument() -> tuple[int, int, int]:
    root_cases = 0
    zero_sum_cases = 0
    maximum_root_multiplicity = 0
    for p in MOMENT_PRIMES:
        h = (p + 1) // 2
        roots = []
        total_multiplicity = 0
        for value in range(1, p):
            polynomial = (
                pow(value, h, p) - pow(value, h - 1, p) + h - 1
            ) % p
            if polynomial != 0:
                continue
            roots.append(value)
            require(value in (h % p, (2 - h) % p), (p, h, value))
            derivative = (
                h * pow(value, h - 1, p)
                - (h - 1) * pow(value, h - 2, p)
            ) % p
            multiplicity = 1 if derivative else 2
            if multiplicity == 2:
                second_derivative = (
                    h * (h - 1) * pow(value, h - 2, p)
                    - (h - 1) * (h - 2) * pow(value, h - 3, p)
                ) % p
                require(second_derivative != 0, (p, "triple root", value))
            total_multiplicity += multiplicity
            maximum_root_multiplicity = max(maximum_root_multiplicity, multiplicity)
        require(len(roots) <= 2, (p, roots))
        require(total_multiplicity <= 3 < h, (p, h, roots, total_multiplicity))
        root_cases += 1

        for constant in range(1, p):
            roots_zero_sum = tuple(
                value for value in range(1, p)
                if pow(value, h, p) == constant
            )
            require(len(roots_zero_sum) <= 2 < h, (p, constant, roots_zero_sum))
            zero_sum_cases += 1
    return root_cases, zero_sum_cases, maximum_root_multiplicity


def audit_composite_hostile() -> tuple[int, int]:
    n = 14
    vector = (0, 1) + (0,) * 11
    covered = sum(
        blocks_candidate(n, vector, a, shift)
        for a in range(n)
        for shift in range(n)
    )
    require(covered == n * n, covered)
    return covered, n * n


def main() -> None:
    source_bytes = Path(__file__).read_bytes()
    source_lf = source_bytes.replace(b"\r\n", b"\n")
    require(b"\r" not in source_lf, "source contains a bare carriage return")

    convention_cases = audit_boundary_convention()
    exhaustive_rows = exhaustive_small_primes()
    dense_rows = deterministic_dense_controls()
    gauge_cases = audit_scalar_gauge()
    moment_rows = audit_critical_moment_argument()
    composite = audit_composite_hostile()

    print("Prime right-boundary micro-staircase rigidity audit")
    print("scope=PROVED_elementary_theorem_plus_FINITE-EXACT_controls;not_LRC")
    print(f"right_boundary_pattern_cases={convention_cases};status=PASS")
    print(f"normalized_exhaustive_rows={exhaustive_rows}")
    print(f"constructive_dense_rows={dense_rows}")
    print(f"scalar_gauge_candidate_equivalences={gauge_cases};status=PASS")
    print(
        "critical_support_moment_audit="
        f"primes={moment_rows[0]},zero_sum_constants={moment_rows[1]},"
        f"maximum_root_multiplicity={moment_rows[2]};status=PASS"
    )
    print(
        f"composite_hostile=n14_(0,1,0^11);boundary_covered={composite[0]}/"
        f"{composite[1]};prime_field_step_is_essential"
    )
    print(f"source_lf_sha256={sha256(source_lf).hexdigest()}")
    print("ALL_CHECKS_PASSED")


if __name__ == "__main__":
    main()
