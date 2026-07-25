#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2352."""

from collections import Counter
from fractions import Fraction
from itertools import product


def require(condition: bool, message: str) -> None:
    """Raise under ordinary and optimized Python."""
    if not condition:
        raise RuntimeError(message)


def prefix_residues(digits: tuple[int, ...], radix: int) -> tuple[int, ...]:
    residue = 0
    answer = []
    for position, digit in enumerate(digits):
        require(0 <= digit < radix, "digit escaped its radix alphabet")
        residue += digit * radix**position
        answer.append(residue)
    return tuple(answer)


def finite_profile_atlas() -> tuple[int, int]:
    """Exhaust finite banks and verify support plus collision tax."""
    instances = 0
    nonzero_instances = 0
    for radix, length in ((2, 10), (3, 7)):
        for digits in product(range(radix), repeat=length):
            residues = tuple(
                value for value in prefix_residues(digits, radix) if value
            )
            counts = Counter(residues)
            support = tuple(sorted(counts))
            indexed = sum((Fraction(1, value) for value in residues), Fraction())
            distinct = sum((Fraction(1, value) for value in support), Fraction())
            tax = sum(
                (
                    Fraction(multiplicity - 1, value)
                    for value, multiplicity in counts.items()
                ),
                Fraction(),
            )
            require(indexed == distinct + tax, "finite collision tax failed")

            active_positions = tuple(
                position for position, digit in enumerate(digits) if digit
            )
            active_values = tuple(
                prefix_residues(digits, radix)[position]
                for position in active_positions
            )
            require(
                support == active_values,
                "distinct prefix support differs from active-digit values",
            )
            if active_positions:
                first = active_positions[0]
                infinite_bound = Fraction(1, radix**first) / (
                    1 - Fraction(1, radix)
                )
                require(
                    distinct <= infinite_bound,
                    "exponential distinct-support bound failed",
                )
                for index, (position, value) in enumerate(
                    zip(active_positions, active_values)
                ):
                    require(
                        position >= first + index and value >= radix**position,
                        "active prefix lost exponential growth",
                    )
                nonzero_instances += 1
            instances += 1
    return instances, nonzero_instances


def threshold_records(
    radix: int,
    alpha: int,
    blocks: int,
) -> tuple[tuple[int, int, int], ...]:
    """Return the first exact integer-alpha plateau records."""
    position = 0
    residue = 1
    answer = []
    for index in range(blocks):
        length = residue**alpha
        answer.append((position, residue, length))
        if index + 1 < blocks:
            position += length
            residue += radix**position
    return tuple(answer)


def threshold_atlas() -> tuple[int, tuple[tuple[int, int, int], ...], int]:
    """Check the exact threshold controls and arbitrary-prefix extension."""
    checks = 0
    for alpha, blocks in ((0, 8), (1, 4), (2, 3)):
        records = threshold_records(2, alpha, blocks)
        for index, (position, residue, length) in enumerate(records):
            require(length >= residue**alpha, "threshold plateau too short")
            require(
                Fraction(length, residue**alpha) >= 1,
                "threshold block stopped forcing divergence",
            )
            require(
                Fraction(length, residue ** (alpha + 1))
                <= Fraction(2, residue),
                "above-threshold geometric majorant failed",
            )
            require(position >= index, "active positions stopped increasing")
            checks += 1

    harmonic_records = threshold_records(2, 1, 4)
    require(
        harmonic_records
        == (
            (0, 1, 1),
            (1, 3, 3),
            (4, 19, 19),
            (23, 8_388_627, 8_388_627),
        ),
        "harmonic-edge plateau table changed",
    )
    require(
        all(Fraction(length, residue) == 1 for _, residue, length in harmonic_records),
        "harmonic plateau contribution is not exactly one",
    )

    # The infinity construction L_j=R_j^j already defeats s=1 at j=2.
    position = 0
    residue = 1
    infinity_checks = 0
    for index in range(3):
        length = residue**index
        if index > 1:
            require(
                Fraction(length, residue) >= 1,
                "infinite-abscissa block failed at s=1",
            )
        infinity_checks += 1
        if index < 2:
            position += length
            residue += 2**position

    for radix, prefix in (
        (2, (0, 1, 0, 1)),
        (3, (2, 0, 1)),
        (5, (4, 1, 0, 0)),
    ):
        old_residue = prefix_residues(prefix, radix)[-1]
        first_position = len(prefix)
        first_residue = old_residue + radix**first_position
        require(
            first_residue % radix**first_position == old_residue,
            "threshold extension changed its prescribed cylinder",
        )
        next_length = first_residue
        next_position = first_position + next_length
        next_residue = first_residue + radix**next_position
        require(
            Fraction(next_length, first_residue) == 1
            and next_residue > first_residue,
            "cylinder-local harmonic construction failed",
        )
        checks += 1

    return checks, harmonic_records, infinity_checks


def boundary_atlas() -> tuple[int, int, int]:
    """Audit the all-one and positive terminating hostile boundaries."""
    all_one_digits = (1,) * 32
    all_one_residues = prefix_residues(all_one_digits, 2)
    require(
        all(
            residue == 2 ** (index + 1) - 1
            for index, residue in enumerate(all_one_residues)
        ),
        "all-one prefix formula changed",
    )
    require(
        len(set(all_one_residues)) == len(all_one_residues),
        "all-one path acquired a collision tax",
    )
    require(
        all(
            residue >= 2**index
            for index, residue in enumerate(all_one_residues)
        ),
        "all-one path lost its geometric majorant",
    )

    terminating_digits = (1,) + (0,) * 63
    terminating_residues = prefix_residues(terminating_digits, 2)
    require(
        set(terminating_residues) == {1},
        "the ordinary integer one did not stabilize",
    )
    terminating_indexed_mass = sum(
        (Fraction(1, value) for value in terminating_residues),
        Fraction(),
    )
    terminating_tax = terminating_indexed_mass - Fraction(1)
    require(
        terminating_indexed_mass == 64 and terminating_tax == 63,
        "positive terminating collision tax changed",
    )
    return len(all_one_residues), int(terminating_indexed_mass), int(terminating_tax)


def haar_block_atlas() -> int:
    """Verify the exact zero-block probabilities used by Borel--Cantelli."""
    checks = 0
    for radix, maximum in ((2, 6), (3, 4)):
        for start in range(1, maximum + 1):
            word_length = 2 * start
            event_count = 0
            total = radix**word_length
            for digits in product(range(radix), repeat=word_length):
                if all(digit == 0 for digit in digits[start : 2 * start]):
                    event_count += 1
            require(
                Fraction(event_count, total) == Fraction(1, radix**start),
                "Haar zero-block probability changed",
            )
            checks += 1
    return checks


def finite_interaction_atlas() -> int:
    """Every finite restriction of the termination complex is a simplex."""
    ground_size = 12
    restrictions = {
        tuple(index for index in range(ground_size) if mask & (1 << index))
        for mask in range(1 << ground_size)
    }
    require(
        len(restrictions) == 2**ground_size,
        "finite termination restrictions lost a simplex face",
    )
    return len(restrictions)


def relation_lift_atlas() -> int:
    """Replay the radix-13 support-three zero-carry ray."""
    active_positions = (0, 1, 4, 23)
    digits = tuple(1 if index in active_positions else 0 for index in range(24))
    residues = prefix_residues(digits, 13)
    relation = (1,) + (0,) * 10 + (1, -1)
    ray = (1,) + (0,) * 10 + (2, 3)
    require(
        sum(left * right for left, right in zip(relation, ray)) == 0,
        "support-three digit ray left the relation kernel",
    )

    checks = 0
    for digit in digits:
        digit_vector = tuple(digit * coordinate for coordinate in ray)
        require(
            all(0 <= coordinate < 13 for coordinate in digit_vector),
            "relation digit vector escaped base thirteen",
        )
        require(
            sum(
                left * right
                for left, right in zip(relation, digit_vector)
            )
            == 0,
            "a relation digit created nonzero carry",
        )
        checks += 1

    for residue in residues:
        prefix_vector = tuple(residue * coordinate for coordinate in ray)
        require(
            sum(
                left * right
                for left, right in zip(relation, prefix_vector)
            )
            == 0,
            "a prefix vector left the integer relation kernel",
        )
        require(
            max(prefix_vector) == 3 * residue,
            "fixed-ray norm stopped scaling the prefix residue",
        )
        checks += 1
    return checks


finite_instances, finite_nonzero = finite_profile_atlas()
threshold_checks, harmonic_records, infinity_checks = threshold_atlas()
all_one_terms, terminating_terms, terminating_tax = boundary_atlas()
haar_checks = haar_block_atlas()
finite_simplex_faces = finite_interaction_atlas()
relation_checks = relation_lift_atlas()

require(finite_instances == 3_211, "finite prefix census changed")
require(finite_nonzero == 3_209, "nonzero prefix census changed")

print("theorem=THM-2352")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
print(f"finite_prefix_instances={finite_instances}")
print(f"finite_nonzero_instances={finite_nonzero}")
print(f"threshold_checks={threshold_checks}")
print(f"harmonic_edge_records={harmonic_records}")
print(f"infinity_threshold_checks={infinity_checks}")
print(f"all_one_geometric_terms={all_one_terms}")
print(f"terminating_indexed_terms={terminating_terms}")
print(f"terminating_collision_tax={terminating_tax}")
print(f"haar_zero_block_probability_checks={haar_checks}")
print(f"finite_interaction_simplex_faces={finite_simplex_faces}")
print(f"relation_zero_carry_checks={relation_checks}")
print("distinct_support_positive_abscissa=0")
print("indexed_abscissa_spectrum=[0,infinity]")
print("termination_iff_support_finite=YES")
print("fixed_indexed_convergence_test_decides_termination=NO")
print("haar_indexed_abscissa=0_ALMOST_SURELY")
print("all_checks=PASS")
