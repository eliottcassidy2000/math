#!/usr/bin/env python3
"""Exact companion for THM-2403.

The script exhausts the finite clean-cover endpoint lemma, including every
three/four-root residual core, every retained two-root word, and every
ordinary blocker gate with zero, one, or two failed shifts.  It also checks
the sharp gap-one word, the exact circle gate which realizes its two killed
shifts, the prime-cyclotomic reduction, the variance bounds, and the target
phase projector.
"""

from fractions import Fraction
from itertools import combinations


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def translated(word: frozenset[int], shift: int) -> frozenset[int]:
    return frozenset((entry + shift) % P for entry in word)


def clean_cover_examples() -> dict[tuple[int, int], tuple[frozenset[int], ...]]:
    # In every example the order is q_*, q_i, q_2, q_3, q_4, guard.
    examples = {
        # The double is q_* cap q_i: |A_0|=1 and |B|=3.
        (1, 3): (
            frozenset((0, 1)),
            frozenset((1, 2)),
            frozenset((3, 4)),
            frozenset((5, 6)),
            frozenset((7, 8)),
            frozenset((9, 10, 11, 12)),
        ),
        # The double is q_i cap q_2: |A_0|=2 and |B|=3.
        (2, 3): (
            frozenset((0, 1)),
            frozenset((2, 3)),
            frozenset((2, 4)),
            frozenset((5, 6)),
            frozenset((7, 8)),
            frozenset((9, 10, 11, 12)),
        ),
        # The double is wholly retained: |A_0|=2 and |B|=4.
        (2, 4): (
            frozenset((0, 1)),
            frozenset((2, 3)),
            frozenset((4, 5)),
            frozenset((4, 6)),
            frozenset((7, 8)),
            frozenset((9, 10, 11, 12)),
        ),
    }
    for key, words in examples.items():
        counts = [sum(root in word for word in words) for root in range(P)]
        require(sorted(counts) == [1] * 12 + [2], "bad clean-cover profile")
        q_star, q_i, *retained = words
        core = frozenset(range(P)) - frozenset().union(*retained)
        require(
            (len(core - q_i), len(core)) == key,
            "bad clean-cover residual type",
        )
        require(len(q_star) == len(q_i) == 2, "bad ordinary word size")
        require(len(words[-1]) == 4, "bad guard word size")
    return examples


def exhaustive_endpoint_lemma():
    minimum: dict[tuple[int, int, int], int] = {}
    configurations = 0
    gated_cases = 0
    unique_mass_vectors: set[tuple[int, ...]] = set()
    type_counts: dict[tuple[int, int], int] = {}

    nonzero_shifts = tuple(range(1, P))
    zero_banks = tuple(
        frozenset(bank)
        for size in range(3)
        for bank in combinations(nonzero_shifts, size)
    )
    require(len(zero_banks) == 79, "wrong blocker-gate bank")

    for core_size in (3, 4):
        for core_tuple in combinations(range(P), core_size):
            core = frozenset(core_tuple)
            for word_tuple in combinations(range(P), 2):
                word = frozenset(word_tuple)
                a0 = len(core - word)
                residual_type = (a0, core_size)
                if residual_type not in ((1, 3), (2, 3), (2, 4)):
                    continue
                configurations += 1
                type_counts[residual_type] = type_counts.get(residual_type, 0) + 1

                ungated = tuple(
                    len(core - translated(word, shift)) for shift in range(P)
                )
                require(sum(ungated) == 11 * core_size, "eleven-cover identity")
                require(ungated[0] == a0, "wrong base deletion support")

                for zero_bank in zero_banks:
                    masses = tuple(
                        0 if shift in zero_bank else ungated[shift]
                        for shift in range(P)
                    )
                    gap = sum(masses) - P * masses[0]
                    lower = 9 * core_size - P * a0
                    require(gap >= lower > 0, "target-axis gap failed")

                    key = (a0, core_size, len(zero_bank))
                    minimum[key] = min(minimum.get(key, gap), gap)
                    unique_mass_vectors.add(masses)
                    gated_cases += 1

    expected_minimum = {
        (1, 3, 0): 20,
        (1, 3, 1): 17,
        (1, 3, 2): 14,
        (2, 3, 0): 7,
        (2, 3, 1): 4,
        (2, 3, 2): 1,
        (2, 4, 0): 18,
        (2, 4, 1): 14,
        (2, 4, 2): 10,
    }
    require(minimum == expected_minimum, "wrong sharp gap table")

    # For degree at most twelve, reduction modulo Phi_13 subtracts the
    # X^12 coefficient from every lower coefficient.  Thus a rational
    # target mode can vanish only for a flat vector.
    for masses in unique_mass_vectors:
        reduced = tuple(masses[j] - masses[12] for j in range(12))
        require(any(reduced), "nonconstant vector reduced to zero")

        total = sum(masses)
        gap = total - P * masses[0]
        variance = sum(
            (P * entry - total) ** 2 for entry in masses
        ) / Fraction(P**3)
        require(
            variance >= Fraction(gap * gap, 2028),
            "variance floor failed",
        )

    return (
        configurations,
        gated_cases,
        len(unique_mass_vectors),
        type_counts,
        minimum,
    )


def sharp_gap_one_control():
    core = frozenset((0, 1, 3))
    word = frozenset((2, 3))
    killed = frozenset((2, 3))
    masses = tuple(
        0
        if shift in killed
        else len(core - translated(word, shift))
        for shift in range(P)
    )
    require(
        masses == (2, 2, 0, 0, 3, 3, 3, 3, 3, 3, 2, 1, 2),
        "wrong sharp mass vector",
    )
    require(sum(masses) - P * masses[0] == 1, "gap one is not sharp")

    # The ordinary strict danger arcs centred at s/13 kill exactly shifts
    # two and three at z=5/26, while the base shift zero is safe.
    z = Fraction(5, 26)
    danger_shifts = tuple(
        shift
        for shift in range(P)
        if circle_norm(z - Fraction(shift, P)) < Fraction(1, 14)
    )
    require(danger_shifts == (2, 3), "wrong physical two-shift gate")
    return masses, danger_shifts


def target_phase_projector():
    # With eta=e_a-e_i, the lawful action shifts the a-factor by -s/13
    # and the i-factor by +s/13.  Multiplying by zeta^(b s) selects
    # n_a-n_i=b modulo thirteen.
    checks = 0
    live = 0
    for n_a in range(P):
        for n_i in range(P):
            for b in range(P):
                residues = tuple(
                    (b + n_i - n_a) * shift % P for shift in range(P)
                )
                selected = all(residue == 0 for residue in residues)
                require(
                    selected == ((n_a - n_i - b) % P == 0),
                    "target phase projector failed",
                )
                live += int(selected)
                checks += 1
    require(live == P * P, "wrong number of selected target phases")
    return checks, live


def exact_floors():
    universal_mass = Fraction(1, 26754)
    common_core_mass = Fraction(66, 4459)
    floors = {
        "universal_fibre_energy": universal_mass**2 / 2028,
        "universal_fibre_max": universal_mass / 156,
        "universal_physical_energy": universal_mass**2 / (2028 * 13**2),
        "universal_physical_max": universal_mass / (156 * 13),
        "common_core_fibre_energy": common_core_mass**2 / 2028,
        "common_core_fibre_max": common_core_mass / 156,
        "common_core_physical_energy": common_core_mass**2 / (2028 * 13**2),
        "common_core_physical_max": common_core_mass / (156 * 13),
    }
    require(
        floors["universal_fibre_energy"] == Fraction(1, 1451594774448),
        "wrong universal energy floor",
    )
    require(
        floors["universal_fibre_max"] == Fraction(1, 4173624),
        "wrong universal max floor",
    )
    require(
        floors["universal_physical_energy"] == Fraction(1, 245319516881712),
        "wrong physical universal energy floor",
    )
    require(
        floors["universal_physical_max"] == Fraction(1, 54257112),
        "wrong physical universal max floor",
    )
    require(
        floors["common_core_fibre_energy"] == Fraction(363, 3360173089),
        "wrong common-core energy floor",
    )
    require(
        floors["common_core_fibre_max"] == Fraction(11, 115934),
        "wrong common-core max floor",
    )
    require(
        floors["common_core_physical_energy"] == Fraction(363, 567869252041),
        "wrong physical common-core energy floor",
    )
    require(
        floors["common_core_physical_max"] == Fraction(11, 1507142),
        "wrong physical common-core max floor",
    )
    return floors


def main() -> None:
    examples = clean_cover_examples()
    (
        configurations,
        gated_cases,
        unique_vectors,
        type_counts,
        minimum,
    ) = exhaustive_endpoint_lemma()
    sharp_masses, danger_shifts = sharp_gap_one_control()
    projector_checks, projector_live = target_phase_projector()
    floors = exact_floors()

    print("THM-2403 clean target-axis imbalance exact companion")
    print(
        f"clean_cover_types={len(examples)} "
        + ",".join(f"{a}:{b}" for a, b in sorted(examples))
    )
    print(
        f"endpoint_configurations={configurations} "
        f"gated_cases={gated_cases} unique_mass_vectors={unique_vectors}"
    )
    print(
        "type_counts="
        + ",".join(
            f"{a}:{b}:{type_counts[a, b]}" for a, b in sorted(type_counts)
        )
    )
    print(
        "sharp_gap_table="
        + ",".join(
            f"{a}:{b}:z{z}:{minimum[a,b,z]}"
            for a, b, z in sorted(minimum)
        )
    )
    print(
        "gap_one_masses="
        + ",".join(map(str, sharp_masses))
        + " danger_shifts="
        + ",".join(map(str, danger_shifts))
    )
    print(
        f"target_phase_checks={projector_checks} "
        f"selected={projector_live}"
    )
    print(
        "floors="
        + ",".join(f"{key}:{value}" for key, value in floors.items())
    )
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
