#!/usr/bin/env python3
"""Dependency-free exact companion for THM-2384.

This companion has two deliberately separate halves.

* The finite-field half exhausts every projective nonzero unit-weight
  tuple in the THM-2309 owner-pivot normal form.  It checks the pure
  primal helper class, its universal failure to annihilate the pivot
  gauge, gauge invariance of that failure, and the true target/helper
  dual dipoles.
* The circle half exhausts the complete rational open-cell refinement
  for both unit and guard failures, every nonzero helper slope modulo
  thirteen, and all compatible anchored phase profiles.  It checks the
  two zero axes, the exact duplicate-probe count, mixed corner, and
  energy floors.

All proof checks use ``require`` rather than ``assert`` and therefore
remain active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import sys
from typing import Iterable, Sequence


P = 13
NONZERO = tuple(range(1, P))
MIXED_COLOURS = (P - 1) ** 2

# Ambient owner-pivot coordinates.
U0 = 0
KA = 1
KB = 2
KPRIME = 3
U4 = 4
U5 = 5
SOURCE = 6
TARGET_A = 7
TARGET_B = 8
AMBIENT_DIMENSION = 9


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dot(left: Sequence[int], right: Sequence[int]) -> int:
    require(
        len(left) == len(right) == AMBIENT_DIMENSION,
        "owner-pivot dot product has the wrong dimension",
    )
    return sum(a * b for a, b in zip(left, right)) % P


def basis(index: int, coefficient: int = 1) -> tuple[int, ...]:
    values = [0] * AMBIENT_DIMENSION
    values[index] = coefficient % P
    return tuple(values)


def add_vectors(*vectors: Sequence[int]) -> tuple[int, ...]:
    require(vectors, "cannot add an empty vector family")
    require(
        all(len(vector) == AMBIENT_DIMENSION for vector in vectors),
        "owner-pivot vector has the wrong dimension",
    )
    return tuple(
        sum(vector[index] for vector in vectors) % P
        for index in range(AMBIENT_DIMENSION)
    )


def scale_vector(scalar: int, vector: Sequence[int]) -> tuple[int, ...]:
    require(
        len(vector) == AMBIENT_DIMENSION,
        "owner-pivot scale has the wrong dimension",
    )
    return tuple(scalar * value % P for value in vector)


def star_row(
    scalar_weights: Sequence[int],
    pivot_coordinate: int,
) -> tuple[int, ...]:
    """The normalized THM-2309 row, including either target graft."""
    omitted_weight = scalar_weights[U0]
    row = [0] * AMBIENT_DIMENSION
    row[U0] = scalar_weights[pivot_coordinate]
    row[pivot_coordinate] = -omitted_weight
    if pivot_coordinate == KA:
        row[TARGET_A] = -omitted_weight
    elif pivot_coordinate == KB:
        row[TARGET_B] = -omitted_weight
    return tuple(value % P for value in row)


def owner_pivot_rows(
    scalar_weights: Sequence[int],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        star_row(scalar_weights, coordinate)
        for coordinate in (KA, KB, KPRIME, U4, U5, SOURCE)
    )


def fractional_part(value: Fraction) -> Fraction:
    return value - value.numerator // value.denominator


def circular_distance(value: Fraction) -> Fraction:
    reduced = fractional_part(value)
    return min(reduced, 1 - reduced)


def danger(width: int, value: Fraction) -> int:
    require(width in (1, 2), "only ordinary and guard widths are in scope")
    return int(circular_distance(value) < Fraction(width, 14))


def target_boundaries(width: int) -> set[Fraction]:
    half_width = Fraction(width, 14)
    return {
        fractional_part(Fraction(shift, P) + sign * half_width)
        for shift in range(P)
        for sign in (-1, 1)
    }


def successor_boundaries(width: int) -> set[Fraction]:
    root_half_width = Fraction(width, 14 * P)
    return {
        fractional_part(Fraction(sheet, P) + sign * root_half_width)
        for sheet in range(P)
        for sign in (-1, 1)
    }


def cell_midpoints(boundaries: Iterable[Fraction]) -> tuple[Fraction, ...]:
    ordered = tuple(sorted(set(boundaries)))
    require(len(ordered) >= 2, "rational refinement has too few boundaries")
    values = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += 1
        values.append(fractional_part((left + right) / 2))
    return tuple(values)


def minus_shift_profile(width: int, phase: Fraction) -> tuple[int, ...]:
    return tuple(
        danger(width, phase - Fraction(shift, P)) for shift in range(P)
    )


def plus_shift_profile(
    width: int, phase: Fraction, slope: int
) -> tuple[int, ...]:
    return tuple(
        danger(width, phase + Fraction(slope * shift, P))
        for shift in range(P)
    )


def rank_one_mixed_energy(
    first_count: int, second_count: int
) -> Fraction:
    """Mixed Fourier energy of two binary one-coordinate profiles."""
    return Fraction(
        first_count
        * (P - first_count)
        * second_count
        * (P - second_count),
        P**4,
    )


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    # Scalar multiplication by the omitted weight does not change any
    # row-space, quotient, or annihilator assertion.  Normalize w_(u0)=1
    # and exhaust all 12^5 projective unit-weight tuples.  This represents
    # all 12^6 raw nonzero tuples exactly.
    normalized_weight_count = 0
    raw_weight_count = len(NONZERO) ** 6
    eta_a = add_vectors(basis(TARGET_A), basis(KA, -1))
    eta_b = add_vectors(basis(TARGET_B), basis(KB, -1))
    for remaining_weights in product(NONZERO, repeat=5):
        scalar_weights = (
            1,
            *remaining_weights,
            0,
            0,
            0,
        )
        require(
            len(scalar_weights) == AMBIENT_DIMENSION,
            "scalar owner-pivot word has the wrong dimension",
        )
        rows = owner_pivot_rows(scalar_weights)
        require(len(rows) == 6, "owner pivot stopped having six rows")
        require(
            all(dot(row, scalar_weights) == 0 for row in rows),
            "owner-pivot row left the scalar relation kernel",
        )

        omitted_weight = scalar_weights[U0]
        helper_weight = scalar_weights[KA]
        h = add_vectors(
            basis(U0, helper_weight),
            basis(KA, -omitted_weight),
        )
        require(
            dot(h, scalar_weights) == 0,
            "balanced helper vector left the primal relation space",
        )

        grafted_a = rows[0]
        primal_residual = add_vectors(h, scale_vector(-1, grafted_a))
        require(
            primal_residual == basis(TARGET_A, omitted_weight),
            "balanced helper stopped representing the pure a primal class",
        )
        helper_inverse = pow(helper_weight, -1, P)
        normalized_h = scale_vector(helper_inverse, h)
        require(
            add_vectors(
                normalized_h,
                scale_vector(-helper_inverse, grafted_a),
            )
            == basis(TARGET_A, helper_inverse * omitted_weight),
            "normalized pure primal class changed",
        )

        witness_row = rows[2]
        witness_value = dot(normalized_h, witness_row)
        require(
            witness_value == scalar_weights[KPRIME] != 0,
            "ungrafted star row stopped detecting the false co-shift",
        )
        require(
            dot(scalar_weights, witness_row) == 0,
            "scalar gauge stopped annihilating the witness row",
        )
        # Bilinearity now checks all thirteen gauges at once:
        # (g+alpha*w).ell=g.ell+alpha*(w.ell)=g.ell.
        require(
            (
                witness_value
                + (P - 1) * dot(scalar_weights, witness_row)
            )
            % P
            == witness_value,
            "scalar gauge unexpectedly repaired the false co-shift",
        )
        require(
            dot(basis(U0), witness_row) == scalar_weights[KPRIME] != 0,
            "omitted-coordinate shift unexpectedly became lawful",
        )

        require(
            all(dot(eta_a, row) == 0 for row in rows),
            "true a-target dipole left the annihilator",
        )
        require(
            all(dot(eta_b, row) == 0 for row in rows),
            "true b-target dipole left the annihilator",
        )
        require(
            dot(eta_a, basis(TARGET_A)) == 1
            and dot(eta_a, basis(TARGET_B)) == 0,
            "a-target dipole lost its pure dual pairing",
        )
        require(
            dot(eta_b, basis(TARGET_B)) == 1
            and dot(eta_b, basis(TARGET_A)) == 0,
            "b-target dipole lost its pure dual pairing",
        )
        require(
            eta_a[TARGET_A] != 0
            and eta_b[TARGET_B] != 0
            and scalar_weights[TARGET_A] == scalar_weights[TARGET_B] == 0,
            "a true target dipole collapsed into the scalar gauge",
        )
        normalized_weight_count += 1

    require(
        normalized_weight_count == len(NONZERO) ** 5 == 248832,
        "projective owner-weight exhaustion changed",
    )
    require(
        normalized_weight_count * (P - 1) == raw_weight_count == 2985984,
        "projective normalization stopped covering all raw weights",
    )

    # The width-one and width-two target/successor boundary sets coincide
    # separately and are disjoint from one another.  Their union is the
    # complete 52-cell open refinement used for every phase below.
    boundaries = {}
    for width in (1, 2):
        target = target_boundaries(width)
        successor = successor_boundaries(width)
        require(len(target) == 26, f"L={width}: boundary count changed")
        require(
            target == successor,
            f"L={width}: shift and successor boundaries diverged",
        )
        boundaries[width] = target
    require(
        boundaries[1].isdisjoint(boundaries[2]),
        "unit and guard boundaries unexpectedly overlap",
    )
    samples = cell_midpoints(boundaries[1] | boundaries[2])
    require(len(samples) == 52, "combined rational refinement changed")

    for phase in samples:
        for width in (1, 2):
            shifted_count = sum(minus_shift_profile(width, phase))
            require(
                shifted_count == 2 * width - danger(width, P * phase),
                f"L={width}: exact shift identity failed at {phase}",
            )

    deep_profiles = tuple(
        minus_shift_profile(1, phase)
        for phase in samples
        if danger(1, phase) == 0
    )
    require(len(deep_profiles) == 45, "deep-safe cell count changed")
    require(
        all(profile[0] == 0 for profile in deep_profiles),
        "deep-safe profile lost its zero axis",
    )
    deep_count_multiplicities = {
        count: sum(1 for profile in deep_profiles if sum(profile) == count)
        for count in (1, 2)
    }
    require(
        sum(deep_count_multiplicities.values()) == len(deep_profiles)
        and all(
            multiplicity > 0
            for multiplicity in deep_count_multiplicities.values()
        ),
        "deep-safe count orbit decomposition changed",
    )

    phase_profile_counts = {}
    total_minima = {}
    actual_energy_minima = {}
    for width in (1, 2):
        failed_phases = tuple(
            phase for phase in samples if danger(width, phase) == 1
        )
        expected_failed_cells = 7 if width == 1 else 13
        require(
            len(failed_phases) == expected_failed_cells,
            f"L={width}: failed-anchor cell count changed",
        )

        profile_count = 0
        minimum_total = None
        minimum_energy = None
        universal_slots = 11 - 2 * width
        for failed_phase in failed_phases:
            failed_danger = minus_shift_profile(width, failed_phase)
            repair_safe = tuple(1 - bit for bit in failed_danger)
            require(
                repair_safe[0] == 0,
                f"L={width}: failed factor lost its zero repair anchor",
            )
            failed_successor = danger(width, P * failed_phase)
            for helper_phase in samples:
                helper_successor = danger(1, P * helper_phase)
                for slope in NONZERO:
                    helper_danger = plus_shift_profile(
                        1, helper_phase, slope
                    )
                    helper_safe = tuple(1 - bit for bit in helper_danger)
                    require(
                        sum(helper_danger)
                        == 2 - helper_successor,
                        "nonzero helper slope stopped permuting all shifts",
                    )
                    overlap = sum(
                        failed * helper
                        for failed, helper in zip(
                            failed_danger, helper_danger
                        )
                    )
                    joint_safe = sum(
                        repair * helper
                        for repair, helper in zip(
                            repair_safe, helper_safe
                        )
                    )
                    exact_joint_safe = (
                        11
                        - 2 * width
                        + failed_successor
                        + helper_successor
                        + overlap
                    )
                    require(
                        joint_safe == exact_joint_safe >= universal_slots,
                        "exact two-probe safe-slot identity changed",
                    )
                    require(
                        repair_safe[0] * helper_safe[0] == 0,
                        "duplicate-probe role lost its zero t-axis",
                    )

                    # Corner and energy depend on the deep profile only
                    # through its exact count.  The 45 profiles were
                    # exhausted above; check each of their two count orbits
                    # once here and retain the full profile multiplicity.
                    for (
                        deep_count,
                        deep_multiplicity,
                    ) in deep_count_multiplicities.items():
                        total = deep_count * joint_safe
                        corner = Fraction(total, P * P)
                        require(
                            corner >= Fraction(universal_slots, P * P),
                            "mixed duplicate-probe corner floor failed",
                        )
                        energy = rank_one_mixed_energy(
                            deep_count, joint_safe
                        )
                        require(
                            energy >= corner * corner / MIXED_COLOURS,
                            "144-colour Cauchy energy floor failed",
                        )
                        require(
                            energy
                            >= Fraction(
                                universal_slots**2,
                                P**4 * MIXED_COLOURS,
                            ),
                            "uniform duplicate-probe energy floor failed",
                        )
                        profile_count += deep_multiplicity
                        if minimum_total is None or total < minimum_total:
                            minimum_total = total
                        if minimum_energy is None or energy < minimum_energy:
                            minimum_energy = energy

        expected_profiles = (
            len(deep_profiles)
            * len(failed_phases)
            * len(samples)
            * len(NONZERO)
        )
        require(
            profile_count == expected_profiles,
            f"L={width}: duplicate-probe profile count changed",
        )
        require(
            minimum_total == universal_slots,
            f"L={width}: exact total-mass minimum changed",
        )
        phase_profile_counts[width] = profile_count
        total_minima[width] = minimum_total
        actual_energy_minima[width] = minimum_energy

    unit_coefficient = Fraction(9, P * P * MIXED_COLOURS)
    guard_coefficient = Fraction(7, P * P * MIXED_COLOURS)
    require(unit_coefficient == Fraction(1, 2704), "unit constant changed")
    require(
        guard_coefficient == Fraction(7, 24336),
        "guard constant changed",
    )
    unit_energy = Fraction(9**2, P**4 * MIXED_COLOURS)
    guard_energy = Fraction(7**2, P**4 * MIXED_COLOURS)
    require(
        unit_energy == Fraction(9, 456976),
        "unit energy constant changed",
    )
    require(
        guard_energy == Fraction(49, 4112784),
        "guard energy constant changed",
    )

    print("THM-2384 owner-pivot primal/dual repair exact companion")
    print(
        "owner weights: "
        f"projective={normalized_weight_count}; raw={raw_weight_count}; "
        "primal_pure=yes; false_dual_witness=universal"
    )
    print(
        "true dual dipoles: eta_a=e_a-e_ka; eta_b=e_b-e_kb; "
        "all normalized owner weights checked"
    )
    print(
        "rational refinement: "
        f"L1 boundaries={len(boundaries[1])}; "
        f"L2 boundaries={len(boundaries[2])}; "
        f"combined cells={len(samples)}; slopes={len(NONZERO)}"
    )
    for width in (1, 2):
        print(
            f"L={width}: profiles={phase_profile_counts[width]}; "
            f"total_min={total_minima[width]}; "
            f"actual_cell_energy_min={actual_energy_minima[width]}"
        )
    print(
        "coefficient floors: "
        f"unit={unit_coefficient}; guard={guard_coefficient}"
    )
    print(
        "energy floors: "
        f"unit={unit_energy}; guard={guard_energy}"
    )
    print(
        "VERDICT: primal purity does not make a target covector; "
        "the two-probe corner survives without target typing"
    )


if __name__ == "__main__":
    main()
