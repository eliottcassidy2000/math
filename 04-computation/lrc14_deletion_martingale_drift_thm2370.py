#!/usr/bin/env python3
"""Exact controls for THM-2370's deletion-martingale theorem.

The companion independently reconstructs the THM-2367 Boolean hostile,
computes its target-drift energy, and then verifies the sharp clone
construction for several deletion depths.  It also computes the complete
Boolean Walsh spectrum of every tested clone cube.  All assertions remain
active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256


P = 13
GRID = 16562
ZERO = Fraction(0)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


EXCLUDED = (
    (16555, 16562),
    (0, 13),
    (1625, 1651),
    (2457, 2463),
    (3263, 3289),
    (4907, 4927),
    (7449, 7455),
    (9087, 9113),
    (10725, 10751),
    (12363, 12389),
)


def grid_danger(
    cell: int,
    speed: int,
    shift: int,
    length_units: int = 1,
) -> bool:
    """Test ||speed*x+shift/13||<length_units/14 at a GRID-cell midpoint."""
    modulus = 2 * GRID * P
    numerator = speed * (2 * cell + 1) * P + shift * 2 * GRID
    distance = min(numerator % modulus, (-numerator) % modulus)
    return 14 * distance < length_units * modulus


def matrix_subtract(
    left: tuple[tuple[Fraction, ...], ...],
    right: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(left[row][column] - right[row][column] for column in range(P))
        for row in range(P)
    )


def matrix_scale(
    matrix: tuple[tuple[Fraction, ...], ...],
    scalar: Fraction,
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(scalar * matrix[row][column] for column in range(P))
        for row in range(P)
    )


def matrix_add(
    left: tuple[tuple[Fraction, ...], ...],
    right: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(left[row][column] + right[row][column] for column in range(P))
        for row in range(P)
    )


def circulant_projection(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[tuple[Fraction, ...], ...]:
    """Average the simultaneous (r,t)-translation orbit."""
    profile = []
    for difference in range(P):
        profile.append(
            sum(
                matrix[(endpoint + difference) % P][endpoint]
                for endpoint in range(P)
            )
            / P
        )
    return tuple(
        tuple(profile[(root - endpoint) % P] for endpoint in range(P))
        for root in range(P)
    )


def drift(
    matrix: tuple[tuple[Fraction, ...], ...],
) -> Fraction:
    projected = circulant_projection(matrix)
    return (
        sum(
            (
                matrix[root][endpoint]
                - projected[root][endpoint]
            )
            ** 2
            for root in range(P)
            for endpoint in range(P)
        )
        / (P * P)
    )


def normalized_walsh(
    values: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    """Return coefficients for chi_A(x)=(-1)^(sum_(j in A)x_j)."""
    size = len(values)
    require(size > 0 and size & (size - 1) == 0, "Walsh size is not dyadic")
    coefficients = list(values)
    half_block = 1
    while half_block < size:
        block = 2 * half_block
        for start in range(0, size, block):
            for offset in range(half_block):
                low = coefficients[start + offset]
                high = coefficients[start + half_block + offset]
                coefficients[start + offset] = low + high
                coefficients[start + half_block + offset] = low - high
        half_block = block
    return tuple(value / size for value in coefficients)


def main() -> None:
    excluded_cells = sum(right - left for left, right in EXCLUDED)
    require(excluded_cells == 182, "hostile deletion count changed")
    require(Fraction(excluded_cells, GRID) == Fraction(1, 91), "mask mass changed")

    segment_string = (
        "D=16562;excluded=[16555,16562);[0,13);[1625,1651);"
        "[2457,2463);[3263,3289);[4907,4927);[7449,7455);"
        "[9087,9113);[10725,10751);[12363,12389)"
    )
    segment_hash = sha256(segment_string.encode()).hexdigest()
    require(
        segment_hash
        == "4e458401005911de96cd61a26638cee0c2c75b8aa033df0c0298a3438c0514eb",
        "hostile segment hash changed",
    )

    all_bits = (1 << GRID) - 1
    mask_bits = 0
    for cell in range(GRID):
        if not any(left <= cell < right for left, right in EXCLUDED):
            mask_bits |= 1 << cell
    require(mask_bits.bit_count() == GRID - excluded_cells, "Boolean mask changed")

    deep_danger_bits: list[int] = []
    deep_safe_bits: list[int] = []
    graft_safe_bits: list[int] = []
    for shift in range(P):
        deep = 0
        graft = 0
        for cell in range(GRID):
            if grid_danger(cell, 13, -shift):
                deep |= 1 << cell
            if grid_danger(cell, 7, shift):
                graft |= 1 << cell
        deep_danger_bits.append(deep)
        deep_safe_bits.append(all_bits ^ deep)
        graft_safe_bits.append(all_bits ^ graft)

    raw_counts: list[list[int]] = []
    masked_counts: list[list[int]] = []
    expected_profile = (0, 1078) + (2002,) * 10 + (1078,)
    for root in range(P):
        raw_row = []
        masked_row = []
        for endpoint in range(P):
            base = (
                deep_danger_bits[root]
                & deep_safe_bits[endpoint]
                & graft_safe_bits[endpoint]
            )
            raw_row.append(base.bit_count())
            masked_row.append((base & mask_bits).bit_count())
        raw_counts.append(raw_row)
        masked_counts.append(masked_row)

    for endpoint in range(P):
        by_difference = tuple(
            masked_counts[(endpoint + difference) % P][endpoint]
            for difference in range(P)
        )
        require(by_difference == expected_profile, "masked table is not circulant")

    raw = tuple(
        tuple(Fraction(raw_counts[root][endpoint], GRID) for endpoint in range(P))
        for root in range(P)
    )
    masked = tuple(
        tuple(
            Fraction(masked_counts[root][endpoint], GRID)
            for endpoint in range(P)
        )
        for root in range(P)
    )
    deleted = matrix_subtract(raw, masked)

    require(
        all(
            raw[root][endpoint] >= masked[root][endpoint] >= 0
            for root in range(P)
            for endpoint in range(P)
        ),
        "positive deletion ordering changed",
    )
    require(
        all(raw[index][index] == masked[index][index] == 0 for index in range(P)),
        "diagonal zero changed",
    )
    require(circulant_projection(masked) == masked, "masked table lost circulation")

    raw_drift = drift(raw)
    masked_drift = drift(masked)
    deleted_drift = drift(deleted)
    require(masked_drift == 0, "terminal table still drifts")
    require(
        raw_drift == deleted_drift == Fraction(852, 11589168409),
        "exact hostile drift changed",
    )
    require(
        matrix_subtract(
            matrix_subtract(raw, circulant_projection(raw)),
            matrix_subtract(deleted, circulant_projection(deleted)),
        )
        == tuple(tuple(ZERO for _ in range(P)) for _ in range(P)),
        "projected deletion no longer equals initial drift",
    )

    clone_depths = (1, 2, 3, 5, 7, 13)
    clone_rows = []
    walsh_rows = []
    for depth in clone_depths:
        layer = matrix_scale(deleted, Fraction(1, depth))
        layer_drift = drift(layer)
        require(
            layer_drift == raw_drift / (depth * depth),
            "single-layer n^-2 equality changed",
        )
        require(
            depth * layer_drift == raw_drift / depth,
            "summed n^-1 equality changed",
        )
        require(
            Fraction(1) - Fraction(1, 91 * depth)
            == Fraction(91 * depth - 1, 91 * depth),
            "clone-mask mass changed",
        )

        for stage in range(depth + 1):
            table = matrix_add(
                masked,
                matrix_scale(deleted, Fraction(depth - stage, depth)),
            )
            require(
                drift(table)
                == Fraction((depth - stage) ** 2, depth * depth) * raw_drift,
                "stagewise quadratic drift changed",
            )
            if stage < depth:
                next_table = matrix_add(
                    masked,
                    matrix_scale(
                        deleted,
                        Fraction(depth - stage - 1, depth),
                    ),
                )
                require(
                    matrix_subtract(table, next_table) == layer,
                    "clone deletion layer changed",
                )

        clone_rows.append((depth, layer_drift))

        # Full Boolean cube: applying the masks indexed by x leaves
        # H_x=C+(1-|x|/n)(T-C).  Under chi_A(x)=(-1)^(A.x), its target
        # coefficient is QT/2 at A=empty, QT/(2n) at singletons, and zero
        # in every higher degree.
        cube_values = tuple(
            Fraction(depth - mask.bit_count(), depth)
            for mask in range(1 << depth)
        )
        coefficients = normalized_walsh(cube_values)
        require(coefficients[0] == Fraction(1, 2), "clone-cube mean changed")
        require(
            all(
                coefficients[1 << coordinate] == Fraction(1, 2 * depth)
                for coordinate in range(depth)
            ),
            "clone-cube singleton coefficient changed",
        )
        require(
            all(
                coefficient == 0
                for mask, coefficient in enumerate(coefficients)
                if mask != 0 and mask.bit_count() >= 2
            ),
            "clone cube acquired a mixed Walsh component",
        )

        coefficient_energies = tuple(value * value for value in coefficients)
        nonconstant_vector_energies = tuple(
            (mask, value)
            for mask, value in enumerate(coefficient_energies)
            if mask != 0 and value != 0
        )
        subcube_dirichlet = []
        for coordinate_set in range(1 << depth):
            value = 2 * sum(
                energy
                for mask, energy in nonconstant_vector_energies
                if mask & coordinate_set
            )
            require(
                value
                == Fraction(coordinate_set.bit_count(), 2 * depth * depth),
                "clone complete-subcube Dirichlet law changed",
            )
            subcube_dirichlet.append(value)

        full_mask = (1 << depth) - 1
        # Since every checked D_U is |U|/(2n^2), additivity of cardinality
        # makes every disjoint pair defect zero.  The following exact count
        # records the entire disjoint-pair universe without a slow 3^n loop.
        disjoint_pair_count = sum(
            1 << (depth - first.bit_count())
            for first in range(1 << depth)
        )
        require(
            disjoint_pair_count == 3**depth,
            "disjoint subcube-pair count changed",
        )

        # Global Boolean complementation preserves the entire squared
        # Dirichlet bank but swaps the terminal values C and T.
        reversed_values = tuple(
            Fraction(mask.bit_count(), depth)
            for mask in range(1 << depth)
        )
        reversed_coefficients = normalized_walsh(reversed_values)
        require(
            all(
                reversed_coefficients[mask] ** 2
                == coefficients[mask] ** 2
                for mask in range(1 << depth)
            ),
            "global clone reversal changed the squared Walsh bank",
        )
        require(
            cube_values[full_mask] == 0
            and reversed_values[full_mask] == 1,
            "global clone reversal stopped swapping terminal drift",
        )

        # Do not confuse the Hilbert-valued packet cube with the scalar
        # quadratic x -> D(H_x)=D(T)(1-|x|/n)^2.  The latter has genuine
        # degree-two Walsh terms and positive pair defects.
        scalar_energy_values = tuple(
            raw_drift * value * value
            for value in cube_values
        )
        scalar_energy_coefficients = normalized_walsh(scalar_energy_values)
        require(
            scalar_energy_coefficients[0]
            == raw_drift * Fraction(depth + 1, 4 * depth),
            "scalar drift-cube mean changed",
        )
        require(
            all(
                scalar_energy_coefficients[1 << coordinate]
                == raw_drift * Fraction(1, 2 * depth)
                for coordinate in range(depth)
            ),
            "scalar drift-cube singleton coefficient changed",
        )
        require(
            all(
                scalar_energy_coefficients[
                    (1 << first) | (1 << second)
                ]
                == raw_drift * Fraction(1, 2 * depth * depth)
                for first in range(depth)
                for second in range(first + 1, depth)
            ),
            "scalar drift-cube pair coefficient changed",
        )
        require(
            all(
                coefficient == 0
                for mask, coefficient in enumerate(scalar_energy_coefficients)
                if mask.bit_count() >= 3
            ),
            "scalar drift cube acquired degree at least three",
        )
        scalar_energy_squares = tuple(
            value * value
            for value in scalar_energy_coefficients
        )
        nonconstant_scalar_energies = tuple(
            (mask, value)
            for mask, value in enumerate(scalar_energy_squares)
            if mask != 0 and value != 0
        )
        scalar_subcube_dirichlet = []
        for coordinate_set in range(1 << depth):
            value = 2 * sum(
                energy
                for mask, energy in nonconstant_scalar_energies
                if mask & coordinate_set
            )
            size = coordinate_set.bit_count()
            require(
                value
                == raw_drift**2
                * Fraction(
                    size * (2 * depth * depth + 2 * depth - size - 1),
                    4 * depth**4,
                ),
                "scalar drift complete-subcube law changed",
            )
            scalar_subcube_dirichlet.append(value)
        if depth >= 2:
            require(
                scalar_subcube_dirichlet[1]
                + scalar_subcube_dirichlet[2]
                - scalar_subcube_dirichlet[3]
                == raw_drift**2 * Fraction(1, 2 * depth**4),
                "scalar drift pair defect changed",
            )

        walsh_rows.append(
            (
                depth,
                coefficients[0],
                coefficients[1],
                subcube_dirichlet[full_mask],
                disjoint_pair_count,
            )
        )

    # Unequal clone weights make both universal bounds strict.
    unequal_weights = (Fraction(1, 2), Fraction(1, 3), Fraction(1, 6))
    require(sum(unequal_weights) == 1, "unequal weights stopped partitioning")
    unequal_drifts = tuple(weight * weight * raw_drift for weight in unequal_weights)
    require(
        max(unequal_drifts) > raw_drift / 9
        and sum(unequal_drifts) > raw_drift / 3,
        "unequal-weight hostile stopped being strict",
    )

    # Omitting the final clone leaves exactly one n^-2 unit of drift.
    for depth in (2, 3, 5, 13):
        incomplete = matrix_add(
            masked,
            matrix_scale(deleted, Fraction(1, depth)),
        )
        require(
            drift(incomplete) == raw_drift / (depth * depth),
            "omitted-final-mask control changed",
        )

    # Sharp principal angle for a table whose whole base slice vanishes.
    def anchored_value(root: int, target: int, endpoint: int) -> Fraction:
        if target == 0 and endpoint == 0:
            return ZERO
        return Fraction(root != endpoint)

    anchored_norm = (
        sum(
            anchored_value(root, target, endpoint) ** 2
            for root in range(P)
            for target in range(P)
            for endpoint in range(P)
        )
        / P**3
    )
    anchored_projection: dict[tuple[int, int, int], Fraction] = {}
    for root in range(P):
        for target in range(P):
            for endpoint in range(P):
                anchored_projection[(root, target, endpoint)] = (
                    sum(
                        anchored_value(
                            (root + diagonal_shift) % P,
                            (target + first_shift) % P,
                            (endpoint + diagonal_shift) % P,
                        )
                        for first_shift in range(P)
                        for diagonal_shift in range(P)
                    )
                    / P**2
                )
    anchored_drift = (
        sum(
            (
                anchored_value(root, target, endpoint)
                - anchored_projection[(root, target, endpoint)]
            )
            ** 2
            for root in range(P)
            for target in range(P)
            for endpoint in range(P)
        )
        / P**3
    )
    anchored_scale = Fraction(P * P - 1, P * P)
    require(
        all(
            anchored_projection[(root, target, endpoint)]
            == anchored_scale * Fraction(root != endpoint)
            for root in range(P)
            for target in range(P)
            for endpoint in range(P)
        ),
        "anchored-zero projection changed",
    )
    require(
        anchored_norm == Fraction(2016, 2197)
        and anchored_drift == Fraction(2016, 371293)
        and anchored_drift == anchored_norm / P**2,
        "sharp anchored-zero drift ratio changed",
    )

    eligible_colours = (P - 1) * (P * P - 1)
    require(
        eligible_colours == 2016
        and P * eligible_colours == 26208,
        "THM-2365 coefficient count changed",
    )

    matrix_string = ",".join(
        str(masked_counts[(endpoint + difference) % P][endpoint])
        for endpoint in range(P)
        for difference in range(P)
    )
    matrix_hash = sha256(matrix_string.encode()).hexdigest()
    require(
        matrix_hash
        == "e79dc643a13e9f4aafecd3f6b007952885902449c7aa8e5825f7f26bae4c7825",
        "masked matrix hash changed",
    )

    print("THM-2370 deletion-martingale exact referee")
    print(f"grid/deleted cells: {GRID}/{excluded_cells}")
    print(f"terminal Boolean mask mass: {Fraction(GRID-excluded_cells, GRID)}")
    print(f"segment sha256: {segment_hash}")
    print(f"masked matrix sha256: {matrix_hash}")
    print(f"initial/terminal drift: {raw_drift}/{masked_drift}")
    print(f"deleted tensor drift: {deleted_drift}")
    for depth, layer_drift in clone_rows:
        print(
            f"n={depth}: mask_mass={Fraction(91*depth-1,91*depth)}; "
            f"one_layer={layer_drift}; sum={depth*layer_drift}"
        )
    for depth, mean, singleton, full_dirichlet, pair_count in walsh_rows:
        print(
            f"cube n={depth}: mean={mean}; singleton={singleton}; "
            f"mixed=0; D_full/D(T)={full_dirichlet}; "
            f"zero_ordered_disjoint_pair_defects={pair_count}"
        )
    print(
        "anchored-zero: "
        f"norm={anchored_norm}; drift={anchored_drift}; "
        f"ratio={anchored_drift/anchored_norm}"
    )
    print(f"eligible deep-target coefficients: {eligible_colours}")
    print("VERDICT: n^-2 maximum and n^-1 summed deletion floors are sharp")


if __name__ == "__main__":
    main()
