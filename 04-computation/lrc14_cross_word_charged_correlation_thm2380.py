#!/usr/bin/env python3
"""Dependency-free exact companion for the THM-2380 proof candidate.

The calculation is over the target group F_13^2.  Cyclotomic values are
represented in Q[zeta_13] using the exact relation

    1 + zeta + ... + zeta^12 = 0.

Every validity check uses ``require`` rather than ``assert`` so that the exact
audit is unchanged under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
import sys


P = 13
N = P * P
ZERO = Fraction(0)
ONE = Fraction(1)

Point = tuple[int, int]
Cyclo = tuple[Fraction, ...]
Gaussian = tuple[Fraction, Fraction]
CycloI = tuple[Cyclo, Cyclo]


def require(condition: bool, message: str) -> None:
    """Raise on a failed exact check, including under ``python -O``."""
    if not condition:
        raise RuntimeError(message)


def canonical(values: list[Fraction]) -> Cyclo:
    """Return canonical Q[zeta_13] coordinates of length 13.

    Subtracting the coefficient of zeta^12 implements the sole rational
    relation among the thirteen powers of a primitive thirteenth root.
    """
    require(len(values) == P, "cyclotomic vector has the wrong length")
    anchor = values[-1]
    return tuple(value - anchor for value in values)


def cyclo_zero() -> Cyclo:
    return tuple(ZERO for _ in range(P))


def cyclo_rational(value: Fraction | int) -> Cyclo:
    entries = [ZERO for _ in range(P)]
    entries[0] = Fraction(value)
    return canonical(entries)


def cyclo_add(left: Cyclo, right: Cyclo) -> Cyclo:
    return canonical(
        [left[index] + right[index] for index in range(P)]
    )


def cyclo_subtract(left: Cyclo, right: Cyclo) -> Cyclo:
    return canonical(
        [left[index] - right[index] for index in range(P)]
    )


def cyclo_scale(value: Cyclo, scalar: Fraction | int) -> Cyclo:
    factor = Fraction(scalar)
    return canonical([factor * entry for entry in value])


def cyclo_phase(value: Cyclo, exponent: int) -> Cyclo:
    """Multiply by zeta**exponent."""
    entries = [ZERO for _ in range(P)]
    shift = exponent % P
    for index, entry in enumerate(value):
        entries[(index + shift) % P] += entry
    return canonical(entries)


def cyclo_multiply(left: Cyclo, right: Cyclo) -> Cyclo:
    entries = [ZERO for _ in range(P)]
    for left_index, left_entry in enumerate(left):
        if left_entry == 0:
            continue
        for right_index, right_entry in enumerate(right):
            if right_entry == 0:
                continue
            entries[(left_index + right_index) % P] += (
                left_entry * right_entry
            )
    return canonical(entries)


def cyclo_conjugate(value: Cyclo) -> Cyclo:
    entries = [ZERO for _ in range(P)]
    for index, entry in enumerate(value):
        entries[(-index) % P] += entry
    return canonical(entries)


def cyclo_norm_squared(value: Cyclo) -> Cyclo:
    return cyclo_multiply(value, cyclo_conjugate(value))


def cyclo_sum(values: list[Cyclo] | tuple[Cyclo, ...]) -> Cyclo:
    total = cyclo_zero()
    for value in values:
        total = cyclo_add(total, value)
    return total


def cyclo_as_rational(value: Cyclo, message: str) -> Fraction:
    require(
        all(entry == 0 for entry in value[1:]),
        message,
    )
    return value[0]


def root(exponent: int) -> Cyclo:
    return cyclo_phase(cyclo_rational(1), exponent)


def gaussian_add(left: Gaussian, right: Gaussian) -> Gaussian:
    return (left[0] + right[0], left[1] + right[1])


def gaussian_conjugate(value: Gaussian) -> Gaussian:
    return (value[0], -value[1])


def gaussian_multiply(left: Gaussian, right: Gaussian) -> Gaussian:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def gaussian_norm_squared(value: Gaussian) -> Fraction:
    return value[0] * value[0] + value[1] * value[1]


def cycloi_zero() -> CycloI:
    return (cyclo_zero(), cyclo_zero())


def cycloi_embed(value: Cyclo) -> CycloI:
    return (value, cyclo_zero())


def cycloi_gaussian(value: Gaussian) -> CycloI:
    return (cyclo_rational(value[0]), cyclo_rational(value[1]))


def cycloi_add(left: CycloI, right: CycloI) -> CycloI:
    return (
        cyclo_add(left[0], right[0]),
        cyclo_add(left[1], right[1]),
    )


def cycloi_subtract(left: CycloI, right: CycloI) -> CycloI:
    return (
        cyclo_subtract(left[0], right[0]),
        cyclo_subtract(left[1], right[1]),
    )


def cycloi_scale(
    value: CycloI,
    scalar: Fraction | int,
) -> CycloI:
    return (
        cyclo_scale(value[0], scalar),
        cyclo_scale(value[1], scalar),
    )


def cycloi_multiply(left: CycloI, right: CycloI) -> CycloI:
    return (
        cyclo_subtract(
            cyclo_multiply(left[0], right[0]),
            cyclo_multiply(left[1], right[1]),
        ),
        cyclo_add(
            cyclo_multiply(left[0], right[1]),
            cyclo_multiply(left[1], right[0]),
        ),
    )


def cycloi_conjugate_i(value: CycloI) -> CycloI:
    """Conjugate only the adjoined i-coordinate."""
    return (value[0], cyclo_scale(value[1], -1))


def cycloi_conjugate(value: CycloI) -> CycloI:
    """Apply complex conjugation to both zeta_13 and i."""
    return (
        cyclo_conjugate(value[0]),
        cyclo_scale(cyclo_conjugate(value[1]), -1),
    )


def cycloi_phase(value: CycloI, exponent: int) -> CycloI:
    """Multiply a Q(zeta_13,i) value by zeta_13**exponent."""
    return (
        cyclo_phase(value[0], exponent),
        cyclo_phase(value[1], exponent),
    )


def cycloi_sum(values: list[CycloI] | tuple[CycloI, ...]) -> CycloI:
    total = cycloi_zero()
    for value in values:
        total = cycloi_add(total, value)
    return total


def dot(left: Point, right: Point) -> int:
    return (left[0] * right[0] + left[1] * right[1]) % P


def add_group(left: Point, right: Point) -> Point:
    return (
        (left[0] + right[0]) % P,
        (left[1] + right[1]) % P,
    )


def character(frequency: Point, target: Point) -> Cyclo:
    return root(dot(frequency, target))


def fourier_transform(
    group: tuple[Point, ...],
    target_values: dict[Point, Fraction],
) -> dict[Point, Cyclo]:
    transformed: dict[Point, Cyclo] = {}
    support = tuple(
        (target, value)
        for target, value in target_values.items()
        if value != 0
    )
    for frequency in group:
        total = cyclo_zero()
        for target, value in support:
            total = cyclo_add(
                total,
                cyclo_scale(character(frequency, target), value),
            )
        transformed[frequency] = total
    return transformed


def cross_correlation_at_shift(
    group: tuple[Point, ...],
    shift: Point,
    deletion_fourier: dict[Point, Cyclo],
    retained_fourier: dict[Point, Cyclo],
) -> Cyclo:
    total = cyclo_zero()
    for frequency in group:
        total = cyclo_add(
            total,
            cyclo_multiply(
                deletion_fourier[add_group(frequency, shift)],
                cyclo_conjugate(retained_fourier[frequency]),
            ),
        )
    return cyclo_scale(total, Fraction(1, N))


def character_orthogonality(
    group: tuple[Point, ...],
) -> dict[Point, Cyclo]:
    """Compute N^-1 sum_l chi_l(q) at every q, with no shortcut."""
    kernel: dict[Point, Cyclo] = {}
    for target in group:
        total = cyclo_zero()
        for frequency in group:
            total = cyclo_add(
                total,
                character(frequency, target),
            )
        kernel[target] = cyclo_scale(total, Fraction(1, N))
        require(
            kernel[target]
            == cyclo_rational(1 if target == (0, 0) else 0),
            "target character orthogonality changed",
        )
    return kernel


def cross_correlation_from_kernel(
    group: tuple[Point, ...],
    deletion: dict[Point, Fraction],
    retained: dict[Point, Fraction],
    orthogonality: dict[Point, Cyclo],
) -> dict[Point, Cyclo]:
    """Expand the frequency correlation through the exact DFT kernel."""
    deletion_support = tuple(
        (target, value)
        for target, value in deletion.items()
        if value != 0
    )
    retained_support = tuple(
        (target, value)
        for target, value in retained.items()
        if value != 0
    )
    correlation: dict[Point, Cyclo] = {}
    for shift in group:
        total = cyclo_zero()
        for deletion_target, deletion_value in deletion_support:
            for retained_target, retained_value in retained_support:
                difference = (
                    (deletion_target[0] - retained_target[0]) % P,
                    (deletion_target[1] - retained_target[1]) % P,
                )
                term = cyclo_multiply(
                    character(shift, deletion_target),
                    orthogonality[difference],
                )
                total = cyclo_add(
                    total,
                    cyclo_scale(
                        term,
                        deletion_value * retained_value,
                    ),
                )
        correlation[shift] = total
    return correlation


def target_pair_energy(
    group: tuple[Point, ...],
    shift: Point,
    twist: int,
    deletion: dict[Point, Fraction],
    retained: dict[Point, Fraction],
) -> Cyclo:
    total = cyclo_zero()
    for target in group:
        deletion_value = deletion[target]
        retained_value = retained[target]
        if deletion_value == 0 and retained_value == 0:
            continue
        total = cyclo_add(
            total,
            cyclo_rational(
                deletion_value * deletion_value
                + retained_value * retained_value
            ),
        )
        charged_product = deletion_value * retained_value
        if charged_product != 0:
            exponent = dot(shift, target) - twist
            total = cyclo_add(
                total,
                cyclo_scale(
                    cyclo_add(root(exponent), root(-exponent)),
                    charged_product,
                ),
            )
    return total


def fourier_pair_energy(
    group: tuple[Point, ...],
    shift: Point,
    twist: int,
    deletion_fourier: dict[Point, Cyclo],
    retained_fourier: dict[Point, Cyclo],
) -> Cyclo:
    total = cyclo_zero()
    for frequency in group:
        value = cyclo_add(
            deletion_fourier[add_group(frequency, shift)],
            cyclo_phase(retained_fourier[frequency], twist),
        )
        total = cyclo_add(total, cyclo_norm_squared(value))
    return cyclo_scale(total, Fraction(1, N))


def main() -> None:
    # Keep the stored exact transcript byte-identical on Windows and POSIX.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    group = tuple(
        (first, second)
        for first in range(P)
        for second in range(P)
    )
    require(len(group) == N, "target group size changed")

    # Rational target-side data.  The four shared cells have products
    # 6, 4, -10, -3.  The final two cells test one-sided support.
    deletion = {target: ZERO for target in group}
    retained = {target: ZERO for target in group}
    deletion[(0, 0)] = Fraction(2)
    retained[(0, 0)] = Fraction(3)
    deletion[(1, 0)] = Fraction(1)
    retained[(1, 0)] = Fraction(4)
    deletion[(0, 2)] = Fraction(-2)
    retained[(0, 2)] = Fraction(5)
    deletion[(3, 4)] = Fraction(3)
    retained[(3, 4)] = Fraction(-1)
    deletion[(5, 6)] = Fraction(7)
    retained[(7, 8)] = Fraction(2)

    shared_product = {
        target: deletion[target] * retained[target]
        for target in group
    }
    expected_cross_energy = sum(
        value * value
        for target, value in shared_product.items()
        if target != (0, 0)
    )
    require(
        expected_cross_energy == 125,
        "shared target energy fixture changed",
    )

    deletion_fourier = fourier_transform(group, deletion)
    retained_fourier = fourier_transform(group, retained)
    orthogonality = character_orthogonality(group)
    correlation = cross_correlation_from_kernel(
        group,
        deletion,
        retained,
        orthogonality,
    )

    # Cross-correlation is the Fourier transform of the pointwise charged
    # product d(q) conjugate(w(q)).  The exact kernel expansion checks every
    # shift; four raw frequency-side controls independently check the route
    # through D(l+delta) conjugate(W(l)).
    cross_correlation_checks = 0
    for shift in group:
        expected = cyclo_zero()
        for target in group:
            expected = cyclo_add(
                expected,
                cyclo_scale(
                    character(shift, target),
                    shared_product[target],
                ),
            )
        require(
            correlation[shift] == expected,
            "cross-correlation normalization or sign changed",
        )
        cross_correlation_checks += 1

    direct_cross_shifts = ((0, 0), (1, 0), (0, 1), (3, 5))
    direct_cross_checks = 0
    for shift in direct_cross_shifts:
        require(
            cross_correlation_at_shift(
                group,
                shift,
                deletion_fourier,
                retained_fourier,
            )
            == correlation[shift],
            "raw frequency cross-correlation control changed",
        )
        direct_cross_checks += 1

    # Inverse target DFT recovers every charged product cell.
    inverse_target_checks = 0
    for target in group:
        total = cyclo_zero()
        for shift in group:
            total = cyclo_add(
                total,
                cyclo_phase(
                    correlation[shift],
                    -dot(shift, target),
                ),
            )
        recovered = cyclo_scale(total, Fraction(1, N))
        require(
            recovered == cyclo_rational(shared_product[target]),
            "inverse target DFT changed",
        )
        inverse_target_checks += 1

    # The cyclic pair twist extracts R at t-character -1 and conjugate(R)
    # at t-character +1.  Pair energies are first calculated target-side;
    # selected independent frequency-side Parseval calculations follow.
    pair_energies: dict[tuple[Point, int], Cyclo] = {}
    pair_twist_checks = 0
    twist_weighted_sums: dict[Point, Cyclo] = {}
    for shift in group:
        energies = []
        for twist in range(P):
            energy = target_pair_energy(
                group,
                shift,
                twist,
                deletion,
                retained,
            )
            pair_energies[(shift, twist)] = energy
            energies.append(energy)

        minus_one = cyclo_scale(
            cyclo_sum(
                [
                    cyclo_phase(energies[twist], twist)
                    for twist in range(P)
                ]
            ),
            Fraction(1, P),
        )
        plus_one = cyclo_scale(
            cyclo_sum(
                [
                    cyclo_phase(energies[twist], -twist)
                    for twist in range(P)
                ]
            ),
            Fraction(1, P),
        )
        require(
            minus_one == correlation[shift]
            and plus_one == cyclo_conjugate(correlation[shift]),
            "pair-twist character extraction changed",
        )
        twist_weighted_sums[shift] = cyclo_scale(minus_one, P)
        pair_twist_checks += P

    parseval_controls = (
        ((0, 0), 0),
        ((0, 0), 4),
        ((1, 0), 0),
        ((1, 0), 4),
        ((0, 1), 0),
        ((0, 1), 4),
    )
    direct_parseval_checks = 0
    for shift, twist in parseval_controls:
        direct_energy = fourier_pair_energy(
            group,
            shift,
            twist,
            deletion_fourier,
            retained_fourier,
        )
        require(
            direct_energy == pair_energies[(shift, twist)],
            "pair-energy Parseval normalization changed",
        )
        direct_parseval_checks += 1

    # With unnormalized frequency energy sums, the target projector has
    # denominator p*N^2 = 371293.  Collapse the t-sum first, but compare the
    # resulting raw numerator at all 169 targets.
    target_projector_denominator = P * N * N
    require(
        target_projector_denominator == 371293,
        "raw target projector denominator changed",
    )
    raw_projector_checks = 0
    for target in group:
        raw_numerator = cyclo_zero()
        for shift in group:
            raw_t_sum = cyclo_scale(
                twist_weighted_sums[shift],
                N,
            )
            raw_numerator = cyclo_add(
                raw_numerator,
                cyclo_phase(raw_t_sum, -dot(shift, target)),
            )
        require(
            raw_numerator
            == cyclo_rational(
                target_projector_denominator
                * shared_product[target]
            ),
            "raw pair-twist target projector changed",
        )
        raw_projector_checks += 1

    # Parseval converts a fixed direction into the exact character-distance
    # multiplier 2-zeta^<e,q>-zeta^-<e,q>.  We compute all 169 directions and
    # independently replay three of them from the actual correlation array.
    direction_averages: dict[Point, Cyclo] = {}
    for direction in group:
        average = cyclo_zero()
        for target, value in shared_product.items():
            if target == (0, 0) or value == 0:
                continue
            exponent = dot(direction, target)
            distance = cyclo_subtract(
                cyclo_rational(2),
                cyclo_add(root(exponent), root(-exponent)),
            )
            average = cyclo_add(
                average,
                cyclo_scale(distance, value * value),
            )
        direction_averages[direction] = average

    direct_direction_checks = 0
    for direction in ((1, 0), (0, 1), (1, 1)):
        raw_total = cyclo_zero()
        for shift in group:
            difference = cyclo_subtract(
                correlation[add_group(shift, direction)],
                correlation[shift],
            )
            raw_total = cyclo_add(
                raw_total,
                cyclo_norm_squared(difference),
            )
        require(
            cyclo_scale(raw_total, Fraction(1, N))
            == direction_averages[direction],
            "fixed-direction charged Parseval control changed",
        )
        direct_direction_checks += 1

    # Complete-Cayley and coordinate-line overlap identities now use all
    # exact direction multipliers.  direction_totals retain the raw sum over
    # the 169 base shifts, matching the displayed denominators.
    direction_totals: dict[Point, Cyclo] = {}
    for direction in group:
        direction_totals[direction] = cyclo_scale(
            direction_averages[direction],
            N,
        )

    complete_cayley_numerator = cyclo_sum(
        [
            direction_totals[direction]
            for direction in group
            if direction != (0, 0)
        ]
    )
    complete_cayley_denominator = N * 338
    require(
        complete_cayley_denominator == 57122,
        "complete-Cayley denominator changed",
    )
    complete_cayley_energy = cyclo_scale(
        complete_cayley_numerator,
        Fraction(1, complete_cayley_denominator),
    )
    require(
        complete_cayley_energy
        == cyclo_rational(expected_cross_energy),
        "complete-Cayley charged overlap identity changed",
    )

    line_1_numerator = cyclo_sum(
        [direction_totals[(step, 0)] for step in range(P)]
    )
    line_2_numerator = cyclo_sum(
        [direction_totals[(0, step)] for step in range(P)]
    )
    line_denominator = P * N
    line_1 = cyclo_as_rational(
        cyclo_scale(line_1_numerator, Fraction(1, line_denominator)),
        "first coordinate-line energy stopped being rational",
    )
    line_2 = cyclo_as_rational(
        cyclo_scale(line_2_numerator, Fraction(1, line_denominator)),
        "second coordinate-line energy stopped being rational",
    )
    require(
        line_1 == 50 and line_2 == 218,
        "coordinate-line charged overlap identity changed",
    )

    energy_10 = expected_cross_energy - line_2 / 2
    energy_01 = expected_cross_energy - line_1 / 2
    energy_11 = (
        (line_1 + line_2) / 2 - expected_cross_energy
    )
    require(
        (energy_10, energy_01, energy_11)
        == (Fraction(16), Fraction(100), Fraction(9))
        and line_1 == 2 * (energy_10 + energy_11)
        and line_2 == 2 * (energy_01 + energy_11),
        "coordinate-mask reconstruction changed",
    )

    # One real quadrature at t=0 and one at t=1 recover a Gaussian-rational
    # current.  Work in the exact quadratic extension Q(zeta_13, i):
    #
    # U_0 = 2*x,
    # U_1 = zeta^-1*(x+i*y) + zeta*(x-i*y)
    #     = cos(2*pi/13)*U_0 + 2*sin(2*pi/13)*y.
    probe: Gaussian = (Fraction(3, 5), Fraction(-7, 4))
    probe_value = cycloi_gaussian(probe)
    probe_conjugate = cycloi_conjugate_i(probe_value)
    zeta = cycloi_embed(root(1))
    zeta_inverse = cycloi_embed(root(-1))
    u_zero = cycloi_add(probe_value, probe_conjugate)
    u_one = cycloi_add(
        cycloi_multiply(zeta_inverse, probe_value),
        cycloi_multiply(zeta, probe_conjugate),
    )
    cosine = cycloi_embed(
        cyclo_scale(cyclo_add(root(1), root(-1)), Fraction(1, 2))
    )
    sine = (
        cyclo_zero(),
        cyclo_scale(
            cyclo_subtract(root(1), root(-1)),
            Fraction(-1, 2),
        ),
    )
    reconstruction_residual = cycloi_subtract(
        u_one,
        cycloi_multiply(cosine, u_zero),
    )
    require(
        u_zero
        == cycloi_gaussian((2 * probe[0], ZERO))
        and reconstruction_residual
        == cycloi_scale(sine, 2 * probe[1])
        and sine != cycloi_zero(),
        "one-quadrature reconstruction identity changed",
    )

    # Ordinary affine translations recover the charged product on a
    # nonzero target line.  Here b=3 and z has three Gaussian-rational
    # entries.  The real translated-union bank is
    #
    # U(s,t)=A_s*zeta^(-tb)+conjugate(A_s*zeta^(-tb)).
    #
    # Its two-variable DFT is z(b,r)+conjugate(z(-b,-r)).
    affine_charge = 3

    # Four raw frequency-side controls independently check the relative
    # shift sign in U(s,t).  This real fixture has z=d*conjugate(w)
    # supported on the same affine line b=3.
    affine_real_deletion = {target: ZERO for target in group}
    affine_real_retained = {target: ZERO for target in group}
    for residue, value in (
        (0, Fraction(1)),
        (4, Fraction(-3)),
        (9, Fraction(2)),
    ):
        affine_real_deletion[(affine_charge, residue)] = value
        affine_real_retained[(affine_charge, residue)] = ONE
    affine_real_dft = fourier_transform(group, affine_real_deletion)
    affine_real_wft = fourier_transform(group, affine_real_retained)
    affine_real_norms = sum(
        value * value
        for value in affine_real_deletion.values()
    ) + sum(
        value * value
        for value in affine_real_retained.values()
    )
    affine_raw_controls = 0
    for shift_s, shift_t in ((0, 0), (1, 0), (0, 1), (4, 7)):
        raw_energy = cyclo_zero()
        for frequency in group:
            shifted_d = affine_real_dft[
                add_group(frequency, (0, shift_s))
            ]
            shifted_w = affine_real_wft[
                add_group(frequency, (shift_t, 0))
            ]
            raw_energy = cyclo_add(
                raw_energy,
                cyclo_norm_squared(cyclo_add(shifted_d, shifted_w)),
            )
        raw_union = cyclo_subtract(
            cyclo_scale(raw_energy, Fraction(1, N)),
            cyclo_rational(affine_real_norms),
        )
        line_value = cyclo_zero()
        for residue in range(P):
            line_value = cyclo_add(
                line_value,
                cyclo_scale(
                    root(shift_s * residue),
                    affine_real_deletion[(affine_charge, residue)]
                    * affine_real_retained[(affine_charge, residue)],
                ),
            )
        oriented = cyclo_phase(
            line_value,
            -shift_t * affine_charge,
        )
        require(
            raw_union
            == cyclo_add(oriented, cyclo_conjugate(oriented)),
            "raw affine translated-union sign or normalization changed",
        )
        affine_raw_controls += 1

    affine_values = {
        target: cycloi_zero()
        for target in group
    }
    affine_values[(affine_charge, 0)] = cycloi_gaussian(
        (Fraction(1), Fraction(2))
    )
    affine_values[(affine_charge, 4)] = cycloi_gaussian(
        (Fraction(-3), Fraction(1))
    )
    affine_values[(affine_charge, 9)] = cycloi_gaussian(
        (Fraction(2), Fraction(-5))
    )

    affine_line_transform: dict[int, CycloI] = {}
    affine_unions: dict[tuple[int, int], CycloI] = {}
    for shift_s in range(P):
        line_value = cycloi_zero()
        for residue in range(P):
            line_value = cycloi_add(
                line_value,
                cycloi_phase(
                    affine_values[(affine_charge, residue)],
                    shift_s * residue,
                ),
            )
        affine_line_transform[shift_s] = line_value
        for shift_t in range(P):
            oriented = cycloi_phase(
                line_value,
                -shift_t * affine_charge,
            )
            affine_unions[(shift_s, shift_t)] = cycloi_add(
                oriented,
                cycloi_conjugate(oriented),
            )

    affine_projector_checks = 0
    for target_charge, target_residue in group:
        raw = cycloi_zero()
        for shift_s in range(P):
            for shift_t in range(P):
                raw = cycloi_add(
                    raw,
                    cycloi_phase(
                        affine_unions[(shift_s, shift_t)],
                        -shift_s * target_residue
                        + shift_t * target_charge,
                    ),
                )
        recovered = cycloi_scale(raw, Fraction(1, N))
        reflected = cycloi_conjugate(
            affine_values[
                ((-target_charge) % P, (-target_residue) % P)
            ]
        )
        require(
            recovered
            == cycloi_add(
                affine_values[(target_charge, target_residue)],
                reflected,
            ),
            "affine symmetrized target projector changed",
        )
        affine_projector_checks += 1

    # Only t=0 and t=b^-1 are needed on a nonzero line.  Cross-multiply
    # the 2x2 inversion to avoid any numerical division by zeta^-1-zeta,
    # then invert the remaining one-dimensional DFT.
    inverse_charge = pow(affine_charge, -1, P)
    determinant = cyclo_subtract(root(-1), root(1))
    require(
        inverse_charge == 9 and determinant != cyclo_zero(),
        "affine two-shift determinant fixture changed",
    )
    scaled_line_transform: dict[int, CycloI] = {}
    two_shift_checks = 0
    for shift_s in range(P):
        left = cycloi_multiply(
            cycloi_embed(determinant),
            affine_line_transform[shift_s],
        )
        right = cycloi_subtract(
            affine_unions[(shift_s, inverse_charge)],
            cycloi_phase(affine_unions[(shift_s, 0)], 1),
        )
        require(
            left == right,
            "two ordinary shifts stopped determining the affine line",
        )
        scaled_line_transform[shift_s] = right
        two_shift_checks += 1

    affine_line_inverse_checks = 0
    for residue in range(P):
        raw = cycloi_zero()
        for shift_s in range(P):
            raw = cycloi_add(
                raw,
                cycloi_phase(
                    scaled_line_transform[shift_s],
                    -shift_s * residue,
                ),
            )
        require(
            cycloi_scale(raw, Fraction(1, P))
            == cycloi_multiply(
                cycloi_embed(determinant),
                affine_values[(affine_charge, residue)],
            ),
            "affine line inverse DFT changed",
        )
        affine_line_inverse_checks += 1

    # The zero line is exactly I+J.  A J-odd pair is invisible to every
    # ordinary affine union, whereas a real-even J-even pair is doubled.
    zero_hostile = {
        target: cycloi_zero()
        for target in group
    }
    zero_hostile[(0, 2)] = cycloi_gaussian((ZERO, ONE))
    zero_hostile[(0, -2 % P)] = cycloi_gaussian((ZERO, ONE))
    zero_hostile_unions = []
    for shift_s in range(P):
        line_value = cycloi_sum(
            [
                cycloi_phase(
                    zero_hostile[(0, residue)],
                    shift_s * residue,
                )
                for residue in range(P)
            ]
        )
        zero_hostile_unions.append(
            cycloi_add(line_value, cycloi_conjugate(line_value))
        )
    require(
        all(value == cycloi_zero() for value in zero_hostile_unions),
        "zero-line I+J kernel hostile stopped being dark",
    )

    zero_even = {
        target: cycloi_zero()
        for target in group
    }
    zero_even[(0, 2)] = cycloi_gaussian((-ONE, ZERO))
    zero_even[(0, -2 % P)] = cycloi_gaussian((-ONE, ZERO))
    for residue in (2, -2 % P):
        symmetrized = cycloi_add(
            zero_even[(0, residue)],
            cycloi_conjugate(zero_even[(0, -residue % P)]),
        )
        require(
            symmetrized == cycloi_gaussian((Fraction(-2), ZERO)),
            "zero-line J-even positive control changed",
        )

    # Hostile 1: deletion and retained target energies can both survive while
    # their charged overlap is identically zero.
    hostile_deletion = {target: ZERO for target in group}
    hostile_retained = {target: ZERO for target in group}
    hostile_deletion[(1, 0)] = ONE
    hostile_retained[(0, 0)] = ONE
    hostile_deletion_fourier = fourier_transform(
        group,
        hostile_deletion,
    )
    hostile_retained_fourier = fourier_transform(
        group,
        hostile_retained,
    )
    hostile_correlation = cross_correlation_from_kernel(
        group,
        hostile_deletion,
        hostile_retained,
        orthogonality,
    )
    require(
        cross_correlation_at_shift(
            group,
            (4, 7),
            hostile_deletion_fourier,
            hostile_retained_fourier,
        )
        == cyclo_zero(),
        "raw disjoint-support hostile correlation changed",
    )
    hostile_deletion_energy = sum(
        value * value
        for target, value in hostile_deletion.items()
        if target != (0, 0)
    )
    hostile_retained_energy = sum(
        value * value
        for target, value in hostile_retained.items()
        if target != (0, 0)
    )
    hostile_overlap_energy = sum(
        (
            hostile_deletion[target]
            * hostile_retained[target]
        )
        ** 2
        for target in group
        if target != (0, 0)
    )
    require(
        hostile_deletion_energy == 1
        and hostile_retained_energy == 0
        and hostile_overlap_energy == 0
        and all(
            value == cyclo_zero()
            for value in hostile_correlation.values()
        ),
        "disjoint-support hostile stopped being dark",
    )

    # Hostile 2: distinct nonzero targets make both singleton spectra live,
    # but their charged product and the entire cross-correlation stay dark.
    split_deletion = {target: ZERO for target in group}
    split_retained = {target: ZERO for target in group}
    split_deletion[(1, 0)] = ONE
    split_retained[(0, 1)] = ONE
    split_deletion_fourier = fourier_transform(group, split_deletion)
    split_retained_fourier = fourier_transform(group, split_retained)
    split_correlation = cross_correlation_from_kernel(
        group,
        split_deletion,
        split_retained,
        orthogonality,
    )
    split_deletion_energy = sum(
        value * value
        for target, value in split_deletion.items()
        if target != (0, 0)
    )
    split_retained_energy = sum(
        value * value
        for target, value in split_retained.items()
        if target != (0, 0)
    )
    split_overlap_energy = sum(
        (split_deletion[target] * split_retained[target]) ** 2
        for target in group
        if target != (0, 0)
    )
    require(
        split_deletion_energy == 1
        and split_retained_energy == 1
        and split_overlap_energy == 0
        and all(
            value == cyclo_zero()
            for value in split_correlation.values()
        )
        and cross_correlation_at_shift(
            group,
            (8, 2),
            split_deletion_fourier,
            split_retained_fourier,
        )
        == cyclo_zero(),
        "distinct-nonzero disjoint-support hostile stopped being dark",
    )

    # Hostile 3: singleton magnitudes and one ordinary union magnitude do not
    # determine the oriented charged phase.
    singleton: Gaussian = (ONE, ZERO)
    positive_i: Gaussian = (ZERO, ONE)
    negative_i: Gaussian = (ZERO, -ONE)
    union_positive = gaussian_norm_squared(
        gaussian_add(singleton, positive_i)
    )
    union_negative = gaussian_norm_squared(
        gaussian_add(singleton, negative_i)
    )
    oriented_positive = gaussian_multiply(
        singleton,
        gaussian_conjugate(positive_i),
    )
    oriented_negative = gaussian_multiply(
        singleton,
        gaussian_conjugate(negative_i),
    )
    require(
        gaussian_norm_squared(singleton) == 1
        and gaussian_norm_squared(positive_i) == 1
        and gaussian_norm_squared(negative_i) == 1
        and union_positive == 2
        and union_negative == 2
        and oriented_positive == (ZERO, -ONE)
        and oriented_negative == (ZERO, ONE),
        "ordinary-union phase ambiguity changed",
    )

    print("THM-2380 CROSS-WORD CHARGED CORRELATION - exact audit")
    print(f"group: p={P} N={N}")
    print(
        "cross-correlation/inverse target cells: "
        f"{cross_correlation_checks}/{inverse_target_checks}"
    )
    print(f"direct frequency cross-correlation controls: {direct_cross_checks}")
    print(f"pair-twist orbit cells: {pair_twist_checks}")
    print(f"direct pair-twist Parseval controls: {direct_parseval_checks}")
    print(
        "raw target-projector cells/denominator: "
        f"{raw_projector_checks}/{target_projector_denominator}"
    )
    print(
        "complete-Cayley denominator/shared energy: "
        f"{complete_cayley_denominator}/{expected_cross_energy}"
    )
    print(f"direct direction Parseval controls: {direct_direction_checks}")
    print(f"coordinate-line energies: L1={line_1} L2={line_2}")
    print(
        "coordinate-mask energies: "
        f"E10={energy_10} E01={energy_01} E11={energy_11}"
    )
    print("one-quadrature reconstruction determinant: nonzero")
    print(
        "affine raw/projectors/two-shift/line inverse: "
        f"{affine_raw_controls}/{affine_projector_checks}/{two_shift_checks}/"
        f"{affine_line_inverse_checks}"
    )
    print(
        "zero-line I+J boundary: "
        "J-odd dark, J-even doubled"
    )
    print(
        "disjoint-support hostile: "
        f"deletion={hostile_deletion_energy} "
        f"retained={hostile_retained_energy} "
        f"overlap={hostile_overlap_energy}"
    )
    print(
        "distinct-nonzero hostile: "
        f"deletion={split_deletion_energy} "
        f"retained={split_retained_energy} "
        f"overlap={split_overlap_energy}"
    )
    print(
        "ordinary-union ambiguity: "
        "singletons=(1,1) union=2 oriented=(-i,+i)"
    )
    print("PASS")


if __name__ == "__main__":
    main()
