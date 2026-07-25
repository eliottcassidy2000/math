#!/usr/bin/env python3
"""Finite exact companion for the proposed THM-2286.

This script checks only the finite arithmetic supporting the intended
endpoint-Prony and owner-multiplier theorem:

* centered lift-window lengths and sharp absolute bounds modulo 13 and 7;
* representative exact Vandermonde determinants for consecutive windows;
* the homometric six-point power-spectrum counterexample;
* the proper-arc arithmetic for multipliers 1,...,6 and the hostile
  one-interval construction from multiplier 7 onward;
* the complete THM-2266 primitive-ratio atlas census; and
* the component and primitive-speed arithmetic inherited from
  THM-2199/THM-2203.

The general distributional-derivative identity, the generic Vandermonde
argument, continuity/support arguments for pushforwards, and their LRC
consequences are mathematical arguments in the theorem document.  Finite
enumeration here is deliberately not presented as proving those universal
claims.  Every validity gate uses ``require`` so normal and optimized
Python executions run the same checks.
"""

from collections import Counter
from fractions import Fraction
from math import gcd, prod
import sys


ATLAS_BOUND = 757
MODULUS_ROOT = 13
MODULUS_APEX = 7

HOMOMETRIC_A = (0, 1, 2, 6, 8, 11)
HOMOMETRIC_B = (0, 1, 6, 7, 9, 11)

CANON_SPEED_CEILING = int(
    "1423300902469616940457319158747393032068079495657519836405567800"
    "6806827357976205356500821273826531456252902899449503310895399282"
    "8270745007301563407410881277403096974669379856420053050310053068"
    "80000"
)


def require(condition: bool, message: str) -> None:
    """Raise in both normal and optimized Python modes."""
    if not condition:
        raise RuntimeError(message)


def fraction_determinant(matrix: list[list[Fraction]]) -> Fraction:
    """Exact determinant by fraction-preserving Gaussian elimination."""
    size = len(matrix)
    require(
        all(len(row) == size for row in matrix),
        "determinant matrix must be square",
    )
    work = [row[:] for row in matrix]
    determinant = Fraction(1)

    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if work[row][column] != 0
            ),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant

        pivot_value = work[column][column]
        determinant *= pivot_value
        for row in range(column + 1, size):
            factor = work[row][column] / pivot_value
            work[row][column] = Fraction(0)
            for entry in range(column + 1, size):
                work[row][entry] -= factor * work[column][entry]

    return determinant


def vandermonde_expected(
    nodes: tuple[Fraction, ...],
    first_exponent: int,
) -> Fraction:
    """Closed determinant for rows of consecutive powers."""
    scale = prod(
        (node**first_exponent for node in nodes),
        start=Fraction(1),
    )
    alternating = prod(
        (
            nodes[right] - nodes[left]
            for left in range(len(nodes))
            for right in range(left + 1, len(nodes))
        ),
        start=Fraction(1),
    )
    return scale * alternating


def centered_lift_window(
    modulus: int,
    residue: int,
    radius: int,
) -> tuple[int, ...]:
    """The 2R consecutive lifts k+m h with -R<=h<R."""
    require(modulus >= 2, "modulus must be at least two")
    require(1 <= residue < modulus, "residue must be nonzero")
    require(radius >= 1, "window radius must be positive")
    return tuple(
        residue + modulus * lift
        for lift in range(-radius, radius)
    )


def cyclic_one_runs(mask: int, cells: int) -> int:
    """Number of selected circular cell components in a bit mask."""
    require(cells >= 1, "cyclic cell count must be positive")
    full_mask = (1 << cells) - 1
    if mask == 0:
        return 0
    if mask == full_mask:
        return 1
    return sum(
        1
        for index in range(cells)
        if (mask >> index) & 1
        and not ((mask >> ((index - 1) % cells)) & 1)
    )


def centered_mod_one(value: Fraction) -> Fraction:
    """Return the representative in [-1/2,1/2)."""
    shifted = value + Fraction(1, 2)
    integer_part = shifted.numerator // shifted.denominator
    result = value - integer_part
    require(
        Fraction(-1, 2) <= result < Fraction(1, 2),
        "centered representative escaped its fundamental interval",
    )
    return result


def hostile_root_count(multiplier: int, target: Fraction) -> int:
    """Count roots in [-1/(2d),1/(2d)) for one non-boundary target."""
    half_width = Fraction(1, 2 * multiplier)
    count = 0
    for root_index in range(multiplier):
        root = centered_mod_one(
            Fraction(target + root_index, multiplier)
        )
        if -half_width <= root < half_width:
            count += 1
    return count


def normalized_translate(points: tuple[int, ...]) -> tuple[int, ...]:
    """Translation-normalized integer support."""
    minimum = min(points)
    return tuple(sorted(point - minimum for point in points))


def normalized_reflection(points: tuple[int, ...]) -> tuple[int, ...]:
    """Reflection- and translation-normalized integer support."""
    maximum = max(points)
    return tuple(sorted(maximum - point for point in points))


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    # A union of R proper circular intervals has at most 2R derivative
    # atoms.  The theorem combines coincident exponential nodes first.
    # Here we freeze the centered consecutive-window arithmetic.
    window_summary: dict[int, str] = {}
    for modulus in (MODULUS_ROOT, MODULUS_APEX):
        for radius in range(1, 65):
            union = []
            for residue in range(1, modulus):
                window = centered_lift_window(
                    modulus,
                    residue,
                    radius,
                )
                require(len(window) == 2 * radius, "lift-window length")
                require(
                    window[0] == residue - modulus * radius,
                    "lift-window lower endpoint",
                )
                require(
                    window[-1]
                    == residue + modulus * (radius - 1),
                    "lift-window upper endpoint",
                )
                require(
                    max(abs(value) for value in window)
                    <= modulus * radius - 1,
                    "centered lift-window absolute bound",
                )
                union.extend(window)
            require(
                max(abs(value) for value in union)
                == modulus * radius - 1,
                "centered lift-window bound must be sharp over residues",
            )

        # The all-radius inequalities reduce to these residue-only slacks:
        # (mR-1)-(mR-k)=k-1 and
        # (mR-1)-(mR-m+k)=m-k-1.
        require(
            all(
                residue - 1 >= 0
                and modulus - residue - 1 >= 0
                for residue in range(1, modulus)
            ),
            "symbolic centered-window slack",
        )
        window_summary[modulus] = f"{modulus}R-1"

    # Representative exact consecutive-power determinants.  These checks
    # verify matrix orientation, negative starting exponents, and the
    # 2R-window bookkeeping.  Generic nonvanishing for arbitrary distinct
    # complex nodes is the theorem-side Vandermonde argument.
    vandermonde_cases = 0
    for radius in range(1, 6):
        endpoint_count = 2 * radius
        nodes = tuple(
            Fraction(index + 1)
            for index in range(endpoint_count)
        )
        for first_exponent in (-radius, 0, radius + 1):
            matrix = [
                [
                    node ** (first_exponent + row)
                    for node in nodes
                ]
                for row in range(endpoint_count)
            ]
            exact = fraction_determinant(matrix)
            expected = vandermonde_expected(nodes, first_exponent)
            require(exact == expected, "Vandermonde determinant formula")
            require(exact != 0, "distinct-node Vandermonde must be nonzero")
            vandermonde_cases += 1

    # Exact homometry: equality of the oriented difference multiset is
    # precisely equality of every coefficient in the Laurent expansion of
    # the two pointwise power spectra.
    differences_a = Counter(
        left - right
        for left in HOMOMETRIC_A
        for right in HOMOMETRIC_A
    )
    differences_b = Counter(
        left - right
        for left in HOMOMETRIC_B
        for right in HOMOMETRIC_B
    )
    require(differences_a == differences_b, "homometric difference bank")
    require(sum(differences_a.values()) == 36, "difference multiplicity")
    require(differences_a[0] == 6, "zero-difference multiplicity")
    require(
        normalized_translate(HOMOMETRIC_A)
        != normalized_translate(HOMOMETRIC_B),
        "homometric supports must not be translates",
    )
    require(
        normalized_reflection(HOMOMETRIC_A)
        != normalized_translate(HOMOMETRIC_B),
        "homometric supports must not be reflected translates",
    )

    # For d<=6, the image dI of I=(-1/14,1/14) is a proper circular
    # interval, and its complement contains at least one full 1/7-period.
    proper_arc_table = {}
    for multiplier in range(1, 7):
        image_length = Fraction(multiplier, 7)
        complement_length = 1 - image_length
        require(image_length < 1, "owner image must be a proper arc")
        require(
            complement_length >= Fraction(1, 7),
            "owner image complement must contain a 1/7-period",
        )
        proper_arc_table[multiplier] = (
            image_length,
            complement_length,
        )

    # The sharp hostile interval is F_d=[-1/(2d),1/(2d)).  For d>=7 it
    # lies in I and multiplication by d maps its endpoints to the standard
    # fundamental interval [-1/2,1/2).  The frequency/length product at
    # frequency dn is exactly n, which is the finite arithmetic behind the
    # vanishing sine factor.
    hostile_sample_checks = 0
    for multiplier in range(7, 65):
        half_width = Fraction(1, 2 * multiplier)
        require(
            half_width <= Fraction(1, 14),
            "hostile interval must lie in the owner sector",
        )
        require(
            multiplier * (-half_width) == Fraction(-1, 2)
            and multiplier * half_width == Fraction(1, 2),
            "hostile interval endpoint mapping",
        )
        require(
            2 * half_width == Fraction(1, multiplier),
            "hostile interval length",
        )
        for frequency_index in range(1, 17):
            require(
                multiplier
                * frequency_index
                * Fraction(1, multiplier)
                == frequency_index,
                "hostile Fourier integer-period arithmetic",
            )

        sample_count = 2 * multiplier + 1
        for sample_index in range(sample_count):
            target = (
                Fraction(-1, 2)
                + Fraction(2 * sample_index + 1, 2 * sample_count)
            )
            require(
                hostile_root_count(multiplier, target) == 1,
                "hostile interval must select one root",
            )
            hostile_sample_checks += 1

    # The containment inequality is sharp at d=7.
    for multiplier in range(1, 129):
        require(
            Fraction(1, 2 * multiplier) <= Fraction(1, 14)
            if multiplier >= 7
            else Fraction(1, 2 * multiplier) > Fraction(1, 14),
            "sharp hostile containment boundary",
        )

    # Complete primitive-ratio atlas from THM-2266.
    oriented_atlas = tuple(
        (owner_multiplier, other_factor)
        for owner_multiplier in range(1, ATLAS_BOUND + 1)
        for other_factor in range(
            1,
            ATLAS_BOUND // owner_multiplier + 1,
        )
        if gcd(owner_multiplier, other_factor) == 1
    )
    unordered_atlas = tuple(
        pair
        for pair in oriented_atlas
        if pair[0] <= pair[1]
    )
    small_owner_atlas = tuple(
        pair
        for pair in oriented_atlas
        if pair[0] <= 6
    )
    small_owner_seven_unit_other = tuple(
        pair
        for pair in small_owner_atlas
        if pair[1] % 7 != 0
    )
    require(len(oriented_atlas) == 3_643, "oriented atlas count")
    require(len(unordered_atlas) == 1_822, "unordered atlas count")
    require(len(small_owner_atlas) == 1_372, "small-owner atlas count")
    require(
        len(small_owner_seven_unit_other) == 1_176,
        "small-owner seven-unit-other atlas count",
    )
    owner_histogram = tuple(
        sum(
            1
            for owner_multiplier, _ in small_owner_atlas
            if owner_multiplier == value
        )
        for value in range(1, 7)
    )
    seven_unit_other_histogram = tuple(
        sum(
            1
            for owner_multiplier, _ in small_owner_seven_unit_other
            if owner_multiplier == value
        )
        for value in range(1, 7)
    )
    require(
        owner_histogram == (757, 189, 168, 95, 121, 42),
        "small-owner atlas histogram",
    )
    require(
        seven_unit_other_histogram == (649, 162, 144, 81, 104, 36),
        "seven-unit-other atlas histogram",
    )

    # A Boolean combination cut by at most 2S distinct circular boundary
    # points has at most S components.  Exhaust the circular cell-run
    # arithmetic through S=8; the general cell argument is theorem-side.
    component_maxima = []
    for comb_order in range(1, 9):
        cells = 2 * comb_order
        maximum = max(
            cyclic_one_runs(mask, cells)
            for mask in range(1 << cells)
        )
        require(maximum == comb_order, "cyclic component maximum")
        component_maxima.append(maximum)

    # THM-2199 equation (23), reconstructed from its anisotropic row
    # heights.  THM-2203's fixed section is
    # (8H,16q_1,...,16q_5,16c_1,16c_2,16c_3).
    speed_ceiling = (
        12**6
        * 105**2
        * 178
        * 204
        * 262
        * 450
        * (78 * 7**21)
        * (78 * 182**13) ** 5
    )
    require(
        speed_ceiling == CANON_SPEED_CEILING,
        "THM-2199 anisotropic speed ceiling",
    )
    require(len(str(speed_ceiling)) == 197, "speed-ceiling digit count")
    require(speed_ceiling % 8 == 0, "speed ceiling divisibility by eight")

    scalar_comb_order_bound = 5 * speed_ceiling // 8
    endpoint_lift_abs_bound = (
        MODULUS_ROOT * scalar_comb_order_bound - 1
    )
    require(
        scalar_comb_order_bound
        == speed_ceiling // 8 + 8 * (speed_ceiling // 16),
        "THM-2203 scalar comb-order bound",
    )
    require(
        endpoint_lift_abs_bound
        == 13 * (5 * speed_ceiling // 8) - 1,
        "global endpoint-lift bound",
    )

    difference_signature = tuple(sorted(differences_a.items()))

    print("scope=FINITE-EXACT arithmetic; universal analysis remains theorem-side")
    print("endpoint_atoms_per_R<=2R")
    print(
        "centered_lift_windows="
        f"mod13:{window_summary[13]},mod7:{window_summary[7]}"
    )
    print(
        f"vandermonde_exact_cases={vandermonde_cases} "
        "radii=1..5 starts=(-R,0,R+1)"
    )
    print(f"homometric_A={HOMOMETRIC_A}")
    print(f"homometric_B={HOMOMETRIC_B}")
    print(f"homometric_difference_signature={difference_signature}")
    print("homometric_translate_or_reflection=False")
    print(
        "proper_owner_arc_table="
        + str(
            tuple(
                (
                    multiplier,
                    str(lengths[0]),
                    str(lengths[1]),
                )
                for multiplier, lengths in proper_arc_table.items()
            )
        )
    )
    print("hostile_multiplier_boundary=7")
    print(
        "hostile_interval=[-1/(2d),1/(2d)) "
        "scaled_interval=[-1/2,1/2)"
    )
    print(
        f"hostile_exact_sample_checks={hostile_sample_checks} "
        "d=7..64"
    )
    print(
        f"atlas_bound={ATLAS_BOUND} "
        f"oriented={len(oriented_atlas)} "
        f"unordered={len(unordered_atlas)}"
    )
    print(
        f"owner_multiplier_le_6={len(small_owner_atlas)} "
        f"histogram={owner_histogram}"
    )
    print(
        "owner_multiplier_le_6_and_other_7_unit="
        f"{len(small_owner_seven_unit_other)} "
        f"histogram={seven_unit_other_histogram}"
    )
    print(
        "cyclic_component_exhaustive="
        f"{tuple(component_maxima)} for S=1..8"
    )
    print(f"primitive_speed_ceiling={speed_ceiling}")
    print(f"primitive_speed_ceiling_digits={len(str(speed_ceiling))}")
    print(f"scalar_comb_order_bound={scalar_comb_order_bound}")
    print(f"endpoint_lift_abs_bound={endpoint_lift_abs_bound}")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
