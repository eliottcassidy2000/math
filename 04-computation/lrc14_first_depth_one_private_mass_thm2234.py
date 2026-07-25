#!/usr/bin/env python3
"""Exact referee for THM-2234's depth-one private two-owner remainder.

The proof itself is analytic.  This companion checks:

* the odd-guard charge rows used in the Hunter lower bound;
* the specialized worst row when the sixth tooth is divisible by 13;
* a single compatible coefficient packet attaining every independent
  charge debit in the ledger;
* the three exact mass constants; and
* local sharpness witnesses for the 10-root guard cap and the subsequent
  12-root complement-of-a-unit-comb cap.

All arithmetic is integer or ``Fraction`` exact.  ``require`` remains active
under ``python -O``.
"""

from fractions import Fraction as F
from math import gcd


P = 13


def require(condition, message):
    """Raise on a failed exact check, including under optimized Python."""
    if not condition:
        raise RuntimeError(message)


def reduced_charge(a, b):
    """THM-2080's centered odd-guard charge for a reduced ratio (a,b)."""
    require(a > 0 and b > 0 and gcd(a, b) == 1, "charge row is not reduced")
    require(a % 2 == 1, "guard quotient is not odd")
    x = F(a % 14, 14)
    y = F(b % 7, 7)
    correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
    return F(2, a * b) * correction


def top_negative_rows(product_cutoff, predicate):
    """Enumerate negative charge rows below an analytic product cutoff."""
    rows = []
    for a in range(1, product_cutoff + 1, 2):
        for b in range(1, product_cutoff // a + 1):
            if gcd(a, b) != 1 or not predicate(a, b):
                continue
            charge = reduced_charge(a, b)
            if charge < 0:
                rows.append((-charge, a, b))
    rows.sort(reverse=True)
    return tuple(rows)


def norm_numerator(value, modulus):
    """Return modulus times the circle norm of value/modulus."""
    residue = value % modulus
    return min(residue, modulus - residue)


def in_danger(multiplier, numerator, modulus):
    """Strict membership in D_multiplier on a rational torsion grid."""
    return 14 * norm_numerator(multiplier * numerator, modulus) < modulus


def in_guard(multiplier, numerator, modulus):
    """Strict membership in C_multiplier on a rational torsion grid."""
    return 7 * norm_numerator(multiplier * numerator, modulus) > modulus


def first_fibre_cap_witness():
    """Return the exact ten-root private-remainder fibre witness."""
    modulus = P**3
    phase_modulus = P**2
    phase = 1
    guard = 1
    units = (1, 2, 3, 4, 5)
    shallow_unit = 14
    shallow = P * shallow_unit
    private_roots = []
    guard_roots = []
    unit_root_sheets = []
    for sheet in range(P):
        numerator = phase + sheet * phase_modulus
        if in_guard(guard, numerator, modulus):
            guard_roots.append(sheet)
        if (
            in_guard(guard, numerator, modulus)
            and all(not in_danger(q, numerator, modulus) for q in units)
            and not in_danger(shallow, numerator, modulus)
        ):
            private_roots.append(sheet)
    for q in units:
        unit_root_sheets.append(
            tuple(
                sheet
                for sheet in range(P)
                if in_danger(q, phase + sheet * phase_modulus, modulus)
            )
        )
    require(
        not in_danger(shallow_unit, phase, phase_modulus),
        "the shallow divided tooth is not safe at the hostile phase",
    )
    require(tuple(guard_roots) == tuple(range(2, 12)), "guard fibre changed")
    require(
        tuple(unit_root_sheets) == ((0,),) * 5,
        "five distinct unit masks no longer share the hostile singleton",
    )
    require(tuple(private_roots) == tuple(range(2, 12)), "ten-cap witness lost")
    return phase, tuple(private_roots), tuple(unit_root_sheets)


def second_fibre_cap_witness():
    """Return a torsion witness attaining 12 roots in T(P)."""
    modulus = P**3
    image_modulus = P**2
    base_modulus = P
    guard = 1
    units = (1, 2, 3, 4, 5)
    shallow = P

    private = {
        numerator
        for numerator in range(modulus)
        if in_guard(guard, numerator, modulus)
        and all(not in_danger(q, numerator, modulus) for q in units)
        and not in_danger(shallow, numerator, modulus)
    }
    image = {numerator % image_modulus for numerator in private}
    roots_of_zero = tuple(k * base_modulus for k in range(P))
    occupied = tuple(y for y in roots_of_zero if y in image)
    require(
        occupied == tuple(k * base_modulus for k in range(1, P)),
        "twelve-cap image witness changed",
    )

    first_preimages = []
    for y in occupied:
        witnesses = tuple(
            y + k * image_modulus
            for k in range(P)
            if y + k * image_modulus in private
        )
        require(witnesses, f"image point {y} has no private preimage")
        first_preimages.append(witnesses[0])
    return occupied, tuple(first_preimages)


def main():
    # The fifth unrestricted debit is 4/441.  THM-2080's analytic
    # -e <= 1/(4ab) makes every row with ab >= 28 strictly smaller, so the
    # complete finite search for the top five is ab <= 27.
    unrestricted = top_negative_rows(
        27, lambda a, b: a % P != 0 and b % P != 0
    )
    expected_unit_top = (
        (F(5, 294), 1, 6),
        (F(8, 539), 11, 1),
        (F(3, 245), 3, 5),
        (F(3, 245), 1, 5),
        (F(4, 441), 9, 2),
    )
    require(unrestricted[:5] == expected_unit_top, "top-five unit ledger drift")
    require(F(1, 4 * 28) < F(4, 441), "unit analytic tail does not separate")

    # If 13 divides the terminal coefficient and not the guard, its reduced
    # denominator b is divisible by 13.  The candidate 5/637 dominates the
    # analytic tail ab >= 32, leaving only the rows (1,13) and (1,26).
    deep_rows = top_negative_rows(
        31, lambda a, b: a % P != 0 and b % P == 0
    )
    expected_deep_rows = (
        (F(5, 637), 1, 13),
        (F(3, 1274), 1, 26),
    )
    require(deep_rows == expected_deep_rows, "13-divisible charge census drift")
    require(F(1, 4 * 32) < F(5, 637), "deep analytic tail does not separate")

    # These six reduced rows are simultaneously compatible with one fixed
    # guard H=99.  Hence the independent-charge debit cannot be lowered.
    guard = 99
    unit_coefficients = (594, 9, 165, 495, 22)
    shallow = 1287
    expected_rows = (
        (1, 6),
        (11, 1),
        (3, 5),
        (1, 5),
        (9, 2),
        (1, 13),
    )
    actual_rows = []
    for coefficient in unit_coefficients + (shallow,):
        common = gcd(guard, coefficient)
        actual_rows.append((guard // common, coefficient // common))
    require(tuple(actual_rows) == expected_rows, "compatible charge packet drift")
    require(guard % 2 == 1 and guard % P != 0, "hostile guard type drift")
    require(
        all(q % P != 0 for q in unit_coefficients),
        "hostile unit coefficient became 13-divisible",
    )
    require(shallow % P == 0 and shallow % P**2 != 0, "hostile depth drift")
    require(len(set(unit_coefficients)) == 5, "hostile units are not distinct")

    unit_debit = sum((row[0] for row in expected_unit_top), F(0))
    deep_debit = F(5, 637)
    private_floor = F(5, 49) - unit_debit - deep_debit
    first_image_floor = F(13, 10) * private_floor
    second_image_floor = F(13, 12) * first_image_floor
    require(private_floor == F(2593, 90090), "private mass floor drift")
    require(first_image_floor == F(2593, 69300), "first image floor drift")
    require(second_image_floor == F(33709, 831600), "second image floor drift")

    first_phase, first_sheets, unit_sheets = first_fibre_cap_witness()
    second_roots, second_preimages = second_fibre_cap_witness()

    print("unit_charge_top5:")
    for debit, a, b in expected_unit_top:
        print(f"  (a,b)=({a},{b}) deficit={debit}")
    print("depth_one_13_divisible_charge_rows:")
    for debit, a, b in expected_deep_rows:
        print(f"  (a,b)=({a},{b}) deficit={debit}")
    print(f"compatible_guard={guard}")
    print(f"compatible_units={unit_coefficients}")
    print(f"compatible_depth_one={shallow}")
    print(f"private_mass_floor={private_floor}")
    print(f"first_image_floor={first_image_floor}")
    print(f"second_image_floor={second_image_floor}")
    print(
        f"ten_cap_witness_phase={first_phase}/169 "
        f"private_sheets={first_sheets} unit_sheets={unit_sheets}"
    )
    print(f"twelve_cap_witness_image_roots={second_roots}")
    print(f"twelve_cap_witness_first_preimages={second_preimages}")
    print("strict_13_power_endpoint_control=PASS")
    print("status=THM2234_PRIVATE_TWO_OWNER_MASS")


if __name__ == "__main__":
    main()
