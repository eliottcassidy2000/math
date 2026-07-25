#!/usr/bin/env python3
"""Exact referee for THM-2244's joint two-step private-fibre cap.

THM-2234 bounds the first private image with guard occupancy 10 and the
second with peeled-complement occupancy 12.  Those caps cannot saturate
independently.  On the intermediate thirteen-root fibre, the allowed
peeled-complement roots and the guard roots must overlap.  This script checks
the resulting four-bit table, every displayed rational constant, and an exact
169-root hostile control attaining the joint cap 112 while all five distinct
unit masks remain present.
"""

from fractions import Fraction


P = 13
PRIVATE_FLOOR = Fraction(2593, 90090)
OLD_SECOND_FLOOR = Fraction(33709, 831600)
EXPECTED_JOINT_CAP = 112
EXPECTED_SECOND_FLOOR = Fraction(33709, 776160)
EXPECTED_GAIN = Fraction(33709, 11642400)
EXPECTED_STOPPING_GAP = Fraction(77171, 776160)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle_norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def four_bit_table():
    rows = []
    for danger_at_z in (0, 1):
        for guard_at_z in (0, 1):
            allowed = 11 + danger_at_z
            guard = 10 - guard_at_z
            intersection_floor = allowed + guard - 13
            cap = 10 * allowed - intersection_floor
            rows.append(
                (
                    danger_at_z,
                    guard_at_z,
                    allowed,
                    guard,
                    intersection_floor,
                    cap,
                )
            )
    return tuple(rows)


def hostile_control():
    """Return the exact locally sharp 169-root configuration."""
    guard_speed = 1
    peeled_core = 2
    deep_owner_core = 2
    phase = Fraction(325007, 700000)
    unit_speeds = tuple(
        1 + P**2 * 700000 * index for index in range(1, 6)
    )

    require(guard_speed % 2 == 1, "hostile guard is not odd")
    require(guard_speed % P != 0, "hostile guard is not a thirteen-unit")
    require(peeled_core % P != 0, "peeled core is not a thirteen-unit")
    require(
        circle_norm(deep_owner_core * phase) < Fraction(1, 14),
        "hostile depth-two owner is not active at the terminal phase",
    )
    require(len(set(unit_speeds)) == 5, "unit hostile speeds are not distinct")
    require(
        all(speed > 0 and speed % P != 0 for speed in unit_speeds),
        "unit hostile speed scope drift",
    )
    require(
        all((speed - guard_speed) * phase % 1 == 0 for speed in unit_speeds),
        "hostile unit phases no longer align with the guard speed",
    )

    intermediate = tuple((phase + digit) / P for digit in range(P))
    allowed_intermediate = tuple(
        value
        for value in intermediate
        if circle_norm(peeled_core * value) >= Fraction(1, 14)
    )
    guard_intermediate = tuple(
        value
        for value in intermediate
        if circle_norm(guard_speed * value) > Fraction(1, 7)
    )
    intersection = set(allowed_intermediate).intersection(guard_intermediate)

    grandchildren = tuple(
        (phase + digit) / P**2 for digit in range(P**2)
    )
    eligible = []
    for value in grandchildren:
        guard = circle_norm(guard_speed * value) > Fraction(1, 7)
        peeled_complement = (
            circle_norm(P * peeled_core * value) >= Fraction(1, 14)
        )
        unit_complement = all(
            circle_norm(speed * value) >= Fraction(1, 14)
            for speed in unit_speeds
        )
        if guard and peeled_complement and unit_complement:
            require(
                circle_norm(P**2 * deep_owner_core * value)
                < Fraction(1, 14),
                "hostile eligible root lost its labelled deep owner",
            )
            eligible.append(value)

        # The five distinct masks are exactly aligned with D_H on this fibre.
        guard_danger = circle_norm(guard_speed * value) < Fraction(1, 14)
        for speed in unit_speeds:
            require(
                (circle_norm(speed * value) < Fraction(1, 14))
                == guard_danger,
                "hostile unit-mask alignment drift",
            )

        # Keep the witness off every load-bearing boundary.
        require(
            circle_norm(guard_speed * value) != Fraction(1, 7),
            "hostile guard boundary encountered",
        )
        require(
            circle_norm(P * peeled_core * value) != Fraction(1, 14),
            "hostile peeled boundary encountered",
        )
        require(
            all(
                circle_norm(speed * value) != Fraction(1, 14)
                for speed in unit_speeds
            ),
            "hostile unit boundary encountered",
        )

    return {
        "guard_speed": guard_speed,
        "peeled_core": peeled_core,
        "deep_owner_core": deep_owner_core,
        "phase": phase,
        "unit_speeds": unit_speeds,
        "allowed_count": len(allowed_intermediate),
        "guard_count": len(guard_intermediate),
        "intersection_count": len(intersection),
        "eligible_count": len(eligible),
    }


def main():
    rows = four_bit_table()
    expected_rows = (
        (0, 0, 11, 10, 8, 102),
        (0, 1, 11, 9, 7, 103),
        (1, 0, 12, 10, 9, 111),
        (1, 1, 12, 9, 8, 112),
    )
    require(rows == expected_rows, "four-bit joint-cap table drift")
    require(
        max(row[-1] for row in rows) == EXPECTED_JOINT_CAP,
        "joint cap drift",
    )

    second_floor = Fraction(P**2, EXPECTED_JOINT_CAP) * PRIVATE_FLOOR
    require(second_floor == EXPECTED_SECOND_FLOOR, "new floor drift")
    require(
        second_floor - OLD_SECOND_FLOOR == EXPECTED_GAIN,
        "old/new floor gain drift",
    )
    require(
        Fraction(1, 7) - second_floor == EXPECTED_STOPPING_GAP,
        "first exhausted-owner stopping gap drift",
    )
    require(
        Fraction(120, EXPECTED_JOINT_CAP) == Fraction(15, 14),
        "relative cap improvement drift",
    )

    witness = hostile_control()
    require(
        (
            witness["allowed_count"],
            witness["guard_count"],
            witness["intersection_count"],
            witness["eligible_count"],
        )
        == (12, 9, 8, EXPECTED_JOINT_CAP),
        "hostile sharpness control drift",
    )

    print("four_bit_table=(Da(z),CH(z),A,G,I_floor,joint_cap)")
    for row in rows:
        print(row)
    print(f"joint_two_step_cap={EXPECTED_JOINT_CAP}")
    print(
        f"private_floor={PRIVATE_FLOOR} "
        f"old_second_floor={OLD_SECOND_FLOOR}"
    )
    print(
        f"new_second_floor={second_floor} "
        f"decimal={float(second_floor):.15f}"
    )
    print(
        f"gain={EXPECTED_GAIN} relative_factor=15/14 "
        f"one_seventh_minus={EXPECTED_STOPPING_GAP}"
    )
    print(
        "hostile="
        f"H={witness['guard_speed']} a={witness['peeled_core']} "
        f"deep_owner={witness['deep_owner_core']} "
        f"z={witness['phase']} q={witness['unit_speeds']}"
    )
    print(
        "hostile_counts="
        f"A={witness['allowed_count']} G={witness['guard_count']} "
        f"I={witness['intersection_count']} "
        f"eligible169={witness['eligible_count']}"
    )


if __name__ == "__main__":
    main()
