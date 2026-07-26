#!/usr/bin/env python3
"""Exact referee for THM-2362's thirteen-shift successor statistic.

All interval boundaries and masses use ``Fraction``.  The script exhausts
the 26 open cells cut by the target-coordinate shifts, checks the distinct
inverse-root kernel, verifies every rational Fourier-floor constant, and
audits the pure/fork blocker truth table.  No floating point or Python
``assert`` is used.
"""

from __future__ import annotations

from fractions import Fraction


P = 13
HALF_DANGER = Fraction(1, 14)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fractional_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def centered_distance(x: Fraction) -> Fraction:
    value = fractional_part(x)
    return min(value, 1 - value)


def danger(x: Fraction) -> int:
    return int(centered_distance(x) < HALF_DANGER)


def target_shift_danger(y: Fraction, shift: int) -> int:
    return danger(y + Fraction(shift, P))


def inverse_root_danger(y: Fraction, shift: int) -> int:
    return danger((y + shift) / P)


def main() -> None:
    # The target-shift boundaries are also the boundaries of d(13y).
    boundaries = {
        fractional_part(sign * HALF_DANGER - Fraction(shift, P))
        for shift in range(P)
        for sign in (-1, 1)
    }
    require(len(boundaries) == 2 * P, "wrong target-shift boundary count")
    ordered = sorted(boundaries)
    samples = [
        fractional_part((ordered[index] + ordered[(index + 1) % len(ordered)]) / 2)
        if index + 1 < len(ordered)
        else (ordered[index] + ordered[0] + 1) / 2
        for index in range(len(ordered))
    ]
    samples[-1] = fractional_part(samples[-1])

    target_cells = 0
    inverse_cells = 0
    for y in samples:
        target_d = sum(target_shift_danger(y, shift) for shift in range(P))
        target_g = P - target_d
        require(
            target_d == 2 - danger(P * y),
            f"target danger count failed at y={y}",
        )
        require(
            target_g == 11 + danger(P * y),
            f"target complement count failed at y={y}",
        )
        target_cells += 1

        root_d = sum(inverse_root_danger(y, shift) for shift in range(P))
        root_g = P - root_d
        require(
            root_d == 2 - danger(y),
            f"inverse-root danger count failed at y={y}",
        )
        require(
            root_g == 11 + danger(y),
            f"inverse-root complement count failed at y={y}",
        )
        inverse_cells += 1

    # A literal rational witness separates the two operations.
    witness = Fraction(1, 100)
    require(danger(witness) == 1, "scale witness left the current danger")
    require(danger(P * witness) == 0, "scale witness entered successor danger")
    require(
        sum(target_shift_danger(witness, shift) for shift in range(P)) == 2,
        "target-shift hostile count changed",
    )
    require(
        sum(inverse_root_danger(witness, shift) for shift in range(P)) == 1,
        "inverse-root hostile count changed",
    )

    # Symbolic rational mass and Fourier-mode ledgers.  Sweep exact
    # successor overlaps rho_+/rho from 0 through 1.
    rho = Fraction(17, 31)
    successor_cases = 14
    for numerator in range(successor_cases):
        rho_plus = rho * Fraction(numerator, successor_cases - 1)

        danger_average = (2 * rho - rho_plus) / P
        danger_mode_sum = rho - danger_average
        require(
            danger_mode_sum == (11 * rho + rho_plus) / P,
            "danger nonzero-mode sum changed",
        )
        require(
            danger_mode_sum / (P - 1) >= Fraction(11, 156) * rho,
            "danger real-mode floor failed",
        )
        require(
            danger_mode_sum**2 / (P - 1)
            >= Fraction(121, 2028) * rho**2,
            "danger energy floor failed",
        )

        complement_average = (11 * rho + rho_plus) / P
        complement_mode_sum = rho - complement_average
        require(
            complement_mode_sum == (2 * rho - rho_plus) / P,
            "complement nonzero-mode sum changed",
        )
        require(
            complement_mode_sum / (P - 1) >= Fraction(1, 156) * rho,
            "complement real-mode floor failed",
        )
        require(
            complement_mode_sum**2 / (P - 1)
            >= Fraction(1, 2028) * rho**2,
            "complement energy floor failed",
        )

    # Inverse-root counts determine only the zero character.  A chosen
    # sheet anchor is extra lift data, and the nonzero-mode sum can have
    # the opposite sign from the formerly claimed constants.
    danger_lift = Fraction(99, 100)
    danger_root_values = [
        inverse_root_danger(danger_lift, shift) for shift in range(P)
    ]
    require(
        danger(danger_lift) == 1
        and sum(danger_root_values) == 1
        and danger_root_values[0] == 0
        and Fraction(danger_root_values[0], 1)
        - Fraction(sum(danger_root_values), P)
        == Fraction(-1, 13),
        "inverse-root danger lift hostile changed",
    )

    complement_lift = Fraction(1, 2)
    complement_root_values = [
        1 - inverse_root_danger(complement_lift, shift) for shift in range(P)
    ]
    require(
        danger(complement_lift) == 0
        and sum(complement_root_values) == 11
        and complement_root_values[0] == 0
        and Fraction(complement_root_values[0], 1)
        - Fraction(sum(complement_root_values), P)
        == Fraction(-11, 13),
        "inverse-root complement lift hostile changed",
    )

    # Scalar-cover truth table.  Exclude (0,0,0): at least one blocker is
    # dangerous.  Pure-a other factors g_j g_b force d_a=1; a fork with
    # d_b=1 does not, and deleting a complement from E_j is not redundant.
    cover_rows = [
        (d_j, d_a, d_b)
        for d_j in (0, 1)
        for d_a in (0, 1)
        for d_b in (0, 1)
        if d_j or d_a or d_b
    ]
    require(len(cover_rows) == 7, "wrong blocker truth-table size")
    pure_a_rows = [row for row in cover_rows if row[0] == 0 and row[2] == 0]
    require(
        pure_a_rows == [(0, 1, 0)],
        "pure-a danger factor stopped being redundant",
    )
    fork_without_a = [
        row for row in cover_rows if row[0] == 0 and row[2] == 1
    ]
    require(
        {row[1] for row in fork_without_a} == {0, 1},
        "fork falsely forced its a-danger factor",
    )
    base_without_g_a = [
        row for row in cover_rows if row[0] == 1 and row[2] == 0
    ]
    require(
        {row[1] for row in base_without_g_a} == {0, 1},
        "base word falsely made its a-complement redundant",
    )

    print("THM-2362 exact thirteen-shift successor referee")
    print(f"target/inverse open cells checked: {target_cells}/{inverse_cells}")
    print("target counts: sum d(y+s/13)=2-d(13y); sum g=11+d(13y)")
    print("inverse-root counts: sum d((y+s)/13)=2-d(y); sum g=11+d(y)")
    print("scale hostile y=1/100: target count 2, inverse-root count 1")
    print("inverse-root lift hostiles: danger -1/13; complement -11/13")
    print(
        "uniform role floors: danger Re>=11rho/156, "
        "energy>=121rho^2/2028"
    )
    print(
        "uniform role floors: complement Re>=rho/156, "
        "energy>=rho^2/2028"
    )
    print("blocker truth table: pure danger redundant; fork/complement not")
    print("VERDICT: the successor overlap rho_+ is a load-bearing jet sidecar")


if __name__ == "__main__":
    main()
