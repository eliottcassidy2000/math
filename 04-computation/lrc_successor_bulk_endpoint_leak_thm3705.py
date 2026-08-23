#!/usr/bin/env python3
"""Exact companion for THM-3705's successor bulk/endpoint gate.

The theorem's general inequality is elementary.  This script independently
checks the centred danger primitive, finite interval identities, all exact
THM-3672 control data, and the resulting scalar-line calibration.
"""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys


sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
P = 13
CONTROL_DENOMINATOR = 50334435734703120
CONTROL_SUCCESSOR_SPEED = 9653618

PINS = {
    "04-computation/lrc_successor_mass_all_pair_swap_control_thm3672.py":
        "a191f934b20494b98e878fc5504a328c47d4cb6a92ea332246f782fac80c01c8",
    "05-knowledge/results/lrc_successor_mass_all_pair_swap_control_thm3672.out":
        "0d72ff662c908d9d9934976edce67895135c17a1eb0ba409e0abc382563b52c6",
    "04-computation/lrc_radial_successor_mass_star_frame_thm3701.py":
        "fadcdfa6ff2c4c1c4f48434ad6268c440fdacc8c3ea5fdb18f9de6a41ec69f0e",
    "05-knowledge/results/lrc_radial_successor_mass_star_frame_thm3701.out":
        "7b0043197502a6b38795ff4cab502f85acc4009bd5e04f4157bf8dc95b383d87",
}

# (successor numerator, marked numerator, danger-intersection numerator,
# component count), ordered S0,A0,...,A5,B0,...,B5.
CONTROL = (
    (110219232915792, 60084076348296, 9948919780800, 188056),
    ( 99299969997228, 54135630512964, 8971291028700, 169431),
    ( 99299969997228, 54135630512964, 8971291028700, 169431),
    ( 99299969997228, 54135630512964, 8971291028700, 169431),
    ( 99050735597814, 54000012243177, 8949288888540, 169006),
    ( 96962693739912, 52861450868226, 8760207996540, 165443),
    ( 96754696164126, 52747826144013, 8740956123900, 165088),
    (113880033524816, 61887542465528, 9895051406240, 186674),
    (113880033524816, 61887542465528, 9895051406240, 186674),
    (113880033524816, 61887542465528, 9895051406240, 186674),
    (113525486798282, 61695142630901, 9864798463520, 186093),
    (111210186545380, 60436521934750, 9662857324120, 182297),
    (111203438664916, 60432849886708, 9662261108500, 182285),
)


def require(condition, payload):
    if condition is not True:
        raise RuntimeError(payload)


def fractional_part(value: Fraction) -> Fraction:
    return value - value.numerator // value.denominator


def centred_danger_primitive(value: Fraction) -> Fraction:
    """The one-period primitive P of 1_(||x||<1/14)-1/7."""
    y = fractional_part(value)
    if y <= Fraction(1, 14):
        return Fraction(6, 7) * y
    if y <= Fraction(13, 14):
        return Fraction(1, 14) - y / 7
    return Fraction(6, 7) * (y - 1)


def intersection_length(left, right, other_left, other_right):
    return max(Fraction(), min(right, other_right) - max(left, other_left))


def brute_danger_integral(left: Fraction, right: Fraction, frequency: int) -> Fraction:
    """Directly intersect [left,right) with all danger windows in [0,1)."""
    total = Fraction()
    for integer in range(frequency):
        first_left = Fraction(integer, frequency)
        first_right = Fraction(14 * integer + 1, 14 * frequency)
        second_left = Fraction(14 * integer + 13, 14 * frequency)
        second_right = Fraction(integer + 1, frequency)
        total += intersection_length(left, right, first_left, first_right)
        total += intersection_length(left, right, second_left, second_right)
    return total


def endpoint_residual(intervals, frequency):
    return sum(
        centred_danger_primitive(frequency * right)
        - centred_danger_primitive(frequency * left)
        for left, right in intervals
    ) / frequency


def distance_to_scalars_squared(vector):
    mean = sum(vector, Fraction()) / len(vector)
    return sum((value - mean) ** 2 for value in vector)


def decimal_sqrt(value: Fraction) -> Decimal:
    return (Decimal(value.numerator) / Decimal(value.denominator)).sqrt()


def main():
    for relative, expected in PINS.items():
        actual = sha256((ROOT / relative).read_bytes()).hexdigest()
        require(actual == expected, ("parent hash drift", relative, actual))

    # Exact primitive range and an independent intersection calculation on
    # several nonaligned rational interval unions and frequencies.
    primitive_controls = 0
    for denominator in (14, 28, 91, 182, 273):
        for numerator in range(denominator + 1):
            value = centred_danger_primitive(Fraction(numerator, denominator))
            require(abs(value) <= Fraction(3, 49), (denominator, numerator, value))
            primitive_controls += 1
    require(centred_danger_primitive(Fraction(1, 14)) == Fraction(3, 49), "positive peak")
    require(centred_danger_primitive(Fraction(13, 14)) == Fraction(-3, 49), "negative peak")
    require(centred_danger_primitive(Fraction()) == centred_danger_primitive(Fraction(1)), "period")

    interval_banks = (
        ((Fraction(1, 31), Fraction(7, 29)),),
        ((Fraction(0), Fraction(3, 17)), (Fraction(8, 19), Fraction(17, 23))),
        ((Fraction(2, 21), Fraction(5, 18)), (Fraction(7, 13), Fraction(11, 14))),
    )
    interval_controls = 0
    for frequency in (1, 2, 5, 13, 26):
        for intervals in interval_banks:
            mass = sum(right - left for left, right in intervals)
            danger = sum(
                brute_danger_integral(left, right, frequency)
                for left, right in intervals
            )
            residual = endpoint_residual(intervals, frequency)
            require(danger == mass / 7 + residual, (frequency, intervals, danger, residual))
            require(abs(residual) <= Fraction(6 * len(intervals), 49 * frequency), "endpoint bound")
            successor = 2 * mass - danger
            require(successor == Fraction(13, 7) * mass - residual, "successor split")
            interval_controls += 1

    successor = tuple(Fraction(row[0], CONTROL_DENOMINATOR) for row in CONTROL)
    marked = tuple(Fraction(row[1], CONTROL_DENOMINATOR) for row in CONTROL)
    endpoint = tuple(
        (Fraction(row[2]) - Fraction(row[1], 7)) / CONTROL_DENOMINATOR
        for row in CONTROL
    )
    require(
        successor == tuple(Fraction(13, 7) * mass - leak for mass, leak in zip(marked, endpoint)),
        "control vector split",
    )

    marked_distance_squared = distance_to_scalars_squared(marked)
    endpoint_distance_squared = distance_to_scalars_squared(endpoint)
    successor_distance_squared = distance_to_scalars_squared(successor)
    require(
        marked_distance_squared
        == Fraction(
            7267086384084568099667227,
            97444439258883015987041328374400,
        ),
        marked_distance_squared,
    )
    require(
        endpoint_distance_squared
        == Fraction(
            31351115173959958275463,
            530530835965029753707225010038400,
        ),
        endpoint_distance_squared,
    )
    require(
        successor_distance_squared
        == Fraction(
            6394584339188742740200237,
            24361109814720753996760332093600,
        ),
        successor_distance_squared,
    )
    require(
        Fraction(169, 49) * marked_distance_squared > endpoint_distance_squared,
        "positive scalar-line gate",
    )

    # The canonical forward chart is retained as a sign check on the same
    # decomposition.  B2 is control index 9; A1 is index 2.
    chart = (CONTROL[0], CONTROL[2], CONTROL[9])
    marked_defect = chart[0][1] + chart[1][1] - 2 * chart[2][1]
    danger_defect = chart[0][2] + chart[1][2] - 2 * chart[2][2]
    successor_defect = chart[0][0] + chart[1][0] - 2 * chart[2][0]
    bulk_defect = Fraction(13 * marked_defect, 7)
    endpoint_defect = Fraction(danger_defect) - Fraction(marked_defect, 7)
    require(bulk_defect == Fraction(-124219914907348, 7), bulk_defect)
    require(endpoint_defect == Fraction(3466134048936, 7), endpoint_defect)
    require(successor_defect == -18240864136612, successor_defect)
    require(Fraction(successor_defect) == bulk_defect - endpoint_defect, "chart split")

    getcontext().prec = 60
    marked_distance = decimal_sqrt(marked_distance_squared)
    endpoint_distance = decimal_sqrt(endpoint_distance_squared)
    successor_distance = decimal_sqrt(successor_distance_squared)
    scaled_bulk_distance = Decimal(13) * marked_distance / Decimal(7)
    gamma = scaled_bulk_distance - endpoint_distance
    require(gamma > 0, gamma)
    require(successor_distance >= gamma, (successor_distance, gamma))
    drift_floor = gamma * gamma / Decimal(474552)
    target_floor = gamma * gamma / Decimal(6169176)

    component_bounds = tuple(
        Fraction(6 * row[3], 49 * CONTROL_SUCCESSOR_SPEED)
        for row in CONTROL
    )
    component_l2_bound_squared = sum(bound * bound for bound in component_bounds)
    component_l2_bound = decimal_sqrt(component_l2_bound_squared)

    print(f"primitive_grid_controls={primitive_controls}")
    print(f"independent_interval_controls={interval_controls}")
    print("primitive_range=[-3/49,3/49]")
    print("endpoint_bound=6*r/(49*C)")
    print("successor_vector=(13/7)*marked_vector-endpoint_vector")
    print(f"control_dist_marked_squared={marked_distance_squared}")
    print(f"control_dist_endpoint_squared={endpoint_distance_squared}")
    print(f"control_dist_successor_squared={successor_distance_squared}")
    print(f"control_dist_marked={marked_distance}")
    print(f"control_scaled_bulk_distance={scaled_bulk_distance}")
    print(f"control_dist_endpoint={endpoint_distance}")
    print(f"control_gamma={gamma}")
    print(f"control_dist_successor={successor_distance}")
    print(f"control_gamma_over_actual={gamma / successor_distance}")
    print(f"control_D_floor={drift_floor}")
    print(f"control_E_dt_floor={target_floor}")
    print(f"control_component_l2_absolute_bound={component_l2_bound}")
    print(f"control_absolute_component_bound_fires={scaled_bulk_distance > component_l2_bound}")
    print(f"canonical_bulk_defect_numerator={bulk_defect}")
    print(f"canonical_endpoint_defect_numerator={endpoint_defect}")
    print(f"canonical_successor_defect_numerator={successor_defect}")
    print("PASS")


if __name__ == "__main__":
    main()
