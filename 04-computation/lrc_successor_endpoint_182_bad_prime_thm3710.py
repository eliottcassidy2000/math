#!/usr/bin/env python3
"""Finite-exact THM-3672 endpoint compression and bad-prime audit.

This companion reconstructs all thirteen marked sets and all thirty lawful
successor chart endpoint ledgers.  It then separates the exact-denominator
182 stratum and audits the characteristic-13 radical filtration on the full
two-coset 182 skeleton.
"""

from __future__ import annotations

from collections import defaultdict
from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
import importlib.util
from math import comb
from pathlib import Path
import sys


sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation/lrc_successor_mass_all_pair_swap_control_thm3672.py"
P = 13

PINS = {
    "04-computation/lrc_successor_mass_all_pair_swap_control_thm3672.py":
        "a191f934b20494b98e878fc5504a328c47d4cb6a92ea332246f782fac80c01c8",
    "05-knowledge/results/lrc_successor_mass_all_pair_swap_control_thm3672.out":
        "0d72ff662c908d9d9934976edce67895135c17a1eb0ba409e0abc382563b52c6",
    "04-computation/lrc14_169_twist_two_twist_referee_thm2334.py":
        "0e4a9e181263647e13d2a6738b6996c45df901d9d2b37d4d589dfddfbdd91480",
    "04-computation/lrc_successor_bulk_endpoint_leak_thm3705.py":
        "82faa434621b965a4d714b932c078698ec0bff00171f43d9172e25615c1534cf",
    "05-knowledge/results/lrc_successor_bulk_endpoint_leak_thm3705.out":
        "035538645e5184ad0bce3625a1d286611968ecd0c901bcd5b85b43e99cc3033b",
}


def require(condition, payload):
    if condition is not True:
        raise RuntimeError(payload)


def load_source():
    spec = importlib.util.spec_from_file_location("thm3672_source", SOURCE)
    require(spec is not None and spec.loader is not None, "cannot load THM-3672")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def add_endpoints(counter, intervals, coefficient, period):
    for left, right in intervals:
        counter[left % period] -= coefficient
        counter[right % period] += coefficient


def primitive(residue, period):
    """Periodic primitive P of d-1/7, normalized by P(0)=0."""
    window = period // 14
    require(period == 14 * window, period)
    if residue <= window:
        return Fraction(6 * residue, 7 * period)
    if residue <= period - window:
        return Fraction(1, 14) - Fraction(residue, 7 * period)
    return Fraction(6 * (residue - period), 7 * period)


def endpoint_residual(counter, period, successor_speed):
    return sum(
        coefficient * primitive(residue, period)
        for residue, coefficient in counter.items()
        if coefficient
    ) / successor_speed


def fractional_part(value):
    return value - value.numerator // value.denominator


def scalar_distance_squared(vector):
    mean = sum(vector, Fraction()) / len(vector)
    return sum((value - mean) ** 2 for value in vector)


def kgrid_fiber(counter, period, sign):
    """Coefficients on {(14*j+sign)/182:j in F_13}."""
    require(period % 182 == 0, period)
    scale = period // 182
    return tuple(
        counter.get(((14 * index + sign) % 182) * scale, 0)
        for index in range(P)
    )


def radical_order(coefficients):
    """(z-1)-adic order in F_13[z]/((z-1)^13), by Hasse jets."""
    require(len(coefficients) == P, len(coefficients))
    for degree in range(P):
        hasse_value = sum(
            coefficient * comb(index, degree)
            for index, coefficient in enumerate(coefficients)
        ) % P
        if hasse_value:
            return degree
    return P


def decimal_sqrt(value):
    return (Decimal(value.numerator) / Decimal(value.denominator)).sqrt()


def main():
    for relative, expected in PINS.items():
        actual = sha256((ROOT / relative).read_bytes()).hexdigest()
        require(actual == expected, ("parent hash drift", relative, actual))

    source = load_source()
    ref = source.load_referee()
    require(ref.W == (1, 14, 27, 40, 53, 66, 13, 2197, 742586), ref.W)
    require(ref.R == 169 and ref.W[ref.TARGET_B] == 2 * 13**5, (ref.R, ref.W))
    successor_speed = P * ref.W[ref.TARGET_B]
    period = ref.NN // successor_speed
    require((successor_speed, period) == (9653618, 5214048840), (successor_speed, period))

    # A boundary of the far blocker inside the delayed word satisfies
    # w_B*u=q+/-1/14 and u=R*x.  Its successor phase is therefore
    # {13*w_B*x}={(14q+/-1)/182}.  This gives two full 13-cosets.
    far_boundary_phases = set()
    for residue in range(P):
        for sign in (-1, 1):
            phase = fractional_part(
                Fraction(P, ref.R) * (residue + Fraction(sign, 14))
            )
            expected = fractional_part(Fraction(14 * residue + sign, 182))
            require(phase == expected, (residue, sign, phase, expected))
            far_boundary_phases.add(phase)
    require(len(far_boundary_phases) == 26, far_boundary_phases)
    require(
        sorted(value.denominator for value in far_boundary_phases).count(182) == 24,
        far_boundary_phases,
    )
    require(
        {value for value in far_boundary_phases if value.denominator == 14}
        == {Fraction(1, 14), Fraction(13, 14)},
        far_boundary_phases,
    )

    zero = (0,) * 9
    word = ref.build_boolean_set(ref.PATTERN_QA, zero)
    interval_sets = {
        "zero": source.marked_intervals(
            ref,
            ref.build_boolean_set(ref.PATTERN_E, zero),
            word,
        )
    }
    for graft in source.UNITS:
        interval_sets[f"a{graft}"] = source.marked_intervals(
            ref,
            ref.build_boolean_set(
                ref.PATTERN_E,
                source.negative_target_dipole(ref.TARGET_A, graft),
            ),
            word,
        )
        interval_sets[f"b{graft}"] = source.marked_intervals(
            ref,
            ref.build_boolean_set(
                ref.PATTERN_E,
                source.negative_target_dipole(ref.TARGET_B, graft),
            ),
            word,
        )

    labels = (
        ("zero",)
        + tuple(f"a{k}" for k in source.UNITS)
        + tuple(f"b{k}" for k in source.UNITS)
    )
    expected_ledger = dict(source.EXPECTED_VALUES)
    counters = {}
    residuals = {}
    exact182_residuals = {}
    for label in labels:
        intervals = interval_sets[label]
        mass = sum(right - left for left, right in intervals)
        window = ref.NN // (14 * successor_speed)
        danger = sum(
            source.danger_prefix(right, period, window)
            - source.danger_prefix(left, period, window)
            for left, right in intervals
        )
        successor = 2 * mass - danger
        require(
            (successor, mass, danger, len(intervals)) == expected_ledger[label],
            (label, successor, mass, danger, len(intervals)),
        )
        counter = defaultdict(int)
        add_endpoints(counter, intervals, 1, period)
        counter = {residue: coefficient for residue, coefficient in counter.items() if coefficient}
        endpoint = endpoint_residual(counter, period, successor_speed)
        direct = (Fraction(danger) - Fraction(mass, 7)) / ref.NN
        require(endpoint == direct, (label, endpoint, direct))
        counters[label] = counter
        residuals[label] = endpoint
        exact182_residuals[label] = sum(
            coefficient * primitive(residue, period)
            for residue, coefficient in counter.items()
            if Fraction(residue, period).denominator == 182
        ) / successor_speed

    expected_primitive_numerators = {
        (14 * index + sign) % 182
        for index in range(P)
        for sign in (-1, 1)
    } - {13, 169}
    require(len(expected_primitive_numerators) == 24, expected_primitive_numerators)

    records = []
    primitive_profiles = []
    full_profiles = []
    exact_only_radical_orders = []
    for k in source.UNITS:
        for ell in source.UNITS:
            if k == ell:
                continue
            counter = defaultdict(int)
            for residue, coefficient in counters["zero"].items():
                counter[residue] += coefficient
            for residue, coefficient in counters[f"a{k}"].items():
                counter[residue] += coefficient
            for residue, coefficient in counters[f"b{ell}"].items():
                counter[residue] -= 2 * coefficient
            counter = {residue: coefficient for residue, coefficient in counter.items() if coefficient}

            endpoint = endpoint_residual(counter, period, successor_speed)
            expected_endpoint = residuals["zero"] + residuals[f"a{k}"] - 2 * residuals[f"b{ell}"]
            require(endpoint == expected_endpoint, (k, ell, endpoint, expected_endpoint))
            exact182 = sum(
                coefficient * primitive(residue, period)
                for residue, coefficient in counter.items()
                if Fraction(residue, period).denominator == 182
            ) / successor_speed
            nonexact182 = endpoint - exact182
            primitive_numerators = {
                Fraction(residue, period).numerator
                for residue in counter
                if Fraction(residue, period).denominator == 182
            }
            require(primitive_numerators == expected_primitive_numerators, (k, ell, primitive_numerators))

            plus = kgrid_fiber(counter, period, 1)
            minus = kgrid_fiber(counter, period, -1)
            require(all(coefficient != 0 for coefficient in plus + minus), (k, ell, plus, minus))
            require((radical_order(plus), radical_order(minus)) == (0, 0), (k, ell, plus, minus))
            # The exact-denominator-182 fibers delete the denominator-14
            # corners: index 12 in the plus fiber and index 1 in the minus
            # fiber.  Keep this attribution separate from the full grid.
            exact_plus = tuple(0 if index == 12 else value for index, value in enumerate(plus))
            exact_minus = tuple(0 if index == 1 else value for index, value in enumerate(minus))
            exact_only_radical_orders.extend(
                (
                    (k, ell, "plus", radical_order(exact_plus)),
                    (k, ell, "minus", radical_order(exact_minus)),
                )
            )
            primitive_profiles.append(tuple(sorted(primitive_numerators)))
            full_profiles.append(
                (
                    k,
                    ell,
                    tuple(value % P for value in plus),
                    tuple(value % P for value in minus),
                )
            )

            effective_l1 = sum(abs(coefficient) for coefficient in counter.values())
            naive_l1 = 2 * (
                len(interval_sets["zero"])
                + len(interval_sets[f"a{k}"])
                + 2 * len(interval_sets[f"b{ell}"])
            )
            records.append(
                (k, ell, len(counter), effective_l1, naive_l1, endpoint, exact182, nonexact182)
            )

    require(len(records) == 30, len(records))
    require(len(set(primitive_profiles)) == 1, len(set(primitive_profiles)))
    require(len({row[2:] for row in full_profiles}) == 13, len({row[2:] for row in full_profiles}))
    exact_only_order_histogram = {
        order: sum(1 for *_, value in exact_only_radical_orders if value == order)
        for order in sorted({value for *_, value in exact_only_radical_orders})
    }
    exact_only_exceptions = tuple(
        row for row in exact_only_radical_orders if row[3] != 0
    )
    require(exact_only_order_histogram == {0: 57, 2: 3}, exact_only_order_histogram)
    require(
        exact_only_exceptions
        == ((4, 0, "minus", 2), (4, 1, "minus", 2), (4, 2, "minus", 2)),
        exact_only_exceptions,
    )
    require(all(record[5] > 0 for record in records), "endpoint sign")
    require(all(record[6] > 0 for record in records), "exact-182 sign")
    require(all(record[7] < 0 for record in records), "nonexact-182 sign")
    total_range = (min(record[5] for record in records), max(record[5] for record in records))
    exact_range = (min(record[6] for record in records), max(record[6] for record in records))
    rest_range = (min(record[7] for record in records), max(record[7] for record in records))
    require(
        total_range
        == (
            Fraction(83117797033, 9034385901100560),
            Fraction(12230170063, 1129298237637570),
        ),
        total_range,
    )
    require(
        exact_range
        == (
            Fraction(78691, 878479238),
            Fraction(598097, 6149354666),
        ),
        exact_range,
    )
    require(
        rest_range
        == (
            Fraction(-61198645499, 694952761623120),
            Fraction(-3420318877, 43434547601445),
        ),
        rest_range,
    )

    canonical = next(record for record in records if record[:2] == (1, 2))
    require(canonical[2:5] == (42, 59408, 1461670), canonical)
    require(canonical[5] == Fraction(11109404003, 1129298237637570), canonical[5])
    require(canonical[6] == Fraction(586695, 6149354666), canonical[6])
    require(canonical[7] == Fraction(-3716699972, 43434547601445), canonical[7])
    effective_bound = Fraction(3 * canonical[3], 49 * successor_speed)
    naive_bound = Fraction(3 * canonical[4], 49 * successor_speed)
    require(effective_bound == Fraction(89112, 236513641), effective_bound)

    _, _, canonical_plus, canonical_minus = next(
        row for row in full_profiles if row[:2] == (1, 2)
    )
    require(
        canonical_plus == (6, 0, 0, 0, 4, 10, 10, 12, 5, 12, 12, 12, 12),
        canonical_plus,
    )
    require(
        canonical_minus == (7, 6, 0, 0, 9, 3, 3, 1, 8, 1, 1, 1, 7),
        canonical_minus,
    )

    # In characteristic thirteen, z^13-1=(z-1)^13.  A uniform 13-fiber is
    # the top radical vector (z-1)^12.  Removing either lower-period corner
    # leaves augmentation 12=-1, so the primitive 12-point fiber has order 0.
    uniform = (1,) * P
    primitive_plus = tuple(0 if index == 12 else 1 for index in range(P))
    primitive_minus = tuple(0 if index == 1 else 1 for index in range(P))
    require(radical_order(uniform) == 12, "uniform radical order")
    require(
        (radical_order(primitive_plus), radical_order(primitive_minus)) == (0, 0),
        (primitive_plus, primitive_minus),
    )

    full_vector = tuple(residuals[label] for label in labels)
    exact_vector = tuple(exact182_residuals[label] for label in labels)
    rest_vector = tuple(full - exact for full, exact in zip(full_vector, exact_vector))
    require(full_vector == tuple(exact + rest for exact, rest in zip(exact_vector, rest_vector)), "vector split")
    full_distance_squared = scalar_distance_squared(full_vector)
    exact_distance_squared = scalar_distance_squared(exact_vector)
    rest_distance_squared = scalar_distance_squared(rest_vector)
    require(
        full_distance_squared
        == Fraction(31351115173959958275463, 530530835965029753707225010038400),
        full_distance_squared,
    )
    require(
        exact_distance_squared == Fraction(2004390109503, 245794658253663815114),
        exact_distance_squared,
    )
    require(
        rest_distance_squared
        == Fraction(21868283125702350430567, 3139235715769406826669970473600),
        rest_distance_squared,
    )

    getcontext().prec = 50
    full_distance = decimal_sqrt(full_distance_squared)
    exact_distance = decimal_sqrt(exact_distance_squared)
    rest_distance = decimal_sqrt(rest_distance_squared)
    inner = (full_distance_squared - exact_distance_squared - rest_distance_squared) / 2
    require(
        inner == Fraction(-418382806563352927, 55555643095377343171212960),
        inner,
    )
    cosine = (
        Decimal(inner.numerator) / Decimal(inner.denominator)
    ) / (exact_distance * rest_distance)

    print("chart_count=30")
    print(f"phase_support_range={min(row[2] for row in records)}..{max(row[2] for row in records)}")
    print(f"effective_endpoint_l1_range={min(row[3] for row in records)}..{max(row[3] for row in records)}")
    print(f"naive_endpoint_l1_range={min(row[4] for row in records)}..{max(row[4] for row in records)}")
    print(f"canonical_support_l1_naive={(canonical[2], canonical[3], canonical[4])}")
    print(f"canonical_endpoint={canonical[5]}")
    print(f"canonical_exact182={canonical[6]}")
    print(f"canonical_nonexact182={canonical[7]}")
    print(f"canonical_effective_bound_over_actual={Decimal(effective_bound.numerator) / Decimal(effective_bound.denominator) / (Decimal(canonical[5].numerator) / Decimal(canonical[5].denominator))}")
    print(f"canonical_naive_bound_over_actual={Decimal(naive_bound.numerator) / Decimal(naive_bound.denominator) / (Decimal(canonical[5].numerator) / Decimal(canonical[5].denominator))}")
    print("full_182_grid_skeleton=26")
    print("primitive_exact_denominator_182_stratum=24")
    print(f"primitive_182_numerators={sorted(expected_primitive_numerators)}")
    print(f"actual_plus_fiber_mod13={canonical_plus}")
    print(f"actual_minus_fiber_mod13={canonical_minus}")
    print("actual_coefficient_profile_count=13")
    print("uniform_13_fiber_radical_order=12")
    print("primitive_corner_deleted_radical_orders=(0,0)")
    print("actual_full_grid_all_chart_radical_orders=(0,0)")
    print(f"actual_exact182_radical_order_histogram={exact_only_order_histogram}")
    print(f"actual_exact182_positive_depth_fibers={exact_only_exceptions}")
    print("exact182_all_positive=True")
    print("nonexact182_all_negative=True")
    print(f"endpoint_range={total_range[0]}..{total_range[1]}")
    print(f"exact182_range={exact_range[0]}..{exact_range[1]}")
    print(f"nonexact182_range={rest_range[0]}..{rest_range[1]}")
    print(f"endpoint_scalar_distance_squared={full_distance_squared}")
    print(f"exact182_scalar_distance_squared={exact_distance_squared}")
    print(f"nonexact182_scalar_distance_squared={rest_distance_squared}")
    print(f"endpoint_scalar_distance={full_distance}")
    print(f"exact182_scalar_distance={exact_distance}")
    print(f"nonexact182_scalar_distance={rest_distance}")
    print(f"exact182_nonexact_scalar_inner={inner}")
    print(f"exact182_nonexact_scalar_cosine={cosine}")
    print("PASS")


if __name__ == "__main__":
    main()
