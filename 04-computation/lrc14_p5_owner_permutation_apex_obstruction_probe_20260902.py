#!/usr/bin/env python3
"""Exact audit of the unit-tail owner formula and one p=5 apex escape.

For 2 <= d <= 7 and a d-unit tail w over quotient phase theta, its
strict 1/14-danger mask on the d lifts is empty or a singleton.  When
||w theta|| < d/14, the singleton is

    -round(w theta) * w^(-1)  (mod d).

The symbolic theorem and its proof live in the matching reflection.  This
script exhaustively checks the formula, including equality and d=7
half-integer boundaries, checks the pair-tooth intersection formula, and
freezes one h=420, r_5=5 row.  It is a certificate audit, not LRC(14).
All load-bearing arithmetic is integer or fractions.Fraction arithmetic.
"""

from __future__ import annotations

import ast
from fractions import Fraction as Q
from itertools import combinations
from math import floor, gcd, isqrt
from pathlib import Path


DELTA = Q(1, 14)
CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def distance(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Q) -> Q:
    return min(distance(speed * phase) for speed in speeds)


def nearest_integer_or_tie(value: Q) -> int | None:
    lower = floor(value)
    offset = value - lower
    if offset < Q(1, 2):
        return lower
    if offset > Q(1, 2):
        return lower + 1
    return None


def literal_bad_sheets(d: int, speed: int, theta: Q) -> tuple[int, ...]:
    return tuple(
        sheet
        for sheet in range(d)
        if distance(Q(speed) * (theta + sheet) / d) < DELTA
    )


def predicted_bad_sheets(d: int, speed: int, theta: Q) -> tuple[int, ...]:
    require(2 <= d <= 7, ("d-range", d))
    require(gcd(d, speed) == 1, ("unit", d, speed))
    if distance(speed * theta) >= Q(d, 14):
        return ()
    nearest = nearest_integer_or_tie(speed * theta)
    require(nearest is not None, ("active-nearest-unique", d, speed, theta))
    owner = (-nearest * pow(speed, -1, d)) % d
    return (owner,)


def owner_data(d: int, speed: int, theta: Q) -> tuple[int, Q, int]:
    value = speed * theta
    nearest = nearest_integer_or_tie(value)
    require(nearest is not None, ("owner-nearest", d, speed, theta))
    error = abs(value - nearest)
    require(error < Q(d, 14), ("owner-active", d, speed, theta, error))
    owner = (-nearest * pow(speed, -1, d)) % d
    return nearest, error, owner


def tooth_intersection_length(
    d: int, u: int, n_u: int, v: int, n_v: int
) -> tuple[Q, Q]:
    centre_u = Q(n_u, u)
    centre_v = Q(n_v, v)
    radius_u = Q(d, 14 * u)
    radius_v = Q(d, 14 * v)
    direct = max(
        Q(0),
        min(centre_u + radius_u, centre_v + radius_v)
        - max(centre_u - radius_u, centre_v - radius_v),
    )
    determinant = abs(n_u * v - n_v * u)
    formula = max(
        Q(0),
        min(
            2 * radius_u,
            2 * radius_v,
            radius_u + radius_v - Q(determinant, u * v),
        ),
    )
    return direct, formula


def factorization(value: int) -> tuple[tuple[int, int], ...]:
    result: list[tuple[int, int]] = []
    divisor = 2
    while divisor * divisor <= value:
        exponent = 0
        while value % divisor == 0:
            value //= divisor
            exponent += 1
        if exponent:
            result.append((divisor, exponent))
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        result.append((value, 1))
    return tuple(result)


def product_divisors(value: int) -> set[int]:
    result = {1}
    for prime, exponent in factorization(value):
        result = {
            divisor * prime**power
            for divisor in result
            for power in range(exponent + 1)
        }
    return result


def trial_divisors(value: int) -> set[int]:
    result = {1, value}
    for divisor in range(2, isqrt(value) + 1):
        if value % divisor == 0:
            result.add(divisor)
            result.add(value // divisor)
    return result


def capacity(row: tuple[int, ...], scale: int) -> Q:
    total = Q(0)
    for speed in row:
        if speed % scale:
            order = scale // gcd(scale, speed)
            total += Q((order + 6) // 7, order)
    return total


def exact_pair_maximum(a: int, b: int) -> Q:
    candidates = {Q(0)}
    for speed in (a, b):
        candidates.update(Q(k, 2 * speed) for k in range(2 * speed))
    for denominator in {a + b, abs(a - b)} - {0}:
        candidates.update(Q(k, denominator) for k in range(denominator))
    return max(min(distance(a * x), distance(b * x)) for x in candidates)


def exact_formula_audit() -> tuple[int, int, int]:
    formula_checks = 0
    gauge_checks = 0
    pair_checks = 0

    for d in range(2, 8):
        for speed in range(1, 46):
            if gcd(d, speed) != 1:
                continue

            # Both exact activity boundaries are safe.  At d=7 these are
            # half-integers, so nearest-integer ownership is deliberately
            # undefined while the strict mask is still empty.
            for theta in (Q(d, 14 * speed), Q(14 - d, 14 * speed)):
                require(
                    distance(speed * theta) == Q(d, 14),
                    ("activity-equality", d, speed, theta),
                )
                require(
                    literal_bad_sheets(d, speed, theta) == (),
                    ("equality-empty-mask", d, speed, theta),
                )
                require(
                    predicted_bad_sheets(d, speed, theta) == (),
                    ("equality-prediction", d, speed, theta),
                )
                if d == 7:
                    require(
                        nearest_integer_or_tie(speed * theta) is None,
                        ("half-integer-tie", speed, theta),
                    )

            for denominator in range(2, 30):
                for numerator in range(denominator):
                    theta = Q(numerator, denominator)
                    literal = literal_bad_sheets(d, speed, theta)
                    predicted = predicted_bad_sheets(d, speed, theta)
                    require(literal == predicted, ("owner-formula", d, speed, theta))
                    require(len(literal) <= 1, ("singleton-cap", d, speed, theta))
                    formula_checks += 1
                    if literal:
                        shifted = literal_bad_sheets(d, speed, theta + 1)
                        require(
                            shifted == tuple((sheet - 1) % d for sheet in literal),
                            ("gauge-shift", d, speed, theta),
                        )
                        gauge_checks += 1

    # A second exact path audits the interval formula and determinant residue.
    # It uses literal common active points, then rebuilds both open teeth from
    # their independently rounded integer addresses.
    for d in range(2, 8):
        units = tuple(speed for speed in range(1, 24) if gcd(speed, d) == 1)
        for u, v in combinations(units, 2):
            for denominator in range(3, 19):
                for numerator in range(denominator):
                    theta = Q(numerator, denominator)
                    mask_u = literal_bad_sheets(d, u, theta)
                    mask_v = literal_bad_sheets(d, v, theta)
                    if not mask_u or not mask_v or mask_u == mask_v:
                        continue
                    n_u, _, owner_u = owner_data(d, u, theta)
                    n_v, _, owner_v = owner_data(d, v, theta)
                    require(owner_u != owner_v, ("owner-distinct", d, u, v, theta))
                    determinant = n_u * v - n_v * u
                    require(determinant % d != 0, ("determinant-residue", d, u, v, theta))
                    require(
                        determinant % gcd(u, v) == 0,
                        ("determinant-gcd", d, u, v, theta),
                    )
                    direct, formula = tooth_intersection_length(d, u, n_u, v, n_v)
                    require(direct == formula > 0, ("intersection-formula", d, u, v, theta))
                    pair_checks += 1

    return formula_checks, gauge_checks, pair_checks


def main() -> None:
    source = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(source)), "no-assert")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(source)
        ),
        "no-float",
    )

    formula_checks, gauge_checks, pair_checks = exact_formula_audit()

    d = 5
    h = 420
    anchor = 2 * h
    core = (13, 168, 349, 375, 711, 737, 1073, 1099)
    tails = (1, 21, 23, 327, 689)
    row = tuple(d * speed for speed in core) + tails
    odd_relatives = tuple(speed for speed in row if speed != anchor)

    require(len(row) == len(set(row)) == 13, "row-cardinality")
    require(anchor in row, "anchor")
    require(len(odd_relatives) == 12 and all(speed % 2 for speed in odd_relatives), "odd-W")
    require(gcd(*row) == 1, "primitive")
    require(h % 420 == 0, "420-wall")
    require(sum(speed % 5 != 0 for speed in odd_relatives) == 5, "r5-equality")
    require(
        tuple(speed // 5 for speed in row if speed % 5 == 0) == core,
        "p5-core",
    )

    theta_zero = Q(7, 181)
    mu = Q(90, 181)
    core_max = max(core)
    rho = (mu - DELTA) / core_max
    require(exact_pair_maximum(13, 168) == mu, "pair-upper-bound")
    require(clearance(core, theta_zero) == mu, "core-attainment")

    apex = tuple((speed,) + owner_data(d, speed, theta_zero) for speed in tails)
    owners = tuple(record[3] for record in apex)
    require(owners == (0, 4, 3, 1, 2), ("apex-owners", owners))
    require(set(owners) == set(range(d)), "apex-permutation")
    moments = tuple(sum(owner**power for owner in owners) % d for power in range(1, d))
    require(moments == (0, 0, 0, 4), ("f5-moments", moments))

    apex_slacks = tuple(
        (speed, Q(d, 14) - error - speed * rho)
        for speed, _, error, _ in apex
    )
    require(dict(apex_slacks)[327] < 0, "tail-327-crosses-core-margin")
    require(dict(apex_slacks)[689] < 0, "tail-689-crosses-core-margin")
    generic_height = Q(d * core_max, 14 * (mu - DELTA))
    require(max(tails) < generic_height, "old-safe-band-does-not-fire")

    theta_endpoint = Q(59, 1526)
    require(theta_zero - theta_endpoint == Q(3, 276206), "endpoint-displacement")
    require(theta_zero - theta_endpoint <= rho, "endpoint-inside-core-margin")
    require(clearance(core, theta_endpoint) == Q(741, 1526), "endpoint-core-clearance")
    endpoint_masks = tuple(
        (speed, literal_bad_sheets(d, speed, theta_endpoint)) for speed in tails
    )
    require(
        endpoint_masks == ((1, (0,)), (21, (4,)), (23, (3,)), (327, ()), (689, ())),
        ("endpoint-masks", endpoint_masks),
    )
    witness_sheet = 2
    witness = (theta_endpoint + witness_sheet) / d
    witness_clearance = clearance(row, witness)
    require(witness == Q(3111, 7630), "witness-phase")
    require(witness_clearance == Q(551, 7630) > DELTA, "witness-clearance")

    missing_denominators = tuple(
        modulus
        for modulus in range(2, 15)
        if not any(speed % modulus == 0 for speed in row)
    )
    require(not missing_denominators, "thm366-complete")

    bank_rows: list[tuple[int, bool, tuple[int, ...]]] = []
    for p in range(8, 15):
        modulus = 2 * p
        represented = {
            min(speed % modulus, modulus - speed % modulus)
            for speed in odd_relatives
            if gcd(speed, modulus) == 1
        }
        target = {
            min(unit, modulus - unit)
            for unit in range(1, modulus)
            if gcd(unit, modulus) == 1
        }
        anchor_killed = h % p == 0
        if not anchor_killed:
            require(represented == target, ("unit-bank-complete", p, represented, target))
        bank_rows.append((p, anchor_killed, tuple(sorted(represented))))

    divisors_trial = {
        divisor
        for speed in row
        for divisor in trial_divisors(speed)
        if divisor >= 2
    }
    divisors_product = {
        divisor
        for speed in row
        for divisor in product_divisors(speed)
        if divisor >= 2
    }
    require(divisors_trial == divisors_product, "divisor-paths")
    capacities = {scale: capacity(row, scale) for scale in sorted(divisors_trial)}
    capacity_minimum = min(capacities.values())
    capacity_minimizers = tuple(
        scale for scale, value in capacities.items() if value == capacity_minimum
    )
    capacity_closures = tuple(scale for scale, value in capacities.items() if value < 1)
    require(len(capacities) == 70, "capacity-scale-count")
    require(capacity_minimum == 1 and capacity_minimizers == (5,), "capacity-minimum")
    require(not capacity_closures, "capacity-has-no-closure")

    half_turn_modulus = 28 * h
    half_turn_phases = (
        Q(1, 2) - Q(1, half_turn_modulus),
        Q(1, 2) + Q(1, half_turn_modulus),
    )
    half_turn_clearances = tuple(clearance(row, phase) for phase in half_turn_phases)
    require(half_turn_clearances == (Q(11, 336), Q(11, 336)), "half-turn-killed")

    grid_survivors = 0
    for shift in range(half_turn_modulus):
        phase = Q(1, 2) + Q(shift, half_turn_modulus)
        if clearance(row, phase) >= DELTA:
            grid_survivors += 1
    require(grid_survivors == 1416, "full-grid-scope-control")

    print("LRC14 P5 OWNER-PERMUTATION APEX OBSTRUCTION -- EXACT AUDIT")
    print("STATUS=SYMBOLIC_LEMMA_AUDIT_PLUS_ONE_VERIFIED_EXACT_ROW_NOT_LRC14")
    print(
        f"FORMULA_CHECKS={formula_checks} GAUGE_CHECKS={gauge_checks} "
        f"PAIR_INTERSECTION_CHECKS={pair_checks} TOTAL_REQUIRE_CALLS={CHECKS}"
    )
    print("BOUNDARY=activity equality gives empty strict mask; d=7 half-integer ownership left undefined")
    print(f"h={h} d={d} anchor={anchor} primitive=1 r5=5 missing_thm366={missing_denominators}")
    print("core=" + ",".join(map(str, core)))
    print("tails=" + ",".join(map(str, tails)))
    print(f"core_maximum={mu} theta0={theta_zero} R={core_max} rho={rho}")
    print("apex=" + ";".join(f"w{w}:n{n}:err{error}:owner{owner}" for w, n, error, owner in apex))
    print(f"owner_moments_mod5={moments} exact_permutation={int(set(owners) == set(range(d)))}")
    print("apex_slacks=" + ";".join(f"w{w}:{slack}" for w, slack in apex_slacks))
    print(f"generic_MCb_threshold={generic_height} max_tail={max(tails)} fires=0")
    for p, anchor_killed, represented in bank_rows:
        print(f"unit_bank_p={p} anchor_killed={int(anchor_killed)} represented={represented}")
    print(
        f"adaptive_scales={len(capacities)} minimum={capacity_minimum} "
        f"minimizers={capacity_minimizers} closures={capacity_closures}"
    )
    print(f"half_turn_clearances={half_turn_clearances} full_grid_survivors={grid_survivors}")
    print(f"endpoint_theta={theta_endpoint} endpoint_masks={endpoint_masks}")
    print(
        f"witness_sheet={witness_sheet} witness={witness} "
        f"clearance={witness_clearance} margin={witness_clearance - DELTA}"
    )
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
