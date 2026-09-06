#!/usr/bin/env python3
"""Exact residual audit for 3*C plus arbitrary ternary-unit triples.

The parity-free width lemma reduces each fixed ten-body C to tails whose
largest entry c satisfies 7*c*L(C)<3.  This script exhausts that finite
residual for every C in binom([13],10).  The quotient lift-mask calculation
is checked against a literal merge of the danger intervals of the thirteen
physical speeds.  All arithmetic is Fraction arithmetic.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from multiprocessing import Pool
from lrc14_bounded_ten_body_geometry_thm4442 import (
    component_mask_pieces,
    exact_safe_geometry,
    first_safe_lift_cell,
    interval_measure,
    safe_at,
    safe_components_by_union,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def wall_witness(speeds: tuple[int, ...]) -> Q | None:
    """Independent endpoint fallback if the positive safe measure vanishes."""
    points = {Q(0)}
    for speed in speeds:
        for k in range(speed + 1):
            for sign in (-1, 1):
                point = Q(14 * k + sign, 14 * speed)
                if 0 <= point < 1:
                    points.add(point)
    return next((point for point in sorted(points) if safe_at(speeds, point)), None)


def audit_body(
    payload: tuple[tuple[int, ...], tuple[tuple[Q, Q], ...], Q, Q]
) -> tuple[
    int,
    int,
    int,
    int,
    int,
    int,
    dict[int, int],
    dict[str, int],
    tuple[Q, tuple[int, ...], tuple[int, int, int], int] | None,
    tuple[tuple[int, ...], tuple[int, int, int]] | None,
    tuple[int, tuple[int, ...], Q, int],
]:
    """One independently checkable body shard; safe for Windows spawn."""
    body, components, _measure, longest = payload
    cutoff_numerator = 3 * longest.denominator
    cutoff_denominator = 7 * longest.numerator
    exclusive_ceiling = (
        cutoff_numerator + cutoff_denominator - 1
    ) // cutoff_denominator
    values = tuple(
        value for value in range(1, exclusive_ceiling)
        if value % 3 != 0 and 7 * value * longest < 3
    )
    require(not values or max(values) < Q(3, 7) / longest, "residual bound")
    rows = positive = isolated = direct_count = piece_count = safe_piece_count = 0
    parity_counts: Counter[int] = Counter()
    circuit_counts: Counter[str] = Counter()
    minimum: tuple[Q, tuple[int, ...], tuple[int, int, int], int] | None = None
    failure: tuple[tuple[int, ...], tuple[int, int, int]] | None = None
    for tails in combinations(values, 3):
        rows += 1
        parity_mask = sum((tail & 1) << i for i, tail in enumerate(tails))
        parity_counts[parity_mask] += 1
        a, b, c = tails
        if a + b == c:
            circuit_counts["111"] += 1
        if c == 2 * a + b or c == a + 2 * b or 2 * b == a + c:
            circuit_counts["112"] += 1
        if (
            a + 2 * b == 2 * c
            or a + 2 * c == 2 * b
            or b + 2 * c == 2 * a
            or 2 * a + 2 * b == c
            or 2 * a + 2 * c == b
            or 2 * b + 2 * c == a
            or 2 * a == 2 * b + c
            or 2 * a == b + 2 * c
            or 2 * b == 2 * a + c
            or 2 * b == a + 2 * c
            or 2 * c == 2 * a + b
            or 2 * c == a + 2 * b
        ):
            circuit_counts["122"] += 1

        certificate = first_safe_lift_cell(components, tails)
        full_speeds = tuple(sorted(tuple(3 * value for value in body) + tails))
        require(len(full_speeds) == 13 and len(set(full_speeds)) == 13, "labels")
        pieces = component_mask_pieces(components, tails)
        mask_measure = sum(
            (right - left) * mask.bit_count() / 3
            for _, left, right, mask in pieces
        )
        piece_count += len(pieces)
        safe_piece_count += sum(mask != 0 for _, _, _, mask in pieces)

        direct_components = safe_components_by_union(full_speeds)
        direct_measure = interval_measure(direct_components)
        require(mask_measure == direct_measure, f"mass mismatch {body} {tails}")
        direct_count += len(direct_components)
        if direct_measure > 0:
            positive += 1
            require(certificate is not None, f"missing mask cell {body} {tails}")
            _, left, right, mask = certificate
            sheet = min(i for i in range(3) if mask & (1 << i))
            witness = ((left + right) / 2 + sheet) / 3
            require(safe_at(full_speeds, witness), f"bad witness {body} {tails}")
            record = (direct_measure, body, tails, len(direct_components))
            if minimum is None or record < minimum:
                minimum = record
        else:
            witness = wall_witness(full_speeds)
            if witness is not None:
                isolated += 1
            else:
                failure = (body, tails)
                break
    return (
        rows,
        positive,
        isolated,
        direct_count,
        piece_count,
        safe_piece_count,
        dict(parity_counts),
        dict(circuit_counts),
        minimum,
        failure,
        (rows, body, longest, max(values) if values else 0),
    )


def main() -> None:
    print("LRC14_TEN_BODY_ALL_PARITY_TERNARY_UNIT_RESIDUAL")
    print("STATUS=FINITE_EXACT_CANONICAL_AUDIT")
    print("body_universe=C in binom([13],10)")
    print("tail_universe=distinct positive T={a<b<c}, 3_not_dividing_abc")
    print("width_rule=7*max(T)*longest_component(G_C)>=3_implies_completion")

    bodies: list[tuple[tuple[int, ...], tuple[tuple[Q, Q], ...], Q, Q]] = []
    for body in combinations(range(1, 14), 10):
        components, measure = exact_safe_geometry(body)
        longest = max(right - left for left, right in components)
        bodies.append((body, components, measure, longest))
    require(len(bodies) == 286, "body count")
    least_longest, least_body = min((row[3], row[0]) for row in bodies)
    require(least_longest == Q(9, 1232), "global longest-component floor")
    require(Q(3, 7) / least_longest == Q(176, 3), "global tail cutoff")
    print(
        f"least_longest_component={least_longest} least_body={least_body} "
        f"global_tail_cutoff={Q(176, 3)}"
    )

    total_rows = 0
    positive_rows = 0
    isolated_only_rows = 0
    total_direct_components = 0
    total_mask_pieces = 0
    total_safe_mask_pieces = 0
    parity_counts: Counter[int] = Counter()
    circuit_counts: Counter[str] = Counter()
    per_body_counts: list[tuple[int, tuple[int, ...], Q, int]] = []
    minimum: tuple[Q, tuple[int, ...], tuple[int, int, int], int] | None = None
    first_failure: tuple[tuple[int, ...], tuple[int, int, int]] | None = None

    with Pool(processes=8) as pool:
        results = pool.map(audit_body, bodies)
    for result in results:
        (
            rows,
            positive,
            isolated,
            direct_count,
            piece_count,
            safe_piece_count,
            parity_part,
            circuit_part,
            body_minimum,
            failure,
            body_record,
        ) = result
        total_rows += rows
        positive_rows += positive
        isolated_only_rows += isolated
        total_direct_components += direct_count
        total_mask_pieces += piece_count
        total_safe_mask_pieces += safe_piece_count
        parity_counts.update(parity_part)
        circuit_counts.update(circuit_part)
        per_body_counts.append(body_record)
        if body_minimum is not None and (minimum is None or body_minimum < minimum):
            minimum = body_minimum
        if failure is not None and first_failure is None:
            first_failure = failure

    print(f"bodies={len(bodies)}")
    print(f"finite_residual_rows={total_rows}")
    print(f"positive_measure_rows={positive_rows}")
    print(f"isolated_only_rows={isolated_only_rows}")
    print(f"total_direct_positive_components={total_direct_components}")
    print(f"total_quotient_mask_pieces={total_mask_pieces}")
    print(f"total_safe_mask_pieces={total_safe_mask_pieces}")
    print(f"parity_mask_counts={sorted(parity_counts.items())}")
    print(f"low_circuit_occurrences={sorted(circuit_counts.items())}")
    print(f"maximum_body_residual={max(per_body_counts)}")
    print(f"finite_residual_minimum={minimum}")
    print(f"first_counterexample={first_failure}")
    require(total_rows == 174045, "finite residual row count")
    require(positive_rows == total_rows, "every residual row has positive measure")
    require(isolated_only_rows == 0, "no isolated-only row")
    require(first_failure is None, "counterexample")
    require(minimum is not None, "minimum exists")
    print("CONCLUSION=EVERY_C_IN_BINOM_[13]_10_ACCEPTS_EVERY_DISTINCT_POSITIVE_TERNARY_UNIT_TRIPLE_AFTER_SCALE_THREE")


if __name__ == "__main__":
    main()
