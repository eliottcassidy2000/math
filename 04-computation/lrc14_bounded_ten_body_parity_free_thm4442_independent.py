#!/usr/bin/env python3
"""Clean-room exact audit of the parity-free bounded ten-body theorem.

No candidate module is imported.  Every literal thirteen-speed row is rebuilt
on the integer grid 14*lcm(speeds); open cells are tested using integer modular
arithmetic.  Workers shard only by the 286 ten-bodies.
"""

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import lcm
from multiprocessing import Pool


EXPECTED_MINIMUM = (
    Q(10517879, 643242600),
    (1, 2, 3, 5, 6, 7, 8, 9, 11, 13),
    (8, 34, 50),
)


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def midpoint_is_safe(speeds: tuple[int, ...], numerator: int, denominator: int) -> bool:
    for speed in speeds:
        residue = (speed * numerator) % denominator
        if 14 * min(residue, denominator - residue) < denominator:
            return False
    return True


def geometry(speeds: tuple[int, ...]) -> tuple[tuple[tuple[Q, Q], ...], Q]:
    common = lcm(*speeds)
    denominator = 14 * common
    cuts = {0, denominator}
    for speed in speeds:
        unit = common // speed
        for k in range(speed + 1):
            base = 14 * k * unit
            for shift in (-unit, unit):
                point = base + shift
                if 0 <= point <= denominator:
                    cuts.add(point)
    ordered = sorted(cuts)
    raw: list[list[int]] = []
    for left, right in zip(ordered, ordered[1:]):
        if midpoint_is_safe(speeds, left + right, 2 * denominator):
            if raw and raw[-1][1] == left:
                raw[-1][1] = right
            else:
                raw.append([left, right])
    for left, right in raw:
        for endpoint in (left, right):
            require(
                midpoint_is_safe(speeds, endpoint, denominator),
                f"non-safe weak endpoint {speeds}",
            )
    components = tuple((Q(left, denominator), Q(right, denominator)) for left, right in raw)
    measure = Q(sum(right - left for left, right in raw), denominator)
    return components, measure


def fraction_is_safe(speeds: tuple[int, ...], x: Q) -> bool:
    for speed in speeds:
        residue = (speed * x) % 1
        if min(residue, 1 - residue) < Q(1, 14):
            return False
    return True


def residual_values(longest: Q) -> tuple[int, ...]:
    numerator = 3 * longest.denominator
    denominator = 7 * longest.numerator
    ceiling = (numerator + denominator - 1) // denominator
    return tuple(
        value for value in range(1, ceiling)
        if value % 3 and 7 * value * longest < 3
    )


def audit_body(payload: tuple[tuple[int, ...], Q]) -> tuple:
    body, longest = payload
    values = residual_values(longest)
    rows = component_count = 0
    parity = Counter()
    minimum = None
    for tails in combinations(values, 3):
        full = tuple(sorted(tuple(3 * speed for speed in body) + tails))
        require(len(full) == len(set(full)) == 13, "physical labels")
        components, measure = geometry(full)
        require(measure > 0 and components, f"counterexample {body} {tails}")
        left, right = components[0]
        x = (left + right) / 2
        three_x = 3 * x
        sheet = three_x.numerator // three_x.denominator
        y = three_x % 1
        require(sheet in (0, 1, 2) and x == (y + sheet) / 3, "reverse lift identity")
        require(fraction_is_safe(body, y), "reverse-lift body failure")
        require(fraction_is_safe(tails, x), "reverse-lift tail failure")
        require(fraction_is_safe(full, x), "literal witness failure")
        rows += 1
        component_count += len(components)
        parity[sum((tail & 1) << i for i, tail in enumerate(tails))] += 1
        record = (measure, body, tails, len(components))
        if minimum is None or record < minimum:
            minimum = record
    return body, longest, values[-1], rows, component_count, dict(parity), minimum


def main() -> None:
    body_payloads = []
    for body in combinations(range(1, 14), 10):
        components, _ = geometry(body)
        require(components, "body has no positive component")
        body_payloads.append((body, max(right - left for left, right in components)))
    require(len(body_payloads) == 286, "body count")
    least_longest, least_body = min((longest, body) for body, longest in body_payloads)
    require(least_longest == Q(9, 1232), "least longest component")
    require(Q(3, 7) / least_longest == Q(176, 3), "global cutoff")

    with Pool(processes=8) as pool:
        results = pool.map(audit_body, body_payloads)

    total_rows = 0
    total_components = 0
    parity = Counter()
    minimum = None
    maximum_tail = 0
    maximum_body_residual = None
    for body, longest, residual_max, rows, component_count, parity_part, local_minimum in results:
        total_rows += rows
        total_components += component_count
        parity.update(parity_part)
        maximum_tail = max(maximum_tail, residual_max)
        body_record = (rows, body, longest, residual_max)
        if maximum_body_residual is None or body_record > maximum_body_residual:
            maximum_body_residual = body_record
        if minimum is None or local_minimum < minimum:
            minimum = local_minimum

    require(total_rows == 174045, "residual row count")
    require(maximum_tail == 58, "residual height")
    require(minimum is not None and minimum[:3] == EXPECTED_MINIMUM, "residual minimum")
    print("INDEPENDENT_TEN_BODY_ALL_PARITY_AUDIT=PASS")
    print("engine=integer_grid_14_lcm_plus_modular_midpoint_classification")
    print(f"bodies={len(body_payloads)} least_longest={least_longest} least_body={least_body} global_cutoff={Q(176,3)}")
    print(f"residual_rows={total_rows} positive_rows={total_rows} residual_components={total_components} max_tail={maximum_tail}")
    print(f"parity_mask_counts={sorted(parity.items())}")
    print(f"maximum_body_residual={maximum_body_residual}")
    print(f"residual_minimum={minimum}")
    print("endpoint_policy=strict_open_danger;all_reported_components_have_exact_weak_safe_endpoints")
    print("lift_policy=literal_physical_midpoint_then_exact_reverse_lift_y=3x_mod_1")


if __name__ == "__main__":
    main()
