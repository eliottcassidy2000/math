#!/usr/bin/env python3
"""Clean-room independent audit of provisional THM-4112.

The implementation imports no repository module.  Its literal topology
engine builds every antipodal tooth wall as an exact ``Fraction``, classifies
the intervening open cells by direct two-phase inequalities, and reconstructs
circular components by connectivity across genuinely dangerous endpoints.
Endpoint contact is therefore kept separate from open overlap.

The audit checks the algebraic ancestry recursion and its two scale regimes,
the adaptive AP7 and parity-free families, the non-AP D0 interval in both
phases, explicit endpoint clocks, the two direct thirteen-speed suppliers,
and both parent-gap hostiles.  It does not import or inspect the primary
THM-4112 referee.
"""

from __future__ import annotations

import json
from fractions import Fraction
from hashlib import sha256
from itertools import combinations


DELTA = Fraction(1, 14)
HALF = Fraction(1, 2)
EXPECTED_SEMANTIC_SHA256 = "a716c2a1656d979a087674ebd3e92fb3a691e02a87f4bc867092d6b1dffcec6e"


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def frow(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def mod_one(value: Fraction) -> Fraction:
    return value % 1


def circle_norm(value: Fraction) -> Fraction:
    residue = mod_one(value)
    return min(residue, 1 - residue)


def phase_danger(speed: int, theta: Fraction) -> bool:
    return (
        circle_norm(speed * theta) < DELTA
        or circle_norm(speed * (theta + HALF)) < DELTA
    )


def two_phase_safe(speed: int, theta: Fraction) -> bool:
    return not phase_danger(speed, theta)


def safe_for_all(speeds: tuple[int, ...], theta: Fraction) -> bool:
    return all(two_phase_safe(speed, theta) for speed in speeds)


def parity_weight(speed: int) -> int:
    return 1 if speed % 2 == 0 else 2


def tooth_length(speed: int) -> Fraction:
    return Fraction(1, 7 * speed)


def tooth_period(speed: int) -> Fraction:
    return Fraction(1, parity_weight(speed) * speed)


def safe_gap(speed: int) -> Fraction:
    return tooth_period(speed) - tooth_length(speed)


def literal_walls(speed: int) -> set[Fraction]:
    """Solve both strict phase equalities directly, modulo one."""
    walls: set[Fraction] = set()
    for phase in (Fraction(0), HALF):
        for integer in range(speed):
            walls.add(mod_one((Fraction(integer) - DELTA) / speed - phase))
            walls.add(mod_one((Fraction(integer) + DELTA) / speed - phase))
    return walls


def danger_components(
    speeds: tuple[int, ...],
) -> list[tuple[Fraction, Fraction, Fraction, frozenset[int]]]:
    """Reconstruct literal open circular components from atomic wall cells."""
    require(bool(speeds), "nonempty component speed bank")
    points = sorted(
        {Fraction(0), Fraction(1)}
        | set().union(*(literal_walls(speed) for speed in speeds))
    )
    segment_count = len(points) - 1
    active = [
        any(
            phase_danger(speed, (points[index] + points[index + 1]) / 2)
            for speed in speeds
        )
        for index in range(segment_count)
    ]

    parent = list(range(segment_count))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def join(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    # Boundary j is points[j], with point 0 also representing point 1.
    for boundary in range(segment_count):
        left = (boundary - 1) % segment_count
        right = boundary
        boundary_danger = any(
            phase_danger(speed, points[boundary]) for speed in speeds
        )
        if boundary_danger:
            require(active[left] and active[right],
                    "strictly dangerous wall has a dangerous neighbourhood")
            join(left, right)

    groups: dict[int, set[int]] = {}
    for index, is_active in enumerate(active):
        if is_active:
            groups.setdefault(find(index), set()).add(index)

    components = []
    for members in groups.values():
        span = sum(
            (points[index + 1] - points[index] for index in members),
            Fraction(0),
        )
        if len(members) == segment_count:
            start = Fraction(0)
        else:
            start_index = next(
                index for index in members
                if (index - 1) % segment_count not in members
            )
            start = points[start_index]
        components.append(
            (start, mod_one(start + span), span, frozenset(members))
        )
    return sorted(components, key=lambda row: (row[0], row[2], row[1]))


def component_row(
    component: tuple[Fraction, Fraction, Fraction, frozenset[int]]
) -> list[list[int]]:
    return [frow(component[0]), frow(component[1]), frow(component[2])]


def component_root_teeth(
    component: tuple[Fraction, Fraction, Fraction, frozenset[int]],
    root: int,
) -> list[int]:
    start, _, span, _ = component
    end = start + span
    count = parity_weight(root) * root
    radius = Fraction(1, 14 * root)
    met: set[int] = set()
    for tooth in range(count):
        centre = Fraction(tooth, count)
        for shift in (-1, 0, 1, 2):
            left = centre - radius + shift
            right = centre + radius + shift
            if max(start, left) < min(end, right):
                met.add(tooth)
                break
    return sorted(met)


def has_consecutive_teeth(indices: list[int], count: int) -> bool:
    index_set = set(indices)
    return any((index + 1) % count in index_set for index in index_set)


def envelopes(speeds: tuple[int, ...]) -> list[Fraction]:
    require(bool(speeds), "nonempty ancestry chain")
    values = [Fraction(0)] * len(speeds)
    values[-1] = tooth_length(speeds[-1])
    for index in range(len(speeds) - 2, -1, -1):
        values[index] = tooth_length(speeds[index]) + 2 * values[index + 1]
    return values


def closed_envelope(speeds: tuple[int, ...], index: int) -> Fraction:
    return sum(
        (
            Fraction(2 ** (position - index), 7 * speeds[position])
            for position in range(index, len(speeds))
        ),
        Fraction(0),
    )


def ancestry_gates(speeds: tuple[int, ...]) -> list[bool]:
    values = envelopes(speeds)
    return [
        values[index + 1] <= safe_gap(speeds[index])
        for index in range(len(speeds) - 1)
    ]


def ceiling_ratio(value: int, ratio: Fraction) -> int:
    return (
        ratio.numerator * value + ratio.denominator - 1
    ) // ratio.denominator


def ratio_chain(first: int, ratio: Fraction, depth: int) -> tuple[int, ...]:
    values = [first]
    while len(values) < depth:
        values.append(ceiling_ratio(values[-1], ratio))
    return tuple(values)


def pair_envelope(y: int, z: int) -> Fraction:
    require(y < z, "ordered pair envelope")
    return min(Fraction(2, 7 * y), Fraction(1, 7 * y) + Fraction(2, 7 * z))


def triple_envelope(x: int, y: int, z: int) -> Fraction:
    require(x < y < z, "ordered triple envelope")
    return Fraction(1, 7 * x) + 2 * pair_envelope(y, z)


def four_root_row(
    row: tuple[int, int, int, int], root: int, core_length: Fraction
) -> dict[str, object]:
    require(root in row and len(set(row)) == 4, "four-root row")
    x, y, z = sorted(value for value in row if value != root)
    residual = triple_envelope(x, y, z)
    span = tooth_length(root) + 2 * residual
    return {
        "root": root,
        "parent_slack": frow(safe_gap(root) - residual),
        "survival_slack": frow(core_length - span),
        "parent_ok": residual <= safe_gap(root),
        "survival_ok": span <= core_length,
    }


def norm_range(
    speed: int,
    interval: tuple[Fraction, Fraction],
    phase: Fraction = Fraction(0),
) -> tuple[Fraction, Fraction]:
    left, right = interval
    breakpoints = {left, right}
    scaled_left = 2 * speed * (left + phase)
    scaled_right = 2 * speed * (right + phase)
    first = scaled_left.numerator // scaled_left.denominator - 2
    last = scaled_right.numerator // scaled_right.denominator + 3
    for integer in range(first, last + 1):
        theta = Fraction(integer, 2 * speed) - phase
        if left <= theta <= right:
            breakpoints.add(theta)
    values = [circle_norm(speed * (theta + phase)) for theta in breakpoints]
    return min(values), max(values)


def interval_is_two_phase_safe(
    speeds: tuple[int, ...], interval: tuple[Fraction, Fraction]
) -> bool:
    for speed in speeds:
        for phase in (Fraction(0), HALF):
            minimum, _ = norm_range(speed, interval, phase)
            if minimum < DELTA:
                return False
    return True


def integer_distance(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def safe_clock_labels(speeds: tuple[int, ...], modulus: int) -> list[int]:
    return [
        label
        for label in range(modulus)
        if all(
            14 * integer_distance(speed * label, modulus) >= modulus
            for speed in speeds
        )
    ]


def eligible_odd_classes(
    speeds: tuple[int, ...], modulus: int
) -> tuple[list[int], list[int]]:
    labels = safe_clock_labels(speeds, modulus)
    eligible = [
        tail
        for tail in range(1, 2 * modulus, 2)
        if all(
            7 * integer_distance(tail * label, modulus) < modulus
            for label in labels
        )
    ]
    return labels, eligible


def endpoint_clock_row(
    name: str,
    speeds: tuple[int, ...],
    theta: Fraction,
    modulus: int,
    bound: int,
) -> dict[str, object]:
    scaled = theta * modulus
    require(scaled.denominator == 1 and modulus % 2 == 0,
            "even rational clock presentation")
    label = scaled.numerator % modulus
    labels, eligible = eligible_odd_classes(speeds, modulus)
    require(label in labels and (label + modulus // 2) % modulus in labels,
            "antipodal endpoint labels are safe")
    require(not eligible, "endpoint pair empties eligible odd classes")
    require(modulus <= bound, "adaptive endpoint clock bound")
    for tail in range(1, 2 * modulus, 2):
        first = integer_distance(tail * label, modulus)
        second = integer_distance(
            tail * (label + modulus // 2), modulus
        )
        require(first + second == modulus // 2,
                "odd-tail endpoint complement identity")
    return {
        "name": name,
        "theta": frow(theta),
        "clock": modulus,
        "label": label,
        "safe_label_count": len(labels),
        "eligible_odd_classes": len(eligible),
        "bound": bound,
    }


def main() -> None:
    # Exact one-comb base: odd speeds have twice as many teeth, but every
    # literal open component has the same exact span a_v.
    single_comb_rows = []
    for speed in range(1, 25):
        components = danger_components((speed,))
        require(len(components) == parity_weight(speed) * speed,
                "single-comb component count")
        require(all(component[2] == tooth_length(speed) for component in components),
                "single-comb exact-span base")
        if speed in {1, 2, 3, 4, 7, 14, 23, 24}:
            single_comb_rows.append(
                [speed, parity_weight(speed), len(components), frow(components[0][2])]
            )

    # Reverse recursion equals the closed 1,2,4,... formula through depth 12.
    recursion_rows = []
    for ratio in (Fraction(2), Fraction(12, 5), Fraction(3)):
        speeds = ratio_chain(5, ratio, 12)
        values = envelopes(speeds)
        require(
            all(values[index] == closed_envelope(speeds, index)
                for index in range(len(speeds))),
            "closed ancestry formula through depth twelve",
        )
        recursion_rows.append(
            [frow(ratio), speeds[0], speeds[-1], frow(values[0]), sum(ancestry_gates(speeds))]
        )

    # Literal suffix controls at the equality-sensitive finite-depth boundary
    # and in an arbitrary-depth admissible ratio-three chain.
    literal_chain_rows = []
    for speeds in (
        (1, 2, 4, 8, 16, 32),
        (1, 3, 9, 27, 81, 243, 729),
    ):
        values = envelopes(speeds)
        require(all(ancestry_gates(speeds)), "literal chain parent gaps")
        suffix_rows = []
        for index in range(len(speeds)):
            components = danger_components(speeds[index:])
            maximum = max(component[2] for component in components)
            if index == len(speeds) - 1:
                require(maximum == values[index], "literal exact-span recursion base")
            else:
                require(maximum < values[index], "literal strict ancestry envelope")
            suffix_rows.append([index, len(components), frow(maximum), frow(values[index])])
        literal_chain_rows.append([list(speeds), suffix_rows])
    require(
        envelopes((1, 2, 4, 8, 16, 32))[1] == safe_gap(1),
        "ratio-two depth-six odd-parent equality gate",
    )

    # Finite doubling and arbitrary-depth geometric majorants.
    finite_thresholds = [
        ["AP7", frow(Fraction(9, 490)), 4, 32],
        ["AP8", frow(Fraction(3, 392)), 5, 94],
        ["D0", frow(Fraction(3, 56)), 6, 16],
    ]
    for _, length_row, depth, threshold in finite_thresholds:
        length = Fraction(*length_row)
        require(Fraction(depth, 7 * threshold) <= length,
                "finite doubling threshold succeeds")
        require(Fraction(depth, 7 * (threshold - 1)) > length,
                "finite doubling threshold is least integer")
        chain = tuple(threshold * (2**index) for index in range(depth))
        require(all(ancestry_gates(chain)), "finite doubling parent gates")
        require(envelopes(chain)[0] <= length, "finite doubling top envelope")

    geometric_thresholds = [
        ["AP7", frow(Fraction(9, 490)), frow(Fraction(12, 5)), 47],
        ["AP8", frow(Fraction(3, 392)), frow(Fraction(12, 5)), 112],
        ["D0", frow(Fraction(3, 56)), frow(Fraction(12, 5)), 16],
        ["AP8", frow(Fraction(3, 392)), frow(Fraction(3)), 56],
    ]
    for _, length_row, ratio_row, threshold in geometric_thresholds:
        length = Fraction(*length_row)
        ratio = Fraction(*ratio_row)
        required = ratio / (7 * length * (ratio - 2))
        require(threshold - 1 < required <= threshold,
                "arbitrary-depth first-speed threshold")
        chain = ratio_chain(threshold, ratio, 12)
        values = envelopes(chain)
        require(all(ancestry_gates(chain)), "arbitrary-depth parent gaps")
        require(values[0] <= length, "arbitrary-depth finite top envelope")
        for index in range(len(chain) - 1):
            require(
                values[index + 1]
                <= Fraction(1, 7 * chain[index] * (ratio - 2)),
                "geometric suffix majorant",
            )
            require(
                Fraction(1, 7 * chain[index] * (ratio - 2))
                <= safe_gap(chain[index]),
                "lambda at least twelve-fifths gap gate",
            )

    ap7 = tuple(range(1, 8))
    ap7_interval = (Fraction(4, 35), Fraction(13, 98))
    ap7_length = Fraction(9, 490)
    require(interval_is_two_phase_safe(ap7, ap7_interval),
            "inherited AP7 interval independently replayed")
    require(
        7 * ap7_length
        - (Fraction(1, 84) + Fraction(2, 85) + Fraction(8, 86))
        == Fraction(1, 8772),
        "AP7 adaptive case-one boundary",
    )
    require(
        7 * ap7_length
        - (Fraction(1, 86) + Fraction(2, 85) + Fraction(8, 87))
        == Fraction(650, 445179),
        "AP7 adaptive case-two boundary",
    )

    adaptive_counts = {"a_even": 0, "b_even_close": 0, "b_even_far": 0}
    adaptive_min: dict[str, tuple[Fraction, tuple[int, int, int, int]]] = {}
    for row in combinations(range(84, 106), 4):
        a, b, _, _ = row
        if a % 2 == 0:
            case = "a_even"
            root = a
        elif b % 2 == 0 and b < 2 * a:
            case = "b_even_close"
            root = b
        else:
            continue
        certificate = four_root_row(row, root, ap7_length)
        require(certificate["parent_ok"] and certificate["survival_ok"],
                "adaptive-root AP7 finite hostile probe")
        adaptive_counts[case] += 1
        slack = Fraction(*certificate["survival_slack"])
        if case not in adaptive_min or slack < adaptive_min[case][0]:
            adaptive_min[case] = (slack, row)

    for a in range(85, 101, 2):
        for offset in range(4):
            b = 2 * a + 2 * offset
            row = (a, b, b + 1, b + 2)
            certificate = four_root_row(row, a, ap7_length)
            require(certificate["parent_ok"] and certificate["survival_ok"],
                    "adaptive-root AP7 far-even case")
            adaptive_counts["b_even_far"] += 1

    require(adaptive_min["a_even"] == (Fraction(1, 61404), (84, 85, 86, 87)),
            "AP7 case-one exact worst row")
    require(
        adaptive_min["b_even_close"]
        == (Fraction(650, 3116253), (85, 86, 87, 88)),
        "AP7 case-two exact worst row",
    )

    parity_free_ap7_checks = 0
    for a in range(47, 61):
        for b in range(2 * a, 2 * a + 5):
            for c in range(b + 1, b + 4):
                for d in range(c + 1, c + 4):
                    certificate = four_root_row((a, b, c, d), a, ap7_length)
                    require(certificate["parent_ok"] and certificate["survival_ok"],
                            "parity-free AP7 family")
                    parity_free_ap7_checks += 1
    require(Fraction(6, 47) < 7 * ap7_length,
            "parity-free AP7 analytic threshold")

    ap7_control_rows = [
        ((85, 86, 91, 101), Fraction(81, 707), 1_414),
        ((47, 95, 97, 99), Fraction(4, 35), 70),
    ]
    for outliers, theta, _ in ap7_control_rows:
        require(safe_for_all(ap7 + outliers, theta), "AP7 exact control survivor")
    require([value % 2 for value in ap7_control_rows[0][0]].count(1) == 3,
            "AP7 parity-weight-seven control")
    require(
        [
            ap7_control_rows[0][0][index + 1] - ap7_control_rows[0][0][index]
            for index in range(3)
        ]
        == [1, 5, 10],
        "AP7 control avoids selected gaps four eight twelve",
    )

    d0 = (3, 4, 5, 6, 8, 10, 12)
    d0_interval = (Fraction(1, 42), Fraction(13, 168))
    d0_length = Fraction(3, 56)
    expected_d0_zero = {
        3: (Fraction(1, 14), Fraction(13, 56)),
        4: (Fraction(2, 21), Fraction(13, 42)),
        5: (Fraction(5, 42), Fraction(65, 168)),
        6: (Fraction(1, 7), Fraction(13, 28)),
        8: (Fraction(4, 21), Fraction(1, 2)),
        10: (Fraction(19, 84), Fraction(1, 2)),
        12: (Fraction(1, 14), Fraction(1, 2)),
    }
    d0_range_rows = []
    for speed in d0:
        zero_range = norm_range(speed, d0_interval)
        half_range = norm_range(speed, d0_interval, HALF)
        require(zero_range == expected_d0_zero[speed], "D0 zero-phase affine range")
        require(zero_range[0] >= DELTA and half_range[0] >= DELTA,
                "D0 both phases weak-safe")
        if speed % 2:
            require(
                half_range == (HALF - zero_range[1], HALF - zero_range[0]),
                "D0 odd half-phase complement range",
            )
        else:
            require(half_range == zero_range, "D0 even repeated phase")
        d0_range_rows.append(
            [speed, frow(zero_range[0]), frow(zero_range[1]),
             frow(half_range[0]), frow(half_range[1])]
        )
    require(interval_is_two_phase_safe(d0, d0_interval), "D0 interval inclusion")

    parity_free_d0_checks = 0
    for a in range(16, 31):
        for b in range(2 * a, 2 * a + 4):
            for c in range(b + 1, b + 3):
                for d in range(c + 1, c + 3):
                    certificate = four_root_row((a, b, c, d), a, d0_length)
                    require(certificate["parent_ok"] and certificate["survival_ok"],
                            "parity-free D0 family")
                    parity_free_d0_checks += 1
    require(Fraction(6, 16) == 7 * d0_length,
            "D0 equality threshold")

    d0_crossing_rows = []
    expected_crossings = {
        (15, 31, 33, 35): (15, Fraction(19, 95480), True),
        (28, 29, 30, 31): (28, Fraction(89, 170520), True),
        (26, 27, 28, 29): (26, Fraction(-457, 137592), False),
    }
    for row, (expected_root, expected_slack, expected_success) in expected_crossings.items():
        candidates = [
            four_root_row(row, root, d0_length)
            for root in row
        ]
        parent_valid = [candidate for candidate in candidates if candidate["parent_ok"]]
        require(bool(parent_valid), "D0 crossing has a parent-valid root")
        best = max(
            parent_valid,
            key=lambda candidate: Fraction(*candidate["survival_slack"]),
        )
        require(best["root"] == expected_root, "D0 adaptive crossing root")
        require(Fraction(*best["survival_slack"]) == expected_slack,
                "D0 exact crossing slack")
        require(best["survival_ok"] == expected_success,
                "D0 crossing sufficient-criterion status")
        d0_crossing_rows.append([list(row), best])
    require(safe_for_all(d0 + (15, 31, 33, 35), Fraction(1, 42)),
            "D0 all-odd crossing survivor")

    # Direct thirteen-speed suppliers: component counts come from the literal
    # endpoint-cell engine, independently of the ancestry bound.
    direct_rows = []
    direct_data = [
        (
            "D0+six",
            d0,
            d0_interval,
            (16, 32, 64, 128, 256, 512),
            [Fraction(3, 56), Fraction(5, 224), Fraction(1, 112),
             Fraction(3, 896), Fraction(1, 896), Fraction(1, 3584)],
            416,
            Fraction(1, 112),
            Fraction(43, 1792),
        ),
        (
            "AP8+five",
            tuple(range(1, 9)),
            (Fraction(11, 49), Fraction(13, 56)),
            (94, 188, 376, 752, 1504),
            [Fraction(5, 658), Fraction(1, 329), Fraction(3, 2632),
             Fraction(1, 2632), Fraction(1, 10528)],
            1_316,
            Fraction(1, 658),
            Fraction(11, 49),
        ),
    ]
    for name, core, interval, outliers, expected_levels, expected_count, expected_max, theta in direct_data:
        require(interval_is_two_phase_safe(core, interval),
                "direct supplier inherited core interval")
        levels = envelopes(outliers)
        require(levels == expected_levels, "direct supplier ancestry levels")
        require(all(ancestry_gates(outliers)), "direct supplier parent gaps")
        require(levels[0] <= interval[1] - interval[0],
                "direct supplier top envelope")
        components = danger_components(outliers)
        maximum = max(component[2] for component in components)
        require(len(components) == expected_count, "direct supplier component count")
        require(maximum == expected_max, "direct supplier actual maximum span")
        require(safe_for_all(core + outliers, theta), "direct supplier exact survivor")
        direct_rows.append(
            {
                "name": name,
                "outliers": list(outliers),
                "levels": [frow(value) for value in levels],
                "components": len(components),
                "maximum_span": frow(maximum),
                "survivor": frow(theta),
            }
        )

    # Endpoint-clock scope: only eleven-speed cores compile the dyadic seam.
    clock_rows = [
        endpoint_clock_row(
            "AP7-adaptive",
            ap7 + ap7_control_rows[0][0],
            ap7_control_rows[0][1],
            ap7_control_rows[0][2],
            14 * max(ap7_control_rows[0][0]),
        ),
        endpoint_clock_row(
            "AP7-parity-free",
            ap7 + ap7_control_rows[1][0],
            ap7_control_rows[1][1],
            ap7_control_rows[1][2],
            14 * max(ap7_control_rows[1][0]),
        ),
        endpoint_clock_row(
            "D0-all-odd",
            d0 + (15, 31, 33, 35),
            Fraction(1, 42),
            42,
            14 * 35,
        ),
    ]
    require(all(len(row["name"]) > 0 for row in clock_rows), "clock row labels")

    # First parent-gap hostile: exact component, three root teeth, and excess
    # over the false no-gap envelope.
    hostile_one = (85, 87, 89, 91)
    residual_q = triple_envelope(87, 89, 91)
    require(residual_q == Fraction(437, 54201), "first hostile residual envelope")
    require(residual_q > safe_gap(85), "first hostile parent-gap failure")
    hostile_components = danger_components(hostile_one)
    target = [
        component for component in hostile_components
        if component[0] == Fraction(69, 1274)
        and component[1] == Fraction(46, 637)
    ]
    require(len(target) == 1 and target[0][2] == Fraction(23, 1274),
            "first hostile exact component")
    met_teeth = component_root_teeth(target[0], 85)
    require(met_teeth == [10, 11, 12], "first hostile meets three root teeth")
    false_envelope = tooth_length(85) + 2 * residual_q
    require(false_envelope == Fraction(11719, 658155),
            "first hostile false envelope")
    require(
        target[0][2] - false_envelope == Fraction(207559, 838489470),
        "first hostile exact excess",
    )
    require(not any(phase_danger(speed, target[0][0]) for speed in hostile_one),
            "first hostile left endpoint remains a safe separator")
    require(not any(phase_danger(speed, target[0][1]) for speed in hostile_one),
            "first hostile right endpoint remains a safe separator")

    # All-root hostile: every possible root loses the one-parent gate.
    hostile_two = (43, 45, 47, 49)
    all_root_rows = []
    for root in hostile_two:
        residual = tuple(value for value in hostile_two if value != root)
        components = danger_components(residual)
        witnesses = []
        for component in components:
            teeth = component_root_teeth(component, root)
            if has_consecutive_teeth(teeth, parity_weight(root) * root):
                witnesses.append((component, teeth))
        require(bool(witnesses), "all-root hostile has a two-tooth residual component")
        component, teeth = witnesses[0]
        all_root_rows.append(
            [root, len(components), len(witnesses), component_row(component), teeth]
        )
    require(safe_for_all(ap7 + hostile_two, Fraction(4, 35)),
            "all-root hostile remains antipodal-safe")

    ledger = {
        "single_comb_rows": single_comb_rows,
        "recursion_rows": recursion_rows,
        "literal_chain_rows": literal_chain_rows,
        "finite_thresholds": finite_thresholds,
        "geometric_thresholds": geometric_thresholds,
        "adaptive_counts": adaptive_counts,
        "adaptive_min": {
            key: [frow(value[0]), list(value[1])]
            for key, value in adaptive_min.items()
        },
        "parity_free_checks": [parity_free_ap7_checks, parity_free_d0_checks],
        "ap7_controls": [
            [list(outliers), frow(theta), clock]
            for outliers, theta, clock in ap7_control_rows
        ],
        "d0_ranges": d0_range_rows,
        "d0_crossings": d0_crossing_rows,
        "direct_rows": direct_rows,
        "clock_rows": clock_rows,
        "parent_gap_hostile": {
            "row": list(hostile_one),
            "residual_Q": frow(residual_q),
            "gap_85": frow(safe_gap(85)),
            "component": component_row(target[0]),
            "root_teeth": met_teeth,
            "false_envelope": frow(false_envelope),
            "excess": frow(target[0][2] - false_envelope),
        },
        "all_root_hostile": {
            "row": list(hostile_two),
            "witnesses": all_root_rows,
            "safe_theta": frow(Fraction(4, 35)),
        },
        "scope": (
            "component ancestry and explicit scale-separated suppliers; "
            "dyadic compilation only for eleven-speed cores; not arbitrary LRC14"
        ),
    }
    semantic = canonical_digest(ledger)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                "independent semantic digest")

    print("THM-4112 ANTIPODAL COMPONENT ANCESTRY INDEPENDENT AUDIT")
    print("verdict=ACCEPT")
    print("method=Fraction-only literal two-phase wall cells and strict endpoint connectivity")
    print("single_comb_exact_span_rows=", single_comb_rows)
    print("recursion_depth12_rows=", recursion_rows)
    print("literal_chain_rows=", literal_chain_rows)
    print("finite_ratio_two_thresholds=", finite_thresholds)
    print("arbitrary_depth_geometric_thresholds=", geometric_thresholds)
    print("AP7_adaptive_case_counts=", adaptive_counts)
    print("AP7_adaptive_exact_minima=", ledger["adaptive_min"])
    print("parity_free_AP7_D0_checks=", ledger["parity_free_checks"])
    print("AP7_control_rows=", ledger["ap7_controls"])
    print("D0_both_phase_ranges=", d0_range_rows)
    print("D0_crossing_rows=", d0_crossing_rows)
    print("direct_thirteen_speed_rows=", direct_rows)
    print("endpoint_clock_rows=", clock_rows)
    print("parent_gap_hostile=", ledger["parent_gap_hostile"])
    print("all_root_hostile=", ledger["all_root_hostile"])
    print("semantic_sha256=", semantic)
    print("smallest_failure=none")
    print("scope=conditional component theorem and explicit suppliers; no arbitrary-core LRC14")
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
