#!/usr/bin/env python3
"""Fraction-exact referee for THM-4100.

The analytic theorem is independent of this finite battery.  The battery
checks literal open antipodal components, every AP8 threshold row in the
displayed finite universe, adaptive even clocks and odd owner classes, and
the two scalar-blindness hostile pairs.  All gates raise through ``require``
and therefore remain active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction as F
from heapq import merge
from itertools import combinations
from math import gcd


DELTA = F(1, 14)
AP8_LEFT = F(11, 49)
AP8_RIGHT = F(13, 56)
AP8_LENGTH = F(3, 392)
AP8_Q = 98

D_STAR = (1, 3, 4, 5, 7, 8, 9, 11, 12)
OBSTRUCTED_ELEVEN = tuple(v for v in range(1, 13) if v != 10)
SAFE_ELEVEN = tuple(v for v in range(1, 13) if v != 4)

AP_COVER14 = tuple(range(1, 13)) + (182,)
NONAP_COVER14 = tuple(v for v in range(1, 14) if v != 3) + (182,)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parity_weight(speed: int) -> int:
    return 1 if speed % 2 == 0 else 2


def circle_norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def literal_safe(speed: int, theta: F) -> bool:
    return all(
        circle_norm(speed * (theta + F(half, 2))) >= DELTA
        for half in (0, 1)
    )


def mod_norm(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def even_presentation_denominator(value: F) -> int:
    return value.denominator if value.denominator % 2 == 0 else 2 * value.denominator


def lifted_teeth(speed: int, left: F, right: F) -> tuple[tuple[F, F], ...]:
    """Literal open teeth meeting a lifted interval, in endpoint order."""
    weight = parity_weight(speed)
    count = weight * speed
    radius = F(1, 14 * speed)
    lower = count * (left - radius)
    upper = count * (right + radius)
    first = lower.numerator // lower.denominator - 2
    last = upper.numerator // upper.denominator + 3
    return tuple(
        (F(index, count) - radius, F(index, count) + radius)
        for index in range(first, last + 1)
        if F(index, count) + radius > left
        and F(index, count) - radius < right
    )


def periodic_teeth(speed: int) -> tuple[tuple[F, F], ...]:
    """Enough ordered teeth to recover every circular component once."""
    count = parity_weight(speed) * speed
    radius = F(1, 14 * speed)
    return tuple(
        (F(index, count) - radius, F(index, count) + radius)
        for index in range(-count - 2, 2 * count + 3)
    )


def circular_component_lengths(
    speeds: tuple[int, ...], tooth_bank: dict[int, tuple[tuple[F, F], ...]]
) -> tuple[tuple[F, F, F], ...]:
    """Return one full lifted representative of each literal circular component."""
    merged: list[tuple[F, F]] = []
    active_left: F | None = None
    active_right: F | None = None
    for left, right in merge(*(tooth_bank[speed] for speed in speeds)):
        if active_left is not None and left < active_right:
            active_right = max(active_right, right)
        else:
            if active_left is not None:
                merged.append((active_left, active_right))
            active_left, active_right = left, right
    if active_left is not None:
        merged.append((active_left, active_right))

    # Component left endpoints repeat with period one.  Selecting starts in
    # [0,1) counts a seam-crossing circular component exactly once and keeps
    # its untruncated lifted length.
    return tuple(
        (left, right, right - left)
        for left, right in merged
        if 0 <= left < 1
    )


def arrangement_witness(outliers: tuple[int, int, int]) -> tuple[F, int, str]:
    candidates: list[tuple[F, int, str]] = [
        (AP8_LEFT, AP8_Q, "core-left"),
        (AP8_RIGHT, 56, "core-right"),
    ]
    for speed in outliers:
        for left, right in lifted_teeth(speed, AP8_LEFT, AP8_RIGHT):
            if AP8_LEFT <= left <= AP8_RIGHT:
                candidates.append((left, 14 * speed, f"{speed}-left"))
            if AP8_LEFT <= right <= AP8_RIGHT:
                candidates.append((right, 14 * speed, f"{speed}-right"))

    body = tuple(range(1, 9)) + outliers
    for theta, natural_clock, source in sorted(set(candidates)):
        if all(literal_safe(speed, theta) for speed in body):
            clock = even_presentation_denominator(theta)
            require(clock <= natural_clock, "reduced even presentation exceeded natural clock")
            return theta, clock, source
    raise RuntimeError(f"no AP8 arrangement witness for {outliers}")


def safe_packet(speeds: tuple[int, ...], clock: int) -> tuple[int, ...]:
    return tuple(
        label
        for label in range(clock)
        if all(14 * mod_norm(speed * label, clock) >= clock for speed in speeds)
    )


def eligible_odd_classes(safe_labels: tuple[int, ...], clock: int) -> tuple[int, ...]:
    return tuple(
        odd_class
        for odd_class in range(1, 2 * clock, 2)
        if all(7 * mod_norm(odd_class * label, clock) < clock for label in safe_labels)
    )


def equality_walls(speeds: tuple[int, ...]) -> tuple[F, ...]:
    walls = {F(0)}
    for speed in speeds:
        weight = parity_weight(speed)
        count = weight * speed
        radius = F(1, 14 * speed)
        for index in range(count):
            walls.add((F(index, count) - radius) % 1)
            walls.add((F(index, count) + radius) % 1)
    return tuple(sorted(walls))


def exhaustive_antipodal_safe_points(speeds: tuple[int, ...]) -> tuple[F, ...]:
    """Test every equality wall and one point in every complementary open cell."""
    walls = equality_walls(speeds)
    probes = list(walls)
    for index, left in enumerate(walls):
        right = walls[(index + 1) % len(walls)]
        if index + 1 == len(walls):
            right += 1
        probes.append(((left + right) / 2) % 1)
    return tuple(theta for theta in probes if all(literal_safe(v, theta) for v in speeds))


def weighted_total(speeds: tuple[int, ...]) -> int:
    return sum(parity_weight(speed) for speed in speeds)


def is_cover14(speeds: tuple[int, ...]) -> bool:
    return all(any(speed % divisor == 0 for speed in speeds) for divisor in range(2, 15))


def is_primitive(speeds: tuple[int, ...]) -> bool:
    common = 0
    for speed in speeds:
        common = gcd(common, speed)
    return common == 1


def is_arithmetic_progression(speeds: tuple[int, ...]) -> bool:
    ordered = tuple(sorted(speeds))
    if len(ordered) <= 2:
        return True
    step = ordered[1] - ordered[0]
    return step > 0 and all(
        ordered[index + 1] - ordered[index] == step
        for index in range(1, len(ordered) - 1)
    )


def lonely_value(speeds: tuple[int, ...], theta: F) -> F:
    return min(circle_norm(speed * theta) for speed in speeds)


def exact_lonely_maximum(speeds: tuple[int, ...]) -> tuple[F, F, tuple[int, ...], int]:
    """Exact maximum from every tent vertex and pairwise affine intersection."""
    candidates = {F(1, 2)}
    for speed in speeds:
        for index in range(speed):
            theta = F(2 * index + 1, 2 * speed)
            if theta <= F(1, 2):
                candidates.add(theta)
    for left, right in combinations(speeds, 2):
        for denominator in (left + right, abs(right - left)):
            if denominator == 0:
                continue
            for numerator in range(1, denominator // 2 + 1):
                candidates.add(F(numerator, denominator))

    best_value = F(0)
    best_theta = F(0)
    for theta in candidates:
        value = lonely_value(speeds, theta)
        if value > best_value or (value == best_value and theta < best_theta):
            best_value = value
            best_theta = theta
    owners = tuple(speed for speed in speeds if circle_norm(speed * best_theta) == best_value)
    return best_value, best_theta, owners, len(candidates)


def transitive_order_fingerprint(speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    ordered = tuple(sorted(speeds))
    return tuple((left, right) for left in range(len(ordered)) for right in range(left + 1, len(ordered)))


def main() -> None:
    require(AP8_RIGHT - AP8_LEFT == AP8_LENGTH, "AP8 interval length changed")
    for speed in range(1, 9):
        require(literal_safe(speed, AP8_LEFT), "AP8 left endpoint lost safety")
        require(literal_safe(speed, AP8_RIGHT), "AP8 right endpoint lost safety")

    # Literal pair and three-comb component universes.  Open intervals merge
    # only under strict overlap; a touching equality wall remains a safe
    # separator.  The pair sweep separately checks both inputs to the minimum
    # used by the strengthened three-comb envelope.
    tooth_bank = {speed: periodic_teeth(speed) for speed in range(1, 51)}
    pair_rows = 0
    pair_component_count = 0
    pair_maximum_ratio = F(0)
    pair_maximum_row: tuple[object, ...] | None = None
    for pair in combinations(range(1, 51), 2):
        lower, upper = pair
        density_bound = F(2, 7 * lower)
        span_bound = F(1, 7 * lower) + F(2, 7 * upper)
        pair_bound = min(density_bound, span_bound)
        for left, right, length in circular_component_lengths(pair, tooth_bank):
            pair_component_count += 1
            require(length < density_bound, f"two-comb density bound failed for {pair}")
            require(length < span_bound, f"two-comb span bound failed for {pair}")
            ratio = length / pair_bound
            if ratio > pair_maximum_ratio:
                pair_maximum_ratio = ratio
                pair_maximum_row = (pair, left, right, length, pair_bound)
        pair_rows += 1

    expected_pair_maximum = (
        (16, 41),
        F(139, 574),
        F(74, 287),
        F(9, 574),
        F(73, 4592),
    )
    require(pair_rows == 1225, "pair census changed")
    require(pair_component_count == 73066, "pair component census changed")
    require(pair_maximum_ratio == F(72, 73), "pair maximum ratio changed")
    require(pair_maximum_row == expected_pair_maximum, "pair maximum row changed")

    component_rows = 0
    component_count = 0
    span_branch_rows = 0
    maximum_ratio = F(0)
    maximum_row: tuple[object, ...] | None = None
    for outliers in combinations(range(1, 51), 3):
        lower, middle, upper = outliers
        density_pair_bound = F(2, 7 * middle)
        span_pair_bound = F(1, 7 * middle) + F(2, 7 * upper)
        pair_bound = min(density_pair_bound, span_pair_bound)
        span_branch_rows += span_pair_bound < density_pair_bound
        bound = F(1, 7 * lower) + 2 * pair_bound
        for left, right, length in circular_component_lengths(outliers, tooth_bank):
            component_count += 1
            require(length < bound, f"component bound failed for {outliers}")
            ratio = length / bound
            if ratio > maximum_ratio:
                maximum_ratio = ratio
                maximum_row = (outliers, left, right, length, bound)
        component_rows += 1

    expected_maximum = (
        (1, 42, 47),
        F(279, 658),
        F(379, 658),
        F(50, 329),
        F(23, 147),
    )
    require(component_rows == 19600, "component triple census changed")
    require(component_count == 1390420, "circular component census changed")
    require(span_branch_rows == 4600, "sharpened span-branch census changed")
    require(maximum_ratio == F(1050, 1081), "maximum component ratio changed")
    require(maximum_row == expected_maximum, "maximum component row changed")

    # Every AP8 row in the threshold bank, with adaptive clock and the two
    # complementary labels replayed against every odd class modulo 2N.
    ap8_rows = 0
    owner_classes = 0
    maximum_clock = 0
    source_left = 0
    source_right = 0
    for outliers in combinations(range(93, 121), 3):
        lower, middle, upper = outliers
        sharpened_reciprocal = F(1, lower) + min(
            F(4, middle), F(2, middle) + F(4, upper)
        )
        require(sharpened_reciprocal <= F(3, 56), "AP8 reciprocal component gate failed")
        theta, clock, source = arrangement_witness(outliers)
        require(clock % 2 == 0 and clock <= 14 * upper, "adaptive clock bound failed")
        label = theta * clock
        require(label.denominator == 1, "clock does not present witness")
        integer_label = int(label)
        half_label = (integer_label + clock // 2) % clock
        body = tuple(range(1, 9)) + outliers
        require(
            all(14 * mod_norm(speed * integer_label, clock) >= clock for speed in body),
            "first manufactured label is unsafe",
        )
        require(
            all(14 * mod_norm(speed * half_label, clock) >= clock for speed in body),
            "half-shifted manufactured label is unsafe",
        )
        for odd_class in range(1, 2 * clock, 2):
            owner_classes += 1
            first_bad = 7 * mod_norm(odd_class * integer_label, clock) < clock
            second_bad = 7 * mod_norm(odd_class * half_label, clock) < clock
            require(not (first_bad and second_bad), "odd class survived complementary labels")
        source_left += source.endswith("left")
        source_right += source.endswith("right")
        maximum_clock = max(maximum_clock, clock)
        ap8_rows += 1

    require(ap8_rows == 3276, "AP8 threshold census changed")
    require(owner_classes == 2227526, "odd owner-class census changed")
    require(maximum_clock == 1680, "maximum adaptive clock changed")
    require((source_left, source_right) == (1771, 1505), "source-side census changed")

    boundary_outliers = (93, 94, 95)
    boundary_body = tuple(range(1, 9)) + boundary_outliers
    boundary_theta, boundary_clock, boundary_source = arrangement_witness(boundary_outliers)
    boundary_label = int(boundary_theta * boundary_clock)
    boundary_packet = safe_packet(boundary_body, boundary_clock)
    boundary_eligible = eligible_odd_classes(boundary_packet, boundary_clock)
    reciprocal_slack = F(3, 392) - F(1, 7 * 93) - F(4, 7 * 94)
    require(boundary_theta == F(11, 49), "boundary witness changed")
    require((boundary_clock, boundary_label) == (98, 22), "boundary clock changed")
    require(boundary_source == "core-left", "boundary source changed")
    require(len(boundary_packet) == 34, "boundary safe-packet size changed")
    require(boundary_eligible == (), "boundary eligible bank is nonempty")
    require(reciprocal_slack == F(65, 1713432), "boundary reciprocal slack changed")

    require(
        F(1, 92) + min(F(4, 93), F(2, 93) + F(4, 94)) - F(3, 56)
        == F(37, 119784),
        "uniform threshold-92 failure boundary changed",
    )
    spread_outliers = (90, 94, 189)
    spread_theta, spread_clock, spread_source = arrangement_witness(spread_outliers)
    old_spread_reciprocal = F(1, 90) + F(4, 94)
    new_spread_reciprocal = F(1, 90) + F(2, 94) + F(4, 189)
    require(old_spread_reciprocal == F(227, 4230), "old spread control changed")
    require(old_spread_reciprocal > F(3, 56), "old spread gate unexpectedly passes")
    require(new_spread_reciprocal == F(4757, 88830), "new spread control changed")
    require(F(3, 56) - new_spread_reciprocal == F(1, 50760), "spread gain changed")
    require(
        (spread_theta, spread_clock, spread_source) == (F(11, 49), 98, "core-left"),
        "spread arrangement witness changed",
    )

    # Same scalar weight and same transitive rank order do not determine
    # antipodal feasibility or an AP maximum-deletion core.
    require(weighted_total(OBSTRUCTED_ELEVEN) == 17, "obstructed weight changed")
    require(weighted_total(SAFE_ELEVEN) == 17, "safe weight changed")
    require(set(D_STAR).issubset(OBSTRUCTED_ELEVEN), "D* containment changed")
    obstructed_points = exhaustive_antipodal_safe_points(OBSTRUCTED_ELEVEN)
    safe_points = exhaustive_antipodal_safe_points(SAFE_ELEVEN)
    require(obstructed_points == (), "claimed eleven-core obstruction has a survivor")
    require(F(18, 77) in safe_points, "strong deletion witness disappeared")
    require(all(literal_safe(v, F(18, 77)) for v in SAFE_ELEVEN), "safe hostile witness failed")

    for family in (AP_COVER14, NONAP_COVER14):
        require(is_primitive(family), "Cover14 hostile lost primitivity")
        require(is_cover14(family), "Cover14 hostile lost divisor coverage")
        require(max(family) == 182, "Cover14 hostile maximum changed")
        require(weighted_total(family) == 19, "Cover14 hostile weight changed")
    require(is_arithmetic_progression(AP_COVER14[:-1]), "AP deletion core changed")
    require(not is_arithmetic_progression(NONAP_COVER14[:-1]), "non-AP deletion core became AP")

    ap_maximum = exact_lonely_maximum(AP_COVER14)
    nonap_maximum = exact_lonely_maximum(NONAP_COVER14)
    require(
        ap_maximum[:3] == (F(14, 183), F(14, 183), (1, 182)),
        "AP Cover14 exact maximum changed",
    )
    require(ap_maximum[3] == 2198, "AP Cover14 candidate census changed")
    require(
        nonap_maximum[:3] == (F(2, 17), F(6, 17), (6, 11)),
        "non-AP Cover14 exact maximum changed",
    )
    require(nonap_maximum[3] == 2187, "non-AP Cover14 candidate census changed")

    rank11 = transitive_order_fingerprint(OBSTRUCTED_ELEVEN)
    rank13 = transitive_order_fingerprint(AP_COVER14)
    require(rank11 == transitive_order_fingerprint(SAFE_ELEVEN), "rank-11 order fibre changed")
    require(rank13 == transitive_order_fingerprint(NONAP_COVER14), "rank-13 order fibre changed")
    require((len(rank11), len(rank13)) == (55, 78), "order edge census changed")

    print("THM-4100 RESIDUAL-COMPONENT THREE-OUTLIER EXACT REFEREE")
    print("status=PASS")
    print(f"pair_rows={pair_rows}")
    print(f"pair_circular_literal_components={pair_component_count}")
    print(
        "maximum_pair_min_bound_ratio="
        f"{pair_maximum_ratio};row=(16,41);component=[139/574,74/287];"
        "length=9/574;bound=73/4592"
    )
    print(f"component_triples={component_rows}")
    print(f"circular_literal_components={component_count}")
    print("seam_double_count_scratch_value=1410020;corrected_by=19600")
    print(f"sharpened_pair_span_branch_rows={span_branch_rows}")
    print(
        "maximum_component_ratio="
        f"{maximum_ratio};row=(1,42,47);component=[279/658,379/658];"
        "length=50/329;bound=23/147"
    )
    print(f"ap8_threshold_rows={ap8_rows}")
    print(f"odd_owner_classes={owner_classes}")
    print(f"maximum_clock={maximum_clock}")
    print(f"witness_source_sides=left:{source_left},right:{source_right}")
    print(
        "boundary_row=(93,94,95);theta=11/49;N=98;r=22;half=71;"
        "A_N=34;E_N=0;reciprocal_slack=65/1713432"
    )
    print("uniform_92_criterion_excess=37/119784")
    print(
        "spread_gain_row=(90,94,189);old_reciprocal=227/4230>3/56;"
        "new_reciprocal=4757/88830;slack=1/50760;theta=11/49;N=98"
    )
    print(
        "height12_scalar_hostile=W17/order55;"
        "obstructed=[12]\\{10};safe=[12]\\{4};safe_theta=18/77"
    )
    print(
        "cover14_scalar_hostile=W19/order78/max182;"
        "AP_M=14/183@14/183 owners(1,182) candidates2198;"
        "nonAP_M=2/17@6/17 owners(6,11) candidates2187"
    )
    print("semantic_source=closed antipodal-safe interval plus three sorted outliers")
    print("semantic_map=group the last two danger combs, then reinsert the slowest tooth comb")
    print("semantic_preserved=weak two-phase safety and complementary half-label packet")
    print("semantic_destroyed_by_scalarization=component ancestry, metric length, endpoint owner, AP support")
    print("semantic_sidecar=ordered residual-component envelope plus parity safe-gap and rational endpoint")
    print("typing_firewall=bare order tournament is transitive scheduling data only (THM-4088)")
    print("scope=structured sufficient compiler; arbitrary-core AP extraction and LRC(14) remain OPEN")


if __name__ == "__main__":
    main()
