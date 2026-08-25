#!/usr/bin/env python3
"""Independent Fraction-exact audit for THM-4100.

This file does not import the primary referee.  Circular open components are
constructed in a different gauge: merge three lifted periods with strict
overlap and retain the unique representative whose midpoint lies in [0,1).
The AP8 clock is recovered from a separately generated arrangement endpoint
bank.  Exact scalar maxima are certified by closed interval covers, rather
than by the primary tent-vertex candidate enumeration.

Every proof gate uses ``require`` and every nonintegral computation uses
``Fraction``.  The transcript is deterministic under normal and optimized
Python.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import gcd


DELTA = F(1, 14)
AP8_LEFT = F(11, 49)
AP8_RIGHT = F(13, 56)
AP8_LENGTH = F(3, 392)

D_STAR = (1, 3, 4, 5, 7, 8, 9, 11, 12)
OBSTRUCTED = tuple(speed for speed in range(1, 13) if speed != 10)
SAFE = tuple(speed for speed in range(1, 13) if speed != 4)
AP_COVER14 = tuple(range(1, 13)) + (182,)
NONAP_COVER14 = tuple(speed for speed in range(1, 14) if speed != 3) + (182,)


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError("FAILED: " + label)


def fraction_text(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def parity_weight(speed: int) -> int:
    return 1 if speed % 2 == 0 else 2


def circle_norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def antipodal_safe(speed: int, theta: F) -> bool:
    return (
        circle_norm(speed * theta) >= DELTA
        and circle_norm(speed * (theta + F(1, 2))) >= DELTA
    )


def mod_norm(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def open_teeth_three_periods(speed: int) -> tuple[tuple[F, F], ...]:
    count = parity_weight(speed) * speed
    radius = F(1, 14 * speed)
    return tuple(
        (F(index, count) - radius, F(index, count) + radius)
        for index in range(-count, 2 * count)
    )


def merge_intervals(
    intervals: list[tuple[F, F]], touching_joins: bool
) -> tuple[tuple[F, F], ...]:
    require(bool(intervals), "nonempty interval family")
    ordered = sorted(intervals)
    merged: list[tuple[F, F]] = []
    active_left, active_right = ordered[0]
    for left, right in ordered[1:]:
        overlaps = left < active_right or (touching_joins and left == active_right)
        if overlaps:
            active_right = max(active_right, right)
        else:
            merged.append((active_left, active_right))
            active_left, active_right = left, right
    merged.append((active_left, active_right))
    return tuple(merged)


def circular_open_components(
    speeds: tuple[int, ...], tooth_bank: dict[int, tuple[tuple[F, F], ...]]
) -> tuple[tuple[F, F, F], ...]:
    """One full lift per circular component, selected by midpoint in [0,1)."""
    intervals = [interval for speed in speeds for interval in tooth_bank[speed]]
    merged = merge_intervals(intervals, touching_joins=False)
    chosen = tuple(
        (left, right, right - left)
        for left, right in merged
        if 0 <= (left + right) / 2 < 1
    )
    require(all(length < 1 for _, _, length in chosen), "proper circular components")
    return chosen


def safe_gap(speed: int) -> F:
    if speed % 2 == 0:
        return F(6, 7 * speed)
    return F(5, 14 * speed)


def pair_envelope(lower: int, upper: int) -> F:
    return min(F(2, 7 * lower), F(1, 7 * lower) + F(2, 7 * upper))


def triple_reciprocal(lower: int, middle: int, upper: int) -> F:
    return F(1, lower) + min(F(4, middle), F(2, middle) + F(4, upper))


def triple_component_bound(lower: int, middle: int, upper: int) -> F:
    return F(1, 7 * lower) + 2 * pair_envelope(middle, upper)


def teeth_meeting_interval(speed: int, left: F, right: F) -> tuple[tuple[F, F], ...]:
    count = parity_weight(speed) * speed
    radius = F(1, 14 * speed)
    scaled_left = count * (left - radius)
    scaled_right = count * (right + radius)
    first = scaled_left.numerator // scaled_left.denominator - 2
    last = scaled_right.numerator // scaled_right.denominator + 2
    return tuple(
        (F(index, count) - radius, F(index, count) + radius)
        for index in range(first, last + 1)
        if F(index, count) + radius > left
        and F(index, count) - radius < right
    )


def arrangement_endpoint(outliers: tuple[int, int, int]) -> F:
    candidates = {AP8_LEFT, AP8_RIGHT}
    for speed in outliers:
        for left, right in teeth_meeting_interval(speed, AP8_LEFT, AP8_RIGHT):
            if AP8_LEFT <= left <= AP8_RIGHT:
                candidates.add(left)
            if AP8_LEFT <= right <= AP8_RIGHT:
                candidates.add(right)
    body = tuple(range(1, 9)) + outliers
    for theta in sorted(candidates):
        if all(antipodal_safe(speed, theta) for speed in body):
            return theta
    raise RuntimeError("no arrangement endpoint for " + str(outliers))


def even_denominator(theta: F) -> int:
    return theta.denominator if theta.denominator % 2 == 0 else 2 * theta.denominator


def integer_safe(speed: int, label: int, clock: int) -> bool:
    return 14 * mod_norm(speed * label, clock) >= clock


def safe_packet(speeds: tuple[int, ...], clock: int) -> tuple[int, ...]:
    return tuple(
        label
        for label in range(clock)
        if all(integer_safe(speed, label, clock) for speed in speeds)
    )


def eligible_odd_classes(packet: tuple[int, ...], clock: int) -> tuple[int, ...]:
    return tuple(
        odd
        for odd in range(1, 2 * clock, 2)
        if all(7 * mod_norm(odd * label, clock) < clock for label in packet)
    )


def weighted_total(speeds: tuple[int, ...]) -> int:
    return sum(parity_weight(speed) for speed in speeds)


def gcd_all(speeds: tuple[int, ...]) -> int:
    common = 0
    for speed in speeds:
        common = gcd(common, speed)
    return common


def cover14(speeds: tuple[int, ...]) -> bool:
    return all(
        any(speed % divisor == 0 for speed in speeds)
        for divisor in range(2, 15)
    )


def is_ap(speeds: tuple[int, ...]) -> bool:
    ordered = tuple(sorted(speeds))
    if len(ordered) <= 2:
        return True
    difference = ordered[1] - ordered[0]
    return difference > 0 and all(
        ordered[index + 1] - ordered[index] == difference
        for index in range(1, len(ordered) - 1)
    )


def open_antipodal_union_covers_circle(speeds: tuple[int, ...]) -> bool:
    intervals = []
    for speed in speeds:
        count = parity_weight(speed) * speed
        radius = F(1, 14 * speed)
        intervals.extend(
            (F(index, count) - radius, F(index, count) + radius)
            for index in range(-count, 2 * count + 1)
        )
    merged = merge_intervals(intervals, touching_joins=False)
    return any(left < 0 and right > 1 for left, right in merged)


def closed_lonely_sublevel_covers_circle(speeds: tuple[int, ...], level: F) -> bool:
    """Certify max_theta min_v ||v theta|| <= level by a closed cover."""
    intervals = []
    for speed in speeds:
        radius = level / speed
        intervals.extend(
            (F(index, speed) - radius, F(index, speed) + radius)
            for index in range(-speed, 2 * speed + 1)
        )
    merged = merge_intervals(intervals, touching_joins=True)
    return any(left <= 0 and right >= 1 for left, right in merged)


def lonely_value(speeds: tuple[int, ...], theta: F) -> F:
    return min(circle_norm(speed * theta) for speed in speeds)


def transformed_support(speeds: tuple[int, ...]) -> set[tuple[int, F]]:
    return {
        (speed if speed % 2 else speed // 2, F(1, 7) if speed % 2 else F(1, 14))
        for speed in speeds
    }


def main() -> None:
    require(AP8_RIGHT - AP8_LEFT == AP8_LENGTH, "AP8 interval length")
    require(all(antipodal_safe(speed, AP8_LEFT) for speed in range(1, 9)),
            "AP8 left endpoint")
    require(all(antipodal_safe(speed, AP8_RIGHT) for speed in range(1, 9)),
            "AP8 right endpoint")

    # A different circular gauge from the primary: midpoint, not left endpoint.
    tooth_bank = {speed: open_teeth_three_periods(speed) for speed in range(1, 51)}
    pair_rows = 0
    pair_components = 0
    pair_seams = 0
    pair_max_ratio = F(0)
    pair_max_row: tuple[object, ...] | None = None
    for lower, upper in combinations(range(1, 51), 2):
        density = F(2, 7 * lower)
        span = F(1, 7 * lower) + F(2, 7 * upper)
        bound = min(density, span)
        require((span < density) == (upper > 2 * lower), "pair branch boundary")
        require(F(1, 7 * upper) < safe_gap(lower), "one fast tooth per slow gap")
        components = circular_open_components((lower, upper), tooth_bank)
        seam_count = sum(left < 0 < right for left, right, _ in components)
        require(seam_count == 1, "one pair seam component")
        pair_seams += seam_count
        pair_components += len(components)
        for left, right, length in components:
            require(length < density, "literal-open pair density bound")
            require(length < span, "literal-open pair protrusion bound")
            ratio = length / bound
            if ratio > pair_max_ratio:
                pair_max_ratio = ratio
                pair_max_row = ((lower, upper), left, right, length, bound)
        pair_rows += 1

    require(pair_rows == 1225, "pair row count")
    require(pair_components == 73066, "pair component count")
    require(pair_seams == pair_rows, "pair seam census")
    require(pair_max_ratio == F(72, 73), "pair maximum ratio")
    require(
        pair_max_row == (
            (16, 41), F(139, 574), F(74, 287), F(9, 574), F(73, 4592)
        ),
        "pair maximum row",
    )

    triple_rows = 0
    triple_components = 0
    triple_seams = 0
    clipped_components = 0
    sharpened_rows = 0
    triple_max_ratio = F(0)
    midpoint_max_row: tuple[object, ...] | None = None
    for lower, middle, upper in combinations(range(1, 51), 3):
        density_pair = F(2, 7 * middle)
        span_pair = F(1, 7 * middle) + F(2, 7 * upper)
        pair_bound = min(density_pair, span_pair)
        sharpened_rows += span_pair < density_pair
        require(pair_bound < safe_gap(lower), "pair component cannot bridge slow teeth")
        bound = F(1, 7 * lower) + 2 * pair_bound
        require(
            bound == F(1, 7) * triple_reciprocal(lower, middle, upper),
            "three-comb envelope algebra",
        )
        components = circular_open_components((lower, middle, upper), tooth_bank)
        seam_count = sum(left < 0 < right for left, right, _ in components)
        require(seam_count == 1, "one triple seam component")
        triple_seams += seam_count
        triple_components += len(components)
        clipped_components += len(components) + seam_count
        for left, right, length in components:
            require(length < bound, "three-comb component bound")
            ratio = length / bound
            if ratio > triple_max_ratio:
                triple_max_ratio = ratio
                midpoint_max_row = (
                    (lower, middle, upper), left, right, length, bound
                )
        triple_rows += 1

    require(triple_rows == 19600, "triple row count")
    require(triple_components == 1390420, "canonical circular component count")
    require(triple_seams == triple_rows, "one seam component per triple")
    require(clipped_components == 1410020, "clipped seam double-count")
    require(sharpened_rows == 4600, "sharpened pair branch rows")
    require(triple_max_ratio == F(1050, 1081), "triple maximum ratio")
    require(
        midpoint_max_row == (
            (1, 42, 47), F(-25, 329), F(25, 329), F(50, 329), F(23, 147)
        ),
        "midpoint-gauge triple maximum",
    )
    primary_tie = (F(279, 658), F(379, 658), F(50, 329))
    require(
        primary_tie in circular_open_components((1, 42, 47), tooth_bank),
        "primary interior maximizer is a tied component",
    )

    # Threshold and arrangement-endpoint/clock audit.
    target = F(3, 56)
    least_uniform = next(
        height
        for height in range(1, 200)
        if triple_reciprocal(height, height + 1, height + 2) <= target
    )
    require(least_uniform == 93, "least uniform criterion height")
    require(
        triple_reciprocal(92, 93, 94) - target == F(37, 119784),
        "height-92 criterion excess",
    )
    require(
        target - (F(1, 93) + F(4, 94)) == F(65, 244776),
        "height-93 reciprocal slack",
    )

    ap8_rows = 0
    odd_classes = 0
    maximum_clock = 0
    distinct_endpoints: set[F] = set()
    for outliers in combinations(range(93, 121), 3):
        require(triple_reciprocal(*outliers) <= target, "AP8 component gate")
        theta = arrangement_endpoint(outliers)
        distinct_endpoints.add(theta)
        clock = even_denominator(theta)
        require(clock % 2 == 0 and clock <= max(98, 14 * outliers[2]),
                "adaptive even clock bound")
        label_value = theta * clock
        require(label_value.denominator == 1, "endpoint clock presentation")
        label = int(label_value)
        half_label = (label + clock // 2) % clock
        body = tuple(range(1, 9)) + outliers
        require(all(integer_safe(speed, label, clock) for speed in body),
                "endpoint label packet")
        require(all(integer_safe(speed, half_label, clock) for speed in body),
                "half endpoint label packet")
        for odd in range(1, 2 * clock, 2):
            first = mod_norm(odd * label, clock)
            second = mod_norm(odd * half_label, clock)
            require(first + second == clock // 2, "odd half-label distance identity")
            require(not (7 * first < clock and 7 * second < clock),
                    "odd class cannot be eligible on complementary labels")
            odd_classes += 1
        maximum_clock = max(maximum_clock, clock)
        ap8_rows += 1

    require(ap8_rows == 3276, "AP8 row count")
    require(odd_classes == 2227526, "odd class count")
    require(maximum_clock == 1680, "maximum adaptive clock")
    require(len(distinct_endpoints) == 19, "arrangement endpoint support")

    boundary = (93, 94, 95)
    boundary_body = tuple(range(1, 9)) + boundary
    boundary_theta = arrangement_endpoint(boundary)
    boundary_clock = even_denominator(boundary_theta)
    boundary_label = int(boundary_theta * boundary_clock)
    boundary_half = (boundary_label + boundary_clock // 2) % boundary_clock
    packet = safe_packet(boundary_body, boundary_clock)
    eligible = eligible_odd_classes(packet, boundary_clock)
    require(boundary_theta == F(11, 49), "boundary endpoint")
    require((boundary_clock, boundary_label, boundary_half) == (98, 22, 71),
            "boundary labels")
    require(len(packet) == 34 and eligible == (), "boundary E_N emptiness")
    require(
        AP8_LENGTH - triple_component_bound(*boundary) == F(65, 1713432),
        "boundary length slack",
    )
    for dilation in range(1, 9):
        transported = tuple(dilation * speed for speed in boundary_body)
        require(all(antipodal_safe(speed, boundary_theta / dilation)
                    for speed in transported), "common dilation transport")
        require((boundary_theta / dilation) * (dilation * boundary_clock)
                == boundary_label, "dilated clock presentation")

    spread = (90, 94, 189)
    old_spread = F(1, 90) + F(4, 94)
    new_spread = triple_reciprocal(*spread)
    spread_theta = arrangement_endpoint(spread)
    require(old_spread == F(227, 4230) and old_spread > target,
            "old spread envelope fails")
    require(new_spread == F(4757, 88830), "new spread envelope")
    require(target - new_spread == F(1, 50760), "spread gain")
    require((spread_theta, even_denominator(spread_theta)) == (F(11, 49), 98),
            "spread endpoint and clock")

    # D*/safe hostile, including the exact transformed-support explanation.
    require(set(D_STAR).issubset(OBSTRUCTED), "D* containment")
    require(open_antipodal_union_covers_circle(D_STAR), "D* exact obstruction")
    require(open_antipodal_union_covers_circle(OBSTRUCTED), "obstructed eleven-core")
    require(not open_antipodal_union_covers_circle(SAFE), "safe hostile not obstructed")
    safe_theta = F(18, 77)
    require(all(antipodal_safe(speed, safe_theta) for speed in SAFE),
            "safe hostile witness")
    require(weighted_total(OBSTRUCTED) == weighted_total(SAFE) == 17,
            "height-12 hostile weight")
    obstructed_support = transformed_support(OBSTRUCTED)
    safe_support = transformed_support(SAFE)
    require((2, F(1, 14)) in obstructed_support
            and (2, F(1, 14)) not in safe_support,
            "new transformed multiplier two")
    require((5, F(1, 7)) in safe_support and (5, F(1, 14)) in safe_support,
            "speed ten is a weaker shadow of odd speed five")

    # Cover14 scalar hostile: exact maxima by closed sublevel covers.
    for family in (AP_COVER14, NONAP_COVER14):
        require(gcd_all(family) == 1, "Cover14 primitivity")
        require(cover14(family), "Cover14 divisor coverage")
        require(max(family) == 182, "Cover14 maximum")
        require(weighted_total(family) == 19, "Cover14 parity weight")
    require(is_ap(AP_COVER14[:-1]), "AP maximum deletion")
    require(not is_ap(NONAP_COVER14[:-1]), "non-AP maximum deletion")

    ap_level, ap_theta = F(14, 183), F(14, 183)
    nonap_level, nonap_theta = F(2, 17), F(6, 17)
    require(closed_lonely_sublevel_covers_circle(AP_COVER14, ap_level),
            "AP exact maximum upper cover")
    require(lonely_value(AP_COVER14, ap_theta) == ap_level,
            "AP exact maximum witness")
    require(tuple(speed for speed in AP_COVER14
                  if circle_norm(speed * ap_theta) == ap_level) == (1, 182),
            "AP maximum owners")
    require(closed_lonely_sublevel_covers_circle(NONAP_COVER14, nonap_level),
            "non-AP exact maximum upper cover")
    require(lonely_value(NONAP_COVER14, nonap_theta) == nonap_level,
            "non-AP exact maximum witness")
    require(tuple(speed for speed in NONAP_COVER14
                  if circle_norm(speed * nonap_theta) == nonap_level) == (6, 11),
            "non-AP maximum owners")
    require(ap_level < F(1, 13) < nonap_level, "deep-residual scope separation")
    require((len(OBSTRUCTED) * (len(OBSTRUCTED) - 1) // 2,
             len(AP_COVER14) * (len(AP_COVER14) - 1) // 2) == (55, 78),
            "transitive rank fingerprints")

    scope = {
        "THM-4100": {
            "base": "AP8",
            "outliers": 3,
            "parity": "arbitrary",
            "uniform_height": 93,
            "carrier": "literal component ancestry",
        },
        "THM-4101": {
            "base": "AP7",
            "outliers": 4,
            "parity": "three odd and one even",
            "resonance": "selected odd gap-four pair",
            "uniform_height": 264,
            "carrier": "second-moment overlap credit",
        },
    }
    require(scope["THM-4100"]["outliers"] != scope["THM-4101"]["outliers"],
            "scope arity separation")
    require(scope["THM-4100"]["parity"] == "arbitrary"
            and "gap-four" in scope["THM-4101"]["resonance"],
            "scope hypothesis separation")

    pair_record = {
        "rows": pair_rows,
        "components": pair_components,
        "seams": pair_seams,
        "max_ratio": fraction_text(pair_max_ratio),
        "max_row": [16, 41],
        "max_component": ["139/574", "74/287"],
        "max_bound": "73/4592",
    }
    triple_record = {
        "rows": triple_rows,
        "components": triple_components,
        "seams": triple_seams,
        "clipped_components": clipped_components,
        "sharpened_rows": sharpened_rows,
        "max_ratio": fraction_text(triple_max_ratio),
        "max_row": [1, 42, 47],
        "midpoint_component": ["-25/329", "25/329"],
        "primary_tied_component": ["279/658", "379/658"],
    }
    ap_record = {
        "least_uniform": least_uniform,
        "rows": ap8_rows,
        "odd_classes": odd_classes,
        "maximum_clock": maximum_clock,
        "distinct_endpoints": len(distinct_endpoints),
        "boundary": ["11/49", 98, 22, 71, len(packet), len(eligible)],
        "spread": [90, 94, 189, "4757/88830", "1/50760", "11/49", 98],
    }
    hostiles = {
        "Dstar_safe": [
            17,
            "18/77",
            [[m, fraction_text(theta)] for m, theta in sorted(obstructed_support)],
            [[m, fraction_text(theta)] for m, theta in sorted(safe_support)],
        ],
        "Cover14": [
            "14/183", "14/183", [1, 182],
            "2/17", "6/17", [6, 11],
        ],
    }
    ledger = {
        "theorem": "THM-4100",
        "component_gauge": "one full lift with midpoint in [0,1)",
        "pair": pair_record,
        "triple": triple_record,
        "AP8": ap_record,
        "hostiles": hostiles,
        "scope": scope,
    }
    semantic = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("THM-4100 INDEPENDENT RESIDUAL-COMPONENT AUDIT")
    print("pair=", pair_record)
    print("triple=", triple_record)
    print("seam_convention=midpoint lift; one seam component per row; clipping adds one")
    print("triple_max_tie=midpoint gauge selects seam [-25/329,25/329];"
          " primary left-end gauge selects interior [279/658,379/658]")
    print("AP8=", ap_record)
    print("Dstar_safe_hostile=weight17;Dstar_and_O_cover;S_theta=18/77")
    print("Cover14_hostile=weight19;max182;AP_M=14/183;nonAP_M=2/17")
    print("scope=THM-4100 AP8+3 arbitrary parity versus THM-4101 AP7+4 resonant weight7")
    print("semantic_sha256=", semantic)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
