#!/usr/bin/env python3
"""Primary Fraction-exact component/owner-clock audit for THM-4087.

The proof carrier is the literal open antipodal danger union.  A speed d has
d teeth when d is even and 2d teeth when d is odd, always of radius 1/(14d).
Touching open teeth are deliberately *not* merged.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction


DELTA = Fraction(1, 14)
CORE = tuple(range(1, 10))
CORE_LEFT = Fraction(4, 49)
CORE_RIGHT = Fraction(3, 35)
CORE_INTERVAL = (CORE_LEFT, CORE_RIGHT)
H8 = (1, 3, 5, 8, 11, 13, 23, 36)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def literal_safe(speed: int, theta: Fraction, threshold: Fraction = DELTA) -> bool:
    return all(
        circle_norm(speed * (theta + Fraction(half, 2))) >= threshold
        for half in (0, 1)
    )


def root_count(speed: int) -> int:
    return speed if speed % 2 == 0 else 2 * speed


def tooth_radius(speed: int) -> Fraction:
    return DELTA / speed


def split_open_teeth(speed: int) -> list[tuple[Fraction, Fraction]]:
    """Open antipodal-danger teeth split at 0 on [0,1]."""
    count = root_count(speed)
    radius = tooth_radius(speed)
    pieces: list[tuple[Fraction, Fraction]] = []
    for index in range(count):
        centre = Fraction(index, count)
        start = centre - radius
        end = centre + radius
        if start < 0:
            pieces.append((Fraction(0), end))
            pieces.append((1 + start, Fraction(1)))
        elif end > 1:
            pieces.append((start, Fraction(1)))
            pieces.append((Fraction(0), end - 1))
        else:
            pieces.append((start, end))
    return pieces


def component_lengths(first: int, second: int) -> tuple[Fraction, ...]:
    pieces = sorted(split_open_teeth(first) + split_open_teeth(second))
    merged: list[tuple[Fraction, Fraction]] = []
    for start, end in pieces:
        # Strict comparison is load-bearing: equality is a safe handoff seam.
        if merged and start < merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    lengths = [end - start for start, end in merged]
    if len(merged) >= 2 and merged[0][0] == 0 and merged[-1][1] == 1:
        lengths = [end - start for start, end in merged[1:-1]] + [
            merged[0][1] + 1 - merged[-1][0]
        ]
    return tuple(sorted(lengths))


def linear_open_components(first: int, second: int) -> tuple[tuple[Fraction, Fraction], ...]:
    """Literal components on [0,1], before identifying the two boundary copies."""
    pieces = sorted(split_open_teeth(first) + split_open_teeth(second))
    merged: list[tuple[Fraction, Fraction]] = []
    for start, end in pieces:
        if merged and start < merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


def local_teeth(speed: int) -> tuple[tuple[Fraction, Fraction, int], ...]:
    left, right = CORE_INTERVAL
    count = root_count(speed)
    radius = tooth_radius(speed)
    scaled_left = count * (left - radius)
    scaled_right = count * (right + radius)
    first = scaled_left.numerator // scaled_left.denominator - 2
    last = scaled_right.numerator // scaled_right.denominator + 3
    return tuple(
        (Fraction(index, count) - radius, Fraction(index, count) + radius, speed)
        for index in range(first, last + 1)
        if Fraction(index, count) + radius > left
        and Fraction(index, count) - radius < right
    )


def arrangement_witness(first: int, second: int) -> tuple[Fraction, int, str]:
    candidates: list[tuple[Fraction, int, str]] = [
        (CORE_LEFT, 98, "core_left"),
        (CORE_RIGHT, 70, "core_right"),
    ]
    for speed, name in ((first, "first"), (second, "second")):
        for start, end, _ in local_teeth(speed):
            if CORE_LEFT <= start <= CORE_RIGHT:
                candidates.append((start, 14 * speed, name))
            if CORE_LEFT <= end <= CORE_RIGHT:
                candidates.append((end, 14 * speed, name))
    for theta, clock, source in candidates:
        if literal_safe(first, theta) and literal_safe(second, theta):
            return theta, clock, source
    raise RuntimeError(f"no arrangement witness for {first},{second}")


def mod_norm(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def transformed_safe(speed: int, time: Fraction) -> bool:
    if speed % 2:
        return circle_norm(speed * time) >= Fraction(1, 7)
    return circle_norm((speed // 2) * time) >= DELTA


def main() -> None:
    require(CORE_RIGHT - CORE_LEFT == Fraction(1, 245), "core interval length changed")
    core_phase_gates = 0
    for theta in (CORE_LEFT, Fraction(1, 12), CORE_RIGHT):
        for speed in CORE:
            for half in (0, 1):
                require(
                    circle_norm(speed * (theta + Fraction(half, 2))) >= DELTA,
                    f"core interval checkpoint failed d={speed},theta={theta},half={half}",
                )
                core_phase_gates += 1
    require(circle_norm(7 * (CORE_LEFT + Fraction(1, 2))) == DELTA, "left owner changed")
    require(circle_norm(5 * (CORE_RIGHT + Fraction(1, 2))) == DELTA, "right owner changed")
    epsilon = Fraction(1, 10**8)
    require(not literal_safe(7, CORE_LEFT - epsilon), "left boundary not maximal")
    require(not literal_safe(5, CORE_RIGHT + epsilon), "right boundary not maximal")

    component_pairs = 0
    component_count = 0
    worst_ratio = Fraction(0)
    worst_pair = (0, 0)
    worst_length = Fraction(0)
    for first in range(1, 81):
        for second in range(first + 1, 81):
            lengths = component_lengths(first, second)
            bound = Fraction(2, 7 * first)
            require(lengths and lengths[-1] < bound, f"component bound failed {first},{second}")
            ratio = lengths[-1] / bound
            if ratio > worst_ratio:
                worst_ratio = ratio
                worst_pair = (first, second)
                worst_length = lengths[-1]
            component_pairs += 1
            component_count += len(lengths)

    arithmetic_pairs = 0
    one_tooth_regime = 0
    odd_multi_regime = 0
    even_multi_regime = 0
    for first in range(1, 1001):
        first_length = Fraction(1, 7 * first)
        first_gap = Fraction(5, 14 * first) if first % 2 else Fraction(6, 7 * first)
        for second in range(first + 1, 1001):
            second_length = Fraction(1, 7 * second)
            require(second_length < first_gap, "one lower-speed tooth can be bridged")
            if second % 2:
                if 2 * second < 5 * first:
                    require(
                        Fraction(1, 2 * second) > first_length + second_length,
                        "odd one-tooth spacing branch failed",
                    )
                    require(first_length + second_length < Fraction(2, 7 * first), "odd single bound failed")
                    one_tooth_regime += 1
                else:
                    require(first_length + 2 * second_length < Fraction(2, 7 * first), "odd multi bound failed")
                    odd_multi_regime += 1
            else:
                if second < 6 * first:
                    require(
                        Fraction(1, second) > first_length + second_length,
                        "even one-tooth spacing branch failed",
                    )
                    require(first_length + second_length < Fraction(2, 7 * first), "even single bound failed")
                    one_tooth_regime += 1
                else:
                    require(first_length + 2 * second_length < Fraction(2, 7 * first), "even multi bound failed")
                    even_multi_regime += 1
            arithmetic_pairs += 1

    family_pairs = 0
    family_phase_gates = 0
    odd_class_gates = 0
    source_histogram = {"core_left": 0, "core_right": 0, "first": 0, "second": 0}
    max_clock = 0
    sample_rows: list[tuple[int, int, str, str, int, int]] = []
    witness_bank: dict[tuple[int, int], tuple[Fraction, int, int]] = {}
    for first in range(70, 181):
        for second in range(first + 1, 181):
            theta, clock, source = arrangement_witness(first, second)
            require(clock % 2 == 0 and clock <= 14 * second, "adaptive clock bound failed")
            scaled = theta * clock
            require(scaled.denominator == 1, "witness is not on its advertised clock")
            label = scaled.numerator % clock
            body = CORE + (first, second)
            for speed in body:
                for half in (0, 1):
                    require(
                        circle_norm(speed * (theta + Fraction(half, 2))) >= DELTA,
                        f"family phase failed B={first},C={second},d={speed},half={half}",
                    )
                    family_phase_gates += 1
            opposite = (label + clock // 2) % clock
            for tail in range(1, 2 * clock, 2):
                require(
                    not (
                        7 * mod_norm(tail * label, clock) < clock
                        and 7 * mod_norm(tail * opposite, clock) < clock
                    ),
                    f"eligible odd class survived B={first},C={second},z={tail}",
                )
                odd_class_gates += 1
            source_histogram[source] += 1
            max_clock = max(max_clock, clock)
            witness_bank[(first, second)] = (theta, clock, label)
            if (first, second) in ((70, 71), (70, 180), (98, 147), (179, 180)):
                sample_rows.append((first, second, str(theta), source, clock, label))
            family_pairs += 1

    dilation_phase_gates = 0
    dilation_clock_gates = 0
    for first, second in ((70, 71), (98, 147), (120, 179), (179, 180)):
        theta, clock, label = witness_bank[(first, second)]
        for dilation in range(1, 16):
            lifted = theta / dilation
            dilated_clock = dilation * clock
            require(dilated_clock % 2 == 0, "dilated clock parity failed")
            require(Fraction(label, dilated_clock) == lifted % 1, "dilated clock phase failed")
            for speed in CORE + (first, second):
                for half in (0, 1):
                    require(
                        circle_norm(dilation * speed * (lifted + Fraction(half, 2))) >= DELTA,
                        f"dilation failed q={dilation},d={speed}",
                    )
                    dilation_phase_gates += 1
            opposite = (label + dilated_clock // 2) % dilated_clock
            for tail in (1, 3, 5, 7, 11, 13):
                require(
                    not (
                        7 * mod_norm(tail * label, dilated_clock) < dilated_clock
                        and 7 * mod_norm(tail * opposite, dilated_clock) < dilated_clock
                    ),
                    "dilated eligibility hostile survived",
                )
                dilation_clock_gates += 1

    hostile_gates = 0
    h8_points = tuple(Fraction(index, 7) for index in range(1, 7))
    for multiple in range(1, 101):
        appended = 7 * multiple
        for time in h8_points:
            require(all(transformed_safe(speed, time) for speed in H8), "H8 endpoint changed")
            require(not transformed_safe(appended, time), "7m failed to kill an H8 endpoint")
            hostile_gates += 1

    hostile_left = Fraction(139, 1722)
    hostile_right = Fraction(74, 861)
    hostile_centre = Fraction(1, 12)
    hostile_component = next(
        (start, end)
        for start, end in linear_open_components(48, 123)
        if start < hostile_centre < end
    )
    require(
        hostile_component == (hostile_left, hostile_right),
        "48/123 literal-open component changed",
    )
    require(hostile_left < CORE_LEFT < CORE_RIGHT < hostile_right, "48/123 local cover ordering changed")
    for theta in (CORE_LEFT, Fraction(1, 12), CORE_RIGHT):
        require(
            not literal_safe(48, theta) or not literal_safe(123, theta),
            "48/123 failed to cover a core checkpoint",
        )

    semantic = {
        "arithmetic_regimes": [one_tooth_regime, odd_multi_regime, even_multi_regime],
        "component_worst": [*worst_pair, str(worst_length), str(worst_ratio)],
        "core_interval": [str(CORE_LEFT), str(CORE_RIGHT)],
        "family_samples": sample_rows,
        "source_histogram": source_histogram,
    }
    digest = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4087 primary literal-open component/owner-clock audit")
    print(f"core_interval=[{CORE_LEFT},{CORE_RIGHT}] length={CORE_RIGHT - CORE_LEFT} core_phase_gates={core_phase_gates}")
    print(f"component_pairs={component_pairs} component_count={component_count} worst_pair={worst_pair[0]},{worst_pair[1]} worst_length={worst_length} bound_ratio={worst_ratio}")
    print(f"arithmetic_pairs={arithmetic_pairs} one_tooth={one_tooth_regime} odd_multi={odd_multi_regime} even_multi={even_multi_regime}")
    print(f"family_pairs={family_pairs} family_phase_gates={family_phase_gates} odd_class_gates={odd_class_gates} max_clock={max_clock}")
    print("witness_sources=" + ",".join(f"{key}:{source_histogram[key]}" for key in sorted(source_histogram)))
    print(f"dilation_phase_gates={dilation_phase_gates} dilation_clock_gates={dilation_clock_gates}")
    print(f"H8_7m_hostile_gates={hostile_gates} local_48_123_component=[{hostile_left},{hostile_right}]")
    print(f"semantic_sha256={digest}")
    print("PASS: two-comb interval compiler, adaptive owner clock, AP9 two-outlier family, and dilation")


if __name__ == "__main__":
    main()
