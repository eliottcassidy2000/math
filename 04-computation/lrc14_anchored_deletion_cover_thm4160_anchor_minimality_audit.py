#!/usr/bin/env python3
"""Exact finite anchor-minimality sidecar for THM-4160."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb, gcd


Q = Fraction
DELTA = Q(1, 14)
THRESHOLD = Q(4, 63)
ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
OPTIONAL = frozenset(POOL) - frozenset(ANCHORS)


def require(predicate: bool, label: object) -> None:
    if not predicate:
        raise RuntimeError(f"requirement failed: {label}")


def qfloor(value: Q) -> int:
    return value.numerator // value.denominator


def intersect_one(intervals: tuple[tuple[Q, Q], ...], speed: int):
    output = []
    for left, right in intervals:
        low = max(0, qfloor(speed * left) - 1)
        high = min(speed - 1, qfloor(speed * right) + 1)
        for tooth in range(low, high + 1):
            a = max(left, Q(14 * tooth + 1, 14 * speed))
            b = min(right, Q(14 * tooth + 13, 14 * speed))
            if a <= b:
                if output and output[-1][1] == a:
                    output[-1] = (output[-1][0], b)
                else:
                    output.append((a, b))
    return tuple(output)


def components(speeds: tuple[int, ...]):
    answer = ((Q(0), Q(1)),)
    for speed in speeds:
        answer = intersect_one(answer, speed)
    return answer


def measure(intervals) -> Q:
    return sum((right - left for left, right in intervals), Q(0))


def main() -> None:
    covers = []
    for triple in combinations(range(1, 144), 3):
        if gcd(gcd(triple[0], triple[1]), triple[2]) != 1:
            continue
        if all(
            any(value % modulus == 0 for value in triple)
            for modulus in range(2, 15)
        ):
            covers.append(triple)
    expected_covers = [
        (70, 72, 143),
        (72, 140, 143),
        (120, 126, 143),
    ]
    require(covers == expected_covers, "primitive divisor-complete triples")
    require(all(triple[-1] == 143 for triple in covers), "maximum 143 is forced")

    geometry = []
    for triple in covers:
        candidate_pool = tuple(sorted(OPTIONAL | frozenset(triple)))
        intervals = components(candidate_pool)
        mass = measure(intervals)
        geometry.append((triple, len(candidate_pool), len(intervals), mass, mass - THRESHOLD))
    expected_geometry = [
        ((70, 72, 143), 30, 102, Q(3242643641, 80005085160),
         Q(3242643641, 80005085160) - THRESHOLD),
        ((72, 140, 143), 30, 106, Q(90567155849, 2280144927060),
         Q(90567155849, 2280144927060) - THRESHOLD),
        ((120, 126, 143), 30, 150, Q(298133356159, 4560289854120),
         Q(298133356159, 4560289854120) - THRESHOLD),
    ]
    require(geometry == expected_geometry, "fixed-optional-pool geometry")
    require(tuple(mass >= THRESHOLD for _, _, _, mass, _ in geometry) ==
            (False, False, True), "unique Haar-compatible anchor triple")

    ownership = {
        modulus: tuple(anchor for anchor in ANCHORS if anchor % modulus == 0)
        for modulus in range(2, 15)
    }
    for anchor in ANCHORS:
        require(any(
            ownership[modulus] == (anchor,) for modulus in ownership
        ), ("anchor inclusion-minimality", anchor))

    semantic = {
        "universe": (1, 143, comb(143, 3)),
        "covers": covers,
        "geometry": [
            (triple, pool_size, component_count, str(mass), str(margin))
            for triple, pool_size, component_count, mass, margin in geometry
        ],
        "ownership": ownership,
    }
    print("LRC14_ANCHORED_DELETION_COVER_THM4160_ANCHOR_AUDIT_20260826")
    print(f"universe=1<=a<b<c<=143;checked={comb(143,3)};"
          "constraints=gcd1_and_each_d_2_through_14_divides_an_anchor")
    print(f"covers={tuple(covers)};all_maxima_143=True")
    print(f"divisor_ownership={ownership};current_anchor_inclusion_minimal=True")
    print("fixed_optional_pool_geometry=(triple,pool_size,components,mass,margin)")
    for row in geometry:
        print(f"  {row}")
    print("haar_threshold_pattern=(False,False,True);"
          "current_anchors_unique_in_complete_height143_cover_universe=True")
    print("semantic_sha256=" + sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest())


if __name__ == "__main__":
    main()
