#!/usr/bin/env python3
"""Exact primary audit for THM-4156.

The certificate deliberately keeps two claims separate:

* the complete closed safe-set arrangement of the 30-label pool has Haar
  measure strictly above THM-4150's sharp ``4/63`` transfer threshold;
* the three mandatory anchors cover every divisor from 2 through 14, so the
  hereditary family lies inside the genuine divisor-complete dyadic seam.

All theorem-bearing arithmetic uses ``Fraction``.  The combination census is
complete over its explicitly printed universe.
"""

from __future__ import annotations

import ast
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

DELTA = Q(1, 14)
THRESHOLD = Q(4, 63)
ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
EXPECTED_MEASURE = Q(298133356159, 4560289854120)
EXPECTED_SURPLUS = Q(8591143199, 4560289854120)
EXPECTED_MAX_WIDTH = Q(37, 25520)
EXPECTED_SEMANTIC = "ed193657974754e64ae4faed03bc63146f453c16290553ef5d88f198bc77bd0e"

CHECKS = 0


def require(predicate: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not predicate:
        raise RuntimeError(f"requirement failed: {label}")


def pair(value: Q) -> tuple[int, int]:
    return value.numerator, value.denominator


def digest(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(raw).hexdigest()


def circle_distance(value: Q) -> Q:
    residue = value % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Q) -> Q:
    return min(circle_distance(speed * phase) for speed in speeds)


def safe_walls(speeds: tuple[int, ...]) -> tuple[Q, ...]:
    walls = {Q(0), Q(1)}
    for speed in speeds:
        for integer in range(speed):
            walls.add((Q(integer) + DELTA) / speed)
            walls.add((Q(integer + 1) - DELTA) / speed)
    return tuple(sorted(walls))


def safe_components(speeds: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    walls = safe_walls(speeds)
    cells: list[tuple[Q, Q]] = []
    for left, right in zip(walls, walls[1:]):
        middle = (left + right) / 2
        if clearance(speeds, middle) < DELTA:
            continue
        require(clearance(speeds, left) >= DELTA, ("safe left", left, right))
        require(clearance(speeds, right) >= DELTA, ("safe right", left, right))
        cells.append((left, right))

    result: list[tuple[Q, Q]] = []
    for left, right in cells:
        if result and result[-1][1] == left:
            result[-1] = (result[-1][0], right)
        else:
            result.append((left, right))
    return tuple(result)


def divisor_audit() -> tuple[tuple[int, int], ...]:
    require(len(ANCHORS) == len(set(ANCHORS)) == 3, "anchor cardinality")
    require(set(ANCHORS) < set(POOL), "anchors are a proper pool subset")
    owners = []
    for modulus in range(2, 15):
        owner = next((anchor for anchor in ANCHORS if anchor % modulus == 0), None)
        require(owner is not None, ("divisor owner", modulus))
        owners.append((modulus, owner))

    # Inclusion-minimality is a useful hostile control for the anchor design.
    for omitted in ANCHORS:
        retained = tuple(value for value in ANCHORS if value != omitted)
        missed = tuple(
            modulus for modulus in range(2, 15)
            if not any(value % modulus == 0 for value in retained)
        )
        require(bool(missed), ("anchor deletion loses a divisor", omitted))

    # The repaired construction must not fall back into MISTAKE-511's
    # no-multiple-of-six branch.
    require(all(anchor % 6 == 0 for anchor in ANCHORS[:2]), "multiple-six anchors")
    require(clearance(tuple(2 * anchor for anchor in ANCHORS), Q(1, 12)) == 0,
            "one-twelfth clock hostile")
    return tuple(owners)


def geometry_audit() -> dict[str, object]:
    require(len(POOL) == len(set(POOL)) == 30, "pool cardinality")
    require(tuple(sorted(POOL)) == POOL, "pool ordering")
    components = safe_components(POOL)
    walls = safe_walls(POOL)
    measure = sum((right - left for left, right in components), Q(0))
    maximum = max(right - left for left, right in components)
    require(len(walls) == 7134, "distinct wall count")
    require(len(components) == 150, "component count")
    require(measure == EXPECTED_MEASURE, "exact Haar measure")
    require(measure - THRESHOLD == EXPECTED_SURPLUS > 0, "strict Haar surplus")
    require(maximum == EXPECTED_MAX_WIDTH, "largest component")
    require(tuple((1 - right, 1 - left) for left, right in reversed(components)) == components,
            "reflection symmetry")

    # Direct positive control on the full pool, stronger than one selected
    # eleven-body.  It is not used in the all-tail proof.
    y = Q(199, 21280)
    x = Q(21479, 42560)
    require(clearance(POOL, y) == Q(199, 2660), "positive quotient control")
    require(clearance(tuple(2 * speed for speed in POOL), x) == Q(199, 2660),
            "positive doubled-pool control")
    require(clearance((1, 3), x) == Q(20683, 42560), "positive tail control")
    return {
        "walls": len(walls),
        "components": len(components),
        "measure": pair(measure),
        "surplus": pair(measure - THRESHOLD),
        "maximum_width": pair(maximum),
        "positive_control": (pair(y), pair(x), pair(Q(199, 2660)), pair(Q(20683, 42560))),
    }


def family_audit() -> dict[str, int]:
    optional = tuple(value for value in POOL if value not in ANCHORS)
    require(len(optional) == 27, "optional-label count")
    total = comb(len(optional), 8)
    require(total == 2_220_075, "hereditary family count")

    thm4148 = 0
    thm4151 = 0
    both = 0
    checked = 0
    for choice in combinations(optional, 8):
        body = tuple(sorted(ANCHORS + choice))
        require(len(body) == len(set(body)) == 11, ("body cardinality", choice))
        minimum = body[0]
        maximum = body[-1]
        gate4148 = 27 * (13 * minimum - maximum) >= 4 * minimum * maximum
        gate4151 = 16 * maximum <= 156 * minimum + 13
        thm4148 += int(gate4148)
        thm4151 += int(gate4151)
        both += int(gate4148 and gate4151)
        checked += 1

    require(checked == total, "complete family census")
    require(thm4148 == 0, "all bodies evade THM-4148 gate")
    require(thm4151 == 344_366, "THM-4151 overlap")
    require(both == 0, "first-window overlap intersection")
    beyond = total - thm4151
    require(beyond == 1_875_709, "beyond both first-window gates")
    return {
        "total": total,
        "thm4148": thm4148,
        "thm4151": thm4151,
        "both": both,
        "beyond_union": beyond,
    }


def main() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "no asserts")
    owners = divisor_audit()
    geometry = geometry_audit()
    families = family_audit()
    semantic = {
        "delta": pair(DELTA),
        "threshold": pair(THRESHOLD),
        "anchors": ANCHORS,
        "pool": POOL,
        "divisor_owners": owners,
        "geometry": geometry,
        "families": families,
    }
    semantic_hash = digest(semantic)
    if EXPECTED_SEMANTIC != "TO_BE_FILLED":
        require(semantic_hash == EXPECTED_SEMANTIC, "frozen semantic digest")

    print("THM4156_DIVISOR_COMPLETE_ANCHOR_POOL_HAAR_TRANSFER_20260825")
    print("status=PASS;scope=exact common-safe-set and complete stated-family census")
    print(f"anchors={ANCHORS};divisor_owners={owners}")
    print(f"pool_size={len(POOL)};pool={POOL}")
    print(f"pool_walls={geometry['walls']};pool_components={geometry['components']}")
    print(f"pool_measure={geometry['measure'][0]}/{geometry['measure'][1]}")
    print(f"threshold=4/63;surplus={geometry['surplus'][0]}/{geometry['surplus'][1]}")
    print(f"pool_max_component={geometry['maximum_width'][0]}/{geometry['maximum_width'][1]}")
    print(f"families={families}")
    print("fixed_clock_hostile=x=1/12 has zero clearance at doubled anchor 120")
    print(f"positive_control={geometry['positive_control']}")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={semantic_hash}")


if __name__ == "__main__":
    main()
