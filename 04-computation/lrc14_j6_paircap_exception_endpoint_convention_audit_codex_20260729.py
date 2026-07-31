#!/usr/bin/env python3
"""Exact hostile audit of THM-2907's endpoint convention.

The child census uses a closed-comb interval engine, while the lonely-runner
danger predicate is strict.  This independent readout verifies that each
surviving family is deleted by the former only at equality and is safe for
the latter, and that the equality point is isolated.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIMARY = (
    ROOT
    / "04-computation/lrc14_j6_paircap_exception_h4_link_child_census_codex_20260729.py"
)
PRIMARY_SHA256 = (
    "88d523ea97235471ecce03c06de5cd1e1ba434ccd41fe0633beadf1017aa8fa3"
)
TARGET = F(1, 14)
TIME = F(1, 28)
EPSILON = F(1, 10**9)

EXPECTED_ROWS = (
    (
        (2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26),
        (),
        (2, 26),
        (2, 26),
        F(249999993, 3500000000),
        F(249999909, 3500000000),
    ),
    (
        (2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 26, 48),
        (),
        (2, 26),
        (2, 26),
        F(249999993, 3500000000),
        F(249999909, 3500000000),
    ),
)
EXPECTED_DIGEST = "b521e84294bc8b7ac9f431ec7cc6841bbc326bddf26ce60fcd546c55c9479c05"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_primary():
    require(file_sha256(PRIMARY) == PRIMARY_SHA256, "primary census changed")
    spec = importlib.util.spec_from_file_location("h4_endpoint_primary", PRIMARY)
    require(spec is not None and spec.loader is not None, "cannot load primary")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_primary()


def circle_distance(speed: int, time: F) -> F:
    residue = (speed * time) % 1
    return min(residue, 1 - residue)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def row_line(row: tuple[object, ...]) -> str:
    family, strict, boundary, closed, left, right = row
    return (
        f"family={family};strict={strict};boundary={boundary};closed={closed};"
        f"left={ftext(left)};right={ftext(right)}\n"
    )


def main() -> None:
    rows = []
    for family in M.EXPECTED_ENDPOINT_FAMILIES:
        distances = {
            speed: circle_distance(speed, TIME)
            for speed in family
        }
        strict = tuple(
            speed for speed, distance in distances.items() if distance < TARGET
        )
        boundary = tuple(
            speed for speed, distance in distances.items() if distance == TARGET
        )
        closed = tuple(
            speed for speed, distance in distances.items() if distance <= TARGET
        )
        good, components, mass = M.R.CORE.good_norm(family)
        require(
            good == [] and components == 0 and mass == 0,
            "closed-comb engine retained an endpoint carrier",
        )
        left = min(
            circle_distance(speed, TIME - EPSILON)
            for speed in family
        )
        right = min(
            circle_distance(speed, TIME + EPSILON)
            for speed in family
        )
        rows.append((family, strict, boundary, closed, left, right))

    rows = tuple(rows)
    require(rows == EXPECTED_ROWS, "endpoint hostile rows changed")
    require(
        all(
            not strict
            and boundary
            and closed == boundary
            and left < TARGET
            and right < TARGET
            for _, strict, boundary, closed, left, right in rows
        ),
        "strict/open endpoint convention failed",
    )

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/paircap-exception/endpoint-convention/v1\n")
    for row in rows:
        digest.update(row_line(row).encode())
    semantic_digest = digest.hexdigest()
    require(semantic_digest == EXPECTED_DIGEST, "semantic digest changed")

    print("LRC14 j6 pair-cap-exception endpoint-convention hostile audit")
    print(f"time={ftext(TIME)}")
    print(f"target={ftext(TARGET)}")
    print(f"epsilon={ftext(EPSILON)}")
    for row in rows:
        print(row_line(row).rstrip())
    print(f"semantic_digest={semantic_digest}")
    print("mode=LOCKED")
    print(
        "scope=two THM2907 survivor families;strict danger versus closed-comb "
        "endpoint convention;not LRC14"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
