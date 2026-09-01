"""Exact labelled-body edit and recurrence probe for THM-4296 boundaries.

The edit metric is replacement distance on nine-subsets:
    d(A,B) = |A triangle B|/2 = 9-|A intersect B|.
Carrier changes are printed as an explicit sidecar; no monotonicity in the
physical endpoint or in the carrier is assumed.
"""

from __future__ import annotations

import csv
import itertools
from collections import Counter
from pathlib import Path


ROOT = Path("05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296")


def read(path: Path, q: int | None = None) -> tuple[int, ...]:
    with path.open(newline="") as handle:
        rows = [
            (int(row["q"]), int(row["r"]), int(row["body_hex"], 16))
            for row in csv.DictReader(handle)
        ]
    bodies = tuple(body for row_q, _, body in rows if q is None or row_q == q)
    assert all(body.bit_count() == 9 for body in bodies)
    assert len(set(bodies)) == len(bodies)
    return bodies


def edit(left: int, right: int) -> int:
    assert left.bit_count() == right.bit_count() == 9
    distance = (left ^ right).bit_count()
    assert distance % 2 == 0
    return distance // 2


def nearest(source: tuple[int, ...], target: tuple[int, ...]):
    rows = []
    for body in source:
        distances = [(edit(body, candidate), candidate) for candidate in target]
        minimum = min(distance for distance, _ in distances)
        witnesses = tuple(candidate for distance, candidate in distances if distance == minimum)
        rows.append((body, minimum, witnesses))
    return rows


def injection(source: tuple[int, ...], target: tuple[int, ...]):
    """Map the small source injectively into target, minimizing (total,max,lex)."""
    assert len(source) <= 4 and len(source) <= len(target)
    best_key = None
    best_map = None
    optimum_count = 0
    for indices in itertools.permutations(range(len(target)), len(source)):
        distances = tuple(edit(source[i], target[j]) for i, j in enumerate(indices))
        mapping = tuple(target[j] for j in indices)
        key = (sum(distances), max(distances), mapping)
        if best_key is None or key < best_key:
            best_key = key
            best_map = tuple(zip(source, mapping, distances, strict=True))
            optimum_count = 1
        elif key[:2] == best_key[:2]:
            optimum_count += 1
    assert best_key is not None and best_map is not None
    return best_key[:2], optimum_count, best_map


def print_transition(name: str, before: tuple[int, ...], after: tuple[int, ...],
                     small_is_after: bool) -> None:
    forward = nearest(before, after)
    backward = nearest(after, before)
    print(
        "TRANSITION",
        name,
        "BEFORE",
        len(before),
        "AFTER",
        len(after),
        "INTERSECTION",
        len(set(before) & set(after)),
    )
    print(
        "DIRECTED BEFORE_TO_AFTER",
        "HIST",
        dict(sorted(Counter(distance for _, distance, _ in forward).items())),
        "HAUSDORFF",
        max(distance for _, distance, _ in forward),
    )
    print(
        "DIRECTED AFTER_TO_BEFORE",
        "HIST",
        dict(sorted(Counter(distance for _, distance, _ in backward).items())),
        "HAUSDORFF",
        max(distance for _, distance, _ in backward),
    )
    small, large = (after, before) if small_is_after else (before, after)
    (total, bottleneck), count, mapping = injection(small, large)
    print(
        "MIN_TOTAL_INJECTION",
        "SOURCE",
        "AFTER" if small_is_after else "BEFORE",
        "TOTAL",
        total,
        "BOTTLENECK",
        bottleneck,
        "OPTIMUM_COUNT",
        count,
    )
    for left, right, distance in mapping:
        print(f"MATCH {left:08x} -> {right:08x} EDIT {distance} XOR {left ^ right:08x}")


def main() -> None:
    boundary = {
        636: read(ROOT / "inputs/endpoint636_failures101.csv", q=100),
        632: read(ROOT / "inputs/post_exchange_r632_failures72.csv", q=100),
        630: read(ROOT / "results/endpoint/pre_r630_failure_boundary.csv", q=100),
        629: read(ROOT / "inputs/r629_failures28.csv", q=100),
        628: read(ROOT / "inputs/r628_failures4.csv", q=100),
        626: read(ROOT / "results/endpoint/final_mixed9019_failures.csv", q=100),
    }
    print("THM4296_ENDPOINT_LABEL_EDIT_PROBE_V1")
    print("METRIC replacement_distance=popcount(xor)/2")
    print("CARRIER_SIDECAR 630:C632 629:C630 628:C629 627:C628 626:C628")
    print_transition("630_TO_629", boundary[630], boundary[629], small_is_after=False)
    print_transition("629_TO_628", boundary[629], boundary[628], small_is_after=True)
    for earlier, later in [(636, 632), (632, 629), (629, 628), (628, 626)]:
        common = sorted(set(boundary[earlier]) & set(boundary[later]))
        print(
            "RECURRENCE",
            f"{earlier}_TO_{later}",
            "COUNT",
            len(common),
            "BODIES",
            ",".join(f"{body:08x}" for body in common) or "-",
        )
    assert set(boundary[628]) <= set(boundary[626])
    print("EXACT_SUBSET F628_IN_F626 4_OF_4 YES")
    print("VERDICT PASS EXACT_LABEL_EDIT_AND_RECURRENCE_SIDECAR")


if __name__ == "__main__":
    main()
