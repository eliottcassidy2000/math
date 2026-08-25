#!/usr/bin/env python3
"""Independent literal-theta wall/cell audit for THM-4087.

This companion never uses the parity-compressed tooth description from the
primary script.  It creates both literal walls
||d theta||=1/14 and ||d(theta+1/2)||=1/14, tests every wall and intervening
cell, and reconstructs open components from strict wall connectivity.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction


DELTA = Fraction(1, 14)
CORE = tuple(range(1, 10))
CORE_LEFT = Fraction(4, 49)
CORE_RIGHT = Fraction(3, 35)
H8 = (1, 3, 5, 8, 11, 13, 23, 36)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def literal_values(speed: int, theta: Fraction) -> tuple[Fraction, Fraction]:
    return (
        circle_norm(speed * theta),
        circle_norm(speed * (theta + Fraction(1, 2))),
    )


def literal_safe(speed: int, theta: Fraction) -> bool:
    first, second = literal_values(speed, theta)
    return first >= DELTA and second >= DELTA


def body_safe(body: tuple[int, ...], theta: Fraction) -> bool:
    return all(literal_safe(speed, theta) for speed in body)


def literal_walls(body: tuple[int, ...]) -> tuple[Fraction, ...]:
    walls = {Fraction(0), Fraction(1)}
    for speed in body:
        for half in (0, 1):
            for integer in range(speed):
                for sign in (-1, 1):
                    walls.add(
                        ((Fraction(integer) + sign * DELTA) / speed - Fraction(half, 2)) % 1
                    )
    return tuple(sorted(walls))


def literal_bad_component_lengths(first: int, second: int) -> tuple[Fraction, ...]:
    body = (first, second)
    walls = literal_walls(body)
    cell_count = len(walls) - 1
    bad_cells = [
        not body_safe(body, (walls[index] + walls[index + 1]) / 2)
        for index in range(cell_count)
    ]
    parents = list(range(cell_count))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parents[root_right] = root_left

    for index in range(cell_count - 1):
        if bad_cells[index] and bad_cells[index + 1] and not body_safe(body, walls[index + 1]):
            union(index, index + 1)
    if bad_cells[-1] and bad_cells[0] and not body_safe(body, Fraction(0)):
        union(cell_count - 1, 0)

    lengths: dict[int, Fraction] = {}
    for index, is_bad in enumerate(bad_cells):
        if is_bad:
            root = find(index)
            lengths[root] = lengths.get(root, Fraction(0)) + walls[index + 1] - walls[index]
    require(lengths, f"literal bad union unexpectedly empty for {first},{second}")
    return tuple(sorted(lengths.values()))


def candidate_walls(first: int, second: int) -> tuple[tuple[Fraction, int, str], ...]:
    candidates: dict[Fraction, tuple[int, str]] = {
        CORE_LEFT: (98, "core_left"),
        CORE_RIGHT: (70, "core_right"),
    }
    for speed, source in ((first, "first"), (second, "second")):
        for half in (0, 1):
            for integer in range(speed):
                for sign in (-1, 1):
                    wall = ((Fraction(integer) + sign * DELTA) / speed - Fraction(half, 2)) % 1
                    if CORE_LEFT <= wall <= CORE_RIGHT:
                        old = candidates.get(wall)
                        if old is None or 14 * speed < old[0]:
                            candidates[wall] = (14 * speed, source)
    return tuple((wall, value[0], value[1]) for wall, value in sorted(candidates.items()))


def arrangement_witness(first: int, second: int) -> tuple[Fraction, int, str]:
    for theta, clock, source in candidate_walls(first, second):
        if literal_safe(first, theta) and literal_safe(second, theta):
            return theta, clock, source
    raise RuntimeError(f"literal wall audit found no witness for {first},{second}")


def mod_norm(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def main() -> None:
    core_walls = literal_walls(CORE)
    require(CORE_LEFT in core_walls and CORE_RIGHT in core_walls, "core endpoints absent from walls")
    interval_walls = tuple(wall for wall in core_walls if CORE_LEFT <= wall <= CORE_RIGHT)
    core_wall_gates = 0
    core_cell_gates = 0
    for wall in interval_walls:
        for speed in CORE:
            require(literal_safe(speed, wall), f"unsafe core wall d={speed},theta={wall}")
            core_wall_gates += 2
    for left, right in zip(interval_walls, interval_walls[1:]):
        midpoint = (left + right) / 2
        for speed in CORE:
            require(literal_safe(speed, midpoint), f"unsafe core cell d={speed},theta={midpoint}")
            core_cell_gates += 2
    left_index = core_walls.index(CORE_LEFT)
    right_index = core_walls.index(CORE_RIGHT)
    require(
        not body_safe(CORE, (core_walls[left_index - 1] + CORE_LEFT) / 2),
        "core safe component extends left",
    )
    require(
        not body_safe(CORE, (CORE_RIGHT + core_walls[right_index + 1]) / 2),
        "core safe component extends right",
    )
    require(CORE_RIGHT - CORE_LEFT == Fraction(1, 245), "core width changed")

    component_pairs = 0
    component_count = 0
    worst_ratio = Fraction(0)
    worst_row = (0, 0, Fraction(0))
    for first in range(1, 51):
        for second in range(first + 1, 51):
            lengths = literal_bad_component_lengths(first, second)
            bound = Fraction(2, 7 * first)
            require(lengths[-1] < bound, f"literal component bound failed {first},{second}")
            ratio = lengths[-1] / bound
            if ratio > worst_ratio:
                worst_ratio = ratio
                worst_row = (first, second, lengths[-1])
            component_pairs += 1
            component_count += len(lengths)

    family_pairs = 0
    family_phase_gates = 0
    owner_class_gates = 0
    source_histogram = {"core_left": 0, "core_right": 0, "first": 0, "second": 0}
    sample_rows: list[tuple[int, int, str, str, int, int]] = []
    witness_bank: dict[tuple[int, int], tuple[Fraction, int, int]] = {}
    for first in range(70, 151):
        for second in range(first + 1, 151):
            theta, clock, source = arrangement_witness(first, second)
            require(clock % 2 == 0 and clock <= 14 * second, "literal adaptive clock bound failed")
            scaled = theta * clock
            require(scaled.denominator == 1, "literal witness missed clock grid")
            label = scaled.numerator % clock
            body = CORE + (first, second)
            for speed in body:
                values = literal_values(speed, theta)
                require(values[0] >= DELTA and values[1] >= DELTA, "literal family phase failed")
                family_phase_gates += 2
            opposite = (label + clock // 2) % clock
            for tail in range(1, 2 * clock, 2):
                require(
                    not (
                        7 * mod_norm(tail * label, clock) < clock
                        and 7 * mod_norm(tail * opposite, clock) < clock
                    ),
                    f"literal eligible class survived B={first},C={second},z={tail}",
                )
                owner_class_gates += 1
            source_histogram[source] += 1
            witness_bank[(first, second)] = (theta, clock, label)
            if (first, second) in ((70, 71), (70, 150), (98, 147), (149, 150)):
                sample_rows.append((first, second, str(theta), source, clock, label))
            family_pairs += 1

    dilation_phase_gates = 0
    dilation_owner_gates = 0
    for first, second in ((70, 71), (98, 147), (120, 149), (149, 150)):
        theta, clock, label = witness_bank[(first, second)]
        for dilation in range(1, 13):
            lifted = theta / dilation
            dilated_clock = dilation * clock
            for speed in CORE + (first, second):
                values = literal_values(dilation * speed, lifted)
                require(values[0] >= DELTA and values[1] >= DELTA, "literal dilation failed")
                dilation_phase_gates += 2
            opposite = (label + dilated_clock // 2) % dilated_clock
            for tail in (1, 3, 5, 7, 11):
                require(
                    not (
                        7 * mod_norm(tail * label, dilated_clock) < dilated_clock
                        and 7 * mod_norm(tail * opposite, dilated_clock) < dilated_clock
                    ),
                    "literal dilated owner hostile survived",
                )
                dilation_owner_gates += 1

    h8_walls = literal_walls(H8)
    h8_circle_walls = h8_walls[:-1]
    h8_safe_walls = tuple(theta for theta in h8_circle_walls if body_safe(H8, theta))
    expected_h8 = tuple(
        Fraction(index, 14) for index in tuple(range(1, 7)) + tuple(range(8, 14))
    )
    require(h8_safe_walls == expected_h8, "H8 literal safe wall set changed")
    h8_safe_cells = 0
    for left, right in zip(h8_walls, h8_walls[1:]):
        if body_safe(H8, (left + right) / 2):
            h8_safe_cells += 1
    require(h8_safe_cells == 0, "H8 acquired a safe open cell")
    h8_append_gates = 0
    for multiple in range(1, 101):
        appended = 7 * multiple
        for theta in expected_h8:
            require(not literal_safe(appended, theta), "literal 7m hostile failed")
            h8_append_gates += 1

    correction_seam = Fraction(1, 14)
    require(circle_norm(correction_seam) == DELTA, "MISTAKE-274 speed-one seam changed")
    require(circle_norm(13 * correction_seam) == DELTA, "MISTAKE-274 speed-thirteen seam changed")
    require(not (circle_norm(correction_seam) < DELTA), "speed-one open comb filled seam")
    require(not (circle_norm(13 * correction_seam) < DELTA), "speed-thirteen open comb filled seam")
    require(Fraction(13, 196) < Fraction(1, 14) < Fraction(15, 196), "{1,14} strict overlap changed")

    hostile_left = Fraction(139, 1722)
    hostile_right = Fraction(74, 861)
    hostile_walls = literal_walls((48, 123))
    require(hostile_left in hostile_walls and hostile_right in hostile_walls, "48/123 boundaries absent")
    require(body_safe((48, 123), hostile_left), "48/123 left boundary is not weak-safe")
    require(body_safe((48, 123), hostile_right), "48/123 right boundary is not weak-safe")
    require(hostile_left < CORE_LEFT < CORE_RIGHT < hostile_right, "48/123 does not contain core interval")
    hostile_segment = tuple(wall for wall in hostile_walls if hostile_left <= wall <= hostile_right)
    hostile_wall_gates = 0
    hostile_cell_gates = 0
    for wall in hostile_segment[1:-1]:
        require(not body_safe((48, 123), wall), "safe internal 48/123 seam found")
        hostile_wall_gates += 1
    for left, right in zip(hostile_segment, hostile_segment[1:]):
        require(not body_safe((48, 123), (left + right) / 2), "safe internal 48/123 cell found")
        hostile_cell_gates += 1
    require(hostile_right - hostile_left == Fraction(3, 574), "48/123 component length changed")

    semantic = {
        "component_worst": [worst_row[0], worst_row[1], str(worst_row[2]), str(worst_ratio)],
        "core_interval": [str(CORE_LEFT), str(CORE_RIGHT), len(interval_walls)],
        "family_samples": sample_rows,
        "hostile_48_123": [str(hostile_left), str(hostile_right)],
        "source_histogram": source_histogram,
    }
    digest = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4087 independent literal-theta wall/cell audit")
    print(f"core_walls_in_component={len(interval_walls)} core_wall_phase_gates={core_wall_gates} core_cell_phase_gates={core_cell_gates}")
    print(f"component_pairs={component_pairs} component_count={component_count} worst_pair={worst_row[0]},{worst_row[1]} worst_length={worst_row[2]} bound_ratio={worst_ratio}")
    print(f"family_pairs={family_pairs} family_phase_gates={family_phase_gates} owner_class_gates={owner_class_gates}")
    print("witness_sources=" + ",".join(f"{key}:{source_histogram[key]}" for key in sorted(source_histogram)))
    print(f"dilation_phase_gates={dilation_phase_gates} dilation_owner_gates={dilation_owner_gates}")
    print(f"H8_safe_walls={len(h8_safe_walls)} H8_safe_cells={h8_safe_cells} H8_7m_gates={h8_append_gates}")
    print(f"mistake274_open_handoff=1/14 local_48_123_component=[{hostile_left},{hostile_right}] hostile_wall_gates={hostile_wall_gates} hostile_cell_gates={hostile_cell_gates}")
    print(f"semantic_sha256={digest}")
    print("PASS: independent literal walls, clocks, dilation, endpoint hostiles, and correction lineage")


if __name__ == "__main__":
    main()
