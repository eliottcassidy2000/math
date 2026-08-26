#!/usr/bin/env python3
"""Exact primary scout for the arbitrary-pair cofinal LRC14 repair deck.

This is a scratch artifact, not canon.  It exhausts all depth-six deletions
of the THM-4156 pool, computes their exact safe-set mass and component count,
forms the two-newcomer discrepancy/Bonferroni activation deck, and proves the
adjacent activation-filter transition 17547 -> 17548 by a path-complete
hitting-set search.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations
from math import comb, gcd


POOL = (
    8,
    10,
    15,
    16,
    20,
    30,
    40,
    42,
    60,
    63,
    80,
    84,
    85,
    88,
    95,
    120,
    126,
    132,
    143,
    145,
    168,
    170,
    176,
    190,
    193,
    240,
    252,
    264,
    286,
    290,
)
DEPTH = 6
BODY_SIZE = 9
Q_PREVIOUS = 17_547
Q_TRANSITION = 17_548
MASK64 = (1 << 64) - 1


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def ceil_div(numerator: int, denominator: int) -> int:
    require(numerator >= 0 and denominator > 0, "invalid ceiling quotient")
    return (numerator + denominator - 1) // denominator


def mask_of(labels: tuple[int, ...] | list[int] | set[int]) -> int:
    index = {label: position for position, label in enumerate(POOL)}
    result = 0
    for label in labels:
        result |= 1 << index[label]
    return result


def labels(mask: int) -> tuple[int, ...]:
    return tuple(POOL[index] for index in range(30) if mask & (1 << index))


class Fnv1a64:
    def __init__(self) -> None:
        self.value = 0xCBF29CE484222325

    def add_u64_le(self, value: int) -> None:
        require(0 <= value <= MASK64, "FNV field outside u64")
        for byte in value.to_bytes(8, "little"):
            self.value ^= byte
            self.value = (self.value * 0x100000001B3) & MASK64


def build_pool_cells() -> tuple[int, list[tuple[int, int]]]:
    common = 1
    for speed in POOL:
        common = lcm(common, 14 * speed)

    walls = {0, common}
    for speed in POOL:
        unit = common // (14 * speed)
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    ordered = sorted(walls)
    require(len(ordered) == 7_134, "pool wall count changed")

    cells: list[tuple[int, int]] = []
    period = 2 * common
    for left, right in zip(ordered, ordered[1:]):
        failed = 0
        midpoint_twice = left + right
        for vertex, speed in enumerate(POOL):
            residue = speed * midpoint_twice % period
            if not (7 * residue >= common and 7 * residue <= 13 * common):
                failed |= 1 << vertex
        cells.append((right - left, failed))
    require(len(cells) == 7_133, "pool cell count changed")
    return common, cells


def build_fixed_stats(
    cells: list[tuple[int, int]],
) -> tuple[dict[int, int], dict[int, int], dict[int, int]]:
    length: dict[int, int] = {}
    cell_count: dict[int, int] = {}
    adjacency_union_count: dict[int, int] = {}

    for width, failed in cells:
        if failed.bit_count() <= DEPTH:
            length[failed] = length.get(failed, 0) + width
            cell_count[failed] = cell_count.get(failed, 0) + 1

    for index, (_, failed) in enumerate(cells):
        joined = failed | cells[(index + 1) % len(cells)][1]
        if joined.bit_count() <= DEPTH:
            adjacency_union_count[joined] = (
                adjacency_union_count.get(joined, 0) + 1
            )
    return length, cell_count, adjacency_union_count


def subset_sum(mask: int, table: dict[int, int]) -> int:
    result = 0
    subset = mask
    while True:
        result += table.get(subset, 0)
        if subset == 0:
            return result
        subset = (subset - 1) & mask


def deletion_geometry(
    deletion: int,
    length: dict[int, int],
    cell_count: dict[int, int],
    adjacency_union_count: dict[int, int],
) -> tuple[int, int]:
    mass = subset_sum(deletion, length)
    safe_cells = subset_sum(deletion, cell_count)
    safe_adjacencies = subset_sum(deletion, adjacency_union_count)
    components = safe_cells - safe_adjacencies
    require(mass > 0 and components > 0, "invalid deletion geometry")
    return mass, components


def enumerate_rows(
    common: int,
    length: dict[int, int],
    cell_count: dict[int, int],
    adjacency_union_count: dict[int, int],
) -> tuple[list[tuple[int, int, int, int, int]], int, int]:
    rows: list[tuple[int, int, int, int, int]] = []
    nonpositive = 0
    equalities = 0

    for choice in combinations(range(30), DEPTH):
        deletion = sum(1 << vertex for vertex in choice)
        mass, components = deletion_geometry(
            deletion, length, cell_count, adjacency_union_count
        )
        surplus_numerator = 45 * mass - 4 * common
        if surplus_numerator == 0:
            equalities += 1
        if surplus_numerator <= 0:
            nonpositive += 1
            continue

        activation = ceil_div(
            108 * components * common, 7 * surplus_numerator
        )
        rows.append(
            (deletion, activation, mass, components, surplus_numerator)
        )

    return rows, nonpositive, equalities


def row_fingerprint(rows: list[tuple[int, int, int, int, int]]) -> int:
    digest = Fnv1a64()
    for deletion, activation, mass, components, surplus in sorted(rows):
        digest.add_u64_le(deletion)
        digest.add_u64_le(activation)
        digest.add_u64_le(mass)
        digest.add_u64_le(components)
        digest.add_u64_le(surplus)
    return digest.value


def active_rows(
    rows: list[tuple[int, int, int, int, int]], cutoff: int
) -> list[tuple[int, int, int, int, int]]:
    return [row for row in rows if row[1] <= cutoff]


def path_complete_cover_search(edges: list[int]) -> tuple[int, int, int]:
    """Return (cover-or-minus-one, visited states, dead states).

    At every nonterminal state the recursion chooses the first uncovered edge
    and branches on every one of its vertices.  Any cover extending the chosen
    mask must contain one of those vertices, so the recursion is path-complete.
    """

    frequencies = [sum((edge >> vertex) & 1 for edge in edges) for vertex in range(30)]
    ordered_edges = sorted(
        edges,
        key=lambda edge: (
            sum(
                frequencies[vertex]
                for vertex in range(30)
                if edge & (1 << vertex)
            ),
            edge,
        ),
    )
    nodes = 0
    dead = 0

    @lru_cache(maxsize=None)
    def search(chosen: int, remaining: int) -> int:
        nonlocal nodes, dead
        nodes += 1
        uncovered = next((edge for edge in ordered_edges if not edge & chosen), 0)
        if uncovered == 0:
            return chosen
        if remaining == 0:
            dead += 1
            return -1

        vertices = [
            vertex for vertex in range(30) if uncovered & (1 << vertex)
        ]
        vertices.sort(key=lambda vertex: (-frequencies[vertex], vertex))
        for vertex in vertices:
            answer = search(chosen | (1 << vertex), remaining - 1)
            if answer != -1:
                return answer
        dead += 1
        return -1

    answer = search(0, BODY_SIZE)
    return answer, nodes, dead


def direct_safe_geometry(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    """Independent literal wall scan for one displayed boundary row."""

    common = 1
    for speed in speeds:
        common = lcm(common, 14 * speed)
    walls = {0, common}
    for speed in speeds:
        unit = common // (14 * speed)
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    ordered = sorted(walls)
    mass = 0
    components = 0
    previous_safe = False
    period = 2 * common
    for left, right in zip(ordered, ordered[1:]):
        midpoint_twice = left + right
        safe = True
        for speed in speeds:
            residue = speed * midpoint_twice % period
            if not (7 * residue >= common and 7 * residue <= 13 * common):
                safe = False
                break
        if safe:
            mass += right - left
            if not previous_safe:
                components += 1
        previous_safe = safe
    return common, mass, components


def reduced(numerator: int, denominator: int) -> tuple[int, int]:
    divisor = gcd(abs(numerator), denominator)
    return numerator // divisor, denominator // divisor


def main() -> None:
    common, cells = build_pool_cells()
    length, cell_count, adjacency_union_count = build_fixed_stats(cells)

    empty_mass, empty_components = deletion_geometry(
        0, length, cell_count, adjacency_union_count
    )
    require(
        empty_mass == 1_192_533_424_636 and empty_components == 150,
        "THM-4156 full-pool control changed",
    )
    triple_control = mask_of((8, 84, 252))
    triple_mass, triple_components = deletion_geometry(
        triple_control, length, cell_count, adjacency_union_count
    )
    require(
        triple_mass == 1_394_198_307_034 and triple_components == 148,
        "THM-4170 component control changed",
    )

    rows, nonpositive, equalities = enumerate_rows(
        common, length, cell_count, adjacency_union_count
    )
    require(len(rows) + nonpositive == comb(30, 6), "depth-six universe changed")
    require(len(rows) == 140_082 and nonpositive == 453_693, "strict census changed")
    require(equalities == 0, "strict-limit equality appeared")

    minimum = min(rows, key=lambda row: (row[1], row[0]))
    maximum = max(rows, key=lambda row: (row[1], -row[0]))
    previous = active_rows(rows, Q_PREVIOUS)
    transition = active_rows(rows, Q_TRANSITION)
    previous_masks = [row[0] for row in previous]
    transition_masks = [row[0] for row in transition]

    witness = mask_of((85, 88, 143, 168, 193, 240, 252, 264, 290))
    require(witness.bit_count() == 9, "previous witness size changed")
    require(
        all(edge & witness for edge in previous_masks),
        "previous witness stopped covering the filtered deck",
    )
    new_breakers = [row for row in transition if not row[0] & witness]
    require(new_breakers, "transition did not break the previous cover")
    breaker = min(new_breakers, key=lambda row: row[0])
    require(
        breaker[0] == mask_of((8, 16, 42, 95, 132, 145))
        and breaker[1] == Q_TRANSITION,
        "first displayed breaker changed",
    )

    answer, nodes, dead = path_complete_cover_search(transition_masks)
    require(answer == -1, "transition deck acquired a nine-cover")

    # Literal controls at the adjacent filtered boundary.  These show that a
    # sufficient-deck cover or activation failure is not literal danger.
    repair_labels = tuple(
        sorted((set(POOL) - set(labels(breaker[0]))) | {17_547, 17_548})
    )
    body_labels = tuple(sorted(set(labels(witness)) | {17_547, 17_548}))
    repair_common, repair_mass, repair_components = direct_safe_geometry(repair_labels)
    body_common, body_mass, body_components = direct_safe_geometry(body_labels)
    repair_delta = 63 * repair_mass - 4 * repair_common
    body_delta = 63 * body_mass - 4 * body_common
    require(repair_delta > 0 and body_delta > 0, "literal boundary hostile lost positivity")

    full_fingerprint = row_fingerprint(rows)
    previous_fingerprint = row_fingerprint(previous)
    transition_fingerprint = row_fingerprint(transition)
    repair_fraction = reduced(repair_mass, repair_common)
    repair_delta_fraction = reduced(repair_delta, repair_common)
    body_fraction = reduced(body_mass, body_common)
    body_delta_fraction = reduced(body_delta, body_common)

    print("LRC14_ARBITRARY_PAIR_COFINAL_D6_EXACT_SCOUT_THM4231")
    print(
        f"POOL COMMON {common} WALLS {len(cells)+1} CELLS {len(cells)} "
        f"FULL_MASS_TICKS {empty_mass} FULL_COMPONENTS {empty_components}"
    )
    print(
        f"THM4170_CONTROL DELETION {labels(triple_control)} MASS_TICKS {triple_mass} "
        f"COMPONENTS {triple_components}"
    )
    print(
        f"D6_UNIVERSE {comb(30,6)} STRICT {len(rows)} NONPOSITIVE {nonpositive} "
        f"EQUALITIES {equalities}"
    )
    print(
        f"ACTIVATION_MIN {minimum[1]} EDGE {labels(minimum[0])} "
        f"MASS_TICKS {minimum[2]} COMPONENTS {minimum[3]} SURPLUS45 {minimum[4]}"
    )
    print(
        f"ACTIVATION_MAX {maximum[1]} EDGE {labels(maximum[0])} "
        f"MASS_TICKS {maximum[2]} COMPONENTS {maximum[3]} SURPLUS45 {maximum[4]}"
    )
    print(
        f"FILTER_Q{Q_PREVIOUS} EDGES {len(previous)} COVER9 {labels(witness)} "
        f"FNV1A64_LE {previous_fingerprint:016x}"
    )
    print(
        f"FILTER_Q{Q_TRANSITION} EDGES {len(transition)} NEW_EDGES "
        f"{len(transition)-len(previous)} BREAKER {labels(breaker[0])} "
        f"BREAKER_MASS_TICKS {breaker[2]} BREAKER_COMPONENTS {breaker[3]} "
        f"BREAKER_SURPLUS45 {breaker[4]} FNV1A64_LE {transition_fingerprint:016x}"
    )
    print(
        f"PATH_COMPLETE_NO_COVER BUDGET {BODY_SIZE} NODES {nodes} DEAD {dead} "
        "COVER NONE"
    )
    print(
        "LITERAL_PREVIOUS_REPAIR "
        f"PAIR (17547,17548) EDGE {labels(breaker[0])} COMPONENTS {repair_components} "
        f"MASS {repair_fraction[0]}/{repair_fraction[1]} "
        f"DELTA63 {repair_delta_fraction[0]}/{repair_delta_fraction[1]} POSITIVE YES"
    )
    print(
        "LITERAL_PREVIOUS_BODY "
        f"PAIR (17547,17548) BODY {labels(witness)} COMPONENTS {body_components} "
        f"MASS {body_fraction[0]}/{body_fraction[1]} "
        f"DELTA63 {body_delta_fraction[0]}/{body_delta_fraction[1]} POSITIVE YES"
    )
    print(
        f"FULL_STRICT_LEDGER_FNV1A64_LE {full_fingerprint:016x} "
        "PAIR_BODY_COUNT 14307150 COMPLETE_32_CHART_COUNT 129024480"
    )
    print(
        "VERDICT PASS Q17548_IS_MINIMAL_FOR_THE_COMPLETE_D6_ACTIVATION_FILTER "
        "EVERY_DISTINCT_Q_R_GE_17548_HAS_ALL_C30_9_TWO_NEWCOMER_BODIES_SAFE "
        "LITERAL_MINIMAL_ENTRY_AND_LRC14_REMAIN_OPEN"
    )


if __name__ == "__main__":
    main()
