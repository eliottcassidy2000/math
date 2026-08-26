#!/usr/bin/env python3
"""Independent exact referee for THM-4231.

Differences from primary.py:

* constructs the pool cells by endpoint toggles rather than midpoint tests;
* scatters mass and signed component-boundary coefficients upward to every
  depth-six deletion rather than summing submasks for each deletion;
* proves the no-cover consequence by scanning all C(30,9) bodies against a
  SplitMix-ordered edge deck rather than using a hitting-set DFS;
* uses the same endpoint-toggle method for the literal boundary controls.

This remains a scratch audit until reviewed and promoted.
"""

from __future__ import annotations

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
Q_PREVIOUS = 17_547
Q_TRANSITION = 17_548
MASK64 = (1 << 64) - 1


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def ceil_div(numerator: int, denominator: int) -> int:
    return (numerator + denominator - 1) // denominator


def labels(mask: int) -> tuple[int, ...]:
    return tuple(POOL[index] for index in range(30) if mask & (1 << index))


def mask_of(values: tuple[int, ...] | set[int]) -> int:
    index = {label: position for position, label in enumerate(POOL)}
    return sum(1 << index[value] for value in values)


class Fnv1a64:
    def __init__(self) -> None:
        self.value = 0xCBF29CE484222325

    def add_u64_le(self, value: int) -> None:
        require(0 <= value <= MASK64, "FNV field outside u64")
        for byte in value.to_bytes(8, "little"):
            self.value ^= byte
            self.value = (self.value * 0x100000001B3) & MASK64


def build_toggle_cells() -> tuple[int, list[tuple[int, int]]]:
    common = 1
    for speed in POOL:
        common = lcm(common, 14 * speed)

    # At (14k+1)/(14p) the p-danger tooth ends; at
    # (14k+13)/(14p) the next p-danger tooth starts.
    leave_at: dict[int, int] = {}
    enter_at: dict[int, int] = {}
    for vertex, speed in enumerate(POOL):
        bit = 1 << vertex
        unit = common // (14 * speed)
        for tooth in range(speed):
            leave = (14 * tooth + 1) * unit
            enter = (14 * tooth + 13) * unit
            leave_at[leave] = leave_at.get(leave, 0) | bit
            enter_at[enter] = enter_at.get(enter, 0) | bit

    event_ticks = sorted(set(leave_at) | set(enter_at))
    require(len(event_ticks) + 2 == 7_134, "toggle wall count changed")
    failed = (1 << 30) - 1
    previous = 0
    cells: list[tuple[int, int]] = []
    for tick in event_ticks:
        require(tick > previous, "non-increasing toggle event")
        cells.append((tick - previous, failed))
        failed &= ~leave_at.get(tick, 0)
        failed |= enter_at.get(tick, 0)
        previous = tick
    cells.append((common - previous, failed))
    require(failed == (1 << 30) - 1, "toggle sweep did not close")
    require(len(cells) == 7_133, "toggle cell count changed")
    return common, cells


def aggregate_coefficients(
    cells: list[tuple[int, int]],
) -> tuple[dict[int, int], dict[int, int]]:
    mass_coeff: dict[int, int] = {}
    component_coeff: dict[int, int] = {}

    for width, failed in cells:
        if failed.bit_count() <= DEPTH:
            mass_coeff[failed] = mass_coeff.get(failed, 0) + width
            component_coeff[failed] = component_coeff.get(failed, 0) + 1

    for index, (_, failed) in enumerate(cells):
        previous_failed = cells[(index - 1) % len(cells)][1]
        joined = previous_failed | failed
        if joined.bit_count() <= DEPTH:
            component_coeff[joined] = component_coeff.get(joined, 0) - 1

    component_coeff = {
        mask: coefficient
        for mask, coefficient in component_coeff.items()
        if coefficient != 0
    }
    return mass_coeff, component_coeff


def scatter_to_depth_six(
    masks: list[int],
    rank: dict[int, int],
    coefficients: dict[int, int],
) -> list[int]:
    result = [0] * len(masks)
    full_vertices = tuple(range(30))

    for base, value in coefficients.items():
        size = base.bit_count()
        require(size <= DEPTH, "oversized scatter base")
        available = tuple(
            vertex for vertex in full_vertices if not base & (1 << vertex)
        )
        for extension in combinations(available, DEPTH - size):
            target = base
            for vertex in extension:
                target |= 1 << vertex
            result[rank[target]] += value
    return result


def row_fingerprint(rows: list[tuple[int, int, int, int, int]]) -> int:
    digest = Fnv1a64()
    for deletion, activation, mass, components, surplus in sorted(rows):
        digest.add_u64_le(deletion)
        digest.add_u64_le(activation)
        digest.add_u64_le(mass)
        digest.add_u64_le(components)
        digest.add_u64_le(surplus)
    return digest.value


def splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return value ^ (value >> 31)


def exhaustive_body_scan(
    edges: list[int],
) -> tuple[int, int, int, int, int, int, int]:
    ordered = sorted(edges, key=lambda edge: (splitmix64(edge), edge))
    bodies = 0
    edge_checks = 0
    maximum_checks = 0
    closest_body = 0
    missed_edge = 0
    cover_count = 0
    first_cover = 0

    for choice in combinations(range(30), 9):
        body = sum(1 << vertex for vertex in choice)
        checked = 0
        missed = 0
        for edge in ordered:
            checked += 1
            if not edge & body:
                missed = edge
                break
        bodies += 1
        edge_checks += checked
        if checked > maximum_checks:
            maximum_checks = checked
            closest_body = body
            missed_edge = missed
        if missed == 0:
            cover_count += 1
            if first_cover == 0:
                first_cover = body
    return (
        bodies,
        edge_checks,
        maximum_checks,
        closest_body,
        missed_edge,
        cover_count,
        first_cover,
    )


def literal_toggle_geometry(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    common = 1
    for speed in speeds:
        common = lcm(common, 14 * speed)
    delta_at: dict[int, int] = {}
    for speed in speeds:
        unit = common // (14 * speed)
        for tooth in range(speed):
            leave = (14 * tooth + 1) * unit
            enter = (14 * tooth + 13) * unit
            delta_at[leave] = delta_at.get(leave, 0) - 1
            delta_at[enter] = delta_at.get(enter, 0) + 1

    failed = len(speeds)
    previous = 0
    previous_safe = False
    mass = 0
    components = 0
    for tick in sorted(delta_at):
        safe = failed == 0
        if safe:
            mass += tick - previous
            if not previous_safe:
                components += 1
        previous_safe = safe
        failed += delta_at[tick]
        require(failed >= 0, "negative literal failure count")
        previous = tick
    safe = failed == 0
    if safe:
        mass += common - previous
        if not previous_safe:
            components += 1
    require(failed == len(speeds), "literal toggle sweep did not close")
    return common, mass, components


def reduced(numerator: int, denominator: int) -> tuple[int, int]:
    divisor = gcd(abs(numerator), denominator)
    return numerator // divisor, denominator // divisor


def main() -> None:
    common, cells = build_toggle_cells()
    mass_coeff, component_coeff = aggregate_coefficients(cells)

    deletion_masks = [
        sum(1 << vertex for vertex in choice)
        for choice in combinations(range(30), DEPTH)
    ]
    require(len(deletion_masks) == comb(30, 6), "deletion mask universe changed")
    rank = {mask: index for index, mask in enumerate(deletion_masks)}
    require(len(rank) == len(deletion_masks), "deletion rank collision")

    masses = scatter_to_depth_six(deletion_masks, rank, mass_coeff)
    components = scatter_to_depth_six(
        deletion_masks, rank, component_coeff
    )

    rows: list[tuple[int, int, int, int, int]] = []
    nonpositive = 0
    equalities = 0
    ceiling_checks = 0
    for deletion, mass, component_count in zip(
        deletion_masks, masses, components
    ):
        require(mass > 0 and component_count > 0, "invalid scattered geometry")
        surplus = 45 * mass - 4 * common
        if surplus == 0:
            equalities += 1
        if surplus <= 0:
            nonpositive += 1
            continue
        numerator = 108 * component_count * common
        denominator = 7 * surplus
        activation = ceil_div(numerator, denominator)
        require(
            denominator * activation >= numerator,
            "ceiling upper inequality failed",
        )
        require(
            denominator * (activation - 1) < numerator,
            "ceiling lower inequality failed",
        )
        ceiling_checks += 1
        rows.append((deletion, activation, mass, component_count, surplus))

    require(len(rows) == 140_082, "strict reverse-scatter census changed")
    require(nonpositive == 453_693 and equalities == 0, "reverse-scatter signs changed")
    previous = [row for row in rows if row[1] <= Q_PREVIOUS]
    transition = [row for row in rows if row[1] <= Q_TRANSITION]
    require(len(previous) == 54_563 and len(transition) == 54_566, "filter count changed")
    require(row_fingerprint(rows) == 0xA8B79AD77AD91A62, "full ledger mismatch")
    require(row_fingerprint(previous) == 0x476FEF92619D2C0B, "previous ledger mismatch")
    require(row_fingerprint(transition) == 0xD20636ACE1522A29, "transition ledger mismatch")

    previous_by_inequality = [
        row
        for row in rows
        if 7 * row[4] * Q_PREVIOUS >= 108 * row[3] * common
    ]
    transition_by_inequality = [
        row
        for row in rows
        if 7 * row[4] * Q_TRANSITION >= 108 * row[3] * common
    ]
    require(previous_by_inequality == previous, "Q17547 rational filter mismatch")
    require(transition_by_inequality == transition, "Q17548 rational filter mismatch")
    previous_boundary_equalities = sum(
        7 * row[4] * Q_PREVIOUS == 108 * row[3] * common for row in rows
    )
    transition_boundary_equalities = sum(
        7 * row[4] * Q_TRANSITION == 108 * row[3] * common for row in rows
    )

    previous_scan = exhaustive_body_scan([row[0] for row in previous])
    (
        previous_bodies,
        previous_edge_checks,
        previous_maximum_checks,
        previous_closest,
        previous_missed,
        previous_cover_count,
        previous_first_cover,
    ) = previous_scan
    transition_scan = exhaustive_body_scan([row[0] for row in transition])
    (
        bodies,
        edge_checks,
        maximum_checks,
        closest,
        missed,
        cover_count,
        first_cover,
    ) = transition_scan
    require(previous_bodies == comb(30, 9), "previous body universe changed")
    require(bodies == comb(30, 9), "transition body universe changed")
    require(previous_edge_checks == 233_058_925, "previous body incidence total changed")
    require(previous_maximum_checks == 54_563, "maximum previous body trace changed")
    require(previous_cover_count == 1, "previous cover census changed")
    require(
        previous_closest
        == previous_first_cover
        == mask_of((85, 88, 143, 168, 193, 240, 252, 264, 290))
        and previous_missed == 0,
        "unique previous cover changed",
    )
    require(cover_count == 0 and first_cover == 0, "referee found a transition cover")
    require(edge_checks == 233_056_301, "transition body incidence total changed")
    require(maximum_checks == 51_934, "maximum transition body trace changed")
    require(
        closest == mask_of((85, 88, 143, 168, 193, 240, 252, 264, 290))
        and missed == mask_of((8, 16, 42, 95, 132, 145)),
        "closest body or missed edge changed",
    )

    zero_outsider = comb(30, 11)
    one_outsider_each = comb(30, 10)
    two_outsider = comb(30, 9)
    lifted_total = zero_outsider + 2 * one_outsider_each + two_outsider
    require(lifted_total == comb(32, 11), "Pascal outsider partition failed")

    repair_labels = tuple(
        sorted((set(POOL) - set(labels(missed))) | {17_547, 17_548})
    )
    body_labels = tuple(sorted(set(labels(closest)) | {17_547, 17_548}))
    repair_common, repair_mass, repair_components = literal_toggle_geometry(repair_labels)
    body_common, body_mass, body_components = literal_toggle_geometry(body_labels)
    repair_delta = 63 * repair_mass - 4 * repair_common
    body_delta = 63 * body_mass - 4 * body_common
    require(repair_delta > 0 and body_delta > 0, "literal toggle hostile changed")
    repair_fraction = reduced(repair_mass, repair_common)
    repair_delta_fraction = reduced(repair_delta, repair_common)
    body_fraction = reduced(body_mass, body_common)
    body_delta_fraction = reduced(body_delta, body_common)

    print("LRC14_ARBITRARY_PAIR_COFINAL_D6_INDEPENDENT_REFEREE_THM4231")
    print(
        f"TOGGLE_GEOMETRY COMMON {common} WALLS {len(cells)+1} CELLS {len(cells)} "
        f"MASS_COEFF {len(mass_coeff)} COMPONENT_COEFF {len(component_coeff)}"
    )
    print(
        f"REVERSE_SCATTER D6_UNIVERSE {len(deletion_masks)} STRICT {len(rows)} "
        f"NONPOSITIVE {nonpositive} EQUALITIES {equalities}"
    )
    print(
        f"FILTER_COUNTS Q17547 {len(previous)} Q17548 {len(transition)} "
        "FINGERPRINTS a8b79ad77ad91a62,476fef92619d2c0b,d20636ace1522a29"
    )
    print(
        f"RATIONAL_CEILING_AUDIT ROWS {ceiling_checks} "
        f"Q17547_MATCH {len(previous_by_inequality)} EQUALITIES {previous_boundary_equalities} "
        f"Q17548_MATCH {len(transition_by_inequality)} EQUALITIES {transition_boundary_equalities}"
    )
    print(
        f"EXHAUSTIVE_Q17547 BODIES {previous_bodies} "
        f"EDGE_CHECKS {previous_edge_checks} MAX_CHECKS {previous_maximum_checks} "
        f"CLOSEST {labels(previous_closest)} MISSED "
        f"{labels(previous_missed) if previous_missed else 'NONE'} "
        f"COVERS {previous_cover_count} FIRST {labels(previous_first_cover)}"
    )
    print(
        f"EXHAUSTIVE_Q17548 BODIES {bodies} EDGE_CHECKS {edge_checks} "
        f"MAX_CHECKS {maximum_checks} CLOSEST {labels(closest)} "
        f"MISSED {labels(missed)} COVERS {cover_count}"
    )
    print(
        f"PASCAL_OUTSIDER_BLOCKS ZERO {zero_outsider} "
        f"ONE_EACH {one_outsider_each} ONE_LABELS 2 "
        f"TWO {two_outsider} TOTAL {lifted_total} C32_11 {comb(32, 11)}"
    )
    print(
        f"LITERAL_TOGGLE_REPAIR COMPONENTS {repair_components} "
        f"MASS {repair_fraction[0]}/{repair_fraction[1]} "
        f"DELTA63 {repair_delta_fraction[0]}/{repair_delta_fraction[1]} POSITIVE YES"
    )
    print(
        f"LITERAL_TOGGLE_BODY COMPONENTS {body_components} "
        f"MASS {body_fraction[0]}/{body_fraction[1]} "
        f"DELTA63 {body_delta_fraction[0]}/{body_delta_fraction[1]} POSITIVE YES"
    )
    print(
        "VERDICT ACCEPT COMPLETE_REVERSE_SCATTER_AND_BODY_SCAN_AGREE "
        "Q17548_ACTIVATION_FILTER_HAS_TAU_GT_9 LRC14_REMAINS_OPEN"
    )


if __name__ == "__main__":
    main()
