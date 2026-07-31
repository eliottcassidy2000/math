#!/usr/bin/env python3
"""Exact semantic-cospan audit of the THM-2707 common packet fibre.

THM-2707 proves that 3,346 fixed-skeleton packet addresses share one open
base cylinder.  This companion tests every point of every such cylinder
against the literal THM-2305 exclusive-source and prescribed-clock terminal
word predicates.  It is a semantic audit of that proved support SCC, not a
new packet census and not an LRC(14) row exclusion.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_full_physical_lift_fibre_thm2707 as fibre


old = fibre.old
m = fibre.m
P = 13
R = P**6
GUARD = 1
UNITS = (14, 27, 40, 53, 66)
BLOCKERS = (P, P**3, 2 * P**5)
CLOCKS = (2, 4, 6)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def centered(value):
    residue = frac(value)
    return residue if residue < Fraction(1, 2) else residue - 1


def danger(speed, value, denominator=14):
    phase = centered(speed * value)
    radius = Fraction(1, denominator)
    return -radius <= phase < radius


def common_safe(value):
    return (
        not danger(GUARD, value, 7)
        and all(not danger(speed, value) for speed in UNITS)
    )


def exclusive_source(value):
    if not common_safe(value):
        return 0
    hits = tuple(danger(speed, value) for speed in BLOCKERS)
    if sum(hits) != 1:
        return 0
    return hits.index(True) + 1


def terminal_word(value, source):
    require(source in (1, 2, 3), "terminal word needs a source label")
    source_index = source - 1
    if not common_safe(value) or danger(BLOCKERS[source_index], value):
        return None
    word = tuple(
        index + 1
        for index, speed in enumerate(BLOCKERS)
        if index != source_index and danger(speed, value)
    )
    return word or None


def semantic_record(value):
    source = exclusive_source(value)
    if not source:
        return 0, None
    endpoint = frac(P ** CLOCKS[source - 1] * value)
    return source, terminal_word(endpoint, source)


def floor_fraction(value):
    return value.numerator // value.denominator


def semantic_partition(value, radius):
    """Partition one open packet interval at every semantic predicate wall."""
    left = value - radius
    right = value + radius
    walls = {left, right}
    banks = [(1, GUARD, 7)]
    banks.extend((1, speed, 14) for speed in UNITS + BLOCKERS)
    for clock in CLOCKS:
        expansion = P**clock
        banks.append((expansion, GUARD, 7))
        banks.extend(
            (expansion, speed, 14) for speed in UNITS + BLOCKERS
        )
    for expansion, speed, denominator in banks:
        scale = expansion * speed
        danger_radius = Fraction(1, denominator)
        scaled_left = scale * left
        scaled_right = scale * right
        for sign in (-1, 1):
            shift = sign * danger_radius
            start = floor_fraction(scaled_left - shift) - 1
            stop = floor_fraction(scaled_right - shift) + 2
            for integer in range(start, stop + 1):
                wall = Fraction(integer, scale) + shift / scale
                if left < wall < right:
                    walls.add(wall)
    ordered = tuple(sorted(walls))
    cells = []
    for cell_left, cell_right in zip(ordered, ordered[1:]):
        if cell_left == cell_right:
            continue
        midpoint = (cell_left + cell_right) / 2
        source, word = semantic_record(midpoint)
        cells.append((cell_left, cell_right, source, word))
    require(sum(cell_right - cell_left
                for cell_left, cell_right, _source, _word in cells)
            == 2 * radius,
            "semantic cell partition lost interval length")
    return tuple(cells)


def main():
    x = Fraction(649039434905733, 1304692766858936)
    z = Fraction(46873542509301, 100360982066072)
    base_radius = Fraction(1, 1304692766858936)
    physical_radius = P * base_radius
    interval = (x - base_radius, x + base_radius)
    require(R == m.R and frac(P * x) == z,
            "THM-2707 base constants changed")

    module, _prefixes, _, _, rails, present, _starts = (
        m.core.build_carrier_data()
    )
    rows = m.shard((0, 1))[6][0]
    rail_support = old.merge_support(rails[0][3])
    present_support = tuple(present[1, (-6) % P])

    good_residues = tuple(
        residue for residue in range(P)
        if (
            (root := (6 + residue) % P)
            and m.is_unit(
                rows[0][0][(2 + 7 * residue) % P][1][6], root, 26
            )
        )
    )
    require(good_residues == (0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12),
            "THM-2707 unit residue bank changed")

    denominator = (z * m.T).denominator
    point = (z * m.T).numerator
    modulus = m.T * denominator
    step = (7 * m.T // R) * denominator
    scaled_rail = tuple(
        (left * denominator, right * denominator)
        for left, right in rail_support
    )
    scaled_present = tuple(
        (left * denominator, right * denominator)
        for left, right in present_support
    )
    rail_starts = tuple(left for left, _ in scaled_rail)
    present_starts = tuple(left for left, _ in scaled_present)

    good = []
    for address in range(R):
        if (
            address % P in good_residues
            and fibre.strict_interval_index(
                point, rail_starts, scaled_rail
            )
            and fibre.strict_interval_index(
                point, present_starts, scaled_present
            )
        ):
            good.append(address)
        point = (point + step) % modulus
    require(len(good) == 3346 and 110 in good,
            "THM-2707 full packet census changed")

    source_counts = Counter()
    semantic_counts = Counter()
    source_residue_counts = Counter()
    e2_rows = []
    cut_addresses = 0
    maximum_cells = 0
    cell_semantic_counts = Counter()
    cell_semantic_lengths = Counter()
    positive_cospan_cells = []
    zero_blocker_nonexclusive = 0
    near_q2_pure_one = []
    for address in good:
        value, _carry, _root, unit = fibre.packet_address_data(
            address, z, rows
        )
        require(unit and 0 < value - physical_radius
                < value + physical_radius < 1,
                "a packet cylinder wrapped or lost its unit")
        source, word = semantic_record(value)
        source_counts[source] += 1
        semantic_counts[(source, word)] += 1
        source_residue_counts[(source, address % P)] += 1
        blocker_hits = tuple(danger(speed, value) for speed in BLOCKERS)
        if source == 0:
            require(common_safe(value) and sum(blocker_hits) == 0,
                    "a nonexclusive address is not the zero-blocker A0 case")
            zero_blocker_nonexclusive += 1
        if source == 2:
            endpoint = frac(P**4 * value)
            e2_rows.append((address, value, endpoint))
            if (
                common_safe(endpoint)
                and danger(BLOCKERS[0], endpoint)
                and not danger(BLOCKERS[2], endpoint)
            ):
                near_q2_pure_one.append(address)

        cells = semantic_partition(value, physical_radius)
        cut_addresses += len(cells) > 1
        maximum_cells = max(maximum_cells, len(cells))
        for cell_left, cell_right, cell_source, cell_word in cells:
            key = (cell_source, cell_word)
            cell_semantic_counts[key] += 1
            cell_semantic_lengths[key] += cell_right - cell_left
            if cell_source and cell_word:
                positive_cospan_cells.append(
                    (address, cell_left, cell_right,
                     cell_source, cell_word)
                )

    require(source_counts == Counter({0: 3081, 2: 265}),
            "full-fibre exclusive-source census changed")
    require(semantic_counts == Counter({(0, None): 3081, (2, None): 265}),
            "midpoint semantic census changed")
    require(len(e2_rows) == 265, "E2 source bank changed")
    require(zero_blocker_nonexclusive == 3081,
            "zero-blocker nonexclusive control changed")
    require(len(near_q2_pure_one) == 16
            and all(address % P != 110 % P
                    for address in near_q2_pure_one),
            "sharp one-bit Q_(2,{1}) near-cospan control changed")

    # For q_n=z+7n/R, the source-c2 phase after its prescribed four-step
    # clock is independent of n:
    #
    #   c2 D^4(q_n) = 13^7 z + 91 n  (mod 1).
    #
    # It remains strictly dangerous on the whole open q_n(I) cylinder.
    invariant_phase = frac(P**7 * z)
    invariant_distance = abs(centered(invariant_phase))
    endpoint_phase_radius = P**7 * physical_radius
    invariant_slack = (
        Fraction(1, 14) - invariant_distance - endpoint_phase_radius
    )
    require(invariant_distance == Fraction(675, 1599416),
            "fixed endpoint c2 phase changed")
    require(invariant_slack > 0,
            "fixed endpoint c2 danger stopped being whole-cylinder strict")
    for address, value, endpoint in e2_rows:
        require(
            frac(BLOCKERS[1] * endpoint) == invariant_phase
            and frac(P**7 * value) == invariant_phase
            and danger(BLOCKERS[1], endpoint),
            f"E2 endpoint invariant failed at address {address}",
        )
    require(not positive_cospan_cells,
            "a positive THM-2305 source-to-word subinterval appeared")

    print("LRC14 THM-2707 FULL-FIBRE SEMANTIC COSPAN AUDIT")
    print("status=FINITE-EXACT WHOLE-OPEN-CYLINDER NO-GO")
    print(f"packet_addresses={len(good)} common_base_I={interval}")
    print(f"physical_q_radius={physical_radius}")
    print(f"source_counts={tuple(sorted(source_counts.items()))}")
    print(
        "source_residue_counts="
        f"{tuple(sorted(source_residue_counts.items()))}"
    )
    print(
        "semantic_counts="
        f"{tuple(sorted(semantic_counts.items(), key=repr))}"
    )
    print("nonexclusive_A0_zero_blockers=3081")
    print("near_Q2_{1}_except_c2_safe=16 all_nonzero_lifts=True")
    print(f"semantic_cut_addresses={cut_addresses} "
          f"maximum_open_cells={maximum_cells}")
    print("open_cell_semantic_counts="
          f"{tuple(sorted(cell_semantic_counts.items(), key=repr))}")
    print("open_cell_semantic_lengths="
          f"{tuple(sorted(cell_semantic_lengths.items(), key=repr))}")
    print(
        "E2_prescribed_endpoint_identity="
        "c2*D^4(q_n)=13^7*z+91*n mod1"
    )
    print(
        f"fixed_c2_endpoint_phase={invariant_phase} "
        f"distance={invariant_distance} "
        f"whole_cylinder_danger_slack={invariant_slack}"
    )
    print(
        "verdict=no packet cylinder contains a positive open literal "
        "THM-2305 exclusive-source-to-terminal-word cospan"
    )
    print(
        "first_failed_predicate=all 265 exclusive E2 sources retain "
        "source blocker c2 at D^4; remaining 3081 addresses are not "
        "exclusive sources"
    )
    print(
        "scope=THM-2707 fixed packet skeleton only; no row exclusion; "
        "LRC14 remains open"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
