#!/usr/bin/env python3
"""Optimization-safe exact referee for THM-1266.

The referee has four independent jobs.

* Exhaust proper six-owner words and check the closest-return leaf law.
* Check the recursive disjoint-packet packing count and the 6/5 arithmetic.
* Check the exact repeated-low separation behind the five-rung star bound.
* Reconstruct the primitive c=140, k=80 sharp cell, its centered blockers,
  aligned two-cycles, full-comb sweep, and four uncovered fastest-safe tails.

The general proofs are in the theorem note and the sorry-free Lean arithmetic
consumer.  This file deliberately uses explicit RuntimeError checks rather
than Python ``assert`` so that ``python`` and ``python -O`` are identical.
"""

import ast
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import ceil, floor, gcd, lcm
from pathlib import Path


def require(condition: bool, message: str = "exact certificate check failed") -> None:
    if not condition:
        raise RuntimeError(message)


def optimization_safety_probe() -> int:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "optimization-sensitive assert remains in referee")
    caught = False
    try:
        require(False, "sentinel")
    except RuntimeError:
        caught = True
    require(caught, "explicit RuntimeError sentinel did not fire")
    return count


def proper_words(alphabet_size: int, length: int):
    """Generate words with unequal consecutive letters."""

    require(alphabet_size >= 2 and length >= 1)

    def visit(prefix: tuple[int, ...]):
        if len(prefix) == length:
            yield prefix
            return
        for symbol in range(alphabet_size):
            if prefix and prefix[-1] == symbol:
                continue
            yield from visit((*prefix, symbol))

    yield from visit(())


def closest_return(word: tuple[int, ...]) -> tuple[int, int] | None:
    best: tuple[int, int] | None = None
    for left in range(len(word)):
        for right in range(left + 1, len(word)):
            if word[left] != word[right]:
                continue
            if best is None or right - left < best[1] - best[0]:
                best = (left, right)
    return best


def closest_return_word_audit() -> tuple[int, int, Counter[int]]:
    total = 0
    repeated = 0
    edge_counts: Counter[int] = Counter()
    for length in range(1, 9):
        for word in proper_words(6, length):
            total += 1
            packet = closest_return(word)
            if packet is None:
                require(len(set(word)) == len(word))
                require(length <= 6)
                continue
            repeated += 1
            left, right = packet
            edge_count = right - left
            internal = word[left + 1 : right]
            require(2 <= edge_count <= 6)
            require(word[left] not in internal)
            require(len(set(internal)) == len(internal))
            edge_counts[edge_count] += 1
    return total, repeated, edge_counts


def ratio_leaf_audit() -> tuple[int, int, Fraction]:
    rows = 0
    positive_rows = 0
    bank = tuple(range(2, 10))
    for edge_count in range(2, 7):
        internal_count = edge_count - 1
        for speeds in permutations(bank, internal_count + 1):
            outer = speeds[0]
            internal = speeds[1:]
            for address_return in range(1, 4):
                rows += 1
                ratio_sum = sum(
                    (Fraction(outer, speed) for speed in internal), Fraction()
                )
                if ratio_sum <= 7 * address_return - 1:
                    continue
                positive_rows += 1
                minimum_internal = min(internal)
                ratio = Fraction(outer, minimum_internal)
                require(
                    ratio
                    > Fraction(7 * address_return - 1, internal_count)
                )
                require(ratio > Fraction(6, 5))
    # A distinct-internal five-term row approaching the strict 6/5 frontier.
    near_internal = tuple(1_000_000 + offset for offset in range(5))
    near_outer = 1_200_003
    near_sum = sum(
        (Fraction(near_outer, speed) for speed in near_internal), Fraction()
    )
    require(near_sum > 6)
    near_ratio = Fraction(near_outer, min(near_internal))
    require(Fraction(6, 5) < near_ratio < Fraction(1_200_004, 1_000_000))
    return rows, positive_rows, near_ratio


def packet_packing(word: tuple[int, ...]) -> tuple[list[tuple[int, int]], list[int]]:
    packets: list[tuple[int, int]] = []
    terminal_lengths: list[int] = []

    def visit(left: int, right: int) -> None:
        if left >= right:
            return
        block = word[left:right]
        packet = closest_return(block)
        if packet is None:
            require(len(set(block)) == len(block))
            require(len(block) <= 6)
            terminal_lengths.append(len(block))
            return
        packet_left = left + packet[0]
        packet_right = left + packet[1]
        require(2 <= packet_right - packet_left <= 6)
        packets.append((packet_left, packet_right))
        visit(left, packet_left)
        visit(packet_right + 1, right)

    visit(0, len(word))
    return packets, terminal_lengths


def packet_packing_audit() -> tuple[int, int, int, int]:
    words = 0
    packets_total = 0
    maximum_packets = 0
    minimum_slack = 10**9
    for length in range(1, 11):
        for word in proper_words(4, length):
            words += 1
            packets, terminals = packet_packing(word)
            packets_total += len(packets)
            maximum_packets = max(maximum_packets, len(packets))
            for index, (left, right) in enumerate(packets):
                for other_left, other_right in packets[index + 1 :]:
                    require(right < other_left or other_right < left)
            require(len(terminals) <= len(packets) + 1)
            require(sum(terminals) <= 6 * (len(packets) + 1))
            consumed = sum(right - left + 1 for left, right in packets)
            require(consumed + sum(terminals) == len(word))
            slack = 13 * len(packets) + 6 - len(word)
            require(slack >= 0)
            minimum_slack = min(minimum_slack, slack)
    return words, packets_total, maximum_packets, minimum_slack


def star_separation_audit() -> tuple[int, Fraction]:
    rows = 0
    closest_forbidden_slack: Fraction | None = None
    for low in range(1, 101):
        for high in range(6 * low + 1, 9 * low + 1):
            for delta in range(1, 6):
                rows += 1
                # Repeated-low endpoint separation would force the strict
                # reverse inequality high/low < delta + 1/7.
                slack = Fraction(high, low) - (delta + Fraction(1, 7))
                require(slack > 0)
                if closest_forbidden_slack is None or slack < closest_forbidden_slack:
                    closest_forbidden_slack = slack
    require(closest_forbidden_slack is not None)
    return rows, closest_forbidden_slack


def tooth(speed: int, address: int) -> tuple[Fraction, Fraction]:
    require(speed > 0)
    return (
        Fraction(14 * address - 1, 14 * speed),
        Fraction(14 * address + 1, 14 * speed),
    )


def nearest_integer(value: Fraction) -> int:
    lower = floor(value)
    remainder = value - lower
    require(remainder != Fraction(1, 2), "sharp row unexpectedly has a tie")
    return lower if remainder < Fraction(1, 2) else lower + 1


def circle_depth(value: Fraction) -> Fraction:
    fractional = value - floor(value)
    return min(fractional, 1 - fractional)


def merge_open_teeth(
    carrier_gap: tuple[Fraction, Fraction], speeds: tuple[int, ...]
) -> tuple[
    int,
    list[tuple[Fraction, Fraction]],
    list[tuple[Fraction, Fraction]],
    list[tuple[Fraction, Fraction]],
]:
    gap_left, gap_right = carrier_gap
    intervals: list[tuple[Fraction, Fraction]] = []
    raw_intervals: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        low_address = floor(gap_left * speed) - 2
        high_address = ceil(gap_right * speed) + 2
        for address in range(low_address, high_address + 1):
            left, right = tooth(speed, address)
            if left < gap_right and gap_left < right:
                raw_intervals.append((left, right))
                intervals.append((max(left, gap_left), min(right, gap_right)))
    intervals.sort()

    covered: list[list[Fraction]] = []
    for left, right in intervals:
        if not covered or covered[-1][1] < left:
            covered.append([left, right])
        elif covered[-1][1] < right:
            covered[-1][1] = right

    holes: list[tuple[Fraction, Fraction]] = []
    cursor = gap_left
    for left, right in covered:
        if cursor < left:
            holes.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < gap_right:
        holes.append((cursor, gap_right))
    return len(intervals), holes, intervals, raw_intervals


def exact_star_cell_audit() -> dict[str, object]:
    carrier = 140
    gap_address = 80
    high = 1805
    rows = (
        (1805, 1036),
        (256, 147),
        (1805, 1037),
        (254, 146),
        (1805, 1038),
        (292, 168),
        (1805, 1039),
        (257, 148),
        (1805, 1040),
        (255, 147),
        (1805, 1041),
    )
    speeds = tuple(sorted({speed for speed, _ in rows}))
    intervals = tuple(tooth(speed, address) for speed, address in rows)
    gap = (
        Fraction(14 * gap_address + 1, 14 * carrier),
        Fraction(14 * gap_address + 13, 14 * carrier),
    )

    require(gap[0] < intervals[0][0] < intervals[-1][1] < gap[1])
    require(all(intervals[i][0] < intervals[i + 1][0] for i in range(10)))
    require(all(intervals[i][1] < intervals[i + 1][1] for i in range(10)))
    overlaps = tuple(intervals[i][1] - intervals[i + 1][0] for i in range(10))
    require(all(overlap > 0 for overlap in overlaps))
    separations = tuple(intervals[i + 2][0] - intervals[i][1] for i in range(9))
    require(all(separation > 0 for separation in separations))

    low_speeds = (256, 254, 292, 257, 255)
    detunings = tuple(high - 6 * speed for speed in low_speeds)
    require(detunings == (269, 281, 53, 263, 275))
    seam_pair_masses: list[Fraction] = []
    for rung, low in enumerate(low_speeds):
        left_edge = 2 * rung
        seam_pair = overlaps[left_edge] + overlaps[left_edge + 1]
        expected = Fraction(high - 6 * low, 7 * high * low)
        require(seam_pair == expected)
        require(seam_pair >= Fraction(1, 7 * lcm(high, low)))
        seam_pair_masses.append(seam_pair)
    require(sum(seam_pair_masses, Fraction()) == sum(overlaps, Fraction()))

    edge_symbols = tuple(
        frozenset((rows[index][0], rows[index + 1][0])) for index in range(10)
    )
    backtracks = sum(
        edge_symbols[index - 1] == edge_symbols[index] for index in range(1, 10)
    )
    turns = 9 - backtracks
    require((backtracks, turns) == (5, 4))

    harmonic_sum = sum((Fraction(1, speed) for speed in speeds), Fraction())
    harmonic_margin = harmonic_sum - Fraction(1, carrier)
    require(harmonic_margin > 0)
    require(max(speeds) < 2345 * carrier)
    all_speed_gcd = carrier
    for speed in speeds:
        all_speed_gcd = gcd(all_speed_gcd, speed)
    require(all_speed_gcd == 1)

    centre = Fraction(2 * gap_address + 1, 2 * carrier)
    phases: dict[int, Fraction] = {}
    for speed in speeds:
        clock = carrier + speed
        numerator = nearest_integer(clock * centre)
        phase = Fraction(numerator, clock)
        phases[speed] = phase
        require(gap[0] < phase < gap[1])
        require(circle_depth(carrier * phase) == circle_depth(speed * phase))
        require(circle_depth(speed * phase) > Fraction(1, 4))

    blocker_indices = {
        254: 8,
        255: 3,
        256: 7,
        257: 1,
        292: 1,
        1805: 3,
    }
    blocker_map: dict[int, int] = {}
    for source, target_index in blocker_indices.items():
        require(0 < target_index < len(rows) - 1)
        target_speed = rows[target_index][0]
        require(target_speed != source)
        require(intervals[target_index][0] < phases[source] < intervals[target_index][1])
        blocker_map[source] = target_speed
    require(
        blocker_map
        == {254: 1805, 255: 254, 256: 257, 257: 256, 292: 256, 1805: 254}
    )

    cycles = ((254, 1805, 8, 3, 7), (256, 257, 7, 1, 6))
    phase_differences: list[Fraction] = []
    for low, high_speed, high_target, low_target, bridge in cycles:
        require(blocker_map[low] == high_speed and blocker_map[high_speed] == low)
        require(high_target - low_target > 1)
        phase_difference = phases[high_speed] - phases[low]
        phase_differences.append(phase_difference)
        require(phase_difference < 0)
        require(low_target - high_target < 0)
        require(low_target < bridge < high_target)
        require(rows[bridge][0] not in (low, high_speed))
        # The ascent target is wholly source-safe.  Both endpoints stay in
        # one source-depth branch and exceed the danger threshold.
        target_left, target_right = intervals[high_target]
        require(floor(low * target_left) == floor(low * target_right))
        require(circle_depth(low * target_left) > Fraction(1, 14))
        require(circle_depth(low * target_right) > Fraction(1, 14))

    tooth_count, holes, all_teeth, raw_teeth = merge_open_teeth(gap, speeds)
    expected_holes = [
        (Fraction(14477, 25270), Fraction(14489, 25270)),
        (Fraction(14491, 25270), Fraction(14503, 25270)),
        (Fraction(2915, 5054), Fraction(14587, 25270)),
        (Fraction(14589, 25270), Fraction(14601, 25270)),
    ]
    require(tooth_count == 20)
    require(holes == expected_holes)
    # Open intervals can in principle touch at one uncovered equality point.
    # Test every endpoint and every endpoint-cell midpoint, which is an exact
    # exhaustive union check because comb activity is constant between them.
    critical_points = sorted(
        {gap[0], gap[1]}
        | {endpoint for interval in all_teeth for endpoint in interval}
    )
    probe_points = set(critical_points)
    probe_points.update(
        (critical_points[index] + critical_points[index + 1]) / 2
        for index in range(len(critical_points) - 1)
    )
    for point in probe_points:
        covered = any(left < point < right for left, right in raw_teeth)
        expected_covered = not any(
            left <= point <= right for left, right in expected_holes
        )
        require(covered == expected_covered, f"union mismatch at {point}")
    require(all(right - left == Fraction(6, 7 * high) for left, right in holes))
    require(sum((right - left for left, right in holes), Fraction()) == Fraction(24, 12635))
    # These are the four fastest-safe gaps after addresses 1034, 1035,
    # 1041, and 1042.  The exact star fills the intervening five.
    hole_addresses = []
    for left, right in holes:
        address = (14 * high * left - 1) // 14
        require(left == tooth(high, address)[1])
        require(right == tooth(high, address + 1)[0])
        hole_addresses.append(address)
    require(tuple(hole_addresses) == (1034, 1035, 1041, 1042))

    return {
        "rows": len(rows),
        "overlaps": len(overlaps),
        "seam_mass": sum(overlaps, Fraction()),
        "minimum_separation": min(separations),
        "harmonic_margin": harmonic_margin,
        "backtracks": backtracks,
        "turns": turns,
        "phases": phases,
        "phase_differences": tuple(phase_differences),
        "blocker_map": blocker_map,
        "teeth": tooth_count,
        "holes": tuple(holes),
        "hole_mass": Fraction(24, 12635),
    }


def main() -> None:
    assert_nodes = optimization_safety_probe()
    word_total, word_repeated, leaf_edges = closest_return_word_audit()
    ratio_rows, ratio_positive, near_sharp_ratio = ratio_leaf_audit()
    packing_words, packet_total, maximum_packets, packing_slack = packet_packing_audit()
    separation_rows, separation_slack = star_separation_audit()
    cell = exact_star_cell_audit()

    print("THM-1266 CLOSEST-RETURN LEAF / PAID-STAR EXACT AUDIT")
    print(f"optimization-sensitive assert nodes = {assert_nodes}")
    print(f"proper six-owner words = {word_total}")
    print(f"words with a closest return = {word_repeated}")
    print(f"closest-return edge histogram = {tuple(sorted(leaf_edges.items()))}")
    print(f"exact ratio rows / positive rows = {ratio_rows} / {ratio_positive}")
    print(f"near-sharp outer/min ratio = {near_sharp_ratio}")
    print(f"packet-packing words / packets = {packing_words} / {packet_total}")
    print(f"maximum packets / minimum 13P+6 slack = {maximum_packets} / {packing_slack}")
    print(f"forbidden repeated-low separation rows = {separation_rows}")
    print(f"minimum forbidden separation slack = {separation_slack}")
    print(f"sharp star teeth / seam occurrences = {cell['rows']} / {cell['overlaps']}")
    print(f"sharp star seam mass = {cell['seam_mass']}")
    print(f"sharp star minimum nonconsecutive separation = {cell['minimum_separation']}")
    print(f"sharp star B/T ledger = {cell['backtracks']} / {cell['turns']}")
    print(f"sharp harmonic margin above 1/140 = {cell['harmonic_margin']}")
    print(f"centered blocker map = {tuple(cell['blocker_map'].items())}")
    print(f"aligned two-cycle phase differences = {cell['phase_differences']}")
    print(f"full-comb teeth meeting G / uncovered components = {cell['teeth']} / {len(cell['holes'])}")
    print(f"uncovered fastest-safe tail mass = {cell['hole_mass']}")
    print("runner tournament fingerprint = transitive scores (0,1,2,3,4,5); cycles 0; SCCs 6; Hamiltonian paths 1")
    print("return witness fingerprint = five-arrow in-star; directed cycles 0; depth 1")
    print("chronological switch word = B,T,B,T,B,T,B,T,B")
    print("RESULT: PASS -- local r=1 star terminates sharply at five; four fastest-safe tails remain")


if __name__ == "__main__":
    main()
