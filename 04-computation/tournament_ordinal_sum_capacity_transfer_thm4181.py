#!/usr/bin/env python3
"""Exact audit for THM-4181 ordinal-sum capacity transfer.

The large census uses one ``gentourng`` representative per unlabelled
tournament class.  It counts ordered block-class presentations (A, B), not
orbits or isomorphism classes of the ordinal-sum child.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from dataclasses import dataclass


EXPECTED_CLASSES = {2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
EXPECTED_NO_SINK = {3: 1, 4: 2, 5: 8, 6: 44, 7: 400}


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def find_gentourng() -> str:
    configured = os.environ.get("THM4181_GENTOURNG")
    candidates = (
        configured,
        shutil.which("gentourng"),
        "/opt/homebrew/bin/gentourng",
        "/usr/local/bin/gentourng",
    )
    for candidate in candidates:
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    raise RuntimeError("gentourng executable not found")


def parse(bits: str) -> tuple[int, ...]:
    order = 0
    while order * (order - 1) // 2 < len(bits):
        order += 1
    need(order * (order - 1) // 2 == len(bits), "invalid pair-bit label")
    out = [0] * order
    cursor = 0
    for left in range(order):
        for right in range(left + 1, order):
            if bits[cursor] == "1":
                out[left] |= 1 << right
            else:
                out[right] |= 1 << left
            cursor += 1
    return tuple(out)


def label(out: tuple[int, ...]) -> str:
    return "".join(
        "1" if out[left] & (1 << right) else "0"
        for left in range(len(out))
        for right in range(left + 1, len(out))
    )


def transitive(order: int) -> tuple[int, ...]:
    return tuple(sum(1 << right for right in range(left + 1, order))
                 for left in range(order))


def has_sink(out: tuple[int, ...]) -> bool:
    return any(row == 0 for row in out)


def classes(executable: str, order: int) -> tuple[tuple[int, ...], ...]:
    process = subprocess.run(
        [executable, "-q", str(order)],
        check=False,
        capture_output=True,
        text=True,
    )
    need(process.returncode == 0, f"gentourng exit at order {order}")
    need(not process.stderr.strip(), f"gentourng stderr at order {order}")
    records = tuple(line.strip() for line in process.stdout.splitlines() if line.strip())
    need(len(records) == EXPECTED_CLASSES[order], f"class count at order {order}")
    need(len(set(records)) == len(records), f"duplicate class at order {order}")
    return tuple(parse(record) for record in records)


@dataclass(frozen=True)
class TournamentData:
    out: tuple[int, ...]
    hamilton: int
    capacities: tuple[tuple[int, ...], ...]
    starts: tuple[tuple[int, int], ...]
    ends: tuple[tuple[int, int], ...]
    subset_hamilton: tuple[int, ...]
    mass: int
    gate: int


def gate(out: tuple[int, ...], capacities: tuple[tuple[int, ...], ...]) -> int:
    order = len(out)
    degrees = [sum(row) for row in capacities]
    currents = [
        sum(
            capacities[vertex][other]
            if out[vertex] & (1 << other)
            else -capacities[vertex][other]
            for other in range(order)
            if other != vertex
        )
        for vertex in range(order)
    ]
    current = sum(degrees[i] * currents[i] for i in range(order))
    mass = sum(capacities[i][j] for i in range(order) for j in range(i + 1, order))
    squares = sum(capacities[i][j] ** 2
                  for i in range(order) for j in range(i + 1, order))
    disjoint = (mass * mass + squares - sum(value * value for value in degrees)) // 2
    return disjoint + 2 * current


def tournament_data(out: tuple[int, ...]) -> TournamentData:
    order = len(out)
    size = 1 << order
    full = size - 1
    start = [[0] * order for _ in range(size)]
    end = [[0] * order for _ in range(size)]
    for vertex in range(order):
        start[1 << vertex][vertex] = 1
        end[1 << vertex][vertex] = 1
    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        todo = mask
        while todo:
            bit = todo & -todo
            todo ^= bit
            vertex = bit.bit_length() - 1
            rest = mask ^ bit
            candidates = rest
            while candidates:
                other_bit = candidates & -candidates
                candidates ^= other_bit
                other = other_bit.bit_length() - 1
                if out[vertex] & other_bit:
                    start[mask][vertex] += start[rest][other]
                if out[other] & bit:
                    end[mask][vertex] += end[rest][other]
    subset_hamilton_list = [sum(end[mask]) for mask in range(size)]
    # The empty complement carries the unique empty Hamilton word.
    subset_hamilton_list[0] = 1
    subset_hamilton = tuple(subset_hamilton_list)
    exposed = [[0] * order for _ in range(order)]
    for left_mask in range(1, full):
        right_mask = full ^ left_mask
        left_vertices = left_mask
        while left_vertices:
            left_bit = left_vertices & -left_vertices
            left_vertices ^= left_bit
            left = left_bit.bit_length() - 1
            left_count = end[left_mask][left]
            if not left_count:
                continue
            right_vertices = right_mask
            while right_vertices:
                right_bit = right_vertices & -right_vertices
                right_vertices ^= right_bit
                right = right_bit.bit_length() - 1
                exposed[left][right] += left_count * start[right_mask][right]
    capacities = [[0] * order for _ in range(order)]
    for left in range(order):
        for right in range(left + 1, order):
            value = exposed[left][right] + exposed[right][left]
            capacities[left][right] = capacities[right][left] = value
    rooted_start = [[0, 0] for _ in range(order)]
    rooted_end = [[0, 0] for _ in range(order)]
    for mask in range(1, size):
        parity = mask.bit_count() & 1
        complement_h = subset_hamilton[full ^ mask]
        if not complement_h:
            continue
        for vertex in range(order):
            rooted_start[vertex][parity] += start[mask][vertex] * complement_h
            rooted_end[vertex][parity] += end[mask][vertex] * complement_h
    capacity_tuple = tuple(tuple(row) for row in capacities)
    mass = sum(capacities[i][j] for i in range(order) for j in range(i + 1, order))
    return TournamentData(
        out=out,
        hamilton=subset_hamilton[full],
        capacities=capacity_tuple,
        starts=tuple(tuple(row) for row in rooted_start),
        ends=tuple(tuple(row) for row in rooted_end),
        subset_hamilton=subset_hamilton,
        mass=mass,
        gate=gate(out, capacity_tuple),
    )


def ordinal_out(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    nleft = len(left)
    nright = len(right)
    answer = [0] * (nleft + nright)
    right_mask = ((1 << nright) - 1) << nleft
    for vertex in range(nleft):
        answer[vertex] = left[vertex] | right_mask
    for vertex in range(nright):
        answer[nleft + vertex] = right[vertex] << nleft
    return tuple(answer)


@dataclass(frozen=True)
class OrdinalPacket:
    out: tuple[int, ...]
    capacities: tuple[tuple[int, ...], ...]
    total_gate: int
    left_gate: int
    right_gate: int
    cross_self: int
    left_right: int
    left_cross: int
    right_cross: int
    remainder: int


def ordinal_packet(left: TournamentData, right: TournamentData) -> OrdinalPacket:
    nleft = len(left.out)
    nright = len(right.out)
    order = nleft + nright
    capacities = [[0] * order for _ in range(order)]
    for i in range(nleft):
        for j in range(i + 1, nleft):
            value = right.hamilton * left.capacities[i][j]
            capacities[i][j] = capacities[j][i] = value
    for i in range(nright):
        for j in range(i + 1, nright):
            value = left.hamilton * right.capacities[i][j]
            capacities[nleft + i][nleft + j] = capacities[nleft + j][nleft + i] = value
    cross = [[0] * nright for _ in range(nleft)]
    for i in range(nleft):
        for j in range(nright):
            value = 2 * (
                left.starts[i][0] * right.ends[j][0]
                + left.starts[i][1] * right.ends[j][1]
            )
            cross[i][j] = value
            capacities[i][nleft + j] = capacities[nleft + j][i] = value
    out = ordinal_out(left.out, right.out)
    capacity_tuple = tuple(tuple(row) for row in capacities)
    total_gate = gate(out, capacity_tuple)
    left_gate = right.hamilton ** 2 * left.gate
    right_gate = left.hamilton ** 2 * right.gate
    wx = right.hamilton * left.mass
    wy = left.hamilton * right.mass
    p = [sum(row) for row in cross]
    q = [sum(cross[i][j] for i in range(nleft)) for j in range(nright)]
    total_cross = sum(p)
    cross_squares = sum(value * value for row in cross for value in row)
    numerator = (
        total_cross * total_cross
        + cross_squares
        + 3 * sum(value * value for value in p)
        - 5 * sum(value * value for value in q)
    )
    need(numerator % 2 == 0, "cross self-gate parity")
    cross_self = numerator // 2
    left_right = wx * wy
    left_cross = 0
    for i in range(nleft):
        degree = right.hamilton * sum(left.capacities[i])
        outgoing = right.hamilton * sum(
            left.capacities[i][j]
            for j in range(nleft)
            if left.out[i] & (1 << j)
        )
        left_cross += p[i] * (wx - degree + 4 * outgoing)
    right_cross = 0
    for j in range(nright):
        degree = left.hamilton * sum(right.capacities[j])
        incoming = left.hamilton * sum(
            right.capacities[u][j]
            for u in range(nright)
            if right.out[u] & (1 << j)
        )
        right_cross += q[j] * (wy - degree - 4 * incoming)
    formula_gate = (
        left_gate + right_gate + cross_self
        + left_right + left_cross + right_cross
    )
    need(formula_gate == total_gate, "generic block-gate formula")
    remainder = formula_gate - left_gate - right_gate
    return OrdinalPacket(
        out=out,
        capacities=capacity_tuple,
        total_gate=total_gate,
        left_gate=left_gate,
        right_gate=right_gate,
        cross_self=cross_self,
        left_right=left_right,
        left_cross=left_cross,
        right_cross=right_cross,
        remainder=remainder,
    )


def odd_path_layers(data: TournamentData) -> dict[int, tuple[tuple[int, ...], ...]]:
    order = len(data.out)
    full = (1 << order) - 1
    layers = {
        length: [[0] * order for _ in range(order)]
        for length in range(1, order, 2)
    }

    def visit(start: int, last: int, mask: int) -> None:
        length = mask.bit_count() - 1
        if length & 1:
            value = 2 * data.subset_hamilton[full ^ mask]
            left, right = sorted((start, last))
            layers[length][left][right] += value
            layers[length][right][left] += value
        candidates = data.out[last] & ~mask & full
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            visit(start, bit.bit_length() - 1, mask | bit)

    for vertex in range(order):
        visit(vertex, vertex, 1 << vertex)
    return {length: tuple(tuple(row) for row in tensor)
            for length, tensor in layers.items()}


def dominance_rows(left: TournamentData, right: TournamentData,
                   packet: OrdinalPacket):
    nleft = len(left.out)
    left_rows = []
    right_rows = []
    for tail in range(nleft):
        for head in range(nleft):
            if left.out[tail] & (1 << head):
                for vertex in range(len(right.out)):
                    old = packet.capacities[tail][head]
                    cross = packet.capacities[head][nleft + vertex]
                    left_rows.append((old - cross, tail, head, vertex, old, cross))
    for tail in range(len(right.out)):
        for head in range(len(right.out)):
            if right.out[tail] & (1 << head):
                for vertex in range(nleft):
                    cross = packet.capacities[vertex][nleft + tail]
                    old = packet.capacities[nleft + tail][nleft + head]
                    right_rows.append((cross - old, vertex, tail, head, cross, old))
    return left_rows, right_rows


def main() -> None:
    executable = find_gentourng()
    tool_digest = hashlib.sha256()
    with open(executable, "rb") as tool_stream:
        for block in iter(lambda: tool_stream.read(1 << 20), b""):
            tool_digest.update(block)
    banks = {1: (transitive(1),)}
    banks.update({order: classes(executable, order) for order in range(2, 8)})
    no_sink = {
        order: tuple(out for out in banks[order] if not has_sink(out))
        for order in range(3, 8)
    }
    for order, expected in EXPECTED_NO_SINK.items():
        need(len(no_sink[order]) == expected, f"no-sink class count at order {order}")
    data_banks = {
        order: tuple(tournament_data(out) for out in bank)
        for order, bank in banks.items()
    }
    no_sink_data = {
        order: tuple(tournament_data(out) for out in no_sink[order])
        for order in range(3, 8)
    }

    sidecar_checks = 0
    layer_coordinates = 0
    exchange_mass_checks = 0
    for order, bank in data_banks.items():
        for data in bank:
            start_totals = [sum(row[parity] for row in data.starts) for parity in (0, 1)]
            end_totals = [sum(row[parity] for row in data.ends) for parity in (0, 1)]
            need(sum(start_totals) == data.mass + data.hamilton, "start cover mass")
            need(sum(end_totals) == data.mass + data.hamilton, "end cover mass")
            if order & 1:
                need(start_totals[1] - start_totals[0] == data.hamilton,
                     "odd start component exchange")
                need(end_totals[1] - end_totals[0] == data.hamilton,
                     "odd end component exchange")
            layers = odd_path_layers(data)
            rebuilt = [[0] * order for _ in range(order)]
            for tensor in layers.values():
                for i in range(order):
                    for j in range(order):
                        rebuilt[i][j] += tensor[i][j]
                        layer_coordinates += i < j
            need(tuple(tuple(row) for row in rebuilt) == data.capacities,
                 "odd-path layer reconstruction")
            if order % 2 == 0:
                for length in range(1, order - 2, 2):
                    partner = order - length - 2
                    left_mass = sum(layers[length][i][j]
                                    for i in range(order) for j in range(i + 1, order))
                    right_mass = sum(layers[partner][i][j]
                                     for i in range(order) for j in range(i + 1, order))
                    need(left_mass == right_mass, "even component-exchange mass")
                    exchange_mass_checks += 1
            sidecar_checks += 1

    direct_transfer_checks = 0
    ordinal_presentations = 0
    negative_remainders = 0
    zero_remainders = 0
    minimum = None
    first_left = None
    first_right = None
    first_right_with_left_at_least = {2: None, 3: None}
    left_block_sharp = {}
    right_block_sharp = {}
    sharp_left = None
    sharp_right = None
    stress_left = None
    stress_right = None
    digest = hashlib.sha256()
    for nleft in range(1, 8):
        for nright in range(3, 8):
            for left in data_banks[nleft]:
                for right in no_sink_data[nright]:
                    packet = ordinal_packet(left, right)
                    witness = (nleft, nright, label(left.out), label(right.out), packet.remainder)
                    if minimum is None or (packet.remainder, witness) < (minimum[0], minimum[1]):
                        minimum = (packet.remainder, witness)
                    negative_remainders += packet.remainder < 0
                    zero_remainders += packet.remainder == 0
                    digest.update(("|".join(map(str, witness)) + "\n").encode("ascii"))
                    if nleft <= 4 and nright <= 4:
                        actual = tournament_data(packet.out).capacities
                        need(actual == packet.capacities, "direct ordinal capacity transfer")
                        direct_transfer_checks += 1
                    left_rows, right_rows = dominance_rows(left, right, packet)
                    block_key = (nleft + nright, nleft, nright)
                    left_label = label(left.out)
                    right_label = label(right.out)
                    for row in left_rows:
                        candidate = (row[0], left_label, right_label, row[1:4], row)
                        if block_key not in left_block_sharp or candidate < left_block_sharp[block_key]:
                            left_block_sharp[block_key] = candidate
                        if row[0] < 0:
                            key = block_key + (left_label, right_label, row[1:4])
                            first_candidate = (key, row, left_label, right_label)
                            if first_left is None or first_candidate < first_left:
                                first_left = first_candidate
                        if sharp_left is None or candidate < sharp_left:
                            sharp_left = candidate
                        if left_label == "1100111111" and right_label == "1100111111":
                            if stress_left is None or candidate < stress_left:
                                stress_left = candidate
                    for row in right_rows:
                        candidate = (row[0], left_label, right_label, row[1:4], row)
                        if block_key not in right_block_sharp or candidate < right_block_sharp[block_key]:
                            right_block_sharp[block_key] = candidate
                        if row[0] < 0:
                            key = block_key + (left_label, right_label, row[1:4])
                            first_candidate = (key, row, left_label, right_label)
                            if first_right is None or first_candidate < first_right:
                                first_right = first_candidate
                            for threshold in (2, 3):
                                if nleft >= threshold:
                                    old = first_right_with_left_at_least[threshold]
                                    if old is None or first_candidate < old:
                                        first_right_with_left_at_least[threshold] = first_candidate
                        if sharp_right is None or candidate < sharp_right:
                            sharp_right = candidate
                        if left_label == "1110101111" and right_label == "1111100111":
                            if stress_right is None or candidate < stress_right:
                                stress_right = candidate
                    ordinal_presentations += 1

    need(ordinal_presentations == 242_060, "ordinal presentation count")
    need(negative_remainders == 0 and zero_remainders == 0,
         "strict finite ordinal remainder")
    need(minimum == (120, (1, 3, "", "101", 120)), "minimum remainder")
    need(first_left is not None and first_left[0][:3] == (6, 3, 3),
         "first left-domination obstruction")
    need(first_left[1] == (-4, 0, 1, 0, 6, 10), "left obstruction row")
    need(first_right is not None and first_right[0][:3] == (7, 1, 6),
         "first right-domination obstruction")
    need(first_right[1] == (-6, 0, 0, 5, 22, 28), "first right obstruction row")
    first_left_block_sharp = left_block_sharp[first_left[0][:3]]
    first_right_block_sharp = right_block_sharp[first_right[0][:3]]
    need(first_left_block_sharp[4] == (-4, 0, 1, 0, 6, 10),
         "first-block left sharp obstruction")
    need(first_right_block_sharp[0] <= -6, "first-block right sharp obstruction")
    need(first_right_block_sharp[1:3] == ("", "111111101111110")
         and first_right_block_sharp[4] == (-8, 0, 0, 5, 18, 26),
         "lex-first sharp row in first right-failing block")
    need(first_right_with_left_at_least[2] is not None
         and first_right_with_left_at_least[2][0][:3] == (8, 2, 6),
         "first right obstruction with left order at least two")
    need(first_right_with_left_at_least[3] is not None
         and first_right_with_left_at_least[3][0][:3] == (8, 3, 5),
         "first right obstruction with left order at least three")
    need(stress_left is not None and stress_left[0] == -642, "5x5 left stress control")
    need(stress_right is not None and stress_right[0] == -60, "5x5 right stress control")

    print("THM4181_PRIMARY_EXACT_AUDIT")
    print("gentourng_sha256", tool_digest.hexdigest())
    print("class_counts", " ".join(f"q{n}={len(banks[n])}" for n in range(1, 8)))
    print("no_sink_counts", " ".join(f"q{n}={len(no_sink[n])}" for n in range(3, 8)))
    print("sidecar_checks", sidecar_checks)
    print("odd_path_layer_coordinates", layer_coordinates)
    print("even_component_exchange_mass_checks", exchange_mass_checks)
    print("direct_transfer_checks", direct_transfer_checks)
    print("ordinal_block_class_presentations", ordinal_presentations)
    print("remainder_negative", negative_remainders, "zero", zero_remainders)
    print("remainder_minimum", minimum)
    print("first_left_obstruction", first_left)
    print("first_right_obstruction", first_right)
    print("first_right_with_left_order_at_least_2", first_right_with_left_at_least[2])
    print("first_right_with_left_order_at_least_3", first_right_with_left_at_least[3])
    print("first_left_block_order_sharp", first_left_block_sharp)
    print("first_right_block_order_sharp", first_right_block_sharp)
    print("global_left_sharp", sharp_left)
    print("global_right_sharp", sharp_right)
    print("five_by_five_left_stress", stress_left)
    print("five_by_five_right_stress", stress_right)
    print("presentation_stream_sha256", digest.hexdigest())
    print("PASS")


if __name__ == "__main__":
    main()
