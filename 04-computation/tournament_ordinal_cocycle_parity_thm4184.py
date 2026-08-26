#!/usr/bin/env python3
"""Exact primary audit for THM-4184 ordinal parity and the R_+ cocycle.

The class census uses one nauty ``gentourng`` representative per unlabelled
tournament class.  Its rows are ordered factor-class presentations, not
isomorphism classes of the ordinal-sum child.
"""

from __future__ import annotations

import hashlib
import itertools
import os
import shutil
import subprocess
from dataclasses import dataclass
from fractions import Fraction


EXPECTED_CLASSES = {2: 1, 3: 2, 4: 4, 5: 12, 6: 56}
EXPECTED_NO_SINK = {3: 1, 4: 2, 5: 8, 6: 44}


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def find_gentourng() -> str:
    configured = os.environ.get("ORDINAL_COCYCLE_GENTOURNG")
    for candidate in (
        configured,
        shutil.which("gentourng"),
        "/opt/homebrew/bin/gentourng",
        "/usr/local/bin/gentourng",
    ):
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    raise RuntimeError("gentourng executable not found")


def parse(bits: str, order: int) -> tuple[int, ...]:
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
    return tuple(
        sum(1 << right for right in range(left + 1, order))
        for left in range(order)
    )


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
    return tuple(parse(record, order) for record in records)


def parity_convolution(
    left: tuple[int, int], right: tuple[int, int]
) -> tuple[int, int]:
    return (
        left[0] * right[0] + left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


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
    mass = sum(
        capacities[i][j]
        for i in range(order)
        for j in range(i + 1, order)
    )
    squares = sum(
        capacities[i][j] ** 2
        for i in range(order)
        for j in range(i + 1, order)
    )
    disjoint = (mass * mass + squares - sum(value * value for value in degrees)) // 2
    return disjoint + 2 * current


@dataclass(frozen=True)
class TournamentData:
    out: tuple[int, ...]
    hamilton: int
    capacities: tuple[tuple[int, ...], ...]
    starts: tuple[tuple[int, int], ...]
    ends: tuple[tuple[int, int], ...]
    mass: int
    gate: int

    @property
    def optional_state(self) -> tuple[int, int]:
        return (
            self.hamilton + sum(row[0] for row in self.starts),
            sum(row[1] for row in self.starts),
        )


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
        for vertex in range(order):
            bit = 1 << vertex
            if not mask & bit:
                continue
            rest = mask ^ bit
            for other in range(order):
                other_bit = 1 << other
                if not rest & other_bit:
                    continue
                if out[vertex] & other_bit:
                    start[mask][vertex] += start[rest][other]
                if out[other] & bit:
                    end[mask][vertex] += end[rest][other]
    subset_hamilton = [sum(end[mask]) for mask in range(size)]
    subset_hamilton[0] = 1
    exposed = [[0] * order for _ in range(order)]
    for left_mask in range(1, full):
        right_mask = full ^ left_mask
        for left in range(order):
            if not left_mask & (1 << left):
                continue
            left_count = end[left_mask][left]
            if not left_count:
                continue
            for right in range(order):
                if right_mask & (1 << right):
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
        for vertex in range(order):
            rooted_start[vertex][parity] += start[mask][vertex] * complement_h
            rooted_end[vertex][parity] += end[mask][vertex] * complement_h
    capacity_tuple = tuple(tuple(row) for row in capacities)
    mass = sum(
        capacities[i][j]
        for i in range(order)
        for j in range(i + 1, order)
    )
    return TournamentData(
        out=out,
        hamilton=subset_hamilton[full],
        capacities=capacity_tuple,
        starts=tuple(tuple(row) for row in rooted_start),
        ends=tuple(tuple(row) for row in rooted_end),
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


def ordinal_data(left: TournamentData, right: TournamentData) -> TournamentData:
    """Compose exact capacities and parity sidecars without child subset DP."""
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
    for i in range(nleft):
        for j in range(nright):
            value = 2 * (
                left.starts[i][0] * right.ends[j][0]
                + left.starts[i][1] * right.ends[j][1]
            )
            capacities[i][nleft + j] = capacities[nleft + j][i] = value
    left_optional = left.optional_state
    right_optional = right.optional_state
    starts = tuple(
        parity_convolution(row, right_optional) for row in left.starts
    ) + tuple(
        (left.hamilton * row[0], left.hamilton * row[1]) for row in right.starts
    )
    ends = tuple(
        (right.hamilton * row[0], right.hamilton * row[1]) for row in left.ends
    ) + tuple(
        parity_convolution(left_optional, row) for row in right.ends
    )
    capacity_tuple = tuple(tuple(row) for row in capacities)
    mass = sum(
        capacities[i][j]
        for i in range(order)
        for j in range(i + 1, order)
    )
    out = ordinal_out(left.out, right.out)
    return TournamentData(
        out=out,
        hamilton=left.hamilton * right.hamilton,
        capacities=capacity_tuple,
        starts=starts,
        ends=ends,
        mass=mass,
        gate=gate(out, capacity_tuple),
    )


def remainder(left: TournamentData, right: TournamentData) -> int:
    child = ordinal_data(left, right)
    return (
        child.gate
        - right.hamilton**2 * left.gate
        - left.hamilton**2 * right.gate
    )


def path_valid(out: tuple[int, ...], path: tuple[int, ...]) -> bool:
    return all(out[path[i]] & (1 << path[i + 1]) for i in range(len(path) - 1))


def terminal_move(
    out: tuple[int, ...], pair: tuple[tuple[int, ...], tuple[int, ...]]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """The all-order sign-reversing involution on ordered two-path covers."""
    left, right = pair
    if not left:
        return ((right[-1],), right[:-1])
    if not right:
        return (left[:-1], (left[-1],))
    if out[left[-1]] & (1 << right[-1]):
        return (left + (right[-1],), right[:-1])
    return (left[:-1], right + (left[-1],))


def involution_audit(data: TournamentData) -> int:
    order = len(data.out)
    covers: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    for word in itertools.permutations(range(order)):
        for cut in range(order + 1):
            pair = (word[:cut], word[cut:])
            if path_valid(data.out, pair[0]) and path_valid(data.out, pair[1]):
                covers.add(pair)
    signed = 0
    for pair in covers:
        image = terminal_move(data.out, pair)
        need(image in covers, "terminal move leaves two-path covers")
        need(terminal_move(data.out, image) == pair, "terminal move involution")
        need((len(pair[0]) ^ len(image[0])) & 1, "terminal move parity toggle")
        signed += -1 if len(pair[0]) & 1 else 1
    need(signed == 0, "two-path-cover signed balance")
    return len(covers)


def direct_three_block_capacity(
    left: TournamentData, middle: TournamentData, right: TournamentData
) -> tuple[tuple[int, ...], ...]:
    """Formula using the retained optional-path state of the middle block."""
    first = ordinal_data(left, middle)
    formula = ordinal_data(first, right).capacities
    second = ordinal_data(middle, right)
    alternate = ordinal_data(left, second).capacities
    need(formula == alternate, "sidecar associativity")
    return formula


def transitive_gate_formula(order: int) -> int:
    if order == 1:
        return 0
    return (
        4 * 4**order
        - 9 * (order + 2) * 2**order
        + 24 * order
        + 32
    ) // 18


def lollipop_gate_formula(prefix_order: int) -> int:
    if prefix_order == 0:
        return 0
    return 2 * (
        37 * 4**prefix_order
        - 9 * (prefix_order + 4) * 2**prefix_order
        + 6 * prefix_order
        - 4
    )


def transitive_capacity_data(order: int) -> tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    out = transitive(order)
    capacities = [[0] * order for _ in range(order)]
    for left in range(order):
        for right in range(left + 1, order):
            distance = right - left
            value = 2 if distance == 1 else 2 ** (distance - 1)
            capacities[left][right] = capacities[right][left] = value
    return out, tuple(tuple(row) for row in capacities)


def lollipop_capacity_data(prefix_order: int) -> tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    cycle = parse("101", 3)
    if prefix_order == 0:
        data = tournament_data(cycle)
        return data.out, data.capacities
    prefix_out, prefix_capacities = transitive_capacity_data(prefix_order)
    out = ordinal_out(prefix_out, cycle)
    order = prefix_order + 3
    capacities = [[0] * order for _ in range(order)]
    for i in range(prefix_order):
        for j in range(prefix_order):
            capacities[i][j] = 3 * prefix_capacities[i][j]
    for i in range(3):
        for j in range(i + 1, 3):
            capacities[prefix_order + i][prefix_order + j] = 2
            capacities[prefix_order + j][prefix_order + i] = 2
    for i in range(prefix_order):
        value = 4 if i == prefix_order - 1 else 3 * 2 ** (prefix_order - i - 1)
        for j in range(3):
            capacities[i][prefix_order + j] = value
            capacities[prefix_order + j][i] = value
    return out, tuple(tuple(row) for row in capacities)


def main() -> None:
    executable = find_gentourng()
    tool_digest = hashlib.sha256()
    with open(executable, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            tool_digest.update(block)

    raw_banks: dict[int, tuple[tuple[int, ...], ...]] = {1: (transitive(1),)}
    raw_banks.update({order: classes(executable, order) for order in range(2, 7)})
    banks = {
        order: tuple(tournament_data(out) for out in raw_banks[order])
        for order in range(1, 7)
    }
    no_sink = {
        order: tuple(data for data in banks[order] if not has_sink(data.out))
        for order in range(3, 7)
    }
    for order, expected in EXPECTED_NO_SINK.items():
        need(len(no_sink[order]) == expected, f"no-sink count at order {order}")

    parity_checks = 0
    involution_objects = 0
    for bank in banks.values():
        for data in bank:
            start_totals = tuple(sum(row[p] for row in data.starts) for p in (0, 1))
            end_totals = tuple(sum(row[p] for row in data.ends) for p in (0, 1))
            need(start_totals == end_totals, "rooted start/end parity totals")
            need(start_totals[1] - start_totals[0] == data.hamilton,
                 "all-order rooted parity balance")
            need(data.optional_state[0] == data.optional_state[1],
                 "balanced optional-path state")
            need(2 * data.optional_state[0] == data.mass + 2 * data.hamilton,
                 "K=(W+2H)/2")
            involution_objects += involution_audit(data)
            parity_checks += 1

    all_data = tuple(data for order in range(1, 7) for data in banks[order])
    all_no_sink = tuple(data for order in range(3, 7) for data in no_sink[order])
    need(len(all_data) == 76, "factor class total")
    need(len(all_no_sink) == 55, "no-sink factor class total")

    pair_direct_checks = 0
    pair_balance_checks = 0
    for left in all_data:
        for right in all_data:
            child = ordinal_data(left, right)
            need(
                child.optional_state
                == parity_convolution(left.optional_state, right.optional_state),
                "optional-state multiplicativity",
            )
            expected_k = 2 * left.optional_state[0] * right.optional_state[0]
            need(child.optional_state == (expected_k, expected_k), "K product law")
            pair_balance_checks += 1
            if len(left.out) + len(right.out) <= 8:
                direct = tournament_data(child.out)
                need(
                    (child.capacities, child.starts, child.ends, child.gate)
                    == (direct.capacities, direct.starts, direct.ends, direct.gate),
                    "direct pair transfer",
                )
                pair_direct_checks += 1

    triple_direct_checks = 0
    cocycle_direct_checks = 0
    for left in all_data:
        for middle in all_data:
            if len(left.out) + len(middle.out) >= 8:
                continue
            left_middle = ordinal_data(left, middle)
            for right in all_data:
                total_order = len(left.out) + len(middle.out) + len(right.out)
                if total_order > 8:
                    continue
                formula = direct_three_block_capacity(left, middle, right)
                direct = tournament_data(ordinal_out(left_middle.out, right.out))
                need(formula == direct.capacities, "direct three-block capacity")
                r_left_middle = remainder(left, middle)
                r_middle_right = remainder(middle, right)
                r_lm_right = remainder(left_middle, right)
                middle_right = ordinal_data(middle, right)
                r_left_mr = remainder(left, middle_right)
                need(
                    r_lm_right + right.hamilton**2 * r_left_middle
                    == r_left_mr + left.hamilton**2 * r_middle_right,
                    "weighted remainder cocycle",
                )
                triple_direct_checks += 1
                cocycle_direct_checks += 1

    # Complete class-presentation census for the local interaction
    # Theta(A,B,C)=R(A▷B,C)-H(A)^2 R(B,C), with C intrinsically no-sink.
    right_remainders = {
        (id(middle), id(right)): remainder(middle, right)
        for middle in all_data
        for right in all_no_sink
    }
    theta_checks = 0
    theta_negative = 0
    theta_zero = 0
    theta_minimum: tuple[Fraction, tuple[object, ...]] | None = None
    theta_digest = hashlib.sha256()
    for left in all_data:
        for middle in all_data:
            left_middle = ordinal_data(left, middle)
            for right in all_no_sink:
                theta = (
                    remainder(left_middle, right)
                    - left.hamilton**2 * right_remainders[(id(middle), id(right))]
                )
                denominator = (
                    left.hamilton**2
                    * middle.hamilton**2
                    * right.hamilton**2
                )
                normalized = Fraction(theta, denominator)
                witness = (
                    len(left.out),
                    len(middle.out),
                    len(right.out),
                    label(left.out),
                    label(middle.out),
                    label(right.out),
                    theta,
                    denominator,
                )
                candidate = (normalized, witness)
                if theta_minimum is None or candidate < theta_minimum:
                    theta_minimum = candidate
                theta_negative += theta < 0
                theta_zero += theta == 0
                theta_digest.update(("|".join(map(str, witness)) + "\n").encode("ascii"))
                theta_checks += 1
    need(theta_checks == 317_680, "no-sink-third-factor presentation count")
    need(theta_negative == 0 and theta_zero == 0, "strict finite Theta census")

    # Attack the overstrong all-third-factor sign through factor order five.
    small = tuple(data for order in range(1, 6) for data in banks[order])
    first_negative: tuple[object, ...] | None = None
    negative_count = 0
    zero_count = 0
    hostile_checks = 0
    for left in small:
        for middle in small:
            left_middle = ordinal_data(left, middle)
            for right in small:
                theta = remainder(left_middle, right) - left.hamilton**2 * remainder(middle, right)
                record = (
                    len(left.out) + len(middle.out) + len(right.out),
                    len(left.out),
                    len(middle.out),
                    len(right.out),
                    label(left.out),
                    label(middle.out),
                    label(right.out),
                    theta,
                )
                if theta < 0:
                    negative_count += 1
                    if first_negative is None or record < first_negative:
                        first_negative = record
                elif theta == 0:
                    zero_count += 1
                hostile_checks += 1
    need(hostile_checks == 20**3, "hostile triple presentation count")
    need(
        first_negative == (5, 3, 1, 1, "101", "", "", -216),
        "minimal overstrong-Theta hostile",
    )

    # Closed-form and convexity controls for P_n and L_n=P_n▷C3.
    formula_checks = 0
    for order in range(1, 65):
        out, capacities = transitive_capacity_data(order)
        need(gate(out, capacities) == transitive_gate_formula(order), "P_n gate formula")
        lollipop_out, lollipop_capacities = lollipop_capacity_data(order)
        need(
            gate(lollipop_out, lollipop_capacities) == lollipop_gate_formula(order),
            "L_n gate formula",
        )
        formula_checks += 2
    need(lollipop_gate_formula(0) == 0, "L_0 boundary")
    increments = [lollipop_gate_formula(n + 1) - lollipop_gate_formula(n) for n in range(64)]
    need(all(value > 0 for value in increments), "positive lollipop increments")
    need(all(increments[i + 1] > increments[i] for i in range(63)),
         "strict lollipop discrete convexity")
    lollipop_remainder_minimum = None
    lollipop_checks = 0
    for left_order in range(1, 33):
        for right_prefix in range(0, 33):
            value = (
                lollipop_gate_formula(left_order + right_prefix)
                - 9 * transitive_gate_formula(left_order)
                - lollipop_gate_formula(right_prefix)
            )
            need(value > 0, "lollipop OS+ remainder")
            record = (value, left_order, right_prefix)
            if lollipop_remainder_minimum is None or record < lollipop_remainder_minimum:
                lollipop_remainder_minimum = record
            lollipop_checks += 1
    need(lollipop_remainder_minimum == (120, 1, 0), "lollipop minimum")

    need(theta_minimum is not None, "Theta minimum exists")
    print("TOURNAMENT_ORDINAL_COCYCLE_PARITY_PRIMARY")
    print("gentourng_sha256", tool_digest.hexdigest())
    print("class_counts", " ".join(f"q{n}={len(banks[n])}" for n in range(1, 7)))
    print("no_sink_counts", " ".join(f"q{n}={len(no_sink[n])}" for n in range(3, 7)))
    print("all_order_parity_class_checks", parity_checks)
    print("terminal_involution_cover_objects", involution_objects)
    print("pair_balance_checks", pair_balance_checks)
    print("pair_direct_transfer_checks", pair_direct_checks)
    print("triple_direct_transfer_checks", triple_direct_checks)
    print("weighted_cocycle_direct_checks", cocycle_direct_checks)
    print("no_sink_third_factor_presentations", theta_checks)
    print("theta_negative", theta_negative, "zero", theta_zero)
    print("theta_normalized_minimum", theta_minimum)
    print("theta_presentation_sha256", theta_digest.hexdigest())
    print("all_third_factor_hostile_checks", hostile_checks)
    print("all_third_factor_negative", negative_count, "zero", zero_count)
    print("first_overstrong_theta_hostile", first_negative)
    print("closed_form_gate_checks", formula_checks)
    print("lollipop_box_checks", lollipop_checks)
    print("lollipop_remainder_minimum", lollipop_remainder_minimum)
    print("PASS")


if __name__ == "__main__":
    main()
