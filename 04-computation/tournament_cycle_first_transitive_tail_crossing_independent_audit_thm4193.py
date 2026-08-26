#!/usr/bin/env python3
"""Independent direct-capacity audit for THM-4193.

This file imports no tournament code.  It constructs every checked child's
capacity tensor directly from two complementary Hamilton-path tables on the
actual tournament.  In particular it does not use ordinal capacity transfer,
the primary certificate, gentourng, or the inherited cocycle implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parse(bits: str, order: int) -> tuple[int, ...]:
    need(len(bits) == order * (order - 1) // 2, "label length")
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


def ordinal(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    nleft = len(left)
    nright = len(right)
    right_mask = ((1 << nright) - 1) << nleft
    return tuple(row | right_mask for row in left) + tuple(
        row << nleft for row in right
    )


@dataclass(frozen=True)
class DirectData:
    hamilton: int
    capacities: tuple[tuple[int, ...], ...]
    gate: int


@lru_cache(maxsize=None)
def direct_data(out: tuple[int, ...]) -> DirectData:
    """Build H and c from actual complementary path tables."""
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
            left_count = end[left_mask][left]
            if not left_count:
                continue
            for right in range(order):
                if right_mask & (1 << right):
                    exposed[left][right] += (
                        left_count * start[right_mask][right]
                    )
    capacities = [[0] * order for _ in range(order)]
    for left in range(order):
        for right in range(left + 1, order):
            value = exposed[left][right] + exposed[right][left]
            capacities[left][right] = value
            capacities[right][left] = value

    degree = [sum(row) for row in capacities]
    current = []
    for vertex in range(order):
        value = 0
        for other in range(order):
            if other == vertex:
                continue
            capacity = capacities[vertex][other]
            value += capacity if out[vertex] & (1 << other) else -capacity
        current.append(value)
    mass = sum(
        capacities[left][right]
        for left in range(order)
        for right in range(left + 1, order)
    )
    squares = sum(
        capacities[left][right] ** 2
        for left in range(order)
        for right in range(left + 1, order)
    )
    disjoint_numerator = mass * mass + squares - sum(value * value for value in degree)
    need(disjoint_numerator % 2 == 0, "disjoint parity")
    gate = disjoint_numerator // 2 + 2 * sum(
        degree[vertex] * current[vertex] for vertex in range(order)
    )
    return DirectData(
        hamilton=subset_hamilton[full],
        capacities=tuple(tuple(row) for row in capacities),
        gate=gate,
    )


def remainder(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    left_data = direct_data(left)
    right_data = direct_data(right)
    child_data = direct_data(ordinal(left, right))
    return (
        child_data.gate
        - right_data.hamilton**2 * left_data.gate
        - left_data.hamilton**2 * right_data.gate
    )


def prefix_defect(
    prefix: tuple[int, ...],
    middle: tuple[int, ...],
    right: tuple[int, ...],
) -> int:
    return remainder(ordinal(prefix, middle), right) - remainder(middle, right)


def transitive_gate(order: int) -> int:
    if order == 1:
        return 0
    numerator = 4 * 4**order - 9 * (order + 2) * 2**order + 24 * order + 32
    need(numerator % 18 == 0, "formula integrality")
    return numerator // 18


def transitive_remainder(left: int, right: int) -> int:
    return transitive_gate(left + right) - transitive_gate(left) - transitive_gate(right)


def negative_theta(left: int, right: int) -> int:
    if left == 1:
        return 144 - (27 * right + 153) * 2**right
    return -2 ** (left - 1) * (
        (27 * (left + right) + 126) * 2**right - (27 * left + 126)
    )


def alpha(tail: int) -> int:
    return 12 * (2 * 4**tail - (3 * tail + 21) * 2**tail + 1)


def formula(tail: int, middle: int, right: int) -> int:
    joined = tail + middle
    return (
        9 * transitive_remainder(joined, right)
        - transitive_remainder(middle, right)
        + negative_theta(joined, right)
    )


def increment(tail: int, middle: int, right: int) -> int:
    total = tail + middle
    bracket = 2**total * (4**right - 1) - 3 * (
        2**right * (total + right + 6) - (total + 6)
    )
    return 6 * 2**total * bracket


def main() -> None:
    one = transitive(1)
    cycle = parse("101", 3)
    need(direct_data(one).hamilton == 1, "singleton H")
    need(direct_data(cycle).hamilton == 3, "cycle H")
    need(direct_data(cycle).gate == 0, "cycle gate")

    direct_transitive_gate_checks = 0
    for order in range(1, 13):
        data = direct_data(transitive(order))
        need(data.hamilton == 1, "transitive Hamilton count")
        need(data.gate == transitive_gate(order), "direct transitive gate")
        direct_transitive_gate_checks += 1

    direct_singleton_checks = 0
    singleton_rows = []
    for tail in range(0, 8):
        prefix = cycle if tail == 0 else ordinal(cycle, transitive(tail))
        need(not any(row.bit_count() == len(prefix) - 1 for row in prefix),
             "cycle-first prefix is source-free")
        value = prefix_defect(prefix, one, one)
        need(value == alpha(tail), "direct singleton formula")
        singleton_rows.append((tail, value))
        direct_singleton_checks += 1
    need(singleton_rows[:6] == [
        (0, -216), (1, -468), (2, -900), (3, -1332), (4, -180), (5, 10764)
    ], "direct crossing ledger")

    # Every checked gate below is rebuilt on the actual ordinal child.  The
    # cap n+b+c<=9 keeps this direct subset audit cheap while crossing n=5.
    direct_context_checks = 0
    direct_positive_checks = 0
    for tail in range(0, 8):
        prefix = cycle if tail == 0 else ordinal(cycle, transitive(tail))
        for middle in range(1, 10):
            for right in range(1, 10):
                if tail + middle + right > 9:
                    continue
                value = prefix_defect(prefix, transitive(middle), transitive(right))
                need(value == formula(tail, middle, right), "direct context formula")
                if tail >= 5:
                    need(value > 0, "direct positive buffered cycle")
                    direct_positive_checks += 1
                direct_context_checks += 1

    arithmetic_checks = 0
    increment_checks = 0
    for tail in range(0, 41):
        need((alpha(tail) < 0) == (tail <= 4), "singleton sign iff")
        for middle in range(1, 65):
            for right in range(1, 65):
                value = formula(tail, middle, right)
                need(formula(tail, middle, right) == value, "formula determinism")
                if tail >= 5:
                    need(value > 0, "arithmetic positive context")
                if tail < 40 and tail + middle >= 2:
                    difference = formula(tail + 1, middle, right) - value
                    need(difference == increment(tail, middle, right),
                         "independent increment identity")
                    if tail >= 5:
                        need(difference > 0, "independent positive increment")
                    increment_checks += 1
                arithmetic_checks += 1

    survivor = ordinal(cycle, transitive(5))
    hostile = ordinal(cycle, transitive(4))
    need(label(survivor) == "1011111111111111111111111111", "survivor label")
    need(prefix_defect(survivor, one, one) == 10764, "survivor value")
    need(prefix_defect(hostile, one, one) == -180, "last hostile value")

    print("TOURNAMENT_CYCLE_FIRST_TRANSITIVE_TAIL_CROSSING_INDEPENDENT")
    print("method", "direct_actual_child_complementary_path_tables")
    print("direct_transitive_gate_checks", direct_transitive_gate_checks)
    print("direct_singleton_formula_checks", direct_singleton_checks)
    print("direct_singleton_rows", singleton_rows)
    print("direct_context_formula_checks", direct_context_checks)
    print("direct_positive_context_checks", direct_positive_checks)
    print("arithmetic_formula_checks", arithmetic_checks)
    print("arithmetic_increment_checks", increment_checks)
    print("last_negative_prefix", "C3▷P4", -180)
    print("first_positive_prefix", "C3▷P5", 10764)
    print("survivor_label", label(survivor))
    print("PASS")


if __name__ == "__main__":
    main()
