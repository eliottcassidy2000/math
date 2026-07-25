#!/usr/bin/env python3
"""Exact finite referee for THM-2253's online dyadic-contrast extractor.

The mathematical proof is an involution on terminal cylinders.  This
companion exhausts binary and ternary controls, checks that the first
contrast node is unique and is preserved by swapping its homogeneous
children, verifies output reversal and fixed-composition bisection, and
tests the critical-run deadline on every hostile continuation in the
declared finite universes.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation_two(n: int) -> int:
    require(n >= 1, "the two-adic valuation requires n>=1")
    exponent = 0
    while n % 2 == 0:
        exponent += 1
        n //= 2
    return exponent


def contrast_nodes_ending(
    word: tuple[int, ...], time: int
) -> list[tuple[int, int, int, int]]:
    """Return (start, half_length, left_symbol, right_symbol) at `time`."""
    nodes: list[tuple[int, int, int, int]] = []
    half = 1
    while 2 * half <= time:
        if time % (2 * half) == 0:
            start = time - 2 * half
            left = word[start : start + half]
            right = word[start + half : time]
            if (
                len(left) == half
                and len(right) == half
                and len(set(left)) == 1
                and len(set(right)) == 1
                and left[0] != right[0]
            ):
                nodes.append((start, half, left[0], right[0]))
        half *= 2
    return nodes


def first_contrast(
    word: tuple[int, ...],
) -> tuple[int, int, int, int, int] | None:
    """Return (time,start,half,left_symbol,right_symbol), if already seen."""
    for time in range(2, len(word) + 1):
        nodes = contrast_nodes_ending(word, time)
        require(
            len(nodes) <= 1,
            f"two contrast nodes end simultaneously: {word=} {time=} {nodes=}",
        )
        if nodes:
            start, half, left, right = nodes[0]
            return time, start, half, left, right
    return None


def swap_selected(
    word: tuple[int, ...],
    node: tuple[int, int, int, int, int],
) -> tuple[int, ...]:
    time, start, half, _left, _right = node
    require(time == len(word), "only a terminal prefix may be swapped")
    return (
        word[:start]
        + word[start + half : time]
        + word[start : start + half]
    )


def binary_orientation(left: int, right: int) -> int:
    require(left != right and left in (0, 1) and right in (0, 1), "bad bit edge")
    return 1 if (left, right) == (0, 1) else -1


TERNARY_EDGES = {(0, 1), (1, 2), (2, 0)}


def ternary_orientation(left: int, right: int) -> int:
    require(left != right and left in range(3) and right in range(3), "bad ternary edge")
    return 1 if (left, right) in TERNARY_EDGES else -1


def counts(word: tuple[int, ...], alphabet_size: int) -> tuple[int, ...]:
    return tuple(word.count(symbol) for symbol in range(alphabet_size))


def audit_terminal_universe(
    alphabet_size: int,
    maximum_length: int,
) -> tuple[int, int]:
    orientation = binary_orientation if alphabet_size == 2 else ternary_orientation
    terminal_count = 0
    shell_balance: dict[tuple[int, tuple[int, ...]], list[int]] = defaultdict(
        lambda: [0, 0]
    )
    for time in range(2, maximum_length + 1):
        for word in product(range(alphabet_size), repeat=time):
            node = first_contrast(word)
            if node is None or node[0] != time:
                continue
            terminal_count += 1
            swapped = swap_selected(word, node)
            swapped_node = first_contrast(swapped)
            require(
                swapped_node is not None
                and swapped_node[:3] == node[:3]
                and swapped_node[3:] == (node[4], node[3]),
                f"terminal node not preserved by swap: {word=} {node=} {swapped_node=}",
            )
            require(
                swap_selected(swapped, swapped_node) == word,
                f"terminal swap is not involutive: {word=}",
            )
            sign = orientation(node[3], node[4])
            swapped_sign = orientation(swapped_node[3], swapped_node[4])
            require(sign == -swapped_sign, f"orientation did not reverse: {word=}")
            require(
                counts(word, alphabet_size) == counts(swapped, alphabet_size),
                f"composition changed: {word=} {swapped=}",
            )
            bucket = shell_balance[(time, counts(word, alphabet_size))]
            bucket[0 if sign == 1 else 1] += 1
    require(
        all(heads == tails for heads, tails in shell_balance.values()),
        f"terminal composition shell failed to bisect for alphabet {alphabet_size}",
    )
    return terminal_count, len(shell_balance)


def audit_deadline(
    alphabet_size: int,
    maximum_n: int,
) -> tuple[int, int]:
    tested = 0
    equality_witnesses = 0
    for n in range(1, maximum_n + 1):
        half = 1 << valuation_two(n)
        deadline = n + half
        for first in range(alphabet_size):
            for changed in range(alphabet_size):
                if first == changed:
                    continue
                for suffix in product(
                    range(alphabet_size), repeat=half - 1
                ):
                    word = (first,) * n + (changed,) + suffix
                    node = first_contrast(word)
                    require(
                        node is not None and node[0] <= deadline,
                        f"deadline failure: {alphabet_size=} {n=} {word=} {node=}",
                    )
                    tested += 1
                constant_right = (first,) * n + (changed,) * half
                node = first_contrast(constant_right)
                require(
                    node is not None and node[0] == deadline,
                    f"claimed sharp continuation is not sharp: {n=} {node=}",
                )
                equality_witnesses += 1
        if n & (n - 1):
            require(
                3 * half <= n,
                f"nonpower critical run violates half<=n/3: {n=} {half=}",
            )
            require(
                3 * deadline <= 4 * n,
                f"nonpower 4n/3 deadline failure: {n=} {deadline=}",
            )
    return tested, equality_witnesses


def main() -> None:
    binary_terminals, binary_shells = audit_terminal_universe(2, 16)
    ternary_terminals, ternary_shells = audit_terminal_universe(3, 10)
    binary_deadlines, binary_equalities = audit_deadline(2, 16)
    ternary_deadlines, ternary_equalities = audit_deadline(3, 8)
    print(
        "binary_terminal_audit=PASS "
        f"max_length=16 terminals={binary_terminals} shells={binary_shells}"
    )
    print(
        "ternary_tournament_terminal_audit=PASS "
        f"max_length=10 terminals={ternary_terminals} shells={ternary_shells}"
    )
    print(
        "binary_deadline_audit=PASS "
        f"max_n=16 hostile_continuations={binary_deadlines} "
        f"sharp_constant_right={binary_equalities}"
    )
    print(
        "ternary_deadline_audit=PASS "
        f"max_n=8 hostile_continuations={ternary_deadlines} "
        f"sharp_constant_right={ternary_equalities}"
    )
    print("deadline_formula=tau<=n+2^nu2(n)")
    print("nonpower_corollary=tau<=4n/3")
    print("status=THM-2253_PROVED_VERIFIED_EXACT")


if __name__ == "__main__":
    main()
