#!/usr/bin/env python3
"""Exact arithmetic referee for THM-1287.

The paper provides the selected-prefix survivor/component argument and the
endpoint placement of the globally unique terminal owner.  This referee
checks every recurrence constant, strict tooth-count ceiling, endpoint-word
alternation consequence, private-count conversion, and THM-1275 nested
ceiling used by the finite selected-word residual.
"""

from __future__ import annotations

import ast
import hashlib
from fractions import Fraction
from itertools import product
from pathlib import Path


F = Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ceil_fraction(value: F) -> int:
    return -(-value.numerator // value.denominator)


def strict_tooth_cap(ratio_upper: F) -> int:
    """Cap ceil((6x+1)/7) under the strict hypothesis x<ratio_upper."""
    return ceil_fraction((6 * ratio_upper + 1) / 7)


PREFIX_COEFFICIENTS = {
    1: F(91, 29),
    2: F(91, 22),
    3: F(49, 15),
    4: F(21, 8),
    5: F(7),
}


EXPECTED_EFFECTIVE_CAPS = {
    2: (2, 1, 15, 54, 165, 237),
    3: (2, 7, 1, 31, 95, 136),
    4: (2, 7, 7, 1, 41, 58),
    5: (2, 7, 7, 7, 1, 24),
    6: (2, 7, 7, 7, 7, 1),
}


EXPECTED_RAW_FASTEST_CAPS = {2: 1429, 3: 823, 4: 355, 5: 151, 6: 1}
EXPECTED_PREFIX_H_BOUNDS = {
    2: F(1666),
    3: F(959),
    4: F(413),
    5: F(175),
    6: F(7),
}
EXPECTED_WORD_CAPS = {2: 474, 3: 272, 4: 116, 5: 48, 6: 31}
EXPECTED_FINAL_H_BOUNDS = {2: 1659, 3: 952, 4: 406, 5: 168, 6: 7}


def selected_prefix_table_audit() -> tuple[
    tuple[int, tuple[int, ...], int, int, int, int], ...
]:
    rows: list[tuple[int, tuple[int, ...], int, int, int, int]] = []
    for terminal_rank in range(2, 7):
        selected_caps = [0] * 7
        selected_caps[1] = 2

        # d_i<d_r<=7c before the terminal rank.  The terminal selected count
        # itself is one, while every earlier selected count is at most the
        # corresponding number (at most seven) of teeth meeting G.
        for rank in range(2, terminal_rank):
            selected_caps[rank] = 7
        selected_caps[terminal_rank] = 1

        prefix_h_bound = F(7)
        for prefix_rank in range(terminal_rank, 6):
            next_rank = prefix_rank + 1
            component_bound = PREFIX_COEFFICIENTS[prefix_rank] * (
                1 + sum(selected_caps[1 : prefix_rank + 1])
            )
            if next_rank == 3:
                component_bound = min(component_bound, F(84, 5))
            selected_caps[next_rank] = strict_tooth_cap(component_bound)
            if next_rank == 6:
                prefix_h_bound = component_bound

        raw_fastest_cap = selected_caps[6]
        if terminal_rank < 6:
            low_total = sum(selected_caps[1:6])
            effective_fastest_cap = low_total
            word_cap = 2 * low_total
            final_h_bound = 7 * low_total
            selected_caps[6] = effective_fastest_cap
        else:
            low_total = sum(selected_caps[1:6])
            effective_fastest_cap = 1
            word_cap = low_total + 1
            final_h_bound = 7

        effective_caps = tuple(selected_caps[1:7])
        require(
            effective_caps == EXPECTED_EFFECTIVE_CAPS[terminal_rank],
            f"selected cap row changed at terminal rank {terminal_rank}",
        )
        require(
            raw_fastest_cap == EXPECTED_RAW_FASTEST_CAPS[terminal_rank],
            f"raw fastest cap changed at terminal rank {terminal_rank}",
        )
        require(
            prefix_h_bound == EXPECTED_PREFIX_H_BOUNDS[terminal_rank],
            f"prefix h bound changed at terminal rank {terminal_rank}",
        )
        require(
            word_cap == EXPECTED_WORD_CAPS[terminal_rank],
            f"word cap changed at terminal rank {terminal_rank}",
        )
        require(
            final_h_bound == EXPECTED_FINAL_H_BOUNDS[terminal_rank],
            f"final h/c bound changed at terminal rank {terminal_rank}",
        )
        require(
            effective_fastest_cap <= raw_fastest_cap,
            "endpoint alternation did not improve the raw fastest cap",
        )
        rows.append(
            (
                terminal_rank,
                effective_caps,
                raw_fastest_cap,
                word_cap,
                int(prefix_h_bound),
                final_h_bound,
            )
        )
    return tuple(rows)


def endpoint_alternation_audit() -> tuple[int, int, int]:
    """Exhaust binary owner shadows in both endpoint orientations."""
    rows = 0
    left_rows = 0
    right_rows = 0
    equality_rows = 0
    for length in range(1, 18):
        for word in product((0, 1), repeat=length):
            # 1 means fastest; adjacent selected occurrences cannot both have
            # owner h.  Low-low adjacency is allowed in this binary quotient.
            if any(word[index] == word[index + 1] == 1 for index in range(length - 1)):
                continue
            for endpoint_side in (-1, 1):
                terminal = word[0] if endpoint_side < 0 else word[-1]
                if terminal != 0:
                    continue
                fastest = sum(word)
                low = length - fastest
                require(fastest <= low, "endpoint alternation injection failed")
                equality_rows += fastest == low
                left_rows += endpoint_side < 0
                right_rows += endpoint_side > 0
                rows += 1
    require(left_rows == right_rows > 0, "endpoint reflection ledger changed")
    require(equality_rows > 0, "alternation cap was never sharp")
    return rows, left_rows, equality_rows


def private_count_conversion_audit() -> int:
    rows = 0
    for terminal_rank in range(2, 7):
        fastest_cap = EXPECTED_EFFECTIVE_CAPS[terminal_rank][5]
        final_bound = EXPECTED_FINAL_H_BOUNDS[terminal_rank]
        for carrier in range(1, 51):
            boundary = final_bound * carrier
            for fastest_speed in range(max(1, boundary - 3), boundary + 4):
                private_floor = ceil_fraction(F(fastest_speed, 7 * carrier))
                if private_floor <= fastest_cap:
                    require(
                        fastest_speed <= final_bound * carrier,
                        "private count did not imply the h/c cap",
                    )
                if fastest_speed > final_bound * carrier:
                    require(
                        private_floor > fastest_cap,
                        "speed above final cap survived the private count",
                    )
                rows += 1
    return rows


def tail_exception_count_audit() -> tuple[int, int]:
    rows = 0
    positive_worst_rows = 0
    for carrier in range(1, 31):
        for fastest_speed in range(carrier + 1, 180 * carrier + 1):
            private_floor = ceil_fraction(F(fastest_speed, 7 * carrier))
            for eligible in range(6):
                exact_from_private = ceil_fraction(F(private_floor, eligible + 1)) - 1
                direct = ceil_fraction(
                    F(fastest_speed, 7 * carrier * (eligible + 1))
                ) - 1
                require(
                    exact_from_private == direct,
                    "nested THM-1275 ceiling identity failed",
                )
                worst = ceil_fraction(F(fastest_speed, 42 * carrier)) - 1
                require(direct >= worst, "e<=5 worst-case tail bound failed")
                positive_worst_rows += worst > 0
                rows += 1
    require(positive_worst_rows > 0, "42c tail bound never became positive")
    return rows, positive_worst_rows


def source_has_no_assert_nodes() -> int:
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "referee contains optimization-sensitive assert nodes")
    return count


def main() -> None:
    no_asserts = source_has_no_assert_nodes()
    cap_rows = selected_prefix_table_audit()
    alternation_rows = endpoint_alternation_audit()
    private_rows = private_count_conversion_audit()
    tail_rows = tail_exception_count_audit()

    print("THM-1287 SELECTED-PREFIX TERMINAL-WORD EXACT AUDIT")
    print(f"Python assert nodes = {no_asserts}")
    for row in cap_rows:
        print(
            "rank/caps/rawK/word/prefixH/finalH = "
            f"{row[0]}/{row[1]}/{row[2]}/{row[3]}/{row[4]}/{row[5]}"
        )
    print(
        "endpoint-alternation rows/one-side/equality = "
        f"{alternation_rows}"
    )
    print(f"private-count boundary rows = {private_rows}")
    print(f"tail nested-ceiling rows/positive-worst = {tail_rows}")
    print("selected_prefix=survivor_mass_preserved; components<=1+sum(selected_n_i)")
    print("endpoint_injection=K<=L; word=K+L<=2L")
    print("private_count=ceil(h/(7c))<=K")
    print("tail_count=X+T>=ceil(K/(e+1))-1>=ceil(h/(42c))-1")
    print("finite=selected_owner_word_only; carrier_and_phase_stalk_unbounded")
    print("STATUS=PASS")
    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(f"source_sha256={digest}")


if __name__ == "__main__":
    main()
