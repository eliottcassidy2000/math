#!/usr/bin/env python3
"""Exact finite audit for the Rule 30 rectangular macroblock compiler.

The universal statements audited here are elementary locality and
left-permutivity lemmas.  The finite checks deliberately use independent
list, centered-word, and right-edge packed implementations.  They do not
imply a lower bound for the distinguished one-seed center or solve any Rule
30 prize problem.
"""

from __future__ import annotations


DIRECT_HORIZON = 192
LOCAL_WORD_WIDTH_CAP = 10
TABLE_CASES = ((1, 1), (2, 1), (3, 2), (4, 2), (4, 3), (5, 2))
TRIANGULAR_CASES = ((1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (4, 2))
COCYCLE_CASES = ((1, 1, 1), (2, 1, 2), (3, 2, 1), (4, 2, 2))
MACRO_CASES = ((2, 1), (3, 2), (4, 2), (5, 3))
RECTANGULAR_COVERAGE_EXPECTED = (
    (1, 2, 4, 11),
    (2, 4, 8, 152),
    (3, 6, 12, 1161),
    (4, 8, 16, 6144),
)
TERNARY_BLOCK_COVERAGE_EXPECTED = (
    (1, 1, 3, 3),
    (2, 2, 6, 52),
    (3, 3, 9, 210),
    (4, 4, 12, 980),
    (5, 5, 15, 3740),
    (6, 6, 18, 14718),
)


def check(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def rule30(left: int, center: int, right: int) -> int:
    return left ^ center ^ right ^ (center & right)


def list_step(row: list[int]) -> list[int]:
    padded = [0, *row, 0]
    return [
        rule30(padded[index], padded[index + 1], padded[index + 2])
        for index in range(len(row))
    ]


def word_step(row: int, width: int) -> int:
    mask = (1 << width) - 1
    left = (row << 1) & mask
    center = row
    right = row >> 1
    return (left ^ center ^ right ^ (center & right)) & mask


def word_to_list(word: int, width: int) -> list[int]:
    return [(word >> index) & 1 for index in range(width)]


def bits_to_word(bits: list[int]) -> int:
    return sum(bit << index for index, bit in enumerate(bits))


def audit_local_implementations() -> None:
    for width in range(1, LOCAL_WORD_WIDTH_CAP + 1):
        for word in range(1 << width):
            observed = word_step(word, width)
            expected = bits_to_word(list_step(word_to_list(word, width)))
            check(observed == expected, f"local width={width} word={word}")


def block_map_list(word: int, output_width: int, height: int) -> int:
    input_width = output_width + 2 * height
    row = word_to_list(word, input_width)
    for _ in range(height):
        row = list_step(row)
    return bits_to_word(row[height : height + output_width])


def block_map_word(word: int, output_width: int, height: int) -> int:
    input_width = output_width + 2 * height
    row = word
    for _ in range(height):
        row = word_step(row, input_width)
    return (row >> height) & ((1 << output_width) - 1)


def block_map_with_trace(word: int, output_width: int, height: int) -> tuple[int, int]:
    """Return the final middle block and its first site's vertical trace."""

    input_width = output_width + 2 * height
    row = word
    trace = 0
    for time in range(height):
        row = word_step(row, input_width)
        trace |= ((row >> height) & 1) << time
    block = (row >> height) & ((1 << output_width) - 1)
    return block, trace


def build_table(output_width: int, height: int) -> list[int]:
    input_width = output_width + 2 * height
    return [
        block_map_word(word, output_width, height)
        for word in range(1 << input_width)
    ]


def build_trace_table(output_width: int, height: int) -> list[tuple[int, int]]:
    input_width = output_width + 2 * height
    return [
        block_map_with_trace(word, output_width, height)
        for word in range(1 << input_width)
    ]


def audit_block_tables() -> None:
    for output_width, height in TABLE_CASES:
        input_width = output_width + 2 * height
        for word in range(1 << input_width):
            check(
                block_map_word(word, output_width, height)
                == block_map_list(word, output_width, height),
                f"block s={output_width} h={height} word={word}",
            )


def audit_triangular_fibres() -> None:
    for output_width, height in TRIANGULAR_CASES:
        output_mask = (1 << output_width) - 1
        for tail in range(1 << (2 * height)):
            image = {
                block_map_word(prefix | (tail << output_width), output_width, height)
                for prefix in range(1 << output_width)
            }
            check(
                image == set(range(1 << output_width)),
                f"triangular fibre s={output_width} h={height} tail={tail}",
            )

            # The r-th output has Boolean derivative one in owner bit r and
            # derivative zero in every earlier owner bit.
            for prefix in range(1 << output_width):
                word = prefix | (tail << output_width)
                output = block_map_word(word, output_width, height)
                for row in range(output_width):
                    toggled = block_map_word(
                        word ^ (1 << row), output_width, height
                    )
                    check(
                        ((output ^ toggled) >> row) & 1 == 1,
                        f"unit diagonal s={output_width} h={height} r={row}",
                    )
                    for earlier in range(row):
                        toggled_earlier = block_map_word(
                            word ^ (1 << earlier), output_width, height
                        )
                        check(
                            ((output ^ toggled_earlier) >> row) & 1 == 0,
                            (
                                "zero below diagonal "
                                f"s={output_width} h={height} r={row} q={earlier}"
                            ),
                        )


def audit_width_changing_cocycle() -> None:
    for output_width, first_height, second_height in COCYCLE_CASES:
        input_width = output_width + 2 * (first_height + second_height)
        for word in range(1 << input_width):
            direct = block_map_word(
                word, output_width, first_height + second_height
            )
            intermediate = block_map_word(
                word, output_width + 2 * second_height, first_height
            )
            composed = block_map_word(
                intermediate, output_width, second_height
            )
            check(
                direct == composed,
                (
                    f"cocycle s={output_width} h={first_height} "
                    f"k={second_height} word={word}"
                ),
            )


def extract_window(row: int, start: int, width: int) -> int:
    mask = (1 << width) - 1
    if start >= 0:
        return (row >> start) & mask
    return (row << (-start)) & mask


def macro_prefix(horizon: int, output_width: int, height: int) -> bytearray:
    """Compute the one-seed center prefix using only charged block tables."""

    table = build_trace_table(output_width, height)
    block_mask = (1 << output_width) - 1
    input_width = output_width + 2 * height
    input_mask = (1 << input_width) - 1

    side_blocks = (horizon + height + output_width - 1) // output_width + 3
    row_width = (2 * side_blocks + 1) * output_width
    center_index = side_blocks * output_width
    row_mask = (1 << row_width) - 1
    row = 1 << center_index
    prefix = bytearray([1])

    macrosteps, remainder = divmod(horizon, height)
    for _ in range(macrosteps):
        following = 0
        center_trace = None
        for block_start in range(0, row_width, output_width):
            address = extract_window(row, block_start - height, input_width)
            address &= input_mask
            block, trace = table[address]
            following |= (block & block_mask) << block_start
            if block_start == center_index:
                center_trace = trace
        check(center_trace is not None, "central macroblock is aligned")
        prefix.extend((center_trace >> time) & 1 for time in range(height))
        row = following & row_mask

    for _ in range(remainder):
        row = word_step(row, row_width)
        prefix.append((row >> center_index) & 1)

    check(len(prefix) == horizon + 1, "macro prefix length")
    return prefix


def direct_center_prefix(horizon: int) -> bytearray:
    """Independent centered-set implementation on a dynamically grown cone."""

    row = frozenset({0})
    prefix = bytearray([1])
    for time in range(horizon):
        following = {
            site
            for site in range(-time - 1, time + 2)
            if rule30(
                int(site - 1 in row),
                int(site in row),
                int(site + 1 in row),
            )
        }
        row = frozenset(following)
        prefix.append(int(0 in row))
    return prefix


def audit_macro_compiler() -> None:
    direct = direct_center_prefix(DIRECT_HORIZON)
    for output_width, height in MACRO_CASES:
        observed = macro_prefix(DIRECT_HORIZON, output_width, height)
        check(
            observed == direct,
            f"macro prefix s={output_width} h={height} n={DIRECT_HORIZON}",
        )


def ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def reverse_low_bits(word: int, width: int) -> int:
    check(0 <= width <= 32, "32-bit reversal width")
    word &= (1 << width) - 1
    word = ((word & 0x55555555) << 1) | ((word >> 1) & 0x55555555)
    word = ((word & 0x33333333) << 2) | ((word >> 2) & 0x33333333)
    word = ((word & 0x0F0F0F0F) << 4) | ((word >> 4) & 0x0F0F0F0F)
    word = ((word & 0x00FF00FF) << 8) | ((word >> 8) & 0x00FF00FF)
    word = ((word & 0x0000FFFF) << 16) | ((word >> 16) & 0x0000FFFF)
    return word >> (32 - width) if width else 0


def edge_window(row: int, time: int, start: int, width: int) -> int:
    """Decode an ascending spatial window from right-edge packed coordinates."""

    low_edge_index = time - (start + width - 1)
    mask = (1 << width) - 1
    if low_edge_index >= 0:
        segment = (row >> low_edge_index) & mask
    else:
        segment = (row << (-low_edge_index)) & mask
    return reverse_low_bits(segment, width)


def aligned_coverage_centered(
    output_width: int, height: int, horizon: int
) -> tuple[int | None, int]:
    window_width = output_width + 2 * height
    universe_size = 1 << window_width
    seen = bytearray(universe_size)
    seen[0] = 1
    distinct = 1

    offset = horizon + window_width + 2
    total_width = 2 * offset + 1
    row_mask = (1 << total_width) - 1
    row = 1 << offset

    for time in range(0, horizon + 1, height):
        first_block = ceil_div(
            -time - (output_width + height - 1), output_width
        )
        last_block = (time + height) // output_width
        for block in range(first_block, last_block + 1):
            start = block * output_width - height + offset
            word = extract_window(row, start, window_width)
            if not seen[word]:
                seen[word] = 1
                distinct += 1
        if distinct == universe_size:
            return time, distinct
        for _ in range(height):
            row = word_step(row, total_width) & row_mask
    return None, distinct


def aligned_coverage_right_edge(
    output_width: int, height: int, horizon: int
) -> tuple[int | None, int]:
    window_width = output_width + 2 * height
    universe_size = 1 << window_width
    seen = bytearray(universe_size)
    seen[0] = 1
    distinct = 1
    row = 1

    for time in range(horizon + 1):
        if time % height == 0:
            first_block = ceil_div(
                -time - (output_width + height - 1), output_width
            )
            last_block = (time + height) // output_width
            for block in range(first_block, last_block + 1):
                start = block * output_width - height
                word = edge_window(row, time, start, window_width)
                if not seen[word]:
                    seen[word] = 1
                    distinct += 1
            if distinct == universe_size:
                return time, distinct
        row ^= (row << 1) | (row << 2)
    return None, distinct


def audit_fixed_seed_address_saturation(
    expected_rows: tuple[tuple[int, int, int, int], ...]
) -> tuple[tuple[int, int, int, int], ...]:
    observed_rows = []
    for height, output_width, window_width, expected_time in expected_rows:
        check(
            window_width == output_width + 2 * height,
            f"coverage width invoice h={height} s={output_width}",
        )
        centered = aligned_coverage_centered(
            output_width, height, expected_time
        )
        right_edge = aligned_coverage_right_edge(
            output_width, height, expected_time
        )
        expected = (expected_time, 1 << window_width)
        check(centered == expected, f"centered coverage h={height}")
        check(right_edge == expected, f"right-edge coverage h={height}")
        observed_rows.append((height, output_width, window_width, expected_time))
    return tuple(observed_rows)


def main() -> None:
    audit_local_implementations()
    audit_block_tables()
    audit_triangular_fibres()
    audit_width_changing_cocycle()
    audit_macro_compiler()
    rectangular_coverage = audit_fixed_seed_address_saturation(
        RECTANGULAR_COVERAGE_EXPECTED
    )
    ternary_coverage = audit_fixed_seed_address_saturation(
        TERNARY_BLOCK_COVERAGE_EXPECTED
    )

    print("RULE 30 RECTANGULAR MACROBLOCK EXACT AUDIT")
    print(f"local_list_word_widths=1..{LOCAL_WORD_WIDTH_CAP}")
    print(f"block_table_cases_s_h={list(TABLE_CASES)}")
    print(f"triangular_fibre_cases_s_h={list(TRIANGULAR_CASES)}")
    print(f"width_changing_cocycle_cases_s_h_k={list(COCYCLE_CASES)}")
    print(
        "macro_center_prefix_cases_s_h="
        f"{list(MACRO_CASES)} n=0..{DIRECT_HORIZON}"
    )
    print(
        "aligned_rectangular_cumulative_full_address_coverage_h_s_L_t="
        f"{list(rectangular_coverage)}"
    )
    print(
        "aligned_three_block_cumulative_full_address_coverage_h_s_L_t="
        f"{list(ternary_coverage)}"
    )
    print(
        "scope=FINITE-EXACT audit of all-n locality/permutivity proofs; "
        "no fixed-seed lower bound and no Rule 30 prize conclusion"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
