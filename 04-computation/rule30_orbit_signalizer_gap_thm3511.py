#!/usr/bin/env python3
"""Finite exact companion for proved THM-3511."""

from __future__ import annotations

import hashlib


A, B, C = 0, 1, 2
LETTERS = "ABC"
RENORMALIZATION_MAX_M = 22
DIRECT_SECTION_MAX_M = 16
PHASE_OWNER_MAX_M = 7
TAIL_BITS = 40
PACKED_WIDTH = 96

EXPECTED_V_0_TO_23 = (
    1,
    3,
    4,
    6,
    7,
    9,
    15,
    16,
    24,
    25,
    27,
    29,
    34,
    36,
    37,
    39,
    41,
    43,
    48,
    49,
    51,
    54,
    55,
    58,
)
EXPECTED_TRANSCRIPT_SHA256 = "d06c6d7291bf043d028760f0dc1af50a225aa103a50a14084715a0ba9672e5ca"
EXPECTED_DIRECT_SHA256 = "e6f413b60231d610f54308302bab7ec3bd87b506930a0d17e9f5223beb049d96"
EXPECTED_S7_WORD_SHA256 = "c72aefdab8539df7d930ac9c4d503ccfd036757fb0f634ea74a261aa89483be9"
EXPECTED_S17_WORD_SHA256 = "0d92bb368f0091b40feb4afa8d9b9bdfca3ac82aca8d394776173d263a529540"
EXPECTED_S7_PERM4 = "7e5cb290f6d43a18"
EXPECTED_S17_PERM4 = "f6543a987edcb210"


def check(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"check failed: {label}")


def nu2(x: int) -> int:
    check(x != 0, "nu2 zero")
    return (x & -x).bit_length() - 1


def phi(x: int, width: int) -> int:
    return (x ^ ((x << 1) | (x << 2))) & ((1 << width) - 1)


def activity(word: bytes | bytearray) -> int:
    return sum(letter != A for letter in word) & 1


def tau(word: bytes | bytearray, start: int) -> tuple[bytearray, int]:
    """Right-action root section; letters are applied from left to right."""
    out = bytearray(len(word))
    control = start
    for index, letter in enumerate(word):
        if control:
            out[index] = B
        else:
            out[index] = C if letter == B else A
        if letter != A:
            control ^= 1
    return out, control


def section_at_prefix(word: bytes | bytearray, prefix: int, length: int) -> bytearray:
    section = bytearray(word)
    for depth in range(length):
        section, _ = tau(section, (prefix >> depth) & 1)
    return section


def thm3463_display_section(word: bytes | bytearray, start: int) -> bytearray:
    """Section in THM-3463's opposite factor-display convention."""
    incoming = [0] * len(word)
    control = start
    for index in range(len(word) - 1, -1, -1):
        incoming[index] = control
        if word[index] != A:
            control ^= 1
    out = bytearray(len(word))
    for index, (letter, bit) in enumerate(zip(word, incoming)):
        if bit:
            out[index] = B
        else:
            out[index] = C if letter == B else A
    return out


def renormalize(word: bytes | bytearray) -> tuple[bytearray, int]:
    check(activity(word) == 1, ("inactive signalizer", len(word)))
    left, left_end = tau(word, 0)
    right, right_end = tau(word, 1)
    check(left_end == 1 and right_end == 0, "active child routing")
    current = left + right
    gap = 1
    while activity(current) == 0:
        current, end = tau(current, 0)
        check(end == 0, ("fixed zero descent", gap))
        gap += 1
        check(gap <= 64, "finite probe overflow")
    return current, gap


def apply_generator(letter: int, value: int, width: int) -> int:
    """Apply one Mealy state to a finite low-to-high binary word."""
    state = letter
    out = 0
    for depth in range(width):
        bit = (value >> depth) & 1
        if state == A:
            output = bit
            state = A if bit == 0 else B
        elif state == B:
            output = 1 - bit
            state = C if bit == 0 else B
        else:
            output = 1 - bit
            state = A if bit == 0 else B
        out |= output << depth
    return out


def apply_word(word: bytes | bytearray, value: int, width: int) -> int:
    out = value
    for letter in word:
        out = apply_generator(letter, out, width)
    return out


def borrow_decode(x: int, y: int, bits: int) -> int:
    out = 0
    borrow = 0
    for depth in range(bits):
        xb = (x >> depth) & 1
        yb = (y >> depth) & 1
        xor_bit = xb ^ yb
        out |= (xor_bit ^ borrow) << depth
        borrow = ((1 - yb) & xb) | (borrow & (1 - xor_bit))
    return out


def generator_actions(width: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(apply_generator(letter, value, width) for value in range(1 << width))
        for letter in (A, B, C)
    )


def permutation(
    word: bytes | bytearray,
    actions: tuple[tuple[int, ...], ...],
    block_cache: dict[bytes, tuple[int, ...]],
    block_size: int = 16,
) -> tuple[int, ...]:
    size = len(actions[0])
    result = list(range(size))
    for offset in range(0, len(word), block_size):
        block = bytes(word[offset : offset + block_size])
        block_permutation = block_cache.get(block)
        if block_permutation is None:
            values = list(range(size))
            for letter in block:
                values = [actions[letter][value] for value in values]
            block_permutation = tuple(values)
            block_cache[block] = block_permutation
        result = [block_permutation[value] for value in result]
    return tuple(result)


def packed_power_rows() -> dict[int, int]:
    requested = {1 << m for m in range(DIRECT_SECTION_MAX_M + 2)}
    rows: dict[int, int] = {0: 1}
    row = 1
    for time in range(1, (1 << (DIRECT_SECTION_MAX_M + 1)) + 1):
        row = phi(row, PACKED_WIDTH)
        if time in requested:
            rows[time] = row
    check(requested <= rows.keys(), "packed power samples")
    return rows


def signalizer_words() -> tuple[list[bytearray], list[int], str, str]:
    actions4 = generator_actions(4)
    block_cache: dict[bytes, tuple[int, ...]] = {}
    words: list[bytearray] = []
    gaps: list[int] = []
    transcript_lines: list[str] = []
    direct_lines: list[str] = []
    word = bytearray((B,))

    for m in range(RENORMALIZATION_MAX_M + 1):
        v = EXPECTED_V_0_TO_23[m]
        next_word, gap = renormalize(word)
        check(gap == EXPECTED_V_0_TO_23[m + 1] - v, ("gap", m, gap))
        check(len(word) == 1 << m, ("word length", m))
        check(activity(word) == 1, ("word activity", m))

        perm4 = permutation(word, actions4, block_cache)
        perm4_text = "".join(format(value, "x") for value in perm4)
        word_digest = hashlib.sha256(word).hexdigest()
        transcript_lines.append(
            f"m={m:02d} v={v:02d} d={gap} len={len(word)} "
            f"perm4={perm4_text} word_sha256={word_digest}"
        )

        if m <= DIRECT_SECTION_MAX_M:
            direct = bytearray((B,)) * (1 << m)
            for _ in range(v - 1):
                direct, end = tau(direct, 0)
                check(end == 0, ("direct fixed prefix", m))
            direct_ok = direct == word
            check(direct_ok, ("direct section word", m))
            direct_lines.append(f"{m}:{v}:{len(word)}:{int(direct_ok)}:{word_digest}")

        words.append(bytearray(word))
        gaps.append(gap)
        word = next_word

    transcript = "RULE30_ORBIT_SIGNALIZER_PROBE_V1\n" + "\n".join(transcript_lines) + "\n"
    transcript_digest = hashlib.sha256(transcript.encode()).hexdigest()
    direct_digest = hashlib.sha256(("\n".join(direct_lines) + "\n").encode()).hexdigest()
    check(transcript_digest == EXPECTED_TRANSCRIPT_SHA256, "transcript digest")
    check(direct_digest == EXPECTED_DIRECT_SHA256, "direct-section digest")
    return words, gaps, transcript_digest, direct_digest


def audit_right_action_order(words: list[bytearray]) -> tuple[str, str, int]:
    square = bytearray((B, B))
    section0, end = tau(square, 0)
    check(end == 0, "B square fixes root")
    right_action = "".join(LETTERS[letter] for letter in section0)
    check(right_action == "CB", "B square section order")
    next_section, end = tau(section0, 0)
    check(end == 0 and "".join(LETTERS[letter] for letter in next_section) == "AB", "B square second section")
    conversion_checks = 0
    for word in words:
        display = word[::-1]
        for start in (0, 1):
            right_section, _ = tau(word, start)
            display_section = thm3463_display_section(display, start)
            check(display_section == right_section[::-1], ("display convention", len(word), start))
            conversion_checks += 1
    return right_action, right_action[::-1], conversion_checks


def audit_marked_units(words: list[bytearray], gaps: list[int]) -> tuple[int, tuple[int, int, int]]:
    rows = packed_power_rows()
    modulus = 1 << TAIL_BITS
    mask = modulus - 1
    checks = 0
    q4 = (0, 0, 0)

    for m in range(DIRECT_SECTION_MAX_M + 1):
        q = 1 << m
        v = EXPECTED_V_0_TO_23[m]
        gap = gaps[m]
        word = words[m]
        next_word = words[m + 1] if m < len(words) - 1 else renormalize(word)[0]

        first = apply_word(word, 0, TAIL_BITS)
        twice = apply_word(word, first, TAIL_BITS)
        next_unit = apply_word(next_word, 0, TAIL_BITS - gap)
        direct_first = ((rows[q] - 1) >> v) & mask
        direct_twice = ((rows[2 * q] - 1) >> v) & mask
        direct_sibling = ((rows[2 * q] - rows[q]) >> v) & mask

        check(first == direct_first, ("marked first unit", m))
        check(twice == direct_twice, ("marked square tail", m))
        check((twice - first) & mask == direct_sibling, ("marked sibling unit", m))
        check(nu2(twice) == gap, ("marked gap", m))
        check((twice >> gap) == next_unit, ("marked next unit", m))
        check(borrow_decode(first, twice, TAIL_BITS) == ((twice - first) & mask), ("borrow", m))
        checks += 1

        if m == 2:
            q4 = (first, (twice - first) & mask, twice >> gap)
            check(q4 == (25, 6403, 1607), "q4 units")
    return checks, q4


def finite_rows(width: int, stop: int) -> list[int]:
    row = 1
    rows = [row]
    for _ in range(stop):
        row = phi(row, width)
        rows.append(row)
    return rows


def audit_phase_owners(words: list[bytearray], gaps: list[int]) -> tuple[int, int, int, int]:
    phase_checks = ratio_checks = transition_checks = conjugacy_checks = 0
    phase_bits = 18
    tail_mask = (1 << phase_bits) - 1

    for m in range(PHASE_OWNER_MAX_M + 1):
        q = 1 << m
        v = EXPECTED_V_0_TO_23[m]
        gap = gaps[m]
        width = v + phase_bits
        rows = finite_rows(width, 3 * q)
        power_word = bytearray((A,)) * q
        doubled_word = bytearray((A,)) * (2 * q)
        marked = words[m]

        for phase in range(q):
            prefix = rows[phase] & ((1 << v) - 1)
            tail = (rows[phase] >> v) & tail_mask
            owner = section_at_prefix(power_word, prefix, v)
            first = apply_word(owner, tail, phase_bits)
            twice = apply_word(owner, first, phase_bits)
            check(first == ((rows[phase + q] >> v) & tail_mask), ("phase first", m, phase))
            check(twice == ((rows[phase + 2 * q] >> v) & tail_mask), ("phase twice", m, phase))

            unit = (first - tail) & tail_mask
            sibling = (twice - first) & tail_mask
            total = (twice - tail) & tail_mask
            check(unit & 1 and sibling & 1, ("phase odd units", m, phase))
            check(nu2(total) == gap, ("phase common gap", m, phase))
            ratio = (-sibling * pow(unit, -1, 1 << phase_bits)) & tail_mask
            check(nu2((ratio - 1) & tail_mask) == gap, ("principal ratio", m, phase))
            phase_checks += 1
            ratio_checks += 1

            extension = tail & ((1 << gap) - 1)
            induced_next = section_at_prefix(owner + owner, extension, gap)
            long_prefix = rows[phase] & ((1 << (v + gap)) - 1)
            direct_next = section_at_prefix(doubled_word, long_prefix, v + gap)
            check(induced_next == direct_next, ("phase induced next", m, phase))
            transition_checks += 1

            transport_word = bytearray((A,)) * phase
            transport = section_at_prefix(transport_word, 1, v)
            check(apply_word(transport, 0, phase_bits) == tail, ("phase transport tail", m, phase))
            left_word = transport + owner
            right_word = marked + transport
            for test_tail in range(16):
                left_value = apply_word(left_word, test_tail, 4)
                right_value = apply_word(right_word, test_tail, 4)
                check(left_value == right_value, ("owner conjugacy portrait", m, phase, test_tail))
                conjugacy_checks += 1

    return phase_checks, ratio_checks, transition_checks, conjugacy_checks


def audit_shallow_hostile(words: list[bytearray], gaps: list[int]) -> tuple[str, str, int]:
    actions3 = generator_actions(3)
    actions4 = generator_actions(4)
    perm3_7 = permutation(words[7], actions3, {})
    perm3_17 = permutation(words[17], actions3, {})
    check(perm3_7 == tuple(range(7, -1, -1)), "s7 depth-three portrait")
    check(perm3_17 == perm3_7, "s7/s17 shared depth-three portrait")
    check(gaps[7] == 8 and gaps[17] == 5, "hostile costs")

    perm4_texts: list[str] = []
    cache: dict[bytes, tuple[int, ...]] = {}
    for word in words:
        perm4 = permutation(word, actions4, cache)
        perm4_texts.append("".join(format(value, "x") for value in perm4))
    check(len(set(perm4_texts)) == len(words), "distinct depth-four portraits")
    check(perm4_texts[7] == EXPECTED_S7_PERM4, "s7 depth-four portrait")
    check(perm4_texts[17] == EXPECTED_S17_PERM4, "s17 depth-four portrait")
    check(hashlib.sha256(words[7]).hexdigest() == EXPECTED_S7_WORD_SHA256, "s7 word digest")
    check(hashlib.sha256(words[17]).hexdigest() == EXPECTED_S17_WORD_SHA256, "s17 word digest")
    return perm4_texts[7], perm4_texts[17], len(set(perm4_texts))


def main() -> None:
    words, gaps, transcript_digest, direct_digest = signalizer_words()
    right_order, thm3463_display, display_checks = audit_right_action_order(words)
    marked_checks, q4_units = audit_marked_units(words, gaps)
    phase_checks, ratio_checks, transition_checks, conjugacy_checks = audit_phase_owners(words, gaps)
    s7_perm4, s17_perm4, distinct4 = audit_shallow_hostile(words, gaps)

    print("RULE30_ORBIT_SIGNALIZER_GAP_THM3511")
    print("status PROVED_VERIFIED_EXACT_INDEPENDENTLY_AUDITED")
    print(
        "right_action_B2_section0",
        right_order,
        "thm3463_factor_display",
        thm3463_display,
        "display_conversion_checks",
        display_checks,
    )
    print(
        "renormalization_scales",
        len(gaps),
        "m_range",
        (0, RENORMALIZATION_MAX_M),
        "terminal_v23",
        EXPECTED_V_0_TO_23[-1],
    )
    print("gap_prefix", tuple(gaps))
    print("direct_section_checks", DIRECT_SECTION_MAX_M + 1, "digest", direct_digest)
    print("marked_unit_borrow_checks", marked_checks, "q4_units", q4_units)
    print(
        "phase_owner_scales",
        (0, PHASE_OWNER_MAX_M),
        "phase_checks",
        phase_checks,
        "ratio_checks",
        ratio_checks,
        "transition_checks",
        transition_checks,
    )
    print("owner_conjugacy_depth4_checks", conjugacy_checks)
    print(
        "shallow_hostile",
        "scales",
        (7, 17),
        "shared_perm3",
        "76543210",
        "costs",
        (gaps[7], gaps[17]),
    )
    print("hostile_perm4", s7_perm4, s17_perm4, "distinct_depth4", distinct4)
    print("transcript_sha256", transcript_digest)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
