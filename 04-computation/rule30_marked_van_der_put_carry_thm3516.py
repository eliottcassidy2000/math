#!/usr/bin/env python3
"""Finite exact controls for proved THM-3516."""

from __future__ import annotations

from fractions import Fraction


A, B, C = 0, 1, 2
PATH_MAX_N = 32
POWER_MAX_M = 12
UNARY_Q_MAX = 256
MAX_TIME = 1 << POWER_MAX_M
EXPECTED_V = (1, 3, 4, 6, 7, 9, 15, 16, 24, 25, 27, 29, 34)
EXPECTED_POWER_CENTERS = (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1)


def check(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"check failed: {label}")


def nu2(value: int) -> int:
    check(value != 0, "nu2 zero")
    return (value & -value).bit_length() - 1


def phi(value: int) -> int:
    return value ^ ((value << 1) | (value << 2))


def packed_rows() -> list[int]:
    rows = [1]
    for _ in range(MAX_TIME):
        rows.append(phi(rows[-1]))
    return rows


def activity(word: bytes | bytearray) -> int:
    return sum(letter != A for letter in word) & 1


def tau(word: bytes | bytearray, start: int) -> tuple[bytearray, int]:
    """Right-action root section, with factors applied left to right."""
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


def renormalize(word: bytes | bytearray) -> tuple[bytearray, int]:
    check(activity(word) == 1, ("inactive signalizer", len(word)))
    left, left_end = tau(word, 0)
    right, right_end = tau(word, 1)
    check((left_end, right_end) == (1, 0), "active child routing")
    current = left + right
    gap = 1
    while activity(current) == 0:
        current, end = tau(current, 0)
        check(end == 0, ("fixed zero descent", gap))
        gap += 1
        check(gap <= 64, "gap probe overflow")
    return current, gap


def shell_unit(rows: list[int], valuations: tuple[int, ...], m: int, phase: int) -> int:
    q = 1 << m
    numerator = rows[phase + q] - rows[phase]
    divisor = 1 << valuations[m]
    check(numerator % divisor == 0, ("shell divisibility", m, phase))
    unit = numerator // divisor
    check(unit & 1 == 1, ("shell oddness", m, phase))
    return unit


def audit_periods_and_paths(rows: list[int], valuations: tuple[int, ...]) -> tuple[int, int, int, int]:
    period_checks = 0
    normalized_checks = 0
    shell_terms = 0
    block_terms = 0
    first_carry_hostile = None

    for n in range(1, PATH_MAX_N + 1):
        e = sum(value <= n for value in valuations)
        check(e < len(valuations), ("valuation sentinel", n, e))
        check(valuations[e] > n, ("next innovation", n, e))
        period = 1 << e
        mask = (1 << (n + 1)) - 1
        orbit = [rows[t] & mask for t in range(period)]
        check(len(set(orbit)) == period, ("distinct seed orbit", n, period))
        check((rows[period] & mask) == 1, ("seed return", n, period))
        period_checks += 1

        owner = n % period
        check((rows[n] - rows[owner]) % (1 << (n + 1)) == 0, ("owner fiber", n))
        j_n = (rows[n] - 1) // 2
        j_owner = (rows[owner] - 1) // 2
        check((j_n - j_owner) % (1 << n) == 0, ("desuspended owner fiber", n))

        owner_sum = 0
        owner_terms: list[tuple[int, int, int]] = []
        for m in range(e):
            check((m < e) == (valuations[m] <= n), ("active shell equivalence", n, m))
            if (owner >> m) & 1:
                phase = owner & ((1 << m) - 1)
                unit = shell_unit(rows, valuations, m, phase)
                term = unit << (valuations[m] - 1)
                owner_sum += term
                owner_terms.append((m, unit, term))
                shell_terms += 1
        check(owner_sum == j_owner, ("van der Put owner path", n, owner))
        center = (rows[n] >> n) & 1
        check(center == ((owner_sum >> (n - 1)) & 1), ("marked owner bit", n))

        alpha = nu2(n)
        v_alpha = valuations[alpha]
        divisor = 1 << v_alpha
        check((rows[n] - 1) % divisor == 0, ("all-n normalization", n))
        normalized = (rows[n] - 1) // divisor
        check(normalized & 1 == 1, ("normalized odd", n))

        q = 1 << alpha
        block_sum = 0
        for phase in range(0, n, q):
            block_sum += shell_unit(rows, valuations, alpha, phase)
            block_terms += 1
        check(block_sum == normalized, ("block telescoping", n))

        shell_sum = 0
        active_normalized_terms: list[int] = []
        for m in range(n.bit_length()):
            if (n >> m) & 1:
                phase = n & ((1 << m) - 1)
                unit = shell_unit(rows, valuations, m, phase)
                term = (1 << (valuations[m] - v_alpha)) * unit
                shell_sum += term
                if valuations[m] <= n:
                    precision = n - valuations[m] + 1
                    active_normalized_terms.append(term % (1 << (n - v_alpha + 1)))
        check(shell_sum == normalized, ("normalized shell factorization", n))

        if v_alpha > n:
            check(center == 0, ("pre-live center zero", n, v_alpha))
        else:
            target = n - v_alpha
            check(center == ((normalized >> target) & 1), ("normalized target bit", n, target))
            modulus = 1 << (target + 1)
            top = 1 << target
            truncated = [term % modulus for term in active_normalized_terms]
            direct_xor = sum((term >> target) & 1 for term in truncated) & 1
            lower_sum = sum(term % top for term in truncated)
            carry_parity = (lower_sum // top) & 1
            check(center == (direct_xor ^ carry_parity), ("carry decomposition", n))
            if carry_parity and first_carry_hostile is None:
                first_carry_hostile = n
        normalized_checks += 1

    check(first_carry_hostile == 6, ("first physical carry hostile", first_carry_hostile))
    return period_checks, normalized_checks, shell_terms, block_terms


def audit_n6(rows: list[int], valuations: tuple[int, ...]) -> tuple[int, ...]:
    u_1_0 = shell_unit(rows, valuations, 1, 0)
    u_2_2 = shell_unit(rows, valuations, 2, 2)
    j_terms = (u_1_0 << (valuations[1] - 1), u_2_2 << (valuations[2] - 1))
    normalized_terms = (u_1_0, u_2_2 << (valuations[2] - valuations[1]))
    j_6 = (rows[6] - 1) // 2
    w_6 = (rows[6] - 1) >> valuations[1]
    check(rows[6] == 6409 and j_6 == 3204 and w_6 == 801, "n6 packed constants")
    check((u_1_0, u_2_2) == (3, 399), "n6 shell units")
    check(j_terms == (12, 3192) and sum(j_terms) == j_6, "n6 J terms")
    check(normalized_terms == (3, 798) and sum(normalized_terms) == w_6, "n6 W terms")
    j_residues = tuple(term % 64 for term in j_terms)
    w_residues = tuple(term % 16 for term in normalized_terms)
    direct_xor = ((w_residues[0] >> 3) ^ (w_residues[1] >> 3)) & 1
    lower_carry = ((w_residues[0] % 8 + w_residues[1] % 8) // 8) & 1
    center = (rows[6] >> 6) & 1
    check(j_residues == (12, 56), "n6 J residues")
    check(w_residues == (3, 14), "n6 W residues")
    check((direct_xor, lower_carry, center) == (1, 1, 0), "n6 carry toggle")
    return (rows[6], j_6, w_6, *j_terms, *j_residues, direct_xor, lower_carry, center)


def audit_projective_origin(rows: list[int], valuations: tuple[int, ...]) -> int:
    units = [shell_unit(rows, valuations, m, 0) for m in range(POWER_MAX_M + 1)]
    check(units[0] == 3, "origin gauge")
    product = Fraction(units[0], 1)
    checks = 0
    for m in range(POWER_MAX_M):
        shifted = shell_unit(rows, valuations, m, 1 << m)
        g = -Fraction(shifted, units[m])
        gap = valuations[m + 1] - valuations[m]
        z_from_g = (1 - g) / (1 << gap)
        z_direct = Fraction(units[m + 1], units[m])
        check(z_from_g == z_direct, ("projective vertical unit", m))
        product *= z_from_g
        check(product == units[m + 1], ("projective product", m))
        checks += 1
    return checks


def audit_signalizer_power_sections(
    rows: list[int], valuations: tuple[int, ...]
) -> tuple[tuple[int | None, ...], tuple[int, ...], int]:
    word = bytearray((B,))
    section_depths: list[int | None] = []
    centers: list[int] = []
    universal_checks = 0
    for m in range(POWER_MAX_M + 1):
        q = 1 << m
        v = valuations[m]
        center = (rows[q] >> q) & 1

        direct_signalizer = bytearray((B,)) * q
        for _ in range(v - 1):
            direct_signalizer, _ = tau(direct_signalizer, 0)
        check(direct_signalizer == word, ("direct marked signalizer", m, q, v))

        universal_section = bytearray((B,)) * q
        for _ in range(q - 1):
            universal_section, _ = tau(universal_section, 0)
        universal_bit = activity(universal_section)
        check(universal_bit == center, ("universal B-power section", m, q))
        universal_checks += 1

        if v > q:
            depth = None
            section_bit = 0
        else:
            depth = q - v
            section = bytearray(word)
            for _ in range(depth):
                section, _ = tau(section, 0)
            section_bit = activity(section)
            check(section == universal_section, ("section-depth cancellation", m, q, v))
        check(section_bit == center, ("power section center", m, q, v, depth))
        if v == q:
            check(depth == 0 and center == 1, ("equality edge", m))
        section_depths.append(depth)
        centers.append(center)
        if m < POWER_MAX_M:
            next_word, gap = renormalize(word)
            check(gap == valuations[m + 1] - valuations[m], ("signalizer gap", m))
            word = next_word
    check(tuple(centers) == EXPECTED_POWER_CENTERS, "power center prefix")
    return tuple(section_depths), tuple(centers), universal_checks


def audit_unary_future_output(rows: list[int]) -> int:
    total_classes = 0
    for q in range(1, UNARY_Q_MAX + 1):
        word = bytearray((B,)) * q
        trace: list[int] = []
        raw_words: set[bytes] = set()
        for depth in range(2 * q + 1):
            trace.append(activity(word))
            check(trace[-1] == ((rows[q] >> (depth + 1)) & 1), ("unary output", q, depth))
            frozen = bytes(word)
            check(frozen not in raw_words, ("unary raw repeat", q, depth))
            raw_words.add(frozen)
            if depth < 2 * q:
                word, _ = tau(word, 0)
        check(trace[2 * q - 1] == 1 and trace[2 * q] == 0, ("unary terminal trace", q))
        check(word == bytearray((A,)) * q, ("unary terminal word", q))
        fixed, _ = tau(word, 0)
        check(fixed == word, ("unary terminal fixed", q))
        suffixes = {
            tuple(trace[index:] + [0] * index) for index in range(2 * q + 1)
        }
        check(len(suffixes) == 2 * q + 1, ("unary future classes", q))
        total_classes += len(suffixes)
    return total_classes


def main() -> None:
    rows = packed_rows()
    valuations = tuple(nu2(rows[1 << m] - 1) for m in range(POWER_MAX_M + 1))
    check(valuations == EXPECTED_V, "valuation prefix")
    period_checks, normalized_checks, shell_terms, block_terms = audit_periods_and_paths(
        rows, valuations
    )
    n6 = audit_n6(rows, valuations)
    projective_checks = audit_projective_origin(rows, valuations)
    section_depths, centers, universal_checks = audit_signalizer_power_sections(rows, valuations)
    unary_classes = audit_unary_future_output(rows)

    print("RULE30_MARKED_VAN_DER_PUT_CARRY_THM3516")
    print("status PROVED_VERIFIED_EXACT_INDEPENDENTLY_AUDITED")
    print("packed_universe", (0, MAX_TIME), "exact_untruncated_rows", len(rows))
    print("valuation_scales", (0, POWER_MAX_M), "v_prefix", valuations)
    print(
        "path_universe",
        (1, PATH_MAX_N),
        "period_checks",
        period_checks,
        "normalized_checks",
        normalized_checks,
        "van_der_put_terms",
        shell_terms,
        "block_terms",
        block_terms,
    )
    print(
        "n6_carry",
        "row_J_W",
        n6[:3],
        "J_terms",
        n6[3:5],
        "J_mod64",
        n6[5:7],
        "direct_xor_lowercarry_center",
        n6[7:],
    )
    print("projective_origin_checks", projective_checks, "calibrated_gauge", 3)
    print("power_section_depths", section_depths)
    print("universal_B_power_sections", universal_checks, "input_depth_rule", "q_minus_1")
    print("power_center_prefix", centers)
    print(
        "unary_future_output_q_range",
        (1, UNARY_Q_MAX),
        "total_classes",
        unary_classes,
        "class_rule",
        "2q_plus_1",
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
