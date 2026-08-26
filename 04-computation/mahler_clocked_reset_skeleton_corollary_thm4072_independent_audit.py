#!/usr/bin/env python3
"""Clean-room exact audit for the THM-4072 clocked-reset corollary.

This audit imports no code from the primary probe.  The all-depth statements
follow algebraically from the affine ceiling recurrence; this script supplies
independent finite hostile controls from that defining map.
"""

from __future__ import annotations

import hashlib
import sys
from fractions import Fraction
from itertools import product


sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def T(value: int) -> int:
    require(value >= 0, "ceiling orbit left the nonnegative integers")
    return (3 * value + (value & 1)) // 2


def orbit(value: int, steps: int) -> int:
    for _ in range(steps):
        value = T(value)
    return value


def carries(value: int, steps: int) -> tuple[int, ...]:
    word = []
    for _ in range(steps):
        word.append(value & 1)
        value = T(value)
    return tuple(word)


def greedy_word(length: int) -> tuple[int, ...]:
    remainder = Fraction(1)
    result = []
    for _ in range(length):
        scaled = 3 * remainder / 2
        digit = scaled.numerator // scaled.denominator
        require(digit in (0, 1), "greedy digit left the binary alphabet")
        result.append(digit)
        remainder = scaled - digit
    return tuple(result)


GREEDY = greedy_word(160)


def follower(word: tuple[int, ...], start: int = 0) -> tuple[bool, int, int]:
    state = start
    for index, digit in enumerate(word):
        boundary = GREEDY[state]
        if digit == boundary:
            state += 1
        elif digit == 0 and boundary == 1:
            state = 0
        else:
            return False, state, index
    return True, state, len(word)


def canonical_residue(word: tuple[int, ...]) -> int:
    depth = len(word)
    if depth == 0:
        return 0
    modulus = 1 << depth
    carry_sum = sum(
        digit * (1 << index) * 3 ** (depth - 1 - index)
        for index, digit in enumerate(word)
    )
    return (-pow(3**depth, -1, modulus) * carry_sum) % modulus


def reset_block(length: int) -> tuple[int, ...]:
    require(length >= 1 and GREEDY[length - 1] == 1,
            "reset length is not a greedy-one position")
    return GREEDY[: length - 1] + (0,)


def append_reset(
    residue: int,
    depth: int,
    final_state: int,
    length: int,
) -> tuple[int, int, int]:
    """Append one reset block to an affine cone.

    Input means T^depth(residue+2^depth*k)=final_state+3^depth*k.
    """

    block_residue = canonical_residue(reset_block(length))
    modulus = 1 << length
    quotient_class = (
        (block_residue - final_state) * pow(3**depth, -1, modulus)
    ) % modulus
    new_residue = residue + (1 << depth) * quotient_class
    launch = final_state + 3**depth * quotient_class
    require(carries(launch, length) == reset_block(length),
            "appended launch missed its reset block")
    return new_residue, depth + length, orbit(launch, length)


def skeletons(total_cap: int, lengths: tuple[int, ...]):
    stack = [((), 0)]
    while stack:
        skeleton, total = stack.pop()
        if skeleton:
            yield skeleton
        for length in reversed(lengths):
            if total + length <= total_cap:
                stack.append((skeleton + (length,), total + length))


def decision_after_terminalization(value: int) -> tuple[int, bool, int, int]:
    depth = value.bit_length()
    word = carries(value, depth)
    safe, state, _ = follower(word)
    require(safe, "requested terminal prefix is not follower-safe")
    terminal_state = orbit(value, depth)
    next_digit = terminal_state & 1
    survives, _, _ = follower((next_digit,), state)
    return state, survives, terminal_state, next_digit


def first_rejection(value: int, cap: int = 256) -> int:
    state = 0
    for index in range(cap):
        digit = value & 1
        survives, state, _ = follower((digit,), state)
        if not survives:
            return index
        value = T(value)
    raise RuntimeError("hostile did not reject inside the exact cap")


def main() -> None:
    # 1. Reconstruct every carry/native cylinder through depth ten and audit
    # the affine completion identity and synchronized transition.
    cone_gates = 0
    transition_gates = 0
    for depth in range(0, 11):
        for word in product((0, 1), repeat=depth):
            residue = canonical_residue(word)
            native_state = orbit(residue, depth)
            require(carries(residue, depth) == word,
                    "canonical residue missed its carry word")
            for quotient in range(32):
                launch = residue + (1 << depth) * quotient
                require(carries(launch, depth) == word,
                        "completion left its carry cylinder")
                require(
                    orbit(launch, depth) == native_state + 3**depth * quotient,
                    "clocked affine completion identity failed",
                )
                cone_gates += 2

                native_digit = quotient & 1
                carry_digit = (native_state + native_digit) & 1
                next_quotient = (quotient - native_digit) // 2
                next_residue = residue + (1 << depth) * native_digit
                next_native_state = (
                    3 * native_state
                    + carry_digit
                    + 3 ** (depth + 1) * native_digit
                ) // 2
                require(
                    orbit(next_residue, depth + 1) == next_native_state,
                    "synchronized native-state update failed",
                )
                require(
                    orbit(launch, depth + 1)
                    == next_native_state + 3 ** (depth + 1) * next_quotient,
                    "synchronized completion transition failed",
                )
                require(
                    next_quotient.bit_length()
                    == max(0, quotient.bit_length() - 1),
                    "radix clock failed to descend",
                )
                transition_gates += 3

    # 2. Reconstruct the reset blocks and their induced affine maps without
    # importing the proposed rho/nu tables.
    reset_lengths = tuple(
        index + 1 for index, digit in enumerate(GREEDY[:128]) if digit == 1
    )
    reset_map_gates = 0
    xor_gates = 0
    for length in reset_lengths:
        block = reset_block(length)
        rho = canonical_residue(block)
        nu = orbit(rho, length)
        equality_residue = canonical_residue(GREEDY[:length])
        require(rho ^ equality_residue == 1 << (length - 1),
                "reset/equality native residues did not differ in the last bit")
        xor_gates += 1
        for quotient in range(64):
            launch = rho + (1 << length) * quotient
            require(carries(launch, length) == block,
                    "single-reset cone missed its block")
            require(orbit(launch, length) == nu + 3**length * quotient,
                    "single-reset affine map failed")
            reset_map_gates += 2

    # The second block condition has exactly one quotient class modulo 2^j.
    congruence_gates = 0
    short_lengths = tuple(length for length in reset_lengths if length <= 12)
    for first in short_lengths:
        rho_first = canonical_residue(reset_block(first))
        nu_first = orbit(rho_first, first)
        for second in short_lengths:
            rho_second = canonical_residue(reset_block(second))
            modulus = 1 << second
            predicted = (
                (rho_second - nu_first) * pow(3**first, -1, modulus)
            ) % modulus
            actual = tuple(
                quotient
                for quotient in range(modulus)
                if carries(nu_first + 3**first * quotient, second)
                == reset_block(second)
            )
            require(actual == (predicted,),
                    "second reset block did not have one quotient class")
            congruence_gates += modulus + 1

    # 3. Compose every complete reset skeleton of total depth at most 18.
    # The count 2290 is independently recovered from the greedy word.
    lengths18 = tuple(length for length in reset_lengths if length <= 18)
    all_skeletons = tuple(skeletons(18, lengths18))
    require(len(all_skeletons) == 2290, "renewal skeleton count changed")
    skeleton_gates = 0
    for skeleton in all_skeletons:
        residue = depth = final_state = 0
        word: tuple[int, ...] = ()
        for length in skeleton:
            residue, depth, final_state = append_reset(
                residue, depth, final_state, length
            )
            word += reset_block(length)
        require(depth == len(word) == sum(skeleton),
                "skeleton depth ledger failed")
        require(residue == canonical_residue(word),
                "composed cone disagreed with direct residue")
        require(orbit(residue, depth) == final_state,
                "composed cone final state failed")
        require(follower(word) == (True, 0, depth),
                "reset skeleton did not return to follower state zero")
        for quotient in range(8):
            launch = residue + (1 << depth) * quotient
            require(carries(launch, depth) == word,
                    "skeleton cone missed its carry word")
            require(
                orbit(launch, depth) == final_state + 3**depth * quotient,
                "skeleton affine final-state map failed",
            )
            skeleton_gates += 2

    # 4. Independently locate the first clock-only quotient collision.
    first_collision = None
    for depth in range(1, 9):
        groups: dict[int, dict[bool, list[int]]] = {}
        for value in range(1 << (depth - 1), 1 << depth):
            word = carries(value, depth)
            safe, state, _ = follower(word)
            if not safe:
                continue
            terminal_state = orbit(value, depth)
            next_survives, _, _ = follower((terminal_state & 1,), state)
            groups.setdefault(state, {}).setdefault(next_survives, []).append(value)
        candidates = []
        for state, by_decision in groups.items():
            if set(by_decision) == {False, True}:
                candidates.append(
                    (min(by_decision[False]), min(by_decision[True]), state)
                )
        if candidates:
            first_collision = (depth,) + min(candidates)
            break
    require(first_collision == (4, 8, 13, 1),
            "minimal clock-only quotient hostile changed")
    state8, survives8, terminal8, next8 = decision_after_terminalization(8)
    state13, survives13, terminal13, next13 = decision_after_terminalization(13)
    require((state8, survives8, terminal8, next8) == (1, False, 41, 1),
            "A=8 hostile ledger changed")
    require((state13, survives13, terminal13, next13) == (1, True, 68, 0),
            "A=13 hostile ledger changed")
    require((first_rejection(8), first_rejection(13)) == (4, 10),
            "hostile rejection indices changed")

    semantic_rows = (
        "universe=T(a)=(3a+(a mod2))/2;nonnegative integers",
        "cone=A=r_m+2^m*k=>T^m(A)=u_m+3^m*k",
        "transition=e=kmod2;c=(u_m+e)mod2;k'=floor(k/2)",
        "reset=b_ell=d[:ell-1]0;rho xor equality=2^(ell-1)",
        "reset_map=rho+2^ell*k=>nu+3^ell*k",
        "finite_skeleton=unique residue modulo 2^sum(ell_i)",
        "hostile=depth4;A8/q1/u41/reject4;A13/q1/u68/reject10",
        "scope=finite exact controls;no infinite survivor or rejection theorem",
    )
    semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("MAHLER CLOCKED-RESET COROLLARY THM-4072 INDEPENDENT AUDIT")
    print(f"checks={CHECKS}")
    print(f"cone_gates={cone_gates};transition_gates={transition_gates}")
    print(f"reset_lengths_through_128={len(reset_lengths)}")
    print(f"reset_map_gates={reset_map_gates};xor_gates={xor_gates}")
    print(f"second_reset_congruence_gates={congruence_gates}")
    print(f"skeletons_total_depth_le_18={len(all_skeletons)}")
    print(f"skeleton_affine_gates={skeleton_gates}")
    print("minimal_clock_quotient_hostile=depth4;A8_vs_A13;q=1")
    print("A8=T4=41;next1;reject_index4")
    print("A13=T4=68;next0;reject_index10")
    print(f"semantic_sha256={semantic}")
    print("verdict=EXACT_IDENTITIES_CONFIRMED;NOVELTY_LIMITED")


if __name__ == "__main__":
    main()
