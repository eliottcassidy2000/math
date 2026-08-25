#!/usr/bin/env python3
"""Independent native-first exact audit for THM-4072."""

from fractions import Fraction
from hashlib import sha256


MAX_EXHAUSTIVE_DEPTH = 16
GREEDY_DEPTH = 640


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def greedy_digits(length):
    state = Fraction(1)
    digits = []
    for index in range(length):
        scaled_numerator = 3 * state.numerator
        scaled_denominator = 2 * state.denominator
        digit = scaled_numerator // scaled_denominator
        require(digit in (0, 1), f"greedy alphabet drift at index={index}")
        state = Fraction(scaled_numerator, scaled_denominator) - digit
        require(Fraction(0) < state < Fraction(1),
                f"greedy state drift at index={index}")
        digits.append(digit)
    return tuple(digits)


def integer_itinerary(start, length):
    """Recover carry letters directly from the ordinary ceiling orbit."""
    value = start
    word = []
    for _ in range(length):
        carry = value & 1
        word.append(carry)
        value = (3 * value + carry) // 2
    return tuple(word), value


def carry_polynomial(word):
    length = len(word)
    return sum(
        digit * 2**index * 3 ** (length - index - 1)
        for index, digit in enumerate(word)
    )


def lex_safe(word, greedy):
    """Independent lexicographic form of all finite suffix gates."""
    return all(
        tuple(word[start:]) <= greedy[: len(word) - start]
        for start in range(len(word))
    )


def renewal_parse(word, greedy):
    """Parse first downward disagreements instead of stepping a transducer."""
    position = 0
    blocks = []
    while position < len(word):
        matched = 0
        while (
            position + matched < len(word)
            and word[position + matched] == greedy[matched]
        ):
            matched += 1
        if position + matched == len(word):
            return True, tuple(blocks), matched
        require(
            greedy[matched] == 1 and word[position + matched] == 0,
            f"renewal parser met an upward disagreement at position={position}",
        )
        blocks.append(matched + 1)
        position += matched + 1
    return True, tuple(blocks), 0


def advance_from_state(state, extension, greedy):
    for digit in extension:
        expected = greedy[state]
        if digit == expected:
            state += 1
        elif expected == 1 and digit == 0:
            state = 0
        else:
            return False, state
    return True, state


def rational_periodic_tail(block, phase=0):
    rotated = tuple(block[phase:]) + tuple(block[:phase])
    weighted = Fraction(0)
    power = Fraction(2, 3)
    for digit in rotated:
        weighted += digit * power
        power *= Fraction(2, 3)
    return weighted / (1 - Fraction(2, 3) ** len(block))


def two_adic_rational_bits(numerator, denominator, length):
    modulus = 2**length
    residue = (numerator * pow(denominator, -1, modulus)) % modulus
    return tuple((residue >> index) & 1 for index in range(length))


def truncated_rational_tail(word):
    value = Fraction(0)
    weight = Fraction(2, 3)
    for digit in word:
        value += digit * weight
        weight *= Fraction(2, 3)
    return value


def fractional_part(value):
    return value - value.numerator // value.denominator


def audit_native_first(greedy):
    total = 0
    safe_counts = [1]
    state_counts = [1]
    terminal_extension_counts = [1]
    safe_level_16 = []
    clock_gates = 0

    for length in range(1, MAX_EXHAUSTIVE_DEPTH + 1):
        carry_words = set()
        safe_words = set()
        terminal_safe_prefixes = 0
        final_states = set()
        modulus = 2**length
        for native_value in range(modulus):
            total += 1
            word, terminal_orbit_value = integer_itinerary(native_value, length)
            require(word not in carry_words,
                    f"native-first carry collision at length={length}")
            carry_words.add(word)

            polynomial = carry_polynomial(word)
            affine_numerator = 3**length * native_value + polynomial
            require(affine_numerator % modulus == 0,
                    f"native-first affine divisibility drift at value={native_value}")
            require(affine_numerator // modulus == terminal_orbit_value,
                    f"native-first orbit identity drift at value={native_value}")

            safe = lex_safe(word, greedy)
            if safe:
                parsed, blocks, final_state = renewal_parse(word, greedy)
                require(parsed, f"renewal parse drift at word={word}")
                require(sum(blocks) + final_state == length,
                        f"renewal length drift at word={word}")
                safe_words.add(word)
                final_states.add(final_state)
                terminal_safe_prefixes += 1
                if length == MAX_EXHAUSTIVE_DEPTH:
                    safe_level_16.append(native_value)

            # An all-zero native prefix receives a 1 next; every other prefix
            # already has a positive zero-tail completion.
            positive_completion = native_value if native_value else modulus
            completion_word, _ = integer_itinerary(positive_completion, length)
            require(completion_word == word,
                    f"positive terminal cylinder extension drift at value={native_value}")

            # Audit an independently chosen nonzero remaining completion k.
            k = (native_value % 31) + 1
            completed_integer = native_value + modulus * k
            require(completed_integer % modulus == native_value,
                    f"ordinary prefix congruence drift at value={native_value}")
            quotient = completed_integer >> length
            require(quotient == k,
                    f"ordinary quotient drift at value={native_value}")
            height = quotient.bit_length()
            while height:
                next_quotient = quotient >> 1
                require(next_quotient.bit_length() == height - 1,
                        f"native clock decrement drift at value={native_value}")
                quotient = next_quotient
                height -= 1
                clock_gates += 1
            require(quotient == 0,
                    f"native clock terminal drift at value={native_value}")

        require(len(carry_words) == modulus,
                f"native-first level bijection drift at length={length}")
        require(final_states == set(range(length + 1)),
                f"native-first follower-state scout drift at length={length}")
        safe_counts.append(len(safe_words))
        state_counts.append(len(final_states))
        terminal_extension_counts.append(terminal_safe_prefixes)

    expected = (1, 2, 3, 5, 8, 12, 18, 27, 40, 60, 90, 134,
                201, 302, 452, 678, 1018)
    require(tuple(safe_counts) == expected, "native-first safe counts drift")
    require(tuple(terminal_extension_counts) == expected,
            "terminal-prefix no-pruning counts drift")
    require(total == 131070, "native-first exhaustive count drift")
    digest_input = ",".join(map(str, safe_level_16)).encode("ascii")
    return (
        total,
        tuple(safe_counts),
        tuple(state_counts),
        tuple(terminal_extension_counts),
        clock_gates,
        sha256(digest_input).hexdigest(),
    )


def audit_follower_separation(greedy):
    depth = MAX_EXHAUSTIVE_DEPTH
    for state in range(depth + 1):
        witness = (0,) * (depth - state) + greedy[:state]
        require(lex_safe(witness, greedy),
                f"state witness unsafe at state={state}")
        _, _, final_state = renewal_parse(witness, greedy)
        require(final_state == state,
                f"state witness ended at {final_state}, expected {state}")

    separators = []
    for left in range(depth + 1):
        for right in range(left + 1, depth + 1):
            offset = 0
            while greedy[left + offset] == greedy[right + offset]:
                offset += 1
                require(offset < GREEDY_DEPTH - depth,
                        f"follower separator search escaped at {left},{right}")
            if greedy[left + offset] == 1:
                larger = left
                smaller = right
            else:
                larger = right
                smaller = left
            extension = greedy[larger: larger + offset + 1]
            larger_ok, _ = advance_from_state(larger, extension, greedy)
            smaller_ok, _ = advance_from_state(smaller, extension, greedy)
            require(larger_ok and not smaller_ok,
                    f"follower separator failed at states={left},{right}")
            separators.append(offset + 1)

    boundary = greedy[:256]
    _, boundary_blocks, boundary_state = renewal_parse(boundary, greedy)
    require(not boundary_blocks and boundary_state == 256,
            "greedy equality ray acquired a false reset")
    zero_predecessor = (0,) + greedy[:255]
    require(lex_safe(zero_predecessor, greedy),
            "zero predecessor of equality ray left K")
    _, predecessor_blocks, predecessor_state = renewal_parse(
        zero_predecessor, greedy
    )
    require(predecessor_blocks == (1,) and predecessor_state == 255,
            "zero predecessor boundary reset drift")
    shifted = greedy[1:257]
    require(lex_safe(shifted, greedy), "shifted greedy control left K")
    _, shifted_blocks, _ = renewal_parse(shifted, greedy)
    require(len(shifted_blocks) > 100,
            "shifted greedy strict control lost recurrent resets")
    return len(separators), max(separators), len(shifted_blocks)


def audit_hostiles(greedy):
    equality_prefix = "".join(map(str, greedy[:64]))

    tails_100 = tuple(rational_periodic_tail((1, 0, 0), phase)
                      for phase in range(3))
    require(tails_100 == (Fraction(18, 19), Fraction(8, 19), Fraction(12, 19)),
            "independent (100)^infinity tail drift")
    native_100 = "".join(map(str, two_adic_rational_bits(-9, 19, 18)))
    require(native_100 == "101100001010011110",
            "independent denominator-19 2-adic expansion drift")

    a1_word, _ = integer_itinerary(1, 8)
    require(a1_word[:4] == (1, 0, 1, 1), "independent A=1 carry drift")
    require(lex_safe(a1_word[:3], greedy) and not lex_safe(a1_word[:4], greedy),
            "independent A=1 boundary drift")
    a1_tail = truncated_rational_tail(a1_word[:4])
    require(a1_tail == Fraction(94, 81), "independent A=1 tail drift")

    native_01 = "".join(map(str, two_adic_rational_bits(-2, 5, 16)))
    require(native_01 == "0110" * 4,
            "independent (01)^infinity 2-adic expansion drift")
    first_01_failure = next(
        length for length in range(1, 17)
        if not lex_safe(((0, 1) * 8)[:length], greedy)
    )
    require(first_01_failure == 6, "independent (01) first failure drift")
    hostile_01 = truncated_rational_tail(((0, 1) * 3)[1:6])
    require(hostile_01 == Fraction(266, 243),
            "independent (01) hostile suffix drift")

    tower_rows = []
    formal_carry = (1, 0, 0) * 24
    formal_native = tuple(map(int, "101100001010011110")) * 4
    for horizon in (18, 36, 54, 72):
        integer_start = 9 * (2**horizon - 1) // 19
        recovered, _ = integer_itinerary(integer_start, horizon)
        require(recovered == formal_carry[:horizon],
                f"independent tower carry drift at horizon={horizon}")
        native = tuple((integer_start >> index) & 1 for index in range(horizon))
        require(native == formal_native[:horizon],
                f"independent tower native drift at horizon={horizon}")
        require(lex_safe(recovered, greedy),
                f"independent tower safe-prefix drift at horizon={horizon}")
        alpha = Fraction(9 * 2**horizon, 19)
        phase_h = fractional_part(alpha * Fraction(3, 2) ** horizon)
        phase_next = fractional_part(alpha * Fraction(3, 2) ** (horizon + 1))
        require(phase_h == Fraction(9, 19),
                f"independent tower last safe phase drift at horizon={horizon}")
        require(phase_next == Fraction(27, 38),
                f"independent tower hostile phase drift at horizon={horizon}")
        for index in range(horizon + 1):
            phase = fractional_part(alpha * Fraction(3, 2) ** index)
            require(phase < Fraction(1, 2),
                    f"independent tower premature failure at horizon={horizon}")
        tower_rows.append((horizon, integer_start.bit_length(), phase_next))

    return (
        equality_prefix,
        tails_100,
        native_100,
        a1_tail,
        native_01[:4],
        hostile_01,
        tuple(tower_rows),
    )


def main():
    greedy = greedy_digits(GREEDY_DEPTH)
    (
        total,
        safe_counts,
        state_counts,
        terminal_counts,
        clock_gates,
        safe_digest,
    ) = audit_native_first(greedy)
    separator_count, max_separator, shifted_resets = audit_follower_separation(greedy)
    (
        equality_prefix,
        tails_100,
        native_100,
        a1_tail,
        native_01,
        hostile_01,
        tower_rows,
    ) = audit_hostiles(greedy)

    print("THM-4072 independent native-first exact audit")
    print(f"exhaustive_depth={MAX_EXHAUSTIVE_DEPTH}")
    print(f"ordinary_prefixes_checked={total}")
    print("lex_safe_counts_0_to_16=" + ",".join(map(str, safe_counts)))
    print("renewal_states_0_to_16=" + ",".join(map(str, state_counts)))
    print("terminal_extendible_safe_prefixes_0_to_16="
          + ",".join(map(str, terminal_counts)))
    print(f"native_clock_descent_gates={clock_gates}")
    print(f"safe_native_level16_sha256={safe_digest}")
    print(f"pairwise_follower_separators={separator_count}")
    print(f"maximum_scout_separator_length={max_separator}")
    print(f"shifted_greedy_resets_in_256={shifted_resets}")
    print(f"greedy_equality_prefix_64={equality_prefix}")
    print("periodic_100_tails=" + ",".join(map(str, tails_100)))
    print(f"periodic_100_native_block_18={native_100}")
    print(f"terminal_A1_truncated_tail={a1_tail}")
    print(f"centered_01_native_block_4={native_01}")
    print(f"centered_01_hostile_suffix={hostile_01}")
    print("tower_horizon_bitlength_hostile=" + ";".join(
        f"{horizon}:{bits}:{hostile}"
        for horizon, bits, hostile in tower_rows
    ))
    print("result=PASS")


if __name__ == "__main__":
    main()
