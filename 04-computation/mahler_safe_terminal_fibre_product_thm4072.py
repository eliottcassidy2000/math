#!/usr/bin/env python3
"""Primary exact audit for THM-4072's Mahler fibre product."""

from fractions import Fraction
from itertools import product


MAX_EXHAUSTIVE_DEPTH = 16
GREEDY_AUDIT_DEPTH = 256
PERIODIC_NATIVE_BLOCK = "101100001010011110"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def greedy_boundary(length):
    """Return the first ``length`` greedy digits of the equality-one tail."""
    state = Fraction(1)
    digits = []
    partial = Fraction(0)
    for index in range(length):
        scaled = Fraction(3, 2) * state
        digit = scaled.numerator // scaled.denominator
        require(digit in (0, 1), f"greedy alphabet drift at index={index}")
        digits.append(digit)
        partial += digit * Fraction(2, 3) ** (index + 1)
        state = scaled - digit
        require(Fraction(0) < state < Fraction(1),
                f"greedy state escaped at index={index}")
        require(
            partial + Fraction(2, 3) ** (index + 1) * state == 1,
            f"greedy identity drift at index={index}",
        )
    return tuple(digits)


def carry_integer(word):
    """C_m=sum c_j 2^j 3^(m-1-j)."""
    length = len(word)
    return sum(
        digit * 2**index * 3 ** (length - 1 - index)
        for index, digit in enumerate(word)
    )


def canonical_residue(word):
    length = len(word)
    if length == 0:
        return 0
    modulus = 2**length
    return (-pow(3, -length, modulus) * carry_integer(word)) % modulus


def truncated_tail(word):
    return sum(
        digit * Fraction(2, 3) ** (index + 1)
        for index, digit in enumerate(word)
    )


def direct_finite_safe(word):
    """Test every strict truncated suffix inequality defining P_m."""
    return all(
        2 * carry_integer(word[start:]) < 3 ** (len(word) - start)
        for start in range(len(word))
    )


def follower_run(word, greedy):
    """Run the countable renewal follower graph from state zero."""
    state = 0
    resets = []
    for index, digit in enumerate(word):
        expected = greedy[state]
        if digit == expected:
            state += 1
        elif expected == 1 and digit == 0:
            state = 0
            resets.append(index)
        else:
            return False, state, tuple(resets), index
    return True, state, tuple(resets), None


def integer_step(value):
    carry = value & 1
    return (3 * value + carry) // 2, carry


def recover_carry_word(value, length):
    word = []
    for _ in range(length):
        value, carry = integer_step(value)
        word.append(carry)
    return tuple(word)


def carry_to_native(word):
    """Use the new u-sidecar transition to compute theta_m(word)."""
    native = []
    residue = 0
    carry_sum = 0
    u_state = 0
    history = [(residue, u_state)]
    for index, carry in enumerate(word):
        native_digit = carry ^ (u_state & 1)
        require(native_digit == (u_state + carry) % 2,
                f"native parity drift at index={index}")
        next_residue = residue + native_digit * 2**index
        next_carry_sum = 3 * carry_sum + carry * 2**index
        numerator = (
            3 * u_state
            + carry
            + 3 ** (index + 1) * native_digit
        )
        require(numerator % 2 == 0,
                f"u transition parity drift at index={index}")
        next_u = numerator // 2
        require(
            3 ** (index + 1) * next_residue + next_carry_sum
            == 2 ** (index + 1) * next_u,
            f"affine u identity drift at index={index}",
        )
        require(
            next_residue == canonical_residue(word[: index + 1]),
            f"canonical residue drift at index={index}",
        )
        check_u = next_residue
        for _ in range(index + 1):
            check_u, _ = integer_step(check_u)
        require(check_u == next_u, f"T-iterate u drift at index={index}")
        native.append(native_digit)
        residue = next_residue
        carry_sum = next_carry_sum
        u_state = next_u
        history.append((residue, u_state))
    return tuple(native), tuple(history)


def periodic_tail(block, phase=0):
    rotated = tuple(block[phase:]) + tuple(block[:phase])
    numerator = sum(
        digit * Fraction(2, 3) ** (index + 1)
        for index, digit in enumerate(rotated)
    )
    return numerator / (1 - Fraction(2, 3) ** len(block))


def fractional_part(value):
    return value - value.numerator // value.denominator


def binary_word(value, length):
    return tuple((value >> index) & 1 for index in range(length))


def maximum_zero_run_periodic(block):
    require(any(block), "zero periodic word has no bounded zero plateau")
    doubled = tuple(block) + tuple(block)
    best = 0
    current = 0
    for digit in doubled:
        if digit == 0:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return min(best, len(block) - 1)


def audit_exhaustive(greedy):
    total_words = 0
    safe_counts = [1]
    state_counts = [1]
    pair_graph_counts = [1]
    clock_gates = 0
    for length in range(1, MAX_EXHAUSTIVE_DEPTH + 1):
        native_images = set()
        safe_count = 0
        safe_states = set()
        for word in product((0, 1), repeat=length):
            total_words += 1
            accepted, state, _, rejection = follower_run(word, greedy)
            direct = direct_finite_safe(word)
            require(accepted == direct,
                    f"follower/direct safe mismatch at word={word}")
            require((rejection is None) == accepted,
                    f"rejection marker drift at word={word}")

            native, history = carry_to_native(word)
            residue = sum(digit * 2**index for index, digit in enumerate(native))
            require(residue == canonical_residue(word),
                    f"theta residue drift at word={word}")
            require(recover_carry_word(residue, length) == word,
                    f"theta inverse drift at word={word}")
            require(native not in native_images,
                    f"theta level collision at length={length}")
            require(len(history) == length + 1,
                    f"theta history length drift at word={word}")
            native_images.add(native)

            # Every native prefix extends to a positive eventually-zero word.
            completion = residue if residue > 0 else 2**length
            require(recover_carry_word(completion, length) == word,
                    f"positive terminal extension drift at word={word}")

            # k is the remaining ordinary completion; bit_length is the exact
            # number of further radix divisions until its owner mask empties.
            for k in (0, 1, 2, 3, 7, 8, 15, 16):
                if residue == 0 and k == 0:
                    continue
                quotient = k
                height = quotient.bit_length()
                while height > 0:
                    digit = quotient & 1
                    next_quotient = (quotient - digit) // 2
                    require(next_quotient.bit_length() == height - 1,
                            f"radix clock drift at length={length}, k={k}")
                    quotient = next_quotient
                    height -= 1
                    clock_gates += 1
                require(quotient == 0,
                        f"radix clock terminal drift at length={length}, k={k}")

            if direct:
                safe_count += 1
                safe_states.add(state)
                # Appending only zeroes leaves all old truncated sums fixed
                # and gives strict zero suffixes, so this is a strict lift.
                require(all(truncated_tail(word[start:]) < 1
                            for start in range(length)),
                        f"strict zero extension drift at word={word}")

        require(len(native_images) == 2**length,
                f"theta level surjectivity drift at length={length}")
        require(safe_states == set(range(length + 1)),
                f"follower state scout drift at length={length}")
        safe_counts.append(safe_count)
        state_counts.append(len(safe_states))
        pair_graph_counts.append(safe_count)

    require(total_words == 131070, "exhaustive word count drift")
    require(
        tuple(safe_counts)
        == (1, 2, 3, 5, 8, 12, 18, 27, 40, 60, 90, 134,
            201, 302, 452, 678, 1018),
        "safe-count sequence drift",
    )
    return (
        total_words,
        tuple(safe_counts),
        tuple(state_counts),
        tuple(pair_graph_counts),
        clock_gates,
    )


def audit_greedy_boundary(greedy):
    for length in range(1, 129):
        prefix = greedy[:length]
        require(direct_finite_safe(prefix),
                f"greedy prefix left P_m at length={length}")
        accepted, state, resets, rejection = follower_run(prefix, greedy)
        require(accepted and rejection is None and not resets and state == length,
                f"greedy follower boundary drift at length={length}")
        require(truncated_tail(prefix) < 1,
                f"finite greedy prefix ceased to be strict at length={length}")
        require(
            truncated_tail(prefix)
            + Fraction(2, 3) ** length
            * truncated_tail(greedy[length:GREEDY_AUDIT_DEPTH])
            < 1,
            f"finite greedy decomposition overshot at length={length}",
        )
    return "".join(map(str, greedy[:64]))


def audit_periodic_safe(greedy):
    block = (1, 0, 0)
    tails = tuple(periodic_tail(block, phase) for phase in range(3))
    phases = tuple(tail / 2 for tail in tails)
    require(tails == (Fraction(18, 19), Fraction(8, 19), Fraction(12, 19)),
            "(100)^infinity tail drift")
    require(phases == (Fraction(9, 19), Fraction(4, 19), Fraction(6, 19)),
            "(100)^infinity oriented-phase drift")
    word = block * 6
    native, _ = carry_to_native(word)
    native_text = "".join(map(str, native))
    require(native_text == PERIODIC_NATIVE_BLOCK,
            "denominator-19 native block drift")
    accepted, _, resets, _ = follower_run(block * 20, greedy)
    require(accepted and resets == tuple(range(2, 60, 3)),
            "periodic reset schedule drift")
    require(maximum_zero_run_periodic(native) == 4,
            "periodic native plateau bound drift")
    for length in range(1, 73):
        formal = -Fraction(9, 19)
        residue = (formal.numerator * pow(formal.denominator, -1, 2**length)) % 2**length
        require(canonical_residue((block * 24)[:length]) == residue,
                f"periodic formal residue drift at length={length}")
    return tails, phases, native_text


def audit_terminal_unsafe(greedy):
    carry_word = recover_carry_word(1, 16)
    require(carry_word[:4] == (1, 0, 1, 1), "A=1 carry prefix drift")
    require(binary_word(1, 16) == (1,) + (0,) * 15,
            "A=1 native terminal word drift")
    accepted3, _, _, rejection3 = follower_run(carry_word[:3], greedy)
    accepted4, _, _, rejection4 = follower_run(carry_word[:4], greedy)
    require(accepted3 and rejection3 is None,
            "A=1 rejected before the hostile boundary")
    require(not accepted4 and rejection4 == 3,
            "A=1 first rejection drift")
    partial = truncated_tail(carry_word[:4])
    require(partial == Fraction(94, 81) > 1,
            "A=1 tail obstruction drift")
    mapped, _ = carry_to_native(carry_word)
    require(mapped == binary_word(1, 16), "A=1 theta drift")
    return carry_word[:4], partial


def audit_centered_trap(greedy):
    word = (0, 1) * 8
    require(periodic_tail((0, 1), 0) == Fraction(4, 5),
            "(01)^infinity even tail drift")
    require(periodic_tail((0, 1), 1) == Fraction(6, 5),
            "(01)^infinity odd tail drift")
    formal = -Fraction(2, 5)
    residue = (formal.numerator * pow(formal.denominator, -1, 2**16)) % 2**16
    native_text = "".join(map(str, binary_word(residue, 16)))
    require(native_text == "0110" * 4, "(01)^infinity native block drift")
    accepted5, _, _, rejection5 = follower_run(word[:5], greedy)
    accepted6, _, _, rejection6 = follower_run(word[:6], greedy)
    require(accepted5 and rejection5 is None,
            "(01)^infinity rejected too early")
    require(not accepted6 and rejection6 == 5,
            "(01)^infinity first rejection drift")
    hostile_suffix = truncated_tail(word[1:6])
    require(hostile_suffix == Fraction(266, 243) > 1,
            "(01)^infinity truncated hostile drift")
    return native_text[:4], hostile_suffix


def audit_denominator_19_tower(greedy):
    block = (1, 0, 0)
    periodic_native = tuple(map(int, PERIODIC_NATIVE_BLOCK))
    rows = []
    for horizon in (18, 36, 54, 72):
        require(pow(3, horizon, 38) == 1,
                f"mod-38 tower phase drift at horizon={horizon}")
        integer_start = 9 * (2**horizon - 1) // 19
        require(19 * integer_start == 9 * (2**horizon - 1),
                f"tower integrality drift at horizon={horizon}")
        formal_word = (block * ((horizon + 2) // 3))[:horizon]
        require(recover_carry_word(integer_start, horizon) == formal_word,
                f"tower carry-prefix drift at horizon={horizon}")
        native_prefix = binary_word(integer_start, horizon)
        expected_native = (
            periodic_native * ((horizon + 17) // 18)
        )[:horizon]
        require(native_prefix == expected_native,
                f"tower joint native-prefix drift at horizon={horizon}")
        require(direct_finite_safe(formal_word),
                f"tower formal prefix left P_m at horizon={horizon}")
        alpha = Fraction(9 * 2**horizon, 19)
        phases = tuple(
            fractional_part(alpha * Fraction(3, 2) ** index)
            for index in range(horizon + 2)
        )
        require(all(phase < Fraction(1, 2) for phase in phases[: horizon + 1]),
                f"tower failed before horizon={horizon}")
        require(phases[horizon] == Fraction(9, 19),
                f"tower terminal safe phase drift at horizon={horizon}")
        require(phases[horizon + 1] == Fraction(27, 38) > Fraction(1, 2),
                f"tower first hostile phase drift at horizon={horizon}")
        require(binary_word(integer_start, horizon + 18)[horizon:] == (0,) * 18,
                f"tower ordinary terminal tail drift at horizon={horizon}")
        rows.append((horizon, integer_start.bit_length(), phases[horizon],
                     phases[horizon + 1]))
    return tuple(rows)


def main():
    greedy = greedy_boundary(GREEDY_AUDIT_DEPTH + MAX_EXHAUSTIVE_DEPTH + 4)
    (
        total_words,
        safe_counts,
        state_counts,
        pair_graph_counts,
        clock_gates,
    ) = audit_exhaustive(greedy)
    greedy_prefix = audit_greedy_boundary(greedy)
    periodic_tails, periodic_phases, native_block = audit_periodic_safe(greedy)
    a1_prefix, a1_partial = audit_terminal_unsafe(greedy)
    centered_native, centered_hostile = audit_centered_trap(greedy)
    tower_rows = audit_denominator_19_tower(greedy)

    print("THM-4072 primary exact fibre-product audit")
    print(f"exhaustive_depth={MAX_EXHAUSTIVE_DEPTH}")
    print(f"nonempty_binary_words_checked={total_words}")
    print("safe_counts_0_to_16=" + ",".join(map(str, safe_counts)))
    print("follower_states_0_to_16=" + ",".join(map(str, state_counts)))
    print("naive_pair_nodes_0_to_16=" + ",".join(map(str, pair_graph_counts)))
    print(f"radix_clock_descent_gates={clock_gates}")
    print(f"greedy_equality_prefix_64={greedy_prefix}")
    print("periodic_100_tails=" + ",".join(map(str, periodic_tails)))
    print("periodic_100_oriented_phases=" + ",".join(map(str, periodic_phases)))
    print(f"periodic_100_native_block_18={native_block}")
    print(f"terminal_A1_carry_prefix={''.join(map(str, a1_prefix))}")
    print(f"terminal_A1_truncated_tail={a1_partial}")
    print(f"centered_01_native_block_4={centered_native}")
    print(f"centered_01_hostile_suffix={centered_hostile}")
    print("tower_horizon_bitlength_safe_hostile=" + ";".join(
        f"{horizon}:{bits}:{safe}:{hostile}"
        for horizon, bits, safe, hostile in tower_rows
    ))
    print("result=PASS")


if __name__ == "__main__":
    main()
