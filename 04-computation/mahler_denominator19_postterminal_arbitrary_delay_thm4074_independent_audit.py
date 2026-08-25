#!/usr/bin/env python3
"""Independent native-first audit for THM-4074."""

from fractions import Fraction
from itertools import product


MAX_S = 6
MAX_H = 12
DIRECT_MAX_S = 10


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def greedy_boundary(length):
    state = Fraction(1)
    digits = []
    for index in range(length):
        scaled = Fraction(3, 2) * state
        digit = scaled.numerator // scaled.denominator
        require(digit in (0, 1), f"greedy digit escaped at index={index}")
        digits.append(digit)
        state = scaled - digit
        require(Fraction(0) < state < Fraction(1),
                f"greedy remainder escaped at index={index}")
    return tuple(digits)


def carry_integer(word):
    length = len(word)
    return sum(
        digit * 2 ** index * 3 ** (length - 1 - index)
        for index, digit in enumerate(word)
    )


def direct_truncated_safe(word):
    return all(
        2 * carry_integer(word[start:]) < 3 ** (len(word) - start)
        for start in range(len(word))
    )


def native_orbit(start, length):
    value = start
    carries = []
    for _ in range(length):
        carry = value % 2
        carries.append(carry)
        value = (3 * value + carry) // 2
    return tuple(carries), value


def follower(word, boundary, initial_state=0):
    state = initial_state
    resets = 0
    for index, carry in enumerate(word):
        wanted = boundary[state]
        if carry > wanted:
            return False, state, resets, index
        if carry < wanted:
            state = 0
            resets += 1
        else:
            state += 1
    return True, state, resets, len(word)


def normalized_residue(scale, parameter, height):
    """Compute W_s(t) mod 2^h directly from 3^m modulo 2^(h+k)."""
    zero_count = scale + 3
    depth = 18 * 2 ** scale * parameter
    large_modulus = 2 ** (height + zero_count)
    power_residue = pow(3, depth, large_modulus)
    require((power_residue - 1) % 2 ** zero_count == 0,
            f"direct modular division failed at s={scale}, t={parameter}")
    quotient = (power_residue - 1) // 2 ** zero_count
    modulus = 2 ** height
    return (
        9 * 3 ** zero_count * quotient * pow(19, -1, modulus)
    ) % modulus


def lift_parameter(scale, target):
    """Lift t one bit at a time, without building the geometric map."""
    require(target and target[0] == 1, "target must begin with one")
    parameter = 1
    lift_gates = 0
    for height in range(2, len(target) + 1):
        candidates = (parameter, parameter + 2 ** (height - 1))
        survivors = []
        for candidate in candidates:
            residue = normalized_residue(scale, candidate, height)
            word, _ = native_orbit(residue, height)
            lift_gates += 1
            if word == target[:height]:
                survivors.append(candidate)
        require(len(survivors) == 1,
                f"Hensel lift not unique at s={scale}, h={height}, target={target}")
        parameter = survivors[0]
    residue = normalized_residue(scale, parameter, len(target))
    word, _ = native_orbit(residue, len(target))
    require(word == target, f"final lifted target failed at s={scale}")
    return parameter, lift_gates


def trace_to_rejection(state, value, boundary, limit):
    resets = 0
    maximum_state = state
    for index in range(limit):
        carry = value % 2
        wanted = boundary[state]
        if carry > wanted:
            return index, resets, maximum_state
        if carry < wanted:
            state = 0
            resets += 1
        else:
            state += 1
            maximum_state = max(maximum_state, state)
        value = (3 * value + carry) // 2
    raise RuntimeError(f"native hostile survived {limit} steps")


def main():
    boundary = greedy_boundary(256)
    require(boundary[:3] == (1, 0, 1), "reset loop changed")

    target_count = 0
    lift_gate_count = 0
    safe_target_count = 0
    safe_counts = []
    greedy_rows = []

    for scale in range(MAX_S + 1):
        local_safe_counts = []
        for height in range(1, MAX_H + 1):
            used_parameters = set()
            local_safe = 0
            for suffix in product((0, 1), repeat=height - 1):
                target = (1,) + suffix
                parameter, gates = lift_parameter(scale, target)
                require(parameter % 2 == 1 and 0 < parameter < 2 ** height,
                        f"lifted parameter escaped at s={scale}, h={height}")
                require(parameter not in used_parameters,
                        f"two targets share a parameter at s={scale}, h={height}")
                used_parameters.add(parameter)
                target_count += 1
                lift_gate_count += gates

                safe = direct_truncated_safe(target)
                accepted, _, _, _ = follower(target, boundary)
                require(safe == accepted,
                        f"suffix/follower disagreement at h={height}, target={target}")
                if safe:
                    local_safe += 1
                    safe_target_count += 1

            require(used_parameters == set(range(1, 2 ** height, 2)),
                    f"Hensel parameters incomplete at s={scale}, h={height}")
            local_safe_counts.append(local_safe)

            if height in (4, 8, 12):
                parameter, _ = lift_parameter(scale, boundary[:height])
                zero_count = scale + 3
                residue = normalized_residue(scale, parameter, height)
                word, _ = native_orbit(residue, height)
                accepted, end_state, resets, _ = follower(
                    (0,) * zero_count + word, boundary
                )
                require(accepted and end_state == height and resets == zero_count,
                        f"greedy follower lift failed at s={scale}, h={height}")
                greedy_rows.append((scale, height, parameter))

        if scale == 0:
            safe_counts = local_safe_counts
        else:
            require(local_safe_counts == safe_counts,
                    f"direct-safe counts depend on scale={scale}")

    require(target_count == 28665, "odd target count drift")
    require(safe_target_count == 1400, "safe target count drift")
    require(safe_counts == [1, 1, 2, 3, 4, 6, 9, 13, 20, 30, 44, 67],
            "direct safe-count row drift")

    direct_steps = 0
    native_controls = []
    for scale in range(DIRECT_MAX_S + 1):
        depth = 18 * 2 ** scale
        start = 9 * (2 ** depth - 1) // 19
        require(19 * start == 9 * (2 ** depth - 1),
                f"start identity failed at s={scale}")
        require(start.bit_length() == depth - 1,
                f"native last-one depth failed at s={scale}")

        prefix, state_at_depth = native_orbit(start, depth)
        require(prefix == (1, 0, 0) * (depth // 3),
                f"native-first prefix failed at s={scale}")
        require(19 * state_at_depth == 9 * (3 ** depth - 1),
                f"native-first endpoint failed at s={scale}")
        accepted, end_state, _, _ = follower(prefix, boundary)
        require(accepted and end_state == 0,
                f"native-first safe prefix failed at s={scale}")

        _, immediate_state = native_orbit(start, depth - 1)
        require(19 * immediate_state == 6 * (3 ** depth - 1),
                f"immediate post-terminal identity failed at s={scale}")
        zero_count = scale + 4
        zeros, odd_state = native_orbit(immediate_state, zero_count)
        require(zeros == (0,) * zero_count and odd_state % 2 == 1,
                f"native reset run failed at s={scale}")
        accepted, end_state, resets, _ = follower(zeros, boundary, 2)
        require(accepted and end_state == 0 and resets == zero_count,
                f"native follower reset run failed at s={scale}")

        rejection, rejection_resets, maximum_state = trace_to_rejection(
            2, immediate_state, boundary, 512
        )
        native_controls.append(
            (scale, depth, zero_count, rejection, rejection_resets, maximum_state)
        )
        direct_steps += depth

    require(direct_steps == 36846, "native direct-step count drift")

    for scale in range(MAX_S + 1):
        even_parameter = 2
        depth = 18 * 2 ** scale * even_parameter
        state = 9 * (3 ** depth - 1) // 19
        zero_count = 0
        while state % 2 == 0:
            state //= 2
            zero_count += 1
        require(zero_count == scale + 4,
                f"even-parameter hostile failed at s={scale}")

    unsafe_parameter, _ = lift_parameter(0, (1, 1))
    unsafe_residue = normalized_residue(0, unsafe_parameter, 2)
    unsafe_word, _ = native_orbit(unsafe_residue, 2)
    require(unsafe_word == (1, 1) and not direct_truncated_safe(unsafe_word),
            "unsafe-cylinder hostile failed")

    print("THM-4074 independent native-first exact audit")
    print(f"greedy_prefix_30={''.join(map(str, boundary[:30]))}")
    print(f"hensel_box=s:0..{MAX_S};h:1..{MAX_H};all_odd_start_words")
    print(f"hensel_target_classes={target_count}")
    print(f"hensel_candidate_gates={lift_gate_count}")
    print(f"direct_suffix_safe_counts={','.join(map(str, safe_counts))}")
    print(f"safe_target_lifts={safe_target_count}")
    print(f"native_direct_prefix_steps={direct_steps}")
    print("greedy_lifts=" + ";".join(
        f"s{s}:h{h}:t{t}" for s, h, t in greedy_rows
    ))
    print("native_t1_controls=" + ";".join(
        f"s{s}:m{m}:zeros{zeros}:reject{reject}:resets{resets}:maxq{maxq}"
        for s, m, zeros, reject, resets, maxq in native_controls
    ))
    print("hostiles=even_t_has_one_extra_v2;word_11_lifts_but_is_not_safe")
    print("scope=finite_exact_layers_and_all_parameter_proof_identities;Mahler_remains_open")
    print("status=PASS")


if __name__ == "__main__":
    main()
