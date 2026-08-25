#!/usr/bin/env python3
"""Primary exact audit for THM-4074's post-terminal delay theorem."""

from itertools import product


MAX_S = 6
MAX_H = 12
GREEDY_CONTROL_HEIGHTS = (1, 4, 8, 12)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuation_two(value):
    require(value != 0, "v2(0) is undefined")
    value = abs(value)
    exponent = 0
    while value % 2 == 0:
        value //= 2
        exponent += 1
    return exponent


def greedy_boundary(length):
    """Generate d from its odd-numerator dyadic recursion."""
    numerator = 1
    digits = []
    for state in range(length):
        digit = int(3 * numerator >= 2 ** (state + 1))
        digits.append(digit)
        numerator = 3 * numerator - digit * 2 ** (state + 1)
        require(0 < numerator < 2 ** (state + 1),
                f"greedy state escaped at depth {state + 1}")
        require(numerator % 2 == 1,
                f"greedy numerator lost oddness at depth {state + 1}")
    return tuple(digits)


def follower(word, boundary):
    state = 0
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


def integer_itinerary(start, length):
    value = start
    digits = []
    for _ in range(length):
        carry = value % 2
        digits.append(carry)
        value = (3 * value + carry) // 2
    return tuple(digits), value


def trace_to_rejection(state, value, boundary, limit):
    resets = 0
    maximum_state = state
    for index in range(limit):
        carry = value % 2
        wanted = boundary[state]
        if carry > wanted:
            return index, resets, maximum_state, state
        if carry < wanted:
            state = 0
            resets += 1
        else:
            state += 1
            maximum_state = max(maximum_state, state)
        value = (3 * value + carry) // 2
    raise RuntimeError(f"control did not reject within {limit} steps")


def layer_map(scale, height):
    """Map odd t modulo 2^h to the normalized odd post-zero state."""
    modulus = 2 ** height
    zero_count = scale + 3
    base_exponent = 18 * 2 ** scale
    base = 3 ** base_exponent
    divisor = 19 * 2 ** zero_count
    require((base - 1) % divisor == 0,
            f"normalizing divisor failed at s={scale}, h={height}")
    unit = 9 * 3 ** zero_count * ((base - 1) // divisor)
    require(unit % 2 == 1, f"normalizing unit is even at s={scale}")

    image = {}
    bracket = 0
    power = 1
    reduced_base = base % modulus
    for parameter in range(modulus):
        if parameter % 2 == 1:
            residue = (unit * bracket) % modulus
            require(residue % 2 == 1,
                    f"odd parameter reached even residue at s={scale}, h={height}")
            require(residue not in image,
                    f"collision at s={scale}, h={height}, residue={residue}")
            image[residue] = parameter
        bracket = (bracket + power) % modulus
        power = (power * reduced_base) % modulus

    odd_residues = set(range(1, modulus, 2))
    require(set(image) == odd_residues,
            f"odd-residue image incomplete at s={scale}, h={height}")
    return image


def safe_odd_start_words(height, boundary):
    words = []
    for suffix in product((0, 1), repeat=height - 1):
        word = (1,) + suffix
        accepted, _, _, _ = follower(word, boundary)
        if accepted:
            words.append(word)
    return tuple(words)


def launch(scale, parameter):
    require(parameter > 0 and parameter % 2 == 1,
            "launch parameter must be positive and odd")
    depth = 18 * 2 ** scale * parameter
    start = 9 * (2 ** depth - 1) // 19
    state_before_depth = 6 * (3 ** depth - 1) // 19
    state_at_depth = 9 * (3 ** depth - 1) // 19
    return depth, start, state_before_depth, state_at_depth


def main():
    boundary = greedy_boundary(256)
    require(boundary[:3] == (1, 0, 1), "denominator-19 loop lost its reset")

    layer_count = 0
    odd_parameter_classes = 0
    safe_realizations = 0
    safe_counts = []
    greedy_controls = []

    for scale in range(MAX_S + 1):
        scale_safe_counts = []
        for height in range(1, MAX_H + 1):
            image = layer_map(scale, height)
            layer_count += 1
            odd_parameter_classes += len(image)

            word_to_parameter = {}
            for residue, parameter in image.items():
                word, _ = integer_itinerary(residue, height)
                require(word[0] == 1,
                        f"odd residue has wrong first carry at h={height}")
                require(word not in word_to_parameter,
                        f"itinerary collision at s={scale}, h={height}")
                word_to_parameter[word] = parameter
            require(len(word_to_parameter) == 2 ** (height - 1),
                    f"odd itinerary image incomplete at s={scale}, h={height}")

            safe_words = safe_odd_start_words(height, boundary)
            scale_safe_counts.append(len(safe_words))
            for word in safe_words:
                require(word in word_to_parameter,
                        f"safe word missing at s={scale}, h={height}")
                safe_realizations += 1

            if height in GREEDY_CONTROL_HEIGHTS:
                target = boundary[:height]
                parameter = word_to_parameter[target]
                depth, start, _, state_at_depth = launch(scale, parameter)
                require(start.bit_length() == depth - 1,
                        f"native terminal depth failed at s={scale}, h={height}")
                zero_count = scale + 3
                require(valuation_two(state_at_depth) == zero_count,
                        f"post-reset valuation failed at s={scale}, h={height}")
                odd_state = 3 ** zero_count * state_at_depth // 2 ** zero_count
                word, _ = integer_itinerary(odd_state, height)
                require(word == target,
                        f"greedy target failed at s={scale}, h={height}")
                accepted, end_state, resets, _ = follower(
                    (0,) * zero_count + word, boundary
                )
                require(accepted and end_state == height and resets == zero_count,
                        f"greedy follower control failed at s={scale}, h={height}")
                greedy_controls.append((scale, height, parameter, depth))

        if scale == 0:
            safe_counts = scale_safe_counts
        else:
            require(scale_safe_counts == safe_counts,
                    f"safe counts depend on scale={scale}")

    require(layer_count == 84, "unexpected isometry-layer count")
    require(odd_parameter_classes == 28665,
            "unexpected odd-parameter class count")
    require(safe_realizations == 1400,
            "unexpected safe-target realization count")
    require(safe_counts == [1, 1, 2, 3, 4, 6, 9, 13, 20, 30, 44, 67],
            "safe odd-start count drift")

    direct_prefix_steps = 0
    rejection_controls = []
    for scale in range(MAX_S + 1):
        depth, start, state_before_depth, state_at_depth = launch(scale, 1)
        word, reached = integer_itinerary(start, depth)
        require(word == (1, 0, 0) * (depth // 3),
                f"direct denominator-19 prefix failed at s={scale}")
        require(reached == state_at_depth,
                f"closed orbit state failed at s={scale}")
        _, reached_before = integer_itinerary(start, depth - 1)
        require(reached_before == state_before_depth,
                f"immediate terminal state failed at s={scale}")
        require(start.bit_length() == depth - 1,
                f"last native one location failed at s={scale}")
        accepted, end_state, _, _ = follower(word, boundary)
        require(accepted and end_state == 0,
                f"periodic safe prefix failed at s={scale}")

        immediate_zero_count = scale + 4
        require(valuation_two(state_before_depth) == immediate_zero_count,
                f"immediate zero-run valuation failed at s={scale}")
        zeros, odd_state = integer_itinerary(state_before_depth,
                                             immediate_zero_count)
        require(zeros == (0,) * immediate_zero_count and odd_state % 2 == 1,
                f"immediate reset run failed at s={scale}")
        # The first zero is read at q=2; all later zeros are read at q=0.
        q = 2
        reset_count = 0
        for carry in zeros:
            require(carry < boundary[q],
                    f"zero stopped being a reset at s={scale}, q={q}")
            q = 0
            reset_count += 1
        require(q == 0 and reset_count == immediate_zero_count,
                f"reset count failed at s={scale}")

        rejection = trace_to_rejection(2, state_before_depth, boundary, 256)
        rejection_controls.append((scale,) + rejection[:3])
        direct_prefix_steps += depth

    require(direct_prefix_steps == 2286, "direct-prefix step count drift")

    even_parameter_controls = []
    for scale in range(MAX_S + 1):
        depth = 18 * 2 ** scale * 2
        state_at_depth = 9 * (3 ** depth - 1) // 19
        even_parameter_controls.append(valuation_two(state_at_depth))
        require(valuation_two(state_at_depth) == scale + 4,
                f"even-parameter valuation hostile failed at s={scale}")

    unsafe_word = (1, 1)
    require(integer_itinerary(3, 2)[0] == unsafe_word,
            "unsafe odd cylinder was not realized")
    accepted, _, _, reject_index = follower(unsafe_word, boundary)
    require(not accepted and reject_index == 1,
            "unsafe-cylinder hostile did not reject")
    require(trace_to_rejection(0, 1, boundary, 32)[0] == 3,
            "A=1 rejection hostile drift")

    print("THM-4074 primary exact audit")
    print(f"greedy_prefix_30={''.join(map(str, boundary[:30]))}")
    print(f"parameter_box=s:0..{MAX_S};h:1..{MAX_H};odd_t_mod_2^h")
    print(f"isometry_layers={layer_count}")
    print(f"odd_parameter_classes={odd_parameter_classes}")
    print(f"odd_residue_classes={odd_parameter_classes}")
    print(f"safe_odd_start_counts={','.join(map(str, safe_counts))}")
    print(f"safe_target_realizations={safe_realizations}")
    print(f"direct_denominator19_prefix_steps={direct_prefix_steps}")
    print("greedy_controls_h4_h8_h12=" + ";".join(
        f"s{s}:h{h}:t{t}:m{m}"
        for s, h, t, m in greedy_controls if h in (4, 8, 12)
    ))
    print("t1_rejection_controls=" + ";".join(
        f"s{s}:reject{reject}:resets{resets}:maxq{maxq}"
        for s, reject, resets, maxq in rejection_controls
    ))
    print("even_t2_v2=" + ",".join(map(str, even_parameter_controls)))
    print("hostiles=even_t_changes_zero_count;word_11_realized_but_rejects;A1_rejects_at_index3")
    print("scope=finite_postterminal_programming_only;no_infinite_survivor_claim")
    print("status=PASS")


if __name__ == "__main__":
    main()
