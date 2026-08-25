#!/usr/bin/env python3
"""Primary exact audit for THM-4082.

This path deliberately avoids 2-adic logarithm and exponential series.  It
constructs the renormalized maps from the affine/geometric recurrence, obtains
the limiting tangent from the cross-scale secants, and checks the exact bit
defect, carry first divergence, safe-prefix transport, inverse fibres, and the
ordinary vertical divisibility sidecar.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path


PRECISION = 48
S_MAX = 12
X_MAX = 4095
STABILIZATION_STEPS = PRECISION + 4
INVERSE_TARGET_NAMES = (
    "one",
    "minus_one",
    "periodic_100",
    "three",
    "five",
    "seven",
    "eleven",
    "two",
    "four",
    "twelve",
    "forty",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation_two_integer(value: int) -> int:
    require(value != 0, "valuation of zero requested")
    value = abs(value)
    return (value & -value).bit_length() - 1


def valuation_two_mod(value: int, precision: int = PRECISION) -> int:
    modulus = 1 << precision
    value %= modulus
    if value == 0:
        return precision
    return (value & -value).bit_length() - 1


def geometric_sum_mod(base: int, exponent: int, modulus: int) -> int:
    """Return 1+base+...+base^(exponent-1) modulo modulus."""
    require(exponent >= 0, "negative geometric exponent")
    power = 1 % modulus
    total = 0
    block_power = base % modulus
    block_sum = 1 % modulus
    remaining = exponent
    while remaining:
        if remaining & 1:
            total = (total + power * block_sum) % modulus
            power = power * block_power % modulus
        block_sum = block_sum * (1 + block_power) % modulus
        block_power = block_power * block_power % modulus
        remaining //= 2
    return total


def build_cross_scale_rows() -> list[tuple[int, int]]:
    modulus = 1 << PRECISION
    lifted_modulus = 1 << (PRECISION + 1)
    g_value = pow(3, 18, lifted_modulus)
    g_zero = 3**18
    u_zero = 9 * 3**3 * (g_zero - 1) // (19 * 2**3)
    v_value = u_zero % modulus
    rows: list[tuple[int, int]] = []
    for _ in range(STABILIZATION_STEPS + 1):
        rows.append((g_value % modulus, v_value))
        factor = ((g_value + 1) // 2) % modulus
        v_value = v_value * factor % modulus
        g_value = g_value * g_value % lifted_modulus
    return rows


CROSS_SCALE_ROWS = build_cross_scale_rows()
LAMBDA_MOD = CROSS_SCALE_ROWS[STABILIZATION_STEPS][1]


def h_mod(scale: int, value: int, precision: int = PRECISION) -> int:
    require(0 <= scale <= S_MAX, "scale outside audited range")
    modulus = 1 << precision
    g_value, v_value = CROSS_SCALE_ROWS[scale]
    g_value %= modulus
    v_value %= modulus
    u_value = pow(3, scale, modulus) * v_value % modulus
    inverse_scale = pow(pow(3, scale, modulus), -1, modulus)
    parameter = inverse_scale * (value % modulus) % modulus
    return u_value * geometric_sum_mod(g_value, parameter, modulus) % modulus


def carry_prefix(state: int, height: int) -> tuple[int, ...]:
    digits: list[int] = []
    current = state
    for _ in range(height):
        digit = current & 1
        digits.append(digit)
        current = (3 * current + digit) // 2
    return tuple(digits)


def greedy_boundary(height: int) -> tuple[int, ...]:
    numerator = 1
    digits: list[int] = []
    for index in range(height):
        digit = int(3 * numerator >= 1 << (index + 1))
        digits.append(digit)
        numerator = 3 * numerator - digit * (1 << (index + 1))
        require(0 < numerator < 1 << (index + 1), "greedy state escaped")
        require(numerator & 1 == 1, "greedy numerator lost oddness")
    return tuple(digits)


GREEDY = greedy_boundary(PRECISION)


def safe_follower(word: tuple[int, ...]) -> bool:
    state = 0
    for digit in word:
        expected = GREEDY[state]
        if digit > expected:
            return False
        state = state + 1 if digit == expected else 0
    return True


def inverse_isometry(scale: int, target: int) -> int:
    representative = 0
    for height in range(1, PRECISION + 1):
        modulus = 1 << height
        upper = representative + (1 << (height - 1))
        wanted = target % modulus
        lower_good = h_mod(scale, representative) % modulus == wanted
        upper_good = h_mod(scale, upper) % modulus == wanted
        require(lower_good != upper_good, "one-bit inverse lift was not unique")
        if upper_good:
            representative = upper
    return representative


def inverse_targets() -> dict[str, int]:
    modulus = 1 << PRECISION
    return {
        "one": 1,
        "minus_one": (-1) % modulus,
        "periodic_100": (-9 * pow(19, -1, modulus)) % modulus,
        "three": 3,
        "five": 5,
        "seven": 7,
        "eleven": 11,
        "two": 2,
        "four": 4,
        "twelve": 12,
        "forty": 40,
    }


def exact_parameters(scale: int) -> tuple[int, int, int, int]:
    length = 18 * 2**scale
    k_value = scale + 3
    g_value = 3**length
    numerator = 9 * 3**k_value * (g_value - 1)
    denominator = 19 * 2**k_value
    require(numerator % denominator == 0, "exact U integrality failed")
    u_value = numerator // denominator
    return length, k_value, g_value, u_value


def source_shape() -> tuple[int, int]:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    return assert_nodes, float_literals


def main() -> None:
    modulus = 1 << PRECISION
    require(LAMBDA_MOD & 1 == 1, "limiting multiplier is not a unit")

    secant_rows: list[tuple[int, int, int]] = []
    secant_gates = 0
    for scale in range(S_MAX + 1):
        _, v_value = CROSS_SCALE_ROWS[scale]
        defect = valuation_two_mod(LAMBDA_MOD - v_value)
        require(defect == scale + 2, "secant convergence valuation failed")
        secant_rows.append((scale, v_value % 65536, defect))
        secant_gates += 1

    valuation_gates = 0
    carry_divergence_gates = 0
    for scale in range(S_MAX + 1):
        for value in range(1, X_MAX + 1):
            input_valuation = valuation_two_integer(value)
            expected = scale + 2 + 2 * input_valuation
            nonlinear = h_mod(scale, value)
            linear = LAMBDA_MOD * value % modulus
            observed = valuation_two_mod(nonlinear - linear)
            require(observed == expected, "sharp nonlinear defect failed")
            nonlinear_word = carry_prefix(nonlinear, expected + 1)
            linear_word = carry_prefix(linear, expected + 1)
            require(
                nonlinear_word[:expected] == linear_word[:expected],
                "carry common prefix was too short",
            )
            require(
                nonlinear_word[expected] != linear_word[expected],
                "carry first-divergence bit did not flip",
            )
            valuation_gates += 1
            carry_divergence_gates += 1
        require(h_mod(scale, 0) == 0, "zero boundary moved")

    safe_prefix_gates = 0
    safe_counts: list[int] = []
    for scale in range(S_MAX + 1):
        height = scale + 2
        local_modulus = 1 << height
        safe_count = 0
        for value in range(1, local_modulus, 2):
            nonlinear_word = carry_prefix(h_mod(scale, value) % local_modulus, height)
            linear_word = carry_prefix(LAMBDA_MOD * value % local_modulus, height)
            require(nonlinear_word == linear_word, "safe-prefix word transport failed")
            nonlinear_safe = safe_follower(nonlinear_word)
            linear_safe = safe_follower(linear_word)
            require(nonlinear_safe == linear_safe, "safe-prefix predicate drifted")
            safe_count += int(nonlinear_safe)
            safe_prefix_gates += 1
        safe_counts.append(safe_count)

    target_values = inverse_targets()
    inverse_fibre_gates = 0
    inverse_zero_gates = 0
    inverse_rows: list[tuple[int, str, int]] = []
    inverse_lambda = pow(LAMBDA_MOD, -1, modulus)
    for scale in range(S_MAX + 1):
        for name in INVERSE_TARGET_NAMES:
            target = target_values[name]
            nonlinear_preimage = inverse_isometry(scale, target)
            linear_preimage = inverse_lambda * target % modulus
            target_valuation = valuation_two_mod(target)
            expected = scale + 2 + 2 * target_valuation
            observed = valuation_two_mod(nonlinear_preimage - linear_preimage)
            require(observed == expected, "inverse-fibre convergence failed")
            require(h_mod(scale, nonlinear_preimage) == target, "inverse target missed")
            inverse_rows.append((scale, name, observed))
            inverse_fibre_gates += 1
        require(inverse_isometry(scale, 0) == 0, "zero inverse fibre moved")
        inverse_zero_gates += 1

    vertical_divisibility_gates = 0
    for scale in range(5):
        _, _, g_value, u_value = exact_parameters(scale)
        _, _, next_g, next_u = exact_parameters(scale + 1)
        require(next_g == g_value * g_value, "vertical base squaring failed")
        for parameter in range(1, 20, 2):
            current = u_value * (g_value**parameter - 1) // (g_value - 1)
            following = next_u * (next_g**parameter - 1) // (next_g - 1)
            multiplier = 3 * (g_value**parameter + 1) // 2
            require(multiplier & 1 == 1, "vertical multiplier lost oddness")
            require(following == multiplier * current, "vertical divisibility failed")
            vertical_divisibility_gates += 1

    periodic_target = target_values["periodic_100"]
    periodic_word = carry_prefix(periodic_target, 18)
    require(periodic_word == (1, 0, 0) * 6, "periodic safe hostile drifted")
    require(safe_follower(periodic_word), "periodic safe hostile rejected")
    one_word = carry_prefix(1, 4)
    require(one_word == (1, 0, 1, 1), "state-one hostile drifted")
    require(not safe_follower(one_word), "state-one hostile became safe")
    require(
        valuation_two_mod(1 - periodic_target) == 2,
        "safe/output hostile separation changed",
    )

    small_scale_controls = [
        (scale, value, valuation_two_mod(h_mod(scale, value) - LAMBDA_MOD * value))
        for scale in (0, 1)
        for value in (1, 2)
    ]
    require(
        small_scale_controls == [(0, 1, 2), (0, 2, 4), (1, 1, 3), (1, 2, 5)],
        "small-scale/even boundary controls changed",
    )

    assert_nodes, float_literals = source_shape()
    require(assert_nodes == 0, "assert statement found in source")
    require(float_literals == 0, "floating literal found in source")

    digest_payload = {
        "secants": secant_rows,
        "safe_counts": safe_counts,
        "inverse_rows": inverse_rows,
        "small_controls": small_scale_controls,
        "lambda_mod_65536": LAMBDA_MOD % 65536,
    }
    semantic_digest = hashlib.sha256(
        json.dumps(digest_payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("THM-4082 primary affine/geometric audit")
    print(f"precision_bits={PRECISION} s_range=0..{S_MAX} x_range=1..{X_MAX}")
    print(f"secant_gates={secant_gates} zero_forward_gates={S_MAX + 1}")
    print(
        f"valuation_gates={valuation_gates} "
        f"carry_first_divergence_gates={carry_divergence_gates}"
    )
    print(f"safe_triangular_gates={safe_prefix_gates} safe_counts={safe_counts}")
    print(
        f"inverse_fibre_gates={inverse_fibre_gates} "
        f"inverse_zero_gates={inverse_zero_gates}"
    )
    print(f"vertical_divisibility_gates={vertical_divisibility_gates}")
    print(f"lambda_mod_65536={LAMBDA_MOD % 65536} secant_rows={secant_rows}")
    print(f"small_scale_even_controls={small_scale_controls}")
    print(f"source_assert_nodes={assert_nodes} source_float_literals={float_literals}")
    print(f"semantic_sha256={semantic_digest}")
    print("PASS: sharp tangent, exact bit defect, fibre transport, and boundaries")


if __name__ == "__main__":
    main()
