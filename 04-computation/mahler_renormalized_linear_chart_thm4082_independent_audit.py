#!/usr/bin/env python3
"""Independent analytic audit for THM-4082.

The load-bearing path constructs the 2-adic logarithm and exponential term by
term, including the even denominators through valuation/odd-part splitting.
It does not use the primary audit's cross-scale secant recurrence.  A direct
closed-form geometric sum is retained only as an independent comparison with
THM-4077's original map.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path


PRECISION = 72
INVERSE_PRECISION = 48
S_MAX = 12
R_MAX = 16
SIGNED_UNITS = (1, -1, 3, -3, 5, -5, 7, -7)
INVERSE_TARGET_LABELS = (
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


def valuation_two_mod(value: int, precision: int) -> int:
    modulus = 1 << precision
    value %= modulus
    if value == 0:
        return precision
    return (value & -value).bit_length() - 1


def logarithm_one_plus(argument: int, precision: int) -> tuple[int, int]:
    """Compute log(1+argument) modulo 2^precision from its defining series."""
    require(valuation_two_integer(argument) >= 2, "logarithm argument outside 1+4Z2")
    modulus = 1 << precision
    argument_valuation = valuation_two_integer(argument)
    total = 0
    retained_terms = 0
    for index in range(1, precision + 1):
        denominator_valuation = valuation_two_integer(index)
        term_valuation = index * argument_valuation - denominator_valuation
        if term_valuation >= precision:
            continue
        odd_denominator = index >> denominator_valuation
        lifted_modulus = modulus << denominator_valuation
        numerator = pow(argument, index, lifted_modulus)
        require(
            numerator % (1 << denominator_valuation) == 0,
            "logarithm term did not clear its dyadic denominator",
        )
        term = (
            (numerator >> denominator_valuation)
            * pow(odd_denominator, -1, modulus)
        ) % modulus
        total = (total + term) % modulus if index & 1 else (total - term) % modulus
        retained_terms += 1
    return total, retained_terms


def exponential(argument: int, precision: int) -> tuple[int, int]:
    """Compute exp(argument) modulo 2^precision from its defining series."""
    require(valuation_two_integer(argument) >= 2, "exponential argument outside 4Z2")
    modulus = 1 << precision
    argument_valuation = valuation_two_integer(argument)
    total = 1 % modulus
    factorial = 1
    retained_terms = 0
    for index in range(1, precision + 1):
        factorial *= index
        denominator_valuation = valuation_two_integer(factorial)
        term_valuation = index * argument_valuation - denominator_valuation
        if term_valuation >= precision:
            continue
        odd_denominator = factorial >> denominator_valuation
        lifted_modulus = modulus << denominator_valuation
        numerator = pow(argument, index, lifted_modulus)
        require(
            numerator % (1 << denominator_valuation) == 0,
            "exponential term did not clear its dyadic denominator",
        )
        term = (
            (numerator >> denominator_valuation)
            * pow(odd_denominator, -1, modulus)
        ) % modulus
        total = (total + term) % modulus
        retained_terms += 1
    return total, retained_terms


MAX_LOG_PRECISION = PRECISION + S_MAX + 6
ELL_MOD, ELL_TERMS = logarithm_one_plus(3**18 - 1, MAX_LOG_PRECISION)
LOG_MINUS_THREE_MOD, LOG_MINUS_THREE_TERMS = logarithm_one_plus(
    -4, MAX_LOG_PRECISION
)


def lambda_mod(precision: int) -> int:
    require(MAX_LOG_PRECISION >= precision + 3, "insufficient logarithm precision")
    modulus = 1 << precision
    ell = ELL_MOD % (1 << (precision + 3))
    require(ell % 8 == 0, "ell did not clear the denominator eight")
    return 243 * pow(19, -1, modulus) * (ell // 8) % modulus


def analytic_h(scale: int, value: int, precision: int) -> tuple[int, int, int]:
    require(0 <= scale <= S_MAX, "scale outside audited range")
    output_modulus = 1 << precision
    working_precision = precision + scale + 3
    working_modulus = 1 << working_precision
    ell = ELL_MOD % working_modulus
    scaled_exponent = (
        pow(2, scale, working_modulus)
        * pow(pow(3, scale, working_modulus), -1, working_modulus)
        * (value % working_modulus)
        * ell
    ) % working_modulus
    if value == 0:
        return 0, 0, scaled_exponent
    exponential_value, retained_terms = exponential(scaled_exponent, working_precision)
    numerator = (exponential_value - 1) % working_modulus
    dyadic_denominator = 1 << (scale + 3)
    require(numerator % dyadic_denominator == 0, "analytic H denominator did not clear")
    odd_coefficient = (
        9 * pow(3, scale + 3, output_modulus) * pow(19, -1, output_modulus)
    ) % output_modulus
    result = odd_coefficient * (numerator // dyadic_denominator) % output_modulus
    return result, retained_terms, scaled_exponent


def geometric_sum_mod(base: int, exponent: int, modulus: int) -> int:
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


EXACT_PARAMETERS: list[tuple[int, int]] = []
for SCALE in range(S_MAX + 1):
    LENGTH = 18 * 2**SCALE
    K_VALUE = SCALE + 3
    G_VALUE = 3**LENGTH
    NUMERATOR = 9 * 3**K_VALUE * (G_VALUE - 1)
    DENOMINATOR = 19 * 2**K_VALUE
    require(NUMERATOR % DENOMINATOR == 0, "independent exact U integrality failed")
    EXACT_PARAMETERS.append((G_VALUE, NUMERATOR // DENOMINATOR))


def direct_original_h(scale: int, value: int, precision: int) -> int:
    modulus = 1 << precision
    g_value, u_value = EXACT_PARAMETERS[scale]
    parameter = pow(pow(3, scale, modulus), -1, modulus) * (value % modulus) % modulus
    return u_value * geometric_sum_mod(g_value, parameter, modulus) % modulus


def inverse_analytic_isometry(scale: int, target: int) -> int:
    representative = 0
    for height in range(1, INVERSE_PRECISION + 1):
        modulus = 1 << height
        upper = representative + (1 << (height - 1))
        wanted = target % modulus
        lower_value = analytic_h(scale, representative, INVERSE_PRECISION)[0] % modulus
        upper_value = analytic_h(scale, upper, INVERSE_PRECISION)[0] % modulus
        lower_good = lower_value == wanted
        upper_good = upper_value == wanted
        require(lower_good != upper_good, "analytic inverse bit was not unique")
        if upper_good:
            representative = upper
    return representative


def inverse_targets(precision: int) -> dict[str, int]:
    modulus = 1 << precision
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
    require(
        (ELL_MOD - 18 * LOG_MINUS_THREE_MOD) % (1 << MAX_LOG_PRECISION) == 0,
        "log(3^18)=18 log(-3) normalization failed",
    )
    require(
        valuation_two_mod(ELL_MOD, MAX_LOG_PRECISION) == 3,
        "ell valuation changed",
    )
    limiting_multiplier = lambda_mod(PRECISION)
    require(limiting_multiplier & 1 == 1, "analytic tangent is not a unit")

    analytic_valuation_gates = 0
    original_map_comparison_gates = 0
    exponential_remainder_gates = 0
    maximum_exponential_terms = 0
    valuation_rows: list[tuple[int, int, int, int]] = []
    for scale in range(S_MAX + 1):
        for input_valuation in range(R_MAX + 1):
            for unit in SIGNED_UNITS:
                value = (1 << input_valuation) * unit
                analytic_value, retained_terms, exponent = analytic_h(
                    scale, value, PRECISION
                )
                direct_value = direct_original_h(scale, value, PRECISION)
                require(analytic_value == direct_value, "analytic/original map mismatch")
                expected = scale + 2 + 2 * input_valuation
                observed = valuation_two_mod(
                    analytic_value - limiting_multiplier * (value % modulus),
                    PRECISION,
                )
                require(observed == expected, "independent exact tangent valuation failed")

                working_precision = PRECISION + scale + 3
                working_modulus = 1 << working_precision
                exponential_value, _ = exponential(exponent, working_precision)
                remainder = (exponential_value - 1 - exponent) % working_modulus
                expected_remainder = 2 * valuation_two_mod(exponent, working_precision) - 1
                require(
                    valuation_two_mod(remainder, working_precision) == expected_remainder,
                    "exponential quadratic term was not uniquely minimal",
                )
                analytic_valuation_gates += 1
                original_map_comparison_gates += 1
                exponential_remainder_gates += 1
                maximum_exponential_terms = max(maximum_exponential_terms, retained_terms)
                if unit == 1:
                    valuation_rows.append((scale, input_valuation, expected, observed))

        zero_value, _, _ = analytic_h(scale, 0, PRECISION)
        require(zero_value == 0, "analytic zero boundary moved")

    inverse_modulus = 1 << INVERSE_PRECISION
    inverse_limiting_multiplier = lambda_mod(INVERSE_PRECISION)
    inverse_linear = pow(inverse_limiting_multiplier, -1, inverse_modulus)
    target_values = inverse_targets(INVERSE_PRECISION)
    inverse_fibre_gates = 0
    inverse_zero_gates = 0
    inverse_rows: list[tuple[int, str, int]] = []
    for scale in range(S_MAX + 1):
        for label in INVERSE_TARGET_LABELS:
            target = target_values[label]
            nonlinear_preimage = inverse_analytic_isometry(scale, target)
            linear_preimage = inverse_linear * target % inverse_modulus
            target_valuation = valuation_two_mod(target, INVERSE_PRECISION)
            expected = scale + 2 + 2 * target_valuation
            observed = valuation_two_mod(
                nonlinear_preimage - linear_preimage, INVERSE_PRECISION
            )
            require(observed == expected, "analytic inverse-fibre valuation failed")
            require(
                analytic_h(scale, nonlinear_preimage, INVERSE_PRECISION)[0] == target,
                "analytic inverse fibre missed its target",
            )
            inverse_rows.append((scale, label, observed))
            inverse_fibre_gates += 1
        require(inverse_analytic_isometry(scale, 0) == 0, "analytic zero inverse moved")
        inverse_zero_gates += 1

    periodic_target = target_values["periodic_100"]
    require(
        valuation_two_mod(1 - periodic_target, INVERSE_PRECISION) == 2,
        "safe/output target separation changed",
    )

    assert_nodes, float_literals = source_shape()
    require(assert_nodes == 0, "assert statement found in independent source")
    require(float_literals == 0, "floating literal found in independent source")

    digest_payload = {
        "ell_mod_65536": ELL_MOD % 65536,
        "lambda_mod_65536": limiting_multiplier % 65536,
        "valuation_rows": valuation_rows,
        "inverse_rows": inverse_rows,
        "max_exp_terms": maximum_exponential_terms,
    }
    semantic_digest = hashlib.sha256(
        json.dumps(digest_payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("THM-4082 independent logarithm/exponential audit")
    print(
        f"precision_bits={PRECISION} inverse_precision_bits={INVERSE_PRECISION} "
        f"s_range=0..{S_MAX} valuation_range=0..{R_MAX}"
    )
    print(
        f"log_terms_ell={ELL_TERMS} log_terms_minus_three={LOG_MINUS_THREE_TERMS} "
        f"ell_v2={valuation_two_mod(ELL_MOD, MAX_LOG_PRECISION)}"
    )
    print(
        f"analytic_valuation_gates={analytic_valuation_gates} "
        f"original_map_comparison_gates={original_map_comparison_gates} "
        f"exponential_remainder_gates={exponential_remainder_gates}"
    )
    print(
        f"inverse_fibre_gates={inverse_fibre_gates} "
        f"inverse_zero_gates={inverse_zero_gates}"
    )
    print(f"maximum_exponential_terms={maximum_exponential_terms}")
    print(
        f"ell_mod_65536={ELL_MOD % 65536} "
        f"lambda_mod_65536={limiting_multiplier % 65536}"
    )
    print(f"source_assert_nodes={assert_nodes} source_float_literals={float_literals}")
    print(f"semantic_sha256={semantic_digest}")
    print("PASS: independent analytic normal form, exact remainder, and inverse fibres")


if __name__ == "__main__":
    main()
