"""Exact hostile controls for the September 21 Collatz energy blueprint.

Universal results are proved in the companion note. Integer/Fraction checks
here cover declared finite universes; floating diagnostics are not proofs.
Uses only the Python standard library and retains checks under python -O.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import cos, floor, isqrt, log, log2, pi, sin, sqrt
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def step(n: int) -> tuple[int, int]:
    require(n > 0 and n % 2 == 1, "positive odd input required")
    value = 3 * n + 1
    k = 0
    while value % 2 == 0:
        value //= 2
        k += 1
    return value, k


def rational_step(n: Fraction) -> tuple[Fraction, int]:
    require(n.numerator % 2 == n.denominator % 2 == 1, "odd unit required")
    value = 3 * n + 1
    k = 0
    numerator = value.numerator
    while numerator % 2 == 0:
        numerator //= 2
        k += 1
    return Fraction(numerator, value.denominator), k


def energy(x: float) -> float:
    u = log2(x) - floor(log2(x))
    return log(x) * (1.5 + 0.5 * cos(pi * cos(pi * u / 2) ** 2))


def energy_derivative(x: float) -> float:
    u = log2(x) - floor(log2(x))
    a = pi * cos(pi * u / 2) ** 2
    profile = 1.5 + 0.5 * cos(a)
    profile_derivative = pi**2 * sin(a) * sin(pi * u) / 4
    return (profile + log2(x) * profile_derivative) / x


def factor(n: int) -> list[tuple[int, int]]:
    result = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            exponent += 1
            n //= p
        if exponent:
            result.append((p, exponent))
        p += 1
    if n > 1:
        result.append((n, 1))
    return result


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def run() -> dict:
    families = []
    step_count = 0
    for length in range(1, 65):
        c = 2 ** (length + 1) - 1
        for m in sorted({2 * length + 4, 2 * length + 5, 3 * length + 10}):
            start = 2**m + c
            current = start
            for j in range(1, length + 1):
                current, exponent = step(current)
                expected = 3**j * (2 ** (m - j) + 2 ** (length + 1 - j)) - 1
                require(exponent == 1, "growth family exponent failure")
                require(current == expected, "growth family affine formula failure")
                require(Fraction(current, start) > Fraction(3, 2) ** j,
                        "strict multiplicative growth failure")
                step_count += 1
            rational_error_bound = 512 * (m + 1) * Fraction(c, 2**m) ** 4
            require(rational_error_bound < Fraction(7, 128), "error bound failure")
            families.append({"L": length, "m": m, "start": start,
                             "last": current,
                             "rational_error_upper_bound": fraction_text(rational_error_bound)})

    same_band = []
    for m in range(65):
        start = 20 * 2**m - 1
        target, exponent = step(start)
        require(exponent == 1 and target == 30 * 2**m - 1, "same-band map failure")
        require(16 * 2**m < start < target < 32 * 2**m, "same-band containment failure")
        same_band.append({"m": m, "start": start, "target": target})

    supplied_gaps = [2, 4, 2, 4, 4, 2, 4, 2, 2]
    binary = [(x - 2) // 2 for x in supplied_gaps]
    sums = [sum(binary[j:j + 2]) for j in range(len(binary) - 1)]
    require(0 in sums and 2 in sums, "nonbalanced witness missing")
    wythoff = [(k + isqrt(5 * k * k)) // 2 for k in range(10)]
    true_gaps = [2 * (b - a) for a, b in zip(wythoff, wythoff[1:])]
    require(true_gaps == [2, 4, 2, 4, 4, 2, 4, 2, 4], "Wythoff calculation failure")

    rational_cycle = [Fraction(5, 7), Fraction(11, 7)]
    rational_exponents = []
    for index, value in enumerate(rational_cycle):
        target, exponent = rational_step(value)
        require(target == rational_cycle[(index + 1) % 2], "rational cycle failure")
        rational_exponents.append(exponent)
    require(rational_exponents == [1, 3], "rational valuations failure")

    squarefree_path = [3**j * 2 ** (11 - j) - 1 for j in range(11)]
    factorizations = []
    for index, n in enumerate(squarefree_path):
        factors = factor(n)
        require(all(exponent == 1 for _, exponent in factors), "squarefree witness failure")
        require(all(n % (d * d) for d in range(2, isqrt(n) + 1)),
                "independent trial-square failure")
        factorizations.append({"n": n, "factorization": factors})
        if index < 10:
            require(step(n) == (squarefree_path[index + 1], 1), "squarefree path failure")

    short_interval_witnesses = []
    for length in [1, 2, 3, 6, 10]:
        modulus = 2 ** (length + 1)
        m = max(16, 2 * length + 4)
        # n=Mq-1 in [2^m,2^m(1+1/m)], with integer endpoints exact.
        first_q = (2**m + 1 + modulus - 1) // modulus
        last_q = ((m + 1) * 2**m + m) // (m * modulus)
        witness = None
        checked = 0
        for q in range(first_q, last_q + 1):
            checked += 1
            nodes = [3**j * 2 ** (length + 1 - j) * q - 1
                     for j in range(length + 1)]
            if all(all(exponent == 1 for _, exponent in factor(n)) for n in nodes):
                witness = (q, nodes)
                break
        require(witness is not None, "no squarefree short-interval witness found")
        q, nodes = witness
        start = nodes[0]
        require(all(all(n % (d * d) for d in range(2, isqrt(n) + 1)) for n in nodes),
                "independent short-interval trial-square test failed")
        require(2**m <= start and m * start <= (m + 1) * 2**m,
                "short-interval containment failed")
        for n, target in zip(nodes, nodes[1:]):
            require(step(n) == (target, 1), "short-interval path failed")
        excess_bound = 512 * (m + 1) * Fraction(start - 2**m, 2**m) ** 4
        require(excess_bound < Fraction(1, 3), "short-interval energy bound failed")
        short_interval_witnesses.append({"L": length, "m": m, "q": q,
                                         "q_interval": [first_q, last_q],
                                         "q_checked_before_witness_inclusive": checked,
                                         "nodes": nodes,
                                         "rational_energy_excess_upper_bound": fraction_text(excess_bound)})

    diagnostics = []
    for length in [1, 2, 3, 6, 10]:
        m = 2 * length + 4
        start = 2**m + 2 ** (length + 1) - 1
        current = start
        differences = []
        for _ in range(length):
            current, _ = step(current)
            differences.append(energy(float(current)) - energy(float(start)))
        diagnostics.append({"L": length, "m": m, "start": start,
                            "initial_energy_excess_over_log": energy(float(start)) - log(start),
                            "energy_differences_from_start": differences})
    x = sqrt(2)
    h = 1e-6
    return {
        "scope": "Finite exact controls; the companion note proves universal statements.",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "exact_growth_families": {"L_range": [1, 64],
                                  "m_choices": ["2L+4", "2L+5", "3L+10"],
                                  "family_count": len(families), "steps_checked": step_count,
                                  "cases": families},
        "exact_same_dyadic_band": {"m_range": [0, 64], "cases": same_band},
        "exact_sturmian_hostile": {"supplied_gaps": supplied_gaps,
                                   "binary_word": binary, "length_two_sums": sums,
                                   "lower_wythoff_k_0_through_9": wythoff,
                                   "doubled_wythoff_gaps": true_gaps},
        "exact_rational_cycle": {"nodes": list(map(fraction_text, rational_cycle)),
                                  "halving_exponents": rational_exponents,
                                  "positive_integer_counterexample": False},
        "exact_squarefree_growth": {"L": 10, "q": 1, "nodes": factorizations},
        "exact_squarefree_short_interval_witnesses": short_interval_witnesses,
        "floating_diagnostics_NOT_PROOFS": {
            "at_sqrt_2": {"actual_derivative": energy_derivative(x),
                          "closed_expression": (1.5 + pi**2 / 8) / sqrt(2),
                          "blueprint_claimed_derivative": 1 / sqrt(2),
                          "central_difference": (energy(x + h) - energy(x - h)) / (2 * h)},
            "nineteen_to_twenty_nine": {"energy_19": energy(19), "energy_29": energy(29)},
            "minimum_profile_families": diagnostics,
        },
        "checks": "PASS (explicit RuntimeError checks remain active under python -O)",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"output": str(args.output), "checks": result["checks"],
                      "families": result["exact_growth_families"]["family_count"],
                      "steps": result["exact_growth_families"]["steps_checked"]}, sort_keys=True))
