#!/usr/bin/env python3
"""Primary finite-level audit for THM-4077.

The implementation builds the normalized denominator-19 map by its affine
recurrence, proves its one-bit lift law at every tested residue, and composes
the resulting permutations with the carry-cylinder coordinate of THM-2228.
"""

from __future__ import annotations

import hashlib
import json


S_MAX = 6
H_MAX = 12
CARRY_H_MAX = 10
TARGET_NAMES = ("periodic_100", "state_1", "all_1")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation_two(n: int) -> int:
    require(n != 0, "valuation of zero requested")
    n = abs(n)
    value = 0
    while n % 2 == 0:
        n //= 2
        value += 1
    return value


def parameters(s: int) -> tuple[int, int, int, int]:
    length = 18 * 2**s
    k = s + 3
    g = 3**length
    numerator = 9 * 3**k * (g - 1)
    denominator = 19 * 2**k
    require(numerator % denominator == 0, f"U integrality failed at s={s}")
    u = numerator // denominator
    require(valuation_two(g - 1) == k, f"LTE valuation failed at s={s}")
    require(u % 2 == 1 and u > 1, f"U oddness/positivity failed at s={s}")
    return length, k, g, u


def carry_prefix(state: int, height: int) -> tuple[int, ...]:
    digits: list[int] = []
    current = state
    for _ in range(height):
        digit = current % 2
        digits.append(digit)
        current = (3 * current + digit) // 2
    return tuple(digits)


def phi_residue(word: tuple[int, ...]) -> int:
    height = len(word)
    modulus = 2**height
    carry_sum = sum(
        digit * 2**index * 3 ** (height - 1 - index)
        for index, digit in enumerate(word)
    )
    return (-pow(3, -height, modulus) * carry_sum) % modulus


def target_residue(name: str, height: int) -> int:
    modulus = 2**height
    if name == "periodic_100":
        return (-9 * pow(19, -1, modulus)) % modulus
    if name == "state_1":
        return 1 % modulus
    if name == "all_1":
        return phi_residue((1,) * height)
    raise RuntimeError(f"unknown target {name}")


def main() -> None:
    layers = 0
    permutation_gates = 0
    one_bit_gates = 0
    carry_cylinder_gates = 0
    rows: list[dict[str, object]] = []

    for s in range(S_MAX + 1):
        length, k, g, u = parameters(s)
        previous_values: list[int] | None = None
        target_representatives = {name: [] for name in TARGET_NAMES}

        for height in range(1, H_MAX + 1):
            modulus = 2**height
            values: list[int] = []
            current = 0
            for t in range(modulus):
                values.append(current)
                require(current % 2 == t % 2, f"parity isometry failed at s={s},h={height},t={t}")
                current = (g * current + u) % modulus
                permutation_gates += 1

            require(len(set(values)) == modulus, f"finite permutation failed at s={s},h={height}")
            require(current == 0, f"affine cycle failed at s={s},h={height}")
            if previous_values is not None:
                half = modulus // 2
                for t in range(half):
                    require(values[t] % half == previous_values[t], f"truncation failed at s={s},h={height},t={t}")
                    require(
                        (values[t + half] - values[t]) % modulus == half,
                        f"one-bit isometry failed at s={s},h={height},t={t}",
                    )
                    one_bit_gates += 1
            previous_values = values
            inverse = {value: t for t, value in enumerate(values)}

            for name in TARGET_NAMES:
                representative = inverse[target_residue(name, height)]
                require(representative % 2 == 1, f"odd target preimage failed at s={s},h={height}")
                history = target_representatives[name]
                if history:
                    require(representative % (modulus // 2) == history[-1], f"target coherence failed for {name},s={s},h={height}")
                history.append(representative)

            if height <= CARRY_H_MAX:
                words = set()
                for t in range(1, modulus, 2):
                    word = carry_prefix(values[t], height)
                    require(word[0] == 1, f"odd carry cylinder failed at s={s},h={height},t={t}")
                    require(phi_residue(word) == values[t], f"Phi inverse failed at s={s},h={height},t={t}")
                    words.add(word)
                    carry_cylinder_gates += 1
                require(len(words) == 2 ** (height - 1), f"odd full cylinder failed at s={s},h={height}")

            layers += 1

        rows.append(
            {
                "s": s,
                "L": length,
                "k": k,
                "U_mod_4096": u % 4096,
                "targets_h12": {
                    name: target_representatives[name][-1]
                    for name in TARGET_NAMES
                },
            }
        )

    require(layers == 84, "layer aggregate changed")
    require(permutation_gates == 57330, "permutation aggregate changed")
    require(one_bit_gates == 28658, "one-bit aggregate changed")
    require(carry_cylinder_gates == 7161, "carry-cylinder aggregate changed")

    digest = hashlib.sha256(
        json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4077 primary inverse-limit audit")
    print(f"s_range=0..{S_MAX} heights=1..{H_MAX} finite_permutation_layers={layers}")
    print(f"permutation_parity_gates={permutation_gates} one_bit_lift_gates={one_bit_gates}")
    print(f"odd_carry_cylinder_gates_through_h{CARRY_H_MAX}={carry_cylinder_gates}")
    print(f"h12_target_representatives={json.dumps(rows, sort_keys=True, separators=(',', ':'))}")
    print(f"semantic_sha256={digest}")
    print("PASS: compatible permutations, one-bit isometry, and full odd carry cylinders")


if __name__ == "__main__":
    main()
