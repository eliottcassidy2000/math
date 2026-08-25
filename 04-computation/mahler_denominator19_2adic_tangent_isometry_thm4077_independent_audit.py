#!/usr/bin/env python3
"""Independent Hensel and native-orbit audit for THM-4077.

This path evaluates geometric sums by binary doubling, lifts three coherent
targets one bit at a time through height 64, and directly replays selected
ordinary denominator-19 launches.  It imports no primary implementation.
"""

from __future__ import annotations

import hashlib
import json


S_MAX = 6
H_MAX = 64
REPLAY_S_MAX = 2
REPLAY_HEIGHTS = (4, 6, 8, 10)
TARGET_NAMES = ("periodic_100", "state_1", "all_1")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def make_parameters(s: int) -> tuple[int, int, int, int]:
    length = 18 << s
    k = s + 3
    base = pow(3, length)
    top = 9 * pow(3, k) * (base - 1)
    bottom = 19 << k
    require(top % bottom == 0, f"normalizer integrality failed at s={s}")
    normalizer = top // bottom
    require(normalizer & 1 == 1, f"normalizer parity failed at s={s}")
    return length, k, base, normalizer


def geometric_mod(base: int, exponent: int, modulus: int) -> int:
    power_result = 1
    sum_result = 0
    block_power = base % modulus
    block_sum = 1 % modulus
    remaining = exponent
    while remaining:
        if remaining & 1:
            sum_result = (sum_result + power_result * block_sum) % modulus
            power_result = power_result * block_power % modulus
        block_sum = block_sum * (1 + block_power) % modulus
        block_power = block_power * block_power % modulus
        remaining >>= 1
    return sum_result


def evaluate(normalizer: int, base: int, parameter: int, height: int) -> int:
    modulus = 1 << height
    return normalizer * geometric_mod(base, parameter, modulus) % modulus


def phi_prefix_residue(word: tuple[int, ...]) -> int:
    height = len(word)
    modulus = 1 << height
    accumulator = 0
    power_two = 1
    power_three = pow(3, height - 1)
    for index, digit in enumerate(word):
        if digit:
            accumulator += power_two * power_three
        power_two <<= 1
        if index + 1 < height:
            power_three //= 3
    return -pow(pow(3, height), -1, modulus) * accumulator % modulus


def target_at_height(name: str, height: int) -> int:
    modulus = 1 << height
    if name == "periodic_100":
        return -9 * pow(19, -1, modulus) % modulus
    if name == "state_1":
        return 1 % modulus
    if name == "all_1":
        return phi_prefix_residue((1,) * height)
    raise RuntimeError(f"unknown target {name}")


def lift_target(normalizer: int, base: int, name: str) -> tuple[list[int], int]:
    representatives = [1]
    candidate_gates = 0
    require(evaluate(normalizer, base, 1, 1) == target_at_height(name, 1), f"height-one target failed for {name}")
    for height in range(2, H_MAX + 1):
        old = representatives[-1]
        candidates = (old, old + (1 << (height - 1)))
        target = target_at_height(name, height)
        valid: list[int] = []
        for candidate in candidates:
            if evaluate(normalizer, base, candidate, height) == target:
                valid.append(candidate)
            candidate_gates += 1
        require(len(valid) == 1, f"unique Hensel lift failed for {name},h={height}")
        representatives.append(valid[0])
    return representatives, candidate_gates


def carry_prefix(state: int, height: int) -> tuple[int, ...]:
    word: list[int] = []
    current = state
    for _ in range(height):
        digit = current & 1
        word.append(digit)
        current = (3 * current + digit) // 2
    return tuple(word)


def expected_word(name: str, height: int) -> tuple[int, ...]:
    if name == "periodic_100":
        return tuple((1, 0, 0)[index % 3] for index in range(height))
    if name == "state_1":
        return carry_prefix(1, height)
    if name == "all_1":
        return (1,) * height
    raise RuntimeError(f"unknown target {name}")


def direct_launch_replay(s: int, height: int, parameter: int, name: str) -> int:
    length, k, base, normalizer = make_parameters(s)
    m = length * parameter
    require(pow(2, m, 19) == 1, f"launch divisibility failed at s={s},h={height}")
    launch = 9 * ((1 << m) - 1) // 19
    current = launch
    steps = 0
    for _ in range(m + k):
        digit = current & 1
        current = (3 * current + digit) // 2
        steps += 1
    expected_state = normalizer * ((pow(base, parameter) - 1) // (base - 1))
    require(current == expected_state, f"moving observation identity failed at s={s},h={height},{name}")
    require(carry_prefix(current, height) == expected_word(name, height), f"launch prefix failed at s={s},h={height},{name}")
    return steps + height


def low_bits(value: int, height: int) -> str:
    return "".join(str((value >> index) & 1) for index in range(height))


def main() -> None:
    all_representatives: dict[int, dict[str, list[int]]] = {}
    hensel_candidate_gates = 0
    logarithm_identity_gates = 0

    for s in range(S_MAX + 1):
        _, k, base, normalizer = make_parameters(s)
        by_target: dict[str, list[int]] = {}
        for name in TARGET_NAMES:
            representatives, gates = lift_target(normalizer, base, name)
            by_target[name] = representatives
            hensel_candidate_gates += gates
            final = representatives[-1]
            require(evaluate(normalizer, base, final, H_MAX) == target_at_height(name, H_MAX), f"final target failed for {name},s={s}")

        for height in range(1, H_MAX + 1):
            modulus = 1 << height
            parameter = by_target["periodic_100"][height - 1]
            left = pow(base, parameter, modulus)
            right = (1 - pow(2 * pow(3, -1, modulus), k, modulus)) % modulus
            require(left == right, f"safe logarithm identity failed at s={s},h={height}")
            logarithm_identity_gates += 1
        all_representatives[s] = by_target

    orbit_steps = 0
    replay_rows: list[dict[str, int | str]] = []
    for s in range(REPLAY_S_MAX + 1):
        for height in REPLAY_HEIGHTS:
            for name in TARGET_NAMES:
                parameter = all_representatives[s][name][height - 1]
                orbit_steps += direct_launch_replay(s, height, parameter, name)
                replay_rows.append({"s": s, "h": height, "target": name, "t": parameter})

    require(hensel_candidate_gates == 2646, "Hensel aggregate changed")
    require(logarithm_identity_gates == 448, "logarithm aggregate changed")

    bit_rows = {
        name: low_bits(all_representatives[0][name][-1], H_MAX)
        for name in TARGET_NAMES
    }
    require(
        bit_rows["periodic_100"] == "1010100100000000011010001100100001100110011100011100101000101100",
        "periodic target bit control changed",
    )
    require(
        bit_rows["state_1"] == "1001100000110001011111101010100011100110110110001010000100000011",
        "state-one target bit control changed",
    )
    require(
        bit_rows["all_1"] == "1111111011011011010110010111001011010010101110011010000000001111",
        "all-one target bit control changed",
    )

    semantic = {
        "bits_s0": bit_rows,
        "replays": replay_rows,
        "orbit_steps": orbit_steps,
    }
    digest = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    print("THM-4077 independent Hensel/native audit")
    print(f"s_range=0..{S_MAX} target_heights=1..{H_MAX} hensel_candidate_gates={hensel_candidate_gates}")
    print(f"safe_logarithm_identity_gates={logarithm_identity_gates}")
    print(f"direct_launch_replays={len(replay_rows)} ordinary_orbit_steps={orbit_steps}")
    print(f"s0_parameter_bits_lsb_first={json.dumps(bit_rows, sort_keys=True, separators=(',', ':'))}")
    print(f"semantic_sha256={digest}")
    print("PASS: unique 64-bit lifts, tangent identity, and direct moving-time launches")


if __name__ == "__main__":
    main()
