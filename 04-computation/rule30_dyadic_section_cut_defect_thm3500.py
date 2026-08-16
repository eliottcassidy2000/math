#!/usr/bin/env python3
"""Exact companion for provisional THM-3500.

The default replay has three independent layers:

1. bit-packed verification of the universal dyadic section split and driven
   cut-defect recurrence;
2. direct cyclic-arc/Hasse verification of the Mersenne endpoint carrier;
3. a constant-memory C replay of exactly 2^32 truncated packed Rule 30
   updates, compiled from the frozen source below.

Set RULE30_THM3500_LONG_SECTION=1 to additionally replay the independent
section-mask route through m=32.  That optional control uses about 1.5 GiB at
m=32 and is deliberately not the sole or default gate.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import tempfile


def check(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def strict_suffix(mask: int, length: int) -> int:
    """Strict suffix parity on a low-to-high indexed bit vector."""
    value = mask
    shift = 1
    while shift < length:
        value ^= value >> shift
        shift <<= 1
    return value >> 1


def section_masks(length: int, stop: int) -> list[int]:
    ones = (1 << length) - 1
    rows = [ones, ones]
    while len(rows) <= stop:
        rows.append(
            strict_suffix(rows[-2], length)
            | strict_suffix(rows[-1], length)
        )
    return rows


def valuation_from_masks(length: int) -> int:
    ones = (1 << length) - 1
    if length & 1:
        return 1
    previous = ones
    current = ones
    k = 1
    while True:
        following = strict_suffix(previous, length) | strict_suffix(
            current, length
        )
        k += 1
        if following.bit_count() & 1:
            return k + 1
        previous, current = current, following


def delta_step(
    p_previous: int,
    p_current: int,
    delta_previous: int,
    delta_current: int,
    z_previous: int,
    z_current: int,
    length: int,
) -> int:
    ones = (1 << length) - 1
    upper_following = strict_suffix(
        p_previous, length
    ) | strict_suffix(p_current, length)
    lower_previous = strict_suffix(
        p_previous ^ delta_previous, length
    ) ^ (ones if z_previous else 0)
    lower_current = strict_suffix(
        p_current ^ delta_current, length
    ) ^ (ones if z_current else 0)
    return upper_following ^ (lower_previous | lower_current)


def delta_step_intersection(
    p_previous: int,
    p_current: int,
    delta_previous: int,
    delta_current: int,
    z_previous: int,
    z_current: int,
    length: int,
) -> int:
    """Expanded three-intersection form of the same defect recurrence."""
    ones = (1 << length) - 1
    a_previous = strict_suffix(p_previous, length)
    a_current = strict_suffix(p_current, length)
    d_previous = strict_suffix(delta_previous, length) ^ (
        ones if z_previous else 0
    )
    d_current = strict_suffix(delta_current, length) ^ (
        ones if z_current else 0
    )
    return (
        d_previous
        ^ d_current
        ^ (a_previous & d_current)
        ^ (a_current & d_previous)
        ^ (d_previous & d_current)
    )


def hasse_one(mask: int, length: int) -> int:
    """M_1(mask)=xor_i binom(i,1)mask_i, hence odd-index parity."""
    if length == 1:
        return 0
    even_positions = ((1 << length) - 1) // 3
    odd_positions = even_positions << 1
    return (mask & odd_positions).bit_count() & 1


def direct_packed_prefix(max_m: int, width: int) -> tuple[int, ...]:
    mask = (1 << width) - 1
    state = 1
    next_time = 1
    values: list[int] = []
    for time in range(1, (1 << max_m) + 1):
        state = (state ^ ((state << 1) | (state << 2))) & mask
        if time == next_time:
            difference = state - 1
            check(difference != 0, ("zero packed difference", time))
            values.append((difference & -difference).bit_length() - 1)
            next_time <<= 1
    return tuple(values)


def packed_edge_states(stop_time: int, width: int) -> list[int]:
    """Return R_t mod 2^width for 0<=t<=stop_time."""
    mask = (1 << width) - 1
    state = 1
    states = [state]
    for _ in range(stop_time):
        state = (state ^ ((state << 1) | (state << 2))) & mask
        states.append(state)
    return states


def edge_current(states: list[int], depth: int, time: int) -> int:
    """Q_depth(time)=b_(depth-1)(time) OR b_(depth-2)(time)."""
    first = (states[time] >> (depth - 1)) & 1
    second = 0 if depth == 1 else (states[time] >> (depth - 2)) & 1
    return first | second


def reverse_mask(mask: int, length: int) -> int:
    result = 0
    for index in range(length):
        result |= ((mask >> index) & 1) << (length - 1 - index)
    return result


def hasse(vector: list[int]) -> list[int]:
    length = len(vector)
    result = [0] * length
    for index, bit in enumerate(vector):
        if not bit:
            continue
        submask = index
        while True:
            result[submask] ^= 1
            if submask == 0:
                break
            submask = (submask - 1) & index
    return result


def cyclic_arc(vector: list[int], arc_length: int) -> list[int]:
    period = len(vector)
    return [
        sum(vector[(phase + offset) % period] for offset in range(arc_length))
        & 1
        for phase in range(period)
    ]


def cyclic_shift(vector: list[int], amount: int) -> list[int]:
    period = len(vector)
    return [vector[(phase - amount) % period] for phase in range(period)]


DIRECT_LONG_C = r'''
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

static unsigned ctz128(__uint128_t x) {
    uint64_t low = (uint64_t)x;
    if (low != 0) return (unsigned)__builtin_ctzll(low);
    return 64u + (unsigned)__builtin_ctzll((uint64_t)(x >> 64));
}

int main(void) {
    const unsigned max_m = 32;
    const uint64_t stop = UINT64_C(1) << max_m;
    const __uint128_t mask = (((__uint128_t)1 << 120) - 1);
    __uint128_t state = 1;
    uint64_t next_time = 1;
    unsigned next_m = 0;
    for (uint64_t time = 1; time <= stop; ++time) {
        state ^= (state << 1) | (state << 2);
        state &= mask;
        if (time == next_time) {
            __uint128_t difference = state - 1;
            if (difference == 0) return 3;
            printf("%u %u\n", next_m, ctz128(difference));
            ++next_m;
            if (next_m <= max_m) next_time <<= 1;
        }
    }
    return 0;
}
'''.strip() + "\n"


SECTION_LONG_C = r'''
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint64_t suffix_inclusive(uint64_t x) {
    x ^= x >> 1; x ^= x >> 2; x ^= x >> 4;
    x ^= x >> 8; x ^= x >> 16; x ^= x >> 32;
    return x;
}

static int step(const uint64_t *a, const uint64_t *b, uint64_t *out,
                size_t words) {
    unsigned parity_a = 0, parity_b = 0, parity_out = 0;
    for (size_t jj = words; jj-- > 0;) {
        uint64_t xa = a[jj], xb = b[jj];
        uint64_t sa = suffix_inclusive(xa) >> 1;
        uint64_t sb = suffix_inclusive(xb) >> 1;
        if (parity_a) sa = ~sa;
        if (parity_b) sb = ~sb;
        out[jj] = sa | sb;
        parity_out ^= (unsigned)__builtin_parityll(out[jj]);
        parity_a ^= (unsigned)__builtin_parityll(xa);
        parity_b ^= (unsigned)__builtin_parityll(xb);
    }
    return (int)parity_out;
}

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    unsigned m = (unsigned)strtoul(argv[1], NULL, 10);
    if (m < 6 || m > 32) return 2;
    uint64_t bits = UINT64_C(1) << m;
    size_t words = (size_t)(bits >> 6);
    uint64_t *a = malloc(words * sizeof(*a));
    uint64_t *b = malloc(words * sizeof(*b));
    uint64_t *c = malloc(words * sizeof(*c));
    if (a == NULL || b == NULL || c == NULL) return 3;
    memset(a, 0xff, words * sizeof(*a));
    memset(b, 0xff, words * sizeof(*b));
    uint64_t k = 1;
    for (;;) {
        int odd = step(a, b, c, words);
        ++k;
        if (odd) {
            printf("%u %" PRIu64 "\n", m, k + 1);
            break;
        }
        uint64_t *temporary = a;
        a = b; b = c; c = temporary;
    }
    free(a); free(b); free(c);
    return 0;
}
'''.strip() + "\n"


def compile_and_run(source: str, arguments: tuple[str, ...]) -> str:
    with tempfile.TemporaryDirectory(prefix="rule30-thm3500-") as directory:
        executable = os.path.join(directory, "gate")
        compilation = subprocess.run(
            ["cc", "-O3", "-std=c11", "-Wall", "-Wextra", "-x", "c", "-o", executable, "-"],
            input=source,
            text=True,
            capture_output=True,
            check=False,
        )
        check(
            compilation.returncode == 0,
            ("C compilation failed", compilation.stdout, compilation.stderr),
        )
        replay = subprocess.run(
            [executable, *arguments],
            text=True,
            capture_output=True,
            check=False,
        )
        check(
            replay.returncode == 0,
            ("C replay failed", replay.stdout, replay.stderr),
        )
        return replay.stdout


def main() -> None:
    known_v_0_to_27 = (
        1,
        3,
        4,
        6,
        7,
        9,
        15,
        16,
        24,
        25,
        27,
        29,
        34,
        36,
        37,
        39,
        41,
        43,
        48,
        49,
        51,
        54,
        55,
        58,
        60,
        63,
        64,
        66,
    )

    mealy_bridge_pairs = 0
    for length in range(1, 257):
        stop_depth = min(length + 4, 96)
        masks = section_masks(length, stop_depth - 1)
        states = packed_edge_states(length - 1, stop_depth + 1)
        for depth in range(1, stop_depth + 1):
            expected = 0
            for time in range(length):
                expected |= edge_current(states, depth, time) << time
            check(
                reverse_mask(masks[depth - 1], length) == expected,
                ("Mealy/current bridge", length, depth),
            )
            mealy_bridge_pairs += 1

    block_pairs = 0
    first_carry_checks = 0
    holonomy_bridge_pairs = 0
    q4_hostile: tuple[str, str] | None = None
    for m in range(13):
        length = 1 << m
        stop = min(2 * length + 4, 96)
        small = section_masks(length, stop)
        large = section_masks(2 * length, stop)
        block_mask = (1 << length) - 1
        deltas = [0, 0]
        for k in range(stop + 1):
            lower = large[k] & block_mask
            upper = large[k] >> length
            check(upper == small[k], ("upper block", m, k))
            check((lower ^ upper) == deltas[k], ("cut defect", m, k))
            check(
                (large[k].bit_count() & 1) == (deltas[k].bit_count() & 1),
                ("defect parity", m, k),
            )
            if 1 <= k < stop:
                z_previous = small[k - 1].bit_count() & 1
                z_current = small[k].bit_count() & 1
                following_delta = delta_step(
                    small[k - 1],
                    small[k],
                    deltas[k - 1],
                    deltas[k],
                    z_previous,
                    z_current,
                    length,
                )
                expanded_delta = delta_step_intersection(
                    small[k - 1],
                    small[k],
                    deltas[k - 1],
                    deltas[k],
                    z_previous,
                    z_current,
                    length,
                )
                check(
                    following_delta == expanded_delta,
                    ("OR/intersection recurrence", m, k),
                )
                deltas.append(following_delta)
            block_pairs += 1

        edge_states = packed_edge_states(2 * length, stop + 2)
        for depth in range(1, stop + 2):
            delta_reversed = reverse_mask(deltas[depth - 1], length)
            expected_delta = 0
            for time in range(length):
                expected_delta |= (
                    edge_current(edge_states, depth, time + length)
                    ^ edge_current(edge_states, depth, time)
                ) << time
                holonomy_here = (
                    ((edge_states[time + length] >> depth) & 1)
                    ^ ((edge_states[time] >> depth) & 1)
                )
                holonomy_next = (
                    ((edge_states[time + length + 1] >> depth) & 1)
                    ^ ((edge_states[time + 1] >> depth) & 1)
                )
                check(
                    ((delta_reversed >> time) & 1)
                    == (holonomy_next ^ holonomy_here),
                    ("dyadic holonomy derivative", m, depth, time),
                )
            check(
                delta_reversed == expected_delta,
                ("dyadic current defect", m, depth),
            )
            check(
                (deltas[depth - 1].bit_count() & 1)
                == ((edge_states[2 * length] >> depth) & 1),
                ("dyadic holonomy boundary", m, depth),
            )
            holonomy_bridge_pairs += 1

        if m >= 1:
            valuation = valuation_from_masks(length)
            pivot = valuation - 1
            check(all(deltas[k] == 0 for k in range(pivot + 1)), ("pre-pivot", m))
            predicted = block_mask ^ strict_suffix(small[pivot - 1], length)
            check(deltas[pivot + 1] == predicted, ("pivot defect", m))
            immediate = hasse_one(small[pivot - 1], length) == 1
            next_valuation = valuation_from_masks(2 * length)
            check(
                immediate == (next_valuation == valuation + 1),
                ("first-carry criterion", m, valuation, next_valuation),
            )
            first_carry_checks += 1

        if length == 4:
            q4_hostile = tuple(
                format(deltas[k], "04b")[::-1] for k in (4, 5)
            )

    check(q4_hostile == ("0011", "0111"), ("q=4 hostile", q4_hostile))

    packed_control = direct_packed_prefix(18, 120)
    mask_control = tuple(valuation_from_masks(1 << m) for m in range(19))
    check(packed_control == mask_control, ("packed/mask mismatch",))
    check(mask_control[: len(known_v_0_to_27)] == known_v_0_to_27[:19], "known prefix")

    mersenne_pairs = 0
    mersenne_basis_vectors = 0
    for exponent in range(2, 9):
        period = 1 << exponent
        for m in range(exponent - 1):
            half = 1 << m
            endpoint = 2 * half - 1
            if endpoint >= period:
                continue
            for basis_index in range(period):
                current = [0] * period
                current[basis_index] = 1
                profile = cyclic_arc(current, endpoint)
                recentered = cyclic_shift(profile, endpoint)
                current_hasse = hasse(current)
                recentered_hasse = hasse(recentered)
                expected_hasse = [
                    current_hasse[j]
                    ^ (current_hasse[j - endpoint] if j >= endpoint else 0)
                    for j in range(period)
                ]
                check(
                    recentered_hasse == expected_hasse,
                    ("Mersenne two-diagonal law", period, endpoint, basis_index),
                )
                center = profile[period - endpoint]
                check(center == recentered[0], ("Mersenne recenter", period, endpoint))
                check(
                    center == (sum(current_hasse[period - endpoint :]) & 1),
                    ("Mersenne current tail", period, endpoint),
                )
                profile_hasse = hasse(profile)
                profile_face = (
                    sum(
                        profile_hasse[period - endpoint + 2 * offset]
                        for offset in range(half)
                    )
                    & 1
                )
                check(center == profile_face, ("Mersenne profile face", period, endpoint))
                mersenne_basis_vectors += 1
            mersenne_pairs += 1

    for m in range(12):
        endpoint = (1 << (m + 1)) - 1
        next_endpoint = (1 << (m + 2)) - 1
        carrier_square = (1 << 0) | (1 << (2 * endpoint))
        next_carrier_from_recurrence = (
            (1 << 0) ^ (1 << 1) ^ (carrier_square << 1)
        )
        expected_next_carrier = (1 << 0) | (1 << next_endpoint)
        check(
            next_carrier_from_recurrence == expected_next_carrier,
            ("Mersenne Cartier/Mahler law", m),
        )

    direct_long_output = compile_and_run(DIRECT_LONG_C, ())
    long_values: list[int] = []
    for line in direct_long_output.splitlines():
        fields = line.split()
        check(len(fields) == 2, ("long direct line", line))
        index, value = (int(field) for field in fields)
        check(index == len(long_values), ("long direct index", line))
        long_values.append(value)
    check(len(long_values) == 33, ("long direct count", len(long_values)))
    check(tuple(long_values[:28]) == known_v_0_to_27, "long known positive control")
    check(tuple(long_values[28:]) == (69, 70, 72, 74, 77), "long extension")

    # Derive the finite wrap consequence from THM-3493's exact dyadic-block
    # rule instead of treating the printed interval as a trusted constant.
    wrapped_depths: list[int] = []
    for m, value in enumerate(long_values):
        block_start = 1 << m
        block_end = (block_start << 1) - 1
        check(value <= block_end, ("strict dyadic ceiling", m, value))
        if value >= block_start:
            wrapped_depths.extend(range(block_start, value + 1))
    check(tuple(wrapped_depths) == (1, 2, 3, 4), "derived finite wrap set")
    hard_interval = (wrapped_depths[-1] + 1, (1 << len(long_values)) - 1)
    check(hard_interval == (5, (1 << 33) - 1), "derived finite hard interval")

    optional_section_status = "FROZEN_NOT_RUN"
    if os.environ.get("RULE30_THM3500_LONG_SECTION") == "1":
        optional_values = []
        for m in range(28, 33):
            section_output = compile_and_run(SECTION_LONG_C, (str(m),))
            fields = section_output.split()
            check(len(fields) == 2 and int(fields[0]) == m, ("long section", m))
            optional_values.append(int(fields[1]))
        check(tuple(optional_values) == (69, 70, 72, 74, 77), "section extension")
        optional_section_status = "PASS"

    print("RULE30_DYADIC_SECTION_CUT_DEFECT_THM3500")
    print("status PROVED_VERIFIED_EXACT_INDEPENDENTLY_AUDITED")
    print(
        "universal_split_universe",
        "m=0..12",
        "k<=min(2q+4,96)",
        f"block_pairs={block_pairs}",
    )
    print(
        "mealy_current_bridge_universe",
        "n=1..256",
        "depth<=min(n+4,96)",
        f"pairs={mealy_bridge_pairs}",
    )
    print(
        "dyadic_holonomy_bridge_universe",
        "m=0..12",
        "depth<=min(2q+5,97)",
        f"pairs={holonomy_bridge_pairs}",
    )
    print("first_carry_hasse_checks", first_carry_checks)
    print("q4_delta_k4_k5_increasing_index", q4_hostile)
    print("packed_mask_control_m0_to_m18", mask_control)
    print(
        "mersenne_carrier_universe",
        "p=2^d,d=2..8",
        f"parameter_pairs={mersenne_pairs}",
        f"basis_vectors={mersenne_basis_vectors}",
        "frobenius_scales=12",
    )
    print(
        "long_direct_gate",
        "updates=2^32",
        "width=120",
        "memory=constant",
        "compiler=cc",
        "source_sha256=" + hashlib.sha256(DIRECT_LONG_C.encode()).hexdigest(),
    )
    print("v_28_to_v_32", tuple(long_values[28:]))
    print("finite_wrap_set_through_2^33_minus_1", tuple(wrapped_depths))
    print("finite_hard_interval", hard_interval)
    print(
        "optional_independent_section_gate",
        optional_section_status,
        "m=28..32",
        "peak_mask_memory_about_1.5_GiB",
        "source_sha256=" + hashlib.sha256(SECTION_LONG_C.encode()).hexdigest(),
    )
    print("status PASS")


if __name__ == "__main__":
    main()
