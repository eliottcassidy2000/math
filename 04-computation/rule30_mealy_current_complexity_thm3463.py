#!/usr/bin/env python3
"""Exact audit for THM-3463's Rule-30 diagonal-frontier package.

The finite checks cover the edge Mealy codec, its suffix-parity compiler,
the Rule-150 collision-current decomposition, the period-lift restriction,
the arrival-trace atlas, the structured Boolean-cube derivative, and a large
finite 2-kernel hostile.  Universal statements are proved in the theorem;
no finite prefix is extrapolated.
"""

from __future__ import annotations

import hashlib
from itertools import product


def check(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def rule30(left: int, center: int, right: int) -> int:
    return left ^ center ^ right ^ (center & right)


def phi(x: int, width: int | None = None) -> int:
    y = x ^ ((x << 1) | (x << 2))
    if width is not None:
        y &= (1 << width) - 1
    return y


def edge_rows(horizon: int) -> list[int]:
    rows = []
    x = 1
    for _ in range(horizon + 1):
        rows.append(x)
        x = phi(x)
    return rows


def direct_local_rows(horizon: int) -> list[int]:
    row = [1]
    packed = []
    for _ in range(horizon + 1):
        packed.append(sum(bit << k for k, bit in enumerate(reversed(row))))
        padded = [0, 0] + row + [0, 0]
        row = [rule30(padded[i], padded[i + 1], padded[i + 2]) for i in range(len(padded) - 2)]
    return packed


# A state transition is input -> (output, section).
MEALY = {
    "A": {0: (0, "A"), 1: (1, "B")},
    "B": {0: (1, "C"), 1: (0, "B")},
    "C": {0: (1, "A"), 1: (0, "B")},
}
ACTIVITY = {"A": 0, "B": 1, "C": 1}


def transduce(state: str, bits: tuple[int, ...]) -> tuple[tuple[int, ...], str]:
    output = []
    for bit in bits:
        out, state = MEALY[state][bit]
        output.append(out)
    return tuple(output), state


def bits_to_int(bits: tuple[int, ...]) -> int:
    return sum(bit << k for k, bit in enumerate(bits))


def strict_suffix_parity(bits: list[int]) -> list[int]:
    parity = 0
    result = [0] * len(bits)
    for i in range(len(bits) - 1, -1, -1):
        result[i] = parity
        parity ^= bits[i]
    return result


def strict_suffix_parity_int(bits: int, width: int) -> int:
    """Bit-parallel strict suffix XOR on positions 0,...,width-1."""
    result = bits >> 1
    shift = 1
    while shift < width:
        result ^= result >> shift
        shift <<= 1
    return result & ((1 << width) - 1)


def section_word_at_zero(word: list[str]) -> list[str]:
    suffix_activity = strict_suffix_parity([ACTIVITY[state] for state in word])
    return [MEALY[state][incoming][1] for state, incoming in zip(word, suffix_activity)]


def mealy_center(n: int) -> int:
    if n == 0:
        return 1
    word = ["B"] * n
    for _ in range(n - 1):
        word = section_word_at_zero(word)
    return sum(ACTIVITY[state] for state in word) & 1


def suffix_center(n: int) -> int:
    if n == 0:
        return 1
    previous = (1 << n) - 1
    current = previous
    if n == 1:
        return 1
    for _ in range(1, n - 1):
        s_previous = strict_suffix_parity_int(previous, n)
        s_current = strict_suffix_parity_int(current, n)
        following = s_previous | s_current
        previous, current = current, following
    return current.bit_count() & 1


def audit_mealy(rows: list[int]) -> str:
    check(direct_local_rows(128) == rows[:129], "local Rule30 versus packed edge")
    for width in range(1, 12):
        for x in range(1 << width):
            bits = tuple((x >> k) & 1 for k in range(width))
            output, _ = transduce("A", bits)
            check(bits_to_int(output) == phi(x, width), f"Mealy A width={width} x={x}")

    for n in range(1, 65):
        word = ["B"] * n
        previous = [1] * n
        current = [1] * n
        check([ACTIVITY[state] for state in word] == previous, f"section p0 n={n}")
        for depth in range(n - 1):
            word = section_word_at_zero(word)
            if depth == 0:
                expected = current
            else:
                s_previous = strict_suffix_parity(previous)
                s_current = strict_suffix_parity(current)
                expected = [a | b for a, b in zip(s_previous, s_current)]
                previous, current = current, expected
            check([ACTIVITY[state] for state in word] == expected, f"section recurrence n={n} k={depth + 1}")

    centers = []
    for n in range(1, 513):
        direct = (rows[n] >> n) & 1
        compiled = suffix_center(n)
        check(compiled == direct, f"suffix center n={n}")
        if n <= 80:
            check(mealy_center(n) == direct, f"Mealy center n={n}")
        centers.append(str(direct))
    return "1" + "".join(centers[:63])


def linear_step(x: int) -> int:
    return x ^ (x << 1) ^ (x << 2)


def collision_injection(x: int) -> int:
    return (x & (x >> 1)) << 2


def linear_power(x: int, exponent: int) -> int:
    for _ in range(exponent):
        x = linear_step(x)
    return x


def green(n: int, r: int) -> int:
    if r < 0 or r > 2 * n:
        return 0
    return (linear_power(1, n) >> r) & 1


def audit_collision_current(rows: list[int]) -> tuple[int, int]:
    for t in range(128):
        check(phi(rows[t]) == linear_step(rows[t]) ^ collision_injection(rows[t]), f"one-step current t={t}")

    for t in range(81):
        reconstructed = linear_power(1, t)
        for s in range(t):
            reconstructed ^= linear_power(collision_injection(rows[s]), t - 1 - s)
        check(reconstructed == rows[t], f"Duhamel t={t}")
        check(((linear_power(1, t) >> t) & 1) == 1, f"central Rule150 t={t}")

    for m in range(33):
        for r in range(2 * m + 2):
            check(green(2 * m, 2 * r) == green(m, r), f"Green even-even m={m} r={r}")
            check(green(2 * m, 2 * r + 1) == 0, f"Green even-odd m={m} r={r}")
            check(green(2 * m + 1, 2 * r) == (green(m, r) ^ green(m, r - 1)), f"Green odd-even m={m} r={r}")
            check(green(2 * m + 1, 2 * r + 1) == green(m, r), f"Green odd-odd m={m} r={r}")

    dyadic_checks = 0
    for h in (1, 2, 4, 8, 16, 32, 64):
        for r in range(h):
            defect = 0
            for j in range(h):
                defect ^= linear_power(collision_injection(rows[r + j]), h - 1 - j)
            lhs = ((rows[h + r] >> (h + r)) & 1) ^ ((rows[r] >> r) & 1)
            rhs = (defect >> (h + r)) & 1
            check(lhs == rhs, f"dyadic current h={h} r={r}")
            dyadic_checks += 1
    return 81, dyadic_checks


def seed_period(width: int) -> int:
    x = 1
    period = 0
    while True:
        x = phi(x, width)
        period += 1
        if x == 1:
            return period


def lift_bit(width: int, period: int) -> int:
    x = 1
    for _ in range(period):
        x = phi(x, width + 1)
    return (x >> width) & 1


def pair_counts(width: int, period: int) -> tuple[int, int, int, int]:
    x = 1
    counts = [0, 0, 0, 0]
    for _ in range(period):
        high = (x >> (width - 1)) & 1
        low = (x >> (width - 2)) & 1
        counts[2 * high + low] += 1
        x = phi(x, width)
    return tuple(counts)  # type: ignore[return-value]


def audit_period_lifts() -> tuple[list[int], str, list[int]]:
    periods = [seed_period(width) for width in range(1, 37)]
    eps = [lift_bit(width, periods[width - 1]) for width in range(1, 37)]
    check(eps[:2] == [1, 0], "initial lift word")
    for width in range(1, 36):
        check(periods[width] == periods[width - 1] * (1 << eps[width - 1]), f"period lift width={width}")
    for width in range(2, 37):
        counts = pair_counts(width, periods[width - 1])
        check(eps[width - 1] == (counts[0] & 1), f"zero-pair parity width={width}")
    for width in range(4, 37):
        if eps[width - 3] == eps[width - 2] == 1:
            p = periods[width - 3]
            check(pair_counts(width, periods[width - 1]) == (p, p, p, p), f"cube uniformity width={width}")
            check(eps[width - 1] == 0, f"forbidden triple width={width}")
    word = "".join(map(str, eps))
    check("111" not in word, "finite no-111 audit")
    for width, period in enumerate(periods, 1):
        exponent_cap = (2 * width + 2) // 3 - 1  # ceil(2w/3)-1
        check(period <= 1 << exponent_cap, f"period cap width={width}")
    allowed3 = {"".join(bits) for bits in product("01", repeat=3) if "111" not in "".join(bits)}
    observed3 = {word[i : i + 3] for i in range(len(word) - 2)}
    check(observed3 == allowed3, "all allowed length-three lift blocks")
    return periods, word, [period.bit_length() - 1 for period in periods]


def audit_arrival_atlas(periods: list[int], eps_word: str) -> dict[int, tuple[int, int]]:
    eps = [int(bit) for bit in eps_word]
    horizon = periods[30] + 31
    rows = edge_rows(horizon)
    atlas = {}
    for cap in range(1, 31):
        innovations = [k for k in range(1, cap + 1) if eps[k - 1]]
        period = periods[cap]
        words = {
            tuple((rows[h + k] >> k) & 1 for k in innovations)
            for h in range(period)
        }
        check(len(words) == period == 1 << len(innovations), f"arrival atlas cap={cap}")
        check(words == set(product((0, 1), repeat=len(innovations))), f"arrival cube cap={cap}")
        atlas[cap] = (len(innovations), period)
    return atlas


def center_from_cone(state: int, width: int) -> int:
    while width > 1:
        state = state ^ (state >> 1) ^ (state >> 2) ^ ((state >> 1) & (state >> 2))
        width -= 2
        state &= (1 << width) - 1
    return state & 1


def mobius_anf(table: list[int], variables: int) -> list[int]:
    coefficients = table[:]
    for bit in range(variables):
        stride = 1 << bit
        for mask in range(1 << variables):
            if mask & stride:
                coefficients[mask] ^= coefficients[mask ^ stride]
    return coefficients


def audit_boolean_cube() -> tuple[list[int], list[tuple[int, ...]]]:
    degrees = []
    marker_rows = []
    for n in range(1, 9):
        variables = 2 * n + 1
        table = [center_from_cone(state, variables) for state in range(1 << variables)]
        coefficients = mobius_anf(table, variables)
        marker_positions = [n] + [n - 1 - 2 * s for s in range(n)]
        marker_mask = sum(1 << (position + n) for position in marker_positions)
        check(coefficients[marker_mask] == 1, f"marker coefficient n={n}")
        supersets = [mask for mask, value in enumerate(coefficients) if value and mask & marker_mask == marker_mask]
        check(supersets == [marker_mask], f"marker maximality n={n}")
        degree = max(mask.bit_count() for mask, value in enumerate(coefficients) if value)
        degrees.append(degree)

        cube_xor = 0
        ordered_positions = sorted(marker_positions)
        for toggles in range(1 << (n + 1)):
            state = 1 << n
            for i, position in enumerate(ordered_positions):
                if (toggles >> i) & 1:
                    state ^= 1 << (position + n)
            cube_xor ^= center_from_cone(state, variables)
        check(cube_xor == 1, f"seed-centered cube parity n={n}")
        marker_rows.append(tuple(ordered_positions))
    check(degrees == [2, 3, 5, 7, 9, 11, 13, 15], "finite full-degree prefix")
    return degrees, marker_rows


def audit_scalar_nonclosure() -> tuple[str, str, int, int]:
    p = [1, 1, 1]
    q0 = [1, 1, 1]
    q1 = [0, 1, 0]
    sp = strict_suffix_parity(p)
    sq0 = strict_suffix_parity(q0)
    sq1 = strict_suffix_parity(q1)
    marginals0 = (sum(p) & 1, sum(i * bit for i, bit in enumerate(p)) & 1,
                  sum(q0) & 1, sum(i * bit for i, bit in enumerate(q0)) & 1,
                  sum(sp) & 1, sum(sq0) & 1)
    marginals1 = (sum(p) & 1, sum(i * bit for i, bit in enumerate(p)) & 1,
                  sum(q1) & 1, sum(i * bit for i, bit in enumerate(q1)) & 1,
                  sum(sp) & 1, sum(sq1) & 1)
    check(marginals0 == marginals1 == (1, 1, 1, 1, 1, 1), "scalar marginals")
    intersection0 = sum(a & b for a, b in zip(sp, sq0)) & 1
    intersection1 = sum(a & b for a, b in zip(sp, sq1)) & 1
    or0 = sum(a | b for a, b in zip(sp, sq0)) & 1
    or1 = sum(a | b for a, b in zip(sp, sq1)) & 1
    check((intersection0, intersection1, or0, or1) == (1, 0, 1, 0), "intersection hostile")
    return "".join(map(str, sp)), "".join(map(str, sq1)), intersection0, intersection1


def audit_kernel_hostile() -> tuple[int, int, int, str, str]:
    horizon = 393_216
    center = bytearray(horizon)
    edge = 1
    mirror = 1
    for t in range(horizon):
        center[t] = (edge >> t) & 1
        if t < 131_072:
            check(center[t] == ((mirror >> t) & 1), f"mirror center t={t}")
            mirror = (mirror << 2) ^ (mirror << 1) ^ mirror ^ ((mirror << 1) & mirror)
        edge = phi(edge)

    hash_262144 = hashlib.sha256(center[:262_144]).hexdigest()
    hash_393216 = hashlib.sha256(center).hexdigest()
    check(hash_262144 == "75be8407c265798fa046baa3ba0f9336e4cbe27bff0be21aeba3e87a7681bea4", "2-kernel prefix hash")
    check(hash_393216 == "7a7debfac229950122dd4272ce7e183fb724fa4298ed96fc3b5ca49023346733", "dyadic-universality hash")

    profiles16 = {
        bytes(center[offset + (1 << exponent) * n] for n in range(16))
        for exponent in range(15)
        for offset in range(1 << exponent)
    }
    profiles12 = {
        bytes(center[offset + (1 << exponent) * n] for n in range(12))
        for exponent in range(16)
        for offset in range(1 << exponent)
    }
    positive_profiles15 = {
        bytes(center[offset + (1 << exponent) * n] for n in range(1, 16))
        for exponent in range(15)
        for offset in range(1 << exponent)
    }
    check(len(profiles16) == 25_830, "length-16 kernel profiles")
    check(len(profiles12) == 4_096, "length-12 dyadic universality")
    check(len(positive_profiles15) == 20_794, "positive-index length-15 kernel profiles")
    return len(profiles16), len(positive_profiles15), len(profiles12), hash_262144, hash_393216


def main() -> None:
    rows = edge_rows(1024)
    center_prefix = audit_mealy(rows)
    duhamel_horizon, dyadic_checks = audit_collision_current(rows)
    periods, eps_word, period_exponents = audit_period_lifts()
    atlas = audit_arrival_atlas(periods, eps_word)
    degrees, markers = audit_boolean_cube()
    hostile = audit_scalar_nonclosure()
    kernel16, positive_kernel15, kernel12, hash1, hash2 = audit_kernel_hostile()

    print("THM-3463 EXACT AUDIT")
    print("mealy_states=A:(A,B),B:(C,B)sigma,C:(A,B)sigma")
    print("direct_local_rule_vs_packed_rows_t=0..128")
    print("mealy_vs_phi_all_words_width=1..11")
    print("mealy_sections_n=1..80 suffix_compiler_n=1..512")
    print(f"center_prefix_t0..63={center_prefix}")
    print(f"collision_duhamel_t=0..{duhamel_horizon - 1} dyadic_current_checks={dyadic_checks}")
    print(f"periods_w1..30={periods[:30]}")
    print(f"lift_bits_eps1..36={eps_word}")
    print(f"period_exponents_E1..30={period_exponents[:30]}")
    print(f"arrival_atlas_cap30=(innovations,period)={atlas[30]}")
    print(f"anf_full_degrees_n1..8={degrees}")
    print(f"marker_sets_n1..8={markers}")
    print(f"scalar_nonclosure=(Sp,Sq,intersection_parities)={hostile}")
    print(f"kernel_length16_distinct={kernel16} positive_index_length15_distinct={positive_kernel15}")
    print(f"dyadic_length12_distinct={kernel12}")
    print("independent_mirror_center_agreement_t=0..131071")
    print(f"center_sha256_0..262143={hash1}")
    print(f"center_sha256_0..393215={hash2}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
