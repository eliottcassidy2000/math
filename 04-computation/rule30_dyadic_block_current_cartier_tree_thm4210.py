#!/usr/bin/env python3
"""Exact audit for THM-4210: the Rule-30 dyadic current Cartier tree."""

from __future__ import annotations

from itertools import product


Poly = frozenset[int]


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def add(*polynomials: Poly) -> Poly:
    result: set[int] = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return frozenset(result)


def shift(polynomial: Poly, amount: int) -> Poly:
    return frozenset(exponent + amount for exponent in polynomial)


def multiply_q(polynomial: Poly) -> Poly:
    return add(shift(polynomial, -1), polynomial, shift(polynomial, 1))


def q_power_times(power: int, polynomial: Poly) -> Poly:
    for _ in range(power):
        polynomial = multiply_q(polynomial)
    return polynomial


def cartier(polynomial: Poly, parity: int) -> Poly:
    return frozenset(
        (exponent - parity) // 2
        for exponent in polynomial
        if exponent % 2 == parity
    )


def uncartier(even: Poly, odd: Poly) -> Poly:
    return add(
        frozenset(2 * exponent for exponent in even),
        frozenset(2 * exponent + 1 for exponent in odd),
    )


def rule30_step(row: Poly) -> Poly:
    if not row:
        return frozenset()
    return frozenset(
        site
        for site in range(min(row) - 1, max(row) + 2)
        if ((site - 1 in row) ^ ((site in row) or (site + 1 in row)))
    )


def adjacent_current(row: Poly) -> Poly:
    return frozenset(site for site in row if site + 1 in row)


def sample(polynomial: Poly, scale: int) -> Poly:
    return frozenset(
        exponent // scale
        for exponent in polynomial
        if exponent % scale == 0
    )


MAX_TIME = 160
ROWS = [frozenset({0})]
for _ in range(MAX_TIME):
    ROWS.append(rule30_step(ROWS[-1]))
CURRENTS = [adjacent_current(row) for row in ROWS]


def block_source(scale: int, residue: int, block: int) -> Poly:
    result: Poly = frozenset()
    start = residue + block * scale
    for offset in range(scale):
        result = add(
            result,
            q_power_times(scale - 1 - offset, CURRENTS[start + offset]),
        )
    return result


def block_row(scale: int, residue: int, block: int) -> Poly:
    return sample(ROWS[residue + block * scale], scale)


def block_current(scale: int, residue: int, block: int) -> Poly:
    return sample(block_source(scale, residue, block), scale)


def packed_center(length: int) -> tuple[int, ...]:
    packed = 1
    answer = []
    for time in range(length):
        answer.append((packed >> time) & 1)
        packed ^= (packed << 1) | (packed << 2)
    return tuple(answer)


DIRECT_CENTER = tuple(int(0 in row) for row in ROWS)
check(
    DIRECT_CENTER == packed_center(MAX_TIME + 1),
    "direct and packed center paths disagree",
)
print(f"PASS  independent direct/packed center prefix through time {MAX_TIME}")


coarse_checks = 0
for exponent in range(5):
    scale = 1 << exponent
    for residue in range(scale):
        max_block = (MAX_TIME - residue) // scale
        for block in range(max_block):
            current = block_current(scale, residue, block)
            predicted = add(
                multiply_q(block_row(scale, residue, block)), current
            )
            check(
                predicted == block_row(scale, residue, block + 1),
                f"coarse recurrence failed at {(scale, residue, block)}",
            )
            coarse_checks += 1
print(f"PASS  exact coarse Rule-150/current law on {coarse_checks} blocks")


scale_checks = 0
for exponent in range(4):
    scale = 1 << exponent
    for residue in range(scale):
        max_block = max(
            0,
            (MAX_TIME - residue - (3 * scale - 1)) // (2 * scale) + 1,
        )
        for block in range(max_block):
            parent = [
                block_current(scale, residue, 2 * block + offset)
                for offset in range(3)
            ]
            left = cartier(add(multiply_q(parent[0]), parent[1]), 0)
            right = cartier(add(multiply_q(parent[1]), parent[2]), 0)
            check(
                left == block_current(2 * scale, residue, block),
                f"left scale law failed at {(scale, residue, block)}",
            )
            check(
                right == block_current(2 * scale, residue + scale, block),
                f"right scale law failed at {(scale, residue, block)}",
            )
            check(
                cartier(block_row(scale, residue, 2 * block), 0)
                == block_row(2 * scale, residue, block),
                "left row decimation failed",
            )
            check(
                cartier(block_row(scale, residue, 2 * block + 1), 0)
                == block_row(2 * scale, residue + scale, block),
                "right row decimation failed",
            )
            scale_checks += 2
print(f"PASS  exact two-child scale law on {scale_checks} current comparisons")


lossless_checks = 0
parity_checks = 0
for exponent in range(4):
    scale = 1 << exponent
    for residue in range(scale):
        count = min(12, (MAX_TIME - residue) // scale)
        currents = [block_current(scale, residue, block) for block in range(count)]
        rows = [block_row(scale, residue, block) for block in range(count)]
        even_details: list[Poly] = []
        odd_details: list[Poly] = []
        for block in range(count - 1):
            defect = add(multiply_q(currents[block]), currents[block + 1])
            even = cartier(defect, 0)
            odd = cartier(defect, 1)
            even_details.append(even)
            odd_details.append(odd)
            check(uncartier(even, odd) == defect, "Cartier inversion failed")
            if block + 2 < count:
                check(
                    cartier(rows[block + 2], 0)
                    == add(multiply_q(cartier(rows[block], 0)), even),
                    "even-coset current law failed",
                )
                check(
                    cartier(rows[block + 2], 1)
                    == add(multiply_q(cartier(rows[block], 1)), odd),
                    "odd-coset current law failed",
                )
                parity_checks += 2

        rebuilt = [currents[0]]
        for even, odd in zip(even_details, odd_details):
            defect = uncartier(even, odd)
            rebuilt.append(add(defect, multiply_q(rebuilt[-1])))
        check(rebuilt == currents, "lossless current reconstruction failed")
        lossless_checks += count - 1
print(f"PASS  lossless E/O inversion across {lossless_checks} current steps")
print(f"PASS  even/odd spatial-coset laws on {parity_checks} comparisons")


delta_words = []
shell_checks = 0
for exponent in range(7):
    scale = 1 << exponent
    delta = tuple(
        int(0 in block_current(scale, residue, 0))
        for residue in range(scale)
    )
    check(
        delta
        == tuple(DIRECT_CENTER[scale + residue] ^ DIRECT_CENTER[residue]
                 for residue in range(scale)),
        f"dyadic flip word failed at scale {scale}",
    )
    delta_words.append("".join(map(str, delta)))
    signs = tuple(2 * bit - 1 for bit in DIRECT_CENTER)
    for length in range(scale + 1):
        shell = sum(signs[scale + residue] for residue in range(length))
        twisted = sum(
            signs[residue] * ((-1) ** delta[residue])
            for residue in range(length)
        )
        check(shell == twisted, "signed shell identity failed")
        shell_checks += 1

    if exponent < 6:
        for residue in range(scale):
            parent = [block_current(scale, residue, block) for block in range(3)]
            e0 = cartier(add(multiply_q(parent[0]), parent[1]), 0)
            e1 = cartier(add(multiply_q(parent[1]), parent[2]), 0)
            check(
                int(0 in e0)
                == int(0 in block_current(2 * scale, residue, 0)),
                "left child delta readout failed",
            )
            check(
                int(0 in e1)
                == int(0 in block_current(2 * scale, residue + scale, 0)),
                "right child delta readout failed",
            )
print(f"PASS  dyadic flip words and {shell_checks} signed shell prefixes")
for exponent, word in enumerate(delta_words[:6]):
    print(f"      delta_{exponent}={word}")


# The smallest odd-channel kernel witness: nonzero reconstructed currents,
# but an identically zero marked-center response.
length = 14
even_details = [frozenset() for _ in range(length)]
odd_details = [frozenset({0})] + [frozenset() for _ in range(length - 1)]
ambient_currents = [frozenset()]
for even, odd in zip(even_details, odd_details):
    ambient_currents.append(
        add(uncartier(even, odd), multiply_q(ambient_currents[-1]))
    )
ambient_rows = [frozenset()]
for current in ambient_currents:
    ambient_rows.append(add(multiply_q(ambient_rows[-1]), current))
check(ambient_currents[1] == frozenset({1}), "odd hostile J_1 mismatch")
check(ambient_currents[2] == frozenset({0, 1, 2}), "odd hostile J_2 mismatch")
check(all(0 not in row for row in ambient_rows), "odd detail reached center")
print("PASS  minimal O_0=1 hostile changes current but never the marked center")


# Every finite spatial window of the even channel misses a sufficiently remote
# impulse, which reaches the center by a unique endpoint path.
for cutoff in range(9):
    exponent = cutoff + 1
    impulse = frozenset({2 * exponent})
    response = q_power_times(2 * exponent, impulse)
    check(0 in response, "remote even impulse missed its endpoint arrival")
    check(all(abs(site) > cutoff for site in frozenset({exponent})),
          "even impulse was not outside its declared E window")
print("PASS  bounded E coefficient windows fail for cutoffs 0,...,8")


def residue_signature(polynomial: Poly, modulus: int) -> tuple[int, ...]:
    return tuple(
        sum(exponent % modulus == residue for exponent in polynomial) % 2
        for residue in range(modulus)
    )


for modulus in range(1, 13):
    near = frozenset({0})
    far = frozenset({modulus})
    check(
        residue_signature(near, modulus) == residue_signature(far, modulus),
        "cyclic quotient hostile did not alias",
    )
    check((0 in near) != (0 in far), "marked coefficients did not differ")
print("PASS  fixed cyclic-character banks fail for moduli 1,...,12")


# A physical primitive-cubic alias: the full residue packet agrees, but the
# marked coefficient does not.
physical_a = block_current(2, 1, 1)
physical_b = block_current(2, 1, 3)
check(physical_a == frozenset({-1}), "first cubic hostile polynomial changed")
check(
    physical_b == frozenset({-3, 0, 2}),
    "second cubic hostile polynomial changed",
)
check(
    residue_signature(physical_a, 3) == residue_signature(physical_b, 3),
    "physical cubic packets do not agree",
)
check((0 in physical_a) != (0 in physical_b), "physical marks do not differ")
print("PASS  physical cubic alias J^(2,1)_1 vs J^(2,1)_3")


# Exhaust the arbitrary-trace realization on all short center words.  This is
# an ambient driven-current control, not a physical Rule-30 assertion.
trace_checks = 0
delta_zero = frozenset({0})
for tail in product((0, 1), repeat=8):
    bits = (1,) + tail
    rows = [delta_zero if bit else frozenset() for bit in bits]
    drivers = [
        add(rows[index + 1], multiply_q(rows[index]))
        for index in range(len(rows) - 1)
    ]
    rebuilt = [rows[0]]
    for driver in drivers:
        rebuilt.append(add(multiply_q(rebuilt[-1]), driver))
    check(rebuilt == rows, "arbitrary trace realization failed")
    for index in range(len(drivers) - 1):
        defect = add(multiply_q(drivers[index]), drivers[index + 1])
        check(not cartier(defect, 1), "arbitrary trace created odd detail")
    trace_checks += 1
print(f"PASS  all {trace_checks} binary center words of length nine are realized")


# Finite physical-prefix firewall: two ambient continuations preserve the
# exact physical rows and currents through the declared horizon.
horizon = 12
continuation_length = 64
prefix = ROWS[: horizon + 2]
constant_rows = prefix + [frozenset({0})] * continuation_length
thue_morse_rows = prefix + [
    frozenset({0}) if index.bit_count() % 2 else frozenset()
    for index in range(continuation_length)
]


def drivers(rows: list[Poly]) -> list[Poly]:
    return [
        add(rows[index + 1], multiply_q(rows[index]))
        for index in range(len(rows) - 1)
    ]


constant_drivers = drivers(constant_rows)
thue_morse_drivers = drivers(thue_morse_rows)
check(
    constant_drivers[: horizon + 1] == CURRENTS[: horizon + 1]
    == thue_morse_drivers[: horizon + 1],
    "ambient continuations changed the physical current prefix",
)
check(
    all(0 in row for row in constant_rows[horizon + 2 :]),
    "constant continuation is not constant one",
)
tm_signs = [2 * int(0 in row) - 1 for row in thue_morse_rows[horizon + 2 :]]
check(
    all(sum(tm_signs[start : start + (1 << power)]) == 0
        for power in range(1, 6)
        for start in range(0, continuation_length, 1 << power)),
    "Thue-Morse continuation block balance failed",
)
print(f"PASS  opposite ambient tails share the physical current through t={horizon}")


print("PASS  THM-4210 exact audit complete")
