#!/usr/bin/env python3
"""No-import scalar audit for THM-4210's Cartier-current identities."""

from __future__ import annotations


Poly = dict[int, int]


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def toggle(target: Poly, exponent: int, bit: int = 1) -> None:
    if not bit:
        return
    target[exponent] = target.get(exponent, 0) ^ 1
    if not target[exponent]:
        del target[exponent]


def plus(*items: Poly) -> Poly:
    answer: Poly = {}
    for item in items:
        for exponent, bit in item.items():
            toggle(answer, exponent, bit)
    return answer


def q_times(item: Poly) -> Poly:
    answer: Poly = {}
    for exponent, bit in item.items():
        for displacement in (-1, 0, 1):
            toggle(answer, exponent + displacement, bit)
    return answer


def q_power_times(power: int, item: Poly) -> Poly:
    for _ in range(power):
        item = q_times(item)
    return item


def c_even(item: Poly) -> Poly:
    return {
        exponent // 2: bit
        for exponent, bit in item.items()
        if exponent % 2 == 0
    }


def c_odd(item: Poly) -> Poly:
    return {
        (exponent - 1) // 2: bit
        for exponent, bit in item.items()
        if exponent % 2 == 1
    }


def join(even: Poly, odd: Poly) -> Poly:
    answer: Poly = {}
    for exponent, bit in even.items():
        toggle(answer, 2 * exponent, bit)
    for exponent, bit in odd.items():
        toggle(answer, 2 * exponent + 1, bit)
    return answer


def rule30_step(row: Poly) -> Poly:
    radius = max((abs(exponent) for exponent in row), default=0) + 1
    answer: Poly = {}
    for site in range(-radius, radius + 1):
        left = row.get(site - 1, 0)
        center = row.get(site, 0)
        right = row.get(site + 1, 0)
        value = left ^ (center | right)
        if value:
            answer[site] = 1
    return answer


MAX_TIME = 48
rows: list[Poly] = [{0: 1}]
for _ in range(MAX_TIME):
    rows.append(rule30_step(rows[-1]))
currents = [
    {
        site: 1
        for site in row
        if row.get(site + 1, 0)
    }
    for row in rows
]


def sampled(item: Poly, scale: int) -> Poly:
    return {
        exponent // scale: bit
        for exponent, bit in item.items()
        if exponent % scale == 0
    }


def source(scale: int, residue: int, block: int) -> Poly:
    answer: Poly = {}
    start = residue + block * scale
    for offset in range(scale):
        answer = plus(
            answer,
            q_power_times(scale - 1 - offset, currents[start + offset]),
        )
    return answer


def coarse_row(scale: int, residue: int, block: int) -> Poly:
    return sampled(rows[residue + block * scale], scale)


def coarse_current(scale: int, residue: int, block: int) -> Poly:
    return sampled(source(scale, residue, block), scale)


coarse_count = 0
scale_count = 0
parity_count = 0
for scale in (1, 2, 4):
    for residue in range(scale):
        block = 0
        while residue + (block + 1) * scale <= MAX_TIME:
            y0 = coarse_row(scale, residue, block)
            y1 = coarse_row(scale, residue, block + 1)
            j0 = coarse_current(scale, residue, block)
            check(y1 == plus(q_times(y0), j0), "coarse scalar check failed")
            coarse_count += 1

            if residue + (block + 2) * scale <= MAX_TIME:
                y2 = coarse_row(scale, residue, block + 2)
                j1 = coarse_current(scale, residue, block + 1)
                defect = plus(q_times(j0), j1)
                even = c_even(defect)
                odd = c_odd(defect)
                check(join(even, odd) == defect, "scalar Cartier inverse failed")
                check(
                    c_even(y2) == plus(q_times(c_even(y0)), even),
                    "scalar even-coset law failed",
                )
                check(
                    c_odd(y2) == plus(q_times(c_odd(y0)), odd),
                    "scalar odd-coset law failed",
                )
                parity_count += 2

            if residue + (2 * block + 3) * scale - 1 <= MAX_TIME:
                j_a = coarse_current(scale, residue, 2 * block)
                j_b = coarse_current(scale, residue, 2 * block + 1)
                j_c = coarse_current(scale, residue, 2 * block + 2)
                check(
                    c_even(plus(q_times(j_a), j_b))
                    == coarse_current(2 * scale, residue, block),
                    "scalar left child law failed",
                )
                check(
                    c_even(plus(q_times(j_b), j_c))
                    == coarse_current(2 * scale, residue + scale, block),
                    "scalar right child law failed",
                )
                scale_count += 2
            block += 1

print(f"PASS  scalar coarse law on {coarse_count} blocks")
print(f"PASS  scalar child law on {scale_count} comparisons")
print(f"PASS  scalar parity law on {parity_count} comparisons")


# Physical odd-channel contrast from two different nodes.
def node_packet(scale: int, residue: int) -> tuple[Poly, Poly, Poly, tuple[int, ...]]:
    j0 = coarse_current(scale, residue, 0)
    j1 = coarse_current(scale, residue, 1)
    defect = plus(q_times(j0), j1)
    profile = tuple(
        rows[residue + block * scale].get(0, 0)
        for block in range(6)
    )
    return j0, c_even(defect), c_odd(defect), profile


packet_a = node_packet(2, 0)
packet_b = node_packet(8, 3)
check(packet_a[0] == packet_b[0] == {0: 1}, "physical J_0 contrast changed")
check(packet_a[1] == packet_b[1] == {}, "physical E_0 contrast changed")
check(packet_a[2] == {} and packet_b[2] == {0: 1}, "physical O_0 contrast changed")
check(packet_a[3] == packet_b[3] == (1, 0, 1, 0, 1, 0),
      "physical center-profile control changed")
print("PASS  physical O contrast (2,0) vs (8,3) with center prefix 101010")


# Primitive cubic alias on the physical block current.
cubic_a = coarse_current(2, 1, 1)
cubic_b = coarse_current(2, 1, 3)


def residues(item: Poly, modulus: int) -> tuple[int, ...]:
    answer = [0] * modulus
    for exponent, bit in item.items():
        answer[exponent % modulus] ^= bit
    return tuple(answer)


check(cubic_a == {-1: 1}, "scalar cubic current A changed")
check(cubic_b == {-3: 1, 0: 1, 2: 1}, "scalar cubic current B changed")
check(residues(cubic_a, 3) == residues(cubic_b, 3), "cubic alias failed")
check(cubic_a.get(0, 0) != cubic_b.get(0, 0), "cubic marks agree")
print("PASS  scalar physical cubic alias with opposite marked coefficient")


# Endpoint response and the odd-detail center kernel.
for distance in range(1, 13):
    response = q_power_times(distance, {distance: 1})
    check(response.get(0, 0) == 1, "endpoint response failed")
    for earlier in range(distance):
        check(
            q_power_times(earlier, {distance: 1}).get(0, 0) == 0,
            "remote impulse arrived too early",
        )
print("PASS  distinct endpoint response times through distance 12")


even_details = [{} for _ in range(12)]
odd_details = [{0: 1}] + [{} for _ in range(11)]
j_values: list[Poly] = [{}]
y_values: list[Poly] = [{}]
for even, odd in zip(even_details, odd_details):
    defect = join(even, odd)
    j_values.append(plus(q_times(j_values[-1]), defect))
    y_values.append(plus(q_times(y_values[-1]), j_values[-2]))
check(any(j_values), "odd kernel current stayed zero")
check(all(value.get(0, 0) == 0 for value in y_values), "odd kernel hit center")
print("PASS  scalar odd-detail kernel has identically zero center")
print("PASS  THM-4210 independent audit complete")
