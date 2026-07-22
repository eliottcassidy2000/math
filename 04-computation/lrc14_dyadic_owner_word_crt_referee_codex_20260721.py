#!/usr/bin/env python3
"""Exact bitset referee for THM-2066."""

from __future__ import annotations

from itertools import combinations
from math import comb, gcd
import random


def require(condition: bool) -> None:
    if not condition:
        raise AssertionError("exact referee check failed")


def abs_residue(k: int, modulus: int) -> int:
    r = k % modulus
    return min(r, modulus - r)


def clock_tables(N: int, max_speed: int) -> tuple[list[int], list[tuple[int, int]]]:
    safe_masks = [0] * (max_speed + 1)
    for c in range(1, max_speed + 1):
        mask = 0
        for r in range(N):
            if 14 * abs_residue(c * r, N) >= N:
                mask |= 1 << r
        safe_masks[c] = mask

    tail_tables = []
    for z in range(1, 2 * N, 2):
        eligible = 0
        owner = 0
        for r in range(N):
            if 7 * abs_residue(z * r, N) < N:
                eligible |= 1 << r
                nearest = (2 * z * r + N) // (2 * N)
                if nearest & 1:
                    owner |= 1 << r
        tail_tables.append((eligible, owner))
    return safe_masks, tail_tables


def core_packet(C: tuple[int, ...], N: int, safe_masks: list[int]) -> int:
    packet = (1 << N) - 1
    for c in C:
        packet &= safe_masks[c]
    return packet


def word_certificate(C: tuple[int, ...], N: int, tables: tuple[list[int], list[tuple[int, int]]]) -> tuple[bool, int, int]:
    safe_masks, tail_tables = tables
    packet = core_packet(C, N, safe_masks)
    if packet == 0:
        return False, 0, 0
    words = {
        owner & packet
        for eligible, owner in tail_tables
        if packet & ~eligible == 0
    }
    has_complement = any((packet ^ word) in words for word in words)
    return not has_complement, packet.bit_count(), len(words)


def direct_pair_survives(C: tuple[int, ...], N: int, x: int, y: int) -> bool:
    safe_residues = [
        r for r in range(N)
        if all(14 * abs_residue(c * r, N) >= N for c in C)
    ]
    for r in safe_residues:
        for j in (0, 1):
            numerator = r + j * N
            x_danger = 14 * abs_residue(x * numerator, 2 * N) < 2 * N
            y_danger = 14 * abs_residue(y * numerator, 2 * N) < 2 * N
            if not (x_danger or y_danger):
                return False
    return True


def word_pair_survives(C: tuple[int, ...], N: int, x: int, y: int) -> bool:
    A = [r for r in range(N) if all(14 * abs_residue(c * r, N) >= N for c in C)]
    for r in A:
        if 7 * abs_residue(x * r, N) >= N or 7 * abs_residue(y * r, N) >= N:
            return False
        nx = (2 * x * r + N) // (2 * N)
        ny = (2 * y * r + N) // (2 * N)
        if (nx - ny) % 2 == 0:
            return False
    return True


def check_local_equivalence(trials: int) -> int:
    rng = random.Random(2066)
    for _ in range(trials):
        C = tuple(sorted(rng.sample(range(1, 25), rng.randint(2, 11))))
        N = rng.randint(2, 40)
        x = rng.randrange(1, 2 * N, 2)
        y = rng.randrange(1, 2 * N, 2)
        require(word_pair_survives(C, N, x, y) == direct_pair_survives(C, N, x, y))
        require(word_pair_survives(C, N, x + 2 * N, y - 2 * N) == word_pair_survives(C, N, x, y))
    return trials


def check_divisor_transport(trials: int) -> tuple[int, int]:
    rng = random.Random(62066)
    antecedents = 0
    for _ in range(trials):
        N = rng.randint(2, 15)
        M = N * rng.randint(2, 5)
        C = tuple(sorted(rng.sample(range(1, 25), rng.randint(2, 11))))
        x = rng.randrange(1, 2 * M, 2)
        y = rng.randrange(1, 2 * M, 2)
        if word_pair_survives(C, M, x, y):
            antecedents += 1
            require(word_pair_survives(C, N, x % (2 * N), y % (2 * N)))
    return trials, antecedents


def census() -> tuple[int, int, dict[int, int]]:
    clocks = tuple(range(15, 35))
    tables = {N: clock_tables(N, 24) for N in clocks}
    divisor_masks = [0] * 25
    for c in range(1, 25):
        for d in range(2, 15):
            if c % d == 0:
                divisor_masks[c] |= 1 << (d - 2)
    full_divisor_mask = (1 << 13) - 1

    candidates = 0
    histogram: dict[int, int] = {}
    for C in combinations(range(1, 25), 11):
        mask = 0
        for c in C:
            mask |= divisor_masks[c]
        if mask != full_divisor_mask or gcd(*C) != 1:
            continue
        candidates += 1
        for N in clocks:
            certificate, _packet_size, _word_count = word_certificate(C, N, tables[N])
            if certificate:
                histogram[N] = histogram.get(N, 0) + 1
                break
        else:
            raise AssertionError(f"uncertified quotient core: {C}")

    require(comb(24, 11) == 2_496_144)
    require(candidates == 59_880)
    require(sum(histogram.values()) == candidates)
    expected = {
        16: 3502, 17: 10779, 18: 1131, 19: 16079, 20: 1465,
        21: 7885, 22: 2947, 23: 10082, 24: 422, 25: 4597,
        26: 733, 27: 196, 28: 28, 29: 20, 31: 10, 32: 3, 34: 1,
    }
    require(histogram == expected)
    return comb(24, 11), candidates, histogram


def main() -> None:
    print("LRC14 DYADIC OWNER-WORD CRT AUDIT -- exact bitsets")
    print(f"local word/direct lift equivalences: {check_local_equivalence(12000)} PASS")
    trials, antecedents = check_divisor_transport(12000)
    print(f"divisor transports: {trials} trials, {antecedents} surviving antecedents PASS")
    total, candidates, histogram = census()
    print(f"quotient cores: {total} total, {candidates} primitive divisor-complete, {sum(histogram.values())} certified PASS")
    print("certificate histogram:", ", ".join(f"{N}:{histogram[N]}" for N in sorted(histogram)))
    print("PASS")


if __name__ == "__main__":
    main()
