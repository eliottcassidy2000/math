#!/usr/bin/env python3
"""Exact minimum-clock set-cover referee for THM-2068.

The universe is the 59,880 primitive divisor-complete eleven-cores in
{1,...,24}.  Clock N covers a core exactly when THM-2066's eligible binary
owner words contain no complementary pair.  Python integers store the cover
sets as exact bitsets.  Exhaustion after dominance deletion proves the
minimum bank size; no heuristic choice is used in the certificate.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import comb, gcd


CLOCKS = tuple(range(15, 35))
MAX_SPEED = 24
CORE_SIZE = 11


def require(condition: bool, message: str = "exact referee check failed") -> None:
    if not condition:
        raise AssertionError(message)


def abs_residue(k: int, modulus: int) -> int:
    r = k % modulus
    return min(r, modulus - r)


def clock_tables(N: int) -> tuple[list[int], list[tuple[int, int]]]:
    safe_masks = [0] * (MAX_SPEED + 1)
    for c in range(1, MAX_SPEED + 1):
        mask = 0
        for r in range(N):
            if 14 * abs_residue(c * r, N) >= N:
                mask |= 1 << r
        safe_masks[c] = mask

    tail_tables: list[tuple[int, int]] = []
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


def certifies(
    C: tuple[int, ...], N: int,
    tables: tuple[list[int], list[tuple[int, int]]],
) -> bool:
    safe_masks, tail_tables = tables
    packet = (1 << N) - 1
    for c in C:
        packet &= safe_masks[c]
    if packet == 0:
        return False

    words = {
        owner & packet
        for eligible, owner in tail_tables
        if packet & ~eligible == 0
    }
    return not any((packet ^ word) in words for word in words)


def divisor_complete(C: tuple[int, ...]) -> bool:
    return all(any(c % d == 0 for c in C) for d in range(2, 15))


def build_cover_bitsets() -> tuple[list[tuple[int, ...]], list[int], Counter[int]]:
    tables = {N: clock_tables(N) for N in CLOCKS}
    cores: list[tuple[int, ...]] = []
    covers = [0] * len(CLOCKS)
    pattern_histogram: Counter[int] = Counter()

    for C in combinations(range(1, MAX_SPEED + 1), CORE_SIZE):
        if gcd(*C) != 1 or not divisor_complete(C):
            continue
        core_index = len(cores)
        cores.append(C)
        pattern = 0
        for i, N in enumerate(CLOCKS):
            if certifies(C, N, tables[N]):
                covers[i] |= 1 << core_index
                pattern |= 1 << i
        require(pattern != 0, f"THM-2066 coverage failed for {C}")
        pattern_histogram[pattern] += 1

    require(comb(MAX_SPEED, CORE_SIZE) == 2_496_144)
    require(len(cores) == 59_880)
    return cores, covers, pattern_histogram


def undominated_indices(covers: list[int]) -> tuple[list[int], dict[int, int]]:
    """Keep one representative of each cover and delete strict subsets."""
    keep: list[int] = []
    replacement: dict[int, int] = {}
    for i, cover_i in enumerate(covers):
        dominators = [
            j for j, cover_j in enumerate(covers)
            if j != i and cover_i | cover_j == cover_j
        ]
        if not dominators:
            keep.append(i)
            continue
        # Equal sets dominate each other. Keep the least index in an equality
        # class unless a strict dominator is available.
        strict = [j for j in dominators if covers[j] != cover_i]
        if strict:
            replacement[i] = min(strict, key=lambda j: (-covers[j].bit_count(), j))
        else:
            representative = min([i, *dominators])
            if i == representative:
                keep.append(i)
            else:
                replacement[i] = representative

    # Every deleted cover is contained in its recorded replacement. Hence any
    # bank can be dominance-replaced without increasing cardinality.
    for i, j in replacement.items():
        require(covers[i] | covers[j] == covers[j])
    return keep, replacement


def exact_minimum_bank(covers: list[int], active: list[int], core_count: int) -> tuple[tuple[int, ...], list[int]]:
    full = (1 << core_count) - 1
    tested_by_size: list[int] = []
    for size in range(1, len(active) + 1):
        tested = 0
        for bank in combinations(active, size):
            tested += 1
            covered = 0
            for i in bank:
                covered |= covers[i]
            if covered == full:
                tested_by_size.append(tested)
                return bank, tested_by_size
        tested_by_size.append(tested)
    raise AssertionError("the full THM-2066 bank did not cover the universe")


def private_witnesses(
    bank: tuple[int, ...], covers: list[int], cores: list[tuple[int, ...]]
) -> dict[int, tuple[int, ...]]:
    witnesses: dict[int, tuple[int, ...]] = {}
    for i in bank:
        other_cover = 0
        for j in bank:
            if j != i:
                other_cover |= covers[j]
        private = covers[i] & ~other_cover
        require(private != 0, f"chosen clock {CLOCKS[i]} is redundant")
        least_bit = private & -private
        core_index = least_bit.bit_length() - 1
        witnesses[CLOCKS[i]] = cores[core_index]
    return witnesses


def main() -> None:
    print("LRC14 MINIMAL OWNER-CLOCK BANK AUDIT -- exact bitsets")
    cores, covers, patterns = build_cover_bitsets()
    print(f"universe: {comb(MAX_SPEED, CORE_SIZE)} total cores, {len(cores)} primitive divisor-complete")
    print(f"distinct nonempty clock-certificate patterns: {len(patterns)}")
    print("individual coverage:", ", ".join(
        f"{N}:{covers[i].bit_count()}" for i, N in enumerate(CLOCKS)
    ))

    active, replacement = undominated_indices(covers)
    print("undominated clocks:", ",".join(str(CLOCKS[i]) for i in active))
    print("dominated replacements:", ", ".join(
        f"{CLOCKS[i]}->{CLOCKS[j]}" for i, j in sorted(replacement.items())
    ) or "none")

    bank, tested = exact_minimum_bank(covers, active, len(cores))
    bank_clocks = tuple(CLOCKS[i] for i in bank)
    print("exhausted subsets by size:", ", ".join(
        f"{size + 1}:{count}" for size, count in enumerate(tested)
    ))
    print(f"minimum bank ({len(bank_clocks)} clocks):", ",".join(map(str, bank_clocks)))

    witnesses = private_witnesses(bank, covers, cores)
    for N in bank_clocks:
        print(f"private witness for {N}:", ",".join(map(str, witnesses[N])))

    full = (1 << len(cores)) - 1
    bank_cover = 0
    for i in bank:
        bank_cover |= covers[i]
    require(bank_cover == full)
    require(all(pattern & sum(1 << i for i in bank) for pattern in patterns))
    print("PASS")


if __name__ == "__main__":
    main()
