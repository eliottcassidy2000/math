#!/usr/bin/env python3
"""Exact finite audit for THM-4225's bad-owner successor-rook hierarchy."""

from __future__ import annotations

import hashlib
import itertools
import math


def tournament(n: int, word: int) -> tuple[int, ...]:
    out = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (word >> bit) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            bit += 1
    return tuple(out)


def arc(out: tuple[int, ...], i: int, j: int) -> bool:
    return bool(out[i] & (1 << j))


def bad_owner_counts(out: tuple[int, ...]) -> list[int]:
    n = len(out)
    counts = [0] * (1 << n)
    for tail in itertools.permutations(range(1, n)):
        cycle = (0,) + tail
        mask = 0
        for k, owner in enumerate(cycle):
            successor = cycle[(k + 1) % n]
            if arc(out, successor, owner):
                mask |= 1 << owner
        counts[mask] += 1
    return counts


def has_partial_cycle(domain: tuple[int, ...], image: tuple[int, ...]) -> bool:
    follow = dict(zip(domain, image))
    for start in domain:
        seen: set[int] = set()
        vertex = start
        while vertex in follow:
            if vertex in seen:
                return True
            seen.add(vertex)
            vertex = follow[vertex]
    return False


def successor_rook_count(out: tuple[int, ...], mask: int) -> int:
    n = len(out)
    domain = tuple(i for i in range(n) if (mask >> i) & 1)
    choices = [tuple(j for j in range(n) if arc(out, j, i)) for i in domain]
    count = 0
    for image in itertools.product(*choices):
        if len(set(image)) != len(image):
            continue
        if has_partial_cycle(domain, image):
            continue
        count += 1
    return count


def audit() -> None:
    digest = hashlib.sha256()
    tournament_rows = 0
    subset_rows = 0
    pair_rows = 0
    order_rows: list[tuple[int, int, int]] = []
    controls: dict[str, tuple[int, tuple[int, ...], tuple[int, ...]]] = {}

    for n in range(2, 6):
        words = 1 << (n * (n - 1) // 2)
        local_subsets = 0
        for word in range(words):
            out = tournament(n, word)
            counts = bad_owner_counts(out)
            full = (1 << n) - 1
            assert sum(counts) == math.factorial(n - 1)
            assert counts[0] == counts[full]

            indegree = tuple(n - 1 - out[i].bit_count() for i in range(n))
            layer = tuple(
                sum(counts[mask] for mask in range(1 << n) if mask.bit_count() == k)
                for k in range(n + 1)
            )
            assert layer == tuple(reversed(layer))

            first_moment = sum(mask.bit_count() * counts[mask] for mask in range(1 << n))
            assert first_moment == math.factorial(n) // 2

            for mask in range(1 << n):
                k = mask.bit_count()
                upper = sum(
                    counts[superset]
                    for superset in range(1 << n)
                    if (superset & mask) == mask
                )
                if k < n:
                    rook = successor_rook_count(out, mask)
                    assert upper == math.factorial(n - k - 1) * rook
                subset_rows += 1
                local_subsets += 1

            for i in range(n):
                singleton_upper = sum(
                    counts[mask] for mask in range(1 << n) if (mask >> i) & 1
                )
                assert singleton_upper == math.factorial(n - 2) * indegree[i]

            for i in range(n):
                for j in range(i + 1, n):
                    if n == 2:
                        continue
                    pair_upper = sum(
                        counts[mask]
                        for mask in range(1 << n)
                        if ((mask >> i) & 1) and ((mask >> j) & 1)
                    )
                    common = sum(arc(out, x, i) and arc(out, x, j) for x in range(n))
                    predicted = math.factorial(n - 3) * (
                        indegree[i] * indegree[j] - common
                    )
                    assert pair_upper == predicted
                    pair_rows += 1

            digest.update(f"{n}:{word}:{counts}:{indegree}\n".encode())
            tournament_rows += 1

            if n == 3 and word == 5:  # upper-triangle word 101: C3
                controls["C3"] = (word, tuple(counts), layer)
            if word == words - 1:  # all lower labels beat higher labels
                controls[f"T{n}"] = (word, tuple(counts), layer)

        order_rows.append((n, words, local_subsets))

    print("THM-4225 exact successor-rook audit")
    print("universe=all labelled tournaments of orders 2..5")
    for n, words, local_subsets in order_rows:
        print(f"order={n} tournaments={words} subset_checks={local_subsets}")
    print(
        f"totals tournaments={tournament_rows} upper_zeta_checks={subset_rows} "
        f"pair_formula_checks={pair_rows}"
    )
    for name in ("C3", "T2", "T3", "T4", "T5"):
        word, counts, layer = controls[name]
        width = (len(layer) - 1) * (len(layer) - 2) // 2
        print(
            f"control={name} word={word:0{width}b} "
            f"C_empty={counts[0]} C_full={counts[-1]} layers={layer}"
        )
    print(f"fingerprint_sha256={digest.hexdigest()}")
    print(
        "checks=PASS total_mass,layer_reversal,first_moment," 
        "all_upper_zeta_successor_rook,pair_codegree_formula"
    )


if __name__ == "__main__":
    audit()
