#!/usr/bin/env python3
"""Exact repair/audit for THM-309's Paley cycle-design mechanism.

The old proof invoked a nonexistent 2-transitive automorphism group.  This
companion checks the correct sharply transitive action on unordered pairs and
enumerates simple directed cycles modulo cyclic rotation.  No floating point
or third-party package is used.
"""

from __future__ import annotations

import hashlib
import itertools
from collections import Counter


PRIMES = (7, 11, 19, 23)


def residues(prime):
    return frozenset((value * value) % prime for value in range(1, prime))


def is_arc(source, target, prime, squares):
    return (target - source) % prime in squares


def audit_unordered_pair_action(prime, squares):
    source_pair = frozenset((0, 1))
    images = Counter()
    for multiplier in squares:
        for translation in range(prime):
            image = frozenset(
                (multiplier * vertex + translation) % prime
                for vertex in source_pair
            )
            images[image] += 1
    expected_pairs = prime * (prime - 1) // 2
    assert len(images) == expected_pairs
    assert set(images.values()) == {1}
    return len(images)


def cycle_inventory(prime, length, squares):
    pair_counts = Counter()
    total = 0
    vertices = range(prime)
    for support in itertools.combinations(vertices, length):
        anchor = support[0]
        for tail in itertools.permutations(support[1:]):
            cycle = (anchor,) + tail
            if all(
                is_arc(cycle[index], cycle[(index + 1) % length], prime, squares)
                for index in range(length)
            ):
                total += 1
                for pair in itertools.combinations(support, 2):
                    pair_counts[pair] += 1
    values = tuple(pair_counts.get(pair, 0)
                   for pair in itertools.combinations(vertices, 2))
    assert len(set(values)) == 1
    incidence_lambda = values[0]
    assert total * (length * (length - 1) // 2) == (
        prime * (prime - 1) // 2
    ) * incidence_lambda
    return total, incidence_lambda


def main():
    print("THM-309 PALEY ALL-CYCLE-LENGTH DESIGN REPAIR")
    digest = hashlib.sha256()
    for prime in PRIMES:
        squares = residues(prime)
        orbit_size = audit_unordered_pair_action(prime, squares)
        print(f"p={prime};square_affine_group={orbit_size};unordered_pair_orbit=SHARP")
        lengths = range(3, 8) if prime in (7, 11) else (5,)
        for length in lengths:
            total, incidence_lambda = cycle_inventory(prime, length, squares)
            if length == 5:
                expected_total = (
                    prime * (prime * prime - 1) * (prime - 2) * (prime - 3)
                    // 160
                )
                expected_lambda = (
                    (prime + 1) * (prime - 2) * (prime - 3) // 8
                )
                assert (total, incidence_lambda) == (expected_total, expected_lambda)
            row = (prime, length, total, incidence_lambda)
            digest.update(repr(row).encode("ascii"))
            print(
                f"p={prime};k={length};cycles_mod_rotation={total};"
                f"pair_lambda={incidence_lambda};uniform=PASS"
            )
    print(f"semantic_sha256={digest.hexdigest()}")
    print("scope=multiset_of_simple_directed_cycles;unordered_pairs_not_ordered_pairs")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
