#!/usr/bin/env python3
"""Exact signed-C4/H1 sidecar for THM-3537's old-L inertia.

The quartic inertia orbit supplies a cyclic four-sheet carrier.  Its
determinant local system has nontrivial mod-two monodromy; the quadratic
orbit has the same bit, and the two bits cancel in the total degree-nine
permutation representation.  This file compares a one-cut representative of
the quartic class with the independently proved Berggren wall and raw
Fibonacci C4 cochains.  Equality is only in H^1 after vertex switching.
"""

from __future__ import annotations

from hashlib import sha256
import json


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def parity_of_cycle(length: int) -> int:
    return (length - 1) % 2


def coboundary(vertices: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        vertices[index] ^ vertices[(index + 1) % len(vertices)]
        for index in range(len(vertices))
    )


def switch(edges: tuple[int, ...], vertices: tuple[int, ...]) -> tuple[int, ...]:
    delta = coboundary(vertices)
    return tuple(left ^ right for left, right in zip(edges, delta))


def rotate_edge(edge: tuple[int, int], amount: int) -> tuple[int, int]:
    return ((edge[0] + amount) % 4, (edge[1] + amount) % 4)


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def main() -> None:
    # A sign cochain records - as 1 and + as 0.  The Keller quartic orbit has
    # only a canonical total monodromy bit; (1,0,0,0) is the one-cut gauge.
    keller_c4 = (1, 0, 0, 0)
    berggren_wall = (1, 0, 1, 1)      # (-,+,-,-), THM-3536
    fibonacci_raw = (0, 1, 0, 0)      # (+,-,+,+), exact P1(F3) cycle
    require(sum(keller_c4) % 2 == parity_of_cycle(4) == 1, "quartic sign")
    require(sum(berggren_wall) % 2 == 1, "wall H1 class")
    require(sum(fibonacci_raw) % 2 == 1, "Fibonacci H1 class")

    wall_gauge = (0, 0, 0, 1)
    require(switch(keller_c4, wall_gauge) == berggren_wall,
            ("Keller-to-wall switching", switch(keller_c4, wall_gauge)))

    # Reconstruct a switching gauge from raw Fibonacci to the wall rather
    # than importing the earlier bridge's chosen representative.
    fibonacci_to_wall = tuple(left ^ right for left, right in zip(
        fibonacci_raw, berggren_wall
    ))
    gauges = []
    for mask in range(16):
        vertices = tuple((mask >> index) & 1 for index in range(4))
        if coboundary(vertices) == fibonacci_to_wall:
            gauges.append(vertices)
    require(len(gauges) == 2, gauges)  # complementary vertex gauges

    all_cochains = [
        tuple((mask >> index) & 1 for index in range(4)) for mask in range(16)
    ]
    odd_cochains = [edges for edges in all_cochains if sum(edges) % 2 == 1]
    orbit = {switch(keller_c4, vertices) for vertices in all_cochains}
    require(len(odd_cochains) == 8 and orbit == set(odd_cochains),
            (len(odd_cochains), len(orbit)))

    quartic_bit = parity_of_cycle(4)
    quadratic_bit = parity_of_cycle(2)
    fixed_bits = tuple(parity_of_cycle(1) for _ in range(3))
    total_old_l_bit = quartic_bit ^ quadratic_bit
    for bit in fixed_bits:
        total_old_l_bit ^= bit
    require((quartic_bit, quadratic_bit, fixed_bits, total_old_l_bit)
            == (1, 1, (0, 0, 0), 0), "old-L determinant cancellation")
    require(parity_of_cycle(2) == 1, "newest-prime transposition bit")

    successor_edges = {(0, 1), (1, 2), (2, 3), (3, 0)}
    missing_pairs = ({0, 2}, {1, 3})
    completions = []
    invariant = []
    for first in ((0, 2), (2, 0)):
        for second in ((1, 3), (3, 1)):
            arcs = successor_edges | {first, second}
            require(len(arcs) == 6, arcs)
            completions.append((first, second))
            if {rotate_edge(edge, 1) for edge in arcs} == arcs:
                invariant.append((first, second))
    require(len(completions) == 4 and not invariant, (completions, invariant))

    record = {
        "keller_c4_one_cut": keller_c4,
        "berggren_wall": berggren_wall,
        "keller_to_wall_vertex_gauge": wall_gauge,
        "fibonacci_raw": fibonacci_raw,
        "fibonacci_to_wall_gauges": gauges,
        "h1_dimension": 1,
        "nonzero_class_cochain_count": len(odd_cochains),
        "old_L_orbit_bits": (quartic_bit, quadratic_bit, *fixed_bits),
        "old_L_total_bit": total_old_l_bit,
        "newest_prime_bit": parity_of_cycle(2),
        "successor_edges": sorted(successor_edges),
        "missing_pairs": [sorted(pair) for pair in missing_pairs],
        "tournament_completions": completions,
        "rotation_invariant_completions": invariant,
    }
    semantic_sha256 = digest(record)
    print("== Keller inertia / signed C4 / H1 XOR bridge ==")
    print(f"keller_c4_one_cut={keller_c4};berggren_wall={berggren_wall};gauge={wall_gauge}")
    print(f"fibonacci_raw={fibonacci_raw};fibonacci_to_wall_gauges={tuple(gauges)}")
    print(f"H1_C4_F2_dimension=1;nonzero_class_cochains={len(odd_cochains)}")
    print(
        f"old_L_orbit_bits={(quartic_bit, quadratic_bit, *fixed_bits)};"
        f"xor_total={total_old_l_bit};newest_prime_bit={parity_of_cycle(2)}"
    )
    print(
        f"successor_edges={tuple(sorted(successor_edges))};"
        f"missing_pairs={tuple(tuple(sorted(pair)) for pair in missing_pairs)}"
    )
    print(
        f"tournament_completions={tuple(completions)};"
        f"rotation_invariant={tuple(invariant)}"
    )
    print(f"semantic_sha256={semantic_sha256}")
    print(
        "scope=determinant-line H1 cospan only;"
        "no canonical edge gauge,LRC current,D5 flux,or Berggren-Keller identification"
    )
    print("all exact checks passed")


if __name__ == "__main__":
    main()
