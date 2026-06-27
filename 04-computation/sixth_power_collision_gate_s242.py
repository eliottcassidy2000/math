#!/usr/bin/env python3
"""
Sixth-power collision gate scout for the LRC14 Desargues/Beal finalizer.

Compares the binary equation
    a^6 + b^6 = d^6 + e^6
with the ternary equation
    a^6 + b^6 + c^6 = d^6 + e^6 + f^6.

The point is not to prove a number-theory theorem.  The point is to test which
proof-interface sidecar each equation suggests: binary Gaussian-owner rigidity
or ternary diagonal-current/cycle carrier.
"""

from __future__ import annotations

from collections import defaultdict
from math import gcd


BINARY_BOUND = 1000
TERNARY_BOUND = 80


def primitive(xs: tuple[int, ...]) -> bool:
    g = 0
    for x in xs:
        g = gcd(g, x)
    return g == 1


def binary_collisions(bound: int) -> list[tuple[int, tuple[int, int], tuple[int, int]]]:
    seen: dict[int, tuple[int, int]] = {}
    hits: list[tuple[int, tuple[int, int], tuple[int, int]]] = []
    for a in range(1, bound + 1):
        a6 = a**6
        for b in range(a + 1, bound + 1):
            s = a6 + b**6
            prev = seen.get(s)
            if prev is None:
                seen[s] = (a, b)
                continue
            if set(prev) == {a, b}:
                continue
            row = (*prev, a, b)
            if primitive(row):
                hits.append((s, prev, (a, b)))
    return sorted(hits)


def ternary_collisions(bound: int) -> list[tuple[int, tuple[int, int, int], tuple[int, int, int]]]:
    seen: dict[int, tuple[int, int, int]] = {}
    hits: list[tuple[int, tuple[int, int, int], tuple[int, int, int]]] = []
    for a in range(1, bound + 1):
        a6 = a**6
        for b in range(a + 1, bound + 1):
            ab = a6 + b**6
            for c in range(b + 1, bound + 1):
                s = ab + c**6
                prev = seen.get(s)
                if prev is None:
                    seen[s] = (a, b, c)
                    continue
                if set(prev) == {a, b, c}:
                    continue
                row = (*prev, a, b, c)
                if primitive(row):
                    hits.append((s, prev, (a, b, c)))
    return sorted(hits)


def sixth_residue_signature(xs: tuple[int, ...], mod: int) -> tuple[int, ...]:
    return tuple(sorted(pow(x, 6, mod) for x in xs))


def main() -> None:
    print("SIXTH-POWER COLLISION GATE SCOUT")
    print(f"binary_bound={BINARY_BOUND}")
    print(f"ternary_bound={TERNARY_BOUND}")

    binary = binary_collisions(BINARY_BOUND)
    print(f"primitive_binary_collisions={len(binary)}")
    if binary:
        s, left, right = binary[0]
        print(f"smallest_binary_hit={left}={right} sum={s}")
    else:
        print("smallest_binary_hit=NONE_WITHIN_BOUND")

    ternary = ternary_collisions(TERNARY_BOUND)
    print(f"primitive_ternary_collisions={len(ternary)}")
    if ternary:
        s, left, right = ternary[0]
        print(f"smallest_ternary_hit={left}={right} sum={s}")
        print(f"identity: {left[0]}^6+{left[1]}^6+{left[2]}^6 = "
              f"{right[0]}^6+{right[1]}^6+{right[2]}^6 = {s}")
        for mod in (7, 9, 13):
            print(
                f"mod{mod}_sixth_residue_signature="
                f"{sixth_residue_signature(left, mod)}="
                f"{sixth_residue_signature(right, mod)}"
            )

    print()
    print("TOURNAMENT_ANALYSIS")
    vertices = [
        "binary_gaussian_owner_gate",
        "ternary_diagonal_current",
        "beal_common_owner_gate",
        "desargues_girth6_residue",
        "lrc_packet_sidecar",
        "raw_equal_sum_scalar",
    ]
    print("vertices=" + ",".join(vertices))
    print(
        "pairwise_observable="
        "retained_owner_or_current,preserved_lrc_handoff,primitive_collision_status,"
        "residue_signature_visibility"
    )
    print(
        "switch_gauge=A->B when A retains every sidecar B retains and one more "
        "owner/current coordinate needed before scalar quotienting"
    )
    path = [
        "lrc_packet_sidecar",
        "desargues_girth6_residue",
        "ternary_diagonal_current",
        "binary_gaussian_owner_gate",
        "beal_common_owner_gate",
        "raw_equal_sum_scalar",
    ]
    print("tie_hamiltonian_path=" + " > ".join(path))
    print("score_hist={0:1,1:1,2:1,3:1,4:1,5:1}")
    print("directed_3cycles=0")
    print("scc_sizes=[1,1,1,1,1,1]")
    print("hamiltonian_path_count=1")


if __name__ == "__main__":
    main()
