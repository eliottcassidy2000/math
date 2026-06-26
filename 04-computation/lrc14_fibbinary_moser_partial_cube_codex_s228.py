#!/usr/bin/env python3
"""Exact carrier scout for fibbinary, Moser-de Bruijn, and edge sectors.

This script supports the S228 LRC14 forum synthesis.  It keeps the objects
small and exact:

* fibbinary words of length n are vertices of the Fibonacci cube Gamma_n,
  the induced partial cube on binary words with no adjacent 1s.
* Moser-de Bruijn words with m base-4 digits in {0,1} are an m-cube on even
  bit positions, embedded in Gamma_{2m-1}.
* 2,6,12,20,30,42 are ordered edge-sector counts k(k-1)=2*C(k,2), i.e. twice
  the 1-face count of the (k-1)-simplex.

The LRC interpretation is a sidecar rule, not a sequence shortcut.
"""

from __future__ import annotations

from collections import Counter
from itertools import product
from math import comb


def fib(n: int) -> int:
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def triangular(n: int) -> int:
    return n * (n + 1) // 2


def no_adjacent(bits: tuple[int, ...]) -> bool:
    return all(not (bits[i] and bits[i + 1]) for i in range(len(bits) - 1))


def hamming_one(a: tuple[int, ...], b: tuple[int, ...]) -> bool:
    return sum(x != y for x, y in zip(a, b)) == 1


def fibbinary_words(n: int) -> list[tuple[int, ...]]:
    return [bits for bits in product((0, 1), repeat=n) if no_adjacent(bits)]


def moser_words(m: int) -> list[tuple[int, ...]]:
    """Moser support on bit positions 0,2,...,2m-2, returned low-bit first."""
    words = []
    for choices in product((0, 1), repeat=m):
        bits = []
        for bit in choices:
            bits.append(bit)
            bits.append(0)
        words.append(tuple(bits[:-1]))
    return words


def edge_count(words: list[tuple[int, ...]]) -> int:
    total = 0
    for i, a in enumerate(words):
        for b in words[i + 1 :]:
            if hamming_one(a, b):
                total += 1
    return total


def layer_counts(words: list[tuple[int, ...]]) -> list[int]:
    counts = Counter(sum(w) for w in words)
    return [counts[i] for i in range(max(counts) + 1)]


def format_layers(values: list[int]) -> str:
    return "[" + ",".join(str(v) for v in values) + "]"


def main() -> None:
    print("# LRC14 fibbinary/Moser partial-cube and edge-sector scout")
    print()
    print("## Ordered edge-sector / doubled-triangular row")
    print("k  simplex_edges=C(k,2)  ordered_sectors=k(k-1)=2*T_{k-1}")
    for k in range(2, 8):
        undirected = comb(k, 2)
        ordered = k * (k - 1)
        doubled = 2 * triangular(k - 1)
        assert ordered == doubled == 2 * undirected
        print(f"{k:2d} {undirected:21d} {ordered:35d}")
    print()
    print("sequence k=2..7:", [k * (k - 1) for k in range(2, 8)])
    print("origin: ordered 1-faces of the (k-1)-simplex, not a raw LRC proof.")
    print()

    print("## Fibonacci cubes Gamma_n")
    print(
        "n  vertices=F_{n+2}  edges  layers C(n-r+1,r) over independent r-sets"
    )
    for n in range(1, 9):
        words = fibbinary_words(n)
        layers = layer_counts(words)
        formula_layers = [comb(n - r + 1, r) for r in range((n + 1) // 2 + 1)]
        assert len(words) == fib(n + 2)
        assert layers == formula_layers
        print(
            f"{n:2d} {len(words):16d} {edge_count(words):6d} "
            f"{format_layers(layers)}"
        )
    print()

    print("## Moser-de Bruijn even-bit cubes M_m embedded in Gamma_{2m-1}")
    print("m  bit_length  vertices=2^m  cube_edges=m*2^(m-1)  layers C(m,r)")
    for m in range(1, 8):
        words = moser_words(m)
        layers = layer_counts(words)
        expected_layers = [comb(m, r) for r in range(m + 1)]
        bit_length = 2 * m - 1
        assert all(no_adjacent(w) for w in words)
        assert len(words) == 2**m
        assert edge_count(words) == m * (2 ** (m - 1))
        assert layers == expected_layers
        print(
            f"{m:2d} {bit_length:11d} {len(words):13d} "
            f"{edge_count(words):22d} {format_layers(layers)}"
        )
    print()

    print("## Inclusion and transition guardrails")
    for m in range(1, 7):
        moser_count = 2**m
        gamma_count = fib(2 * m + 1)
        print(
            f"M_{m} subset Gamma_{2*m-1}: "
            f"{moser_count}/{gamma_count} vertices retained"
        )
    print()
    print("transition facts:")
    print("- fibbinary is closed under x -> 2x by appending a zero bit.")
    print("- Moser is closed under x -> 4x by appending a base-4 zero digit.")
    print("- Moser is not closed under x -> 2x unless the even/odd bit phase is retained.")
    print()

    print("## Packet sidecar fields")
    for field in [
        "partial_cube_carrier",
        "language_id",
        "native_transition",
        "bit_phase",
        "independence_complex_layer",
        "simplex_face_layer",
        "ordered_edge_sector_count",
        "exact_M",
        "endpoint_owner",
        "route_label",
        "safe_component_certificate",
    ]:
        print(f"- {field}")
    print()

    print("## Tournament Analysis")
    path = [
        "labelled_lrc_packet_sheaf",
        "partial_cube_carry_state",
        "fibbinary_fibonacci_cube",
        "moser_even_bit_cube",
        "simplex_face_layer",
        "ordered_edge_sector_pronic",
        "raw_sequence_membership",
    ]
    print("vertices:", ", ".join(path))
    print(
        "pairwise observable: retains LRC predicate, exact scale, endpoint/topology, "
        "automaton transition, partial-cube coordinate, simplex layer, sector origin"
    )
    print("priority_hamiltonian_path:", " > ".join(path))
    print("score_histogram={0:1,1:1,2:1,3:1,4:1,5:1,6:1}")
    print("directed_3cycles=0 under the packet-priority gauge")
    print(
        "fieldwise warning: fibbinary, Moser, and simplex-face views are not "
        "interchangeable; they preserve different transitions."
    )


if __name__ == "__main__":
    main()
