#!/usr/bin/env python3
"""Independent exact certificate for the exceptional loose leaf of THM-862.

The full scale-three census has two nonempty depth-six leaves and no covering
leaf.  This small referee reconstructs context 888's literal packet without
using the C++ recursion, intersects all twelve strict 1/13-safe combs with
``Fraction`` endpoints, and recomputes its global maximin from the complete
piecewise-linear candidate set.

Tournament Analysis is deliberately diagnostic.  The vertices are the six
labelled insertion obligations and the pair observable is numerical insertion
order.  This produces a transitive tournament, but it cannot certify looseness:
the exact carrier is the terminal packet together with its safe components.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from lrc13_hamming_five_terminal_endpoint_crosscheck_codex_S10 import (
    exact_maximin,
    measure,
    safe_components,
)


CONTEXT = 888
ROOT_MASK = 2830
MISSING_LABELS = (2, 3, 4, 9, 10, 12)
ORDERS = (3, 3, 3, 3, 3, 3)
UNITS = (1, 1, 1, 1, 1, 2)
# (label, speed), in the unique increasing-speed terminal order.
PATH = ((9, 1), (10, 4), (2, 19), (3, 22), (12, 23), (4, 25))
EXPECTED_PACKET = (1, 3, 4, 15, 18, 19, 21, 22, 23, 24, 25, 33)
EXPECTED_MAXIMUM = F(1, 9)
EXPECTED_WITNESSES = (F(4, 27), F(8, 27), F(19, 27), F(23, 27))
EXPECTED_COMPONENTS = 20
EXPECTED_MEASURE = F(5_068_153, 131_231_100)


def tournament_fingerprint(path: tuple[tuple[int, int], ...]) -> tuple:
    """Fingerprint the insertion-order tournament, including its unique path."""
    vertices = tuple(label for label, _ in path)
    speed = dict(path)
    edge = {
        (left, right)
        for left, right in combinations(vertices, 2)
        if speed[left] < speed[right]
    } | {
        (right, left)
        for left, right in combinations(vertices, 2)
        if speed[right] < speed[left]
    }
    scores = tuple(sorted(sum((vertex, other) in edge for other in vertices) for vertex in vertices))
    cycles = sum(
        ((a, b) in edge and (b, c) in edge and (c, a) in edge)
        or ((b, a) in edge and (c, b) in edge and (a, c) in edge)
        for a, b, c in combinations(vertices, 3)
    )
    return scores, cycles, tuple(label for label, _ in path)


def main() -> None:
    retained = tuple(3 * label for label in range(1, 13) if label not in MISSING_LABELS)
    replacements = tuple(speed for _, speed in PATH)
    packet = tuple(sorted(retained + replacements))

    assert packet == EXPECTED_PACKET
    assert len(packet) == len(set(packet)) == 12
    assert all(order == 3 for order in ORDERS)
    for label, speed in PATH:
        position = MISSING_LABELS.index(label)
        # The removed AP runner is 3*label, so proper replacement preserves
        # that residue modulo 13 (not the unscaled label itself).
        assert speed % 13 == (3 * label) % 13
        assert speed % 3 == UNITS[position]

    components = safe_components(packet)
    assert len(components) == EXPECTED_COMPONENTS
    assert measure(components) == EXPECTED_MEASURE
    assert all(F(0) <= lo < hi <= F(1) for lo, hi in components)

    maximum, witnesses = exact_maximin(packet)
    assert maximum == EXPECTED_MAXIMUM
    assert witnesses == EXPECTED_WITNESSES
    assert maximum > F(1, 13)

    component_lines = tuple(f"{lo}|{hi}|{hi-lo}" for lo, hi in components)
    component_digest = sha256("\n".join(component_lines).encode()).hexdigest()
    scores, cycles, path_word = tournament_fingerprint(PATH)

    print("THM862_SCALE_THREE_CONTEXT_888_EXACT_PACKET_CERTIFICATE")
    print(f"context={CONTEXT} root_mask={ROOT_MASK}")
    print(f"missing_labels={','.join(map(str, MISSING_LABELS))}")
    print("path=" + ",".join(f"{label}:{speed}" for label, speed in PATH))
    print("packet=" + ",".join(map(str, packet)))
    print(f"safe_components={len(components)} safe_measure={measure(components)}")
    print(f"maximin={maximum} witnesses={','.join(map(str, witnesses))}")
    print(f"component_manifest_sha256={component_digest}")
    print("TOURNAMENT_ANALYSIS")
    print("vertices=labelled_insertion_obligations pair_observable=increasing_terminal_speed")
    print(f"score_sequence={scores} directed_cycles={cycles} SCCs=singletons")
    print(f"tie_hamiltonian_path={path_word} hamiltonian_path_count=1 edge_flips=0")
    print("preserves=numerical_insertion_order")
    print("destroys=safe_component_endpoints+widths+maximin;_therefore_non_diagnostic")
    print(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")
    print("PASS: context 888 is independently and strictly 1/13-loose")


if __name__ == "__main__":
    main()
