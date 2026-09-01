#!/usr/bin/env python3
"""Build deterministic ledgers for the endpoint-590 critical 41-core."""

from __future__ import annotations

from collections import Counter
import csv
import hashlib
import itertools
from pathlib import Path
import sys


CORE = (
    0, 5, 9, 12, 14, 18, 20, 22, 24, 25, 29, 32, 35, 37, 38, 40,
    43, 47, 48, 53, 55, 57, 63, 65, 68, 71, 73, 76, 77, 79, 80, 83,
    84, 88, 89, 90, 94, 95, 96, 97, 99,
)
DUAL = {
    0: 1, 9: 1, 12: 1, 14: 1, 20: 1, 22: 1, 24: 1, 25: 1,
    29: 1, 32: 1, 35: 2, 37: 1, 40: 1, 43: 1, 73: 1, 76: 1,
    80: 1, 83: 1, 95: 1, 97: 1, 99: 1,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_atlas(path: Path) -> list[int]:
    signatures = []
    with path.open(newline="", encoding="ascii") as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) == 14_368, "response-atlas row count changed")
    for row in rows:
        signature = int(row["w0"], 16) | (int(row["w1"], 16) << 64)
        require(signature.bit_count() == int(row["weight"]),
                "response-atlas weight changed")
        require(bool(row["least8"] or row["least9"]),
                "response signature lost all representatives")
        if row["least8"]:
            require(int(row["least8"], 16).bit_count() == 8,
                    "rank-eight representative changed")
        if row["least9"]:
            require(int(row["least9"], 16).bit_count() == 9,
                    "rank-nine representative changed")
        signatures.append(signature)
    require(len(set(signatures)) == len(signatures),
            "duplicate full response signature")
    return signatures


def maximal_traces(traces: set[int]) -> set[int]:
    kept: list[int] = []
    for trace in sorted(traces, key=lambda value: (-value.bit_count(), value)):
        if not any(trace & ~other == 0 for other in kept):
            kept.append(trace)
    return set(kept)


def incompatibility_adjacency(signatures: list[int], vertices: list[int]) -> list[int]:
    local = {vertex: position for position, vertex in enumerate(vertices)}
    compatible = [0] * len(vertices)
    vertex_mask = sum(1 << vertex for vertex in vertices)
    for signature in signatures:
        members = signature & vertex_mask
        selected = []
        while members:
            least = members & -members
            selected.append(local[least.bit_length() - 1])
            members ^= least
        for left, right in itertools.combinations(selected, 2):
            compatible[left] |= 1 << right
            compatible[right] |= 1 << left
    full = (1 << len(vertices)) - 1
    return [full & ~(1 << vertex) & ~compatible[vertex]
            for vertex in range(len(vertices))]


def clique_census(adjacency: list[int]) -> tuple[int, int, int]:
    maximum = 0
    maximum_count = 0

    def visit(size: int, possible: int, excluded: int) -> None:
        nonlocal maximum, maximum_count
        if size + possible.bit_count() < maximum:
            return
        if possible == 0 and excluded == 0:
            if size > maximum:
                maximum = size
                maximum_count = 1
            elif size == maximum:
                maximum_count += 1
            return
        pivot_pool = possible | excluded
        pivot = -1
        if pivot_pool:
            pivot = max(
                (index for index in range(len(adjacency))
                 if pivot_pool >> index & 1),
                key=lambda index: (possible & adjacency[index]).bit_count(),
            )
        candidates = possible if pivot < 0 else possible & ~adjacency[pivot]
        candidates &= (1 << len(adjacency)) - 1
        while candidates:
            least = candidates & -candidates
            vertex = least.bit_length() - 1
            visit(size + 1, possible & adjacency[vertex],
                  excluded & adjacency[vertex])
            possible ^= least
            excluded |= least
            candidates ^= least

    visit(0, (1 << len(adjacency)) - 1, 0)
    edge_count = sum(bits.bit_count() for bits in adjacency) // 2
    return edge_count, maximum, maximum_count


def main() -> None:
    if len(sys.argv) != 6:
        raise RuntimeError(
            "usage: build FAILURES ATLAS CORE_CSV TRACE_CSV STRUCTURE_OUT"
        )
    failures_path, atlas_path, core_path, trace_path, output_path = map(
        Path, sys.argv[1:]
    )
    with failures_path.open(newline="", encoding="ascii") as handle:
        failure_rows = list(csv.DictReader(handle))
    require(len(failure_rows) == 100, "failure count changed")
    bodies = []
    for row in failure_rows:
        require((int(row["q"]), int(row["r"])) == (210, 590),
                "failure row escaped endpoint 590")
        body = int(row["body_hex"], 16)
        require(body.bit_count() == 9, "failure body rank changed")
        bodies.append(body)
    require(len(set(bodies)) == 100, "failure bodies ceased to be distinct")

    signatures = read_atlas(atlas_path)
    core_mask = sum(1 << ordinal for ordinal in CORE)
    trace_preimages = Counter(signature & core_mask for signature in signatures)
    trace_preimages.pop(0, None)
    require(len(trace_preimages) == 2041, "restricted trace count changed")
    maximal = maximal_traces(set(trace_preimages))
    require(len(maximal) == 395, "maximal trace quotient changed")

    with core_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("ordinal", "body_hex", "layer", "dual_weight"))
        for ordinal in CORE:
            writer.writerow((ordinal, f"{bodies[ordinal]:08x}",
                             "dual" if ordinal in DUAL else "satellite",
                             DUAL.get(ordinal, 0)))

    with trace_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("trace_w0", "trace_w1", "weight",
                         "full_signature_preimages", "maximal"))
        for trace in sorted(trace_preimages):
            writer.writerow((f"{trace & ((1 << 64) - 1):016x}",
                             f"{trace >> 64:016x}", trace.bit_count(),
                             trace_preimages[trace], int(trace in maximal)))

    core_bodies = [bodies[ordinal] for ordinal in CORE]
    label_frequency = [sum(body >> label & 1 for body in core_bodies)
                       for label in range(30)]
    intersections = Counter((left & right).bit_count()
                            for left, right in itertools.combinations(core_bodies, 2))
    load_distribution = Counter(
        sum(DUAL.get(bit, 0) for bit in CORE if trace >> bit & 1)
        for trace in maximal
    )
    require(max(load_distribution) == 3 and sum(DUAL.values()) == 22,
            "dual capacity identity changed")
    full_clique = clique_census(incompatibility_adjacency(
        signatures, list(range(100))))
    core_clique = clique_census(incompatibility_adjacency(
        signatures, list(CORE)))
    require(full_clique == (740, 4, 50), "full clique census changed")
    require(core_clique == (175, 4, 7), "core clique census changed")

    lines = [
        "ENDPOINT590_CORE41_STRUCTURE_V1",
        "CORE_ORDINALS " + ",".join(map(str, CORE)),
        f"DUAL_POINTS {len(DUAL)} DUAL_WEIGHT {sum(DUAL.values())} "
        f"SATELLITES {len(CORE) - len(DUAL)}",
        "LABEL_FREQUENCY " + " ".join(
            f"{label}:{count}" for label, count in enumerate(label_frequency)
        ),
        "PAIR_INTERSECTION " + " ".join(
            f"{size}:{count}" for size, count in sorted(intersections.items())
        ),
        "MAXIMAL_TRACE_DUAL_LOAD " + " ".join(
            f"{load}:{count}" for load, count in sorted(load_distribution.items())
        ),
        f"PAIRWISE_INCOMPATIBILITY_FULL EDGES {full_clique[0]} "
        f"OMEGA {full_clique[1]} MAXIMUM_CLIQUES {full_clique[2]}",
        f"PAIRWISE_INCOMPATIBILITY_CORE EDGES {core_clique[0]} "
        f"OMEGA {core_clique[1]} MAXIMUM_CLIQUES {core_clique[2]}",
        "SATURATION_IDENTITY HYPOTHETICAL_8_COVER "
        "DELTA_PLUS_WEIGHTED_OVERLAP_EQUALS_2",
        "MAP RESTRICT_FULL_RESPONSE_SIGNATURE_TO_CORE41_TRACE_THEN_"
        "QUOTIENT_BY_INCLUSION",
        "PRESERVES EXISTENCE_OF_K_RESPONSE_COVER_ON_CORE41",
        "LOSES OFF_CORE59_BEHAVIOR_MASK_MULTIPLICITY_RANK_GEOMETRIC_MARGIN_"
        "AND_DELETION_SAFETY",
        "SCOPE FINITE_EXACT_FIXED_RESPONSE_HYPERGRAPH_ONLY_NO_PHYSICAL_ENTRY_"
        "NO_LRC14",
        "FAILURES_SEMANTIC_SHA256 " + hashlib.sha256(
            "".join(f"210,590,{body:08x}\n" for body in bodies).encode("ascii")
        ).hexdigest(),
        "ATLAS_SIGNATURE_SEMANTIC_SHA256 " + hashlib.sha256(
            "".join(f"{signature:025x}\n" for signature in signatures).encode("ascii")
        ).hexdigest(),
        "VERDICT PASS",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
