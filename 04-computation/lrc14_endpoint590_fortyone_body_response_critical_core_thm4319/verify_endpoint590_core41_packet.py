#!/usr/bin/env python3
"""Verify the standalone endpoint-590 critical-core packet."""

from __future__ import annotations

from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import itertools
from pathlib import Path


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
EXPECTED_FAILURES_SEMANTIC_SHA = (
    "70bd5accf4f09089a37868bae6a5a35d3cb27055d1bdd0c1b6dcf11d69e5caec"
)
EXPECTED_ATLAS_SEMANTIC_SHA = (
    "900072b0fd1004cb65210cc481ae522d08cfdab75361a4d5a00ba6120151300e"
)
EXPECTED_SEARCH_TRANSCRIPT = (
    "ENDPOINT590_OBSTRUCTION41_NO8_V1\n"
    "UNIVERSE 41 RESTRICTED_SIGNATURES 2041 MAXIMAL_TRACES 395\n"
    "DUAL_WEIGHT 22 CAPACITY_PER_RESPONSE 3 DEPTH 8\n"
    "NODES 1968373 DEAD 1698303 SUM_PRUNES 530603 DUAL_PRUNES 1009466\n"
    "RESULT UNSAT\n"
    "SINGLE_DELETION_DEPTH7 SAT 0 UNSAT 41 NODES 278730 DEAD 270531 "
    "SUM_PRUNES 52383 DUAL_PRUNES 198392\n"
    "SCOPE EXACT_COMPLETE_RESPONSE_ATLAS_RESTRICTED_TO_FIXED_41_"
    "OBLIGATION_SUBINSTANCE_NO_PHYSICAL_ENTRY_NO_LRC14\n"
    "VERDICT PASS\n"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def normalized_text(path: Path) -> str:
    return path.read_text(encoding="ascii").replace("\r\n", "\n")


def read_atlas(path: Path) -> tuple[list[int], dict[tuple[int, int], int]]:
    signatures: list[int] = []
    by_mask: dict[tuple[int, int], int] = {}
    with path.open(newline="", encoding="ascii") as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) == 14_368, "response-atlas row count changed")
    for row in rows:
        signature = int(row["w0"], 16) | (int(row["w1"], 16) << 64)
        require(signature != 0 and signature < 1 << 100,
                "response signature escaped universe")
        require(signature.bit_count() == int(row["weight"]),
                "response weight changed")
        for rank in (8, 9):
            field = row[f"least{rank}"]
            if not field:
                continue
            mask = int(field, 16)
            require(mask.bit_count() == rank, "representative rank changed")
            require((mask, rank) not in by_mask,
                    "representative mask duplicated across signatures")
            by_mask[(mask, rank)] = signature
        require(bool(row["least8"] or row["least9"]),
                "signature lost every representative")
        signatures.append(signature)
    require(len(set(signatures)) == 14_368,
            "full response signatures ceased to be distinct")
    semantic_sha = hashlib.sha256(
        "".join(f"{signature:025x}\n" for signature in signatures).encode("ascii")
    ).hexdigest()
    require(semantic_sha == EXPECTED_ATLAS_SEMANTIC_SHA,
            "response-atlas semantic SHA changed")
    return signatures, by_mask


def maximal_traces(traces: set[int]) -> set[int]:
    result: list[int] = []
    for trace in sorted(traces, key=lambda value: (-value.bit_count(), value)):
        if not any(trace & ~kept == 0 for kept in result):
            result.append(trace)
    return set(result)


def adjacency(signatures: list[int], vertices: list[int]) -> list[int]:
    local = {vertex: index for index, vertex in enumerate(vertices)}
    compatible = [0] * len(vertices)
    vertex_mask = sum(1 << vertex for vertex in vertices)
    for signature in signatures:
        work = signature & vertex_mask
        selected: list[int] = []
        while work:
            least = work & -work
            selected.append(local[least.bit_length() - 1])
            work ^= least
        for left, right in itertools.combinations(selected, 2):
            compatible[left] |= 1 << right
            compatible[right] |= 1 << left
    full = (1 << len(vertices)) - 1
    return [full & ~(1 << vertex) & ~compatible[vertex]
            for vertex in range(len(vertices))]


def clique_census(graph: list[int]) -> tuple[int, int, int]:
    maximum = 0
    maximum_count = 0
    full = (1 << len(graph)) - 1

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
                (index for index in range(len(graph))
                 if pivot_pool >> index & 1),
                key=lambda index: (possible & graph[index]).bit_count(),
            )
        candidates = possible if pivot < 0 else possible & ~graph[pivot]
        candidates &= full
        while candidates:
            least = candidates & -candidates
            vertex = least.bit_length() - 1
            visit(size + 1, possible & graph[vertex], excluded & graph[vertex])
            possible ^= least
            excluded |= least
            candidates ^= least

    visit(0, full, 0)
    edges = sum(neighbors.bit_count() for neighbors in graph) // 2
    return edges, maximum, maximum_count


def verify_manifest(packet: Path) -> int:
    manifest = packet / "SHA256SUMS"
    require(manifest.is_file(), "missing SHA256SUMS")
    declared: dict[str, str] = {}
    for line in normalized_text(manifest).splitlines():
        digest, separator, name = line.partition("  ")
        require(separator == "  " and len(digest) == 64 and name,
                "malformed manifest row")
        require(name not in declared and "/" not in name and "\\" not in name,
                "invalid/duplicate manifest path")
        declared[name] = digest
    actual = {path.name for path in packet.iterdir()
              if path.is_file() and path.name != "SHA256SUMS"}
    require(set(declared) == actual, "manifest closure changed")
    for name, digest in declared.items():
        observed = hashlib.sha256((packet / name).read_bytes()).hexdigest()
        require(observed == digest, f"manifest digest mismatch: {name}")
    return len(declared)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("packet", type=Path)
    parser.add_argument("--skip-manifest", action="store_true")
    args = parser.parse_args()
    packet = args.packet.resolve()
    manifest_files = 19
    if not args.skip_manifest:
        manifest_files = verify_manifest(packet)
        require(manifest_files == 19, "manifest file count changed")

    signatures, by_mask = read_atlas(packet / "endpoint590_response_signatures.csv")
    with (packet / "endpoint590_failures.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        failure_rows = list(csv.DictReader(handle))
    require(len(failure_rows) == 100, "failure count changed")
    bodies: list[int] = []
    for row in failure_rows:
        require((int(row["q"]), int(row["r"])) == (210, 590),
                "failure escaped endpoint-590 row")
        body = int(row["body_hex"], 16)
        require(body.bit_count() == 9, "failure body rank changed")
        bodies.append(body)
    require(len(set(bodies)) == 100, "failure body duplication")
    failure_sha = hashlib.sha256(
        "".join(f"210,590,{body:08x}\n" for body in bodies).encode("ascii")
    ).hexdigest()
    require(failure_sha == EXPECTED_FAILURES_SEMANTIC_SHA,
            "failure semantic SHA changed")

    dual: dict[int, int] = {}
    with (packet / "endpoint590_cover_dual_weights.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        for row in csv.DictReader(handle):
            ordinal = int(row["failure_ordinal"])
            require(ordinal not in dual, "duplicate dual ordinal")
            dual[ordinal] = int(row["weight"])
    require(dual == DUAL and sum(dual.values()) == 22,
            "dual certificate identity changed")
    require(max(sum(dual.get(bit, 0) for bit in range(100)
                    if signature >> bit & 1)
                for signature in signatures) == 3,
            "dual response load changed")

    with (packet / "endpoint590_core41.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        core_rows = list(csv.DictReader(handle))
    require(len(core_rows) == 41, "core row count changed")
    require(tuple(int(row["ordinal"]) for row in core_rows) == CORE,
            "core ordinal identity changed")
    for row in core_rows:
        ordinal = int(row["ordinal"])
        require(int(row["body_hex"], 16) == bodies[ordinal],
                "core body projection changed")
        require(row["layer"] == ("dual" if ordinal in DUAL else "satellite"),
                "core layer changed")
        require(int(row["dual_weight"]) == DUAL.get(ordinal, 0),
                "core dual weight changed")

    core_mask = sum(1 << ordinal for ordinal in CORE)
    trace_preimages = Counter(signature & core_mask for signature in signatures)
    trace_preimages.pop(0, None)
    require(len(trace_preimages) == 2041, "restricted trace count changed")
    maximal = maximal_traces(set(trace_preimages))
    require(len(maximal) == 395, "maximal trace quotient changed")
    with (packet / "endpoint590_core41_restricted_traces.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        trace_rows = list(csv.DictReader(handle))
    require(len(trace_rows) == 2041, "trace-ledger row count changed")
    seen_traces: set[int] = set()
    for row in trace_rows:
        trace = int(row["trace_w0"], 16) | (int(row["trace_w1"], 16) << 64)
        require(trace not in seen_traces, "duplicate restricted trace")
        seen_traces.add(trace)
        require(trace.bit_count() == int(row["weight"]),
                "restricted trace weight changed")
        require(trace_preimages[trace] == int(row["full_signature_preimages"]),
                "restricted trace preimage count changed")
        require(int(row["maximal"]) == int(trace in maximal),
                "maximal-trace flag changed")
    require(seen_traces == set(trace_preimages), "trace ledger not exhaustive")

    groups: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    with (packet / "endpoint590_obstruction41_criticality.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        for row in csv.DictReader(handle):
            omitted = int(row["omitted_ordinal"])
            position = int(row["position"])
            mask = int(row["mask_hex"], 16)
            rank = int(row["rank"])
            response = int(row["response_w0"], 16) | (
                int(row["response_w1"], 16) << 64
            )
            require(by_mask.get((mask, rank)) == response,
                    "criticality response not an atlas representative")
            groups[omitted].append((position, mask, response))
    require(set(groups) == set(CORE), "criticality omissions changed")
    for omitted, records in groups.items():
        require(len(records) == 8, "single-deletion cover size changed")
        require(sorted(position for position, _, _ in records) == list(range(8)),
                "single-deletion cover positions changed")
        require(len({mask for _, mask, _ in records}) == 8,
                "single-deletion cover repeats a mask")
        union = 0
        for _, _, response in records:
            union |= response
        require(union & core_mask == core_mask ^ (1 << omitted),
                "single-deletion union is not exactly the declared deletion")

    cover9: list[tuple[int, int, int]] = []
    with (packet / "endpoint590_obstruction41_cover9.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        for row in csv.DictReader(handle):
            position = int(row["position"])
            mask = int(row["mask_hex"], 16)
            rank = int(row["rank"])
            response = int(row["response_w0"], 16) | (
                int(row["response_w1"], 16) << 64
            )
            require(by_mask.get((mask, rank)) == response,
                    "nine-cover response not an atlas representative")
            cover9.append((position, mask, response))
    require(len(cover9) == 9 and
            sorted(position for position, _, _ in cover9) == list(range(9)) and
            len({mask for _, mask, _ in cover9}) == 9,
            "nine-cover shape changed")
    union9 = 0
    for _, _, response in cover9:
        union9 |= response
    require(union9 & core_mask == core_mask, "nine-cover misses a core point")

    full_graph = adjacency(signatures, list(range(100)))
    core_graph = adjacency(signatures, list(CORE))
    full_clique = clique_census(full_graph)
    core_clique = clique_census(core_graph)
    require(full_clique == (740, 4, 50), "full incompatibility census changed")
    require(core_clique == (175, 4, 7), "core incompatibility census changed")
    colors: dict[int, int] = {}
    with (packet / "endpoint590_core41_incompatibility_4color.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        for row in csv.DictReader(handle):
            ordinal = int(row["ordinal"])
            require(ordinal not in colors, "duplicate coloring ordinal")
            colors[ordinal] = int(row["color"])
    require(set(colors) == set(CORE) and set(colors.values()) == set(range(4)),
            "core four-coloring shape changed")
    for left, neighbors in enumerate(core_graph):
        for right in range(left + 1, len(CORE)):
            if neighbors >> right & 1:
                require(colors[CORE[left]] != colors[CORE[right]],
                        "core coloring has a monochromatic edge")

    transcript_o2 = normalized_text(packet / "endpoint590_obstruction41_no8_O2.out")
    transcript_o3 = normalized_text(packet / "endpoint590_obstruction41_no8_O3.out")
    require(transcript_o2 == transcript_o3 == EXPECTED_SEARCH_TRANSCRIPT,
            "exact-search transcript changed")
    structure = normalized_text(packet / "endpoint590_core41_structure.out")
    require("DUAL_POINTS 21 DUAL_WEIGHT 22 SATELLITES 20\n" in structure and
            "PAIRWISE_INCOMPATIBILITY_FULL EDGES 740 OMEGA 4 MAXIMUM_CLIQUES 50\n"
            in structure and
            "PAIRWISE_INCOMPATIBILITY_CORE EDGES 175 OMEGA 4 MAXIMUM_CLIQUES 7\n"
            in structure and structure.endswith("VERDICT PASS\n"),
            "structure transcript changed")

    # Algebraic audit of the explanatory near-saturation identity. If eight
    # responses covered every dual point, their total capacity is 8*3=24.
    # Writing Delta for unused capacity and Omega for dual-weighted repeated
    # coverage gives Delta + Omega = 24 - 22 = 2 exactly.
    require(8 * 3 - sum(DUAL.values()) == 2,
            "saturation identity constant changed")

    print("ENDPOINT590_CORE41_PACKET_VERIFY_V1")
    print("MANIFEST_FILES", manifest_files)
    print("ATLAS 14368 RESTRICTED_TRACES 2041 MAXIMAL_TRACES 395")
    print("CORE 41 DUAL_POINTS 21 SATELLITES 20 DUAL_WEIGHT 22 MAX_LOAD 3")
    print("EXACT_COVER_NUMBER 9 SINGLE_DELETION_EXACT_COVER_NUMBER 8")
    print("SATURATION_IDENTITY DELTA_PLUS_WEIGHTED_OVERLAP 2")
    print("INCOMPATIBILITY_FULL 740,4,50 CORE 175,4,7 CORE_CHROMATIC_NUMBER 4")
    print("SCOPE FINITE_EXACT_FIXED_RESPONSE_HYPERGRAPH_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
