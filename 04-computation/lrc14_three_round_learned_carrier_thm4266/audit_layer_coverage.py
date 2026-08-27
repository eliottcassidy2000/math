#!/usr/bin/env python3
"""Link THM-4266 layer transcripts to raw and current residual universes."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path


FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1

LAYER_RE = re.compile(
    r"^LAYER (\d+) ROWS (\d+) RESISTANT (\d+)(?: |$)"
)
RESISTANT_RE = re.compile(r"^RESISTANT PAIR (\d+),(\d+)(?: |$)")
UPPER_RE = re.compile(
    r"^CASCADE_UNION_BOUNDARY_V1 PAIR (\d+),(\d+) .* "
    r"VERDICT (CLOSED|RESISTANT)$"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fnv_add(state: int, word: int) -> int:
    for shift in range(0, 64, 8):
        state ^= (word >> shift) & 0xFF
        state = (state * FNV_PRIME) & MASK64
    return state


def edge_fnv(edges: list[tuple[int, int]]) -> int:
    state = FNV_OFFSET
    for q, r in edges:
        state = fnv_add(state, q)
        state = fnv_add(state, r)
    return state


def edge_sha(edges: list[tuple[int, int]]) -> str:
    payload = b"".join(f"{q},{r}\n".encode("ascii") for q, r in edges)
    return hashlib.sha256(payload).hexdigest()


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_edges(path: Path) -> list[tuple[int, int]]:
    edges: list[tuple[int, int]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        q_text, r_text = line.split(",")
        edges.append((int(q_text), int(r_text)))
    require(len(edges) == len(set(edges)), f"duplicate edge in {path}")
    return edges


def parse_layer_file(path: Path) -> tuple[
    dict[int, tuple[int, int]], list[tuple[int, int]], str
]:
    text = path.read_text(encoding="ascii")
    layers: dict[int, tuple[int, int]] = {}
    resistant: list[tuple[int, int]] = []
    for line in text.splitlines():
        layer_match = LAYER_RE.match(line)
        if layer_match:
            endpoint, rows, failures = map(int, layer_match.groups())
            require(endpoint not in layers, f"duplicate layer in {path}")
            layers[endpoint] = (rows, failures)
        resistant_match = RESISTANT_RE.match(line)
        if resistant_match:
            resistant.append(tuple(map(int, resistant_match.groups())))
    return layers, resistant, text


def check_layer_block(
    path: Path,
    residual: list[tuple[int, int]],
    endpoints: range,
    expected_resistant: set[tuple[int, int]],
) -> int:
    layers, resistant, _ = parse_layer_file(path)
    expected_endpoints = set(endpoints)
    require(set(layers) == expected_endpoints,
            f"layer coverage mismatch in {path}")
    expected_by_layer = {
        r: sum(1 for _, edge_r in residual if edge_r == r)
        for r in expected_endpoints
    }
    for endpoint, (rows, failure_count) in layers.items():
        require(rows == expected_by_layer[endpoint],
                f"row count mismatch at {endpoint} in {path}")
        require(failure_count == sum(
                    1 for _, r in expected_resistant if r == endpoint),
                f"resistance count mismatch at {endpoint} in {path}")
    require(set(resistant) == expected_resistant,
            f"resistant rows mismatch in {path}")
    return sum(rows for rows, _ in layers.values())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "artifact_root",
        type=Path,
        help=(
            "packet root containing results/, or the canonical result "
            "directory itself"
        ),
    )
    args = parser.parse_args()
    artifact_root = args.artifact_root.resolve()
    packet_results = artifact_root / "results"
    results = packet_results if packet_results.is_dir() else artifact_root
    require(results.is_dir(), f"missing result directory: {results}")

    residual = read_edges(results / "post_thm4256_residual.csv")
    require(len(residual) == 180_991,
            "post-THM-4256 residual count changed")
    require(edge_fnv(residual) == 0x021BF0ED1581657F,
            "post-THM-4256 residual FNV changed")
    require(edge_sha(residual) ==
            "9192c5d73aa5f123ddd10f0115dcaf7231fa518980610042e4cd3f8e73afd44f",
            "post-THM-4256 residual SHA changed")

    high = [edge for edge in residual if edge[1] >= 688]
    require(len(high) == 3_337, "high universe count changed")

    current = read_edges(results / "post_thm4261_4262_residual.csv")
    require(len(current) == 180_622 and
            edge_fnv(current) == 0x0CEF4E2887C8F24E and
            edge_sha(current) ==
            "fa1c5672b0f2cd2490413e9b69a4720bf1dc4eef8aee694c1c73d390aba58e11",
            "post-THM-4261+4262 residual changed")

    upper_windows = [
        (results / "union_boundary_725_734.out", 725, 734),
        (results / "union_boundary_735_744.out", 735, 744),
        (results / "union_boundary_745_749.out", 745, 749),
        (results / "union_boundary_750_754.out", 750, 754),
    ]
    upper_status: dict[tuple[int, int], str] = {}
    for path, lower, upper in upper_windows:
        window_edges: list[tuple[int, int]] = []
        for line in path.read_text(encoding="ascii").splitlines():
            match = UPPER_RE.match(line)
            require(match is not None, f"bad upper row in {path}")
            require(" UNION 4675 " in line and
                    " BODIES 14307150 " in line,
                    f"upper carrier/body universe mismatch in {path}")
            q_text, r_text, verdict = match.groups()
            edge = (int(q_text), int(r_text))
            require(edge not in upper_status, "duplicate upper edge")
            upper_status[edge] = verdict
            window_edges.append(edge)
        require(window_edges == [
                    edge for edge in residual if lower <= edge[1] <= upper
                ], f"upper inherited order mismatch in {path}")
    expected_upper = {edge for edge in residual if 725 <= edge[1] <= 754}
    require(set(upper_status) == expected_upper and len(expected_upper) == 498,
            "upper transcript universe mismatch")
    upper_resistant = {
        edge for edge, verdict in upper_status.items()
        if verdict == "RESISTANT"
    }
    require(upper_resistant == {(542, 732)},
            "upper base resistance changed")

    layer_732, resistant_732, text_732 = parse_layer_file(
        results / "round2_layer_732.out")
    require(layer_732 == {732: (24, 0)} and not resistant_732 and
            "STOP FIRST_RESISTANT_LAYER -1 VERDICT RANGE_CLOSED" in text_732,
            "round-two 732 closure changed")

    gap_rows = check_layer_block(
        results / "round3_coverage_724_705.out", residual,
        range(705, 725), set())
    require(gap_rows == 1_088, "724--705 row count changed")

    base_upper_rows = check_layer_block(
        results / "base_descent_724_718.out", residual,
        range(718, 725), set())
    require(base_upper_rows == 288, "base 724--718 row count changed")
    base_lower_rows = check_layer_block(
        results / "base_descent_717_704.out", residual,
        range(704, 718), {(416, 704)})
    require(base_lower_rows == 884, "base 717--704 row count changed")

    seed_rows = check_layer_block(
        results / "seed_layer_704_prefix_compare.out", residual,
        range(704, 705), set())
    require(seed_rows == 84, "round-one seed layer count changed")

    round1_rows = check_layer_block(
        results / "prefix_descent_703_680.out", residual,
        range(700, 704), {(416, 700), (520, 700)})
    require(round1_rows == 356, "round-one descent row count changed")

    round2_rows = check_layer_block(
        results / "round2_descent_699_650.out", residual,
        range(694, 700), {(384, 694)})
    require(round2_rows == 612, "round-two descent row count changed")
    _, _, round2_text = parse_layer_file(
        results / "round2_descent_699_650.out")
    require("HOSTILE_416_700 OLD_FAILURES 38 NEW_FAILURES 0" in round2_text and
            "HOSTILE_520_700 OLD_FAILURES 1 NEW_FAILURES 0" in round2_text,
            "round-two seed closures changed")

    round3_rows = check_layer_block(
        results / "round3_descent_693_650.out", residual,
        range(688, 694), {(520, 688)})
    require(round3_rows == 699, "round-three descent row count changed")
    _, _, round3_text = parse_layer_file(
        results / "round3_descent_693_650.out")
    require("HOSTILE_384_694 OLD_FAILURES 1 NEW_FAILURES 0" in round3_text,
            "round-three seed closure changed")

    retained = [(520, 688)]
    closed = [edge for edge in high if edge not in set(retained)]
    require(len(closed) == 3_336 and
            edge_fnv(closed) == 0xC1BA162D8A364FB3 and
            edge_sha(closed) ==
            "95a9c0eb185847f1d64d949cef2b1343e85701dcf416c27f3790c31a91c40854",
            "linked carrier deletion ledger changed")

    prior_set = set(residual) - set(current)
    require(len(prior_set) == 369, "proved-prior deletion count changed")
    overlap = [edge for edge in closed if edge in prior_set]
    require(len(overlap) == 299, "raw/prior overlap changed")
    novel = read_edges(results / "thm4266_novel_edges.csv")
    require(novel == [edge for edge in closed if edge not in prior_set] and
            len(novel) == 3_037 and
            edge_fnv(novel) == 0x24B36D7047589076 and
            edge_sha(novel) ==
            "fcfec867819898ec1a0e1072f27747aec29b6785794328a153bc2b85956ba112",
            "THM-4266 novel deletion changed")
    current_high = [edge for edge in current if edge[1] >= 688]
    require(len(current_high) == 3_038 and
            edge_fnv(current_high) == 0x4ED8CEB63E3EDC9E and
            edge_sha(current_high) ==
            "4333d7306bb1f8df464b3bd3261559e9d91e6dc88c835e37518206b8a9d0e643",
            "THM-4266 current high universe changed")
    final = read_edges(results / "thm4266_final_residual.csv")
    require(final == [edge for edge in current if edge not in set(novel)] and
            len(final) == 177_585 and
            edge_fnv(final) == 0x6CE05D05EB01DAED and
            edge_sha(final) ==
            "009614651bb81e9763b2a9ff4b580497bfb6978a6c69d18cf986346e369374d9",
            "THM-4266 final residual changed")

    print("LRC14_THM4266_RAW_LAYER_PARTITION_V1")
    print(
        f"POST_THM4256 COUNT {len(residual)} FNV {edge_fnv(residual):016x} "
        f"SHA256 {edge_sha(residual)}"
    )
    print(
        "BLOCK R_754_725 ROWS 498 BASE_RESISTANT 542,732 "
        "ROUND2_CLOSED 1"
    )
    print(
        f"BLOCK R_724_705 ROWS {gap_rows} RESISTANT 0 OUTPUT_SHA256 "
        f"{file_sha(results / 'round3_coverage_724_705.out')}"
    )
    print(
        f"BASE_DESCENT R_724_718 ROWS {base_upper_rows} RESISTANT 0 "
        f"OUTPUT_SHA256 {file_sha(results / 'base_descent_724_718.out')}"
    )
    print(
        f"BASE_DESCENT R_717_704 ROWS {base_lower_rows} RETAINED 416,704 "
        f"OUTPUT_SHA256 {file_sha(results / 'base_descent_717_704.out')}"
    )
    print(
        f"BLOCK R_704 ROWS {seed_rows} RESISTANT 0 OUTPUT_SHA256 "
        f"{file_sha(results / 'seed_layer_704_prefix_compare.out')}"
    )
    print(
        f"BLOCK R_703_700 ROWS {round1_rows} RETAINED_OLD 416,700 520,700 "
        f"OUTPUT_SHA256 {file_sha(results / 'prefix_descent_703_680.out')}"
    )
    print(
        f"BLOCK R_699_694 ROWS {round2_rows} RETAINED_OLD 384,694 "
        f"OUTPUT_SHA256 {file_sha(results / 'round2_descent_699_650.out')}"
    )
    print(
        f"BLOCK R_693_688 ROWS {round3_rows} RETAINED 520,688 "
        f"OUTPUT_SHA256 {file_sha(results / 'round3_descent_693_650.out')}"
    )
    print(
        f"CARRIER_REMOVED COUNT {len(closed)} FNV {edge_fnv(closed):016x} "
        f"SHA256 {edge_sha(closed)}"
    )
    print(
        f"RAW_OVERLAP_PROVED_PRIOR COUNT {len(overlap)} "
        "THM4261 297 THM4262 2"
    )
    print(
        f"POST_THM4261_4262_HIGH COUNT {len(current_high)} "
        f"FNV {edge_fnv(current_high):016x} SHA256 {edge_sha(current_high)}"
    )
    print(
        f"THM4266_NOVEL_REMOVED COUNT {len(novel)} "
        f"FNV {edge_fnv(novel):016x} SHA256 {edge_sha(novel)}"
    )
    print(
        f"THM4266_FINAL COUNT {len(final)} FNV {edge_fnv(final):016x} "
        f"SHA256 {edge_sha(final)} TOP 520,688"
    )
    print("VERDICT PASS EXACT_THM4266_LAYER_PARTITION LRC14_OPEN")


if __name__ == "__main__":
    main()
