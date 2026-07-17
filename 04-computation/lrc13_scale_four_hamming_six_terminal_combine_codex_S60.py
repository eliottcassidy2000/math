#!/usr/bin/env python3
"""Exact shard combiner for the THM-952 scale-four terminal census.

The combiner requires four contiguous 64-context shards and one fully ungated
single-context replay.  It verifies raw file hashes, every row/summary total,
global context coverage, the recursion tree identity, frozen exact totals,
and equality of proof counts in the gated/ungated replay.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path


EXPECTED_RANGES = ((0, 64), (64, 128), (128, 192), (192, 256))
EXPECTED_SHARD_HASHES = (
    "b931a594fbe88b7654ba3825f5f4b39f9cd25ac0bb9a68422b81af36849ace2b",
    "3f04d203ddd7e6b63aaa9b6dea36715d6f915ff503801a26dfd0ace0ab3bb720",
    "b52f1a926419f22f1089dd68a7fae4e7b9e7371464dc3f6316ee2f7414945e21",
    "b78ad6926545050eb54a353c321341362c443a785ceb4885a83f12651ca513c6",
)
EXPECTED_UNGATED_HASH = (
    "ee3a4241f81a1b97cd450bb35c1dc333c6707fb8e17cb25a11703f5e9b6048e7"
)
EXPECTED_ROW_HASH = (
    "09a46adda0ffce3647075c563ab83e7616a9e892a22d6f172ace591e5ed5b11e"
)
EXPECTED = {
    "nodes": (256, 25_132, 2_577_024, 163_876_444, 496_938, 643, 0),
    "dead": (0, 0, 1, 163_695_372, 496_386, 643, 0),
    "full_tooth": (0, 0, 0, 156_889_649, 0, 0, 0),
    "streaming_cap": (0, 0, 1, 6_805_723, 496_386, 643, 0),
    "candidate_edges": 166_976_181,
    "covers": 0,
    "loose": 0,
    "maximum_cap": 7_665,
    "maximum_candidate_speed": 7_659,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_array(word: str) -> tuple[int, ...]:
    answer = tuple(map(int, word.split(",")))
    require(len(answer) == 7, f"expected seven recursion depths: {word}")
    return answer


@dataclass(frozen=True)
class Stats:
    nodes: tuple[int, ...]
    dead: tuple[int, ...]
    full_tooth: tuple[int, ...]
    streaming_cap: tuple[int, ...]
    candidate_edges: int
    covers: int
    loose: int
    maximum_cap: int
    maximum_candidate_speed: int


@dataclass(frozen=True)
class Row:
    index: int
    root_mask: int
    labels: tuple[int, ...]
    units: tuple[int, ...]
    stats: Stats
    raw: str


@dataclass(frozen=True)
class Shard:
    start: int
    end: int
    early_gates: int
    rows: tuple[Row, ...]
    summary: Stats


def parse_stats_fields(fields: list[str], offset: int) -> Stats:
    return Stats(
        nodes=parse_array(fields[offset]),
        dead=parse_array(fields[offset + 1]),
        full_tooth=parse_array(fields[offset + 2]),
        streaming_cap=parse_array(fields[offset + 3]),
        candidate_edges=int(fields[offset + 4]),
        covers=int(fields[offset + 5]),
        loose=int(fields[offset + 6]),
        maximum_cap=int(fields[offset + 7]),
        maximum_candidate_speed=int(fields[offset + 8]),
    )


def parse_row(line: str) -> Row:
    fields = line.split("|")
    require(len(fields) == 14 and fields[0] == "GENERIC_ROW", "bad row format")
    labels = tuple(map(int, fields[3].split(",")))
    units = tuple(map(int, fields[4].split(",")))
    require(len(labels) == len(units) == 6, "bad row decoration length")
    require(labels == tuple(sorted(labels)) and len(set(labels)) == 6, "bad labels")
    require(set(units) <= {1, 3}, "bad scale-four units")
    expected_mask = sum(1 << (label - 1) for label in labels)
    require(int(fields[2]) == expected_mask, "root-mask / label mismatch")
    return Row(
        index=int(fields[1]),
        root_mask=int(fields[2]),
        labels=labels,
        units=units,
        stats=parse_stats_fields(fields, 5),
        raw=line,
    )


def add_stats(rows: tuple[Row, ...]) -> Stats:
    return Stats(
        nodes=tuple(sum(row.stats.nodes[depth] for row in rows) for depth in range(7)),
        dead=tuple(sum(row.stats.dead[depth] for row in rows) for depth in range(7)),
        full_tooth=tuple(
            sum(row.stats.full_tooth[depth] for row in rows) for depth in range(7)
        ),
        streaming_cap=tuple(
            sum(row.stats.streaming_cap[depth] for row in rows) for depth in range(7)
        ),
        candidate_edges=sum(row.stats.candidate_edges for row in rows),
        covers=sum(row.stats.covers for row in rows),
        loose=sum(row.stats.loose for row in rows),
        maximum_cap=max(row.stats.maximum_cap for row in rows),
        maximum_candidate_speed=max(
            row.stats.maximum_candidate_speed for row in rows
        ),
    )


def parse_summary(line: str) -> Stats:
    fields = line.split("|")
    require(fields[0] == "GENERIC_SHARD_SUMMARY", "bad shard summary tag")
    values = dict(field.split("=", 1) for field in fields[1:])
    expected_keys = {
        "nodes",
        "dead",
        "full_tooth",
        "streaming_cap",
        "candidate_edges",
        "covers",
        "loose",
        "maximum_cap",
        "maximum_candidate_speed",
    }
    require(set(values) == expected_keys, "bad shard summary keys")
    return Stats(
        nodes=parse_array(values["nodes"]),
        dead=parse_array(values["dead"]),
        full_tooth=parse_array(values["full_tooth"]),
        streaming_cap=parse_array(values["streaming_cap"]),
        candidate_edges=int(values["candidate_edges"]),
        covers=int(values["covers"]),
        loose=int(values["loose"]),
        maximum_cap=int(values["maximum_cap"]),
        maximum_candidate_speed=int(values["maximum_candidate_speed"]),
    )


def parse_shard(path: Path) -> Shard:
    lines = path.read_text().splitlines()
    require(
        lines[0] == "SCALE_FOUR_HAMMING_SIX_GENERIC_RECURSION_SHARD",
        f"bad shard heading: {path}",
    )
    arithmetic = dict(field.split("=", 1) for field in lines[1].split())
    require(arithmetic["arithmetic"] == "integer+rational", "nonexact arithmetic")
    require(arithmetic["floating_point"] == "none", "floating point enabled")
    require(arithmetic["height_cutoff"] == "none", "height cutoff enabled")
    require(int(arithmetic["depth_limit"]) == 6, "nonterminal shard")
    meta = dict(field.split("=", 1) for field in lines[2].split())
    start, end = int(meta["context_start"]), int(meta["context_end"])
    require(int(meta["context_count"]) == end - start, "context-count mismatch")
    require(int(meta["all_contexts"]) == 256, "wrong global context count")
    require(int(meta["direct_sheet_checks"]) == 59_136, "wrong sheet-check count")
    require(lines[-1] == "GENERIC_SHARD_DONE", "missing shard done marker")
    summary_lines = [line for line in lines if line.startswith("GENERIC_SHARD_SUMMARY|")]
    require(len(summary_lines) == 1, "wrong shard-summary count")
    rows = tuple(parse_row(line) for line in lines if line.startswith("GENERIC_ROW|"))
    require(tuple(row.index for row in rows) == tuple(range(start, end)), "row gap")
    summary = parse_summary(summary_lines[0])
    require(add_stats(rows) == summary, "row / shard-summary disagreement")
    return Shard(start, end, int(arithmetic["early_gates"]), rows, summary)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("shards", nargs=4, type=Path)
    parser.add_argument("--ungated", required=True, type=Path)
    arguments = parser.parse_args()

    shards = tuple(parse_shard(path) for path in arguments.shards)
    require(
        tuple((shard.start, shard.end) for shard in shards) == EXPECTED_RANGES,
        "shards are not the frozen contiguous partition",
    )
    require(all(shard.early_gates == 1 for shard in shards), "gated shard missing")
    shard_hashes = tuple(digest(path) for path in arguments.shards)
    require(shard_hashes == EXPECTED_SHARD_HASHES, "frozen shard hash mismatch")

    rows = tuple(row for shard in shards for row in shard.rows)
    require(tuple(row.index for row in rows) == tuple(range(256)), "global row gap")
    require(len({(row.labels, row.units) for row in rows}) == 256, "duplicate context")
    row_hash = sha256(("\n".join(row.raw for row in rows) + "\n").encode()).hexdigest()
    require(row_hash == EXPECTED_ROW_HASH, "row payload hash mismatch")
    aggregate = add_stats(rows)
    require(aggregate == Stats(**EXPECTED), "frozen global totals mismatch")
    require(
        aggregate.candidate_edges == sum(aggregate.nodes[1:]),
        "candidate-edge / child-node identity failed",
    )
    require(
        aggregate.dead
        == tuple(
            aggregate.full_tooth[depth] + aggregate.streaming_cap[depth]
            for depth in range(7)
        ),
        "dead nodes are not exactly shortcut-certified nodes",
    )
    require(aggregate.nodes[1] == 25_132, "THM-952 first-layer mismatch")

    ungated_hash = digest(arguments.ungated)
    require(ungated_hash == EXPECTED_UNGATED_HASH, "ungated replay hash mismatch")
    ungated = parse_shard(arguments.ungated)
    require((ungated.start, ungated.end, ungated.early_gates) == (255, 256, 0), "wrong ungated replay")
    gated_255 = rows[255].stats
    ungated_255 = ungated.rows[0].stats
    require(
        (
            gated_255.nodes,
            gated_255.dead,
            gated_255.candidate_edges,
            gated_255.covers,
            gated_255.loose,
            gated_255.maximum_cap,
            gated_255.maximum_candidate_speed,
        )
        == (
            ungated_255.nodes,
            ungated_255.dead,
            ungated_255.candidate_edges,
            ungated_255.covers,
            ungated_255.loose,
            ungated_255.maximum_cap,
            ungated_255.maximum_candidate_speed,
        ),
        "gated / ungated proof counts disagree",
    )
    require(
        ungated_255.full_tooth == (0,) * 7
        and ungated_255.streaming_cap == (0,) * 7,
        "ungated replay unexpectedly used a shortcut",
    )

    repo = Path(__file__).resolve().parents[1]
    terminal_source = repo / "04-computation/lrc13_scale_four_hamming_six_terminal_scout_codex_S60.cpp"
    print("THM-952 SCALE-FOUR HAMMING-SIX EXACT TERMINAL CENSUS")
    print("arithmetic=integer+rational floating_point=none height_cutoff=none")
    print(f"contexts={aggregate.nodes[0]} context_range=[0,256) shards=4")
    print(f"nodes={','.join(map(str, aggregate.nodes))}")
    print(f"dead={','.join(map(str, aggregate.dead))}")
    print(f"full_tooth={','.join(map(str, aggregate.full_tooth))}")
    print(f"streaming_cap={','.join(map(str, aggregate.streaming_cap))}")
    print(f"candidate_edges={aggregate.candidate_edges}")
    print(f"covers={aggregate.covers} loose_terminals={aggregate.loose}")
    print(
        f"maximum_cap={aggregate.maximum_cap} "
        f"maximum_candidate_speed={aggregate.maximum_candidate_speed}"
    )
    print("deepest_nodes=643@depth5 terminal_depth6_nodes=0")
    print("verdict=all_256_context_trees_cap_dead_scale_four_face_closed")
    print(f"terminal_source_sha256={digest(terminal_source)}")
    print(f"combiner_source_sha256={digest(Path(__file__))}")
    for (start, end), shard_hash in zip(EXPECTED_RANGES, shard_hashes):
        print(f"shard_{start}_{end}_sha256={shard_hash}")
    print(f"row_payload_sha256={row_hash}")
    print(f"ungated_context_255_sha256={ungated_hash}")
    print("ungated_context_255_proof_counts_match=yes full_ungated_bank_run=no")


if __name__ == "__main__":
    main()
