#!/usr/bin/env python3
"""Verify and canonically combine sharded H6 closed-danger replay output."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


HEADER = "H6_CLOSED_DANGER_UNION_REPLAY"
METHOD = (
    "method=rebuild_closed_danger_union_at_expanded_prefixes"
    "+subtract_closed_candidate_comb"
)
PASS_LINE = "PASS: independent closed-danger-union H6 replay completed"

EXPECTED_NODES = (924, 83881, 8906315, 559202706, 12671505, 53812, 21)
EXPECTED_DEAD = (0, 0, 0, 555565824, 12638291, 53792, 0)
EXPECTED_FULL = (0, 0, 0, 495797163, 0, 0, 0)
EXPECTED_EARLY = (0, 0, 0, 59768661, 12638291, 53792, 20)
EXPECTED_EDGES = 580918240
EXPECTED_COVERS = 1
EXPECTED_LOOSE = 20
EXPECTED_COVER = (
    "COVER root=286 missing=1,3,5,7,9,11, "
    "ordered=1:14,3:16,5:18,7:20,9:22,11:24,"
)

ROOT_RE = re.compile(
    r"^root=(?P<root>[0-9]+) "
    r"missing=(?P<missing>(?:[0-9]+,){6}) "
    r"nodes=(?P<nodes>(?:[0-9]+,){7}) "
    r"dead=(?P<dead>(?:[0-9]+,){7}) "
    r"full_prunes=(?P<full>(?:[0-9]+,){7}) "
    r"early_prunes=(?P<early>(?:[0-9]+,){7}) "
    r"edges=(?P<edges>[0-9]+) covers=(?P<covers>[0-9]+) "
    r"loose=(?P<loose>[0-9]+) "
    r"certificate_sha256=(?P<certificate>[0-9a-f]{64})$"
)
RANGE_RE = re.compile(r"^root_start=([0-9]+) root_limit=([0-9]+)$")


def parse_csv(word: str, expected_length: int) -> tuple[int, ...]:
    fields = word.split(",")
    if fields[-1] != "":
        raise ValueError(f"missing trailing comma: {word}")
    values = tuple(int(field) for field in fields[:-1])
    if len(values) != expected_length:
        raise ValueError(f"wrong count-vector length: {word}")
    return values


def counts_word(values: Sequence[int]) -> str:
    return "".join(f"{value}," for value in values)


@dataclass(frozen=True)
class RootRow:
    root: int
    missing: tuple[int, ...]
    nodes: tuple[int, ...]
    dead: tuple[int, ...]
    full: tuple[int, ...]
    early: tuple[int, ...]
    edges: int
    covers: int
    loose: int
    certificate: str
    raw: str


def parse_root_row(line: str) -> RootRow:
    match = ROOT_RE.fullmatch(line)
    if match is None:
        raise ValueError(f"malformed root row: {line}")
    return RootRow(
        root=int(match["root"]),
        missing=parse_csv(match["missing"], 6),
        nodes=parse_csv(match["nodes"], 7),
        dead=parse_csv(match["dead"], 7),
        full=parse_csv(match["full"], 7),
        early=parse_csv(match["early"], 7),
        edges=int(match["edges"]),
        covers=int(match["covers"]),
        loose=int(match["loose"]),
        certificate=match["certificate"],
        raw=line,
    )


def sum_vectors(rows: Iterable[RootRow], field: str) -> tuple[int, ...]:
    total = [0] * 7
    for row in rows:
        for index, value in enumerate(getattr(row, field)):
            total[index] += value
    return tuple(total)


def manifest(rows: Sequence[RootRow]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update((row.raw + "\n").encode("ascii"))
    return digest.hexdigest()


def make_total_line(rows: Sequence[RootRow]) -> str:
    return (
        f"TOTAL roots={len(rows)} nodes={counts_word(sum_vectors(rows, 'nodes'))} "
        f"dead={counts_word(sum_vectors(rows, 'dead'))} "
        f"full_prunes={counts_word(sum_vectors(rows, 'full'))} "
        f"early_prunes={counts_word(sum_vectors(rows, 'early'))} "
        f"edges={sum(row.edges for row in rows)} "
        f"covers={sum(row.covers for row in rows)} "
        f"loose={sum(row.loose for row in rows)} "
        f"manifest_sha256={manifest(rows)}"
    )


def verify_shard(path: Path) -> tuple[list[RootRow], list[str]]:
    lines = path.read_text(encoding="ascii").splitlines()
    if len(lines) < 6 or lines[0] != HEADER or lines[1] != METHOD:
        raise ValueError(f"{path}: invalid replay header")
    range_match = RANGE_RE.fullmatch(lines[2])
    if range_match is None:
        raise ValueError(f"{path}: invalid root range header")
    root_start, root_limit = map(int, range_match.groups())

    rows = [parse_root_row(line) for line in lines if line.startswith("root=")]
    covers = [line for line in lines if line.startswith("COVER root=")]
    totals = [line for line in lines if line.startswith("TOTAL roots=")]
    passes = [line for line in lines if line == PASS_LINE]
    if [row.root for row in rows] != list(
        range(root_start, root_start + root_limit)
    ):
        raise ValueError(f"{path}: noncanonical or incomplete root range")
    if totals != [make_total_line(rows)]:
        raise ValueError(f"{path}: shard total or manifest mismatch")
    if passes != [PASS_LINE]:
        raise ValueError(f"{path}: missing or repeated PASS line")
    if len(covers) != sum(row.covers for row in rows):
        raise ValueError(f"{path}: cover-row count mismatch")
    return rows, covers


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("shards", nargs="+", type=Path)
    args = parser.parse_args()

    all_rows: list[RootRow] = []
    all_covers: list[str] = []
    for path in args.shards:
        rows, covers = verify_shard(path)
        all_rows.extend(rows)
        all_covers.extend(covers)
    all_rows.sort(key=lambda row: row.root)

    if [row.root for row in all_rows] != list(range(924)):
        raise ValueError("combined input does not contain each root 0..923 once")
    expected_missing = list(itertools.combinations(range(1, 13), 6))
    for row, missing in zip(all_rows, expected_missing, strict=True):
        if row.missing != missing:
            raise ValueError(f"root {row.root}: deletion-set ordering mismatch")

    actual = (
        sum_vectors(all_rows, "nodes"),
        sum_vectors(all_rows, "dead"),
        sum_vectors(all_rows, "full"),
        sum_vectors(all_rows, "early"),
        sum(row.edges for row in all_rows),
        sum(row.covers for row in all_rows),
        sum(row.loose for row in all_rows),
    )
    expected = (
        EXPECTED_NODES,
        EXPECTED_DEAD,
        EXPECTED_FULL,
        EXPECTED_EARLY,
        EXPECTED_EDGES,
        EXPECTED_COVERS,
        EXPECTED_LOOSE,
    )
    if actual != expected:
        raise ValueError(f"all-root frozen census mismatch: {actual!r}")
    if all_covers != [EXPECTED_COVER]:
        raise ValueError(f"sole-cover row mismatch: {all_covers!r}")

    print(HEADER)
    print(METHOD)
    print("root_start=0 root_limit=924")
    cover_by_root = {286: EXPECTED_COVER}
    for row in all_rows:
        print(row.raw)
        if row.root in cover_by_root:
            print(cover_by_root[row.root])
    print(make_total_line(all_rows))
    print(PASS_LINE)


if __name__ == "__main__":
    main()
