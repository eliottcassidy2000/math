#!/usr/bin/env python3
"""Exact typed row-set join of THM-4296 with the independent THM-4295 nodes.

This consumer deliberately joins consequences only.  It does not merge or
identify any of the underlying carrier or common-deck certificates.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import defaultdict
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OURS = ROOT / "05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296"
INCOMING = ROOT / "05-knowledge/results/lrc14_endpoint636_exchange_recursive_ideals_thm4295"

PATHS = {
    "S": OURS / "results/singletons/success_union1219.csv",
    "J": OURS / "results/j19/ideal36.csv",
    "E": OURS / "results/proof_graph/endpoint_success70.csv",
    "C": INCOMING / "inputs/carrier_layers636_633.csv",
    "H19": INCOMING / "inputs/signature19_fibre36.csv",
    "H294": INCOMING / "inputs/signature294_fibre21.csv",
    "H372": INCOMING / "inputs/signature372_fibre54.csv",
}
LIVE_PATH = OURS / "inputs/current_residual22647.csv"

FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def read_pairs(path: Path) -> tuple[tuple[int, int], ...]:
    rows = tuple(
        tuple(map(int, line.split(",")))
        for line in path.read_text(encoding="ascii").splitlines()
        if line and line != "q,r"
    )
    assert rows == tuple(sorted(set(rows))), path
    assert all(0 < q < r for q, r in rows), path
    return rows


def serialize(rows: list[tuple[int, int]]) -> bytes:
    return b"".join(f"{q},{r}\n".encode("ascii") for q, r in rows)


def fnv(rows: list[tuple[int, int]] | tuple[tuple[int, int], ...]) -> str:
    state = FNV_OFFSET
    for q, r in rows:
        for value in (q, r):
            for byte in value.to_bytes(8, "little"):
                state ^= byte
                state = state * FNV_PRIME & MASK64
    return f"{state:016x}"


def identity(label: str, values: set[tuple[int, int]]) -> None:
    rows = sorted(values)
    payload = serialize(rows)
    print(
        f"{label} rows={len(rows)} fnv={fnv(rows)} "
        f"sha256={hashlib.sha256(payload).hexdigest()}"
    )


def show(values: set[tuple[int, int]]) -> str:
    return " ".join(f"({q},{r})" for q, r in sorted(values)) or "EMPTY"


def write_rows(path: Path, values: set[tuple[int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(serialize(sorted(values)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--union-out", type=Path, required=True)
    parser.add_argument("--residual-out", type=Path, required=True)
    parser.add_argument("--new-out", type=Path, required=True)
    args = parser.parse_args()

    live_rows = read_pairs(LIVE_PATH)
    live = set(live_rows)
    nodes = {name: set(read_pairs(path)) for name, path in PATHS.items()}
    assert len(live) == 22647 and fnv(live_rows) == "df5374d4aca67677"
    assert all(values <= live for values in nodes.values())

    print("THM4295_CROSS_PACKET_TYPED_UNION_AUDIT_V1")
    identity("LIVE", live)
    for name, values in nodes.items():
        identity(name, values)

    print("RELATIONS")
    for left, right in combinations(nodes, 2):
        if nodes[left] == nodes[right]:
            print(f"{left} {right} EQUAL")
        elif nodes[left] < nodes[right]:
            print(f"{left} {right} {left}_PROPER_SUBSET_{right}")
        elif nodes[right] < nodes[left]:
            print(f"{left} {right} {right}_PROPER_SUBSET_{left}")

    print("PAIRWISE_INTERSECTIONS")
    for left, right in combinations(nodes, 2):
        overlap = nodes[left] & nodes[right]
        identity(f"{left}_AND_{right}", overlap)
        if overlap:
            print(f"{left}_AND_{right}_ROWS {show(overlap)}")

    profiles: dict[tuple[str, ...], set[tuple[int, int]]] = defaultdict(set)
    union = set().union(*nodes.values())
    for row in union:
        profile = tuple(name for name, values in nodes.items() if row in values)
        profiles[profile].add(row)
    print("EXACT_MEMBERSHIP_PROFILES")
    for profile in sorted(profiles, key=lambda value: (len(value), value)):
        values = profiles[profile]
        label = "+".join(profile)
        identity(f"PROFILE_{label}", values)
        if len(profile) > 1:
            print(f"PROFILE_{label}_ROWS {show(values)}")

    ours = nodes["S"] | nodes["J"] | nodes["E"]
    incoming = nodes["C"] | nodes["H19"] | nodes["H294"] | nodes["H372"]
    added = incoming - ours
    duplicate = incoming & ours
    residual = live - union
    identity("OURS_UNION", ours)
    identity("INCOMING_UNION", incoming)
    identity("INCOMING_ALREADY_COVERED", duplicate)
    identity("INCOMING_STRICTLY_ADDED", added)
    identity("CROSS_PACKET_UNION", union)
    identity("FINAL_COMPLEMENT", residual)

    assert nodes["C"] < nodes["E"]
    assert nodes["H19"] == nodes["J"]
    assert added == {(147, 294), (147, 590), (372, 619)}
    assert len(union) == 1324 and fnv(sorted(union)) == "f55ee025df29bb65"
    assert len(residual) == 21323 and fnv(sorted(residual)) == "09a0dfc4515d556b"
    assert not union & residual and union | residual == live
    max_r = max(r for _, r in residual)
    top = {(q, r) for q, r in residual if r == max_r}
    assert max_r == 626 and top == {(100, 626), (256, 626)}

    write_rows(args.union_out, union)
    write_rows(args.residual_out, residual)
    write_rows(args.new_out, added)
    print(f"FINAL_MAX {max_r} TOP {show(top)}")
    print("SCOPE FINITE_EXACT_TYPED_ROW_SET_UNION_NO_COMMON_DECK_NO_LRC14")


if __name__ == "__main__":
    main()
