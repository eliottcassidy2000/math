#!/usr/bin/env python3
"""Exact current proof-graph placement of the two THM-4286 surgeries."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
INDEX = 396


def require(condition: bool, message: str) -> None:
    """Optimization-safe audit gate."""
    if not condition:
        raise RuntimeError(message)


def sha_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fnv_words(words: list[int]) -> int:
    state = OFFSET
    for word in words:
        for byte in int(word).to_bytes(8, "little", signed=False):
            state ^= byte
            state = (state * PRIME) & MASK64
    return state


def fnv_rows(rows: list[tuple[int, int]]) -> int:
    return fnv_words([coordinate for row in rows for coordinate in row])


def raw_rows(rows: list[tuple[int, int]]) -> bytes:
    return "".join(f"{q},{r}\n" for q, r in rows).encode("ascii")


def sha_rows(rows: list[tuple[int, int]]) -> str:
    return hashlib.sha256(raw_rows(rows)).hexdigest()


def read_rows(path: Path) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for number, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        fields = line.split(",")
        require(len(fields) == 2, f"bad row at {path}:{number}")
        rows.append((int(fields[0]), int(fields[1])))
    require(rows == sorted(set(rows)), f"rows not unique and sorted in {path}")
    return rows


def read_index_fibre(path: Path) -> list[tuple[int, int]]:
    expected_header = ["q", "r", "inactive_count"] + [f"w{i}" for i in range(7)]
    wanted = [0] * 7
    wanted[INDEX // 64] = 1 << (INDEX % 64)
    fibre: list[tuple[int, int]] = []
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == expected_header, "signature header changed")
        row_count = 0
        for row in reader:
            row_count += 1
            signature = [int(row[f"w{i}"], 16) for i in range(7)]
            require(
                sum(word.bit_count() for word in signature)
                == int(row["inactive_count"]),
                "signature popcount changed",
            )
            if signature == wanted:
                fibre.append((int(row["q"]), int(row["r"])))
    require(row_count == 24_223, "signature universe count changed")
    require(fibre == sorted(set(fibre)), "index-396 fibre changed")
    return fibre


def describe(label: str, rows: list[tuple[int, int]]) -> None:
    print(
        label,
        len(rows),
        f"FNV={fnv_rows(rows):016x}",
        f"SHA256={sha_rows(rows)}",
    )


def emit(path: Path | None, rows: list[tuple[int, int]]) -> None:
    if path is not None:
        path.write_bytes(raw_rows(rows))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--post4282", type=Path, required=True)
    parser.add_argument("--thm4283-union", type=Path, required=True)
    parser.add_argument("--post4283", type=Path, required=True)
    parser.add_argument("--signatures", type=Path, required=True)
    parser.add_argument("--two-mask-net36", type=Path, required=True)
    parser.add_argument("--emit-fibre36", type=Path)
    parser.add_argument("--emit-redundant-union72", type=Path)
    args = parser.parse_args()

    require(
        sha_file(args.post4282)
        == "c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3",
        "post-THM4282 residual bytes changed",
    )
    require(
        sha_file(args.thm4283_union)
        == "c5646e81b3815bdef5168e36bcd76174065ed21339a5d8853d9efddc8fa3efae",
        "current THM-4283 proof-union bytes changed",
    )
    require(
        sha_file(args.post4283)
        == "7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102",
        "current post-THM4283 residual bytes changed",
    )
    require(
        sha_file(args.signatures)
        == "4f0e8da3fdab8bd5a0e14f3b4fa30050602025f63486aa35e0cf03374e6e3832",
        "THM-4282 signature bytes changed",
    )
    require(
        sha_file(args.two_mask_net36)
        == "7a8b049abe73018e650420d5773fae630733ccfb47f7b5d775023affe23220cd",
        "two-mask net36 bytes changed",
    )

    post4282 = read_rows(args.post4282)
    thm4283_union = read_rows(args.thm4283_union)
    post4283 = read_rows(args.post4283)
    fibre36 = read_index_fibre(args.signatures)
    two_mask_net36 = read_rows(args.two_mask_net36)
    require(len(post4282) == 23_373, "post-THM4282 count changed")
    require(len(thm4283_union) == 691, "THM-4283 proof-union count changed")
    require(len(post4283) == 22_682, "post-THM4283 count changed")
    require(len(fibre36) == 36, "index-396 fibre count changed")
    require(len(two_mask_net36) == 36, "two-mask net count changed")
    require(
        set(thm4283_union) <= set(post4282),
        "THM-4283 proof union leaves inherited universe",
    )
    require(set(fibre36) <= set(post4282), "fibre36 leaves inherited universe")
    require(set(two_mask_net36) <= set(post4282),
            "two-mask net36 leaves inherited universe")

    require(
        set(post4282) - set(thm4283_union) == set(post4283),
        "current THM-4283 residual is not the exact set difference",
    )
    require(not set(thm4283_union) & set(post4283), "THM-4283 split overlaps")
    require(not set(fibre36) & set(two_mask_net36), "surgery branches overlap")
    index_overlap = sorted(set(fibre36) & set(thm4283_union))
    two_mask_overlap = sorted(set(two_mask_net36) & set(thm4283_union))
    index_novel = sorted(set(fibre36) - set(thm4283_union))
    two_mask_novel = sorted(set(two_mask_net36) - set(thm4283_union))
    redundant_union72 = sorted(set(fibre36) | set(two_mask_net36))
    require(index_overlap == fibre36, "index-396 branch is not subsumed")
    require(two_mask_overlap == two_mask_net36, "two-mask branch is not subsumed")
    require(not index_novel and not two_mask_novel, "current ledger decrement changed")
    require(len(redundant_union72) == 72, "alternate-node union count changed")
    require(
        fnv_rows(fibre36) == 0x3D92AB45B46A72C0
        and sha_rows(fibre36)
        == "114050c0100fe58793ed32f38f5c9e6bb530854bc74cc46765e4307a877b5fc6",
        "fibre36 identity changed",
    )
    require(
        fnv_rows(thm4283_union) == 0x4B299B49D107A139
        and sha_rows(thm4283_union)
        == "c5646e81b3815bdef5168e36bcd76174065ed21339a5d8853d9efddc8fa3efae",
        "current THM-4283 proof-union identity changed",
    )
    require(
        fnv_rows(post4283) == 0xF7563445F15EFEBF
        and sha_rows(post4283)
        == "7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102",
        "post-THM4283 residual identity changed",
    )
    require(
        fnv_rows(two_mask_net36) == 0x8DA60395E47E11A3
        and sha_rows(two_mask_net36)
        == "7a8b049abe73018e650420d5773fae630733ccfb47f7b5d775023affe23220cd",
        "two-mask net36 identity changed",
    )
    require(
        fnv_rows(redundant_union72) == 0x4F8F4C79707540A6
        and sha_rows(redundant_union72)
        == "253071871ede4041c658ac7705de5283794f1baa230c9df1e22d16e22ac830b3",
        "alternate-node union identity changed",
    )

    maximum = max(r for _, r in post4283)
    top = [row for row in post4283 if row[1] == maximum]
    require(
        maximum == 637 and top == [(100, 637), (294, 637), (520, 637)],
        "current top changed",
    )

    emit(args.emit_fibre36, fibre36)
    emit(args.emit_redundant_union72, redundant_union72)
    print("THM4286_PROOF_GRAPH_CONSEQUENCE_V2")
    describe("FIBRE36", fibre36)
    print("INDEX_OVERLAP_WITH_CURRENT_THM4283", len(index_overlap))
    print("INDEX_NOVEL_OVER_CURRENT_THM4283", len(index_novel))
    describe("TWO_MASK_NET36", two_mask_net36)
    print("TWO_MASK_OVERLAP_WITH_CURRENT_THM4283", len(two_mask_overlap))
    print("TWO_MASK_NOVEL_OVER_CURRENT_THM4283", len(two_mask_novel))
    describe("ALTERNATE_NODE_UNION72", redundant_union72)
    describe("CURRENT_THM4283_PROOF_UNION", thm4283_union)
    describe("CURRENT_RESIDUAL", post4283)
    print("TOP", maximum, top)
    print(
        "VERDICT PASS EXACT SUBSUMPTION; TWO ALTERNATE DECK NODES; "
        "ZERO CURRENT LEDGER DECREMENT; LRC14 OPEN"
    )


if __name__ == "__main__":
    main()
