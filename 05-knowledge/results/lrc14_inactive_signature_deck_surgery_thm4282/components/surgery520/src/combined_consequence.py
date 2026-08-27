#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_edges(path: Path) -> list[tuple[int, int]]:
    edges: list[tuple[int, int]] = []
    for line in path.read_text().splitlines():
        need(bool(line), f"empty edge row in {path}")
        fields = line.split(",")
        need(len(fields) == 2, f"malformed edge row in {path}")
        edge = (int(fields[0]), int(fields[1]))
        need(0 < edge[0] < edge[1], f"invalid edge in {path}")
        edges.append(edge)
    need(edges == sorted(set(edges)), f"edge ledger not ordered/distinct: {path}")
    return edges


def write_edges(path: Path, edges: list[tuple[int, int]]) -> None:
    path.write_text("".join(f"{q},{r}\n" for q, r in edges))


def fnv_edges(edges: list[tuple[int, int]]) -> int:
    state = 0xCBF29CE484222325
    for edge in edges:
        for value in edge:
            for byte in range(8):
                state ^= (value >> (8 * byte)) & 0xFF
                state = (state * 0x100000001B3) & ((1 << 64) - 1)
    return state


def sha_edges(edges: list[tuple[int, int]]) -> str:
    raw = "".join(f"{q},{r}\n" for q, r in edges).encode()
    return hashlib.sha256(raw).hexdigest()


def read_signature_fibre(path: Path) -> list[tuple[int, int]]:
    rows: list[tuple[tuple[int, int], tuple[int, ...]]] = []
    target: tuple[int, ...] | None = None
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        need(reader.fieldnames == ["q", "r", "inactive_count", "w0", "w1",
                                   "w2", "w3", "w4", "w5", "w6"],
             "full-signature header changed")
        for row in reader:
            edge = (int(row["q"]), int(row["r"]))
            words = tuple(int(row[f"w{i}"], 16) for i in range(7))
            need(sum(word.bit_count() for word in words) ==
                 int(row["inactive_count"]), "signature weight mismatch")
            rows.append((edge, words))
            if edge == (366, 663):
                target = words
    need(len(rows) == 24223 and [edge for edge, _ in rows] ==
         sorted({edge for edge, _ in rows}), "signature universe changed")
    need(target is not None, "missing (366,663) signature")
    bits = [64 * word_index + bit
            for word_index, word in enumerate(target)
            for bit in range(64) if (word >> bit) & 1]
    need(bits == [367], "(366,663) is not the index-367 singleton")
    fibre = [edge for edge, signature in rows if signature == target]
    need(len(fibre) == 26, "index-367 signature fibre size changed")
    return fibre


def describe(name: str, edges: list[tuple[int, int]]) -> str:
    maximum = max((r for _, r in edges), default=0)
    top = [edge for edge in edges if edge[1] == maximum]
    return (f"{name} COUNT {len(edges)} FNV {fnv_edges(edges):016x} "
            f"SHA256 {sha_edges(edges)} MAX_ENDPOINT {maximum} "
            f"TOP " + ";".join(f"{q},{r}" for q, r in top))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signatures", type=Path, required=True)
    parser.add_argument("--old-common", type=Path, required=True)
    parser.add_argument("--new-common", type=Path, required=True)
    parser.add_argument("--post4281", type=Path, required=True)
    parser.add_argument("--emit-520-gain", type=Path, required=True)
    parser.add_argument("--emit-post520", type=Path, required=True)
    parser.add_argument("--emit-367-fibre", type=Path, required=True)
    parser.add_argument("--emit-combined-gain", type=Path, required=True)
    parser.add_argument("--emit-final", type=Path, required=True)
    args = parser.parse_args()

    old_common = read_edges(args.old_common)
    new_common = read_edges(args.new_common)
    post4281 = read_edges(args.post4281)
    need(len(old_common) == 148063 and len(new_common) == 145122 and
         len(post4281) == 24223, "input counts changed")
    old_set = set(old_common)
    new_set = set(new_common)
    post_set = set(post4281)

    canonical_gains = sorted(new_set - old_set)
    gain520 = sorted(new_set & post_set)
    carrier_overlap = sorted(set(canonical_gains) - set(gain520))
    need(len(canonical_gains) == 562 and len(carrier_overlap) == 2 and
         len(gain520) == 560, "520 gain partition changed")
    need(carrier_overlap == [(520, 667), (567, 664)],
         "carrier-overlap pair identity changed")

    fibre367 = read_signature_fibre(args.signatures)
    need(set(fibre367) <= post_set, "index-367 fibre not post-THM4281")
    overlap = sorted(set(gain520) & set(fibre367))
    need(not overlap, "520 gain overlaps index-367 fibre")
    post520 = [edge for edge in post4281 if edge not in set(gain520)]
    need(len(post520) == 23663 and
         [edge for edge in post520 if edge[1] == 663] ==
         [(256, 663), (366, 663)], "520-only residual top changed")
    combined = sorted(set(gain520) | set(fibre367))
    final = [edge for edge in post4281 if edge not in set(combined)]
    need(len(combined) == 586 and len(final) == 23637,
         "combined proof-graph accounting changed")
    need(max(r for _, r in final) == 663 and
         [edge for edge in final if edge[1] == 663] == [(256, 663)],
         "combined residual top changed")

    write_edges(args.emit_520_gain, gain520)
    write_edges(args.emit_post520, post520)
    write_edges(args.emit_367_fibre, fibre367)
    write_edges(args.emit_combined_gain, combined)
    write_edges(args.emit_final, final)

    print("THM4281_SIGNATURE_SURGERY_COMBINED_CONSEQUENCE_V1")
    print("NEW_DECK_COMMON", len(new_common),
          "CANONICAL_RESIDUAL_GAINS", len(canonical_gains),
          "CARRIER_OVERLAP", len(carrier_overlap),
          "POST4281_GAINS", len(gain520))
    print("CARRIER_OVERLAP_ROWS", ";".join(f"{q},{r}" for q, r in carrier_overlap))
    print(describe("GAIN520", gain520))
    print(describe("POST520_RESIDUAL", post520))
    print(describe("FIBRE367", fibre367))
    print("GAIN_INTERSECTION", len(overlap))
    print(describe("COMBINED_GAIN", combined))
    print(describe("FINAL_RESIDUAL", final))
    print("VERDICT PASS EXACT_DISJOINT_UNION_AND_RESIDUAL")


if __name__ == "__main__":
    main()
