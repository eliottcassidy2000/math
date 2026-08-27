#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from pathlib import Path


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_edges(path: Path) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for line in path.read_text().splitlines():
        fields = line.split(",")
        need(len(fields) == 2, f"malformed edge row in {path}")
        edge = (int(fields[0]), int(fields[1]))
        need(0 < edge[0] < edge[1], f"invalid edge in {path}")
        rows.append(edge)
    need(rows == sorted(set(rows)), f"edge ledger not ordered/distinct: {path}")
    return rows


def write_edges(path: Path, rows: set[tuple[int, int]]) -> None:
    path.write_text("".join(f"{q},{r}\n" for q, r in sorted(rows)))


def fnv(rows: set[tuple[int, int]]) -> int:
    state = 0xCBF29CE484222325
    for q, r in sorted(rows):
        for value in (q, r):
            for byte in value.to_bytes(8, "little"):
                state ^= byte
                state = (state * 0x100000001B3) & ((1 << 64) - 1)
    return state


def sha(rows: set[tuple[int, int]]) -> str:
    raw = "".join(f"{q},{r}\n" for q, r in sorted(rows)).encode()
    return hashlib.sha256(raw).hexdigest()


def describe(name: str, rows: set[tuple[int, int]]) -> str:
    maximum = max((r for _, r in rows), default=0)
    top = sorted(edge for edge in rows if edge[1] == maximum)
    return (
        f"{name} COUNT {len(rows)} FNV {fnv(rows):016x} SHA256 {sha(rows)} "
        f"MAX_ENDPOINT {maximum} TOP "
        + ";".join(f"{q},{r}" for q, r in top)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--post4281", type=Path, required=True)
    parser.add_argument("--combined586", type=Path, required=True)
    parser.add_argument("--postcombined", type=Path, required=True)
    parser.add_argument("--surgery256", type=Path, required=True)
    parser.add_argument("--emit-carrier", type=Path, required=True)
    parser.add_argument("--emit-surgery256-net", type=Path, required=True)
    parser.add_argument("--emit-union", type=Path, required=True)
    parser.add_argument("--emit-final", type=Path, required=True)
    args = parser.parse_args()

    post4281 = set(read_edges(args.post4281))
    combined586 = set(read_edges(args.combined586))
    postcombined = set(read_edges(args.postcombined))
    surgery256 = set(read_edges(args.surgery256))
    need(len(post4281) == 24223 and len(combined586) == 586 and
         len(postcombined) == 23637 and len(surgery256) == 188,
         "input counts changed")
    need(postcombined == post4281 - combined586,
         "post-combined residual is not the exact complement")
    need(surgery256 <= post4281, "256-surgery family leaves post4281")

    # The detached carrier audit proves every row of this exact retained band.
    carrier = {edge for edge in postcombined if edge[1] >= 645}
    need(len(carrier) == 90 and min(r for _, r in carrier) == 645 and
         max(r for _, r in carrier) == 663,
         "carrier band changed")
    need(len(combined586 & surgery256) == 9 and
         not (combined586 & carrier) and
         len(surgery256 & carrier) == 5 and
         not (combined586 & surgery256 & carrier),
         "certificate-overlap accounting changed")

    surgery256_net = surgery256 - (combined586 | carrier)
    union = combined586 | surgery256 | carrier
    final = post4281 - union
    need(len(surgery256_net) == 174 and len(union) == 850 and
         len(final) == 23373,
         "proof-graph sizes changed")
    need(max(r for _, r in final) == 644 and
         sorted(edge for edge in final if edge[1] == 644) == [
             (220, 644), (256, 644), (258, 644), (294, 644),
             (366, 644), (416, 644), (512, 644)],
         "final residual boundary changed")

    write_edges(args.emit_carrier, carrier)
    write_edges(args.emit_surgery256_net, surgery256_net)
    write_edges(args.emit_union, union)
    write_edges(args.emit_final, final)

    layers = Counter(r for _, r in carrier)
    print("THM4282_PROOF_GRAPH_CONSEQUENCE_V1")
    print(describe("COMBINED_520_367", combined586))
    print(describe("SURGERY_256", surgery256))
    print(describe("CARRIER_POSTCOMBINED_BAND", carrier))
    print("OVERLAPS COMBINED_SURGERY256", len(combined586 & surgery256),
          "COMBINED_CARRIER", len(combined586 & carrier),
          "SURGERY256_CARRIER", len(surgery256 & carrier), "TRIPLE", 0)
    print(describe("SURGERY_256_NET", surgery256_net))
    print(describe("THEOREM_UNION", union))
    print(describe("FINAL_RESIDUAL", final))
    print("CARRIER_LAYERS", " ".join(
        f"{endpoint}:{layers[endpoint]}" for endpoint in sorted(layers,
                                                                 reverse=True)))
    print("VERDICT PASS EXACT_TYPED_SET_UNION_AND_BOUNDARY")


if __name__ == "__main__":
    main()
