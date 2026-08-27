#!/usr/bin/env python3
"""Exact cross-filter of a dormant carrier through the (520,663) response atlas.

This is a reconstruction/audit sidecar only.  It does not assert theorem-level
provenance or body coverage for the dormant carrier.
"""

from __future__ import annotations

import csv
import hashlib
import pathlib
import sys


FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
EXPECTED_CARRIER_SHA = "a5dac3c7e5a2715e4c9ef8bb1b54bc98792e904f7b0b5ef55e4dd4313ebc87f6"
EXPECTED_ACTIVE_SHA = "dbdd1dab8c1e0d4de03866d20155d2a985929bd9e279b3fb15eb643c351a851e"
EXPECTED_LOST_SHA = "9d4a5bff8436e8d922ac4ae685e1340e43a9a143399e2590ef8c1051586b026d"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def fnv_add(state: int, value: int) -> int:
    require(0 <= value < 1 << 64, "FNV input escaped u64")
    for byte in value.to_bytes(8, "little"):
        state ^= byte
        state = (state * FNV_PRIME) & ((1 << 64) - 1)
    return state


def read_carrier(path: pathlib.Path) -> tuple[list[int], str, int]:
    digest = sha256(path)
    require(digest == EXPECTED_CARRIER_SHA, "dormant carrier SHA changed")
    masks: list[int] = []
    state = FNV_OFFSET
    with path.open(encoding="ascii", newline="") as stream:
        for line_number, raw in enumerate(stream, 1):
            token = raw.strip()
            require(token and len(token) == 8, f"bad carrier token at {line_number}")
            mask = int(token, 16)
            require(mask < 1 << 30 and mask.bit_count() == 8,
                    f"carrier token off rank eight at {line_number}")
            masks.append(mask)
            state = fnv_add(state, mask)
    require(len(masks) == 8951 and len(set(masks)) == len(masks),
            "dormant carrier count/distinctness changed")
    return masks, digest, state


def read_obligations(path: pathlib.Path) -> tuple[list[int], str]:
    digest = sha256(path)
    require(digest == EXPECTED_LOST_SHA, "lost-obligation SHA changed")
    bodies: list[int] = []
    with path.open(encoding="ascii", newline="") as stream:
        rows = csv.DictReader(stream)
        require(rows.fieldnames == ["ordinal", "body_hex", "deleted_response_hex"],
                "lost-obligation header changed")
        for row in rows:
            require(int(row["ordinal"]) == len(bodies), "obligation ordinal changed")
            body = int(row["body_hex"], 16)
            require(body.bit_count() == 9 and body < 1 << 30,
                    "obligation body off rank nine")
            bodies.append(body)
    require(len(bodies) == 53 and bodies == sorted(bodies),
            "obligation count/order changed")
    return bodies, digest


def main() -> int:
    require(len(sys.argv) == 6,
            "usage: cross-filter CARRIER ACTIVE_TSV LOST_CSV FILTER_TSV HITS_TSV")
    carrier_path, active_path, lost_path, filter_path, hits_path = map(
        pathlib.Path, sys.argv[1:])
    carrier, carrier_sha, carrier_fnv = read_carrier(carrier_path)
    carrier_set = set(carrier)
    bodies, lost_sha = read_obligations(lost_path)
    active_sha = sha256(active_path)
    require(active_sha == EXPECTED_ACTIVE_SHA, "active-response TSV SHA changed")

    rows_seen = 0
    previous_rank = -1
    intersection: list[tuple[int, int, int, int]] = []
    all_ledger = FNV_OFFSET
    nonempty_ledger = FNV_OFFSET
    union = 0
    incidences = 0
    nonempty = 0
    hit_counts = [0] * len(bodies)
    with active_path.open(encoding="ascii", newline="") as stream:
        header = stream.readline().rstrip("\n")
        require(header == "colex_rank\tmask_hex\tclass\tpattern_hex",
                "active-response header changed")
        for raw in stream:
            fields = raw.rstrip("\n").split("\t")
            require(len(fields) == 4, "malformed active-response row")
            rank = int(fields[0])
            mask = int(fields[1], 16)
            class_id = int(fields[2])
            pattern = int(fields[3], 16)
            require(rank > previous_rank and mask.bit_count() == 8 and
                    mask < 1 << 30 and pattern < 1 << len(bodies),
                    "active-response rank/mask/pattern changed")
            previous_rank = rank
            rows_seen += 1
            if mask not in carrier_set:
                continue
            intersection.append((rank, mask, class_id, pattern))
            all_ledger = fnv_add(fnv_add(all_ledger, mask), pattern)
            if pattern == 0:
                continue
            nonempty += 1
            nonempty_ledger = fnv_add(fnv_add(nonempty_ledger, mask), pattern)
            union |= pattern
            incidences += pattern.bit_count()
            bits = pattern
            while bits:
                bit = (bits & -bits).bit_length() - 1
                hit_counts[bit] += 1
                bits &= bits - 1
    require(rows_seen == 2_879_147, "active-response row count changed")
    require(len(intersection) == 7_220 and nonempty == 3_007 and
            incidences == 6_580 and union == (1 << len(bodies)) - 1,
            "dormant carrier cross-filter totals changed")
    require(min(hit_counts) == 18 and max(hit_counts) == 290 and
            hit_counts.index(min(hit_counts)) == 38,
            "dormant carrier obligation-hit range changed")

    with filter_path.open("w", encoding="ascii", newline="") as stream:
        stream.write("colex_rank\tmask_hex\tclass\tpattern_hex\n")
        for rank, mask, class_id, pattern in intersection:
            stream.write(f"{rank}\t{mask:08x}\t{class_id}\t{pattern:014x}\n")
    with hits_path.open("w", encoding="ascii", newline="") as stream:
        stream.write("obligation\tbody_hex\thits\n")
        for ordinal, (body, hits) in enumerate(zip(bodies, hit_counts)):
            stream.write(f"{ordinal}\t{body:08x}\t{hits}\n")

    print("THM4281_DORMANT_CARRIER_RESPONSE_CROSS_FILTER_V1")
    print(f"CARRIER 8951 SHA256 {carrier_sha} FNV {carrier_fnv:016x}")
    print(f"ACTIVE_RESPONSE_ROWS {rows_seen} SHA256 {active_sha}")
    print(f"OBLIGATIONS {len(bodies)} SHA256 {lost_sha}")
    print(f"ACTIVE_INTERSECTION {len(intersection)} FNV {all_ledger:016x}")
    print(f"NONEMPTY_RESPONSE {nonempty} FNV {nonempty_ledger:016x}")
    print(f"RESPONSE_UNION {union:014x} COVERED {union.bit_count()} OF {len(bodies)}")
    print(f"INCIDENCES {incidences} HIT_MIN {min(hit_counts)} MIN_OBLIGATION "
          f"{hit_counts.index(min(hit_counts))} HIT_MAX {max(hit_counts)}")
    print("SCOPE RECONSTRUCTION_AUDIT_ONLY NO_CARRIER_COVERAGE_OR_PROVENANCE_CLAIM")
    print("VERDICT PASS DORMANT_CARRIER_CROSS_FILTER")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"CROSS_FILTER_ERROR {error}", file=sys.stderr)
        raise SystemExit(1)
