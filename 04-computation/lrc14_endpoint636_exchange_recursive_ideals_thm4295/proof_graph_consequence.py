#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import sys
from pathlib import Path


EMPTY = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def read_pairs(path: Path) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text().splitlines() if line]
    if rows != sorted(set(rows)):
        raise ValueError(f"unordered or duplicate rows: {path}")
    return rows


def raw(rows: list[tuple[int, int]]) -> bytes:
    return b"".join(f"{q},{r}\n".encode() for q, r in rows)


def fnv(rows: list[tuple[int, int]]) -> int:
    state = EMPTY
    for q, r in rows:
        for value in (q, r):
            for byte in value.to_bytes(8, "little"):
                state ^= byte
                state = state * PRIME & MASK64
    return state


def identity(label: str, rows: list[tuple[int, int]]) -> None:
    payload = raw(rows)
    print(f"{label} COUNT {len(rows)} FNV {fnv(rows):016x} "
          f"SHA256 {hashlib.sha256(payload).hexdigest()}")


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit(
            "usage: union LIVE22647 CARRIER10 H19 H294 H372 OUT_DIR")
    live_path, carrier_path, h19_path, h294_path, h372_path = \
        map(Path, sys.argv[1:6])
    out_dir = Path(sys.argv[6])
    out_dir.mkdir(parents=True, exist_ok=True)
    ledgers = {
        "LIVE": read_pairs(live_path),
        "CARRIER10": read_pairs(carrier_path),
        "H19": read_pairs(h19_path),
        "H294": read_pairs(h294_path),
        "H372": read_pairs(h372_path),
    }
    expected = {
        "LIVE": (22647, 0xDF5374D4ACA67677),
        "CARRIER10": (10, 0x9926701692E6F8D4),
        "H19": (36, 0x5C8AF37CF2F002E7),
        "H294": (21, 0xEADEFA2FAE582CA7),
        "H372": (54, 0x47AB2AF18F07FF59),
    }
    for label, rows in ledgers.items():
        if (len(rows), fnv(rows)) != expected[label]:
            raise ValueError(f"{label} identity changed")
    live = set(ledgers["LIVE"])
    nodes = {label: set(rows) for label, rows in ledgers.items()
             if label != "LIVE"}
    for label, rows in nodes.items():
        if not rows <= live:
            raise ValueError(f"{label} escapes live residual")

    print("LRC14_CARRIER10_THREE_IDEAL_PROOF_UNION_V1")
    for label, rows in ledgers.items():
        identity(label, rows)
    labels = list(nodes)
    for i, left in enumerate(labels):
        for right in labels[i + 1:]:
            overlap = sorted(nodes[left] & nodes[right])
            identity(f"OVERLAP_{left}_{right}", overlap)
            if overlap:
                print(f"OVERLAP_ROWS_{left}_{right}",
                      " ".join(f"{q},{r}" for q, r in overlap))

    union = sorted(set().union(*nodes.values()))
    final = sorted(live - set(union))
    if len(union) != 118 or len(final) != 22529:
        raise ValueError("typed proof-union count changed")
    identity("UNION", union)
    identity("FINAL", final)
    max_r = max(r for _, r in final)
    top = [pair for pair in final if pair[1] == max_r]
    print(f"FINAL_MAX {max_r} TOP_COUNT {len(top)} TOP", end="")
    for q, r in top:
        print(f" {q},{r}", end="")
    print("\nSCOPE FINITE_EXACT_TYPED_SET_UNION_SEPARATE_CERTIFICATE_NODES_NO_LRC14")
    (out_dir / "carrier10_three_ideal_union.csv").write_bytes(raw(union))
    (out_dir / "post_carrier10_three_ideal_residual.csv").write_bytes(raw(final))


if __name__ == "__main__":
    main()
