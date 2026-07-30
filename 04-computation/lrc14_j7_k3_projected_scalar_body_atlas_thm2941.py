#!/usr/bin/env python3
"""Lossless body ledger for the projected ``k=3`` scalar atlas.

The canonical one-pass atlas already pins all 6,060 surviving row records by
one global semantic digest, but its transcript expands only the ``z1=364``
slice.  This sidecar reruns that same exact 3,003-body computation once and
emits every row.  It therefore supplies lossless body ledgers for all 179
first drifts without repeating a separate global scan at every lower spike.

Every row still uses exact singleton coverage through 7,000, the proved
THM-1094 omitted-label bound, and the projected high-label obligation.  The
canonical row digest is rechecked before any output is written; a second
digest freezes just the ordered ``(first drift, body)`` incidence ledger.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import multiprocessing as mp
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ATLAS_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_projected_scalar_atlas_thm2941.py"
)
ATLAS_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_projected_scalar_atlas_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
)
EXPECTED_ATLAS_SHA256 = (
    "ddb2d19c02c4d70cfa74141265ceac585685932e85768baf4cb98aeb3e37935b"
)
EXPECTED_ATLAS_OUTPUT_SHA256 = (
    "ce6807de6d6b7022c97839d0bf9fc8ba3b90e7b97bc5b0d4069e88563e232be6"
)
EXPECTED_ROWS = 6_060
EXPECTED_ROW_SEMANTIC_SHA256 = (
    "46535b2b075e15244ec87101f44d74a8034b093244fd8559fe11b5b70425fda8"
)
EXPECTED_BODY_LEDGER_SHA256 = (
    "d75815ba6ea06b826cd7cb4363678b013c984d0f1891882564ceb5b034e94fed"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


require(file_sha256(ATLAS_PATH) == EXPECTED_ATLAS_SHA256, "scalar atlas changed")
require(
    file_sha256(ATLAS_OUTPUT_PATH) == EXPECTED_ATLAS_OUTPUT_SHA256,
    "scalar atlas transcript changed",
)
sys.path.insert(0, str(ROOT / "04-computation"))
atlas = importlib.import_module("lrc14_j7_k3_projected_scalar_atlas_thm2941")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=max(1, mp.cpu_count() // 2))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, "process count must be positive")

    rows = atlas.atlas(args.processes)
    require(len(rows) == EXPECTED_ROWS, ("atlas row count changed", len(rows)))
    counts = Counter(row[5] for row in rows)
    require(
        min(counts) >= atlas.FIRST_MIN and max(counts) <= atlas.FIRST_MAX,
        "first-drift support changed",
    )

    row_digest = hashlib.sha256()
    body_digest = hashlib.sha256()
    row_texts = []
    previous_key = None
    for row in rows:
        text = atlas.row_text(row)
        key = (row[5], row[0])
        require(previous_key is None or previous_key < key, "row order changed")
        previous_key = key
        row_texts.append(text)
        row_digest.update((text + "\n").encode())
        body_digest.update(f"{row[5]}|{row[0]}\n".encode())
    row_sha256 = row_digest.hexdigest()
    body_sha256 = body_digest.hexdigest()
    require(
        row_sha256 == EXPECTED_ROW_SEMANTIC_SHA256,
        "canonical row semantic digest changed",
    )
    if EXPECTED_BODY_LEDGER_SHA256 is not None:
        require(
            body_sha256 == EXPECTED_BODY_LEDGER_SHA256,
            "body incidence ledger changed",
        )

    lines = [
        "LRC14 projected k=3 lossless scalar body atlas",
        f"scalar_atlas_source_sha256={file_sha256(ATLAS_PATH)}",
        f"scalar_atlas_output_sha256={file_sha256(ATLAS_OUTPUT_PATH)}",
        (
            "universe=six_body_roots=3003;first_drift=200..378;"
            "three_distinct_later_drifts;nonzero_mod_L;exact_horizon=7000"
        ),
        (
            "omitted_tail=delta(z)<=6r/(49z);"
            "projected_high_wall=max(15,floor(13L/132)+1)"
        ),
        f"total_projected_scalar_rows={len(rows)}",
        f"counts_by_first={tuple(sorted(counts.items()))}",
        f"body_ledger_sha256={body_sha256}",
        "rows_begin",
        *(f"row={text}" for text in row_texts),
        "rows_end",
        f"row_semantic_sha256={row_sha256}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
