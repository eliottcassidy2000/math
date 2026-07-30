#!/usr/bin/env python3
"""Exact projected ``k=3`` scalar slice at first drift 328.

This fixed-purpose replay imports the pinned global scanner from the adjacent
``z1=336`` package and evaluates all 3,003 literal six-body carriers.  Exactly
nine bodies survive the scalar necessary condition.  The complete 2,944-row
crude ledger, not merely those survivors, is frozen by a semantic digest.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z336_scalar_slice_thm2941.py"
)
BASE_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z336_scalar_slice_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z328_scalar_slice_thm2941.out"
)
EXPECTED_BASE_SHA256 = (
    "8ca30b3f6e8336e87b558c5e24007153aa0f0331e3b00dfe66c93cd4c93dd789"
)
EXPECTED_BASE_OUTPUT_SHA256 = (
    "e3de2a3d6c65c0e8eb62d400ff47deafe48c8aa6a2db45788eab45e8f50892b9"
)
EXPECTED_FIRST = 328
EXPECTED_CRUDE_ROWS = 2_944
EXPECTED_COMPLETE_RECORD_SHA256 = (
    "a276e07fe6c2ef9fc3aeced14aa1e877761d7b0493724afc323eca9476bdf99e"
)
EXPECTED_BODIES = (
    (1, 2, 6, 8, 12, 14),
    (1, 2, 8, 10, 12, 14),
    (1, 4, 6, 8, 12, 14),
    (1, 4, 8, 10, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (2, 4, 6, 8, 12, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 6, 8, 10, 12, 14),
    (2, 8, 10, 11, 12, 14),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


require(file_sha256(BASE_PATH) == EXPECTED_BASE_SHA256, "global scanner changed")
require(
    file_sha256(BASE_OUTPUT_PATH) == EXPECTED_BASE_OUTPUT_SHA256,
    "z1=336 scanner transcript changed",
)
base = load_module("z328_global_scalar_scanner", BASE_PATH)
require(
    base.atlas.FIRST_MIN <= EXPECTED_FIRST <= base.atlas.FIRST_MAX,
    "z1=328 lies outside the canonical atlas",
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=int, default=EXPECTED_FIRST)
    parser.add_argument("--processes", type=int, default=max(1, mp.cpu_count() // 2))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(
        args.first == EXPECTED_FIRST,
        "this fixed-purpose transcript only permits --first 328",
    )
    crude_count, _records, survivors, digest = base.scan_fixed_first(
        args.first, args.processes
    )
    require(crude_count == EXPECTED_CRUDE_ROWS, "z1=328 crude count changed")
    require(
        tuple(record[0] for record in survivors) == EXPECTED_BODIES,
        "z1=328 body ledger changed",
    )
    require(digest == EXPECTED_COMPLETE_RECORD_SHA256, "z1=328 digest changed")

    positive = base.evaluate((base.engine.EXPECTED_POSITIVE[0], 380))
    require(
        positive == (True, base.engine.EXPECTED_POSITIVE),
        "z1=380 control changed",
    )
    negative = base.evaluate((base.engine.EXPECTED_POSITIVE[0], 379))
    require(negative[0] and negative[1][-1] < 0, "z1=379 control changed")

    lines = [
        "LRC14 projected k=3 fixed-first scalar slice",
        f"global_scanner_sha256={file_sha256(BASE_PATH)}",
        f"global_scanner_output_sha256={file_sha256(BASE_OUTPUT_PATH)}",
        f"scalar_engine_sha256={file_sha256(base.ENGINE_PATH)}",
        f"scalar_atlas_sha256={file_sha256(base.ATLAS_PATH)}",
        f"scalar_atlas_output_sha256={file_sha256(base.ATLAS_OUTPUT_PATH)}",
        (
            f"universe=six_body_roots=3003;first_drift={args.first};"
            "three_distinct_later_drifts;nonzero_mod_L;exact_horizon=7000"
        ),
        (
            "omitted_tail=delta(z)<=6r/(49z);"
            "projected_high_wall=max(15,floor(13L/132)+1)"
        ),
        f"crude_all_tail_rows={crude_count}",
        f"scalar_surviving_bodies={len(survivors)}",
        f"complete_record_sha256={digest}",
        "survivor_rows_begin",
        *(base.engine.record_text(record) for record in survivors),
        "survivor_rows_end",
        "positive_control=z1=380:PASS;negative_control=z1=379:PASS",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
