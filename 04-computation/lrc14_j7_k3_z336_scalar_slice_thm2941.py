#!/usr/bin/env python3
"""Exact projected ``k=3`` scalar slice at first drift 336.

The canonical scalar atlas records eight surviving six-body carriers at this
height.  This fixed-purpose replay evaluates all 3,003 literal carriers,
integrates every later label through 7,000 exactly, and bounds every omitted
label by the proved THM-1094 tail estimate ``delta(z) <= 6r/(49z)``.  The
complete 2,872-row crude ledger is pinned by a semantic digest, so the eight
rows handed to the all-label ray/status closure cannot be selected by hand.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENGINE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_next_frontier_scalar_closure_thm2941.py"
)
ATLAS_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_projected_scalar_atlas_thm2941.py"
)
ATLAS_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_projected_scalar_atlas_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z336_scalar_slice_thm2941.out"
)
EXPECTED_ENGINE_SHA256 = (
    "88c563a247d59b2d9feb552935d91a2bbc5018beeed56df74c84a37a1174894b"
)
EXPECTED_ATLAS_SHA256 = (
    "ddb2d19c02c4d70cfa74141265ceac585685932e85768baf4cb98aeb3e37935b"
)
EXPECTED_ATLAS_OUTPUT_SHA256 = (
    "ce6807de6d6b7022c97839d0bf9fc8ba3b90e7b97bc5b0d4069e88563e232be6"
)
EXPECTED_FIRST = 336
EXPECTED_CRUDE_ROWS = 2_872
EXPECTED_COMPLETE_RECORD_SHA256 = (
    "492265682865b03a4a8df4d5bcc488c6c788b4d46955a4fb504025d714f4cb1f"
)
EXPECTED_BODIES = (
    (1, 2, 6, 8, 12, 14),
    (1, 2, 6, 10, 12, 14),
    (1, 4, 6, 8, 12, 14),
    (1, 4, 8, 10, 12, 14),
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


require(file_sha256(ENGINE_PATH) == EXPECTED_ENGINE_SHA256, "scalar engine changed")
require(file_sha256(ATLAS_PATH) == EXPECTED_ATLAS_SHA256, "scalar atlas changed")
require(
    file_sha256(ATLAS_OUTPUT_PATH) == EXPECTED_ATLAS_OUTPUT_SHA256,
    "scalar atlas transcript changed",
)
spec = importlib.util.spec_from_file_location("z336_scalar_engine", ENGINE_PATH)
require(spec is not None and spec.loader is not None, "cannot load scalar engine")
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)
atlas_spec = importlib.util.spec_from_file_location("z336_scalar_atlas", ATLAS_PATH)
require(
    atlas_spec is not None and atlas_spec.loader is not None,
    "cannot load scalar atlas",
)
atlas = importlib.util.module_from_spec(atlas_spec)
atlas_spec.loader.exec_module(atlas)
require(
    dict(atlas.EXPECTED_UPPER_SPIKES)[EXPECTED_FIRST] == len(EXPECTED_BODIES),
    "canonical atlas z1=336 count changed",
)


def evaluate(args):
    return engine.evaluate(args)


def scan_fixed_first(first, processes):
    """Return the complete crude ledger and surviving scalar rows."""
    bodies = tuple(combinations(range(1, 15), 6))
    with mp.Pool(processes) as pool:
        evaluated = tuple(
            pool.imap_unordered(
                evaluate,
                ((body, first) for body in bodies),
                chunksize=8,
            )
        )
    crude_count = sum(passed for passed, _record in evaluated)
    records = tuple(sorted(record for _passed, record in evaluated if record is not None))
    survivors = tuple(record for record in records if record[-1] >= 0)
    require(len(bodies) == 3_003, "body universe changed")
    require(len(records) == crude_count, "complete crude ledger missing")
    digest = hashlib.sha256()
    for record in records:
        digest.update((engine.record_text(record) + "\n").encode())
    return crude_count, records, survivors, digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=int, default=EXPECTED_FIRST)
    parser.add_argument("--processes", type=int, default=max(1, mp.cpu_count() // 2))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(
        args.first == EXPECTED_FIRST,
        "this fixed-purpose transcript only permits --first 336; import "
        "scan_fixed_first for exploratory slices",
    )
    crude_count, _records, survivors, digest = scan_fixed_first(
        args.first, args.processes
    )
    survivor_bodies = tuple(record[0] for record in survivors)
    require(crude_count == EXPECTED_CRUDE_ROWS, "z1=336 crude count changed")
    require(survivor_bodies == EXPECTED_BODIES, "z1=336 body ledger changed")
    require(digest == EXPECTED_COMPLETE_RECORD_SHA256, "z1=336 digest changed")

    positive = evaluate((engine.EXPECTED_POSITIVE[0], 380))
    require(positive == (True, engine.EXPECTED_POSITIVE), "z1=380 control changed")
    negative = evaluate((engine.EXPECTED_POSITIVE[0], 379))
    require(negative[0] and negative[1][-1] < 0, "z1=379 control changed")

    lines = [
        "LRC14 projected k=3 fixed-first scalar slice",
        f"scalar_engine_sha256={file_sha256(ENGINE_PATH)}",
        f"scalar_atlas_sha256={file_sha256(ATLAS_PATH)}",
        f"scalar_atlas_output_sha256={file_sha256(ATLAS_OUTPUT_PATH)}",
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
        *(engine.record_text(record) for record in survivors),
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
