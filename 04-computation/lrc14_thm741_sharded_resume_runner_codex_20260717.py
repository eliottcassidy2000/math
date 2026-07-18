#!/usr/bin/env python3
"""Portable, collision-free sharded runner for THM-741's 2002-body engine.

The validated mathematical kernel is hash-pinned and left unchanged.  This
file owns only orchestration: deterministic sharding, strict JSONL reads,
per-shard append targets, cross-shard deduplication, and final canonical merge.
It replaces the 2026-07-17 inline sharding edit that accidentally moved the
full-run block out of `main` and still wrote every worker to one state file.

Examples:

  python3 ... --self-test
  python3 ... --state-dir /tmp/thm741 --audit-only --num-shards 6
  python3 ... --state-dir /tmp/thm741 --worker-id 0 --num-shards 6 --workers 4
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import pickle
import sys
import tempfile
import time
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
BODIES = tuple(combinations(range(1, 15), 9))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(file_sha256(CORE_PATH) == CORE_SHA256, "THM-741 kernel hash changed")
    spec = importlib.util.spec_from_file_location("thm741_sharded_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 kernel")
    module = importlib.util.module_from_spec(spec)
    # Multiprocessing pickles CORE.body_work by module/name.  Register the
    # hash-pinned dynamic module in both the parent and spawned children.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


CORE = load_core()


def shard_of(index: int, num_shards: int) -> int:
    """Stable lexicographic round-robin partition; no process hash is used."""
    require(num_shards >= 1, "num_shards must be positive")
    return index % num_shards


def state_paths(state_dir: Path) -> tuple[Path, ...]:
    canonical = state_dir / "thm741_results.jsonl"
    shards = tuple(sorted(state_dir.glob("thm741_results.shard*.jsonl")))
    return ((canonical,) if canonical.exists() else ()) + shards


def read_rows(paths: tuple[Path, ...]) -> dict[tuple[int, ...], dict[str, object]]:
    """Read every row strictly; conflicting duplicates are a hard error."""
    rows: dict[tuple[int, ...], dict[str, object]] = {}
    expected = set(BODIES)
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                require(line.endswith("\n"), f"partial final line {path}:{line_number}")
                try:
                    row = json.loads(line)
                    body = tuple(int(value) for value in row["E"])
                except (json.JSONDecodeError, KeyError, TypeError, ValueError) as error:
                    raise RuntimeError(f"malformed row {path}:{line_number}: {error}") from error
                require(body in expected, f"unknown body {body} in {path}:{line_number}")
                if body in rows:
                    require(rows[body] == row,
                            f"conflicting duplicate body {body} in {path}:{line_number}")
                else:
                    rows[body] = row
    return rows


def pending_by_shard(done: set[tuple[int, ...]], num_shards: int):
    pending = [[] for _ in range(num_shards)]
    for index, body in enumerate(BODIES):
        if body not in done:
            pending[shard_of(index, num_shards)].append(body)
    return tuple(tuple(bucket) for bucket in pending)


def atomic_canonical_merge(state_dir: Path, rows: dict[tuple[int, ...], dict[str, object]]) -> Path:
    require(len(rows) == len(BODIES), "canonical merge requires all 2002 bodies")
    destination = state_dir / "thm741_results.jsonl"
    temporary = state_dir / "thm741_results.jsonl.tmp"
    with temporary.open("w", encoding="utf-8") as handle:
        for body in BODIES:
            handle.write(json.dumps(rows[body], sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(destination)
    return destination


def self_test() -> None:
    require(len(BODIES) == 2002 and len(set(BODIES)) == 2002, "body census")
    require(len(pickle.dumps(CORE.body_work)) > 0, "kernel worker is not picklable")
    balance = {}
    for num_shards in range(1, 17):
        buckets = pending_by_shard(set(), num_shards)
        flattened = [body for bucket in buckets for body in bucket]
        require(len(flattened) == 2002 and set(flattened) == set(BODIES),
                f"partition coverage n={num_shards}")
        require(sum(len(bucket) for bucket in buckets) == len(set(flattened)),
                f"partition overlap n={num_shards}")
        sizes = tuple(map(len, buckets))
        require(max(sizes) - min(sizes) <= 1, f"partition imbalance n={num_shards}")
        balance[num_shards] = sizes

    # Strict reader, identical duplicate, conflicting duplicate, and partial
    # final-line behavior are exercised without touching the repository.
    with tempfile.TemporaryDirectory(prefix="thm741-shard-selftest-") as directory:
        root = Path(directory)
        row = {"E": list(BODIES[0]), "marker": 1}
        first = root / "thm741_results.shard0.jsonl"
        second = root / "thm741_results.shard1.jsonl"
        first.write_text(json.dumps(row) + "\n", encoding="utf-8")
        second.write_text(json.dumps(row) + "\n", encoding="utf-8")
        require(len(read_rows((first, second))) == 1, "identical duplicate merge")
        second.write_text(json.dumps({"E": list(BODIES[0]), "marker": 2}) + "\n",
                          encoding="utf-8")
        try:
            read_rows((first, second))
        except RuntimeError as error:
            require("conflicting duplicate" in str(error), "wrong conflict error")
        else:
            raise RuntimeError("conflicting duplicate was accepted")
        second.write_text(json.dumps(row), encoding="utf-8")
        try:
            read_rows((second,))
        except RuntimeError as error:
            require("partial final line" in str(error), "wrong partial-line error")
        else:
            raise RuntimeError("partial final line was accepted")

    print("THM-741 SHARDED RESUME RUNNER SELF-TEST")
    print(f"kernel_sha256={CORE_SHA256}")
    print("bodies=2002; deterministic lex-index partition; shard counts differ by at most one")
    print("partition audit n=1..16: disjoint, exhaustive, balanced; worker picklable")
    print("state audit: identical duplicate coalesced; conflict rejected; partial line rejected")
    print("six-shard sizes=" + ",".join(map(str, balance[6])))
    print(f"source_sha256={file_sha256(Path(__file__).resolve())}")
    print("ALL SHARD ORCHESTRATION CHECKS PASSED")


def audit(state_dir: Path, num_shards: int) -> dict[tuple[int, ...], dict[str, object]]:
    paths = state_paths(state_dir)
    rows = read_rows(paths)
    pending = pending_by_shard(set(rows), num_shards)
    print("THM-741 SHARDED STATE AUDIT")
    print(f"kernel_sha256={CORE_SHA256}")
    print(f"state_dir={state_dir}; files={len(paths)}; completed={len(rows)}/2002")
    print("pending_by_shard=" + ",".join(map(lambda bucket: str(len(bucket)), pending)))
    print("reader policy=strict JSONL; no malformed/partial rows; duplicate rows must agree exactly")
    return rows


def run_worker(state_dir: Path, worker_id: int, num_shards: int, workers: int) -> None:
    require(0 <= worker_id < num_shards, "worker_id outside shard range")
    state_dir.mkdir(parents=True, exist_ok=True)
    rows = audit(state_dir, num_shards)
    pending = pending_by_shard(set(rows), num_shards)[worker_id]
    output = state_dir / f"thm741_results.shard{worker_id}.jsonl"
    print(
        f"worker={worker_id}/{num_shards}; pending={len(pending)}; "
        f"processes={workers}; output={output}",
        flush=True,
    )
    started = time.time()
    with mp.Pool(workers) as pool, output.open("a", encoding="utf-8") as handle:
        for index, row in enumerate(pool.imap_unordered(CORE.body_work, pending, chunksize=1), start=1):
            handle.write(json.dumps(row, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
            if index % 20 == 0 or index == len(pending):
                print(
                    f"worker={worker_id} complete={index}/{len(pending)} "
                    f"elapsed_minutes={(time.time()-started)/60:.1f}",
                    flush=True,
                )

    merged = read_rows(state_paths(state_dir))
    print(f"post-worker globally completed={len(merged)}/2002", flush=True)
    if len(merged) == len(BODIES):
        destination = atomic_canonical_merge(state_dir, merged)
        print(f"canonical merge complete: {destination}", flush=True)
    else:
        print("global run incomplete; canonical theorem report deliberately not emitted", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--state-dir", type=Path)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--worker-id", type=int)
    parser.add_argument("--num-shards", type=int, default=1)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return
    require(args.state_dir is not None, "--state-dir is required outside --self-test")
    require(args.num_shards >= 1 and args.workers >= 1, "positive shard/process counts required")
    if args.audit_only:
        audit(args.state_dir, args.num_shards)
        return
    require(args.worker_id is not None, "--worker-id is required to run a shard")
    run_worker(args.state_dir, args.worker_id, args.num_shards, args.workers)


if __name__ == "__main__":
    main()
