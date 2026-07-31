#!/usr/bin/env python3
"""Verify byte parity of all ordinary and optimized ranked-H1 shards."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCOUT = HERE / "scout.py"
SCOUT_SHA256 = (
    "6abc06972cda64bb4f53db26cca62bbea5362ca178bdbf5bf7398e8b0f28317a"
)
MERGE = HERE / "merge_all_32_shards.py"
MERGE_SHA256 = (
    "0677f03ec4699d2b652d7e118e96553d9dcefd2c9ee0e2109835bf429870a3c9"
)
EXPECTED_FILE_DIGEST = (
    "d686d1afd01b494ca6aa9a4733816c774f7e8e4ed398081f556c44e4af7869a0"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_merge():
    require(sha256(SCOUT) == SCOUT_SHA256, "scout source changed")
    require(sha256(MERGE) == MERGE_SHA256, "merge verifier changed")
    spec = importlib.util.spec_from_file_location("h1_merge_replay", MERGE)
    require(spec is not None and spec.loader is not None, "cannot load merge")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    merge = load_merge()
    require(len(merge.SHARD_HASHES) == 32, "shard hash battery changed")
    file_digest = hashlib.sha256(b"LRC14/j6/H1-shards/v1\n")
    corpus_digest = hashlib.sha256(b"LRC14/j6/H1-ordinary-replay-corpus/v1\n")
    worker_modes = []
    for index, expected_hash in enumerate(merge.SHARD_HASHES):
        optimized = HERE / f"all_s{index:02d}_of32.out"
        ordinary = HERE / f"all_s{index:02d}_of32.ordinary.out"
        optimized_bytes = optimized.read_bytes()
        ordinary_bytes = ordinary.read_bytes()
        require(
            ordinary_bytes == optimized_bytes,
            f"ordinary/optimized byte mismatch at shard {index}",
        )
        actual_hash = hashlib.sha256(ordinary_bytes).hexdigest()
        require(actual_hash == expected_hash, f"shard {index} hash changed")
        require(
            ordinary_bytes.endswith(b"all_exact_controls=PASS\n"),
            f"shard {index} did not finish",
        )
        parameter_line = ordinary_bytes.splitlines()[1].decode()
        expected_workers = 4 if index == 0 else 1
        require(
            parameter_line.startswith(f"parameters=workers:{expected_workers},"),
            f"shard {index} worker mode changed",
        )
        worker_modes.append(expected_workers)
        file_digest.update(f"{index:02d};{actual_hash}\n".encode())
        corpus_digest.update(f"{index:02d};{len(ordinary_bytes)}\n".encode())
        corpus_digest.update(ordinary_bytes)
    require(
        file_digest.hexdigest() == EXPECTED_FILE_DIGEST,
        "aggregate shard-file digest changed",
    )
    require(
        tuple(worker_modes) == (4,) + (1,) * 31,
        "ordinary replay worker-mode schedule changed",
    )
    print("LRC14 ranked-H1 all-32 ordinary replay parity")
    print("shards=32;byte_matches=32;failures=0")
    print("worker_schedule=(shard0:4,shards1-31:1)")
    print(f"shard_file_digest={file_digest.hexdigest()}")
    print(f"ordinary_replay_corpus_digest={corpus_digest.hexdigest()}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
