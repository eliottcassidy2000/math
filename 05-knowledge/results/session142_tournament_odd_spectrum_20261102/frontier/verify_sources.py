#!/usr/bin/env python3
"""Offline check of preserved source bytes against Git blob IDs and SHA256."""
from hashlib import sha1, sha256
import json
from pathlib import Path


def main():
    base = Path(__file__).resolve().parent
    records = json.loads((base / "repository_sources.json").read_text())
    snapshots = set()
    for r in records:
        target = (base / r["local_snapshot"]).resolve()
        if not target.is_relative_to(base / "sources"):
            raise RuntimeError("snapshot outside local sources directory")
        data = target.read_bytes()
        actual_blob = sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
        if (len(data), sha256(data).hexdigest(), actual_blob) != (
            r["bytes"], r["sha256"], r["git_blob"]
        ):
            raise RuntimeError(f"source identity mismatch: {r['local_snapshot']}")
        snapshots.add(r["local_snapshot"])
    print(json.dumps({"source_pins": len(records), "distinct_snapshots": len(snapshots),
                      "raw_bytes_sha256_and_git_blob": "PASS"}, sort_keys=True))


if __name__ == "__main__":
    main()
