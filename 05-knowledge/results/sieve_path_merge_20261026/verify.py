#!/usr/bin/env python3
"""Verify frozen merge sources and replay the two-step decoder control."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

P = Path(__file__).resolve().parent
R = P.parents[2]


def require(value, message):
    if not value:
        raise RuntimeError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    m = json.loads((P / "MANIFEST.json").read_text())
    expected = set(m["files"]) | {"MANIFEST.json"}
    require({p.name for p in P.iterdir()} == expected, "closed file universe")
    for name, digest in m["files"].items():
        p = P / name
        require(p.is_file() and not p.is_symlink(), "regular file: " + name)
        require(sha(p.read_bytes()) == digest, "file hash: " + name)
    for source in m["reviewed_sources"]:
        raw = subprocess.check_output(
            ["git", "show", source["commit"] + ":" + source["path"]],
            cwd=R, env={**os.environ, "GIT_NO_LAZY_FETCH": "1"}, timeout=60)
        require(sha(raw) == source["sha256"], "historical source: " + source["path"])
    dp = R / m["support_source"]["path"]
    require(sha(dp.read_bytes()) == m["support_source"]["sha256"], "frozen DP dependency")
    receipts = []
    for optimized, name in [(False, "normal.json"), (True, "optimized.json")]:
        command = [sys.executable, "-I", "-B"] + (["-O"] if optimized else [])
        command += [str(P / "check_two_step.py"), "--dp-source", str(dp)]
        result = subprocess.run(command, cwd=R, capture_output=True, timeout=120)
        require(result.returncode == 0 and not result.stderr,
                "decoder replay failed: " + result.stderr.decode(errors="replace"))
        expected_receipt = (P / name).read_bytes()
        require(result.stdout == expected_receipt, "receipt mismatch: " + name)
        receipts.append(expected_receipt)
    require(receipts[0] == receipts[1], "normal/optimized mismatch")
    receipt = json.loads(receipts[0])
    require(receipt["status"] == "FINITE-EXACT PASS" and receipt["gates"] == 83912,
            "receipt status and count")
    print(json.dumps({"status": "PASS", "gates_per_mode": 83912,
                      "labelled_tournaments": 33867,
                      "historical_sources": len(m["reviewed_sources"])}, sort_keys=True))


if __name__ == "__main__":
    main()
