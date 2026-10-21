#!/usr/bin/env python3
"""Verify frozen local files and replay exact controls in normal/-O modes."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
manifest = json.loads((HERE / "manifest.json").read_text())
for name, expected in manifest["files"].items():
    actual = hashlib.sha256((HERE / name).read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError(f"File hash mismatch: {name}: {actual} != {expected}")

receipt = (HERE / "receipt.json").read_bytes()
if receipt != (HERE / "receipt.optimized.json").read_bytes():
    raise RuntimeError("Frozen normal and optimized receipts differ")
for flags in ([], ["-O"]):
    result = subprocess.run([sys.executable, *flags, str(HERE / "check.py")],
                            capture_output=True, check=True)
    if result.stdout != receipt:
        raise RuntimeError(f"Replay differs from receipt in mode {flags}")
    if result.stderr:
        raise RuntimeError(f"Unexpected replay stderr in mode {flags}: {result.stderr!r}")
print(json.dumps({"status": "PASS", "files_hashed": len(manifest["files"]),
                  "normal_optimized_replay": "byte-identical",
                  "receipt_sha256": hashlib.sha256(receipt).hexdigest()},
                 indent=2, sort_keys=True))
