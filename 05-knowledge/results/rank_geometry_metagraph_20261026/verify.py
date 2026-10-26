#!/usr/bin/env python3
"""Replay each lane normally and under -O; compare deterministic receipts."""
import json
from pathlib import Path
import subprocess
import sys


root = Path(__file__).resolve().parent
results = {}
for lane in ("information", "enumeration", "unit_distance", "rank"):
    path = root / lane / "check.py"
    print("Checking " + lane, flush=True)
    normal = subprocess.check_output([sys.executable, str(path)], cwd=path.parent)
    optimized = subprocess.check_output([sys.executable, "-O", str(path)], cwd=path.parent)
    a, b = json.loads(normal), json.loads(optimized)
    if a != b:
        raise RuntimeError(lane + ": normal/-O receipts differ")
    results[lane] = a
    print(lane + ": PASS, normal == -O", flush=True)
print(json.dumps({"paired_lanes": list(results), "status": "PASS"}, sort_keys=True))
