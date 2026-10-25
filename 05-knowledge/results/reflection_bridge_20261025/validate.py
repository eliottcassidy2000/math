"""Replay all local exact controls in normal and optimized Python modes."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CHECKS = {
    "check_divisor_ideal_graph.py": "divisor_ideal_graph_receipt.json",
    "graph_bridge_check.py": "graph_bridge_receipt.json",
    "audit_graph_bridge.py": "audit_graph_bridge_receipt.json",
    "seven_check.py": "seven_check.json",
    "geometry_check.py": "geometry_receipt.json",
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    results = []
    for name, receipt in CHECKS.items():
        runs = []
        for mode in ([], ["-O"]):
            proc = subprocess.run([sys.executable, *mode, name], cwd=ROOT,
                                  capture_output=True)
            if proc.returncode:
                raise RuntimeError(f"{name} {mode}: {proc.stderr.decode(errors='replace')}")
            runs.append(((ROOT / receipt).read_bytes(), proc.stdout))
        if runs[0] != runs[1]:
            raise RuntimeError(f"normal/optimized mismatch: {name}")
        results.append({"checker": name, "receipt": receipt,
                        "checker_sha256": digest((ROOT / name).read_bytes()),
                        "receipt_sha256": digest(runs[0][0]),
                        "normal_and_optimized_identical": True})
    output = json.dumps(results, indent=2) + "\n"
    (ROOT / "validation_receipt.json").write_text(output, encoding="utf-8")
    print(f"PASS: {len(CHECKS)} checkers, normal and optimized outputs identical")


if __name__ == "__main__":
    main()
