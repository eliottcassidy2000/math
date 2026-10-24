#!/usr/bin/env python3
"""Copy only manifest-listed bytes; reproduce normal and -O reports in isolation."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "MANIFEST.json").read_text())
    seen = set()
    with tempfile.TemporaryDirectory(prefix="observer-moments-manifest-only-") as temporary:
        closed = Path(temporary)
        for row in manifest["files"]:
            relative = Path(row["path"])
            require(not relative.is_absolute() and ".." not in relative.parts,
                    "manifest path must remain inside the packet")
            require(row["path"] not in seen, "duplicate manifest path")
            seen.add(row["path"])
            source = root / relative
            data = source.read_bytes()
            require(len(data) == row["bytes"], "byte-size mismatch: " + row["path"])
            require(hashlib.sha256(data).hexdigest() == row["sha256"],
                    "SHA256 mismatch: " + row["path"])
            destination = closed / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
        require({"check.py", "normal.json", "optimized.json"} <= seen,
                "missing closed reproduction inputs")
        reports = []
        for flags, expected, generated in [([], "normal.json", "fresh-normal.json"),
                                            (["-O"], "optimized.json", "fresh-optimized.json")]:
            command = [sys.executable, "-I", "-B", *flags, "check.py", "--output", generated]
            result = subprocess.run(command, cwd=closed, text=True, capture_output=True)
            require(result.returncode == 0, "closed reproduction failed: " + result.stderr)
            data = (closed / generated).read_bytes()
            require(data == (closed / expected).read_bytes(), "reproduction differs: " + expected)
            reports.append(data)
        require(reports[0] == reports[1], "normal and optimized reports differ")
        checks = json.loads(reports[0])["checks"]
        print(json.dumps({"status": "PASS", "copied_manifest_files": len(seen),
                          "normal_and_optimized_identical": True,
                          "mathematical_checks_per_run": checks,
                          "report_sha256": hashlib.sha256(reports[0]).hexdigest()},
                         indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
