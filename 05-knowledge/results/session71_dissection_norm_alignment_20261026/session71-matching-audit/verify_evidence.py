#!/usr/bin/env python3
"""Archived byte-integrity verifier, not a proof checker."""
import hashlib
import json
from pathlib import Path
import tempfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(path, expected):
    if digest(path) != expected:
        raise ValueError("SHA256 mismatch: " + path.name)


def main():
    base = Path(__file__).resolve().parent
    manifest = json.loads((base / "EVIDENCE.json").read_text())
    records = []
    for item in manifest["owned_artifacts"]:
        path = base / item["path"]
        verify(path, item["sha256"])
        records.append(item["path"])
    corruptions = []
    with tempfile.TemporaryDirectory(prefix="session71-audit-corruption-") as directory:
        for name in ("check.py", "check_penrose.py", "normal.json", "penrose.normal.json"):
            source = base / name
            expected = digest(source)
            copy = Path(directory) / name
            raw = source.read_bytes()
            if not raw:
                raise ValueError("empty corruption-test source")
            copy.write_bytes(raw)
            verify(copy, expected)
            changed = bytes((raw[0] ^ 1,)) + raw[1:]
            copy.write_bytes(changed)
            try:
                verify(copy, expected)
            except ValueError:
                corruptions.append({"source": name, "one_byte_change_rejected": True})
            else:
                raise RuntimeError("one-byte corruption was accepted")
    print(json.dumps({"status": "PASS", "verified_files": records,
                      "actual_copy_corruptions": corruptions,
                      "scope": "frozen byte integrity, not mathematical proof verification"},
                     sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
