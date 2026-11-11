#!/usr/bin/env python3
"""Verify this published packet's bytes, not the truth of its mathematics."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parent
manifest = json.loads((root / "MANIFEST.json").read_text())
errors = []
for record in manifest["files"]:
    path = root / record["path"]
    if not path.is_file():
        errors.append(f"Missing: {record['path']}")
        continue
    data = path.read_bytes()
    if len(data) != record["bytes"] or hashlib.sha256(data).hexdigest() != record["sha256"]:
        errors.append(f"Changed: {record['path']}")
if errors:
    raise SystemExit("\n".join(errors))
print(f"PASS: {len(manifest['files'])} published files match the SHA-256 manifest.")
