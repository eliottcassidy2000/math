#!/usr/bin/env python3
"""Freeze a closed SHA256 manifest, excluding executables and the manifest."""

from __future__ import annotations

import hashlib
from pathlib import Path


packet = Path(__file__).resolve().parent
files = sorted(path for path in packet.rglob("*") if path.is_file() and
               path.name != "SHA256SUMS" and path.suffix.lower() != ".exe" and
               "__pycache__" not in path.parts)
rows = []
for path in files:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    rows.append(f"{digest}  {path.relative_to(packet).as_posix()}\n")
(packet / "SHA256SUMS").write_bytes("".join(rows).encode("ascii"))
print(f"MANIFEST_FILES {len(files)}")
print(hashlib.sha256((packet / "SHA256SUMS").read_bytes()).hexdigest())
