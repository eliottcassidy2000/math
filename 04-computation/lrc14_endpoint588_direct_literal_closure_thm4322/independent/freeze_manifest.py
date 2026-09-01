#!/usr/bin/env python3
"""Freeze SHA-256 closure for every non-executable packet file."""

from hashlib import sha256
from pathlib import Path

packet = Path(__file__).resolve().parent
paths = sorted(
    path for path in packet.rglob("*")
    if path.is_file() and path.name != "SHA256SUMS" and path.suffix != ".exe"
)
payload = "".join(
    f"{sha256(path.read_bytes()).hexdigest()}  {path.relative_to(packet).as_posix()}\n"
    for path in paths
)
(packet / "SHA256SUMS").write_text(payload, encoding="utf-8", newline="\n")
print(f"MANIFEST_FILES {len(paths)}")
print(f"MANIFEST_SHA256 {sha256(payload.encode()).hexdigest()}")
