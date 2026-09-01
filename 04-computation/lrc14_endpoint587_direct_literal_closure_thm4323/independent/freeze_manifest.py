#!/usr/bin/env python3
"""Freeze all non-executable packet artifacts except verifier transcripts."""

from __future__ import annotations

import hashlib
from pathlib import Path

packet = Path(__file__).resolve().parent
excluded = {
    "SHA256SUMS",
    "verify_endpoint587_cleanroom_packet.out",
    "verify_endpoint587_cleanroom_packet_opt.out",
}
files = sorted(
    path for path in packet.rglob("*")
    if path.is_file() and path.name not in excluded and path.suffix.lower() != ".exe"
)
lines = [
    f"{hashlib.sha256(path.read_bytes()).hexdigest()}  "
    f"{path.relative_to(packet).as_posix()}\n"
    for path in files
]
(packet / "SHA256SUMS").write_text("".join(lines), encoding="ascii", newline="")
print(f"MANIFEST_FILES {len(files)}")
