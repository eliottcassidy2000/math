#!/usr/bin/env python3
"""Regenerate the packet SHA-256 manifest (the manifest excludes itself)."""

from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "SHA256SUMS.txt"


def main() -> int:
    paths = sorted(
        path for path in ROOT.rglob("*")
        if path.is_file() and path != OUTPUT
    )
    lines = []
    for path in paths:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="ascii")
    print(f"MANIFEST_FILES {len(paths)}")
    print(f"SHA256SUMS_SHA256 {hashlib.sha256(OUTPUT.read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
