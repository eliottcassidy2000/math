#!/usr/bin/env python3
"""Freeze every non-executable packet file except SHA256SUMS itself."""

from __future__ import annotations

import hashlib
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    names = sorted(
        str(path.relative_to(root)).replace("\\", "/")
        for path in root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS"
        and not path.name.endswith(".exe")
    )
    lines = [f"{hashlib.sha256((root / name).read_bytes()).hexdigest()}  {name}\n"
             for name in names]
    (root / "SHA256SUMS").write_bytes("".join(lines).encode("ascii"))


if __name__ == "__main__":
    main()
