"""Freeze every non-executable packet file except SHA256SUMS itself."""

from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parent
paths = sorted(
    path for path in ROOT.rglob("*")
    if path.is_file() and path.name != "SHA256SUMS" and path.suffix != ".exe"
)
(ROOT / "SHA256SUMS").write_text(
    "".join(
        f"{sha256(path.read_bytes()).hexdigest()}  {path.relative_to(ROOT).as_posix()}\n"
        for path in paths
    ),
    encoding="ascii",
)
print(f"manifest_files={len(paths)}")
