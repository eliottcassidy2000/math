#!/usr/bin/env python3
"""Closed replay verifier for THM-4332's two exact implementations."""

from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "04-computation/lrc14_fixed_pool_implication_rigidity_thm4332.py"
INDEPENDENT = ROOT / "04-computation/lrc14_fixed_pool_implication_rigidity_thm4332_independent_audit.py"
PRIMARY_OUT = ROOT / "05-knowledge/results/lrc14_fixed_pool_implication_rigidity_thm4332.out"
INDEPENDENT_OUT = ROOT / "05-knowledge/results/lrc14_fixed_pool_implication_rigidity_thm4332_independent_audit.out"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(script: Path, optimized: bool = False) -> bytes:
    command = [sys.executable]
    if optimized:
        command.append("-O")
    command.extend([str(script), "591"])
    completed = subprocess.run(command, cwd=ROOT, check=True, capture_output=True)
    if completed.stderr:
        raise RuntimeError(f"unexpected stderr from {script.name}: {completed.stderr!r}")
    return completed.stdout.replace(b"\r\n", b"\n")


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    expected_primary = PRIMARY_OUT.read_bytes().replace(b"\r\n", b"\n")
    expected_independent = INDEPENDENT_OUT.read_bytes().replace(b"\r\n", b"\n")
    got_primary = run(PRIMARY)
    got_independent = run(INDEPENDENT)
    got_independent_o = run(INDEPENDENT, optimized=True)
    if got_primary != expected_primary:
        raise RuntimeError("primary replay differs from frozen output")
    if got_independent != expected_independent:
        raise RuntimeError("independent replay differs from frozen output")
    if got_independent_o != expected_independent:
        raise RuntimeError("optimized independent replay differs from frozen output")
    for path in (PRIMARY, PRIMARY_OUT, INDEPENDENT, INDEPENDENT_OUT):
        print(f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}")
    print("THM4332_REPLAY PASS")


if __name__ == "__main__":
    main()
