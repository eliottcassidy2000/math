#!/usr/bin/env python3
"""Run the clean-room verifier and proposed primary in four exact modes."""

import hashlib
import os
from pathlib import Path
import subprocess
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PYTHON = Path(sys.executable).resolve()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def run_target(label, source, frozen=None):
    modes = (
        ("normal", [str(PYTHON), "-B", str(source)], None),
        ("optimized", [str(PYTHON), "-B", "-O", str(source)], None),
        ("isolated", [str(PYTHON), "-B", "-I", str(source)], None),
        ("hashseed", [str(PYTHON), "-B", str(source)], "314159265"),
    )
    outputs = {}
    for mode, command, seed in modes:
        environment = os.environ.copy()
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        if seed is not None:
            environment["PYTHONHASHSEED"] = seed
        completed = subprocess.run(
            command,
            cwd=ROOT,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"{label}/{mode} failed with {completed.returncode}: "
                + completed.stderr.decode("utf-8", "replace")
            )
        if completed.stderr:
            raise RuntimeError(
                f"{label}/{mode} wrote stderr: "
                + completed.stderr.decode("utf-8", "replace")
            )
        if b"\r" in completed.stdout:
            raise RuntimeError(f"{label}/{mode} output is not raw LF")
        outputs[mode] = completed.stdout
        (HERE / f"{label}.{mode}.out").write_bytes(completed.stdout)

    normal = outputs["normal"]
    if any(data != normal for data in outputs.values()):
        raise RuntimeError(f"{label} mode outputs differ")
    if frozen is not None and frozen.read_bytes() != normal:
        raise RuntimeError(f"{label} does not match frozen output {frozen}")
    (HERE / f"{label}.frozen.out").write_bytes(normal)
    source_bytes = source.read_bytes()
    if b"\r" in source_bytes:
        raise RuntimeError(f"{label} source is not raw LF")
    return {
        "source": sha256(source_bytes),
        "output": sha256(normal),
        "bytes": len(normal),
    }


def main():
    cleanroom = run_target(
        "cleanroom",
        HERE / "verify_thm4374_cleanroom.py",
    )
    primary = run_target(
        "primary",
        ROOT / "04-computation" /
        "lrc14_seventeen_step_metric_exit_observability_thm4374.py",
        ROOT / "05-knowledge" / "results" /
        "lrc14_seventeen_step_metric_exit_observability_thm4374.out",
    )
    print("STATUS PASS")
    print("MODES normal optimized isolated hashseed=314159265")
    for label, result in (("cleanroom", cleanroom), ("primary", primary)):
        print(
            f"TARGET {label} source_sha256={result['source']} "
            f"output_sha256={result['output']} bytes={result['bytes']} raw_lf=yes"
        )


if __name__ == "__main__":
    main()
