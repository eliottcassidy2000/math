"""Rebuild the library root and reject forbidden proof dependencies."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SOURCE_PATHS = (
    "lean-toolchain",
    "lakefile.toml",
    "CollatzBlueprintAudit.lean",
    "CollatzBlueprintAudit/Basic.lean",
    "AxiomAudit.lean",
    "verify.py",
)
ALLOWED_AXIOMS = {"propext", "Quot.sound"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    (ROOT / "verification.json").write_text('{"status": "RUNNING"}\n', encoding="utf-8")
    transcript: list[str] = []

    def run(args: list[str]) -> str:
        executable = shutil.which(args[0])
        require(executable is not None, f"Missing executable: {args[0]}")
        result = subprocess.run(
            [executable, *args[1:]],
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            errors="replace",
            capture_output=True,
            check=False,
        )
        output = result.stdout + result.stderr
        transcript.append(f"$ {' '.join(args)}\n{output}exit_code={result.returncode}\n")
        (ROOT / "verification.log").write_text("\n".join(transcript), encoding="utf-8")
        require(result.returncode == 0, f"Failed: {' '.join(args)}\n{output}")
        return output

    for name in SOURCE_PATHS:
        if name.endswith(".lean"):
            source = (ROOT / name).read_text(encoding="utf-8")
            require(not re.search(r"\b(sorry|admit|native_decide|sorryAx)\b", source),
                    f"Forbidden proof escape in {name}")
            require(not re.search(r"^\s*axiom\s", source, re.MULTILINE),
                    f"Custom axiom in {name}")

    basic = (ROOT / "CollatzBlueprintAudit/Basic.lean").read_text(encoding="utf-8")
    names = re.findall(r"^theorem\s+(\w+)", basic, re.MULTILINE)
    audit_source = (ROOT / "AxiomAudit.lean").read_text(encoding="utf-8")
    for name in names:
        require(f"#print axioms CollatzBlueprintAudit.{name}" in audit_source,
                f"Theorem is missing from axiom audit: {name}")
    require("import CollatzBlueprintAudit\n" in audit_source,
            "Axiom audit must import the public library root")

    version = run(["lean", "--version"]).strip()
    run(["lake", "clean"])
    run(["lake", "build"])
    output = run(["lake", "env", "lean", "AxiomAudit.lean"])

    axiom_map: dict[str, list[str]] = {}
    for name in names:
        qualified = f"CollatzBlueprintAudit.{name}"
        match = re.search(
            rf"'{re.escape(qualified)}' (does not depend on any axioms|depends on axioms: \[(.*?)\])",
            output,
        )
        require(match is not None, f"Missing audit output: {qualified}")
        axioms = [] if match.group(2) is None else [x.strip() for x in match.group(2).split(",")]
        require(set(axioms) <= ALLOWED_AXIOMS, f"Unexpected axioms for {qualified}: {axioms}")
        axiom_map[qualified] = axioms

    record = {
        "status": "PASS",
        "scope": "Premise refutation, descent equivalence, finite controls; not Collatz convergence",
        "lean_version": version,
        "root_import": "CollatzBlueprintAudit",
        "external_packages": [],
        "theorems_audited": len(names),
        "axioms": axiom_map,
        "sha256": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in SOURCE_PATHS
        },
    }
    (ROOT / "verification.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"PASS: clean root build and axiom audit for {len(names)} theorems.")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        (ROOT / "verification.json").write_text(
            json.dumps({"status": "FAIL", "error": str(error)}, indent=2) + "\n",
            encoding="utf-8",
        )
        raise
