"""Rebuild the package from clean and reject forbidden proof dependencies.

Steps:
  1. scan every Lean source for proof escapes (`sorry`, `admit`, `native_decide`,
     `decide +native`, `ofReduceBool`, custom `axiom` declarations);
  2. check that every `theorem` of the library modules appears in AxiomAudit.lean;
  3. `lake clean`, then build the modules ONE AT A TIME in dependency order
     (so at most one Lean process runs; the lakefile caps each at ~1.4 GB), then
     the library root;
  4. run AxiomAudit.lean (which imports only the root) and require every audited
     theorem to depend on a subset of {propext, Quot.sound};
  5. write verification.json (status, versions, axioms per theorem, sha256 of the
     sources, per-module build times, peak child RSS) and verification.log.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
LIB = "ProcgenSelfieEdim"
MODULES = (
    "ListLemmas",
    "Hypercube",
    "EdimQ6",
    "CountingBound",
    "Stabilizer",
    "KeyBlocks",
    "EdimQ7Data",
    "EdimQ7P0",
    "EdimQ7P1",
    "EdimQ7",
    "EdimQ8Data",
    "EdimQ8P0",
    "EdimQ8P1",
    "EdimQ8P2",
    "EdimQ8P3",
    "EdimQ8",
    "EdimQ9Data",
    "EdimQ9P0",
    "EdimQ9P1",
    "EdimQ9P2",
    "EdimQ9P3",
    "EdimQ9P4",
    "EdimQ9P5",
    "EdimQ9P6",
    "EdimQ9P7",
    "EdimQ9P8",
    "EdimQ9P9",
    "EdimQ9P10",
    "EdimQ9P11",
    "EdimQ9",
    "Tournament",
    "HamPath",
    "ArcParity",
    "AltSum",
    "TournamentCode",
    "SelfieFinite",
    "Shaved",
    "CycleCount",
    "Circulant",
    "HPExist",
    "SwitchSum",
    "Redei",
    "ConstantH",
    "ParityBreak",
    "CollatzDrop",
)
SOURCE_PATHS = (
    "lean-toolchain",
    "lakefile.toml",
    f"{LIB}.lean",
    *(f"{LIB}/{m}.lean" for m in MODULES),
    "AxiomAudit.lean",
    "verify.py",
    "gen_certificates.py",
)
ALLOWED_AXIOMS = {"propext", "Quot.sound"}
EXPECTED_LEAN = "4.30.0"
FORBIDDEN = re.compile(r"\b(sorry|admit|native_decide|sorryAx|ofReduceBool)\b|\+native")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def executable(name: str) -> str:
    elan = Path.home() / ".elan" / "bin" / name
    if elan.exists():
        return str(elan)
    found = shutil.which(name)
    require(found is not None, f"Missing executable: {name}")
    return str(found)


def peak_child_rss_mb() -> float:
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    # ru_maxrss is in bytes on macOS and in kilobytes on Linux.
    return peak / (1024 * 1024) if sys.platform == "darwin" else peak / 1024


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    (ROOT / "verification.json").write_text('{"status": "RUNNING"}\n', encoding="utf-8")
    transcript: list[str] = []

    def run(args: list[str]) -> tuple[str, float]:
        start = time.monotonic()
        result = subprocess.run(
            [executable(args[0]), *args[1:]],
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            errors="replace",
            capture_output=True,
            check=False,
            env={**os.environ, "LEAN_NUM_THREADS": "1"},
        )
        elapsed = time.monotonic() - start
        output = result.stdout + result.stderr
        transcript.append(
            f"$ {' '.join(args)}\n{output}exit_code={result.returncode} seconds={elapsed:.1f}\n"
        )
        (ROOT / "verification.log").write_text("\n".join(transcript), encoding="utf-8")
        require(result.returncode == 0, f"Failed: {' '.join(args)}\n{output}")
        return output, elapsed

    # 1. proof escapes
    for name in SOURCE_PATHS:
        if name.endswith(".lean"):
            source = (ROOT / name).read_text(encoding="utf-8")
            require(not FORBIDDEN.search(source), f"Forbidden proof escape in {name}")
            require(not re.search(r"^\s*axiom\s", source, re.MULTILINE), f"Custom axiom in {name}")

    # 2. audit coverage
    names: list[str] = []
    for module in MODULES:
        source = (ROOT / LIB / f"{module}.lean").read_text(encoding="utf-8")
        names.extend(re.findall(r"^theorem\s+([A-Za-z_][A-Za-z0-9_']*)", source, re.MULTILINE))
    require(len(set(names)) == len(names), "Duplicate theorem names in audited modules")
    audit_source = (ROOT / "AxiomAudit.lean").read_text(encoding="utf-8")
    require(audit_source.startswith(f"import {LIB}\n"), "Axiom audit must import the library root")
    for name in names:
        require(f"#print axioms {LIB}.{name}\n" in audit_source,
                f"Theorem is missing from axiom audit: {name}")

    # 3. clean, sequential build
    version, _ = run(["lean", "--version"])
    require(EXPECTED_LEAN in version, f"Unexpected Lean version: {version.strip()}")
    run(["lake", "clean"])
    build_seconds: dict[str, float] = {}
    for module in MODULES:
        _, seconds = run(["lake", "build", f"{LIB}.{module}"])
        build_seconds[module] = round(seconds, 1)
    _, seconds = run(["lake", "build"])
    build_seconds["(root)"] = round(seconds, 1)

    # 4. axiom audit
    output, audit_seconds = run(["lake", "env", "lean", "AxiomAudit.lean"])
    axiom_map: dict[str, list[str]] = {}
    for name in names:
        qualified = f"{LIB}.{name}"
        match = re.search(
            rf"'{re.escape(qualified)}' (does not depend on any axioms|depends on axioms: \[(.*?)\])",
            output,
        )
        require(match is not None, f"Missing audit output: {qualified}")
        axioms = [] if match.group(2) is None else [x.strip() for x in match.group(2).split(",")]
        require(set(axioms) <= ALLOWED_AXIOMS, f"Unexpected axioms for {qualified}: {axioms}")
        axiom_map[qualified] = axioms

    histogram: dict[str, int] = {}
    for axioms in axiom_map.values():
        key = ", ".join(axioms) if axioms else "(none)"
        histogram[key] = histogram.get(key, 0) + 1

    record = {
        "status": "PASS",
        "scope": (
            "THM-4525 (Q_6 resolving 15-set, L1, L2, L3, L4, explicit sets edim_m(Q_7) <= 19, Q_8 <= 26, Q_9 <= 38), THM-4524 (gauge theorem A1, verified "
            "Hamiltonian-path enumerator, double counting and arc parity, QR_7 / QR_7 - v counts, "
            "no all-odd tournament for N = 3..5, C1 for circulants of odd order, HP existence, no HP-covered tournament with a universal arc, A2 loop-set sum 2 N!, Redei's theorem and the unconditional parity corollaries, A3(e) constant-H classes need N = 2^k, parity-break number 2 for all-even tournaments), THM-4526 (H_4, H_5, #copies(H_n) = H - n hc), S15 Collatz drop "
            "multiplicity and label map; not the lower bound edim_m(Q_6) >= 15"
        ),
        "lean_version": version.strip(),
        "root_import": LIB,
        "external_packages": [],
        "theorems_audited": len(names),
        "axiom_histogram": histogram,
        "axioms": axiom_map,
        "build_seconds": build_seconds,
        "audit_seconds": round(audit_seconds, 1),
        "peak_child_rss_mb": round(peak_child_rss_mb(), 1),
        "sha256": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in SOURCE_PATHS
        },
    }
    (ROOT / "verification.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"PASS: clean sequential build and axiom audit for {len(names)} theorems "
        f"({sum(build_seconds.values()):.0f} s build, peak child RSS {record['peak_child_rss_mb']} MB)."
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        (ROOT / "verification.json").write_text(
            json.dumps({"status": "FAIL", "error": str(error)}, indent=2) + "\n",
            encoding="utf-8",
        )
        raise
