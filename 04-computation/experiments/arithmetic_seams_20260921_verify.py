"""Bounded master replay for four arithmetic-seams lanes and the Lean package.

Usage from any working directory:
  python /path/to/repo/04-computation/experiments/arithmetic_seams_20260921_verify.py
All paths are derived from this file, not from a particular machine.
Checks use exceptions and remain active under python -O.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
STEM = "arithmetic_seams_20260921"
EXPERIMENTS = Path("04-computation/experiments")
RESULTS = Path("05-knowledge/results")
MASTER = EXPERIMENTS / f"{STEM}_verify.py"
MANIFEST = RESULTS / f"{STEM}_manifest.json"
SYNTHESIS = RESULTS / f"{STEM}_synthesis.md"
LEAN = Path("04-computation/lean/ArithmeticSeams")
LANES = ("operations", "tournaments", "primes", "dynamics")
LEAN_INPUTS = (
    "lean-toolchain", "lakefile.toml", "lake-manifest.json", ".gitignore",
    "ArithmeticSeams.lean", "ArithmeticSeams/Basic.lean", "AxiomAudit.lean",
    "verify.py", "README.md",
)
ALLOWED_AXIOMS = {"propext", "Quot.sound"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def relative(path: Path) -> str:
    return path.as_posix()


def raw_hash(path: Path) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def lf_hash(path: Path) -> str:
    raw = (ROOT / path).read_bytes()
    raw.decode("utf-8")  # Every pinned artifact must actually be UTF-8 text.
    return hashlib.sha256(raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")).hexdigest()


def write_manifest(record: dict) -> None:
    (ROOT / MANIFEST).write_text(json.dumps(record, indent=2, sort_keys=True)+"\n",
                                 encoding="utf-8", newline="\n")


def decoded_json(path: Path) -> dict:
    value = json.loads((ROOT / path).read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"Expected JSON object: {path}")
    return value


def marker_present(value: dict) -> bool:
    return value.get("status") == "PASS" or bool(
        re.match(r"^PASS(?:\b|;)", str(value.get("checks", ""))))


def run_python(script: Path, optimized: bool, cwd: Path, record: dict) -> None:
    command = [sys.executable, *(["-O"] if optimized else []), str(ROOT / script)]
    shown = ["python", *(["-O"] if optimized else []), relative(script)]
    print("RUN: " + " ".join(shown), flush=True)
    result = subprocess.run(command, cwd=ROOT / cwd, capture_output=True, text=True,
                            encoding="utf-8", errors="replace", timeout=180, check=False)
    output = result.stdout + result.stderr
    # Printed paths are diagnostic only; the stored command/cwd and artifacts are relative.
    portable_output = output
    root_spellings = {str(ROOT), ROOT.as_posix(), json.dumps(str(ROOT))[1:-1],
                      json.dumps(ROOT.as_posix())[1:-1]}
    for spelling in sorted(root_spellings, key=len, reverse=True):
        portable_output = portable_output.replace(spelling, "<repo>")
    event = {"command": shown, "cwd": relative(cwd), "exit_code": result.returncode,
             "output": portable_output.strip()}
    record["commands"].append(event)
    write_manifest(record)
    require(result.returncode == 0, f"Command failed: {' '.join(shown)}\n{portable_output}")
    require(re.search(r"\bPASS\b", output) is not None,
            f"Missing explicit PASS/check marker: {' '.join(shown)}")


def check_embedded_source_hashes(value: dict, script: Path) -> list[str]:
    checked = []
    for key in ("source_sha256", "script_sha256"):
        if key in value:
            require(value[key] == raw_hash(script), f"Stale embedded source hash: {script}:{key}")
            checked.append(key)
    return checked


def run_lane(lane: str, record: dict) -> list[Path]:
    script = EXPERIMENTS / f"{STEM}_{lane}.py"
    output = script.with_suffix(".json")
    decoded = []
    markers = []
    for optimized in (False, True):
        # A no-op or incomplete script cannot accidentally pass with stale output.
        (ROOT / output).write_text('{"status": "MASTER_REPLAY_RUNNING"}\n', encoding="utf-8")
        run_python(script, optimized, Path("."), record)
        value = decoded_json(output)
        require(marker_present(value), f"Lane JSON lacks a PASS/check marker: {lane}")
        markers.append(check_embedded_source_hashes(value, script))
        decoded.append(value)
    require(decoded[0] == decoded[1], f"Normal/-O decoded JSON mismatch: {lane}")
    value = decoded[1]
    record["lanes"][lane] = {
        "status": "PASS", "normal_optimized_decoded_json_equal": True,
        "script": relative(script), "certificate": relative(output),
        "note": relative(RESULTS / f"{STEM}_{lane}.md"),
        "lane_status": value.get("status"), "lane_checks": value.get("checks"),
        "embedded_raw_source_hash_fields_checked": markers[0],
        "exact_universes": value.get("universes", value.get("universe", {
            "declared_in_certificate": relative(output),
            "sections": sorted(k for k in value if k not in ("source_sha256", "script_sha256")),
            "scope": "Use the finite ranges and full control lists in these certificate sections; no extrapolation.",
        })),
        "decoded_json_sha256": hashlib.sha256(
            json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest(),
    }
    write_manifest(record)
    return [output]


def run_lean(record: dict) -> list[Path]:
    run_python(LEAN / "verify.py", False, LEAN, record)
    report_path = LEAN / "verification.json"
    report = decoded_json(report_path)
    require(report.get("status") == "PASS", "Lean verifier did not finish with PASS")
    require(report.get("root_import") == "ArithmeticSeams", "Unexpected Lean root")
    require(report.get("external_packages") == [], "Lean package gained an external dependency")
    axioms = report.get("axioms")
    require(isinstance(axioms, dict) and axioms, "Missing theorem axiom report")
    require(report.get("theorems_audited") == len(axioms), "Lean theorem count mismatch")
    for theorem, dependencies in axioms.items():
        require(isinstance(dependencies, list) and set(dependencies) <= ALLOWED_AXIOMS,
                f"Unapproved Lean proof dependencies: {theorem}")
    for name, digest in report.get("sha256", {}).items():
        target = Path(name)
        require(not target.is_absolute() and ".." not in target.parts, "Unsafe Lean hash path")
        require(raw_hash(LEAN / target) == digest, f"Stale Lean source hash: {name}")
    require(report.get("sha256"), "Missing Lean source hashes")
    record["lean"] = {
        "status": "PASS", "root_import": "ArithmeticSeams", "external_packages": [],
        "lean_version": report["lean_version"], "theorems_audited": report["theorems_audited"],
        "axioms": axioms, "source_hashes_rechecked": True,
        "report": relative(report_path), "log": relative(LEAN / "verification.log"),
        "scope": report["scope"],
    }
    write_manifest(record)
    return [report_path, LEAN / "verification.log"]


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    record = {
        "status": "RUNNING", "started_utc": now(), "commands": [], "lanes": {},
        "scope": "Finite exact replay of four arithmetic-seams lanes and the shifted-arithmetic Lean root.",
        "global_collatz_convergence_proved": False, "new_sphere_theorem_proved": False,
        "hash_basis": "SHA-256 of UTF-8 file bytes after CRLF and lone CR are normalized to LF.",
        "embedded_source_hash_basis": "Each lane's embedded hash is independently checked against its raw script bytes.",
        "manifest_self_hash": "Excluded to avoid a circular hash; the master verifier source is pinned.",
    }
    write_manifest(record)
    try:
        inputs = [MASTER]
        for lane in LANES:
            inputs.extend([EXPERIMENTS / f"{STEM}_{lane}.py", RESULTS / f"{STEM}_{lane}.md"])
        inputs.extend(LEAN / name for name in LEAN_INPUTS)
        synthesis_present = (ROOT / SYNTHESIS).is_file()
        record["synthesis_present_at_invocation"] = synthesis_present
        if synthesis_present:
            inputs.append(SYNTHESIS)
        for path in inputs:
            require((ROOT / path).is_file(), f"Missing input: {path}")
        frozen_inputs = {relative(path): lf_hash(path) for path in inputs}
        outputs = []
        for lane in LANES:
            outputs.extend(run_lane(lane, record))
        outputs.extend(run_lean(record))
        for path in inputs:
            require(lf_hash(path) == frozen_inputs[relative(path)],
                    f"Input changed during replay; rerun on frozen files: {path}")
        record["artifacts_sha256_lf"] = {
            relative(path): lf_hash(path) for path in sorted(set(inputs + outputs), key=relative)}
        require(len(record["lanes"]) == len(LANES) and record.get("lean", {}).get("status") == "PASS",
                "Missing completed lane or Lean report")
        record["status"] = "PASS"
        record["finished_utc"] = now()
        write_manifest(record)
        print("PASS: four normal/-O lane pairs match; clean Lean root and axiom audit pass; artifact hashes pinned.")
    except Exception as error:
        record["status"] = "FAIL"
        record["finished_utc"] = now()
        record["error"] = f"{type(error).__name__}: {error}"
        write_manifest(record)
        raise


if __name__ == "__main__":
    main()
