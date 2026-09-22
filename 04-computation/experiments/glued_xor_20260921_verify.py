"""Portable, optimized-safe replay of the four glued-number-line lanes.

Finite certificates verify the stated universes; infinite claims depend on
the separately pinned written proofs and explicitly cited sparsity input.
No new general Lean formalization is asserted by this replay.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
STEM = "glued_xor_20260921"
EXPERIMENTS = Path("04-computation/experiments")
RESULTS = Path("05-knowledge/results")
MANIFEST = RESULTS / f"{STEM}_manifest.json"
LANES = ("conjugacies", "tournaments", "density", "blueprint")
SOURCE = Path("05-knowledge/reference/GLUED-AFFINE-BLUEPRINT-2026-09-21-SOURCE.md")
INHERITED_LEAN = Path("04-computation/lean/CollatzBlueprintAudit")
LEAN_INVENTORY = ("AxiomAudit.lean", "CollatzBlueprintAudit.lean",
                  "CollatzBlueprintAudit/Basic.lean",
                  "CollatzBlueprintAudit/ClockAudit.lean", "README.md")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def stamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def digest(path: Path, *, lf: bool = True) -> str:
    data = (ROOT / path).read_bytes()
    data.decode("utf-8")
    if lf:
        data = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


def save(record: dict) -> None:
    (ROOT / MANIFEST).write_text(json.dumps(record, indent=2, sort_keys=True)+"\n",
                                 encoding="utf-8", newline="\n")


def check_marker(value: dict, lane: str) -> None:
    require(value.get("status") == "PASS" or
            str(value.get("status", "")).startswith("FINITE-EXACT"),
            f"Missing final certificate marker: {lane}")
    if lane == "conjugacies":
        require(value["known_cycle_core"]["tagged_states"] == 20 and
                value["known_cycle_core"]["centralizer_size"] == 1568,
                "Conjugacy universe changed")
    elif lane == "tournaments":
        require(value["universe"]["masks"] == [0, 63] and
                value["universe"]["gauge_pair_checks"] == 4096,
                "Tournament universe changed")
    elif lane == "density":
        require(value["CRT_inputs_checked"] == 93026 and
                value["squarefree_prefix"]["X"] == 600000,
                "Density universe changed")
    else:
        require(value["checks_passed"] == 120097 and
                value["orbit_prefixes_checked"] == 18264 and
                value["formal_attribution"]["general_discrepancy_formalized"] is False,
                "Blueprint universe/formal scope changed")


def replay(lane: str, record: dict) -> Path:
    script = EXPERIMENTS / f"{STEM}_{lane}.py"
    output = script.with_suffix(".json")
    values = []
    checked_fields = []
    for optimized in (False, True):
        # Stale saved PASS files cannot substitute for fresh output.
        (ROOT / output).write_text('{"status":"MASTER_REPLAY_RUNNING"}\n',
                                   encoding="utf-8", newline="\n")
        command = [sys.executable, *(["-O"] if optimized else []), str(ROOT / script)]
        if lane == "conjugacies":
            command.extend(["--output", str(ROOT / output)])
        shown = ["python", *(["-O"] if optimized else []), script.as_posix()]
        if lane == "conjugacies":
            shown.extend(["--output", output.as_posix()])
        print("RUN: " + " ".join(shown), flush=True)
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                                encoding="utf-8", errors="replace", timeout=180,
                                check=False)
        diagnostic = result.stdout + result.stderr
        spellings = {str(ROOT), ROOT.as_posix(), json.dumps(str(ROOT))[1:-1]}
        for spelling in sorted(spellings, key=len, reverse=True):
            diagnostic = diagnostic.replace(spelling, "<repo>")
        record["commands"].append({"command": shown, "cwd": ".",
                                    "exit_code": result.returncode,
                                    "output": diagnostic.strip()})
        save(record)
        require(result.returncode == 0, f"Replay failed: {lane}\n{diagnostic}")
        value = json.loads((ROOT / output).read_text(encoding="utf-8"))
        require(isinstance(value, dict), f"Expected object certificate: {lane}")
        check_marker(value, lane)
        fields = []
        for name in ("source_sha256", "script_sha256", "source_sha256_lf"):
            if name in value:
                require(value[name] == digest(script, lf=name.endswith("_lf")),
                        f"Stale embedded source hash: {lane}:{name}")
                fields.append(name)
        require(fields, f"Missing source binding: {lane}")
        checked_fields = fields
        values.append(value)
    require(values[0] == values[1], f"Normal/-O decoded JSON differs: {lane}")
    value = values[1]
    record["lanes"][lane] = {
        "status": "PASS", "normal_optimized_decoded_json_equal": True,
        "script": script.as_posix(), "certificate": output.as_posix(),
        "note": (RESULTS / f"{STEM}_{lane}.md").as_posix(),
        "embedded_source_hash_fields_checked": checked_fields,
        "finite_universe": value.get("universe", {
            "specified_in_certificate": output.as_posix(),
            "sections": sorted(value), "no_extrapolation": True}),
        "decoded_json_sha256": hashlib.sha256(json.dumps(
            value, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest(),
    }
    save(record)
    return output


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    record = {
        "status": "RUNNING", "started_utc": stamp(), "commands": [], "lanes": {},
        "scope": "Four finite exact lanes, written proof notes, and inherited formal inventory.",
        "global_collatz_convergence_proved": False,
        "all_signed_cycles_classified": False,
        "new_general_lean_formalization": False,
        "cited_input": {
            "source": "Garcia--Tal, Acta Arithmetica 90 (1999), 245--250",
            "url": "https://matwbn.icm.edu.pl/ksiazki/aa/aa90/aa9033.pdf",
            "pin": "Proposition 1, Lemma 3, equation (6); Proposition 1 cites Heppner",
            "role": "Quantitative orbit sparsity; reciprocal summability and limits are derived in the written note.",
            "formalized_or_reproved_here": False,
        },
        "hash_basis": "SHA-256 of UTF-8 bytes with CRLF and lone CR normalized to LF; embedded hashes checked in their own declared basis.",
        "manifest_self_hash": "Excluded to avoid a circular hash; verifier source is pinned.",
    }
    save(record)
    try:
        inputs = [EXPERIMENTS / f"{STEM}_verify.py", SOURCE,
                  RESULTS / f"{STEM}_synthesis.md"]
        for lane in LANES:
            inputs.extend([EXPERIMENTS / f"{STEM}_{lane}.py",
                           RESULTS / f"{STEM}_{lane}.md"])
        inputs.extend(INHERITED_LEAN / name for name in LEAN_INVENTORY)
        frozen = {path.as_posix(): digest(path) for path in inputs}
        outputs = [replay(lane, record) for lane in LANES]
        for path in inputs:
            require(digest(path) == frozen[path.as_posix()],
                    f"Input changed during replay; rerun after freezing: {path}")
        for lane, details in record["lanes"].items():
            path = Path(details["certificate"])
            current = json.loads((ROOT / path).read_text(encoding="utf-8"))
            current_hash = hashlib.sha256(json.dumps(
                current, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
            require(current_hash == details["decoded_json_sha256"],
                    f"Output changed after replay; rerun after freezing: {lane}")
        record["artifacts_sha256_lf"] = {
            path.as_posix(): digest(path) for path in sorted(set(inputs + outputs))}
        require(len(record["lanes"]) == 4, "Missing completed lane")
        record["status"] = "PASS"
        record["finished_utc"] = stamp()
        save(record)
        print("PASS: four fresh normal/-O pairs agree; proof notes and source hashes pinned.")
    except Exception as error:
        record.update(status="FAIL", finished_utc=stamp(),
                      error=f"{type(error).__name__}: {error}")
        save(record)
        raise


if __name__ == "__main__":
    main()
