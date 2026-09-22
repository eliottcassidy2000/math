"""Replay the scoped Collatz blueprint audit, never a convergence claim.

Normal and optimized Python controls precede a clean Lean root/axiom audit.
Invalidate stale PASS before running; record FAIL on any failed step.
"""
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "05-knowledge/results"
STEM = "collatz_blueprint_20260921"
MANIFEST = RESULTS / f"{STEM}_manifest.json"
TRANSCRIPT = RESULTS / f"{STEM}_verification.out"
LEAN = ROOT / "04-computation/lean/CollatzBlueprintAudit"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def save_manifest(data):
    MANIFEST.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8", newline="\n")


def normalized_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    report = {"status": "RUNNING", "collatz_convergence_proved": False,
              "scope": "Scoped proofs, finite controls, and Lean reductions only.",
              "hash_basis": "artifact bytes with CRLF normalized to LF", "runs": []}
    save_manifest(report)
    log = []

    def run(args, label):
        proc = subprocess.run(args, cwd=ROOT, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=180)
        log.extend([label, proc.stdout, proc.stderr, f"exit_code={proc.returncode}\n"])
        TRANSCRIPT.write_text("\n".join(log), encoding="utf-8", newline="\n")
        require(proc.returncode == 0, f"command failed: {label}")
        return proc.stdout

    try:
        for lane in ("affine", "energy"):
            script = ROOT / "04-computation/experiments" / f"{STEM}_{lane}.py"
            output = script.with_suffix(".json")
            relative = script.relative_to(ROOT).as_posix()
            stdout = run([sys.executable, str(script)], f"python {relative}")
            if lane == "affine":
                require(json.loads(stdout).get("all_checks_passed") is True,
                        "affine output lacks success flag")
                output.write_text(stdout, encoding="utf-8", newline="\n")
            normal = output.read_bytes()
            stdout = run([sys.executable, "-O", str(script)], f"python -O {relative}")
            if lane == "affine":
                output.write_text(stdout, encoding="utf-8", newline="\n")
            require(normal == output.read_bytes(), f"optimized replay differs: {lane}")
            report["runs"].append({"lane": lane, "ordinary_pass": True,
                                   "optimized_pass": True, "byte_identical": True,
                                   "command": f"python {relative}",
                                   "certificate": output.relative_to(ROOT).as_posix()})

        run([sys.executable, str(LEAN / "verify.py")],
            "python 04-computation/lean/CollatzBlueprintAudit/verify.py")
        lean_report = json.loads((LEAN / "verification.json").read_text(encoding="utf-8"))
        require(lean_report.get("status") == "PASS", "Lean verification did not pass")
        report["lean_verification"] = "04-computation/lean/CollatzBlueprintAudit/verification.json"
        paths = [Path(__file__), TRANSCRIPT,
                 RESULTS / f"{STEM}_synthesis.md",
                 RESULTS / f"{STEM}_affine.md",
                 RESULTS / f"{STEM}_energy.md",
                 ROOT / "05-knowledge/reference/COLLATZ-BLUEPRINT-2026-09-21-SOURCE.md"]
        for lane in ("affine", "energy"):
            paths += [ROOT / "04-computation/experiments" / f"{STEM}_{lane}{ext}"
                      for ext in (".py", ".json")]
        paths += [p for p in LEAN.rglob("*")
                  if p.is_file() and ".lake" not in p.relative_to(LEAN).parts
                  and "__pycache__" not in p.relative_to(LEAN).parts]
        report["artifacts"] = [{"path": p.relative_to(ROOT).as_posix(),
                                 "sha256_lf": normalized_hash(p)}
                                for p in sorted(set(paths))]
        report["status"] = "PASS"
        save_manifest(report)
        print("PASS: both exact lanes match under Python optimization; Lean root build and axiom audit passed.")
        print("Collatz convergence remains OPEN.")
    except Exception as error:
        report["status"] = "FAIL"
        report["error"] = str(error)
        save_manifest(report)
        raise


if __name__ == "__main__":
    main()
