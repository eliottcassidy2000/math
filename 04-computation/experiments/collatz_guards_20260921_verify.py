"""Replay scoped Collatz guard research and the extended Lean package.

The older blueprint runner refreshes its manifest and performs the clean
Lean audit. This runner pins that verification along with the three new
independent mathematical lanes. It never certifies global convergence.
"""
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "05-knowledge/results"
STEM = "collatz_guards_20260921"
MANIFEST = RESULTS / f"{STEM}_manifest.json"
TRANSCRIPT = RESULTS / f"{STEM}_verification.out"
LEAN = ROOT / "04-computation/lean/CollatzBlueprintAudit"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def normalized_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    report = {"status": "RUNNING", "collatz_convergence_proved": False,
              "scope": "Scoped written proofs, finite certificates, and Lean reductions only",
              "hash_basis": "artifact bytes with CRLF normalized to LF", "runs": []}
    log = []

    def save():
        MANIFEST.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8", newline="\n")

    def run(args, label):
        proc = subprocess.run(args, cwd=ROOT, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=240)
        log.extend([label, proc.stdout, proc.stderr, f"exit_code={proc.returncode}\n"])
        TRANSCRIPT.write_text("\n".join(log), encoding="utf-8", newline="\n")
        require(proc.returncode == 0, f"command failed: {label}")
        return proc.stdout

    save()
    try:
        baseline = ROOT / "04-computation/experiments/collatz_blueprint_20260921_verify.py"
        run([sys.executable, str(baseline)], "python " + baseline.relative_to(ROOT).as_posix())
        baseline_manifest = RESULTS / "collatz_blueprint_20260921_manifest.json"
        require(json.loads(baseline_manifest.read_text(encoding="utf-8"))["status"] == "PASS",
                "baseline audit failed")

        for lane in ("valves", "squarefree", "discrepancy"):
            script = ROOT / "04-computation/experiments" / f"{STEM}_{lane}.py"
            output = script.with_suffix(".json")
            relative = script.relative_to(ROOT).as_posix()
            for optimized in (False, True):
                args = [sys.executable] + (["-O"] if optimized else []) + [str(script)]
                stdout = run(args, ("python -O " if optimized else "python ") + relative)
                if lane == "valves":
                    require(json.loads(stdout).get("all_checks_passed") is True,
                            "valve output lacks success flag")
                    output.write_text(stdout, encoding="utf-8", newline="\n")
                data = json.loads(output.read_text(encoding="utf-8"))
                if lane == "squarefree":
                    require(data.get("checks", "").startswith("PASS"), "squarefree checks failed")
                if lane == "discrepancy":
                    require(data.get("status") == "PASS", "discrepancy checks failed")
                if not optimized:
                    ordinary_bytes = output.read_bytes()
                else:
                    require(ordinary_bytes == output.read_bytes(), "optimized output differs: " + lane)
            report["runs"].append({"lane": lane, "ordinary_pass": True,
                                   "optimized_pass": True, "byte_identical": True,
                                   "command": "python " + relative,
                                   "certificate": output.relative_to(ROOT).as_posix()})

        lean = json.loads((LEAN / "verification.json").read_text(encoding="utf-8"))
        require(lean["status"] == "PASS" and lean["theorems_audited"] == 35, "Lean audit scope changed")
        for name, axioms in lean["axioms"].items():
            require(set(axioms) <= {"propext", "Quot.sound"}, "unexpected axiom: " + name)
        report["lean"] = {"root_import": lean["root_import"], "theorems_audited": 35,
                          "verification": "04-computation/lean/CollatzBlueprintAudit/verification.json",
                          "new_written_analytic_and_sieve_proofs_formalized": False}
        paths = [Path(__file__), TRANSCRIPT, baseline, baseline_manifest,
                 RESULTS / "collatz_blueprint_20260921_verification.out",
                 RESULTS / f"{STEM}_synthesis.md",
                 ROOT / "05-knowledge/reference/COLLATZ-GUARDS-2026-09-21-SOURCE.md"]
        for lane in ("valves", "squarefree", "discrepancy"):
            paths.append(RESULTS / f"{STEM}_{lane}.md")
            paths += [ROOT / "04-computation/experiments" / f"{STEM}_{lane}{ext}"
                      for ext in (".py", ".json")]
        paths += [p for p in LEAN.rglob("*") if p.is_file()
                  and ".lake" not in p.relative_to(LEAN).parts
                  and "__pycache__" not in p.relative_to(LEAN).parts]
        report["artifacts"] = [{"path": p.relative_to(ROOT).as_posix(), "sha256_lf": normalized_hash(p)}
                               for p in sorted(set(paths))]
        report["status"] = "PASS"
        save()
        print("PASS: three exact lanes match under -O; prior controls and 35-theorem Lean audit pass.")
        print("Written proofs have independent audits; the new analytic/sieve theorems are not yet in Lean.")
        print("Global Collatz convergence remains OPEN.")
    except Exception as error:
        report["status"] = "FAIL"
        report["error"] = str(error)
        save()
        raise


if __name__ == "__main__":
    main()
