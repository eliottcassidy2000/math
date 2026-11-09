#!/usr/bin/env python3
"""Replay in a temporary copy, never overwrite the frozen research evidence."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_generator(path):
    spec = importlib.util.spec_from_file_location("generator", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true",
                        help="recompute all exhaustive scans, including order six")
    args = parser.parse_args()
    commands = []
    comparisons = []
    with tempfile.TemporaryDirectory(prefix="digraph234-replay-") as temp:
        mirror = Path(temp) / "repo"
        work = mirror / "05-knowledge" / "results" / ROOT.name
        shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns("__pycache__"))
        citation = Path("01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md")
        (mirror / citation).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT.parents[2] / citation, mirror / citation)
        census, audit = work / "census", work / "audit"

        def run(cmd, cwd, output=None):
            result = subprocess.run(cmd, cwd=cwd, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            require(result.returncode == 0,
                    f"Command failed: {cmd}\n{result.stdout}\n{result.stderr}")
            commands.append({"command": cmd, "directory": cwd.name, "status": "PASS"})
            if output:
                (cwd / output).write_text(result.stdout)

        if args.full:
            for name in ("scan.py", "analyze.py", "classify3.py"):
                run([sys.executable, name], census)
            run(["g++", "-O3", "-std=c++17", "scan6.cpp", "-o", "scan6"], census)
            run(["./scan6"], census)
            run([sys.executable, "enrich6.py"], census)
            run(["g++", "-O3", "-std=c++17", "general_scan.cpp", "-o", "general_scan"], census)
            run(["./general_scan"], census)
            run(["g++", "-O3", "-std=c++17", "audit.cpp", "-o", "audit"], audit)
            for base, orders, prefix in [(3, range(1, 7), "semicomplete"),
                                         (2, range(3, 8), "tournament")]:
                for n in orders:
                    run(["./audit", str(base), str(n)], audit, f"{prefix}-{n}.txt")
            run(["g++", "-O3", "-std=c++17", "audit_all.cpp", "-o", "audit_all"], audit)
            for n in range(1, 6):
                run(["./audit_all", "4", str(n)], audit, f"unrestricted-{n}.txt")

        run([sys.executable, "generator.py", "--validate"], census)
        run([sys.executable, "verify.py"], audit)
        run([sys.executable, "verify_all.py"], audit)
        run([sys.executable, "check.py"], work / "final_review", "check-output.txt")

        generator = load_generator(census / "generator.py")
        seeds = [([[int(v) for v in row] for row in rows.split("/")], target)
                 for rows, target in [
                     ("0111/0011/1001/1000", 7),
                     ("01111/00110/00011/10001/11000", 21),
                     ("011110/001111/000111/100011/000001/100000", 21),
                     ("010/101/100", 3)]]
        for matrix, target in seeds:
            require(generator.hp(matrix) == target, "Explicit seed count failed")
            for i in range(len(matrix)):
                matrix[i][i] = 1
            require(generator.hp(matrix) == target, "Diagonal invariance failed")
        # A hostile control must be rejected by the same gate used above.
        rejected = False
        try:
            require(generator.hp(seeds[0][0]) == 21, "wrong target")
        except RuntimeError:
            rejected = True
        require(rejected, "Wrong-target negative control was not rejected")

        paths = ["census/generator_validation.json", "final_review/validation.json",
                 "final_review/check-output.txt"]
        if args.full:
            paths += ["census/" + name for name in [
                "results.json", "strong_h3.json", "classes6.csv", "classes6.json",
                "results6.csv", "general_census.csv", "general_verification.csv"]]
            paths += [str(p.relative_to(ROOT)) for p in (ROOT / "audit").glob("hist-*.csv")]
            paths += [str(p.relative_to(ROOT)) for pattern in
                      ("semicomplete-*.txt", "tournament-*.txt", "unrestricted-*.txt")
                      for p in (ROOT / "audit").glob(pattern)]
        for rel in sorted(paths):
            expected, actual = ROOT / rel, work / rel
            same = (json.loads(expected.read_text()) == json.loads(actual.read_text())
                    if expected.suffix == ".json" else expected.read_bytes() == actual.read_bytes())
            require(same, "Frozen-output mismatch: " + rel)
            comparisons.append({"path": rel, "comparison": "JSON content" if
                                expected.suffix == ".json" else "exact bytes", "status": "PASS"})

    print(json.dumps({"status": "PASS", "mode": "full" if args.full else "quick",
                      "commands": commands, "comparisons": comparisons,
                      "explicit_seed_and_diagonal_checks": 8,
                      "wrong_target_negative_control": "rejected"}, indent=2))


if __name__ == "__main__":
    main()
