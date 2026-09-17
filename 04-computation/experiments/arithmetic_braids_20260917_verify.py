"""Replay the four arithmetic-braids lanes and pin normalized artifact hashes.

Run from any directory with standard-library Python. Ordinary execution is
mandatory for the summand lane, whose controls intentionally use assert.
The other three lanes are additionally replayed with -O and compared.
"""
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
PREFIX = "arithmetic_braids_20260917"


def normalized(payload):
    return payload.replace(b"\r\n", b"\n")


def main():
    runs = []
    for lane in ("collatz", "summand", "divisors", "geometry"):
        script = f"04-computation/experiments/{PREFIX}_{lane}.py"
        command = [sys.executable, script]
        result = subprocess.run(command, cwd=ROOT, capture_output=True, check=True)
        output = normalized(result.stdout)
        if lane == "collatz":
            # The script prints an executor-specific absolute output path.
            # Its JSON result, rather than that cosmetic path, is pinned below.
            record = json.loads(output)
            record["output"] = f"05-knowledge/results/{PREFIX}_collatz.json"
            output = (json.dumps(record, indent=2) + "\n").encode()
        if lane != "summand":
            optimized = subprocess.run([sys.executable, "-O", script], cwd=ROOT,
                                       capture_output=True, check=True)
            if normalized(optimized.stdout) != normalized(result.stdout):
                raise RuntimeError(f"optimized transcript differs: {lane}")
        target = f"05-knowledge/results/{PREFIX}_{lane}.out"
        (ROOT / target).write_bytes(output)
        runs.append({"lane": lane, "command": f"python3 {script}",
                     "ordinary_pass": True, "optimized_pass": lane != "summand",
                     "optimized_note": "explicit checks retained" if lane != "summand"
                     else "not run; assert-based controls require ordinary Python",
                     "transcript": target})
        print(f"{lane}: PASS")
    artifacts = []
    paths = {Path(__file__).relative_to(ROOT).as_posix()}
    for row in runs:
        lane = row["lane"]
        paths.update({f"04-computation/experiments/{PREFIX}_{lane}.py",
                      f"05-knowledge/results/{PREFIX}_{lane}.md", row["transcript"]})
    paths.update({f"05-knowledge/results/{PREFIX}_synthesis.md",
                  f"05-knowledge/results/{PREFIX}_collatz.json",
                  f"04-computation/experiments/{PREFIX}_divisors.json"})
    for relative in sorted(paths):
        payload = normalized((ROOT / relative).read_bytes())
        artifacts.append({"path": relative, "sha256_lf": sha256(payload).hexdigest(),
                          "bytes_lf": len(payload)})
    manifest = {"status": "FINITE-EXACT replay of stated universes; proofs in notes",
                "hash_basis": "UTF-8 file bytes with CRLF normalized to LF",
                "runs": runs, "artifacts": artifacts,
                "scope": "No global Collatz, prime-pair, or Thue classification claim"}
    target = ROOT / f"05-knowledge/results/{PREFIX}_manifest.json"
    target.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"manifest: {target.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
