#!/usr/bin/env python3
"""procgen-bridges 2026-09-23: run the four parts one at a time and write the combined output.

    python3 04-computation/experiments/procgen_bridges_20260923_run.py [--quick]

Writes 05-knowledge/results/procgen_bridges_20260923.out.  Stdout of every part is deterministic; timing goes
to stderr.  One process at a time; peak memory below 200 MB; full run about 5 minutes.
"""
import hashlib
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EXP = REPO / "04-computation" / "experiments"
OUT = REPO / "05-knowledge" / "results" / "procgen_bridges_20260923.out"
PARTS = ["procgen_bridges_20260923_phi.py", "procgen_bridges_20260923_adelic.py",
         "procgen_bridges_20260923_defects.py", "procgen_bridges_20260923_twoblock.py"]


def main():
    quick = "--quick" in sys.argv
    lines = ["procgen-bridges 2026-09-23: LRC(14), AMM 12592 and the Collatz/Periodicity-Conjecture work.",
             "Session collatz-procgen-20260922 (mac-mini), bridges lane.  Mode: " + ("QUICK" if quick else "FULL"),
             "Script SHA-256 (raw bytes):"]
    for p in PARTS + [Path(__file__).name]:
        h = hashlib.sha256((EXP / p).read_bytes()).hexdigest()
        lines.append(f"  {h}  04-computation/experiments/{p}")
    lines.append("")
    body = []
    for p in PARTS:
        t0 = time.time()
        cmd = [sys.executable, str(EXP / p)] + (["--quick"] if quick else [])
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO)
        print(f"[{time.time() - t0:7.1f}s] {p} exit {r.returncode}", file=sys.stderr, flush=True)
        if r.returncode != 0:
            print(r.stderr, file=sys.stderr)
            raise SystemExit(f"{p} failed")
        body.append(r.stdout)
    OUT.write_text("\n".join(lines) + "\n" + "\n".join(body), encoding="utf-8")
    print(f"wrote {OUT.relative_to(REPO)}", file=sys.stderr)


if __name__ == "__main__":
    main()
