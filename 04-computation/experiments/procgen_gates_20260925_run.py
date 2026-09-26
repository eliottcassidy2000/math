#!/usr/bin/env python3
"""procgen_gates_20260925_run.py

Runner for the lane "cycle gates: exponential sums and equidistribution"
(session collatz-procgen-20260922, 2026-09-25).  Runs the four parts one after another (one
process at a time; each stays below 700 MB) and writes
05-knowledge/results/procgen_gates_20260925.out.  Exits nonzero if any part fails.

  python3 04-computation/experiments/procgen_gates_20260925_run.py

Part C reads the eligibility dump that part B writes to scratch/procgen_gates/partB_elig.tsv.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "05-knowledge", "results", "procgen_gates_20260925.out")
PARTS = [
    ("A", "procgen_gates_20260925_expsum.py"),
    ("B", "procgen_gates_20260925_stats.py"),
    ("C", "procgen_gates_20260925_heuristic.py"),
    ("D", "procgen_gates_20260925_switching.py"),
]


def main():
    os.makedirs(os.path.join(ROOT, "scratch", "procgen_gates"), exist_ok=True)
    chunks = []
    status = []
    for tag, fn in PARTS:
        t = time.time()
        env = dict(os.environ)
        env["MallocLargeCache"] = "0"  # macOS: return freed large blocks (keeps RSS below 700 MB)
        p = subprocess.run([sys.executable, os.path.join(HERE, fn)], capture_output=True, text=True, env=env)
        dt = time.time() - t
        memline = [ln for ln in p.stderr.splitlines() if ln.startswith("[mem] end")]
        status.append((tag, fn, p.returncode, dt, memline[-1] if memline else "[mem] n/a"))
        chunks.append("#" * 100 + f"\n# PART {tag}: {fn}\n" + "#" * 100 + "\n" + p.stdout)
        if p.returncode != 0:
            chunks.append("STDERR:\n" + p.stderr)
            break
    head = ["procgen_gates_20260925.out -- lane: cycle gates, exponential sums and equidistribution",
            "session collatz-procgen-20260922, 2026-09-25; produced by procgen_gates_20260925_run.py",
            "parts: A = DP certification, counting identity, census with SHEET/DRIFT controls;",
            "       B = full residue statistics and spectra for p <= 30 (q=3), p <= 24 (q=5);",
            "       C = heuristics (naive vs perigee-refined), bound hierarchy, archimedean limit law;",
            "       D = switching (Weyl-type) bound: carries modulo small primes, certificate at the gate modulus."]
    for tag, fn, rc, dt, mem in status:
        head.append(f"  part {tag}: exit {rc}, {dt:.0f} s, {mem}")
    with open(OUT, "w") as fh:
        fh.write("\n".join(head) + "\n\n" + "\n".join(chunks))
    print("\n".join(head))
    sys.exit(0 if all(s[2] == 0 for s in status) and len(status) == len(PARTS) else 1)


if __name__ == "__main__":
    main()
