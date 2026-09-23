#!/usr/bin/env python3
"""collatz_procgen_20260922_ladder_run.py -- reproduce collatz_procgen_20260922_ladder.out.

Runs, in order (about 12 minutes on an M-series laptop; peak memory about 0.7 GB, one process at a time):
  1. collatz_procgen_20260922_ladder_ballot.py      no-choice counts, q=3 dimension fit, mu_q enclosures
  2. collatz_procgen_20260922_ladder_f2x.c          F_2[x] sibling: identities, census deg<=24, parity bijection
  3. collatz_procgen_20260922_ladder_mahler.py      THM-3848 safe language re-counted three ways
  4. collatz_procgen_20260922_ladder_choice.py      choice ladder table (lane-one DP program) + cross-checks
  5. collatz_procgen_20260922_ladder_choice_mc.py   exact per-class DFS, Monte Carlo to m=64, first-moment indices
usage: python3 collatz_procgen_20260922_ladder_run.py [--quick] [--out PATH]
"""
import os
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
P = "collatz_procgen_20260922_ladder"


def main():
    quick = "--quick" in sys.argv
    out = os.path.join(ROOT, "05-knowledge", "results", P + ".out")
    if "--out" in sys.argv:
        out = sys.argv[sys.argv.index("--out") + 1]
    tmp = tempfile.mkdtemp(prefix="ladder_run_")
    f2x = os.path.join(tmp, "f2x")
    subprocess.run(["clang", "-O3", "-o", f2x, os.path.join(HERE, P + "_f2x.c"), "-lm"], check=True)
    steps = [
        ("ballot", [sys.executable, os.path.join(HERE, P + "_ballot.py")] + (["20000"] if quick else [])),
        ("f2x", [f2x, "20" if quick else "24", "20"]),
        ("mahler", [sys.executable, os.path.join(HERE, P + "_mahler.py")]),
        ("choice", [sys.executable, os.path.join(HERE, P + "_choice.py")]),
        ("choice_mc", [sys.executable, os.path.join(HERE, P + "_choice_mc.py")] + (["--quick"] if quick else [])),
    ]
    status = 0
    with open(out, "w") as fh:
        fh.write(f"# {P}.out -- produced by {P}_run.py{' --quick' if quick else ''}\n")
        fh.write("# sibling dimension ladder: see 05-knowledge/results/collatz_procgen_20260922_sibling_dimension_ladder.md\n\n")
        for name, cmd in steps:
            t0 = time.time()
            r = subprocess.run(cmd, capture_output=True, text=True)
            fh.write(r.stdout)
            if r.returncode != 0:
                status = 1
                fh.write(f"!!! step {name} exited with code {r.returncode}\n{r.stderr}\n")
            fh.write(f"[step {name}: {time.time() - t0:.1f}s, exit {r.returncode}]\n\n")
            fh.flush()
            print(f"step {name}: exit {r.returncode}, {time.time() - t0:.1f}s", flush=True)
        fh.write(f"LADDER RUN TOTAL: {'ALL STEPS EXITED 0' if status == 0 else 'SOME STEP FAILED'}\n")
    subprocess.run(["rm", "-rf", tmp])
    print(f"wrote {out}")
    return status


if __name__ == "__main__":
    sys.exit(main())
