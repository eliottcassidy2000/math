#!/usr/bin/env python3
"""procgen_mlr_20260930_run.py -- runner for the multiplicative lonely runner lane (mlr), 2026-09-30.

Runs the three parts one at a time (each a separate `nice`d process, so peak memory is per part):
  A  procgen_mlr_20260930_lrc.py    LRC side: kappa of boxes, x2x3 lonely spectrum certificate,
                                    discrete loneliness, 6-free counting criterion
  B  procgen_mlr_20260930_crowd.py  crowded side: Triangle/Uniqueness lemmas, coset bound, census
  C  procgen_mlr_20260930_gates.py  Collatz side: gate spectra vs runner numerators, switching bound,
                                    major/minor arcs and the Parseval barrier
Output goes to stdout only (redirect to 05-knowledge/results/procgen_mlr_20260930.out).
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PARTS = [("A", "procgen_mlr_20260930_lrc.py", "PART A CHECKS PASSED"),
         ("B", "procgen_mlr_20260930_crowd.py", "PART B CHECKS PASSED"),
         ("C", "procgen_mlr_20260930_gates.py", "PART C CHECKS PASSED")]


def main():
    t0 = time.time()
    print("procgen_mlr_20260930 -- multiplicative lonely runner (session collatz-procgen-20260922)")
    print("command: python3 04-computation/experiments/procgen_mlr_20260930_run.py")
    sys.stdout.flush()
    env = dict(os.environ)
    env["MallocLargeCache"] = "0"
    for name, script, token in PARTS:
        print(f"\n############ PART {name}: {script} ############")
        sys.stdout.flush()
        cmd = ["nice", "-n", "10", sys.executable, "-u", os.path.join(HERE, script)]
        proc = subprocess.run(cmd, cwd=HERE, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        sys.stdout.write(proc.stdout)
        sys.stdout.flush()
        if proc.returncode != 0 or token not in proc.stdout:
            print(f"\nPART {name} FAILED (exit code {proc.returncode})")
            sys.exit(1)
    print(f"\ntotal wall time {time.time() - t0:.0f}s")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
