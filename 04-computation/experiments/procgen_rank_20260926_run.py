#!/usr/bin/env python3
"""
procgen_rank_20260926_run.py -- the single runner of the rank lane (session collatz-procgen-20260922, 2026-09-26):
what a Lyapunov function (rank) for Collatz must look like; infinite banks of 2-adic valuation counters,
adaptive centers, heights.

    python3 -u 04-computation/experiments/procgen_rank_20260926_run.py \
        > 05-knowledge/results/procgen_rank_20260926.out

Every 'ok:' line is a check(...) that aborts the run on failure; unlabelled lines print tables quoted in the note.
Requirements: numpy, ortools (CP-SAT, 2 workers).  One process; peak RSS is checked to stay below 700 MB.
"""
import hashlib
import os
import platform
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy  # noqa: E402
import ortools  # noqa: E402

import procgen_rank_20260926_cube as cube  # noqa: E402
import procgen_rank_20260926_q1 as q1  # noqa: E402
import procgen_rank_20260926_q2 as q2  # noqa: E402
import procgen_rank_20260926_q3 as q3  # noqa: E402
from procgen_rank_20260926_lib import check, peak_rss_mb, say  # noqa: E402

FILES = ["procgen_rank_20260926_lib.py", "procgen_rank_20260926_q1.py", "procgen_rank_20260926_cube.py",
         "procgen_rank_20260926_q2.py", "procgen_rank_20260926_q3.py", "procgen_rank_20260926_run.py"]


def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    t0 = time.time()
    say("procgen_rank_20260926 -- two-place Lyapunov functions: forced charges, banks, adaptive centers, heights")
    say("python %s, numpy %s, ortools %s, %s" % (platform.python_version(), numpy.__version__, ortools.__version__,
                                                platform.machine()))
    for f in FILES:
        say("sha256 %s  %s" % (sha(os.path.join(HERE, f)), f))
    say()
    for name, fn in (("Q1", q1.run), ("cube", cube.run), ("Q2", q2.run), ("Q3", q3.run)):
        t = time.time()
        fn()
        say("  [%s wall %.1f s]" % (name, time.time() - t))
        say()
    rss = peak_rss_mb()
    check(rss < 700, "peak RSS %.0f MB < 700 MB (getrusage)" % rss)
    say("  [total wall %.1f s]" % (time.time() - t0))


if __name__ == "__main__":
    main()
