#!/usr/bin/env python3
"""
procgen_tension_20260926_run.py -- the single runner of the tension lane (session collatz-procgen-20260922,
2026-09-26).  Produces 05-knowledge/results/procgen_tension_20260926.out:

    python3 04-computation/experiments/procgen_tension_20260926_run.py \
        > 05-knowledge/results/procgen_tension_20260926.out

Every 'ok:' line is a check(...) that aborts the run on failure.  Options:
    --no-level6      skip the level-6 CP-SAT part of Q3 (about 6 minutes)
    --level7-bound   also run the level-7 upper-bound CP-SAT query (time limit 2 h; may end undecided)
    --level8-search  also search level 8 for rho_max = 41/65, then 29/46 (refined CP-SAT, 30 min each)
    --mip6           add the level-6 query "density > 5/8" to the flow-formulation MIP cross-check (slow)
Resources: one process, CP-SAT with 2 workers, peak RSS well below 700 MB.
"""
import gc
import hashlib
from fractions import Fraction
import os
import platform
import resource
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy  # noqa: E402
import ortools  # noqa: E402
import scipy  # noqa: E402

import procgen_tension_20260926_q1 as q1  # noqa: E402
import procgen_tension_20260926_q2 as q2  # noqa: E402
import procgen_tension_20260926_q3 as q3  # noqa: E402
import procgen_tension_20260926_q4 as q4  # noqa: E402
from procgen_tension_20260926_lib import claim  # noqa: E402

FILES = ["procgen_tension_20260926_lib.py", "procgen_tension_20260926_q1.py", "procgen_tension_20260926_q2.py",
         "procgen_tension_20260926_q3.py", "procgen_tension_20260926_q4.py", "procgen_tension_20260926_mip.py",
         "procgen_tension_20260926_run.py"]


def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    args = set(sys.argv[1:])
    t0 = time.time()
    print("procgen_tension_20260926 -- ranks = certificates, Christoffel maximizers, rigidity, negation duality")
    print("python %s, numpy %s, scipy %s, ortools %s, %s" % (platform.python_version(), numpy.__version__,
                                                           scipy.__version__, ortools.__version__, platform.machine()))
    for f in FILES:
        p = os.path.join(HERE, f)
        if os.path.exists(p):
            print("sha256 %s  %s" % (sha(p), f))
    print()
    for name, fn in (("Q1", q1.run), ("Q2", q2.run)):
        t = time.time()
        fn()
        gc.collect()
        print("  [%s wall %.1f s]" % (name, time.time() - t))
        print()
    t = time.time()
    q3.run(level6=("--no-level6" not in args), level7_bound=("--level7-bound" in args),
           level8_search=("--level8-search" in args))
    gc.collect()
    print("  [Q3 wall %.1f s]" % (time.time() - t))
    print()
    import procgen_tension_20260926_mip as mip
    t = time.time()
    lv = [(4, Fraction(1, 2), 'SAT'), (4, Fraction(3, 5), 'UNSAT'), (5, Fraction(3, 5), 'SAT'),
          (5, Fraction(5, 8), 'UNSAT'), (6, Fraction(3, 5), 'SAT')]
    if "--mip6" in args:
        lv.append((6, Fraction(5, 8), 'UNSAT'))
    mip.run(levels=lv)
    print("  [MIP wall %.1f s]" % (time.time() - t))
    print()
    t = time.time()
    q4.run()
    print("  [Q4 wall %.1f s]" % (time.time() - t))
    print()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        rss /= 1024.0 * 1024.0
    else:
        rss /= 1024.0
    claim(rss < 700, "peak RSS of the runner %.0f MB (< 700 MB)" % rss)
    print("total wall %.1f s" % (time.time() - t0))
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
