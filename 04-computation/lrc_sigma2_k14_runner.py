#!/usr/bin/env python3
"""Standalone k=14 runner: writes result to results file. Detached, no timeout."""
import sys, time, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from lrc_sigma2_scan import run

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "05-knowledge", "results", "lrc_sigma2_k14.out")
OUT = os.path.abspath(OUT)

k = 14
t = time.time()
sig, wit, surv, fl, med = run(k)
elapsed = time.time() - t
line = (f"k={k} survivors={surv} time={elapsed:.2f}s sigma2={sig} "
        f"wit={wit} fl={fl} med={med}\n")
with open(OUT, "w") as f:
    f.write(line)
    f.flush()
print(line)
