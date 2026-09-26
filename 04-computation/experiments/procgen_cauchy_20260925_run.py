#!/usr/bin/env python3
"""procgen_cauchy_20260925 -- runner.  Reproduces 05-knowledge/results/procgen_cauchy_20260925.out:
   python3 04-computation/experiments/procgen_cauchy_20260925_run.py > 05-knowledge/results/procgen_cauchy_20260925.out
Parts (each a separate process, run one after the other; the largest peaks at about 500 MB):
  (1) procgen_cauchy_20260925_moments.py   exact M2(L) = E[W_L^2] (L <= 16), the criticality identity, the Cauchy-Schwarz
                                            bound, Holder exponents, distribution-only optimum, tails, SHEET, DRIFT;
  (2) procgen_cauchy_20260925_syracuse.py  the fixed point (Syracuse law pi): level energies, Poisson-kernel identity,
                                            ladder autocorrelations, tail of d pi / d mu;
  (3) procgen_cauchy_20260925_majorant.py  certified growth rate theta_r (r <= 15) of ||Z_k||^2 and the proved exponent;
  (4) procgen_cauchy_20260925_harmonic.py  the harmonic (pointwise covering) bound H_L, exact L <= 13, Monte Carlo L <= 16;
  (5) procgen_cauchy_20260925_meanfield.py mean-field smoothing transform, g_q(s), Zolotarev constants, strips.
Runtime about 1 minute.  Timing and memory lines vary between runs; everything else is deterministic (fixed seeds).
"""
import os, subprocess, sys, time, resource

HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
for script in ("procgen_cauchy_20260925_moments.py", "procgen_cauchy_20260925_syracuse.py",
               "procgen_cauchy_20260925_majorant.py", "procgen_cauchy_20260925_harmonic.py",
               "procgen_cauchy_20260925_meanfield.py"):
    t1 = time.time()
    out = subprocess.run([sys.executable, os.path.join(HERE, script)], capture_output=True, text=True, check=True).stdout
    print(out, end="", flush=True)
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    rss_mb = rss / 2 ** 20 if sys.platform == "darwin" else rss / 2 ** 10
    print("[runner] %s: %.1f s; max RSS of children so far %.0f MB\n" % (script, time.time() - t1, rss_mb), flush=True)
print("[runner] total %.1f s" % (time.time() - t0))
