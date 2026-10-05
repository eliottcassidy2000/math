#!/usr/bin/env python3
"""
E5e (opus, 2026-10-04): is the extra deficit of the real 3^-n digits a Jensen-gap effect (larger fluctuations)?

For the real digits of q^-n (q = 3, 5, 7, 11) and for i.i.d. digits (4 seeds) we save the whole sequence
v_n = 3^(n/2)|f_n(0)| to N and compute over the window 300..N: the least-squares rate, the variance of the
detrended log (log v_n minus the fitted line), the largest upward excursion above the fitted line, and the
Cesaro mean of v_n^2 / (fitted trend)^2.  A lower typical rate with the same mean square must come with larger
log-fluctuations (Jensen); if q = 3 has the largest log-variance among the digit sources, its extra deficit is
the ridge structure, not a different second moment.
Usage: python collatz_fixed_frequency_fluctuations_20261004.py [N=1500] [A=40]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_random_digits_20261004 import run

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"fluctuation structure of 3^(n/2)|f_n(0)|, N={N}, A={A}, window 300..N")
    P("source | rate | log-variance (detrended) | max excursion (log units) at n | Cesaro mean of v^2/trend^2 | n with v >= 3x trend")
    rows = []
    t0 = time.time()
    sources = [("real q=3", "real", 3, 0), ("real q=5", "real", 5, 0), ("real q=7", "real", 7, 0), ("real q=11", "real", 11, 0)] + \
              [(f"iid seed={s}", "iid", 3, 200 + s) for s in range(4)]
    for label, mode, q, seed in sources:
        v = run(N, A, mode, q=q, seed=seed)
        np.save(f"collatz_fixed_frequency_fluct_{label.replace(' ', '_').replace('=', '')}.npy", v)
        ns = np.arange(300, N + 1); y = np.log(np.maximum(v[300:], 1e-300))
        sl, ic = np.polyfit(ns, y, 1)
        resid = y - (sl * ns + ic)
        rate = 3 ** -0.5 * math.exp(sl)
        trend = np.exp(sl * ns + ic)
        ces = float(np.mean((v[300:] / trend) ** 2))
        imax = int(np.argmax(resid))
        big = [int(n) for n, r in zip(ns, resid) if r >= math.log(3.0)]
        P(f"{label:14s} | {rate:.4f} | {float(np.var(resid)):.3f} | {float(resid[imax]):+.2f} at n={300 + imax} | {ces:.2f} | {len(big)} ({big[:8]})")
    P(f"  [{time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
