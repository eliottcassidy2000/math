#!/usr/bin/env python3
"""
E5n (opus, 2026-10-04): large ensembles to settle mean-square versus Jensen for the Collatz cold-rate excess.

Two hundred units u (odd, prime to 3, u <= 1200) of the q = 3 tower against two hundred i.i.d.-digit towers, each to
N = 600 (A = 40).  Reported for both: the ensemble mean square E[3^n |f_n(0)|^2] by 50-blocks over n = 200..600 and
its least-squares rate (the i.i.d. model's true value is exactly 1, E5c: so its fitted rate measures the finite-
ensemble bias), the ensemble typical rate exp(E log), and the per-member typical rates.  Decisive reading: if the
i.i.d. ensemble rms rate converges to 0.5774 while the q = 3 one stays near 0.570, the cube-root coupling produces
a genuine mean-square anticorrelation; if both move together, the excess is a Jensen (typical-value) effect.
Usage: python collatz_fixed_frequency_large_ensemble_20261004.py [N=600] [A=40] [M=200]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_drift_control_20261004 import run_uq
from collatz_fixed_frequency_random_digits_20261004 import run

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    M = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"large ensembles, {M} members each, N={N}, A={A}, window 200..{N}; values 3^(n/2)|f_n(0)|")
    t0 = time.time()
    units = [u for u in range(1, 3000, 2) if u % 3][:M]
    V3 = np.array([run_uq(u, 3, N, A) for u in units]); P(f"  q=3 ensemble done [{time.time()-t0:.0f}s]")
    Vi = np.array([run(N, A, "iid", seed=9000 + s) for s in range(M)]); P(f"  iid ensemble done [{time.time()-t0:.0f}s]")
    ns = np.arange(200, N + 1)
    for label, V in (("q=3, 200 units", V3), ("iid, 200 seeds", Vi)):
        ms = np.mean(V[:, 200:] ** 2, axis=0); tl = np.mean(np.log(np.maximum(V[:, 200:], 1e-300)), axis=0)
        sl_ms = np.polyfit(ns, np.log(ms), 1)[0]; sl_tl = np.polyfit(ns, tl, 1)[0]
        per = [3 ** -0.5 * math.exp(np.polyfit(ns, np.log(np.maximum(V[i, 200:], 1e-300)), 1)[0]) for i in range(M)]
        P(f"{label}: ensemble rms rate {3**-0.5*math.exp(sl_ms/2):.4f}; ensemble typical rate {3**-0.5*math.exp(sl_tl):.4f}; per-member typical {np.mean(per):.4f} +- {np.std(per, ddof=1)/math.sqrt(M):.4f} (sd {np.std(per, ddof=1):.4f})")
        P("   ensemble mean square by 50-blocks (iid exact = 1): " + " ".join(f"{float(np.mean(ms[a:a+50])):.3f}" for a in range(0, N - 200, 50)))
        # also the median-of-members mean square and the trimmed mean (drop the top 5%) to see the tail's role
        sq = V[:, 200:] ** 2
        trimmed = np.array([np.mean(np.sort(sq[:, j])[: int(0.95 * M)]) for j in range(sq.shape[1])])
        P(f"   trimmed-mean (top 5% dropped) rms rate {3**-0.5*math.exp(np.polyfit(ns, np.log(trimmed), 1)[0]/2):.4f}; median-of-members rms rate {3**-0.5*math.exp(np.polyfit(ns, np.log(np.median(sq, axis=0)), 1)[0]/2):.4f}")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
