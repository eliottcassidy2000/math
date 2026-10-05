#!/usr/bin/env python3
"""
E5q (opus, 2026-10-04): Monte-Carlo kurtosis of the random-start q-tower at larger n (Fourier note 4k).

Exact facts: E|f_n(0)|^2 = 3^-n for every random-start q-tower (4k, PROVED).  The normalised fourth moment
K_n = 9^n E|f_n(0)|^4 is a weighted additive energy; at n <= 4 it is ~3x the i.i.d. value for every q.  Here:
S random odd starts u (40-bit), the recursion to N = 40 for q = 3, 5, 7, 11 and i.i.d. digits, and K_n estimated at
n = 10, 20, 30, 40 (heavy-tailed; standard errors reported).  Also the empirical typical rate over the window.
Usage: python collatz_phase_tower_kurtosis_mc_20261004.py [N=40] [A=40] [S=2000]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_drift_control_20261004 import run_uq
from collatz_fixed_frequency_random_digits_20261004 import run

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    S = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"Monte-Carlo kurtosis K_n = 9^n E|f_n(0)|^4 (E|f|^2 = 3^-n exactly), N={N}, A={A}, S={S} random starts")
    rng = np.random.default_rng(1)
    checkpoints = [10, 20, 30, 40]
    t0 = time.time()
    for label, q in (("q=3", 3), ("q=5", 5), ("q=7", 7), ("q=11", 11), ("iid", None)):
        F = np.zeros((S, N + 1))
        for s in range(S):
            if q is None:
                v = run(N, A, "iid", seed=20000 + s)              # already 3^(n/2)|f|
                F[s] = v
            else:
                u = int(rng.integers(1, 1 << 40)) | 1
                while u % q == 0:                                   # PRIMITIVE frequencies only: a multiple of q has
                    u = int(rng.integers(1, 1 << 40)) | 1           # |mu_hat_n(q u')| = |mu_hat_(n-1)(u')| (projectivity),
                v = run_uq(u, q, N, A)                              # which inflated the first version's moments (audit C3)
                F[s] = v * (3.0 / q) ** (np.arange(N + 1) / 2.0)
        line = f"{label:5s}:"
        for n in checkpoints:
            x = F[:, n]
            m2 = float(np.mean(x ** 2)); k = float(np.mean(x ** 4)); se = float(np.std(x ** 4) / math.sqrt(S))
            line += f"  n={n}: 3^n E|f|^2 = {m2:.3f}, K = {k:.2f} +- {se:.2f};"
        tl = np.mean(np.log(np.maximum(F[:, 10:], 1e-300)), axis=0); ns = np.arange(10, N + 1)
        line += f"  typical rate 10..{N}: {3**-0.5 * math.exp(np.polyfit(ns, tl, 1)[0]):.4f}   [{time.time()-t0:.0f}s]"
        P(line)
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
