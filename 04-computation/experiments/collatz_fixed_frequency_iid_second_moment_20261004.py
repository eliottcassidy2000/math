#!/usr/bin/env python3
"""
E5c (opus, 2026-10-04): in the i.i.d.-digit model the phases omega(-d) = e(theta_d), theta_d = frac(R/2^d) =
0.b_d ... b_1, are pairwise uncorrelated: for d < e the top bit b_e enters theta_e - theta_d with coefficient 1/2 and
nothing else, so E e(theta_d - theta_e) = (1 + e(1/2))/2 * (...) = 0.  Hence E|f_n(k)|^2 = sum_a 4^-a E|f_(n-1)(k-a)|^2
and, by induction from f_0 = 1, E|f_n(k)|^2 = 3^-n exactly for every window position k and every n.
Check: average 3^n |f_n(0)|^2 over S seeds at n = 20, 40, 60, 80, 100 (should be 1 +- O(S^-1/2)), and the
typical (geometric-mean) value exp(E log(3^(n/2)|f_n(0)|)) against the rms 1: the Jensen gap.
Also the same statistics with the real digits of 3^-n are a single sample (printed for reference).
Usage: python collatz_fixed_frequency_iid_second_moment_20261004.py [N=100] [A=40] [S=200]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_random_digits_20261004 import run

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    S = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    t0 = time.time()
    checkpoints = [n for n in (10, 20, 40, 60, 80, 100, 150, 200) if n <= N]
    sq = {n: [] for n in checkpoints}; lg = {n: [] for n in checkpoints}
    for s in range(S):
        v = run(N, A, "iid", seed=1000 + s)          # v[n] = 3^(n/2) |f_n(0)|
        for n in checkpoints:
            sq[n].append(v[n] ** 2); lg[n].append(math.log(max(v[n], 1e-300)))
    P(f"iid-digit model, {S} seeds, A = {A}: E[3^n |f_n(0)|^2] (prediction exactly 1) and the typical value exp(E log(3^(n/2)|f_n(0)|))")
    for n in checkpoints:
        m = float(np.mean(sq[n])); se = float(np.std(sq[n]) / math.sqrt(S))
        g = math.exp(float(np.mean(lg[n])))
        P(f"  n={n:3d}: E[3^n|f|^2] = {m:.4f} +- {se:.4f}   typical 3^(n/2)|f| = {g:.4f}   => typical rate = 3^(-1/2) * {g:.4f}^(1/n) = {3**-0.5 * g**(1/n):.5f}")
    v = run(N, A, "real", q=3)
    P("real 3^-n digits (one sample): 3^(n/2)|f_n(0)| at the checkpoints: " + " ".join(f"{v[n]:.4f}" for n in checkpoints))
    P(f"  [{time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
