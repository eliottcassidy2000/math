#!/usr/bin/env python3
"""E5l control (opus, 2026-10-04): the same 16-member ensemble statistics for the i.i.d.-digit model (whose true
mean-square rate is exactly 1/sqrt3, E5c) and for the random-odd-multiplier model, to calibrate the finite-ensemble
bias of the 'ensemble rms rate' estimator used in collatz_fixed_frequency_ensemble_rms_20261004.py.
Usage: python collatz_fixed_frequency_ensemble_rms_control_20261004.py [N=1000] [A=40]"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_random_digits_20261004 import run
if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"ensemble control (16 members), N={N}, A={A}, window 300..{N}; values rescaled by 3^(n/2)")
    t0 = time.time()
    for mode, label in (("iid", "i.i.d. digits (true rms rate exactly 0.5774)"), ("randomunit", "random odd multiplier")):
        V = np.array([run(N, A, mode, seed=500 + s) for s in range(16)])
        ns = np.arange(300, N + 1)
        ms = np.mean(V[:, 300:] ** 2, axis=0); tl = np.mean(np.log(np.maximum(V[:, 300:], 1e-300)), axis=0)
        rms_rate = 3 ** -0.5 * math.exp(np.polyfit(ns, np.log(ms), 1)[0] / 2)
        typ_rate = 3 ** -0.5 * math.exp(np.polyfit(ns, tl, 1)[0])
        per = [3 ** -0.5 * math.exp(np.polyfit(ns, np.log(np.maximum(V[i, 300:], 1e-300)), 1)[0]) for i in range(16)]
        P(f"{label}: ensemble rms rate {rms_rate:.4f}; ensemble typical rate {typ_rate:.4f}; per-member typical {np.mean(per):.4f} +- {np.std(per, ddof=1):.4f}; Jensen gap {math.log(rms_rate/typ_rate):+.4f}   [{time.time()-t0:.0f}s]")
        P("   ensemble mean square by 100-blocks: " + " ".join(f"{float(np.mean(ms[a:a+100])):.2e}" for a in range(0, N - 300, 100)))
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
