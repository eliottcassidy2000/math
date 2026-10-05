#!/usr/bin/env python3
"""
E5l (opus, 2026-10-04): the ensemble mean square of the u-family -- is the q = 3 excess a mean-square effect or a
Jensen (typical-versus-rms) effect?

For q = 3 and q = 5, run the fixed-frequency recursion for 16 units u (odd, prime to 3q) to N and record
v_n(u) = q^(n/2) |mu_hat_n(u)| (rescaled by sqrt q: for q = 3 the incoherent and Parseval rates coincide, so a flat
ensemble mean square means "exactly incoherent"; for q = 5 the ensemble is rescaled by sqrt 5, and the typical
0.573 shows as growth (0.573/0.447)^n -- we therefore report, for q = 5, the slope of the ensemble mean square
in absolute terms).  Reported: the least-squares slope of log E_u[v_n^2] (ensemble mean square, in units of the
rescaling), the least-squares slope of E_u[log v_n] (the typical rate), both over n = 300..N, and the Jensen gap.
If for q = 3 the ensemble mean square is flat (rms rate 1/sqrt3) while the typical rate is 0.5700, the excess is a
Jensen effect of the cube-root coupling; if the mean square itself decays faster than 1/3 per level, the
coupling produces genuine mean-square anticorrelation.
Usage: python collatz_fixed_frequency_ensemble_rms_20261004.py [N=1000] [A=40]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_drift_control_20261004 import run_uq

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"ensemble (16 units) mean square versus typical rate, N={N}, A={A}, window 300..{N}")
    t0 = time.time()
    for q in (3, 5):
        units = [u for u in range(1, 80, 2) if u % 3 and u % q][:16]
        V = np.array([run_uq(u, q, N, A) for u in units])          # shape (16, N+1), rescaled by q^(n/2)
        ns = np.arange(300, N + 1)
        ms = np.mean(V[:, 300:] ** 2, axis=0)                      # ensemble mean square (rescaled)
        tl = np.mean(np.log(np.maximum(V[:, 300:], 1e-300)), axis=0)   # ensemble mean log (typical)
        sl_ms = np.polyfit(ns, np.log(ms), 1)[0]
        sl_tl = np.polyfit(ns, tl, 1)[0]
        rms_rate = q ** -0.5 * math.exp(sl_ms / 2)
        typ_rate = q ** -0.5 * math.exp(sl_tl)
        # per-unit typical rates for the scatter
        per = [q ** -0.5 * math.exp(np.polyfit(ns, np.log(np.maximum(V[i, 300:], 1e-300)), 1)[0]) for i in range(len(units))]
        P(f"q={q}: units {units}")
        P(f"   ensemble rms rate {rms_rate:.4f}   (incoherent 1/sqrt3 = 0.5774; Parseval q^(-1/2) = {q**-0.5:.4f})")
        P(f"   ensemble typical rate {typ_rate:.4f};  per-unit typical rates mean {np.mean(per):.4f} +- {np.std(per, ddof=1):.4f}")
        P(f"   Jensen gap log(rms/typical) per level = {math.log(rms_rate / typ_rate):+.4f}")
        P(f"   ensemble mean square by 100-blocks (rescaled): " + " ".join(f"{float(np.mean(ms[a:a+100])):.2e}" for a in range(0, N - 300, 100)) + f"   [{time.time()-t0:.0f}s]")
        np.save(f"collatz_fixed_frequency_ensemble_q{q}.npy", V)
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
