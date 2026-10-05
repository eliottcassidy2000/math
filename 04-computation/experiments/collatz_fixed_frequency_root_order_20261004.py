#!/usr/bin/env python3
"""E5i (opus, 2026-10-04): the typical cold-frequency rate r(q) for composite and 3-power root orders.
q = 9, 15, 21, 27 (and 3 again as control), four units each, N = 1500: does r(q) track the root order q, or the
3-adic content of q (9, 27 = powers of 3; 15, 21 = 3 x odd)?  Same recursion as the drift control (phases =
binary digits of -u q^-n; rescaled by the incoherent rate sqrt3 so that a flat sequence = exactly incoherent).
Usage: python collatz_fixed_frequency_root_order_20261004.py [N=1500] [A=40]"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_drift_control_20261004 import phases_u_q

def run_uq3(u, q, N, A):
    w = math.sqrt(3.0) * 2.0 ** (-np.arange(1, A + 1))   # rescale by sqrt3, not sqrt q
    prev_lo = -N * A - A
    prev = np.ones(-prev_lo + 1, dtype=np.complex128)
    v0 = np.zeros(N + 1); v0[0] = 1.0
    for n in range(1, N + 1):
        lo = -(N - n) * A
        ph = phases_u_q(u, q, n, lo - A, 0)
        Pp = ph * prev[(lo - A) - prev_lo: (0 - prev_lo) + 1]
        cur = np.zeros(-lo + 1, dtype=np.complex128)
        for a in range(1, A + 1):
            cur += w[a - 1] * Pp[A - a: A - a + (-lo + 1)]
        prev, prev_lo = cur, lo
        v0[n] = abs(cur[-lo])
    return v0

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"root-order dependence of the typical cold-frequency rate, N={N}, A={A}")
    t0 = time.time()
    for q in (9, 15, 21, 27, 3):
        rates = []
        for u in (1, 5, 7, 11) if q % 5 else (1, 7, 11, 13):
            v = run_uq3(u, q, N, A)
            ns = np.arange(300, N + 1); y = np.log(np.maximum(v[300:N + 1], 1e-300)) - ns * math.log(3) / 2
            r = 3 ** -0.5 * math.exp(np.polyfit(ns, y, 1)[0]); rates.append(r)
            P(f"  q={q:2d} u={u:2d}: rate 300..{N} = {r:.4f}   [{time.time()-t0:.0f}s]")
        a = np.array(rates)
        P(f"  q={q:2d}: mean {a.mean():.4f} +- {a.std(ddof=1):.4f}   (i.i.d. 0.5723 +- 0.0008; q=3 six units 0.5700 +- 0.0009)")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
