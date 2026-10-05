#!/usr/bin/env python3
"""
E5g (opus, 2026-10-04): is the typical cold-frequency rate a function of the multiplier q?

E5b/E5f: for q = 3 the six units u = 1, 5, 7, 11, 13, 17 give rates 0.569-0.571 (N = 1500), all below the twelve
i.i.d. runs (0.5723 +- 0.0008), while single units at q = 5, 7, 11, 13 gave 0.5736, 0.5727, 0.5722, 0.5724.  The
level-to-level phase relation is e(theta_(n+1,d))^q = e(theta_(n,d)) (a q-th root with the branch fixed by the
residue), so the x q^-n family (u varying) and the Q^-n family (Q varying) are different random models.
Here: six units each for q = 5 and q = 7 at N = 1500, to see whether r(q) is a smooth function of q.
Usage: python collatz_fixed_frequency_unit_families_20261004.py [N=1500] [A=40]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_drift_control_20261004 import run_uq

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"unit families: typical rates of |mu_hat_n(u)| for the q-adic law of qx+1, N={N}, A={A}, windows 300..N and 750..N")
    t0 = time.time()
    for q, units in ((5, (1, 2, 3, 7, 11, 13)), (7, (1, 2, 3, 5, 11, 13))):
        rates = []
        for u in units:
            v = run_uq(u, q, N, A)
            rs = []
            for a_ in (300, 750):
                ns = np.arange(a_, N + 1); y = np.log(np.maximum(v[a_:N + 1], 1e-300)) - ns * math.log(q) / 2
                rs.append(math.exp(np.polyfit(ns, y, 1)[0]))
            rates.append(rs)
            P(f"  q={q} u={u:2d}: rate 300..{N} = {rs[0]:.4f}   750..{N} = {rs[1]:.4f}   [{time.time()-t0:.0f}s]")
        arr = np.array(rates)
        P(f"  q={q}: mean rate 300..{N} = {arr[:,0].mean():.4f} +- {arr[:,0].std(ddof=1):.4f}; 750..{N} = {arr[:,1].mean():.4f} +- {arr[:,1].std(ddof=1):.4f}   (i.i.d. digits: 0.5723 +- 0.0008; q = 3 units: 0.5700 +- 0.0009)")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
