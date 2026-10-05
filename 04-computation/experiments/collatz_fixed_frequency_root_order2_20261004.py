#!/usr/bin/env python3
"""E5j (opus, 2026-10-04): root orders q = 25, 49, 45, 33, 5 (control) -- prime squares and 3 x odd versus the
Parseval rate q^(-1/2) and the incoherent rate 1/sqrt3.  Three units each, N = 1200."""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_root_order_20261004 import run_uq3
if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"root orders q = 25, 49, 45, 33, 5: typical rate of |mu_hat_n(u)|, N={N}, A={A}")
    t0 = time.time()
    for q in (25, 49, 45, 33, 5):
        rates = []
        for u in (1, 2, 7) if q % 7 else (1, 2, 11):
            v = run_uq3(u, q, N, A)
            ns = np.arange(250, N + 1); y = np.log(np.maximum(v[250:N + 1], 1e-300)) - ns * math.log(3) / 2
            r = 3 ** -0.5 * math.exp(np.polyfit(ns, y, 1)[0]); rates.append(r)
            P(f"  q={q:2d} u={u:2d}: rate 250..{N} = {r:.4f}   (Parseval q^(-1/2) = {q**-0.5:.4f}, incoherent 0.5774)   [{time.time()-t0:.0f}s]")
        a = np.array(rates)
        P(f"  q={q:2d}: mean {a.mean():.4f} +- {a.std(ddof=1):.4f}")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
