#!/usr/bin/env python3
"""
E7 (opus, 2026-10-04): the ridge test (b) of HYP-9176 -- the frequency-one coefficient of the 3-adic Syracuse law
to level 5000 with the rescaled window recursion (mac-mini's n <= 2500 data: largest value of 3^(n/2)|mu_hat_n(1)|
past n = 200 is 0.078 at n = 261; block rates 0.567-0.572).
Reports: 250-block medians of v_n = 3^(n/2)|mu_hat_n(1)|, the local maxima above 0.02 for n > 2000, the ratio of
each local maximum to the local trend (the 'ridge height'), least-squares rates over 300..2500, 2500..5000, and the
sup of |mu_hat_n(1)| / rho^n for rho = 0.585, 0.58, 0.5774 over n >= 200 (H1's constant C).
Truncation A = 60 (A = 40 vs 60 agree to 1e-9 relative to n = 750, E5f).
Usage: python collatz_fixed_frequency_ridges_n5000_20261004.py [N=5000] [A=60]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_rates_20261004 import run_u

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    t0 = time.time()
    v = run_u(1, N, A)
    np.save(__file__.replace(".py", "_v.npy"), v)
    P(f"frequency-one coefficient to N={N}, A={A}: 3^(n/2)|mu_hat_n(1)|   [{time.time()-t0:.0f}s]")
    P("250-block medians: " + " ".join(f"{float(np.median(v[a:a+250])):.2e}" for a in range(1, N, 250)))
    for (a_, b_) in ((300, 2500), (2500, N), (300, N)):
        ns = np.arange(a_, b_ + 1); y = np.log(np.maximum(v[a_:b_ + 1], 1e-300)) - ns * math.log(3) / 2
        P(f"least-squares rate of |mu_hat_n(1)| over {a_}..{b_}: {math.exp(np.polyfit(ns, y, 1)[0]):.4f}")
    ns = np.arange(300, N + 1); y = np.log(np.maximum(v[300:], 1e-300)); sl, ic = np.polyfit(ns, y, 1)
    trend = np.exp(sl * ns + ic)
    resid = y - np.log(trend)
    peaks = [(int(ns[i]), float(v[300 + i]), float(math.exp(resid[i]))) for i in range(1, len(ns) - 1)
             if resid[i] >= resid[i - 1] and resid[i] >= resid[i + 1] and resid[i] >= math.log(8.0)]
    P(f"local maxima of the detrended series at least 8x the trend (n, 3^(n/2)|mu_hat|, height/trend): {len(peaks)}")
    for p in peaks[:40]:
        P(f"   n={p[0]:5d}  v={p[1]:.3e}  height {p[2]:.1f}x")
    big = [(int(ns[i]), float(v[300 + i])) for i in range(len(ns)) if v[300 + i] >= 0.02]
    P(f"values of 3^(n/2)|mu_hat_n(1)| >= 0.02 for n >= 300: {len(big)}; largest: {sorted(big, key=lambda t: -t[1])[:8]}")
    for rho in (0.585, 0.58, 0.5774, 0.575):
        r = v[1:] * (3 ** -0.5 / rho) ** np.arange(1, N + 1)
        P(f"sup_n |mu_hat_n(1)|/{rho}^n over n >= 200: {r[199:].max():.3e} at n={r[199:].argmax()+200}; over n >= 2500: {r[2499:].max():.3e} at n={r[2499:].argmax()+2500}")
    P(f"[total {time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
