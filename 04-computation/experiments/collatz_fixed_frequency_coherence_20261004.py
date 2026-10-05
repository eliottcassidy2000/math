#!/usr/bin/env python3
"""
E5m (opus, 2026-10-04): the per-level coherence ratio of the window recursion, along one run.

At each level, f_(n+1)(0) = sum_a 2^-a omega_(n+1)(-a) f_n(-a).  Define the coherence ratio
    kappa_n = |f_(n+1)(0)|^2 / sum_a 4^-a |f_n(-a)|^2,
the output energy at exponent 0 divided by the incoherent prediction computed from the ACTUAL input vector.  Under
random phases E[kappa_n] = 1 exactly, whatever the input.  Along one run: the arithmetic mean of kappa_n over levels
is the mean-square coherence (1 = incoherent), and exp(E log kappa_n) is the typical coherence; the typical rate of
|f_n(0)| is sqrt(typical coherence) times the rate of the input norm.  Real digits of 3^-n (u = 1, 5, 7), of 5^-n and
7^-n (u = 1), and i.i.d. digits (three seeds), N = 1500, A = 40, window 300..1500.
Usage: python collatz_fixed_frequency_coherence_20261004.py [N=1500] [A=40]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_random_digits_20261004 import phases_from_R

def run_coh(N, A, mode, q=3, u=1, seed=0):
    rng = np.random.default_rng(seed)
    wts = 2.0 ** (-np.arange(1, A + 1))
    w = math.sqrt(3.0) * wts
    prev_lo = -N * A - A
    prev = np.ones(-prev_lo + 1, dtype=np.complex128)
    kap = np.zeros(N + 1)
    for n in range(1, N + 1):
        lo = -(N - n) * A
        M = -(lo - A)
        if mode == "real":
            R = (-u * pow(q, -n, 1 << (M + 1))) % (1 << (M + 1))
        else:
            R = int.from_bytes(rng.bytes((M + 1) // 8 + 1), "little") % (1 << (M + 1)) | 1
        ph = np.empty(M + 1, dtype=np.complex128); ph[:M] = phases_from_R(R, M); ph[M] = 1.0
        # incoherent prediction for exponent 0 from the input: sum_a 4^-a |prev(-a)|^2 (prev indexed from prev_lo)
        inp = prev[(0 - prev_lo) - A: (0 - prev_lo)]          # prev(-A) .. prev(-1)
        inco = float(np.sum((wts[::-1] ** 2) * np.abs(inp) ** 2))   # weights 4^-a matched to prev(-a)
        Pp = ph * prev[(lo - A) - prev_lo: (0 - prev_lo) + 1]
        cur = np.zeros(-lo + 1, dtype=np.complex128)
        for a in range(1, A + 1):
            cur += w[a - 1] * Pp[A - a: A - a + (-lo + 1)]
        out0 = abs(cur[-lo]) ** 2 / 3.0                       # undo the sqrt3 rescaling of this level
        kap[n] = out0 / inco if inco > 0 else float("nan")
        prev, prev_lo = cur, lo
    return kap

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"per-level coherence ratio kappa_n (1 = incoherent), N={N}, A={A}, window 300..{N}")
    P("source | mean kappa (mean-square coherence) | exp(mean log kappa) (typical coherence) | median | P(kappa > 4) | implied typical rate factor sqrt(typ)")
    t0 = time.time()
    runs = [("real q=3 u=1", "real", 3, 1, 0), ("real q=3 u=5", "real", 3, 5, 0), ("real q=3 u=7", "real", 3, 7, 0),
            ("real q=5 u=1", "real", 5, 1, 0), ("real q=7 u=1", "real", 7, 1, 0),
            ("iid seed=0", "iid", 3, 1, 0), ("iid seed=1", "iid", 3, 1, 1), ("iid seed=2", "iid", 3, 1, 2)]
    for label, mode, q, u, seed in runs:
        k = run_coh(N, A, mode, q=q, u=u, seed=seed)[300:]
        k = k[np.isfinite(k)]
        P(f"{label:14s} | {np.mean(k):.4f} | {math.exp(np.mean(np.log(np.maximum(k, 1e-300)))):.4f} | {np.median(k):.4f} | {np.mean(k > 4):.4f} | {math.sqrt(math.exp(np.mean(np.log(np.maximum(k, 1e-300))))):.4f}   [{time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
