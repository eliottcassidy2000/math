#!/usr/bin/env python3
"""
E5d (opus, 2026-10-04): mean-square (window-energy) rate versus typical rate, real digits of 3^-n against i.i.d. digits.

Along the window recursion f_n(k) = sum_a 2^-a omega_n(k-a) f_{n-1}(k-a) (k <= 0), define the window energy
    e_n = sum_{k=-W}^{0} |f_n(k)|^2   (W fixed; these are the coefficients at the 3-adically fixed frequencies 2^k).
In the i.i.d.-digit model E e_n = (W+1) 3^-n exactly (E5c).  Here we print, for the real digits of -3^-n, of -5^-n,
and for i.i.d. digits, the per-level energy ratios e_n/e_{n-1} (median and mean over blocks), their least-squares rate,
and the typical rate of |f_n(0)| for comparison.  If the energy rate is 1/3 per level for the real digits too, the
real-versus-random difference of the typical rates is a difference of fluctuation structure (the Jensen gap), not of
the second moment; THM-4521(3) measured the full-cycle analogue (ratio 1/3 + covariance +0.0002..0.0006, n <= 12).
Usage: python collatz_fixed_frequency_window_energy_20261004.py [N=600] [A=40] [W=100]
"""
import sys, math, time
import numpy as np
sys.path.insert(0, __file__.rsplit("\\", 1)[0] if "\\" in __file__ else ".")
from collatz_fixed_frequency_random_digits_20261004 import phases_from_R

def run_energy(N, A, W, mode, q=3, seed=0):
    rng = np.random.default_rng(seed)
    w = math.sqrt(3.0) * 2.0 ** (-np.arange(1, A + 1))
    prev_lo = -N * A - A
    prev = np.ones(-prev_lo + 1, dtype=np.complex128)
    e = np.zeros(N + 1); v0 = np.zeros(N + 1); e[0] = W + 1; v0[0] = 1.0
    for n in range(1, N + 1):
        lo = -(N - n) * A
        M = -(lo - A)
        if mode == "real":
            R = (-pow(q, -n, 1 << (M + 1))) % (1 << (M + 1))
        else:
            R = int.from_bytes(rng.bytes((M + 1) // 8 + 1), "little") % (1 << (M + 1)) | 1
        ph = np.empty(M + 1, dtype=np.complex128); ph[:M] = phases_from_R(R, M); ph[M] = 1.0
        Pp = ph * prev[(lo - A) - prev_lo: (0 - prev_lo) + 1]
        cur = np.zeros(-lo + 1, dtype=np.complex128)
        for a in range(1, A + 1):
            cur += w[a - 1] * Pp[A - a: A - a + (-lo + 1)]
        prev, prev_lo = cur, lo
        v0[n] = abs(cur[-lo])
        Wn = min(W, -lo)
        e[n] = float(np.sum(np.abs(cur[-lo - Wn: -lo + 1]) ** 2))
    return e, v0

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    W = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"window energy (W={W}) versus typical rate, N={N}, A={A}; values rescaled by 3^(n/2) (energy by 3^n): a flat energy = exactly incoherent")
    def report(label, e, v):
        r = e[1:] / np.maximum(e[:-1], 1e-300)   # per-level energy ratio in rescaled units (incoherent = 1)
        a_, b_ = 100, N
        ns = np.arange(a_, b_ + 1)
        sl_e = np.polyfit(ns, np.log(np.maximum(e[a_:b_ + 1], 1e-300)), 1)[0]
        sl_v = np.polyfit(ns, np.log(np.maximum(v[a_:b_ + 1], 1e-300)), 1)[0]
        P(f"{label}: rescaled energy medians by 100-blocks " + " ".join(f"{float(np.median(e[a:a+100])):.3e}" for a in range(1, N, 100)) +
          f" | energy ratio mean {float(np.mean(r[a_:])):.4f} median {float(np.median(r[a_:])):.4f} (incoherent 1.0000) | mean-square amplitude rate {3**-0.5 * math.exp(sl_e / 2):.5f} | typical rate of |f_n(0)| {3**-0.5 * math.exp(sl_v):.5f}")
    t0 = time.time()
    for q in (3, 5, 7):
        e, v = run_energy(N, A, W, "real", q=q); report(f"real q={q}", e, v)
    for s in range(3):
        e, v = run_energy(N, A, W, "iid", seed=s); report(f"iid seed={s}", e, v)
    P(f"  [{time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
