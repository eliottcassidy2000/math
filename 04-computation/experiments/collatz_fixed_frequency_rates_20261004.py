#!/usr/bin/env python3
"""
Experiment E1 (opus, 2026-10-04): is the H1 decay rate universal across FIXED frequencies?
(Prediction P6 of the coalescence note.)

The 3-adic Syracuse law mu_n (Y_0 = 0, Y_n = 2^-a (3 Y_{n-1} + 1) mod 3^n, a geometric of mean 2) has
    mu_hat_n(t) = E e(t Y_n / 3^n) = sum_{a>=1} 2^-a e(t 2^-a / 3^n) mu_hat_{n-1}(t 2^-a mod 3^{n-1}),
so for a fixed unit u the family f_n(k) = mu_hat_n(u 2^k), k <= 0, is closed:
    f_n(k) = sum_a 2^-a omega^{(u)}_n(k-a) f_{n-1}(k-a),   omega^{(u)}_n(j) = e((u 2^j mod 3^n)/3^n).
(For u = 1 this is the S22/S23 closed family on the powers of two; mac-mini's collatz_h1_20260929_mu1_rescaled.py.)
We run the rescaled window recursion v_n = 3^(n/2) f_n on [-(N-n)A, 0] for u in {1, 5, 7, 11, 13, 17} and report
3^(n/2) |mu_hat_n(u)| (= |v_n(0)|), block decay rates, and the ratio to u = 1.
Control: at n <= 7 the coefficients are compared with a direct enumeration of the law (valuations <= A).
Phases at negative exponents use the 2-adic digits of -u 3^-n (the same device as mac-mini's phases_level):
    (u 2^-M mod 3^n)/3^n  =  frac( R_u / 2^M ) + u 2^-M 3^-n,   R_u = -u 3^-n mod 2^(M+1).
Usage: python collatz_fixed_frequency_rates_20261004.py [N=800] [A=40]
"""
import sys, math, time
import numpy as np

def phases_u(u, n, kmin, kmax):
    """omega^{(u)}_n(k) for k in [kmin, kmax] (kmax <= 0), complex array."""
    mod = 3 ** n
    K = kmax - kmin + 1
    out = np.empty(K, dtype=np.complex128)
    # non-negative k (if any): direct
    if kmax >= 0:
        k0 = max(kmin, 0)
        r = (u * pow(2, k0, mod)) % mod
        vals = []
        for k in range(k0, kmax + 1):
            vals.append(r / mod); r = (2 * r) % mod
        out[k0 - kmin:] = np.exp(2j * np.pi * np.array(vals))
    if kmin < 0:
        M = -kmin                      # need k = -M .. -1 (and below if kmax < 0)
        Mtop = -max(kmin, kmin)        # all negative k in [kmin, min(kmax,-1)]
        R = (-u * pow(3, -n, 1 << (M + 1))) % (1 << (M + 1))
        bstr = bin(R)[2:].zfill(M + 1)[::-1]            # bstr[i] = bit b_(i+1) (2^i coefficient)
        bits = (np.frombuffer(bstr.encode(), dtype=np.uint8) - 48).astype(np.float64)[:M]
        kern = 2.0 ** (-np.arange(1, 54))               # frac(R/2^Mp) ~ sum_{t=0}^{52} b_(Mp-t) 2^-(t+1)
        D = np.convolve(bits, kern)[:M]                 # D[Mp-1] for Mp = 1..M
        Mp = np.arange(1, M + 1, dtype=np.float64)
        corr = (u * np.exp2(-Mp) * (3.0 ** (-n))) if n < 640 else np.zeros(M)
        ph = np.exp(2j * np.pi * (D + corr))             # index Mp-1 <-> k = -Mp
        kneg_max = min(kmax, -1)
        # fill k = kmin .. kneg_max  <-> Mp = M .. -kneg_max
        seg = ph[(-kneg_max) - 1: M][::-1]               # Mp from -kneg_max .. M, reversed to k ascending
        out[: kneg_max - kmin + 1] = seg
    return out

def run_u(u, N, A):
    w = math.sqrt(3.0) * 2.0 ** (-np.arange(1, A + 1))
    prev_lo = -N * A - A
    prev = np.ones(-prev_lo + 1, dtype=np.complex128)   # v_0 = 1 on [prev_lo, 0]
    v0 = np.zeros(N + 1); v0[0] = 1.0
    for n in range(1, N + 1):
        lo = -(N - n) * A
        ph = phases_u(u, n, lo - A, 0)                   # omega on [lo-A, 0]
        P = ph * prev[(lo - A) - prev_lo: (0 - prev_lo) + 1]
        cur = np.zeros(-lo + 1, dtype=np.complex128)
        for a in range(1, A + 1):
            cur += w[a - 1] * P[A - a: A - a + (-lo + 1)]
        prev, prev_lo = cur, lo
        v0[n] = abs(cur[-lo])
    return v0

def direct_law(n, A):
    """exact-ish law of Y_n mod 3^n with valuations <= A each (mass 1 - n 2^-A dropped), as a dict."""
    mod = 3 ** n
    law = {0: 1.0}
    inv2 = pow(2, -1, mod)
    for lvl in range(1, n + 1):
        new = {}
        for y, p in law.items():
            z = (3 * y + 1) % mod
            r = z
            for a in range(1, A + 1):
                r = (r * inv2) % mod
                new[r] = new.get(r, 0.0) + p * 2.0 ** (-a)
        law = new
    return law

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 800
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"fixed-frequency rates: N={N}, A={A}")
    # control at small n
    P("control: recursion vs direct enumeration of the law at n <= 7 (A = 30), |mu_hat_n(u)|:")
    for u in (1, 5, 7):
        v = run_u(u, 7, 30)
        for n in (3, 5, 7):
            law = direct_law(n, 30); mod = 3 ** n
            val = abs(sum(p * np.exp(2j * np.pi * (u * y % mod) / mod) for y, p in law.items()))
            P(f"  u={u:2d} n={n}: recursion {v[n] * 3.0 ** (-n / 2):.10f}  direct {val:.10f}")
    res = {}
    for u in (1, 5, 7, 11, 13, 17):
        t0 = time.time()
        v = run_u(u, N, A)
        res[u] = v
        med = [float(np.median(v[a:a + 200])) for a in range(1, N, 200)]
        P(f"u={u:2d}: 3^(n/2)|mu_hat_n(u)| 200-block medians: " + " ".join(f"{m:.3e}" for m in med) + f"   [{time.time()-t0:.1f}s]")
        for (a_, b_) in ((200, 400), (400, 600), (600, N), (200, N)):
            if b_ <= N and b_ > a_ + 20:
                ns = np.arange(a_, b_ + 1); y = np.log(v[a_:b_ + 1]) - ns * math.log(3) / 2
                sl = np.polyfit(ns, y, 1)[0]
                P(f"     least-squares rate of |mu_hat_n({u})| over n={a_}..{b_}: {math.exp(sl):.4f}   (Parseval 0.5774, critical 0.5850)")
        peaks = [(n, round(float(v[n]), 3)) for n in range(2, N) if v[n] >= 1.0 and v[n] >= v[n - 1] and v[n] >= v[n + 1]]
        P(f"     local maxima of 3^(n/2)|mu_hat_n({u})| >= 1: {peaks[:12]}")
    # cross-frequency comparison
    P("ratio medians |mu_hat_n(u)| / |mu_hat_n(1)| by 200-blocks:")
    for u in (5, 7, 11, 13, 17):
        r = res[u][1:] / np.maximum(res[1][1:], 1e-300)
        P(f"  u={u:2d}: " + " ".join(f"{float(np.median(r[a:a+200])):.3f}" for a in range(0, N, 200)))
    # truncation check at N/2 with A+20
    N2 = N // 2
    for u in (1, 5):
        va = run_u(u, N2, A); vb = run_u(u, N2, A + 20)
        rel = np.max(np.abs(va[1:] - vb[1:]) / np.maximum(np.abs(vb[1:]), 1e-300))
        P(f"truncation control u={u}: A={A} vs A={A+20} to n={N2}: max relative difference {rel:.2e}")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
