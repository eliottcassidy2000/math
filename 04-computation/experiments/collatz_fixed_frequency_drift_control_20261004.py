#!/usr/bin/env python3
"""
Experiment E4 (opus, 2026-10-04): DRIFT control for the cold-frequency rate of E1.

For an odd multiplier q, the q-adic Syracuse law of qx+1 is Y_0 = 0, Y_n = 2^-a (q Y_{n-1} + 1) mod q^n (a geometric of
mean 2).  Its Fourier coefficient at a fixed unit u obeys the same closed window recursion as for q = 3:
    f_n(k) = sum_{a>=1} 2^-a e((u 2^{k-a} mod q^n)/q^n) f_{n-1}(k-a),  f_0 = 1,  mu_hat_n(u) = f_n(0).
The Parseval (random-phase) rate is q^(-1/2).  E1 found 0.568-0.571 < 3^(-1/2) = 0.5774 for q = 3 and six units.
Question: is the sub-Parseval decay a drift phenomenon (q = 3: 3 < 4, Tao/GGM's condition q < p^(p/(p-1))) or generic?
5x+1 (drift positive, divergent orbits expected) and 7x+1 are the controls.  We report the least-squares rates of
|mu_hat_n(1)| over n = 200..N for q in {3, 5, 7, 11}, against q^(-1/2), with the rescaled recursion v_n = q^(n/2) f_n.
Usage: python collatz_fixed_frequency_drift_control_20261004.py [N=600] [A=40]
"""
import sys, math, time
import numpy as np

def phases_u_q(u, q, n, kmin, kmax):
    mod = q ** n
    K = kmax - kmin + 1
    out = np.empty(K, dtype=np.complex128)
    if kmax >= 0:
        k0 = max(kmin, 0)
        r = (u * pow(2, k0, mod)) % mod
        vals = []
        for k in range(k0, kmax + 1):
            vals.append(r / mod); r = (2 * r) % mod
        out[k0 - kmin:] = np.exp(2j * np.pi * np.array(vals))
    if kmin < 0:
        M = -kmin
        R = (-u * pow(q, -n, 1 << (M + 1))) % (1 << (M + 1))
        bstr = bin(R)[2:].zfill(M + 1)[::-1]
        bits = (np.frombuffer(bstr.encode(), dtype=np.uint8) - 48).astype(np.float64)[:M]
        kern = 2.0 ** (-np.arange(1, 54))
        D = np.convolve(bits, kern)[:M]
        Mp = np.arange(1, M + 1, dtype=np.float64)
        corr = (u * np.exp2(-Mp) * (float(q) ** (-n))) if n * math.log(q) < 700 else np.zeros(M)
        ph = np.exp(2j * np.pi * (D + corr))
        kneg_max = min(kmax, -1)
        seg = ph[(-kneg_max) - 1: M][::-1]
        out[: kneg_max - kmin + 1] = seg
    return out

def run_uq(u, q, N, A):
    w = math.sqrt(q) * 2.0 ** (-np.arange(1, A + 1))
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

def direct_law_q(q, n, A):
    mod = q ** n
    law = {0: 1.0}
    inv2 = pow(2, -1, mod)
    for lvl in range(1, n + 1):
        new = {}
        for y, p in law.items():
            z = (q * y + 1) % mod
            r = z
            for a in range(1, A + 1):
                r = (r * inv2) % mod
                new[r] = new.get(r, 0.0) + p * 2.0 ** (-a)
        law = new
    return law

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    A = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"drift control: fixed-frequency rates for qx+1, N={N}, A={A}")
    P("control at n = 4, 5 (direct law, A = 30):")
    for q in (3, 5, 7):
        v = run_uq(1, q, 5, 30)
        for n in (4, 5):
            law = direct_law_q(q, n, 30); mod = q ** n
            val = abs(sum(p * np.exp(2j * np.pi * (y % mod) / mod) for y, p in law.items()))
            P(f"  q={q} n={n}: recursion {v[n] * float(q) ** (-n / 2):.10f}  direct {val:.10f}")
    for q in (3, 5, 7, 11):
        for u in (1, 5 if q != 5 else 2):
            t0 = time.time()
            v = run_uq(u, q, N, A)
            med = [float(np.median(v[a:a + 100])) for a in range(1, N, 100)]
            P(f"q={q:2d} u={u}: q^(n/2)|mu_hat_n(u)| 100-block medians: " + " ".join(f"{m:.2e}" for m in med) + f"  [{time.time()-t0:.1f}s]")
            for (a_, b_) in ((100, 300), (300, N), (200, N)):
                if b_ <= N and b_ > a_ + 20:
                    ns = np.arange(a_, b_ + 1); y = np.log(np.maximum(v[a_:b_ + 1], 1e-300)) - ns * math.log(q) / 2
                    sl = np.polyfit(ns, y, 1)[0]
                    P(f"     rate of |mu_hat_n({u})| over n={a_}..{b_}: {math.exp(sl):.4f}   (Parseval q^(-1/2) = {q**-0.5:.4f}; ratio {math.exp(sl)/q**-0.5:.4f})")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
