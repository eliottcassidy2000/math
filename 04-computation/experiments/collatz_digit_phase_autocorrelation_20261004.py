#!/usr/bin/env python3
"""
E6 (opus, 2026-10-04): autocorrelations of the window phases along the real digit string of -q^-n.

The window phases are theta_d = 0.b_d b_(d-1) ... b_1 (binary), b_i the 2-adic digits of -q^-n, i.e.
theta_d = frac(R / 2^d) with R = -q^-n mod 2^M -- the orbit of the inverse doubling map theta_d = (theta_(d-1) + b_d)/2.
For i.i.d. digits E e(theta_(d+s) - theta_d) = 0 for every lag s >= 1 (the top bit b_(d+s) enters with coefficient 1/2).
Here we measure, for the REAL digit strings (q = 3, 5, 7, 11; n = 60; M = 2^20 digits) and for i.i.d. strings,
    C_s = (1/M) sum_d e(theta_(d+s) - theta_d),  s = 1..40,
their moduli against the i.i.d. noise level M^(-1/2), and the digit balance and the lag-s bit correlations
(1/M) sum_d (-1)^(b_d + b_(d+s)).  Note: the digits of q^-n mod 2^M are those of q^(2^(M-2) - n) mod 2^M, a power of
q, so this is a statement about the low binary digits of powers of q (Mahler's (3/2)^n territory; Dupuy-Weirich's
averaged equidistribution).  A systematic nonzero C_s for q = 3 would be an arithmetic correlation that the window
recursion could turn into the 0.5% extra cancellation of E5b.
Usage: python collatz_digit_phase_autocorrelation_20261004.py [n=60] [logM=20]
"""
import sys, math, time
import numpy as np

def digit_string(q, n, M, mode="real", seed=0):
    if mode == "real":
        R = (-pow(q, -n, 1 << (M + 1))) % (1 << (M + 1))
    else:
        rng = np.random.default_rng(seed)
        R = int.from_bytes(rng.bytes((M + 1) // 8 + 1), "little") % (1 << (M + 1)) | 1
    bstr = bin(R)[2:].zfill(M + 1)[::-1]
    bits = (np.frombuffer(bstr.encode(), dtype=np.uint8) - 48).astype(np.float64)[:M]
    return bits

def phases(bits):
    kern = 2.0 ** (-np.arange(1, 54))
    D = np.convolve(bits, kern)[:len(bits)]      # theta_d for d = 1..M (index d-1)
    return np.exp(2j * np.pi * D)

if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    logM = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    M = 1 << logM
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"phase autocorrelations along the digit string of -q^-n, n={n}, M=2^{logM} digits; i.i.d. noise level M^(-1/2) = {M**-0.5:.2e}")
    lags = list(range(1, 41))
    t0 = time.time()
    for label, q, mode, seed in (("real q=3", 3, "real", 0), ("real q=5", 5, "real", 0), ("real q=7", 7, "real", 0),
                                 ("real q=11", 11, "real", 0), ("iid seed=0", 0, "iid", 0), ("iid seed=1", 0, "iid", 1)):
        bits = digit_string(q, n, M, mode, seed)
        ph = phases(bits)
        bal = float(np.mean(bits))
        C = []
        B = []
        for s in lags:
            C.append(np.mean(ph[s:] * np.conj(ph[:-s])))
            B.append(float(np.mean((1 - 2 * bits[s:]) * (1 - 2 * bits[:-s]))))
        Cm = np.abs(np.array(C))
        P(f"{label:11s}: digit balance {bal:.5f}; |C_s| for s=1..12: " + " ".join(f"{c:.1e}" for c in Cm[:12]) +
          f"; max |C_s| s<=40: {Cm.max():.1e} at s={lags[int(np.argmax(Cm))]}; rms |C_s| {float(np.sqrt(np.mean(Cm**2))):.1e}; bit correlations s=1..8: " +
          " ".join(f"{b:+.1e}" for b in B[:8]))
        # also the mean phase itself (should be ~0) and the pair sums with the 2^-d-2^-e weights of the recursion
        P(f"             mean phase |E e(theta)| = {abs(np.mean(ph)):.1e}; |E e(2 theta)| = {abs(np.mean(ph**2)):.1e}; |E e(theta_d - 2 theta_(d-1))| (the odometer relation, exact 1 if b_d=0 only) = {abs(np.mean(ph[1:] * np.conj(ph[:-1]**2))):.3f}")
    P(f"  [{time.time()-t0:.0f}s]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
