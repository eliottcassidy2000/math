#!/usr/bin/env python3
"""
Experiment E5 (opus, 2026-10-04): the random-digit control for the cold-frequency rate.

The level-n phases of the fixed-frequency recursion are the binary digits of -u q^-n:
    omega_n(-a) = e( 0.b_a b_{a-1} ... b_1 (binary) + u 2^-a q^-n ),  b_i = i-th 2-adic digit of -u q^-n,
so mu_hat_n(u) = f_n(0) is a functional of the 2-adic digit strings of q^-1, ..., q^-n (consecutive strings are
related by the 2-adic multiplication by q^-1, i.e. the carry automaton of x -> q x run backwards).
E1/E4: the decay rate is ~0.57 for every q and every fixed unit, i.e. the incoherent rate sqrt(sum_a 4^-a) = 1/sqrt3 =
0.57735 of the recursion, with a measured deficit of ~1.3% for q = 3 at n <= 2500 (mac-mini) and an unresolved one
for other q at n <= 600.
Control: replace the digit string of -q^-n at every level by (i) i.i.d. random bits (fresh per level), and
(ii) the digit string of a random odd 2-adic unit raised to the n-th power (a random q!), and measure the rate.
If (i) gives 1/sqrt3 within noise while q = 3 gives 0.569, the deficit is arithmetic; if (i) also sits ~0.57, it is
a property of the recursion itself.
Usage: python collatz_fixed_frequency_random_digits_20261004.py [N=1500] [A=40] [seeds=3]
"""
import sys, math, time
import numpy as np

def phases_from_R(R, M, extra=0.0):
    """phases for k = -M..-1 from the integer R (digits b_1..b_(M+1) low to high): e(frac(R/2^Mp))."""
    bstr = bin(R)[2:].zfill(M + 1)[::-1]
    bits = (np.frombuffer(bstr.encode(), dtype=np.uint8) - 48).astype(np.float64)[:M]
    kern = 2.0 ** (-np.arange(1, 54))
    D = np.convolve(bits, kern)[:M]
    ph = np.exp(2j * np.pi * D)
    return ph[::-1]   # k ascending from -M to -1

def run(N, A, mode, q=3, seed=0):
    rng = np.random.default_rng(seed)
    w = math.sqrt(3.0) * 2.0 ** (-np.arange(1, A + 1))   # rescale by the INCOHERENT rate sqrt3 (not q^(1/2))
    prev_lo = -N * A - A
    prev = np.ones(-prev_lo + 1, dtype=np.complex128)
    v0 = np.zeros(N + 1); v0[0] = 1.0
    if mode == "randomunit":
        qq = int(rng.integers(1, 1 << 40)) | 1   # a random odd 'multiplier'
    if mode == "uniformstart":
        # the UNIFORM-START q-tower of the Fourier note, section 4k: the deepest level N has a uniform odd digit
        # string R0 (to the full window depth), and level n < N is R0 * q^(N-n) (the q-th-root / Pascal law).
        Mmax = N * A + A
        R0 = int.from_bytes(rng.bytes((Mmax + 1) // 8 + 1), "little") % (1 << (Mmax + 1)) | 1
    for n in range(1, N + 1):
        lo = -(N - n) * A
        M = -(lo - A)
        if mode == "real":
            R = (-pow(q, -n, 1 << (M + 1))) % (1 << (M + 1))
        elif mode == "iid":
            R = int.from_bytes(rng.bytes((M + 1) // 8 + 1), "little") % (1 << (M + 1)) | 1
        elif mode == "randomunit":
            R = (-pow(qq, -n, 1 << (M + 1))) % (1 << (M + 1))
        elif mode == "uniformstart":
            R = (R0 * pow(q, N - n, 1 << (M + 1))) % (1 << (M + 1))
        ph_neg = phases_from_R(R, M)          # k = -M .. -1
        ph = np.empty(M + 1, dtype=np.complex128); ph[:M] = ph_neg; ph[M] = 1.0   # k = 0: e(u/q^n) ~ 1
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
    S = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    P(f"random-digit control, N={N}, A={A}; values are 3^(n/2)|f_n(0)| (rescaled by the incoherent rate)")
    def report(label, v):
        meds = [float(np.median(v[a:a + 300])) for a in range(1, N, 300)]
        line = f"{label}: 300-block medians " + " ".join(f"{m:.2e}" for m in meds)
        for (a_, b_) in ((300, N), (N // 2, N)):
            ns = np.arange(a_, b_ + 1); y = np.log(np.maximum(v[a_:b_ + 1], 1e-300)) - ns * math.log(3) / 2
            sl = np.polyfit(ns, y, 1)[0]
            line += f" | rate {a_}..{b_}: {math.exp(sl):.4f} (ratio to 1/sqrt3: {math.exp(sl)*math.sqrt(3):.4f})"
        P(line)
    t0 = time.time()
    qs = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [3, 5]
    for q in qs:
        v = run(N, A, "real", q=q); report(f"real q={q}", v)
    P(f"  [{time.time()-t0:.0f}s]")
    for s in range(S):
        v = run(N, A, "iid", seed=s); report(f"iid digits seed={s}", v)
    for s in range(S):
        v = run(N, A, "randomunit", seed=100 + s); report(f"random odd multiplier seed={s}", v)
    P(f"  [{time.time()-t0:.0f}s total]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("DONE")
