#!/usr/bin/env python3
"""
Lyapunov exponent of the step-1 Pascal tower (opus, 2026-10-05).

OBJECT (Fourier note collatz_fourier_experiments_20261004.md, sections 4j-4m): the uniform-start q-tower.
Level-n window recursion
    f_n(k) = sum_{c>=1} 2^-c e(theta_{n,k-c}) f_{n-1}(k-c),   f_0 = 1,   theta_{n,-d} = (x_n mod 2^d)/2^d,
with x_n = q^{N-n} x_N and x_N a Haar-random 2-adic unit R.  Expanding the recursion, f_N(0) is the path sum
    f_N(0) = E_c e( R Phi_N ),   Phi_N = sum_{i=1}^N q^{i-1} 2^{-S_i},   S_i = c_1 + ... + c_i,  c_i iid geometric(1/2),
i.e. the Fourier transform, at the random integer frequency R, of the law of the N-fold random backward Syracuse
iterate y -> (1 + q y)/2^c of 0 (real-valued; the i-th one of a fair coin sequence sits at S_i).

EXACT REFORMULATION USED HERE (twisted Pascal triangle; no window truncation).  Let
    H_m(j) := 2^-m sum_{coin words of length m with j ones} prod_{i<=j} e( R q^{i-1} / 2^{S_i} ).
Pascal's rule with a twist on the up-step:
    H_m(j) = (1/2) [ H_{m-1}(j) + omega_{m,j} H_{m-1}(j-1) ],   omega_{m,j} = e( (R q^{j-1} mod 2^m) / 2^m ),
and
    f_j(0) = hat nu_j(R) = sum_{m>=j} (1/2) omega_{m,j} H_{m-1}(j-1)        (the total flux into column j).
Column j is computed from column j-1 by a one-dimensional recursion in the depth m with a per-column scale
(no over/underflow), to depth M_max = 5N: the neglected flux is at most P(S_j > 5j) ~ e^{-0.964 j} against
|f_j(0)| ~ e^{-0.5625 j}.  Cost 5 N^2 twists per sample (the window recursion costs N^2 A^2 / 2 with A = 40).
One sample gives log|f_j(0)| for EVERY j <= N at the same R (Haar invariance makes each level's law the
uniform-start law of the Fourier note).

TWIST MODELS.  tower q: x_j = R q^j (q = 3 is the step-1 Pascal tower, q = 5 step-2, ...);  iid: independent Haar
columns (the i.i.d.-digit model);  randmult: x_j = R Q^j with Q a random odd 40-bit multiplier;
collatz: the exact coefficient mu_hat_N(u) of the 3-adic Syracuse law, phases (u 2^-m mod 3^n)/3^n at level n,
written 2-adically as frac((-u 3^-n mod 2^m)/2^m) + u 3^-n 2^-m (identity checked exactly in --selftest).

Usage:
  python collatz_pascal_tower_lyapunov_20261005.py --selftest
  python collatz_pascal_tower_lyapunov_20261005.py --model tower --q 3 --N 1500 --seeds 48 --workers 10
  python collatz_pascal_tower_lyapunov_20261005.py --model iid --N 1500 --seeds 48
  python collatz_pascal_tower_lyapunov_20261005.py --model collatz --u 1 --N 131
"""
import sys, os, math, time, argparse
from fractions import Fraction
import numpy as np
import numba as nb


@nb.njit(cache=True)
def columns_lyap(bits, corr, Mmax, N, Lout, pup=0.5):
    """bits: (N, W) uint64, little-endian packed binary digits of the 2-adic unit of column j+1 (row j).
    corr: (N,) float64 additive phase corrections corr[j] * 2^-m (zero for the random models).
    pup: probability of a one (an odd step) per coin; the valuations are then geometric(pup) with weights
    pup (1-pup)^(c-1) (pup = 1/2 is the Collatz/Syracuse weight 2^-c).
    Fills Lout[j-1] = log|f_j(0)| for j = 1..N."""
    P = np.zeros(Mmax + 1, dtype=np.complex128)
    C = np.zeros(Mmax + 1, dtype=np.complex128)
    qdn = 1.0 - pup
    v = 1.0
    for m in range(Mmax + 1):          # column 0: H_m(0) = (1-p)^m
        P[m] = v
        v *= qdn
    s_prev = 0.0
    twopi = 2.0 * np.pi
    for j in range(1, N + 1):
        row = bits[j - 1]
        cj = corr[j - 1]
        theta = 0.0
        pw = 1.0
        flux = 0.0 + 0.0j
        C[j - 1] = 0.0
        for m in range(1, Mmax + 1):
            b = (row[(m - 1) >> 6] >> ((m - 1) & 63)) & 1
            theta = 0.5 * (theta + b)              # (x mod 2^m)/2^m
            pw *= 0.5
            if m >= j:
                ph = theta + cj * pw
                om = np.exp(1j * twopi * ph)
                t = om * P[m - 1]
                flux += pup * t
                C[m] = qdn * C[m - 1] + pup * t
        mx = 0.0
        for m in range(j, Mmax + 1):
            a = abs(C[m])
            if a > mx:
                mx = a
        Lout[j - 1] = s_prev + math.log(abs(flux))
        inv = 1.0 / mx
        for m in range(j, Mmax + 1):
            C[m] *= inv
        s_prev += math.log(mx)
        P, C = C, P
    return 0


@nb.njit(cache=True)
def window_lyap(bits, corr, A, N, Lout, pup=0.5):
    """The level (window) recursion of the Fourier note with valuations truncated to c <= A:
    f_n(k) = sum_{c=1}^A w_c e(theta_{n,c-k}) f_{n-1}(k-c), k in [-(N-n)A, 0], w_c = pup (1-pup)^(c-1),
    theta_{n,d} = (x_n mod 2^d)/2^d from the bits of row N-n (column n <-> row N-n, as in make_bits).
    Fills Lout[n-1] = log|f_n(0)|.  Independent of columns_lyap (no path reindexing)."""
    Wlen = N * A + A + 1
    prev = np.zeros(Wlen, dtype=np.complex128)     # prev[d] = f_{n-1}(-d)
    cur = np.zeros(Wlen, dtype=np.complex128)
    for d in range(Wlen):
        prev[d] = 1.0
    w = np.zeros(A + 1)
    for c in range(1, A + 1):
        w[c] = pup * (1.0 - pup) ** (c - 1)
    scale = 0.0
    twopi = 2.0 * np.pi
    for n in range(1, N + 1):
        row = bits[N - n]
        cj = corr[N - n]
        D = (N - n) * A + A           # phases needed to depth D (k - c >= -(N-n)A - A)
        ph = np.zeros(D + 1, dtype=np.complex128)
        theta = 0.0
        pw = 1.0
        ph[0] = 1.0
        for d in range(1, D + 1):
            b = (row[(d - 1) >> 6] >> ((d - 1) & 63)) & 1
            theta = 0.5 * (theta + b)
            pw *= 0.5
            ph[d] = np.exp(1j * twopi * (theta + cj * pw))
        Kmax = (N - n) * A
        mx = 0.0
        for k in range(0, Kmax + 1):       # k = -depth of the output index
            acc = 0.0 + 0.0j
            for c in range(1, A + 1):
                acc += w[c] * ph[k + c] * prev[k + c]
            cur[k] = acc
            a = abs(acc)
            if a > mx:
                mx = a
        Lout[n - 1] = scale + math.log(abs(cur[0]))
        inv = 1.0 / mx
        for k in range(0, Kmax + 1):
            cur[k] *= inv
        scale += math.log(mx)
        prev, cur = cur, prev
    return 0


def pack(x, W):
    return np.frombuffer(x.to_bytes(W * 8, "little"), dtype=np.uint64).copy()


def make_bits(model, N, Mmax, q=3, u=1, seed=0, Q=None):
    W = (Mmax + 63) // 64
    mod = 1 << Mmax
    rng = np.random.default_rng(seed)
    bits = np.zeros((N, W), dtype=np.uint64)
    corr = np.zeros(N, dtype=np.float64)
    if model in ("tower", "randmult"):
        R = int.from_bytes(rng.bytes(Mmax // 8 + 2), "little") % mod | 1
        mult = q if model == "tower" else (Q if Q is not None else (int(rng.integers(1, 1 << 40)) | 1))
        x = R
        for j in range(N):
            bits[j] = pack(x % mod, W)
            x = (x * mult) % mod
    elif model == "iid":
        for j in range(N):
            x = int.from_bytes(rng.bytes(Mmax // 8 + 2), "little") % mod | 1
            bits[j] = pack(x, W)
    elif model == "collatz":
        # column j+1 <-> level n = N - j; 2-adic unit -u 3^-n, correction u 3^-n
        for j in range(N):
            n = N - j
            un = u % 3 ** n if n < 64 else u          # the phase depends on u mod 3^n only
            x = (-un * pow(3, -n, mod)) % mod
            bits[j] = pack(x, W)
            corr[j] = un * 3.0 ** (-n) if n < 640 else 0.0
    else:
        raise ValueError(model)
    return bits, corr


def one_sample(args):
    model, N, q, u, seed, Mfac = args[:6]
    pup = args[6] if len(args) > 6 else 0.5
    A = args[7] if len(args) > 7 else 0
    L = np.zeros(N, dtype=np.float64)
    if A > 0:
        Mmax = N * A + A + 1
        bits, corr = make_bits(model, N, Mmax, q=q, u=u, seed=seed)
        window_lyap(bits, corr, A, N, L, pup)
    else:
        Mmax = Mfac * N
        bits, corr = make_bits(model, N, Mmax, q=q, u=u, seed=seed)
        columns_lyap(bits, corr, Mmax, N, L, pup)
    return L


def selftest():
    print("selftest 1: 2-adic form of the exact Collatz phase (u 2^-d mod 3^n)/3^n = frac((-u 3^-n mod 2^d)/2^d) + u 3^-n 2^-d")
    bad = 0
    for u in (1, 5, 7, 11, 13):
        for n in range(1, 9):
            for d in range(1, 30):
                un = u % 3 ** n
                lhs = Fraction((u * pow(2, -d, 3 ** n)) % 3 ** n, 3 ** n)
                t = (-un * pow(3, -n, 2 ** d)) % 2 ** d
                rhs = Fraction(t, 2 ** d) + Fraction(un, 2 ** d * 3 ** n)
                if lhs != rhs:
                    bad += 1
    print("   mismatches:", bad)
    print("selftest 2: triangle flux vs brute-force path sum (N = 6, valuations <= 14, random R), tower q = 3 and q = 5")
    rng = np.random.default_rng(1)
    from itertools import product
    for q in (3, 5):
        N = 6; cap = 14; Mmax = N * cap + 2
        R = int(rng.integers(1, 1 << 60)) | 1
        tot = 0.0 + 0.0j
        for cs in product(range(1, cap + 1), repeat=N):
            S = 0; ph = 0.0
            for i, c in enumerate(cs):
                S += c
                ph += ((R * q ** i) % (1 << S)) / (1 << S)
            tot += 2.0 ** (-sum(cs)) * np.exp(2j * np.pi * ph)
        W = (Mmax + 63) // 64
        bits = np.zeros((N, W), dtype=np.uint64); corr = np.zeros(N)
        x = R
        for j in range(N):
            bits[j] = pack(x % (1 << Mmax), W); x = x * q
        L = np.zeros(N); columns_lyap(bits, corr, Mmax, N, L)
        print(f"   q={q}: brute |f_6| = {abs(tot):.12f}   triangle = {math.exp(L[-1]):.12f}   (cap error ~ {N*2.0**-cap:.1e})")
    print("selftest 3: exact Collatz coefficients against the Fourier note E1 (3^(n/2)|mu_hat_n(1)|): n=5 1.341, 17 2.043, 127 3.882, 129 7.221, 131 10.377, 135 2.624; |mu_hat_7(7)| = 0.0129696631")
    for (u, N, ref) in ((1, 5, 1.341), (1, 17, 2.043), (1, 127, 3.882), (1, 129, 7.221), (1, 131, 10.377), (1, 135, 2.624)):
        L = one_sample(("collatz", N, 3, u, 0, 6))
        print(f"   u={u} n={N}: 3^(n/2)|mu_hat| = {math.exp(L[-1] + N * math.log(3) / 2):.4f}   (E1: {ref})")
    L = one_sample(("collatz", 7, 3, 7, 0, 6))
    print(f"   u=7 n=7: |mu_hat_7(7)| = {math.exp(L[-1]):.10f}   (E1: 0.0129696631)")
    print("selftest 4: tower q=3, one seed, N = 300: same R for all levels; rates")
    L = one_sample(("tower", 300, 3, 1, 7, 5))
    j = np.arange(1, 301)
    sl = np.polyfit(j[60:], L[60:], 1)[0]
    print(f"   slope 60..300: {sl:.4f}  rate {math.exp(sl):.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--model", default="tower")
    ap.add_argument("--q", type=int, default=3)
    ap.add_argument("--u", type=int, default=1)
    ap.add_argument("--N", type=int, default=1500)
    ap.add_argument("--seeds", type=int, default=24)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--Mfac", type=int, default=5)
    ap.add_argument("--save", default="")
    ap.add_argument("--p", type=float, default=0.5, help="probability of an odd step per coin (valuation law geometric(p))")
    ap.add_argument("--A", type=int, default=0, help="if > 0: window recursion with valuations truncated to c <= A")
    ap.add_argument("--umax", type=int, default=0, help="collatz mode: run every odd unit u <= umax prime to 3 (exact mu_hat_n(u), all n <= N, needs --A)")
    a = ap.parse_args()
    if a.selftest:
        selftest(); return
    t0 = time.time()
    if a.model == "collatz" and a.umax > 0:
        jobs = [(a.model, a.N, a.q, u, 0, a.Mfac, a.p, a.A) for u in range(1, a.umax + 1, 2) if u % 3 != 0]
    else:
        jobs = [(a.model, a.N, a.q, a.u, a.seed0 + s, a.Mfac, a.p, a.A) for s in range(a.seeds)]
    if a.workers > 1 and a.seeds > 1:
        import multiprocessing as mp
        with mp.Pool(a.workers) as pool:
            Ls = pool.map(one_sample, jobs)
    else:
        Ls = [one_sample(j) for j in jobs]
    Ls = np.array(Ls)                       # (seeds, N)
    N = a.N
    tag = f"{a.model}" + (f" q={a.q}" if a.model in ("tower",) else "") + (f" u={a.u}" if a.model == "collatz" and a.umax == 0 else "") + (f" units u<={a.umax} prime to 3 ({len(jobs)} units)" if a.umax > 0 else "") + (f" p={a.p}" if a.p != 0.5 else "") + (f" A={a.A}" if a.A > 0 else "")
    # reference mean-square rate for these weights: sqrt(sum_c w_c^2)
    wc = [a.p * (1 - a.p) ** (c - 1) for c in range(1, (a.A if a.A > 0 else 400) + 1)]
    rms = math.sqrt(sum(x * x for x in wc)); wsum = sum(wc)
    if a.save:
        np.save(a.save, Ls)
    N0 = N // 5
    lam_end = (Ls[:, N - 1] - Ls[:, N0 - 1]) / (N - N0)
    j = np.arange(1, N + 1)
    slopes = np.array([np.polyfit(j[N0 - 1:], L[N0 - 1:], 1)[0] for L in Ls])
    m = lam_end.mean(); se = lam_end.std(ddof=1) / math.sqrt(len(lam_end)) if len(lam_end) > 1 else 0.0
    ms = slopes.mean(); ses = slopes.std(ddof=1) / math.sqrt(len(slopes)) if len(slopes) > 1 else 0.0
    print(f"{tag}: N={N} Mmax={a.Mfac}N seeds={len(jobs)} [{time.time()-t0:.0f}s]")
    print(f"  lambda (L_N - L_N0)/(N-N0), N0={N0}: {m:.5f} +- {se:.5f}   rate {math.exp(m):.5f} +- {math.exp(m)*se:.5f}")
    print(f"  lambda least-squares slope {N0}..{N}: {ms:.5f} +- {ses:.5f}   rate {math.exp(ms):.5f} +- {math.exp(ms)*ses:.5f}")
    print(f"  per-seed sd of (L_N-L_N0)/(N-N0): {lam_end.std(ddof=1):.5f};  mean-square rate reference sqrt(sum w_c^2) = {rms:.5f} (sum w_c = {wsum:.5f}); Jensen gap log(rate/rms) = {m - math.log(rms):.5f}")
    for D in (1, 10, 100, min(1000, N - N0 - 1)):
        if D <= 0:
            continue
        inc = Ls[:, N0 - 1 + D:N] - Ls[:, N0 - 1:N - D]
        print(f"  D={D}: mean increment/D {inc.mean()/D:.5f}  Var/D {inc.var()/D:.5f}")
    print(f"  mean log|f_N| / N = {Ls[:, N-1].mean()/N:.5f}   (includes the transient)")


if __name__ == "__main__":
    main()
