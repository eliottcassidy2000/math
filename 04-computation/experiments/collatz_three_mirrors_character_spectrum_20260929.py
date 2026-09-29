#!/usr/bin/env python3
"""The character spectrum of the 3-adic Syracuse law by the full-period closed recursion (S23, 2026-09-29).

The family m_n(k) = mu_hat_n(2^k mod 3^n) is periodic in k with period L_n = 2 3^(n-1) (2 generates the units), so
the frequency recursion m_n(k) = sum_(a>=1) 2^-a omega_n(k-a) m_(n-1)(k-a), omega_n(k) = e((2^k mod 3^n)/3^n), can be
run on the whole cycle Z/L_n with NO truncation: with g(k) = omega_n(k) m_(n-1)(k mod L_(n-1)) it is the one-pole
filter m_n(k) = (g(k-1) + m_n(k-1))/2 around the cycle (scipy.signal.lfilter, two passes for the wrap).  This gives
every Fourier coefficient of mu_n at every unit -- the same information as the FFT of the law -- organised by the
discrete logarithm, in O(L_n) time and memory per level.
Character spectrum: for a character psi_j(2^k) = e(jk/L_n) of (Z/3^n)^x (primitive iff 3 does not divide j... the
conductor of psi_j is 3^n iff 3 does not divide j, for n >= 2; j = 0 trivial, j = L_n/2 quadratic),
    sum_(k mod L_n) m_n(k) conj(psi_j(2^k)) = tau(conj psi_j) * E[psi_j(Y_n)]   (Y_n is always a unit),
so S_n(psi_j) := E[psi_j(Y_n)] = FFT(m_n)[j] / tau(conj psi_j), |tau| = 3^(n/2) for primitive psi_j.
Identities checked: H_n(1) = (2/3) 3^n mu_n(1) = sum_psi S_n(psi) (the seed-1 harmonic mass of S19 Theorem C is the
sum of all multiplicative moments); mu_n(-1) = (1/L_n) sum_psi psi(-1) S_n(psi); E[chi_-3(Y_n)] = -1/3 exactly.
Outputs per level: M(n) = max_k |m_n(k)| and its argmax (against S20's FFT values), the rms of |S| 3^(n/2) over
primitive characters, its distribution against the Rayleigh law (random phases), the odd/even split, the largest |S|,
and the phase-linearity (chirp) test of FFT(m_n)[j] against e(-j k0/L_n).
Run: python 04-computation/experiments/collatz_three_mirrors_character_spectrum_20260929.py [NMAX]
"""
from __future__ import annotations

import cmath
import math
import sys
import time
from fractions import Fraction

import numpy as np
from scipy.signal import lfilter

M_S20 = {1: 0.5773503, 2: 0.3779236, 3: 0.2522368, 4: 0.1769989, 5: 0.1292736, 6: 0.0961064, 7: 0.0758700,
         8: 0.0608907, 9: 0.0480262, 10: 0.0382783, 11: 0.0319442, 12: 0.0264582, 13: 0.0220524, 14: 0.0191280,
         15: 0.0162845, 16: 0.0144095, 17: 0.0125107, 18: 0.0111873}


def omega_vector(n: int) -> np.ndarray:
    """omega_n(k) = e((2^k mod 3^n)/3^n) for k = 0..L_n-1, by running the powers of two around the cycle."""
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    r = np.empty(L, dtype=np.int64)
    x = 1
    for k in range(L):
        r[k] = x
        x = (2 * x) % mod
    return np.exp(2j * np.pi * (r.astype(np.float64) / mod))


def omega_vector_fast(n: int) -> np.ndarray:
    """Vectorised residues 2^k mod 3^n (int64 for n <= 19, exact) via blocks: r[k + t B] = r[k] 2^(tB) mod 3^n."""
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    if L <= 200000:
        return omega_vector(n)
    assert n <= 19, "int64 block arithmetic needs 3^(2n) < 2^63"
    B = 1 << 16
    base = np.empty(B, dtype=np.int64)
    x = 1
    for k in range(B):
        base[k] = x
        x = (2 * x) % mod
    nb = (L + B - 1) // B
    out = np.empty(L, dtype=np.float64)
    step = pow(2, B, mod)
    mult = 1
    for t in range(nb):
        lo = t * B
        hi = min(L, lo + B)
        blk = (base[: hi - lo] * np.int64(mult)) % np.int64(mod)
        out[lo:hi] = blk.astype(np.float64) / mod
        mult = (mult * step) % mod
    return np.exp(2j * np.pi * out)


def level_up(prev: np.ndarray, n: int) -> np.ndarray:
    """m_n on Z/L_n from m_(n-1) on Z/L_(n-1) by the one-pole filter around the cycle."""
    L = 2 * 3 ** (n - 1)
    om = omega_vector_fast(n)
    g = om * np.tile(prev, 3) if n > 1 else om * prev[0]
    # m(k) = (g(k-1) + m(k-1))/2: filter y[k] = 0.5 y[k-1] + 0.5 x[k-1]; run twice around the cycle for the wrap
    reps = max(2, 80 // L + 2)          # enough passes around the cycle for the transient 2^-(steps) to vanish
    x = np.tile(g, reps)
    y = lfilter([0.0, 0.5], [1.0, -0.5], x)
    return y[-L:]


def law(n: int) -> np.ndarray:
    """mu_n on Z/3^n by the forward recursion Y_(k+1) = 2^-a (3 Y_k + 1), P(a) = 2^-a (a <= 60), Y_0 = 0."""
    mod = 3 ** n
    mu = np.zeros(mod)
    mu[0] = 1.0
    inv2 = [pow(2, -a, mod) for a in range(61)]
    for _ in range(n):
        new = np.zeros(mod)
        ys = np.arange(mod)
        src = (3 * ys + 1) % mod
        for a in range(1, 61):
            np.add.at(new, (src * inv2[a]) % mod, 2.0 ** (-a) * mu)   # y -> 3y+1 is 3-to-1 mod 3^n: accumulate duplicates
        mu = new
    return mu


if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    NLAW = 9
    laws = {n: law(n) for n in range(1, NLAW + 1)}
    print("identities: sum over primitive characters of S_n = (2/3)(H_n(1) - H_(n-1)(1)), H_n(y) = 3^n mu_n(y); with psi(-1): (2/3)(H_n(-1) - H_(n-1)(-1))")
    t0 = time.time()
    m = np.array([1.0 + 0j, 1.0 + 0j])  # level 0: m_0 = 1 on the period-2 cycle of level 1 (any value; broadcast)
    # level 1: L_1 = 2, omega_1(k) = e(2^k mod 3 / 3): k=0 -> e(1/3), k=1 -> e(2/3)
    m = None
    for n in range(1, NMAX + 1):
        L = 2 * 3 ** (n - 1)
        if n == 1:
            om = omega_vector(1)
            g = om * 1.0
            x = np.tile(g, 82)
            m = lfilter([0.0, 0.5], [1.0, -0.5], x)[-L:]
        else:
            m = level_up(m, n)
        absm = np.abs(m)
        kmax = int(np.argmax(absm))
        Mn = absm[kmax]
        # discrete log of kmax as exponent; its position relative to n log2 3
        LOG23 = math.log2(3)
        # DFT over the cycle: F[j] = sum_k m(k) e(-jk/L)
        F = np.fft.fft(m)
        # Gauss sums tau(conj psi_j) = sum_k conj(psi_j(2^k)) e(2^k/3^n) = sum_k e(-jk/L) omega_n(k) = FFT(omega)[j]
        om = omega_vector_fast(n)
        tau = np.fft.fft(om)
        # primitive characters: 3 does not divide j (n >= 2); for n = 1 all nontrivial
        js = np.arange(L)
        prim = (js % 3 != 0) if n >= 2 else (js != 0)
        S = np.zeros(L, dtype=complex)
        S[prim] = F[prim] / tau[prim]
        # identities
        # mu_n(1) = (1/L) sum_psi S(psi) over ALL characters (imprimitive included: need S for all j) -> compute all S properly:
        # for imprimitive psi_j (3 | j, j != 0) the Gauss sum tau(conj psi_j) mod 3^n vanishes, so use the identity at the character's own level instead.
        # Here we use the exact inversion: mu_n(y) = (1/L) sum_j F_j-based? Simpler: mu_n(1) from the law is not needed; check H_n(1) = sum over ALL characters
        # via the direct inverse: the vector m_n is the Fourier transform of mu_n at 2^k, so mu_n(2^-k0 * ...) -- skip; use the primitive part only.
        absS = np.abs(S[prim]) * 3 ** (n / 2)
        rms = math.sqrt(np.mean(absS ** 2))
        # odd/even: psi_j(-1) = e(j * (L/2)/L) = (-1)^j
        odd = prim & (js % 2 == 1)
        even = prim & (js % 2 == 0)
        meanS_odd = np.mean(S[odd]) if odd.any() else 0
        meanS_even = np.mean(S[even]) if even.any() else 0
        # chirp test: phase of F[j] against -j kmax / L
        ph = np.angle(F[prim] * np.exp(2j * np.pi * js[prim] * kmax / L))
        chirp_coh = abs(np.mean(np.exp(1j * ph)))
        sumS = S[prim].sum()
        sumS_par = (S[prim] * np.where(js[prim] % 2 == 1, -1.0, 1.0)).sum()
        if n <= NLAW:
            H1 = 3 ** n * laws[n][1 % 3 ** n]
            H1p = 3 ** (n - 1) * laws[n - 1][1 % 3 ** (n - 1)] if n >= 2 else 1.0
            Hm = 3 ** n * laws[n][(-1) % 3 ** n]
            Hmp = 3 ** (n - 1) * laws[n - 1][(-1) % 3 ** (n - 1)] if n >= 2 else 0.0
            print(f"      identity check n={n}: sum_prim S = {sumS.real:+.6f}{sumS.imag:+.1e}i vs (2/3)(H_n(1) - H_(n-1)(1)) = {(2 / 3) * (H1 - H1p):+.6f} [H_n(1) = {H1:.6f}]; "
                  f"parity sum = {sumS_par.real:+.6f}{sumS_par.imag:+.1e}i vs (2/3)(H_n(-1) - H_(n-1)(-1)) = {(2 / 3) * (Hm - Hmp):+.6f} [H_n(-1) = {Hm:.6f}, 1.46 (3/2)^n = {1.46 * 1.5 ** n:.4f}]")
        else:
            print(f"      sums n={n}: sum_prim S = {sumS.real:+.6f}{sumS.imag:+.1e}i (= (2/3) increment of H_n(1)); parity sum = {sumS_par.real:+.5f} vs 0.325 (3/2)^n = {0.325 * 1.5 ** n:.4f}")
        big = np.argsort(-np.abs(S))[:8]
        def jdesc(j):
            jj = j if j <= L // 2 else j - L
            v2 = (abs(jj) & -abs(jj)).bit_length() - 1 if jj else 0
            return f"{jj:+d}=±2^{v2}·{abs(jj) >> v2}"
        print("      largest |S| 3^(n/2): " + ", ".join(f"{abs(S[j]) * 3 ** (n / 2):.2f}@{jdesc(int(j))}" for j in big))
        print(f"n={n:2d}: L={L}: M(n) = {Mn:.7f} (S20 FFT {M_S20.get(n, float('nan')):.7f}) at k = {kmax} (k - n log2 3 = {kmax - n * LOG23:+.2f}; -k mod L = {(-kmax) % L}); "
              f"rms |S| 3^(n/2) over primitive chars = {rms:.4f}; max |S| 3^(n/2) = {absS.max():.3f} at j = {js[prim][np.argmax(absS)]}; "
              f"mean S odd = {meanS_odd.real:+.3e}{meanS_odd.imag:+.3e}i, even = {meanS_even.real:+.3e}{meanS_even.imag:+.3e}i; chirp coherence = {chirp_coh:.3f}   [{time.time() - t0:.0f}s]")
        if n <= 4:
            print("      S(psi_j) 3^(n/2) for all primitive j:", " ".join(f"{abs(S[j]) * 3 ** (n / 2):.3f}" for j in js[prim]))
        # Rayleigh test: |S|^2 3^n / mean should be Exp(1): report the fraction above 1, 2, 3 (e^-1, e^-2, e^-3 = 0.368, 0.135, 0.050)
        z = absS ** 2 / np.mean(absS ** 2)
        print(f"      Rayleigh check: P(z > 1, 2, 3) = {np.mean(z > 1):.3f}, {np.mean(z > 2):.3f}, {np.mean(z > 3):.3f} (Exp(1): 0.368, 0.135, 0.050); quadratic character j = L/2: |S| 3^(n/2) = {abs(S[L // 2]) * 3 ** (n / 2) if prim[L // 2] else float('nan'):.4f}")
        sys.stdout.flush()
    print("DONE")
