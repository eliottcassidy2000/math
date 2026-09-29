#!/usr/bin/env python3
"""Independent audit of collatz_three_mirrors_20260929.md (S23), sections 0-5 and 7.  Own code throughout; the repo's
scripts are not imported or rerun.  Written 2026-09-29 by the auditor session.

What is recomputed
  A. the exact law mu_n of Tao's Syracuse variable on Z/3^n, n <= NLAW (default 16), by a forward DP that embeds
     level n-1 into level n (y -> 3y+1 is injective from Z/3^(n-1) into Z/3^n; the mixture over the valuation is a
     gather, no scatter), truncation a <= 60;
  B. Theorem 1: for every level, (i) F_n(j) = sum_k m_n(k) e(-jk/L) against tau(conj psi_j) S_n(psi_j) with
     S_n(psi_j) = E[psi_j(Y_n)] computed DIRECTLY from the law in discrete-log coordinates, and the vanishing of
     both sides at the imprimitive characters; (ii) sum_prim S_n and the parity sum against L mu_n(+-1) -
     (L/3) mu_(n-1)(+-1) = rho_n(+-1) - rho_(n-1)(+-1), the n = 1 case, the accumulated sums against S19's
     independent values (mazur_harmonic_mass_deep18 output), and the spectrum statistics of the note;
  C. Proposition 2: the one-pole recursion implemented three ways -- a plain sequential filter in Python (n <= 10),
     an exact circular-convolution closed form by FFT (n <= NLAW), and a sequential C program from level 0 to
     level 19 (streaming the last level) -- each compared elementwise with the m_n obtained from the law, and the
     maxima against the S20 FFT values quoted in the note;
  D. the 3-adic zero-line census: box A (u <= 1000, Q <= 2000) in Python integers, box B (u <= 5000, Q <= 10^4) in
     C with residues mod 3^30, the Poisson tails, the deepest lines, the S22 seeds, an LTE check, and a variance
     test showing where the geometric law is forced by equidistribution rather than "Poisson-consistent";
  E. the constants theta*, I and the identity e^(-I) 3^(theta*/ln 2) = log_2 3 - 1;
  F. the saddle chain, factorisations, merges, the -17 cycle and its clock;
  G. the basins by an independent bitmask sieve in C on [1, 2^30] (records EVERY target on the orbit, not the
     first, so nestedness is tested rather than assumed) cross-checked by a pure-Python sieve on [1, 2^22], and the
     Krasikov-Lagarias arithmetic.
Run: python3 04-computation/experiments/collatz_three_mirrors_20260929_audit.py [NLAW=16] [NREC=19] [NBASIN=30]
"""
from __future__ import annotations

import cmath
import math
import os
import subprocess
import sys
import tempfile
import time
from fractions import Fraction

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
T0 = time.time()
NLAW = int(sys.argv[1]) if len(sys.argv) > 1 else 16
NREC = int(sys.argv[2]) if len(sys.argv) > 2 else 19
NBASIN = int(sys.argv[3]) if len(sys.argv) > 3 else 30
SCRATCH = os.environ.get("AUDIT_SCRATCH") or os.path.join(tempfile.gettempdir(), "three_mirrors_audit_tmp")   # never inside the repo
os.makedirs(SCRATCH, exist_ok=True)
LOG23 = math.log2(3)


def stamp() -> str:
    return f"[{time.time() - T0:6.0f}s]"


def banner(s: str) -> None:
    print()
    print("=" * 100)
    print(s)
    print("=" * 100)
    sys.stdout.flush()


# ----------------------------------------------------------------------------------------------------------------------
# E. constants (done first: trivial)
# ----------------------------------------------------------------------------------------------------------------------
banner("E. constants of S22 and the identity e^(-I) 3^(theta*/ln 2) = log_2 3 - 1")
m_ = math.log2(3)
p_ = math.log(2, 3)
q_ = 1 - p_
theta = math.log(2 * q_)
I_ = theta * m_ - math.log(m_ - 1)
print(f"m = log_2 3 = {m_:.6f}, p = log_3 2 = {p_:.6f}, q = 1 - p = {q_:.6f}, theta* = ln(2q) = {theta:.6f}, I = theta* m - ln(m-1) = {I_:.6f}")
print(f"e^(-I) = {math.exp(-I_):.6f} (S22: 3^(h*-1) = 0.946505); theta*/ln 2 = {theta / math.log(2):.5f} (note: -0.438)")
lhs = math.exp(-I_) * 3 ** (theta / math.log(2))
print(f"e^(-I) 3^(theta*/ln 2) = {lhs:.15f}; log_2 3 - 1 = {m_ - 1:.15f}; difference {lhs - (m_ - 1):.1e}")
print("algebra: 3^(theta*/ln 2) = e^(theta* ln 3/ln 2) = e^(theta* m); e^(-I) = e^(-theta* m) (m-1); product = m - 1.  The identity is the")
print("         definition of I rewritten (I := theta* m - ln(m-1)); it holds for ANY theta*, not only the tilt -- it carries no information about the law.")
print(f"3^(-1/2) = {3 ** -0.5:.6f}; (m-1) = {m_ - 1:.6f}; ratio 3^(-1/2)/(m-1) = {3 ** -0.5 / (m_ - 1):.5f} -> margin {100 * (1 - 3 ** -0.5 / (m_ - 1)):.2f}% (note: 1.3%)")
print(f"Chernoff extrapolation at u = 3^n: e^(-nI) 3^(n theta*/ln 2) = (m-1)^n = {m_ - 1:.4f}^n; Parseval rms scale 3^(-n/2) = {3 ** -0.5:.4f}^n")
print(f"3^n (m-1)^(2n) = ({3 * (m_ - 1) ** 2:.4f})^n: a family of ~3^n units each of size (m-1)^n would carry a Fourier mass growing like {3 * (m_ - 1) ** 2:.4f}^n -- the note's 'Parseval forbids' reading (heuristic: the families u 2^j overlap)")
sys.stdout.flush()

# ----------------------------------------------------------------------------------------------------------------------
# A. the law by forward DP (level-embedding form)
# ----------------------------------------------------------------------------------------------------------------------
banner(f"A. the law mu_n on Z/3^n by forward DP, n <= {NLAW} (embedding y -> 3y+1 of level n-1 into level n; gather over a <= 60)")


def law_levels(nmax: int, amax: int = 60) -> dict[int, np.ndarray]:
    laws = {}
    mu = np.array([1.0])                     # level 0: point mass on Z/1
    for n in range(1, nmax + 1):
        mod, modp = 3 ** n, 3 ** (n - 1)
        w = np.zeros(mod)
        ys = np.arange(modp, dtype=np.int64)
        w[(3 * ys + 1) % mod] = mu           # injective: no duplicates to accumulate at this step
        new = np.zeros(mod)
        idx = np.arange(mod, dtype=np.int64)
        for a in range(1, amax + 1):
            idx = (idx * 2) % mod            # 2^a z mod 3^n
            new += (2.0 ** -a) * w[idx]      # P(Y_n = z) = sum_a 2^-a P(3 Y_(n-1) + 1 = 2^a z)
        mu = new
        laws[n] = mu
        print(f"   level {n:2d}: 3^n = {mod:9d}, total mass {mu.sum():.15f}, mass on non-units {mu[::3].sum():.1e}, mu_n(1) = {mu[1 % mod]:.10f}, mu_n(-1) = {mu[(-1) % mod]:.10f}   {stamp()}")
        sys.stdout.flush()
    return laws


laws = law_levels(NLAW)
# exact fractions at level 2 (the note: H_2(1) = 8/7, H_2(-1) = 22/7)
mu2_1 = Fraction(1, 3) * Fraction(1, 4) * Fraction(64, 63) + Fraction(2, 3) * Fraction(1, 16) * Fraction(64, 63)
print(f"   exact level 2: mu_2(1) = {mu2_1} -> H_2(1) = 9 mu_2(1) = {9 * mu2_1} (note: 8/7); DP gives H_2(1) = {9 * laws[2][1]:.12f}, H_2(-1) = {9 * laws[2][8]:.12f} (note: 22/7 = {22 / 7:.12f})")
# consistency check: projection of mu_n to Z/3^(n-1) equals mu_(n-1)
for n in range(2, NLAW + 1):
    modp = 3 ** (n - 1)
    proj = np.bincount(np.arange(3 ** n) % modp, weights=laws[n], minlength=modp)
    print(f"   consistency ||proj mu_{n} - mu_{n - 1}||_1 = {np.abs(proj - laws[n - 1]).sum():.1e}")

# ----------------------------------------------------------------------------------------------------------------------
# helpers for the unit cycle
# ----------------------------------------------------------------------------------------------------------------------


def pow2_residues(n: int) -> np.ndarray:
    """r[k] = 2^k mod 3^n for k = 0..L_n - 1 (int64, exact for n <= 19), by a base block and block multipliers."""
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    B = min(L, 1 << 16)
    base = np.empty(B, dtype=np.int64)
    x = 1
    for k in range(B):
        base[k] = x
        x = (2 * x) % mod
    out = np.empty(L, dtype=np.int64)
    step = pow(2, B, mod)
    mult = 1
    lo = 0
    while lo < L:
        hi = min(L, lo + B)
        out[lo:hi] = (base[: hi - lo] * np.int64(mult)) % np.int64(mod)
        mult = (mult * step) % mod
        lo = hi
    assert out[0] == 1 and (out[L - 1] * 2) % mod == 1, "2 must generate the units with period L_n"
    return out


def signed_index(j: int, L: int) -> int:
    return j if j <= L // 2 else j - L


def jdesc(j: int, L: int) -> str:
    jj = signed_index(j, L)
    if jj == 0:
        return "0"
    v2 = (abs(jj) & -abs(jj)).bit_length() - 1
    return f"{jj:+d}" + (f"=±2^{v2}·{abs(jj) >> v2}" if v2 else "")


# ----------------------------------------------------------------------------------------------------------------------
# B. Theorem 1 and the spectrum, level by level, from the law
# ----------------------------------------------------------------------------------------------------------------------
banner("B. Theorem 1 (character-spectrum duality, the two increment identities) and the spectrum statistics, from the law")
S19_RHO1 = {1: 0.66667, 2: 0.76190, 3: 0.66138, 4: 0.61842, 5: 0.64287, 6: 0.63695, 7: 0.57357, 8: 0.51638, 9: 0.46450,
            10: 0.42497, 11: 0.39428, 12: 0.35824, 13: 0.33343, 14: 0.31549, 15: 0.30591, 16: 0.29011, 17: 0.28053, 18: 0.27364}
S19_RHOM = {1: 1.333, 2: 2.095, 3: 3.198, 4: 4.846, 5: 7.319, 6: 11.02, 7: 16.57, 8: 24.89, 9: 37.38, 10: 56.12, 11: 84.23,
            12: 126.4, 13: 189.6, 14: 284.5, 15: 426.8, 16: 640.2, 17: 960.4, 18: 1441.0}
NOTE_RMS = {2: 0.8452, 3: 0.8321, 4: 0.8345, 5: 0.8356, 6: 0.8362, 7: 0.8356, 8: 0.8360, 9: 0.8367, 10: 0.8373, 11: 0.8380,
            12: 0.8387, 13: 0.8393, 14: 0.8399, 15: 0.8404, 16: 0.8408}
NOTE_MAX = {2: 1.00, 3: 1.41, 4: 1.59, 5: 2.02, 6: 2.21, 7: 2.82, 8: 2.95, 9: 3.60, 10: 4.09, 11: 4.35, 12: 5.75, 13: 6.04,
            14: 7.41, 15: 7.79, 16: 9.11}
NOTE_ARGMAX = {3: 2, 4: 2, 5: 8, 6: 2, 7: 80, 8: 278, 9: 521, 10: 1193, 11: 1127, 12: 4219, 13: 40162, 14: 71333, 15: 154112, 16: 1516613}
NOTE_T = {1: -0.5, 2: 0.143, 3: -0.151, 4: -0.064, 5: 0.037, 6: -0.009, 7: -0.095, 8: -0.086, 9: -0.078, 10: -0.059, 11: -0.046,
          12: -0.054, 13: -0.037, 14: -0.027, 15: -0.014, 16: -0.024}

m_law = {}          # m_n(k) = mu_hat_n(2^k) from the law
sumS_prim = {}
sumS_par = {}
acc1 = 2.0 / 3.0    # rho_1(1)
accm = 4.0 / 3.0    # rho_1(-1)
for n in range(1, NLAW + 1):
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    mu = laws[n]
    r = pow2_residues(n)
    x = mu[r]                                    # law on the units in discrete-log order
    S = L * np.fft.ifft(x)                       # S[j] = sum_k x[k] e(+jk/L) = E[psi_j(Y_n)]
    muhat = np.fft.fft(mu)                       # muhat[t] = sum_y mu(y) e(-ty/3^n)
    m = muhat[(-r) % mod]                        # m_n(k) = mu_hat_n(2^k) with the e(+) convention
    del muhat
    m_law[n] = m
    F = np.fft.fft(m)                            # F[j] = sum_k m(k) e(-jk/L)
    omega = np.exp(2j * np.pi * (r.astype(np.float64) / mod))
    tau = np.fft.fft(omega)                      # tau[j] = sum_k conj(psi_j)(2^k) e(2^k/3^n) = tau(conj psi_j)
    del omega
    js = np.arange(L)
    prim = (js % 3 != 0) if n >= 2 else (js != 0)
    imp = ~prim
    scale = 3 ** (n / 2)
    err_i = np.max(np.abs(F[prim] - tau[prim] * S[prim])) / scale
    tau_mod = np.abs(tau[prim]) / scale
    print(f"n={n:2d} L={L:9d}: (i) max|F - tau S|/3^(n/2) over primitive j = {err_i:.1e}; |tau|/3^(n/2) in [{tau_mod.min():.12f}, {tau_mod.max():.12f}]; "
          f"imprimitive j: max|F| = {np.max(np.abs(F[imp])) if imp.any() else 0:.1e}, max|tau| = {np.max(np.abs(tau[imp])) if imp.any() else 0:.1e}, "
          f"max|S| = {np.max(np.abs(S[imp])) if imp.any() else 0:.3f} (S itself does NOT vanish there; it is the lower-level moment; the vanishing of F and tau is for n >= 2)")
    if n >= 2:
        jj = js[imp]
        dev = np.max(np.abs(S[imp] - S_prev[(jj // 3) % len(S_prev)]))
        print(f"      imprimitive consistency: max_j |S_n(psi_3j') - S_(n-1)(psi_j')| = {dev:.1e} (psi_j primitive iff 3 does not divide j: the induced characters give the lower-level moments)")
    S_prev = S.copy()
    # (ii)
    sS = S[prim].sum()
    sP = (S[prim] * np.where(js[prim] % 2 == 1, -1.0, 1.0)).sum()
    sumS_prim[n], sumS_par[n] = sS.real, sP.real
    mu1, mum = mu[1 % mod], mu[(-1) % mod]
    if n >= 2:
        modp = 3 ** (n - 1)
        mu1p = mu[np.arange(mod) % modp == 1 % modp].sum()          # mu_n(y = 1 mod 3^(n-1)) = mu_(n-1)(1) by consistency
        mump = mu[np.arange(mod) % modp == (-1) % modp].sum()
        rhs1 = L * mu1 - (L / 3) * mu1p
        rhsm = L * mum - (L / 3) * mump
        rho_inc1 = (2 / 3) * (3 ** n * mu1 - 3 ** (n - 1) * laws[n - 1][1 % modp])
        rho_incm = (2 / 3) * (3 ** n * mum - 3 ** (n - 1) * laws[n - 1][(-1) % modp])
        print(f"      (ii) sum_prim S = {sS.real:+.9f}{sS.imag:+.1e}i vs L mu_n(1) - (L/3) mu_(n-1)(1) = {rhs1:+.9f} = rho_n(1) - rho_(n-1)(1) = {rho_inc1:+.9f}; "
              f"parity sum = {sP.real:+.9f}{sP.imag:+.1e}i vs {rhsm:+.9f} = {rho_incm:+.9f}")
        acc1 += sS.real
        accm += sP.real
    else:
        print(f"      (ii) n = 1: the only primitive character is chi_-3: S_1(chi_-3) = {sS.real:+.6f} (note: -1/3), parity sum {sP.real:+.6f} (note: +1/3); "
              f"rho_1(1) = {(2 / 3) * 3 * mu1:.6f} = 1 - 1/3, rho_1(-1) = {(2 / 3) * 3 * mum:.6f} = 1 + 1/3; the general formula would give rho_1 - rho_0 = {(2 / 3) * 3 * mu1 - 2 / 3:+.3f}, so n = 1 is indeed a separate case")
    rho1 = (2 / 3) * 3 ** n * mu1
    rhom = (2 / 3) * 3 ** n * mum
    print(f"      accumulated: rho_n(1) = 2/3 + sum_(m=2..n) sum_prim S_m = {acc1:.6f} vs law {rho1:.6f} vs S19 {S19_RHO1.get(n, float('nan')):.5f}; "
          f"rho_n(-1) = 4/3 + sum parity = {accm:.4f} vs law {rhom:.4f} vs S19 {S19_RHOM.get(n, float('nan')):.4g}; "
          f"H_n(1) = (3/2) rho_n(1) = {1.5 * rho1:.6f}; T_n = (3/2) sum_prim S_n = {1.5 * sS.real:+.4f} (note {NOTE_T.get(n, float('nan')):+.3f}); "
          f"0.325 (3/2)^n = {0.325 * 1.5 ** n:.3f}")
    # spectrum statistics over the primitive characters
    absS = np.abs(S[prim]) * scale
    rms = math.sqrt(np.mean(absS ** 2))
    mass_units = float(np.sum(np.abs(m) ** 2))                 # sum over the units of |mu_hat_n|^2 (S20's 'Fourier mass at level n')
    z = absS ** 2 / np.mean(absS ** 2)
    jp = js[prim]
    order = np.argsort(-absS)
    top = order[:8]
    nprim = int(prim.sum())
    ray_all = rms * math.sqrt(max(math.log(L), 0.0))
    ray_prim = rms * math.sqrt(max(math.log(nprim), 0.0))
    ray_pairs = rms * math.sqrt(max(math.log(nprim / 2) + 0.5772, 0.0))
    odd = jp % 2 == 1
    mo = S[prim][odd].mean() if odd.any() else 0j
    me = S[prim][~odd].mean() if (~odd).any() else 0j
    print(f"      spectrum: rms |S| 3^(n/2) = {rms:.4f} (note {NOTE_RMS.get(n, float('nan')):.4f}); rms^2 = {rms ** 2:.4f} = (3/2) sum_units |mu_hat|^2 = {1.5 * mass_units:.4f} "
          f"[sum_units |mu_hat_n|^2 = {mass_units:.4f} = S20's 'mass at level n'; the fullperiod script's 'mass per level' 3^n sum|m|^2/L = {3 ** n * mass_units / L:.4f}]")
    print(f"      P(z > 1, 2, 3) = {np.mean(z > 1):.3f}, {np.mean(z > 2):.3f}, {np.mean(z > 3):.3f} (Exp(1): 0.368, 0.135, 0.050); "
          f"max |S| 3^(n/2) = {absS.max():.3f} (note {NOTE_MAX.get(n, float('nan')):.2f}) at j = {jdesc(int(jp[order[0]]), L)} (note ±{NOTE_ARGMAX.get(n, '-')}); "
          f"top 8: " + ", ".join(f"{absS[i]:.2f}@{jdesc(int(jp[i]), L)}" for i in top))
    print(f"      Rayleigh-max reference: rms sqrt(ln L) = {ray_all:.2f}; rms sqrt(ln #prim) = {ray_prim:.2f} (#prim = {nprim}); rms sqrt(ln(#prim/2) + gamma) = {ray_pairs:.2f} (conjugate pairs |S_j| = |S_-j|); "
          f"max/rms = {absS.max() / rms:.2f}; n/2 = {n / 2:.1f}")
    print(f"      parity: mean S odd = {mo.real:+.4e}{mo.imag:+.1e}i, even = {me.real:+.4e}{me.imag:+.1e}i, sum of means = {(mo + me).real:+.2e}; #odd = {int(odd.sum())}, #even = {int((~odd).sum())}; "
          f"predicted mean_even ~ (parity sum)/(2 #even) = {sP.real / (2 * max(int((~odd).sum()), 1)):+.4e}")
    del F, tau, S, x
    sys.stdout.flush()

print()
print("Summary of Theorem 1(ii) against S19 (deep18 output):")
acc = 2 / 3
accm_ = 4 / 3
for n in range(2, NLAW + 1):
    acc += sumS_prim[n]
    accm_ += sumS_par[n]
    flag1 = abs(acc - S19_RHO1[n]) < 6e-6
    print(f"   n={n:2d}: rho_n(1) accumulated {acc:.6f} vs S19 {S19_RHO1[n]:.5f} [{'5 digits' if flag1 else 'MISMATCH'}]; parity increment {sumS_par[n]:.5f} vs S19 difference "
          f"{S19_RHOM[n] - S19_RHOM[n - 1]:.4g} (S19 prints rho_n(-1) to 4 significant digits: differences known to ~3-4 digits only)")
print("   the note's 'partial sums 1 + sum_(m<=n) T_m' with T_1 = -1/2: at n = 2 this is", 1 - 0.5 + 1.5 * sumS_prim[2], "but H_2(1) = 8/7 =", 8 / 7,
      "; the correct partial sum is 1 + sum_(2<=m<=n) T_m (or 3/2 + sum_(1<=m<=n) T_m):", 1 + 1.5 * sumS_prim[2])
n16 = min(NLAW, 16)
Lm = 2 * 3 ** (n16 - 1)
print(f"   random-phase size of sum_prim S_m at m = {n16}: sqrt(#prim) rms 3^(-m/2) = sqrt((2/3) L_m) 0.84 3^(-m/2) = {math.sqrt(2 / 3 * Lm) * 0.84 * 3 ** (-n16 / 2):.3f} "
      f"(with L_m in place of the primitive count (2/3) L_m one gets {math.sqrt(Lm) * 0.84 * 3 ** (-n16 / 2):.3f}, the note's 0.69); observed |sum| = {abs(sumS_prim[n16]):.4f}")
sys.stdout.flush()

# ----------------------------------------------------------------------------------------------------------------------
# C. Proposition 2: the one-pole recursion, three implementations
# ----------------------------------------------------------------------------------------------------------------------
banner("C. Proposition 2 (full-period one-pole recursion): sequential Python filter (n <= 10), FFT closed form (n <= NLAW), sequential C to level NREC")
S20_M = {1: 0.5773503, 2: 0.3779236, 3: 0.2522368, 4: 0.1769989, 5: 0.1292736, 6: 0.0961064, 7: 0.0758700, 8: 0.0608907,
         9: 0.0480262, 10: 0.0382783, 11: 0.0319442, 12: 0.0264582, 13: 0.0220524, 14: 0.0191280, 15: 0.0162845,
         16: 0.0144095, 17: 0.0125107, 18: 0.0111873}
S20_ARGS = {1: 0, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 8, 8: 9, 9: 10, 10: 12, 11: 13, 12: 14, 13: 16, 14: 17, 15: 18, 16: 20, 17: 21, 18: 23}
print("derivation: mu_hat_n(t) = sum_a 2^-a e((t 2^-a mod 3^n)/3^n) mu_hat_(n-1)(t 2^-a mod 3^(n-1)) [condition on the last valuation a;")
print("  Y_n = 2^-a (3 Y_(n-1) + 1); the inverse of 2 mod 3^n reduces to the inverse mod 3^(n-1)].  With t = 2^k: m_n(k) = sum_(a>=1) 2^-a g(k-a),")
print("  g(k) = omega_n(k) m_(n-1)(k mod L_(n-1)).  Then m_n(k-1) = sum_(a>=1) 2^-a g(k-1-a) = 2 [m_n(k) - g(k-1)/2], i.e. m_n(k) = (g(k-1) + m_n(k-1))/2.")
print("  The homogeneous solutions c 2^k are not L-periodic (2^L != 1), so the periodic solution is unique; from any state the error is 2^-(steps).")


def filter_sequential(prev: np.ndarray, n: int, warm: int = 100) -> np.ndarray:
    """Plain sequential one-pole filter around the cycle Z/L_n from m_(n-1) (Python loop)."""
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    Lp = len(prev)
    out = np.empty(L, dtype=np.complex128)
    k = (L - 1 - warm) % L                  # index of the first g used; after `warm` steps we are at k = L-1 feeding m(0)
    rr = pow(2, k, mod)
    mstate = 0j
    for step in range(warm + L):
        g = cmath.exp(2j * math.pi * rr / mod) * complex(prev[k % Lp])
        mstate = 0.5 * (g + mstate)
        if step >= warm:
            out[(k + 1) % L] = mstate
        k = (k + 1) % L
        rr = (2 * rr) % mod
    return out


def filter_fft(prev: np.ndarray, n: int) -> np.ndarray:
    """Exact closed form: m = circular convolution of g with h_L(a) = 2^-a/(1 - 2^-L) (a = 1..L); transfer (w/2)/(1 - w/2), w = e(-j/L)."""
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    r = pow2_residues(n)
    g = np.exp(2j * np.pi * (r.astype(np.float64) / mod)) * np.tile(prev, L // len(prev))
    w = np.exp(-2j * np.pi * np.arange(L) / L)
    H = (w / 2) / (1 - w / 2)
    return np.fft.ifft(np.fft.fft(g) * H)


m0 = np.array([1.0 + 0j])
prev = m0
for n in range(1, NLAW + 1):
    mfft = filter_fft(prev if n > 1 else m0, n)      # from the LAW's m_(n-1) (one step, no accumulation)
    d_fft = np.max(np.abs(mfft - m_law[n]))
    line = f"n={n:2d}: one recursion step from the law's m_(n-1): FFT closed form max|m_rec - m_law| = {d_fft:.1e}"
    if n <= 10:
        mseq = filter_sequential(prev if n > 1 else m0, n)
        line += f"; sequential filter max|m_rec - m_law| = {np.max(np.abs(mseq - m_law[n])):.1e}"
    absm = np.abs(m_law[n])
    kmax = int(np.argmax(absm))
    Mn = absm[kmax]
    Ln = 2 * 3 ** (n - 1)
    kk = signed_index(kmax, Ln)
    mirror = signed_index((kmax + Ln // 2) % Ln, Ln)
    s20 = S20_ARGS[n]
    match = (s20 % Ln == kmax) or (s20 % Ln == (kmax + Ln // 2) % Ln)
    s_pow = kk if abs(kk - n * LOG23) < abs(mirror - n * LOG23) else mirror
    line += (f"; M(n) = max_k |m_n(k)| = {Mn:.7f} (S20 FFT {S20_M[n]:.7f}, diff {Mn - S20_M[n]:+.1e}) at k = {kmax} (signed {kk}, mirror {mirror}; S20 argmax 2^{s20} is one of ±2^k: {match}); "
             f"nearest-to-n log2 3 exponent s = {s_pow}: s - n log2 3 = {s_pow - n * LOG23:+.2f}; floor(n log2 3) - 6 = {math.floor(n * LOG23) - 6}")
    print(line)
    prev = m_law[n]
    sys.stdout.flush()

# the C program: sequential recursion from level 0, streaming the last level, dumping m_12 and m_16 for elementwise comparison
FP_C = r'''
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
/* usage: fp NMAX DUMPDIR : full-period one-pole recursion m_n(k) = (g(k-1) + m_n(k-1))/2, g(k) = e(2^k mod 3^n / 3^n) m_(n-1)(k mod L_(n-1)),
   from m_0 = 1, warm-up 100 steps then one pass around Z/L_n; the last level is streamed (not stored). */
typedef double complex cd;
static void report(int n, uint64_t L, uint64_t mod, double sumsq, uint64_t *topk, double *topv, int nt, double t) {
  double mass = pow(3.0, n) * sumsq / (double)L;
  printf("n=%2d: L=%llu: M(n) = %.7f at k = %llu; Fourier mass 3^n sum|m|^2/L = %.4f (sum over units of |mu_hat|^2 = %.4f); top-5:", n, (unsigned long long)L, topv[0], (unsigned long long)topk[0], mass, sumsq);
  for (int i = 0; i < nt; i++) {
    long long k = (long long)topk[i]; long long kk = k <= (long long)(L / 2) ? k : k - (long long)L;
    long long mir = (long long)((topk[i] + L / 2) % L); long long mm = mir <= (long long)(L / 2) ? mir : mir - (long long)L;
    printf(" k=%lld(|m|=%.7f; signed %lld, mirror %lld)", k, topv[i], kk, mm);
  }
  printf("   [%.0fs]\n", t); fflush(stdout);
}
int main(int argc, char **argv) {
  int NMAX = atoi(argv[1]); const char *dumpdir = argc > 2 ? argv[2] : NULL;
  cd *prev = malloc(sizeof(cd)); prev[0] = 1.0; uint64_t Lp = 1;
  clock_t c0 = clock();
  for (int n = 1; n <= NMAX; n++) {
    uint64_t mod = 1; for (int i = 0; i < n; i++) mod *= 3;
    uint64_t L = 2 * (mod / 3);
    int store = (n < NMAX);
    cd *cur = store ? malloc(sizeof(cd) * L) : NULL;
    if (store && !cur) { fprintf(stderr, "alloc fail at n=%d\n", n); return 1; }
    const double invmod = 1.0 / (double)mod;
    int warm = 100;
    uint64_t k = (L - 1 - (warm % L) + L) % L;           /* index of the first g used */
    /* r = 2^k mod 3^n */
    uint64_t r = 1; { uint64_t b = 2 % mod, e = k; while (e) { if (e & 1) r = (r * b) % mod; b = (b * b) % mod; e >>= 1; } }
    uint64_t idx = k % Lp;
    cd m = 0;
    double sumsq = 0; uint64_t topk[5] = {0}; double topv[5] = {0};
    for (uint64_t step = 0; step < (uint64_t)warm + L; step++) {
      double ph = 2 * M_PI * (double)r * invmod;
      cd g = (cos(ph) + I * sin(ph)) * prev[idx];
      m = 0.5 * (g + m);
      if (step >= (uint64_t)warm) {
        uint64_t kout = (k + 1) % L;
        if (store) cur[kout] = m;
        double a = cabs(m); sumsq += a * a;
        if (a > topv[4]) { int i = 4; while (i > 0 && topv[i - 1] < a) { topv[i] = topv[i - 1]; topk[i] = topk[i - 1]; i--; } topv[i] = a; topk[i] = kout; }
      }
      k = (k + 1 == L) ? 0 : k + 1;
      r = (2 * r) % mod;
      idx = (idx + 1 == Lp) ? 0 : idx + 1;
    }
    report(n, L, mod, sumsq, topk, topv, 5, (double)(clock() - c0) / CLOCKS_PER_SEC);
    if (store && dumpdir && (n == 12 || n == 16)) {
      char path[1024]; snprintf(path, sizeof path, "%s/m_%d.bin", dumpdir, n);
      FILE *f = fopen(path, "wb"); if (f) { fwrite(cur, sizeof(cd), L, f); fclose(f); }
    }
    free(prev); prev = cur; Lp = L;
    if (!store) break;
  }
  return 0;
}
'''
src = os.path.join(SCRATCH, "fp_audit.c")
exe = os.path.join(SCRATCH, "fp_audit.exe")
with open(src, "w") as f:
    f.write(FP_C)
cc = subprocess.run(["gcc", "-O2", "-o", exe, src, "-lm"], capture_output=True, text=True)
print("gcc (fp_audit):", "ok" if cc.returncode == 0 else cc.stderr)
sys.stdout.flush()
if cc.returncode == 0:
    proc = subprocess.run([exe, str(NREC), SCRATCH.replace("\\", "/")], capture_output=True, text=True)
    print(proc.stdout)
    if proc.stderr:
        print("stderr:", proc.stderr)
    # parse M(n) and compare with S20 and with the law-derived maxima
    for line in proc.stdout.splitlines():
        if line.startswith("n=") and "M(n) = " in line:
            n = int(line[2:4])
            Mn = float(line.split("M(n) = ")[1].split(" ")[0])
            kmax = int(line.split("at k = ")[1].split(";")[0])
            L = 2 * 3 ** (n - 1)
            kk = signed_index(kmax, L)
            mirror = signed_index((kmax + L // 2) % L, L)
            s_pow = kk if abs(kk - n * LOG23) < abs(mirror - n * LOG23) else mirror
            ref = S20_M.get(n)
            extra = f"S20 {ref:.7f} diff {Mn - ref:+.1e}" if ref else "(beyond S20's FFT range)"
            s20 = S20_ARGS.get(n)
            match = "" if s20 is None else f"; S20 argmax 2^{s20} is one of ±2^k: {(s20 % L == kmax) or (s20 % L == (kmax + L // 2) % L)}"
            print(f"   C recursion n={n:2d}: M(n) = {Mn:.7f} {extra}; argmax k = {kmax} (signed {kk}, mirror {mirror}){match}; s = {s_pow}: s - n log2 3 = {s_pow - n * LOG23:+.2f}; floor(n log2 3) - 6 = {math.floor(n * LOG23) - 6}"
                  + (f"; law max {np.max(np.abs(m_law[n])):.7f}" if n in m_law else ""))
    for n in (12, 16):
        path = os.path.join(SCRATCH, f"m_{n}.bin")
        if n in m_law and os.path.exists(path):
            mc = np.fromfile(path, dtype=np.complex128)
            print(f"   elementwise: C recursion (from level 0, {n} levels of float64 accumulation) vs law-derived m_{n}: max|diff| = {np.max(np.abs(mc - m_law[n])):.1e} over L = {len(mc)}")
            del mc
sys.stdout.flush()

# ----------------------------------------------------------------------------------------------------------------------
# D. the zero-line census
# ----------------------------------------------------------------------------------------------------------------------
banner("D. the 3-adic zero-line census: depth d(u, Q) = v_3(u 2^Q -+ 1), box A (Python integers) and box B (C, residues mod 3^30)")


def v3(x: int) -> int:
    if x == 0:
        return 10 ** 9
    d = 0
    while x % 3 == 0:
        x //= 3
        d += 1
    return d


def poisson_tail(lam: float, k: int) -> float:
    s, term = 0.0, math.exp(-lam)
    for i in range(k):
        s += term
        term *= lam / (i + 1)
    return max(0.0, 1.0 - s)


print("null model: u 2^Q is a unit mod 3, so exactly one of u 2^Q - 1, u 2^Q + 1 is 0 mod 3 (never both): depth >= 1 with probability 1 is a")
print("  tautology; conditionally on the class, P(depth >= d) = 3^-(d-1) if u 2^Q mod 3^d is uniform on its class.  For fixed u the map Q -> u 2^Q mod 3^d")
print("  is a bijection onto the units every L_d = 2 3^(d-1) consecutive Q: when W >= L_d each u contributes 2 floor(W/L_d) + O(1) lines of depth >= d")
print("  DETERMINISTICALLY (equidistribution of 2, not randomness); the Poisson model is informative only when L_d > W (d >= 9 for W = 10^4, d >= 8 for W = 2000).")
for d in range(1, 12):
    print(f"   L_{d} = {2 * 3 ** (d - 1)}", end=";")
print()

# exact seed checks
print("exact seed checks (Python integers):")
print(f"   v_3(55 2^423 + 1) = {v3(55 * 2 ** 423 + 1)} (note 15); v_3(55 2^423 - 1) = {v3(55 * 2 ** 423 - 1)}")
print(f"   v_3(13 2^154 - 1) = {v3(13 * 2 ** 154 - 1)} (note 7); v_3(2^486 - 1) = {v3(2 ** 486 - 1)} (note 6; LTE: 1 + v_3(486/2) = {1 + v3(243)}); v_3(2^480 - 1) = {v3(2 ** 480 - 1)} (note 2; LTE: 1 + v_3(240) = {1 + v3(240)})")
print(f"   v_3(1187 2^5031 - 1) = {v3(1187 * 2 ** 5031 - 1)} (note: 1187 2^5031 = 1 mod 3^17); v_3(2441 2^8384 + 1) = {v3(2441 * 2 ** 8384 + 1)} (note: = -1 mod 3^16)")
for (u, Q, sgn) in ((3025, 846, -1), (1685, 8663, -1), (55, 423, +1)):
    print(f"   v_3({u} 2^{Q} {'+' if sgn > 0 else '-'} 1) = {v3(u * 2 ** Q + sgn)} (note: depth 15)")
for (u, Q) in ((4639, 4171), (4573, 2226), (4541, 7132), (3853, 9384), (3037, 5995), (1997, 873), (917, 1557), (521, 1365), (235, 8097)):
    x = u * 2 ** Q
    print(f"   ({u}, {Q}): v_3(u 2^Q - 1) = {v3(x - 1)}, v_3(u 2^Q + 1) = {v3(x + 1)} (note: depth 14)")
sys.stdout.flush()

# box A in Python
U, W = 1000, 2000
CAPA = 25
MA = 3 ** CAPA
us = [u for u in range(1, U + 1, 2) if u % 3]
hist = {}
lines = []
window = {}
for u in us:
    r = u % MA
    for Q in range(1, W + 1):
        r = (2 * r) % MA
        if r % 3 == 1:
            v, sgn = r - 1, "-"
        else:
            v, sgn = (r + 1) % MA, "+"
        d = v3(v) if v else CAPA
        hist[d] = hist.get(d, 0) + 1
        if d >= 11:
            lines.append((d, u, Q, sgn))
        if Q <= 60 and d >= 8:
            window.setdefault(Q, []).append((d, u, sgn))
N = len(us) * W
print(f"box A: {len(us)} multipliers u <= {U}, Q <= {W}: N = {N} pairs")
dmax = max(hist)
for d in range(1, dmax + 1):
    obs = hist.get(d, 0)
    exp_ = N * (2 / 3) * 3 ** (-(d - 1))
    cum = sum(hist.get(e, 0) for e in range(d, dmax + 1))
    cexp = N * 3 ** (-(d - 1))
    print(f"   d={d:2d}: {obs:7d} exp {exp_:10.2f} ratio {obs / exp_:6.3f} | >= d: {cum:7d} exp {cexp:9.3f} Poisson P(X >= obs) = {poisson_tail(cexp, cum) if cum else 1.0:.3g}")
print("   lines of depth >= 11:", sorted(lines, reverse=True))
print("   window Q <= 60, depth >= 8:", {Q: sorted(v, reverse=True) for Q, v in sorted(window.items())})
print(f"   note's box-A claims: 3 lines of depth >= 14 vs 0.42 (P = {poisson_tail(N * 3 ** -13, 3):.4f}); 4 of depth >= 13 vs 1.25 (P = {poisson_tail(N * 3 ** -12, 4):.4f})")
sys.stdout.flush()

# box B in C
ZL_C = r'''
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
/* usage: zl U W CAP : census of depth d(u,Q) = v_3(u 2^Q -+ 1) over odd u <= U prime to 3 and 1 <= Q <= W, residues mod 3^CAP */
int main(int argc, char **argv) {
  int U = atoi(argv[1]), W = atoi(argv[2]), CAP = atoi(argv[3]);
  uint64_t M = 1; for (int i = 0; i < CAP; i++) M *= 3;
  uint64_t hist[64] = {0};
  int nu = 0; for (int u = 1; u <= U; u += 2) if (u % 3) nu++;
  double *cnt6 = calloc(nu, sizeof(double)), *cnt9 = calloc(nu, sizeof(double)), *cnt10 = calloc(nu, sizeof(double));
  int iu = 0;
  printf("deep lines (depth >= 13): depth u Q sign\n");
  for (int u = 1; u <= U; u += 2) {
    if (u % 3 == 0) continue;
    uint64_t r = (uint64_t)u % M;
    for (int Q = 1; Q <= W; Q++) {
      r = (2 * r) % M;
      uint64_t v; char sgn;
      if (r % 3 == 1) { v = r - 1; sgn = '-'; } else { v = (r + 1) % M; sgn = '+'; }
      int d = 0;
      if (v == 0) d = CAP; else while (d < CAP && v % 3 == 0) { v /= 3; d++; }
      hist[d]++;
      if (d >= 6) cnt6[iu] += 1; if (d >= 9) cnt9[iu] += 1; if (d >= 10) cnt10[iu] += 1;
      if (d >= 13) printf("   %2d %5d %6d %c\n", d, u, Q, sgn);
    }
    iu++;
  }
  double N = (double)nu * W;
  printf("box B: %d multipliers u <= %d, Q <= %d: N = %.0f pairs, cap %d\n", nu, U, W, N, CAP);
  printf("depth d: observed, expected N (2/3) 3^-(d-1), ratio | cumulative >= d: observed, expected\n");
  for (int d = 1; d <= CAP; d++) {
    uint64_t cum = 0; for (int e = d; e <= CAP; e++) cum += hist[e];
    double ex = N * (2.0 / 3.0) * pow(3.0, -(d - 1)), cex = N * pow(3.0, -(d - 1));
    printf("   d=%2d: %9llu %13.3f %7.3f | >= d: %9llu %12.4f\n", d, (unsigned long long)hist[d], ex, hist[d] / ex, (unsigned long long)cum, cex);
  }
  for (int which = 0; which < 3; which++) {
    double *c = which == 0 ? cnt6 : which == 1 ? cnt9 : cnt10; int d = which == 0 ? 6 : which == 1 ? 9 : 10;
    double mean = 0, var = 0; for (int i = 0; i < nu; i++) mean += c[i]; mean /= nu; for (int i = 0; i < nu; i++) var += (c[i] - mean) * (c[i] - mean); var /= (nu - 1);
    printf("per-multiplier count of lines of depth >= %d: mean %.3f (W 3^-(d-1) = %.3f), variance %.3f (Poisson would be %.3f; L_d = %d %s W)\n", d, mean, W * pow(3.0, -(d - 1)), var, mean, 2 * (int)pow(3, d - 1), 2 * (int)pow(3, d - 1) <= W ? "<=" : ">");
  }
  return 0;
}
'''
src = os.path.join(SCRATCH, "zl_audit.c")
exe = os.path.join(SCRATCH, "zl_audit.exe")
with open(src, "w") as f:
    f.write(ZL_C)
cc = subprocess.run(["gcc", "-O2", "-o", exe, src, "-lm"], capture_output=True, text=True)
print("gcc (zl_audit):", "ok" if cc.returncode == 0 else cc.stderr)
if cc.returncode == 0:
    proc = subprocess.run([exe, "5000", "10000", "30"], capture_output=True, text=True)
    out = proc.stdout
    print(out)
    # Poisson tails for the cumulative counts, from the parsed histogram
    NB = 1667 * 10000
    cum = {}
    for line in out.splitlines():
        if line.strip().startswith("d=") and ">= d:" in line:
            d = int(line.split("d=")[1].split(":")[0])
            c = int(line.split(">= d:")[1].split()[0])
            cum[d] = c
    print("Poisson tail probabilities P(X >= observed) for the cumulative counts, box B:")
    for d in range(10, 21):
        if d in cum:
            lam = NB * 3 ** (-(d - 1))
            print(f"   d >= {d:2d}: observed {cum[d]:5d}, mean {lam:9.4f}, P = {poisson_tail(lam, cum[d]) if cum[d] else 1.0:.3f}")
    print(f"   maximum depth in a box of N pairs: E[max] ~ log_3 N + O(1) = {math.log(NB, 3):.1f}; observed deepest 17")
sys.stdout.flush()

# ----------------------------------------------------------------------------------------------------------------------
# F. the saddle chain
# ----------------------------------------------------------------------------------------------------------------------
banner("F. the saddle chain of 1 + 4^m, 1729, 17, 137, the -17 cycle")
from sympy import factorint, isprime, cyclotomic_poly, symbols


def T(n: int) -> int:
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def orbit(n: int, stop=1):
    out = [n]
    while n != stop:
        n = T(n)
        out.append(n)
    return out


ok = all(T(T(1 + 4 ** j * t)) == 1 + 3 * 4 ** (j - 1) * t for j in range(1, 8) for t in range(1, 200))
print(f"   T^2(1 + 4^j t) = 1 + 3 4^(j-1) t for j = 1..7, t = 1..199: {ok}  (x = 1 mod 4: T(x) = (3x+1)/2 = 2 + 3 2^(2j-1) t even for j >= 1, T^2 = 1 + 3 4^(j-1) t; induction gives T^(2j)(1 + 4^j t) = 1 + 3^j t)")
ok2 = all(all(T(T(x)) == 1 + 3 ** (j + 1) * 4 ** (mm - j - 1) for j, x in [(j, 1 + 3 ** j * 4 ** (mm - j)) for j in range(mm)]) for mm in range(1, 12))
print(f"   chain 1 + 3^j 4^(m-j) -> 1 + 3^(j+1) 4^(m-j-1) under T^2 for m <= 11: {ok2}")
chain = [1 + 3 ** j * 4 ** (6 - j) for j in range(7)]
it = [4097]
for _ in range(6):
    it.append(T(T(it[-1])))
print(f"   m = 6 chain 1 + 3^j 4^(6-j): {chain}; T^2-iterates of 4097: {it}; equal: {it == chain}")
for c in chain:
    print(f"      {c:5d} = {dict(factorint(c))}, prime: {isprime(c)}; mod 3 = {c % 3}, mod 4 = {c % 4}")
print(f"   1729 = 1 + 12^3: {1 + 12 ** 3 == 1729}; 1 + 3^3 4^3 = {1 + 27 * 64}; 9^3 + 10^3 = {9 ** 3 + 10 ** 3}; 3^6 = {3 ** 6} = 730 - 1: {3 ** 6 == 729}; 12^3 + 1^3 = {12 ** 3 + 1}")
print(f"   1 + 12^m midpoint of the chain of 1 + 4^(2m): " + ", ".join(f"m={mm}: {1 + 12 ** mm} (j={mm} point = {1 + 3 ** mm * 4 ** mm})" for mm in range(1, 8)) + " (identical by definition: 1 + 3^m 4^m = 1 + 12^m)")
print("   T^(2m)(1 + 4^(2m)) = 1 + 12^m for m = 1..8 by iteration:", end=" ")
res = []
for mm in range(1, 9):
    x = 1 + 4 ** (2 * mm)
    for _ in range(mm):
        x = T(T(x))
    res.append(x == 1 + 12 ** mm)
print(res)
xs = symbols('x')
print(f"   Phi_6(x) = {cyclotomic_poly(6, xs)}; Phi_6(12) = {cyclotomic_poly(6, 12)} = {dict(factorint(int(cyclotomic_poly(6, 12))))}; 1729 = 13 * 133: {13 * 133}; primes of 1729 mod 6: {[p % 6 for p in factorint(1729)]}")
print(f"   1297 = 6^4 + 1 = {6 ** 4 + 1}, prime: {isprime(1297)}; 973 = 7 * 139: {7 * 139}; 139 = 3^7 - 2^11: {3 ** 7 - 2 ** 11}; 4 3^5 + 1 = {4 * 3 ** 5 + 1} = 7 (3^7 - 2^11) = {7 * (3 ** 7 - 2 ** 11)}")
print(f"   17 = 1 + 4^2: {1 + 16}; chain of 17 under T^2: {[17, T(T(17)), T(T(T(T(17))))]} = 1 + 3^j 4^(2-j): {[1 + 3 ** j * 4 ** (2 - j) for j in range(3)]}; 10 = 3^2 + 1")
o1729, o27 = orbit(1729), orbit(27)
common = next(v for v in o1729 if v in set(o27))
print(f"   orbit of 1729 under T: {len(o1729) - 1} steps, max {max(o1729)}; orbit of 27: {len(o27) - 1} steps, max {max(o27)} (9232 under the Collatz map C = 2 * 4616)")
print(f"   first common element {common} at position {o1729.index(common)} in the orbit of 1729 and {o27.index(common)} in that of 27; orbit of 1729 to it: {o1729[:o1729.index(common) + 1]}")
print(f"   137 = 1 + 8 * 17: {1 + 8 * 17}; T-orbit from 137: {orbit(137)[:12]} ... (the note's '137 -> 103' skips 206 = T(137))")
neg = [-17]
for _ in range(11):
    neg.append(T(neg[-1]))
odd_steps = sum(1 for v in neg[:-1] if v % 2)
print(f"   T-orbit of -17: {neg}; returns to -17: {neg[-1] == -17}; odd steps {odd_steps}, even steps {11 - odd_steps}; clock 2^11 - 3^7 = {2 ** 11 - 3 ** 7}; ")
# cycle constant: x 2^11 = 3^7 x + c
c_cyc = -17 * (2 ** 11 - 3 ** 7)
print(f"      cycle equation x (2^11 - 3^7) = c gives c = {c_cyc} = 17 * 139: {17 * 139}")
print(f"   the j = 6 row of the extended-Collatz Theorem 2.3 reads 730, 973, 1297, 1729, 2305, 3073, 1024: the last entry 1024 = (3073 - 1)/3 = {(3073 - 1) // 3} is the backward-greedy shrink (C-preimage), not the chain point 4097 = 1 + 4 * 1024 = {1 + 4 * 1024}")
sys.stdout.flush()

# ----------------------------------------------------------------------------------------------------------------------
# G. basins: bitmask sieve in C to 2^NBASIN, Python cross-check to 2^22, Krasikov-Lagarias arithmetic
# ----------------------------------------------------------------------------------------------------------------------
banner(f"G. basins B(a) = {{n : the T-orbit of n passes through a}} on [1, 2^{NBASIN}] by an independent bitmask sieve (every target on the orbit is recorded)")
BS_C = r'''
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
/* usage: bs S : mask[n] = bitmask of the targets on the T-orbit of n, n <= 2^S; bits 0..8 = 4097, 3073, 2305, 1729, 1297, 973, 730, 27, 137 */
static const uint64_t targets[9] = {4097, 3073, 2305, 1729, 1297, 973, 730, 27, 137};
int main(int argc, char **argv) {
  int S = atoi(argv[1]); uint64_t N = 1ULL << S;
  uint16_t *mask = calloc(N + 1, sizeof(uint16_t)); if (!mask) { fprintf(stderr, "alloc fail\n"); return 1; }
  for (uint64_t n = 1; n <= N; n++) {
    u128 v = n; uint16_t acc = 0;
    for (;;) {
      if (v <= 4097) for (int i = 0; i < 9; i++) if (v == targets[i]) acc |= (uint16_t)(1u << i);
      if (v < n) { acc |= mask[(uint64_t)v]; break; }
      if (v == 1) break;
      v = (v & 1) ? (3 * v + 1) / 2 : v >> 1;
    }
    mask[n] = acc;
  }
  printf("range        B(4097)  B(3073)  B(2305)  B(1729)  B(1297)  B(973)   B(730)   B(27)    B(137)  | |B(1729)| vs X^0.84\n");
  uint64_t cum[9] = {0}, viol_chain = 0, viol_137 = 0, viol_27 = 0;
  for (int s = 0; s < S; s++) {
    uint64_t lo = 1ULL << s, hi = (1ULL << (s + 1)) - 1;
    for (uint64_t n = lo; n <= hi; n++) {
      uint16_t mk = mask[n];
      for (int i = 0; i < 9; i++) if (mk & (1u << i)) cum[i]++;
      for (int i = 0; i < 6; i++) if ((mk & (1u << i)) && !(mk & (1u << (i + 1)))) viol_chain++;
      if ((mk & 0x7F) && !(mk & (1u << 8))) viol_137++;
      if ((mk & (1u << 7)) && !(mk & (1u << 8))) viol_27++;
    }
    if (s + 1 >= 20) {
      double X = (double)(hi + 1);
      printf("[1,2^%2d]  ", s + 1);
      for (int i = 0; i < 9; i++) printf(" %.6f", (double)cum[i] / X);
      printf("  | %llu vs %.4e ratio %.4f\n", (unsigned long long)cum[3], pow(X, 0.84), (double)cum[3] / pow(X, 0.84));
    }
  }
  /* also the count at X = 2^S exactly (n <= 2^S): the number 2^S itself is a power of two: not in any basin */
  printf("nestedness violations: chain B(1+3^j 4^(6-j)) not subset of B(next): %llu; chain point without 137: %llu; 27 without 137: %llu\n", (unsigned long long)viol_chain, (unsigned long long)viol_137, (unsigned long long)viol_27);
  printf("B(27) = {");
  int c27 = 0; for (uint64_t n = 1; n <= N; n++) if (mask[n] & (1u << 7)) { printf("%s%llu", c27 ? "," : "", (unsigned long long)n); c27++; }
  printf("} : %d numbers\n", c27);
  printf("top dyadic range [2^%d, 2^%d): densities", S - 1, S);
  { uint64_t lo = 1ULL << (S - 1), hi = (1ULL << S) - 1, c[9] = {0}; for (uint64_t n = lo; n <= hi; n++) for (int i = 0; i < 9; i++) if (mask[n] & (1u << i)) c[i]++;
    for (int i = 0; i < 9; i++) printf(" %.6f", (double)c[i] / (double)(hi - lo + 1)); printf("\n"); }
  /* side entries of 1729 and 730 */
  { uint64_t c1729_not2305 = 0, c730_not973 = 0; for (uint64_t n = 1; n <= N; n++) { uint16_t mk = mask[n]; if ((mk & (1u << 3)) && !(mk & (1u << 2))) c1729_not2305++; if ((mk & (1u << 6)) && !(mk & (1u << 5))) c730_not973++; }
    printf("B(1729) \\ B(2305): %llu (density %.6f); B(730) \\ B(973): %llu (density %.6f)\n", (unsigned long long)c1729_not2305, (double)c1729_not2305 / (double)N, (unsigned long long)c730_not973, (double)c730_not973 / (double)N); }
  free(mask); return 0;
}
'''
src = os.path.join(SCRATCH, "bs_audit.c")
exe = os.path.join(SCRATCH, "bs_audit.exe")
with open(src, "w") as f:
    f.write(BS_C)
cc = subprocess.run(["gcc", "-O2", "-o", exe, src, "-lm"], capture_output=True, text=True)
print("gcc (bs_audit):", "ok" if cc.returncode == 0 else cc.stderr)
sys.stdout.flush()
count1729 = None
if cc.returncode == 0:
    proc = subprocess.run([exe, str(NBASIN)], capture_output=True, text=True)
    print(proc.stdout)
    for line in proc.stdout.splitlines():
        if line.startswith(f"[1,2^{NBASIN:2d}]"):
            count1729 = int(line.split("|")[1].split("vs")[0].strip())
            dens = [float(t) for t in line.split("|")[0].split()[1:]]
    print(f"   {stamp()}")

# Python cross-check at 2^22 (same bitmask bookkeeping, independent implementation)
SP = 22
NP_ = 1 << SP
targets = [4097, 3073, 2305, 1729, 1297, 973, 730, 27, 137]
tbit = {t: 1 << i for i, t in enumerate(targets)}
maskp_l = [0] * (NP_ + 1)
for n in range(1, NP_ + 1):
    v, acc = n, 0
    while True:
        if v <= 4097:
            b = tbit.get(v)
            if b:
                acc |= b
        if v < n:
            acc |= maskp_l[v]
            break
        if v == 1:
            break
        v = (3 * v + 1) // 2 if v & 1 else v // 2
    maskp_l[n] = acc
maskp = np.array(maskp_l, dtype=np.uint16)
print(f"Python bitmask sieve on [1, 2^{SP}]: densities", " ".join(f"{np.count_nonzero(maskp & (1 << i)) / NP_:.6f}" for i in range(9)),
      f"; |B(1729)| = {np.count_nonzero(maskp & 8)}; |B(27)| = {np.count_nonzero(maskp & 128)}   {stamp()}")
print("   (compare with the C sieve's [1,2^22] line above; the repo's basins output gives 0.003362 0.003537 0.004178 0.004421 0.004907 0.005839 0.014542 0.000004 0.300484 at 2^22)")

# Krasikov-Lagarias arithmetic
banner("G'. Krasikov-Lagarias arithmetic at the root 1729")
X = 2.0 ** NBASIN
if count1729 is not None:
    ratio = count1729 / X ** 0.84
    d = count1729 / X
    print(f"   |B(1729) ∩ [1, 2^{NBASIN}]| = {count1729} (note: 4,719,191 at 2^30); X^0.84 = {X ** 0.84:.4e} (note 3.85e7); ratio {ratio:.4f} (note 12.2%)")
    print(f"   per doubling the ratio grows by 2^0.16 = {2 ** 0.16:.4f} if the density is constant; doublings to ratio 1: {-math.log(ratio) / math.log(2 ** 0.16):.2f} -> X = 2^{NBASIN - math.log(ratio) / math.log(2 ** 0.16):.1f}")
    print(f"   equivalently X* = dens^(-1/0.16) = {d:.6f}^(-6.25) = {d ** -6.25:.3e} = 2^{math.log2(d ** -6.25):.2f} (note: 2^49.3 ≈ 7e14)")
    print(f"   2^49.3 = {2 ** 49.3:.3e}; 2^48.9 = {2 ** 48.9:.3e}; 3 2^53 (Oliveira e Silva 1999, exhaustive check at the time of the theorem) = {3 * 2 ** 53:.3e} = 2^{math.log2(3 * 2 ** 53):.1f}; 2^68 (Barina) = {2.0 ** 68:.3e}")
    print(f"   the count fails the inequality at X = 2^{NBASIN}: {count1729} < {X ** 0.84:.4e}: so any X_0(1729) for which |B ∩ [1,X]| >= X^0.84 holds for ALL X >= X_0 satisfies X_0 > 2^{NBASIN}")
print("   the theorem (Krasikov-Lagarias 2003, Thm 1.1 as recalled): for a != 0 mod 3, pi_a(x) := #{n <= x : T^k(n) = a for some k >= 0} >= x^0.84 for all sufficiently large x.")
print("   B(a) here = {n >= 1 : T-orbit passes through a} (a itself included, k = 0): the same set for odd a; for the even chain point 730 the Collatz map C")
print(f"   also reaches 730 as 3 243 + 1 from 243 = 3^5 (C-basin adds the doubling ray of 243: {sum(1 for k in range(64) if 243 * 2 ** k <= 2 ** 30)} numbers below 2^30), a density-zero difference.")
sys.stdout.flush()

banner("DONE")
print(stamp())
