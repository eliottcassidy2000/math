#!/usr/bin/env python3
"""Audit 4 (independent) of section 2d (S22) of 05-knowledge/results/collatz_five_mirrors_20260929.md.

Everything is recomputed from the definitions (Y_0 = 0, Y_{k+1} = 2^-a (3 Y_k + 1) in Z_3, P(a) = 2^-a; mu_hat_n(t) =
E e(t Y_n / 3^n) with Y_n read in [0, 3^n)); none of the repo's S22 scripts is imported.  Sections:
  A  the constants of section 2d (theta*, I, e^-I = 3^(h*-1), 1/(m-1), the Theorem C bracket, the J-window width)
  B  Lemma G: the 2-adic reading as an INTEGER identity, the two-term bound over all units n <= 9, the class
     identification (xi mod 4 against t mod 4 and the parity of n), the resonant units, the digit corollary
  C  the exponent-walk identity, the real-phase identity (A <= s) and the Ramanujan average, by brute force at h = 5
  D  the law mu_h itself on Z/3^h (forward DP, h <= 10): mu_hat_h(2^k) against the closed recursion and the walk DP
  E  the negative family by the closed recursion (window [-60, 0], a <= 50; against a <= 40, 60): Ntilde_n, rates
  F  Lemma R' term by term with an independent walk DP (state = (suffix sum T_j, counter)) at h = 20, 40, 80, 120;
     J_eff, phases, per-mass weights; the floor counter F
  G  the mass law: exact masses against Chernoff, the Chernoff-normalised growth factors
  H  the bound B_h(s) = sum_J mass_h(J) Ntilde_(J-1) at EVERY h <= 81 and EVERY 0 <= s <= floor(h log2 3), against |m_h(s)|
  I  the profile: ceiling side against 3^(-n/2), floor side against 2^(-d) M(n)
  J  a 30-digit mpmath check of the closed recursion on the negative family (float64 rounding, same truncation)
Run: python3 04-computation/experiments/collatz_five_mirrors_20260929_audit4.py   (a few minutes)
"""
from __future__ import annotations

import cmath
import itertools
import math
import time
from fractions import Fraction

import numpy as np

T0 = time.time()
LOG23 = math.log2(3.0)
MM = LOG23                      # m
PP = 1.0 / MM                   # log_3 2
QQ = 1.0 - PP
THETA = math.log(2.0 * QQ)      # theta*
IRATE = THETA * MM - math.log(MM - 1.0)
HSTAR = -PP * math.log2(PP) - QQ * math.log2(QQ)


def stamp() -> str:
    return f"[{time.time() - T0:4.0f}s]"


def ph(r: int, mod: int) -> complex:
    """e(r/mod) with the quotient correctly rounded."""
    return cmath.exp(2j * math.pi * float(Fraction(r, mod)))


def hdr(s: str) -> None:
    print(f"\n== {s} == {stamp()}")


# ----------------------------------------------------------------------------------------------------------------
hdr("A. constants of section 2d")
print(f"   m = log2 3 = {MM:.6f}, p = log_3 2 = {PP:.6f}, q = 1 - p = {QQ:.6f}, 2q = e^theta* = {2 * QQ:.6f} (note 0.73814)")
print(f"   theta* = ln(2q) = {THETA:.6f} (note -0.30362); I = theta* m - ln(m-1) = {IRATE:.6f} (note 0.054979)")
print(f"   e^-I = {math.exp(-IRATE):.6f}; 3^(h*-1) = {3 ** (HSTAR - 1):.6f}, h* = {HSTAR:.6f}; |difference| = {abs(math.exp(-IRATE) - 3 ** (HSTAR - 1)):.1e}")
lhs = -math.log(PP) - (QQ / PP) * math.log(QQ) - math.log(3)
rhs = -(1 / PP) * math.log(2 * QQ) + math.log(QQ) - math.log(PP)
print(f"   ln 3^(h*-1) = -ln p - (q/p) ln q - ln 3 = {lhs:.10f}; -I = -(1/p) ln(2q) + ln q - ln p = {rhs:.10f}; -I direct = {-IRATE:.10f}")
print(f"   Z(theta*) = q/(1-q) = {QQ / (1 - QQ):.6f} = m - 1 = {MM - 1:.6f}; tilted mean 1/(1 - e^theta*/2) = {1 / (1 - QQ):.6f} = m")
print(f"   1/(m-1) = {1 / (MM - 1):.6f} (note 1.70951); 3^(-1/2)/(m-1) = {3 ** -0.5 / (MM - 1):.6f} (note 0.98699)")
rho = 3 ** -0.5
br = (1 / (MM - 1)) / (1 - rho / (MM - 1))
print(f"   Theorem C bracket at rho = 3^(-1/2): 1 + {br:.1f} C (note 1 + 131.4 C); with C = 1.56 (max Ntilde_n 3^(n/2)): {1 + br * 1.56:.0f}")
print(f"   (3/2)/|ln(rho/(m-1))| = {1.5 / abs(math.log(rho / (MM - 1))):.1f} (note ~115); (1/(1 - rho/(m-1)))^2 = {(1 / (1 - rho / (MM - 1))) ** 2:.0f} (note ~6000)")
var = QQ / (1 - QQ) ** 2
print(f"   tilted variance q/p^2 = {var:.4f}, sd {math.sqrt(var):.4f}; J-window width sd/m = {math.sqrt(var) / MM:.3f} sqrt(h) (note 0.61 sqrt(h)); centre delta/m at delta = 6.4: {6.4 / MM:.2f}")

# ----------------------------------------------------------------------------------------------------------------
hdr("B. Lemma G (one-step gap)")


def G_2adic(n: int, t: int, R: int = 60) -> complex:
    mod = 3 ** n
    s = 0j
    for r in range(1, R + 1):
        m_r = (-t * pow(3, -n, 2 ** r)) % 2 ** r
        s += 2.0 ** (-r) * cmath.exp(2j * math.pi * (t / (2 ** r * mod) + m_r / 2 ** r))
    return s


def G_mod(n: int, t: int, R: int = 60) -> complex:
    mod = 3 ** n
    return sum(2.0 ** (-r) * ph((t * pow(2, -r, mod)) % mod, mod) for r in range(1, R + 1))


bad = checked = 0
for n in range(1, 7):
    mod = 3 ** n
    for t in range(1, mod):
        if t % 3 == 0:
            continue
        for r in range(1, 21):
            m_r = (-t * pow(3, -n, 2 ** r)) % 2 ** r
            num = t + m_r * mod
            checked += 1
            if num % 2 ** r != 0 or num // 2 ** r != (t * pow(2, -r, mod)) % mod:
                bad += 1
print(f"   integer identity (t 2^-r mod 3^n) = (t + m_r 3^n)/2^r with m_r = (-t 3^-n) mod 2^r in [0, 2^r): {checked} cases (units, n <= 6, r <= 20), failures {bad}")
for n in range(2, 6):
    mod = 3 ** n
    err = max(abs(G_2adic(n, t) - G_mod(n, t)) for t in range(1, mod) if t % 3)
    print(f"   n={n}: max |G_2adic - G_modular| over units = {err:.1e}")
GAP = (1 + math.sqrt(5)) / 4
print(f"   two-term bound: |G| <= 1/4 + sqrt(5/16 + cos(2 pi phi)/4), phi = (2 xi_2 - xi_1 - theta)/4; gap value (1+sqrt5)/4 = {GAP:.4f}")
for n in range(2, 10):
    mod = 3 ** n
    sup, cnt = {}, {}
    worst = -1.0
    idbad = 0
    resonant = []
    for t in range(1, mod):
        if t % 3 == 0:
            continue
        g = abs(G_2adic(n, t))
        theta = t / mod
        xi4 = (-t * pow(3, -n, 4)) % 4
        xi1, xi2 = xi4 & 1, xi4 >> 1
        v = (t & -t).bit_length() - 1
        if xi1 != (t & 1) or xi4 != ((-1) ** (n + 1) * t) % 4:
            idbad += 1
        phi = (2 * xi2 - xi1 - theta) / 4
        bnd = 0.25 + math.sqrt(5 / 16 + math.cos(2 * math.pi * phi) / 4)
        worst = max(worst, g - bnd)
        if v == 1:
            key = "v=1"
        elif v == 0 and xi4 == 1:
            key = "odd,xi=1(4)"
        elif v == 0:
            key = "odd,xi=3(4),theta<=1/2" if theta <= 0.5 else "odd,xi=3(4),theta>1/2"
        else:
            key = "4|t,theta<=1/2" if theta <= 0.5 else "4|t,theta>1/2"
        sup[key] = max(sup.get(key, 0.0), g)
        cnt[key] = cnt.get(key, 0) + 1
        if g > 0.99:
            resonant.append(t)
    gapc = cnt.get("v=1", 0) + cnt.get("odd,xi=1(4)", 0)
    print(f"   n={n}: max(|G| - two-term bound) = {worst:+.1e} (must be <= 0); class-identification failures {idbad}; gap classes {gapc} of {2 * 3 ** (n - 1)} units (other {2 * 3 ** (n - 1) - gapc})")
    print("      sup|G| by class: " + "; ".join(f"{k}: {sup[k]:.4f} ({cnt[k]})" for k in sorted(sup)))
    if n == 9:
        desc = []
        for t in resonant:
            u = t if t <= mod - t else t - mod
            v2 = (abs(u) & -abs(u)).bit_length() - 1
            desc.append(f"{u:+d}={'+' if u > 0 else '-'}2^{v2}*{abs(u) >> v2}")
        print(f"      |G_9(t)| > 0.99 at {len(resonant)} units (signed representative t or t - 3^9): " + ", ".join(desc))
bad = tot = 0
for n in range(1, 10):
    mod = 3 ** n
    for m in range(1, 31):
        t = pow(2, -m, mod)
        xi4 = (-t * pow(3, -n, 4)) % 4
        eta = (-pow(3, -n, 2 ** (m + 2))) % 2 ** (m + 2)
        tot += 1
        if (xi4 & 1) != (eta >> m) & 1 or (xi4 >> 1) != (eta >> (m + 1)) & 1:
            bad += 1
print(f"   corollary (the digits xi_1, xi_2 of t = 2^-m mod 3^n are the digits m+1, m+2 of -3^-n in Z_2): {tot} cases, failures {bad}")
for n in (20, 40, 80):
    eta = (-pow(3, -n, 2 ** 70)) % 2 ** 70
    gaps = sum(1 for m in range(1, 61) if ((eta >> m) & 1) != ((eta >> (m + 1)) & 1))
    print(f"   n={n}: steps m <= 60 with a Lemma-G gap at (n, m): {gaps} of 60")

# ----------------------------------------------------------------------------------------------------------------
hdr("C. exponent-walk identity, real-phase identity, Ramanujan average (brute force)")
h = 5
mod = 3 ** h
w1 = w2 = 0.0
w3 = 0
for w in itertools.product(range(1, 6), repeat=h):
    A = sum(w)
    Y = 0
    for a in w:
        Y = (pow(2, -a, mod) * (3 * Y + 1)) % mod
    T = [0] * (h + 2)
    for j in range(h, 0, -1):
        T[j] = T[j + 1] + w[j - 1]
    Yr = sum(3 ** (h - j) * pow(2, -T[j], mod) for j in range(1, h + 1)) % mod
    w3 = max(w3, abs(Yr - Y))
    Pp = [0]
    for a in w:
        Pp.append(Pp[-1] + a)
    C = sum(3 ** (h - j) * 2 ** Pp[j - 1] for j in range(1, h + 1))
    for s in range(-8, 21):
        lhs = ph((pow(2, s, mod) * Y) % mod, mod)
        rhs = 1
        for j in range(1, h + 1):
            mj = 3 ** j
            rhs *= ph(pow(2, s - T[j], mj), mj)
        w1 = max(w1, abs(lhs - rhs))
        if A <= s:
            w2 = max(w2, abs(lhs - ph((2 ** (s - A) * C) % mod, mod)))
print(f"   h=5, all 3125 words with letters <= 5, s in [-8, 20]: max |e(2^s Y/3^h) - prod_j omega_j(s - T_j)| = {w1:.1e}; Y_h = sum_j 3^(h-j) 2^(-T_j) mod 3^h: max discrepancy {w3}; real-phase identity for A <= s: max error {w2:.1e}")
for n in range(1, 8):
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    avg = sum(ph(pow(2, k, mod), mod) for k in range(L)) / L
    print(f"   Ramanujan n={n}: (1/L_n) sum_k omega_n(k) = {avg.real:+.6f}{avg.imag:+.6f}i; mu(3^n)/L_n = {(-1 if n == 1 else 0) / L:+.4f}")


# ----------------------------------------------------------------------------------------------------------------
def closed_all(N: int, lo: int, hi: int, amax: int, keep=None):
    """m_n(k) = mu_hat_n(2^k mod 3^n) for k in [lo, hi] at every level n <= N (exact shrinking window; the only
    truncation is a <= amax).  Returns {n: array}, array index i <-> k = lo + i."""
    Wlo = lo - amax * N
    size = hi - Wlo + 1
    cur = np.ones(size, dtype=complex)
    w = 2.0 ** -np.arange(1, amax + 1)
    out = {}
    for n in range(1, N + 1):
        mod = 3 ** n
        lo_n = lo - amax * (N - n)
        klo = lo_n - amax
        inv2 = pow(2, -1, mod)
        r = pow(2, hi, mod)
        vals = [0] * (hi - klo + 1)
        for idx in range(hi - klo, -1, -1):
            vals[idx] = r
            r = (r * inv2) % mod
        phs = np.zeros(size, dtype=complex)
        phs[klo - Wlo: hi - Wlo + 1] = [ph(x, mod) for x in vals]
        pc = phs * cur
        new = np.zeros(size, dtype=complex)
        i0, i1 = lo_n - Wlo, hi - Wlo
        for a in range(1, amax + 1):
            new[i0:i1 + 1] += w[a - 1] * pc[i0 - a:i1 + 1 - a]
        cur = new
        if keep is None or n in keep:
            out[n] = cur[lo - Wlo: hi - Wlo + 1].copy()
    return out


def cost_laws(N: int, tmax: int) -> np.ndarray:
    """P[N', t] = P(S_N' = t), sums of N' i.i.d. geometric(1/2) on {1, 2, ...}, by the halving recursion."""
    P = np.zeros((N + 1, tmax + 1))
    P[0, 0] = 1.0
    for t in range(1, tmax + 1):
        P[1:, t] = 0.5 * (P[:-1, t - 1] + P[1:, t - 1])
    return P


def exact_masses(h: int, s: int, P: np.ndarray) -> np.ndarray:
    """mass_h(J) = P(T_(J+1) <= s < T_J), J = 0..h, with T_(J+1) = S_(h-J) and T_J - T_(J+1) = a_J."""
    out = np.zeros(h + 1)
    ts = np.arange(s + 1)
    wts = 2.0 ** (-(s - ts))
    out[0] = P[h, :s + 1].sum()
    out[1:] = P[:h][::-1][:, :s + 1] @ wts
    return out


def walk_dp(h: int, s: int, counter: str = "J", amax: int = 40, Tmax=None, K: int = 3):
    """Exact DP over the words from the top level down: state (T_j, counter).  Returns (c[.], mass[.]) indexed by the
    counter value: c[J] = E[1_{counter = J} prod_j omega_j(s - T_j)].  counter J: levels with s - T_j < 0;
    counter F: levels with 0 <= s - T_j and s - T_j > j log2 3 - K.  No window in the total cost A."""
    if Tmax is None:
        Tmax = 2 * h + 140
    W = np.zeros((Tmax + 1, h + 2), dtype=complex)
    W[0, 0] = 1.0
    Wm = np.zeros((Tmax + 1, h + 2))
    Wm[0, 0] = 1.0
    w = 2.0 ** -np.arange(1, amax + 1)
    Ts = np.arange(Tmax + 1)
    for j in range(h, 0, -1):
        mod = 3 ** j
        inv2 = pow(2, -1, mod)
        r = pow(2, s, mod)
        phs = np.empty(Tmax + 1, dtype=complex)
        for T in range(Tmax + 1):
            phs[T] = ph(r, mod)
            r = (r * inv2) % mod
        ks = s - Ts
        counted = (ks < 0) if counter == "J" else ((ks >= 0) & (ks > j * LOG23 - K))
        new = np.zeros_like(W)
        newm = np.zeros_like(Wm)
        for a in range(1, min(amax, Tmax) + 1):
            src = W[:Tmax + 1 - a]
            srcm = Wm[:Tmax + 1 - a]
            contrib = (w[a - 1] * phs[a:])[:, None] * src
            contribm = w[a - 1] * srcm
            cnt = counted[a:]
            unc = ~cnt
            tgt = new[a:]
            tgtm = newm[a:]
            tgt[unc] += contrib[unc]
            tgtm[unc] += contribm[unc]
            tgt[cnt, 1:] += contrib[cnt, :-1]
            tgtm[cnt, 1:] += contribm[cnt, :-1]
        W, Wm = new, newm
    return W.sum(axis=0)[:h + 1], Wm.sum(axis=0)[:h + 1]


# ----------------------------------------------------------------------------------------------------------------
hdr("D. the law on Z/3^h by forward DP (h <= 10) against the closed recursion and the walk DP")


def law(h: int, amax: int = 60) -> np.ndarray:
    mod = 3 ** h
    mu = np.zeros(mod)
    mu[0] = 1.0
    y = np.arange(mod, dtype=np.int64)
    z = (3 * y + 1) % mod
    for _ in range(h):
        new = np.zeros(mod)
        for a in range(1, amax + 1):
            idx = (z * pow(2, -a, mod)) % mod
            new += 2.0 ** (-a) * np.bincount(idx, weights=mu, minlength=mod)
        mu = new
    return mu


def fourier(mu: np.ndarray, t: int) -> complex:
    mod = len(mu)
    y = np.arange(mod, dtype=np.int64)
    return complex(np.sum(mu * np.exp(2j * np.pi * ((t * y) % mod) / mod)))


rec10 = closed_all(10, -60, 30, 50)
for h in range(1, 11):
    mu = law(h)
    mod = 3 ** h
    kc = int(math.floor(h * LOG23))
    ks = list(range(-60, kc + 8))
    err = max(abs(fourier(mu, pow(2, k, mod)) - rec10[h][k + 60]) for k in ks)
    Nt_law = sum(2.0 ** (-m) * abs(fourier(mu, pow(2, -m, mod))) for m in range(1, 61))
    Nt_rec = sum(2.0 ** (-m) * abs(rec10[h][60 - m]) for m in range(1, 61))
    line = f"   h={h:2d}: max |mu_hat_h(2^k) [law] - m_h(k) [recursion]| over k in [-60, {kc + 7}] = {err:.1e}; Ntilde_h law {Nt_law:.6e} rec {Nt_rec:.6e}"
    if h in (6, 8, 10):
        s = kc - 2 if h < 9 else kc - 4
        cJ, mJ = walk_dp(h, s, "J", amax=60, Tmax=200)
        direct = fourier(mu, pow(2, s, mod))
        line += f"; walk DP at s={s}: |sum_J c_J - mu_hat| = {abs(cJ.sum() - direct):.1e}, mass {mJ.sum():.12f}"
    print(line)

# ----------------------------------------------------------------------------------------------------------------
hdr("E. the negative family Ntilde_n = sum_(m>=1) 2^-m |m_n(-m)| (window [-60, 0]; m <= 60, tail <= 2^-60)")
NEG = {}
NT = {}
NS = {}
for amax in (40, 50, 60):
    N = 120 if amax == 50 else 80
    res = closed_all(N, -60, 0, amax)
    NEG[amax] = res
    NT[amax] = {0: 1.0}
    NS[amax] = {}
    for n in range(1, N + 1):
        arr = res[n]
        NT[amax][n] = float(sum(2.0 ** (-m) * abs(arr[60 - m]) for m in range(1, 61)))
        NS[amax][n] = float(max(abs(arr[60 - m]) for m in range(1, 61)))
    print(f"   recursion with a <= {amax} to n = {N} done {stamp()}")
for amax in (40, 60):
    rel = max(abs(NT[amax][n] - NT[50][n]) / NT[50][n] for n in range(1, 81))
    relc = max(abs(NEG[amax][n][60 - m] - NEG[50][n][60 - m]) / abs(NEG[50][n][60 - m]) for n in range(1, 81) for m in range(1, 61))
    print(f"   truncation: max_n |Ntilde(a<={amax}) - Ntilde(a<=50)|/Ntilde(a<=50) over n <= 80 = {rel:.1e}; single coefficients (m <= 60): {relc:.1e}")
print(f"   a priori truncation bound per level 2^-50 = {2.0 ** -50:.1e} against Ntilde_80 = {NT[50][80]:.2e}: NOT self-certifying a priori; certified only by the agreement of the three truncations")
print("   n : Ntilde_n : Ntilde_n 3^(n/2) : N_n 3^(n/2) (sup over m <= 60)")
for n in range(1, 121):
    if n <= 20 or n % 4 == 0 or n in (81,):
        print(f"   {n:3d}: {NT[50][n]:.4e}  {NT[50][n] * 3 ** (n / 2):.3f}  {NS[50][n] * 3 ** (n / 2):.3f}")
vals = {n: NT[50][n] * 3 ** (n / 2) for n in range(1, 121)}
for lo, hi in ((1, 80), (20, 80), (13, 19), (81, 120), (20, 120)):
    sub = {n: vals[n] for n in range(lo, hi + 1)}
    nmin = min(sub, key=sub.get)
    nmax = max(sub, key=sub.get)
    print(f"   Ntilde_n 3^(n/2) over n in [{lo}, {hi}]: min {sub[nmin]:.3f} (n={nmin}), max {sub[nmax]:.3f} (n={nmax})")
for n in (76, 77, 78, 79, 80):
    print(f"   Ntilde_n 3^(n/2) at n={n}: {vals[n]:.3f}")
for n in (84, 88, 96, 104, 112, 120):
    arr = NEG[50][n]
    mm = max(range(1, 61), key=lambda m: abs(arr[60 - m]))
    print(f"   n={n}: sup_m |m_n(-m)| 3^(n/2) = {NS[50][n] * 3 ** (n / 2):.2f} at m = {mm}; sup over m <= 20: {max(abs(arr[60 - m]) for m in range(1, 21)) * 3 ** (n / 2):.2f}; 2^-m |m_n(-m)| 3^(n/2) at the argmax: {2.0 ** -mm * abs(arr[60 - mm]) * 3 ** (n / 2):.2e}")
neg60b = closed_all(120, -20, 0, 60, keep=set(range(84, 121)))
relb = max(abs(neg60b[n][20 - m] - NEG[50][n][60 - m]) / abs(neg60b[n][20 - m]) for n in range(84, 121) for m in range(1, 21))
print(f"   truncation at n = 84..120, m <= 20: max relative |a<=50 - a<=60| = {relb:.1e}")
print("   |m_n(-16)| 3^(n/2) for n = 96..120 (the outlying coefficient of the sup at n = 120): " + " ".join(f"{abs(NEG[50][n][60 - 16]) * 3 ** (n / 2):.2f}" for n in range(96, 121, 2)))
svals = {n: NS[50][n] * 3 ** (n / 2) for n in range(1, 81)}
print(f"   N_n 3^(n/2) over n <= 80: min {min(svals.values()):.3f} (n={min(svals, key=svals.get)}), max {max(svals.values()):.3f} (n={max(svals, key=svals.get)})")
for lo, hi in ((20, 80), (40, 80), (20, 120), (60, 120), (81, 120)):
    ns = np.arange(lo, hi + 1)
    slope, icpt = np.polyfit(ns, [math.log(NT[50][n]) for n in ns], 1)
    print(f"   least-squares rate of Ntilde_n over {lo}..{hi}: {math.exp(slope):.4f} per level (3^-1/2 = {3 ** -0.5:.4f}, m - 1 = {MM - 1:.4f}); geometric-mean ratio {(NT[50][hi] / NT[50][lo]) ** (1 / (hi - lo)):.4f}")
NTV = NT[50]

# ----------------------------------------------------------------------------------------------------------------
hdr("F. Lemma R' term by term with the walk DP; J_eff; per-mass weights; the floor counter")
PROF = closed_all(81, -60, int(math.floor(81 * LOG23)) + 60, 50)
KLO = -60


def sstar_of(n: int, upper=None):
    """argmax_k |m_n(k)| over 0 <= k <= upper (default floor(n m) + 60) and the maximum."""
    arr = PROF[n]
    kc = int(math.floor(n * LOG23))
    if upper is None:
        upper = kc + 60
    ks = np.arange(KLO, KLO + len(arr))
    mask = (ks >= 0) & (ks <= upper)
    sub = np.abs(arr[mask])
    return int(ks[mask][np.argmax(sub)]), float(sub.max())


PC = cost_laws(120, 460)
SUMMARY_F = {}
for h in (20, 40, 80, 120):
    if h <= 81:
        s_star, Mh = sstar_of(h)
        full = PROF[h][s_star - KLO]
    else:
        s0 = int(round(h * LOG23)) - 6
        fam = closed_all(h, s0 - 8, s0 + 8, 50, keep={h})[h]
        i = int(np.argmax(np.abs(fam)))
        s_star, full = s0 - 8 + i, fam[i]
        Mh = abs(full)
    delta = h * LOG23 - s_star
    cJ, mJ = walk_dp(h, s_star, "J")
    em = exact_masses(h, s_star, PC)
    NtJ = np.array([1.0] + [NTV[J - 1] + 2.0 ** -60 if J >= 2 else 1.0 for J in range(1, h + 1)])
    rhs = em * NtJ
    ratio = np.where(rhs > 0, np.abs(cJ) / np.maximum(rhs, 1e-300), 0.0)
    ratio_complete = ratio.copy()
    ratio_complete[em < 1e-200] = 0.0
    jw = int(np.argmax(ratio_complete))
    print(f"   h={h}: s* = {s_star} (delta = {delta:.2f}), |m_h(s*)| = {Mh:.4e}; walk DP: |sum_J c_J - m_h(s*)| = {abs(cJ.sum() - full):.2e} ({abs(cJ.sum() - full) / Mh:.1e} relative), DP mass {mJ.sum():.10f}, max_J |DP mass - exact mass| = {np.max(np.abs(mJ - em)):.1e}")
    print(f"      worst |c_J|/(mass_h(J) Ntilde_(J-1)) = {ratio_complete.max():.3f} at J = {jw}; bound sum_J mass Ntilde = {rhs.sum():.4e} = {rhs.sum() / math.exp(-h * IRATE):.4f} e^-hI; |full|/bound = {Mh / rhs.sum():.4f}")
    print("      J : exact mass : |c_J| : ratio : |c_J|/|full| : arg c_J : per-mass x h : remainder |full - cum_J|/|full|")
    cum = np.cumsum(cJ)
    rem = np.abs(full - cum) / Mh
    for J in range(0, min(h, 30) + 1):
        print(f"      {J:2d}: {em[J]:.3e}  {abs(cJ[J]):.3e}  {ratio[J]:.3f}  {abs(cJ[J]) / Mh:.3f}  {cmath.phase(cJ[J]):+.2f}  {abs(cJ[J]) / em[J] * h if em[J] > 0 else 0:.2f}  {rem[J]:.3f}")
    jeff10 = next((J for J in range(h + 1) if rem[J] < 0.10), None)
    jeff10s = next((J for J in range(h + 1) if all(rem[J:] < 0.10)), None)
    jeff5s = next((J for J in range(h + 1) if all(rem[J:] < 0.05)), None)
    print(f"      J_eff: first J with remainder < 10%: {jeff10}; stable (< 10% from J on): {jeff10s}; stable 5%: {jeff5s}; 1.6 sqrt(h) = {1.6 * math.sqrt(h):.1f}, h/8 = {h / 8:.1f}, 8.2 ln h - 20 = {8.2 * math.log(h) - 20:.1f}")
    SUMMARY_F[h] = (s_star, Mh, ratio_complete.max(), rhs.sum() / math.exp(-h * IRATE), jeff10, jeff10s, [abs(cJ[J]) / em[J] * h for J in range(2, 11)])
    if h in (40, 80):
        cF, mF = walk_dp(h, s_star, "F")
        print(f"      floor counter F (K = 3): mass with F >= 1: {1 - mF[0] / mF.sum():.4f} of the DP mass ({mF.sum():.6f}); |c_(F=0)|/|full| = {abs(cF[0]) / Mh:.4f}; vector remainder |full - c_(F=0)|/|full| = {abs(full - cF[0]) / Mh:.4f}; sum_(F>=1) |c_F| / |full| = {np.abs(cF[1:]).sum() / Mh:.4f}")
    print(f"      {stamp()}")
print("   per-mass weights |c_J|/mass_h(J) x h, J = 2..10:")
for h in SUMMARY_F:
    print(f"      h={h:3d}: " + " ".join(f"{v:.2f}" for v in SUMMARY_F[h][6]))

# ----------------------------------------------------------------------------------------------------------------
hdr("G. the mass law: exact masses against Chernoff, normalised growth factors")
for h in (40, 80, 120, 300, 1000):
    s = int(round(h * LOG23)) - 6
    delta = h * LOG23 - s
    P = cost_laws(h, s + 1)
    em = exact_masses(h, s, P)
    chern = np.array([math.exp(-h * IRATE) * math.exp(THETA * delta) * (MM - 1) ** (-J) for J in range(h + 1)])
    tail = np.array([P[h - J, :s + 1].sum() for J in range(h + 1)])
    norm = em / (math.exp(-h * IRATE) * math.exp(THETA * delta))
    ratios = norm[1:] / np.maximum(norm[:-1], 1e-300)
    print(f"   h={h}: s = {s}, delta = {delta:.2f}; max_J (mass - Chernoff) = {np.max(em - chern):+.1e}; max_J (P(T_(J+1) <= s) - Chernoff) = {np.max(tail - chern):+.1e}; sum_J mass = {em.sum():.12f}")
    print(f"      normalised masses J = 0..12: " + " ".join(f"{v:.2e}" for v in norm[:13]))
    print(f"      ratios J/(J-1), J = 1..12: " + " ".join(f"{r:.2f}" for r in ratios[:12]) + f"; argmax_J mass = {int(np.argmax(em))} (0.21 h = {0.21 * h:.0f}); mass-weighted mean J = {(np.arange(h + 1) * em).sum():.1f}")

# ----------------------------------------------------------------------------------------------------------------
hdr("H. B_h(s) = sum_J mass_h(J) Ntilde_(J-1) at every h <= 81 and every 0 <= s <= floor(h log2 3), against |m_h(s)|")
print("   Ntilde_n taken from the a <= 50 recursion plus the tail bound 2^-60.  Columns: h, s*, delta*, M(h), B(s*)/e^-hI, max_s B(s)/e^-hI (at s), B(floor(hm))/e^-hI, max_s |m_h(s)|/B(s) (at s), |m_h(s*)|/B(s*)")
PC81 = cost_laws(81, 140)
worst_all = 0.0
maxB_star = {}
maxB_all = {}
for h in range(1, 82):
    kc = int(math.floor(h * LOG23))
    arr = PROF[h]
    absr = np.abs(arr)
    s_star, Mh = sstar_of(h, kc)
    s_glob, M_glob = sstar_of(h)
    if s_glob != s_star:
        print(f"   h={h:2d}: global argmax over the powers of two is k = {s_glob} > floor(h m) = {kc} (|m| = {M_glob:.3e}); the bound is evaluated for 0 <= s <= floor(h m) only")
    NtJ = np.array([1.0] + [NTV[J - 1] + 2.0 ** -60 if J >= 2 else 1.0 for J in range(1, h + 1)])
    Bs = np.zeros(kc + 1)
    for s in range(kc + 1):
        Bs[s] = float((exact_masses(h, s, PC81) * NtJ).sum())
    rat = absr[np.arange(kc + 1) - KLO] / Bs
    eI = math.exp(-h * IRATE)
    worst_all = max(worst_all, rat.max())
    maxB_star[h] = Bs[s_star] / eI
    maxB_all[h] = Bs.max() / eI
    if h <= 12 or h % 4 == 0 or h in (81, 41, 61):
        print(f"   h={h:2d}: s*={s_star:3d} d*={h * LOG23 - s_star:5.2f} M={Mh:.3e}  B(s*)/e^-hI={Bs[s_star] / eI:.4f}  max_s B/e^-hI={Bs.max() / eI:.4f} (s={int(np.argmax(Bs))}, d={h * LOG23 - int(np.argmax(Bs)):.2f})  B(kc)/e^-hI={Bs[kc] / eI:.4f}  max_s |m|/B={rat.max():.4f} (s={int(np.argmax(rat))})  |m(s*)|/B(s*)={rat[s_star]:.4f}")
print(f"   Lemma R' over ALL h <= 81 and ALL 0 <= s <= floor(h m): max |m_h(s)|/B_h(s) = {worst_all:.4f} (must be <= 1)")
print(f"   B_h(s*)/e^-hI over 40 <= h <= 81: max {max(maxB_star[h] for h in range(40, 82)):.4f} (h={max(range(40, 82), key=lambda h: maxB_star[h])}), min {min(maxB_star[h] for h in range(40, 82)):.4f}; over 20 <= h <= 39: max {max(maxB_star[h] for h in range(20, 40)):.4f}; over h <= 19: max {max(maxB_star[h] for h in range(1, 20)):.4f}")
print(f"   max_s B_h(s)/e^-hI (all 0 <= s <= floor(h m)): over 40 <= h <= 81: max {max(maxB_all[h] for h in range(40, 82)):.4f} (h={max(range(40, 82), key=lambda h: maxB_all[h])}); over 20 <= h <= 39: {max(maxB_all[h] for h in range(20, 40)):.4f}; over 1 <= h <= 81: {max(maxB_all.values()):.4f} (h={max(maxB_all, key=maxB_all.get)})")
print("   unconditional tails: without any Ntilde_n beyond n = N0 the bound carries the extra term sum_(J-1 > N0) mass_h(J) (Ntilde <= 1):")
for h, sh in ((100, 152), (120, 184), (150, 231), (200, 310), (250, 390), (300, 469)):
    P = cost_laws(h, sh + 1)
    em = exact_masses(h, sh, P)
    b = sum(em[J] * (1.0 if J <= 1 else NTV[J - 1] + 2.0 ** -60) for J in range(min(h, 121) + 1))
    t80 = em[82:].sum() if h > 81 else 0.0
    t120 = em[122:].sum() if h > 121 else 0.0
    eI = math.exp(-h * IRATE)
    print(f"      h={h}: bound (Ntilde to 120)/e^-hI = {b / eI:.4f}; tail sum_(J>=82) mass = {t80:.2e} = {t80 / eI:.2e} e^-hI; tail sum_(J>=122) mass = {t120:.2e} = {t120 / eI:.2e} e^-hI")
print("   bound at the argmax s of the level-300 run with Ntilde_n computed to n = 120 and 1.292 * 3^(-n/2) beyond:")
CN = max(NTV[n] * 3 ** (n / 2) for n in range(20, 81))
for h, sh in ((100, 152), (120, 184), (200, 310), (300, 469)):
    P = cost_laws(h, sh + 1)
    em = exact_masses(h, sh, P)
    b_comp = b_ext = b_ext81 = 0.0
    for J in range(h + 1):
        if J <= 1:
            nt = 1.0
        elif J - 1 <= 120:
            nt = NTV[J - 1] + 2.0 ** -60
        else:
            nt = CN * 3 ** (-(J - 1) / 2)
        term = em[J] * nt
        if J - 1 <= 80:
            b_comp += term
        elif J - 1 <= 120:
            b_ext81 += term
        else:
            b_ext += term
    tot = b_comp + b_ext81 + b_ext
    print(f"      h={h}: bound/e^-hI = {tot / math.exp(-h * IRATE):.4f}; share from J-1 <= 80 (computed in the note): {b_comp / tot:.4f}; from 81 <= J-1 <= 120 (computed here): {b_ext81 / tot:.4f}; extrapolated beyond 120: {b_ext / tot:.1e}")

# ----------------------------------------------------------------------------------------------------------------
hdr("H2. the same for 82 <= h <= 150 with Ntilde_n computed to n = 120 (tail J - 1 > 120 bounded by Ntilde <= 1)")
PROF150 = closed_all(150, -60, int(math.floor(150 * LOG23)) + 2, 50, keep=set(range(82, 151)))
PC150 = cost_laws(150, 245)
rows = []
for h in range(82, 151):
    kc = int(math.floor(h * LOG23))
    arr = PROF150[h]
    absr = np.abs(arr)
    ks = np.arange(-60, -60 + len(arr))
    mask = (ks >= 0) & (ks <= kc)
    s_star = int(ks[mask][np.argmax(absr[mask])])
    Mh = float(absr[mask].max())
    NtJ = np.array([1.0] + [NTV[J - 1] + 2.0 ** -60 if 2 <= J <= 121 else 1.0 for J in range(1, h + 1)])
    Bs = np.zeros(kc + 1)
    for s in range(kc + 1):
        Bs[s] = float((exact_masses(h, s, PC150) * NtJ).sum())
    rat = absr[np.arange(kc + 1) + 60] / Bs
    eI = math.exp(-h * IRATE)
    rows.append((h, s_star, h * LOG23 - s_star, Mh, Bs[s_star] / eI, Bs.max() / eI, rat.max(), int(np.argmax(rat)), rat[s_star]))
    if h % 10 == 0 or h in (82, 150):
        print(f"   h={h:3d}: s*={s_star} d*={h * LOG23 - s_star:.2f} M={Mh:.3e}  B(s*)/e^-hI={Bs[s_star] / eI:.4f}  max_s B/e^-hI={Bs.max() / eI:.4f}  max_s |m|/B={rat.max():.4f} (s={int(np.argmax(rat))})  |m(s*)|/B(s*)={rat[s_star]:.4f}")
print(f"   82 <= h <= 150: max_h B_h(s*)/e^-hI = {max(r[4] for r in rows):.4f} (h={max(rows, key=lambda r: r[4])[0]}), min {min(r[4] for r in rows):.4f}; max_h max_s B_h(s)/e^-hI = {max(r[5] for r in rows):.4f}; max_h max_s |m_h(s)|/B_h(s) = {max(r[6] for r in rows):.4f} (h={max(rows, key=lambda r: r[6])[0]}, s={max(rows, key=lambda r: r[6])[7]}); d* range [{min(r[2] for r in rows):.2f}, {max(r[2] for r in rows):.2f}]")

# ----------------------------------------------------------------------------------------------------------------
hdr("I. the profile |m_n(k)| on the powers of two: ceiling side and floor side")
for n in (40, 80):
    arr = PROF[n]
    kc = int(math.floor(n * LOG23))
    s_star, Mn = sstar_of(n)
    seg = np.abs(arr[np.arange(-60, 0) - KLO])
    print(f"   n={n}: M(n) = {Mn:.3e} at k = {s_star} = nm - {n * LOG23 - s_star:.2f}; ceiling k in [-60, -1]: rms {math.sqrt(np.mean(seg ** 2)):.3e}, max {seg.max():.3e} against 3^(-n/2) = {3 ** (-n / 2):.3e}")
    print("      floor side, k = kc + d (kc = floor(nm)); columns: d, k - nm, |m_n(k)|/M(n), 2^-d, ratio")
    for d in range(-2, 13):
        k = kc + d
        v = abs(arr[k - KLO]) / Mn
        print(f"         d={d:3d}  k-nm={k - n * LOG23:+.2f}  {v:.3e}  {2.0 ** -d:.3e}  ratio {v / 2.0 ** -d:.3f}")
    seg2 = np.abs(arr[np.arange(kc + 21, kc + 61) - KLO])
    seg3 = np.abs(arr[np.arange(kc + 1, kc + 21) - KLO])
    print(f"      k in [kc+1, kc+20]: rms {math.sqrt(np.mean(seg3 ** 2)):.3e} max {seg3.max():.3e}; k in [kc+21, kc+60]: rms {math.sqrt(np.mean(seg2 ** 2)):.3e} max {seg2.max():.3e}; against 3^(-n/2) = {3 ** (-n / 2):.1e}")
    lo_side = np.abs(arr[np.arange(s_star - 7, s_star + 1) - KLO])
    print(f"      low side k = s*-7..s*: " + " ".join(f"{v:.2e}" for v in lo_side) + "; step ratios " + " ".join(f"{lo_side[i] / lo_side[i + 1]:.2f}" for i in range(7)) + f" (e^theta* = {math.exp(THETA):.3f})")

# ----------------------------------------------------------------------------------------------------------------
hdr("J. 30-digit mpmath check of the closed recursion on the negative family (same truncation, float64 rounding only)")
import mpmath as mp  # noqa: E402

mp.mp.dps = 30


def closed_mp(N: int, lo: int, hi: int, amax: int):
    Wlo = lo - amax * N
    size = hi - Wlo + 1
    cur = [mp.mpc(1)] * size
    w = [mp.mpf(2) ** (-a) for a in range(1, amax + 1)]
    for n in range(1, N + 1):
        mod = 3 ** n
        lo_n = lo - amax * (N - n)
        klo = lo_n - amax
        inv2 = pow(2, -1, mod)
        r = pow(2, hi, mod)
        vals = [0] * (hi - klo + 1)
        for idx in range(hi - klo, -1, -1):
            vals[idx] = r
            r = (r * inv2) % mod
        pc = [mp.mpc(0)] * size
        for idx, x in enumerate(vals):
            i = klo - Wlo + idx
            pc[i] = mp.expjpi(2 * mp.mpf(x) / mod) * cur[i]
        new = [mp.mpc(0)] * size
        for i in range(lo_n - Wlo, hi - Wlo + 1):
            acc = mp.mpc(0)
            for a in range(1, amax + 1):
                acc += w[a - 1] * pc[i - a]
            new[i] = acc
        cur = new
    return cur[lo - Wlo: hi - Wlo + 1]


for N, lo, amax in ((40, -10, 25), (60, -6, 20)):
    hp = closed_mp(N, lo, 0, amax)
    fl = closed_all(N, lo, 0, amax, keep={N})[N]
    rel = max(abs(complex(hp[i]) - fl[i]) / abs(complex(hp[i])) for i in range(len(hp)))
    print(f"   N={N}, k in [{lo}, 0], a <= {amax}: max relative |float64 - mpmath(30 digits)| = {rel:.1e}; |m_N(-1)| = {abs(complex(hp[-2])):.6e} (mp) vs {abs(fl[-2]):.6e} (float64)  {stamp()}")
print("DONE " + stamp())
