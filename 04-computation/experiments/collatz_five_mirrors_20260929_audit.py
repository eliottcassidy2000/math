#!/usr/bin/env python3
"""Independent audit of `collatz_five_mirrors_20260929.md` (S20): every number and every proof step recomputed
with own code.  No import from the session's scripts.

The 3-adic Syracuse law is built by the FORWARD recursion z = 2^-a (3y + 1) mod 3^n with weights 2^-a
(float64, a <= 64, to level 14; exact rationals to level 6 using the period L = 2 3^(n-1) of 2^-a mod 3^n),
not by the session's discrete-log correlation.

Sections (tags [A]..[H] refer to the output):
  [A] law: exact vs float; consistency (Y_n mod 3^m has the law of Y_m); the frequency recursion.
  [B] Fourier profile: Proposition 1 numerically (same M(h) at n = 10, 12, 14), M(h), argmax, mass per level,
      typical |mu_hat|^2 3^h, M(1) = 1/sqrt 3, |mu_hat_12(2^s)|, Parseval against E_units[rho^2].
  [C] Proposition 2 numerics: the lift has no primitive coefficients; |mu_hat| <= ell^1 distance; d(m,N) vs M(m+1).
  [D] Gauss sums: sup, argmax, mean, share > 0.9, the 2-adic reading; the stated inequality
      |G_j(2^s)| >= 1 - 2 pi 2^-K - 2^-s tested for all (j, s, K) and against the corrected tail 2^(1-s).
  [E] same-length collisions in the tree of 1: the table with valuation cap 24 (the note's), with cap 70
      (certified exact for all sums <= 70 + 3d), the minimal colliding cost sum at EVERY depth <= 5.
  [F] carry reciprocity: symbolic identity (sympy), the (1,2) example, the three cycles, the orbit of -17,
      -13801 on the integer cycle of x -> (3x + 139)/2^v; the 2-adic source class mod 2^A vs 2^(A+1).
  [G] the resonance at h = 10 by total cost and by last valuation (own joint recursion); LD rate function.
  [H] decay models on M(h), h <= 17: local exponents, ratio windows, pure / shifted power law / geometric fits.
Run: python 04-computation/experiments/collatz_five_mirrors_20260929_audit.py   (about 3 minutes, < 2 GB)
"""
from __future__ import annotations

import cmath
import math
import random
import time
from fractions import Fraction

import numpy as np

T0 = time.time()
LOG23 = math.log(3, 2)


def stamp() -> str:
    return f"[{time.time() - T0:.0f}s]"


# ----------------------------------------------------------------------------------------------------------
# [A] the law
# ----------------------------------------------------------------------------------------------------------
def forward_float(N: int, AMAX: int = 64) -> dict[int, np.ndarray]:
    """mu_n on Z/3^n, n = 0..N, by the forward step z = 2^-a (3y + 1) mod 3^n, weight 2^-a (a <= AMAX)."""
    laws = {0: np.array([1.0])}
    mu = laws[0]
    for n in range(1, N + 1):
        mod = 3 ** n
        y = np.arange(3 ** (n - 1), dtype=np.int64)
        base = (3 * y + 1) % mod
        new = np.zeros(mod)
        for a in range(1, AMAX + 1):
            inv = pow(2, -a, mod)
            z = (inv * base) % mod  # y -> 3y+1 mod 3^n is injective, so no duplicate indices
            new[z] += (2.0 ** -a) * mu
        laws[n] = new
        mu = new
    return laws


def forward_exact(N: int) -> dict[int, dict[int, Fraction]]:
    """Exact rational law to level N: the weight of the class a = r mod L is 2^-r / (1 - 2^-L)."""
    laws = {0: {0: Fraction(1)}}
    for n in range(1, N + 1):
        mod, L = 3 ** n, 2 * 3 ** (n - 1)
        cL = 1 / (1 - Fraction(1, 2 ** L))
        w = [None] + [Fraction(1, 2 ** r) * cL for r in range(1, L + 1)]
        invs = [None] + [pow(2, -r, mod) for r in range(1, L + 1)]
        new: dict[int, Fraction] = {}
        for y, p in laws[n - 1].items():
            base = (3 * y + 1) % mod
            for r in range(1, L + 1):
                z = (invs[r] * base) % mod
                new[z] = new.get(z, Fraction(0)) + w[r] * p
        laws[n] = new
    return laws


NF = 14
NE = 6
laws = forward_float(NF)
print(f"[A] float law by the forward recursion to level {NF} {stamp()}")
ex = forward_exact(NE)
print(f"[A] exact rational law to level {NE} {stamp()}")
print(f"[A1] level 1 exact: {dict(sorted(ex[1].items()))}  (mu_1(1) = 1/3, mu_1(2) = 2/3)")
for n in range(1, NE + 1):
    err = max(abs(float(ex[n].get(y, Fraction(0))) - laws[n][y]) for y in range(3 ** n))
    print(f"[A2] n={n}: max |float - exact| = {err:.1e}; total float mass {laws[n].sum():.15f}")
# consistency: Y_n mod 3^m has the law of Y_m (exact), and in float to level 14
ok = True
for n in range(1, NE + 1):
    for m in range(1, n):
        red: dict[int, Fraction] = {}
        for y, p in ex[n].items():
            red[y % 3 ** m] = red.get(y % 3 ** m, Fraction(0)) + p
        ok &= all(red.get(y, Fraction(0)) == ex[m].get(y, Fraction(0)) for y in range(3 ** m))
print(f"[A3] consistency (reduction of mu_n mod 3^m equals mu_m) EXACT for all 1 <= m < n <= {NE}: {ok}")
worst = 0.0
for n in range(2, NF + 1):
    for m in range(1, n):
        red = laws[n].reshape(3 ** (n - m), 3 ** m).sum(axis=0)
        worst = max(worst, float(np.abs(red - laws[m]).max()))
print(f"[A4] consistency in float64, all 1 <= m < n <= {NF}: max deviation {worst:.1e}")
# a Monte-Carlo-free proof check of the consistency: Y_n mod 3^h depends only on the last h valuations
# (Y_(k+1) = 2^-a (3 Y_k + 1): 3 Y_k mod 3^h depends on Y_k mod 3^(h-1) only), verified on random words
rng = random.Random(7)
ok = True
for _ in range(2000):
    n, h = rng.randint(2, 9), rng.randint(1, 8)
    h = min(h, n)
    w = [rng.randint(1, 6) for _ in range(n)]
    w2 = [rng.randint(1, 6) for _ in range(n - h)] + w[n - h:]

    def Y(word: list[int], mod: int) -> int:
        y = 0
        for a in word:
            y = (pow(2, -a, mod) * (3 * y + 1)) % mod
        return y

    ok &= Y(w, 3 ** n) % 3 ** h == Y(w2, 3 ** n) % 3 ** h == Y(w[n - h:], 3 ** h)
print(f"[A5] Y_n mod 3^h is the function F_h of the last h valuations (2000 random word pairs): {ok}")
# the frequency recursion mu_hat_n(t) = sum_a 2^-a e(t 2^-a / 3^n) mu_hat_(n-1)(t 2^-a mod 3^(n-1)), n = 6
n = 6
mod, pm = 3 ** n, 3 ** (n - 1)
mh6, mh5 = np.fft.fft(laws[6]), np.fft.fft(laws[5])
worst = 0.0
for t in (1, 2, 5, 7, 100, 728, 3 ** 5, 2 * 3 ** 4):
    s = 0
    for a in range(1, 80):
        inv = pow(2, -a, mod)
        s += 2.0 ** (-a) * cmath.exp(-2j * math.pi * t * inv / mod) * mh5[(t * inv) % pm]
    worst = max(worst, abs(s - mh6[t]))
print(f"[A6] frequency recursion at n=6, eight t including non-units: max deviation {worst:.1e}")

# ----------------------------------------------------------------------------------------------------------
# [B] Fourier profile
# ----------------------------------------------------------------------------------------------------------
def levels_of(n: int) -> np.ndarray:
    mod = 3 ** n
    t = np.arange(mod, dtype=np.int64)
    v3 = np.zeros(mod, dtype=np.int64)
    m = t.copy()
    for _ in range(n):
        div = (m % 3 == 0) & (m > 0)
        v3[div] += 1
        m[div] //= 3
    return np.where(t == 0, 0, n - v3)


def profile(n: int):
    mod = 3 ** n
    amp = np.abs(np.fft.fft(laws[n]))
    lev = levels_of(n)
    t = np.arange(mod, dtype=np.int64)
    out = {}
    for h in range(1, n + 1):
        sel = np.flatnonzero(lev == h)
        a = amp[sel]
        k = int(np.argmax(a))
        out[h] = (float(a[k]), int(t[sel[k]]) // 3 ** (n - h), float((a ** 2).sum()), float((a ** 2).mean() * 3 ** h), len(sel))
    return amp, lev, out


profs = {}
for n in (10, 12, 14):
    amp, lev, out = profile(n)
    profs[n] = (amp, lev, out)
    print(f"[B] profile at n={n} {stamp()}")
_, _, P14 = profs[14]
print("[B1] h, M(h), argmax u (as 2^s), mass at level h, typical |mu_hat|^2 3^h  (n = 14, own law):")
for h in range(1, 15):
    M, u, mass, typ, cnt = P14[h]
    s = math.log2(u) if u > 0 else float("nan")
    print(f"      h={h:2d}: M={M:.6f}  u={u} (2^{s:.2f})  mass={mass:.5f}  typical={typ:.4f}  (#t={cnt})")
dev = max(abs(profs[n][2][h][0] - P14[h][0]) for n in (10, 12) for h in range(1, n + 1))


def upto_sign(u: int, h: int) -> int:
    """|mu_hat(-u)| = |mu_hat(u)|: the argmax is defined up to sign; normalise to min(u, 3^h - u)."""
    return min(u, 3 ** h - u)


same_arg = all(upto_sign(profs[n][2][h][1], h) == upto_sign(P14[h][1], h) for n in (10, 12) for h in range(1, n + 1))
print(f"[B2] Proposition 1 numerically: max |M_n(h) - M_14(h)| over n = 10, 12 = {dev:.1e}; same argmax u up to sign: {same_arg}")
print(f"[B3] M(1) = {P14[1][0]:.10f}; 1/sqrt(3) = {1 / math.sqrt(3):.10f}; |(1/3) e(1/3) + (2/3) e(2/3)|^2 = 1/9 + 4/9 + (4/9) cos(2 pi/3) = 1/3 exactly")
note_M = [0.577350, 0.377924, 0.252237, 0.176999, 0.129274, 0.096106, 0.075870, 0.060891, 0.048026, 0.038278, 0.031944, 0.026458, 0.022052, 0.019128, 0.016284, 0.014409, 0.012511, 0.011187]  # h = 15..18 from the session's fourier_deep / deep18 outputs
note_u = [1, 4, 8, 16, 32, 64, 256, 512, 1024, 4096, 8192, 16384, 65536, 131072, 262144, 1048576, 2097152, 8388608]
print("[B4] against the note's table (deep .out, h <= 14): max |M diff| = "
      f"{max(abs(P14[h][0] - note_M[h - 1]) for h in range(1, 15)):.1e}; argmax agree up to sign: {all(upto_sign(P14[h][1], h) == upto_sign(note_u[h - 1], h) for h in range(1, 15))}"
      f" (own argmax at h=9 is 18659 = -2^10 mod 3^9: the FFT breaks the |mu_hat(u)| = |mu_hat(-u)| tie arbitrarily)")
amp12 = profs[12][0]
print("[B5] |mu_hat_12(2^s)|, s = 0..40: " + ", ".join(f"{amp12[pow(2, s, 3 ** 12)]:.4f}" for s in range(41)))
# t/3^h of the argmax and the generic size
print("[B6] argmax frequency 2^s / 3^h: " + ", ".join(f"h={h}: {note_u[h - 1] / 3 ** h:.4f}" for h in range(7, 18)))
print("[B6b] s - h for the argmax u = 2^s (own to h = 14, the deep .outs beyond): "
      + ", ".join(f"h={h}: {round(math.log2(upto_sign(P14[h][1], h) if h <= 14 else note_u[h - 1])) - h:+d}" for h in range(2, 19))
      + "  -- the note's 's = h + 3 to h + 5' holds only for h = 13..18; s - h steps up by one every two or three levels")
print(f"[B7] generic coefficient size at level 17 from typical |mu_hat|^2 3^h = 0.708: sqrt(0.708) 3^-8.5 = {math.sqrt(0.708) * 3 ** -8.5:.2e}; "
      f"0.84 * 3^(-17/2) = {0.84 * 3 ** -8.5:.2e}; the note's '3e-4' would be {3e-4 / 3 ** -8.5:.2f} * 3^(-17/2)")
# Parseval: sum_t |mu_hat_n|^2 = 3^n sum mu^2 = (3/2) E_units[rho_n^2]; mass per level vs second-moment increment
print("[B8] Parseval and the second moment: n, sum|mu_hat|^2, 3^n sum mu^2, E_units[rho^2], (3/2)E, increment of E, (3/2) increment, mass at level n")
Eprev = None
for n in range(1, 15):
    mod = 3 ** n
    mu = laws[n]
    units = np.arange(mod) % 3 != 0
    rho = (2 / 3) * mod * mu
    E = float((rho[units] ** 2).mean())
    tot = float((np.abs(np.fft.fft(mu)) ** 2).sum()) if n <= 12 else float(mod * (mu ** 2).sum())
    inc = E - Eprev if Eprev is not None else float("nan")
    massn = P14[n][2]
    print(f"      n={n:2d}: {tot:.5f}  {mod * float((mu ** 2).sum()):.5f}  E={E:.5f}  (3/2)E={1.5 * E:.5f}  dE={inc:.4f}  (3/2)dE={1.5 * inc:.4f}  mass(n)={massn:.5f}")
    Eprev = E

# ----------------------------------------------------------------------------------------------------------
# [C] Proposition 2 numerics
# ----------------------------------------------------------------------------------------------------------
n = 14
amp14, lev14, _ = profs[14]
print("[C1] N=14: m, ||mu_N - lift mu_m||_1, max |lift^hat| at conductor > 3^m, max |mu_hat_N| at conductor > 3^m, M(m+1), ratio d/M(m+1)")
for m in range(1, 8):
    lift = np.tile(laws[m], 3 ** (n - m)) / 3 ** (n - m)
    d1 = float(np.abs(laws[n] - lift).sum())
    lh = np.abs(np.fft.fft(lift))
    above = lev14 > m
    print(f"      m={m}: d={d1:.4f}  lift primitive max={lh[above].max():.1e}  max|mu_hat| above={amp14[above].max():.4f}  M(m+1)={P14[m + 1][0]:.4f}  d/M={d1 / P14[m + 1][0]:.2f}")
# the inequality |mu_hat_h(u)| <= ||mu_h - lift mu_(h-1)||_1 at every level, every primitive u
worst = -1.0
for h in range(2, 13):
    lift = np.tile(laws[h - 1], 3) / 3
    d1 = float(np.abs(laws[h] - lift).sum())
    a = np.abs(np.fft.fft(laws[h]))
    prim = np.arange(3 ** h) % 3 != 0
    worst = max(worst, float(a[prim].max() / d1))
print(f"[C2] max over h <= 12 of M(h) / ||mu_h - lift mu_(h-1)||_1 = {worst:.4f} (must be <= 1)")
# normalisation: Mazur's ||rho_q - rho_m o pi||_q (mean over G_q) = (2/3) ||mu_q - lift mu_m||_1
q, m = 8, 3
rho_q = (2 / 3) * 3 ** q * laws[q]
rho_m = (2 / 3) * 3 ** m * laws[m]
lhs = float(np.abs(rho_q - np.tile(rho_m, 3 ** (q - m))).mean())
rhs = (2 / 3) * float(np.abs(laws[q] - np.tile(laws[m], 3 ** (q - m)) / 3 ** (q - m)).sum())
print(f"[C3] normalisation check q=8, m=3: mean_G_q |rho_q - rho_m o pi| = {lhs:.6f} = (2/3) ||mu_q - lift mu_m||_1 = {rhs:.6f}")

# ----------------------------------------------------------------------------------------------------------
# [D] Gauss sums
# ----------------------------------------------------------------------------------------------------------
def gauss_units(j: int, R: int = 60):
    mod, L = 3 ** j, 2 * 3 ** (j - 1)
    R = min(L, R)
    units = np.array([t for t in range(1, mod) if t % 3], dtype=np.int64)
    G = np.zeros(len(units), dtype=complex)
    for r in range(1, R + 1):
        inv = pow(2, -r, mod)
        G += (2.0 ** -r) * np.exp(2j * math.pi * ((units * inv) % mod) / mod)
    return units, G / (1 - 2.0 ** -L)


def gauss_at(j: int, t: int, R: int = 60) -> complex:
    mod, L = 3 ** j, 2 * 3 ** (j - 1)
    R = min(L, R)
    s = sum((2.0 ** -r) * cmath.exp(2j * math.pi * ((t * pow(2, -r, mod)) % mod) / mod) for r in range(1, R + 1))
    return s / (1 - 2.0 ** -L)


def gauss_2adic(j: int, t: int, R: int = 60) -> complex:
    mod, L = 3 ** j, 2 * 3 ** (j - 1)
    R = min(L, R)
    s = 0
    for r in range(1, R + 1):
        m_r = (-pow(3, -j, 2 ** r)) % (2 ** r)
        s += (2.0 ** -r) * cmath.exp(2j * math.pi * (t * m_r / 2 ** r + t / (2 ** r * mod)))
    return s / (1 - 2.0 ** -L)


print("[D1] j, sup_t |G_j| over units, argmax t (and t mod 3^j as +-2^s), mean |G_j|, share > 0.9, 2-adic reading deviation (10 random units)")
rng = random.Random(3)
for j in range(1, 13):
    units, G = gauss_units(j)
    aG = np.abs(G)
    k = int(np.argmax(aG))
    t = int(units[k])
    mod = 3 ** j
    pw = [s for s in range(0, 2 * 3 ** (j - 1)) if pow(2, s, mod) in (t, mod - t)]
    dev = max(abs(gauss_2adic(j, int(u)) - gauss_at(j, int(u))) for u in rng.sample(list(units), min(10, len(units))))
    print(f"      j={j:2d}: sup={aG.max():.5f} at t={t} (= +-2^{pw[:2]} mod 3^j)  mean={aG.mean():.4f}  share>0.9={float((aG > 0.9).mean()):.4f}  2-adic dev={dev:.1e}")
print("[D2] the stated inequality |G_j(2^s)| >= 1 - 2 pi 2^-K - 2^-s for 2^s <= 3^j / 2^K, tested for all j <= 12, s <= 40, K >= 0;")
print("     and the corrected one with tail 2^(1-s) (head: sum_(r<=s) 2^-r = 1 - 2^-s; tail |sum_(r>s)| <= 2^-s; both lose 2^-s)")
viol_stated, viol_corr, worst = 0, 0, (0.0, None)
by_s: dict[int, int] = {}
for j in range(1, 13):
    mod = 3 ** j
    for s in range(0, 41):
        if 2 ** s > mod:
            break
        g = abs(gauss_at(j, pow(2, s, mod)))
        K = 0
        while 2 ** s <= mod / 2 ** K:
            b_stated = 1 - 2 * math.pi * 2 ** -K - 2 ** -s
            b_corr = 1 - 2 * math.pi * 2 ** -K - 2 ** (1 - s)
            if g < b_stated - 1e-12:
                viol_stated += 1
                by_s[s] = by_s.get(s, 0) + 1
                if b_stated - g > worst[0]:
                    worst = (b_stated - g, (j, s, K, g, b_stated))
            if g < b_corr - 1e-12:
                viol_corr += 1
            K += 1
print(f"      violations of the stated inequality: {viol_stated} (by s: {dict(sorted(by_s.items()))}); of the corrected inequality: {viol_corr}; "
      f"worst stated violation (j, s, K, |G|, bound): {worst[1]}")
print("[D3] |G_j(2)| and |G_j(4)| for j = 4..12 (s = 1, 2: the stated bound allows 1 - 1/2 - eps and 1 - 1/4 - eps for large j): "
      + "; ".join(f"j={j}: {abs(gauss_at(j, 2)):.3f}, {abs(gauss_at(j, 4)):.3f}" for j in range(4, 13)))

# ----------------------------------------------------------------------------------------------------------
# [E] same-length collisions in the tree of 1
# ----------------------------------------------------------------------------------------------------------
DEPTH, NMAXC = 5, 12
MODBIG = 3 ** (NMAXC + DEPTH)  # residues carried mod 3^17, one power lost per level


def layers_mod(cap: int):
    res = np.array([1], dtype=np.int64)
    cost = np.array([0], dtype=np.int64)
    out = {}
    mod = MODBIG
    for d in range(1, DEPTH + 1):
        r3 = res % 3
        rs, cs = [], []
        for a in range(1, cap + 1):
            sel = r3 == (1 if a % 2 == 0 else 2)  # 2^a y = 1 mod 3
            if not sel.any():
                continue
            child = ((pow(2, a, mod) * res[sel] - 1) % mod) // 3  # exact: 3 | 2^a y - 1
            rs.append(child)
            cs.append(cost[sel] + a)
        res, cost = np.concatenate(rs), np.concatenate(cs)
        mod //= 3
        out[d] = (res, cost, mod)
    return out


def min_collision(res: np.ndarray, cost: np.ndarray, n: int):
    r = res % 3 ** n
    order = np.lexsort((cost, r))
    r, c = r[order], cost[order]
    eq = np.flatnonzero(r[1:] == r[:-1])
    if len(eq) == 0:
        return None, 0
    sums = c[eq] + c[eq + 1]
    k = int(np.argmin(sums))
    return int(sums[k]), int(r[eq[k]])


for cap in (24, 70):
    lay = layers_mod(cap)
    print(f"[E] tree of 1 to depth {DEPTH} with valuations <= {cap}: layer sizes {[len(lay[d][0]) for d in range(1, DEPTH + 1)]} {stamp()}")
    print(f"[E{'1' if cap == 24 else '2'}] cap {cap}: n | first colliding depth d, minimal A+A' there (class) | minimal colliding A+A' at each depth 1..{DEPTH} with the bound (n-d) log2 3 + 1")
    for n in range(4, NMAXC + 1):
        per = []
        first = None
        for d in range(1, DEPTH + 1):
            s, cls = min_collision(lay[d][0], lay[d][1], n)
            per.append(s)
            if s is not None and first is None:
                first = (d, s, cls)
        row = " ".join(f"d={d}: {'-' if per[d - 1] is None else per[d - 1]} (bound {(n - d) * LOG23 + 1:.1f})" for d in range(1, DEPTH + 1))
        print(f"      n={n:2d}: first {first} | {row}")
print("[E3] certification: a colliding pair with a valuation > cap has A + A' >= (cap + d) + 2d, so every minimum <= cap + 3d above is exact for that depth")
# how fast do the depth-d layers cover the units mod 3^n?  (the note: "mixing ... depth = n log 3 / log(4/3) = 3.8 n")
print("[E7] share of the 2 3^(n-1) unit classes mod 3^n hit by the depth-d layer (cap 70), d = 1..5; the note's 3.8 n against n log3/log4 = 0.79 n:")
for n in range(4, 9):
    cover = []
    for d in range(1, DEPTH + 1):
        r = np.unique(lay[d][0] % 3 ** n)
        cover.append(float((r % 3 != 0).sum() / (2 * 3 ** (n - 1))))
    print(f"      n={n}: " + ", ".join(f"d={d}: {c:.3f}" for d, c in enumerate(cover, 1)) + f"   (3.8 n = {n * math.log(3) / math.log(4 / 3):.1f}; n log3/log4 = {n * math.log(3) / math.log(4):.1f})")
# the proof's ingredients: C_w odd; C_w <= 2^A (3^d - 1)/2; C_w does not involve a_d; (d, C_w, A) determines w
def carry(word, u=3, v=2, inclusive=False):
    d, s, pref = len(word), 0, 0
    for j in range(1, d + 1):
        if inclusive:
            pref += word[j - 1]
        s += u ** (d - j) * v ** pref
        if not inclusive:
            pref += word[j - 1]
    return s


rng = random.Random(11)
odd_ok = bnd_ok = True
for _ in range(3000):
    d = rng.randint(1, 8)
    w = [rng.randint(1, 7) for _ in range(d)]
    C, A = carry(w), sum(w)
    odd_ok &= C % 2 == 1
    bnd_ok &= C <= 2 ** A * (3 ** d - 1) // 2
    Y1 = (C * pow(2, -A, 3 ** 12)) % 3 ** 12
    y = 0
    for a in w:
        y = (pow(2, -a, 3 ** 12) * (3 * y + 1)) % 3 ** 12
    bnd_ok &= Y1 == y
print(f"[E4] on 3000 random words: C_w odd {odd_ok}; C_w <= 2^A (3^d-1)/2 and Y(w) = C_w 2^-A mod 3^12: {bnd_ok}")
print(f"[E5] C_w does not involve the last valuation: C_(1,2) = {carry([1, 2])}, C_(1,3) = {carry([1, 3])}, C_(1,9) = {carry([1, 9])} -- (d, C_w) alone does NOT determine w; (d, C_w, A) does")


def recover(d: int, C: int, A: int) -> list[int]:
    w = []
    for _ in range(d - 1):
        C -= 3 ** (d - 1 - len(w))
        a = (C & -C).bit_length() - 1
        w.append(a)
        C >>= a
    w.append(A - sum(w))
    return w


rec_ok = all(recover(len(w), carry(w), sum(w)) == w for w in ([rng.randint(1, 7) for _ in range(rng.randint(1, 8))] for _ in range(2000)))
print(f"[E6] recovery a_1 = v_2(C_w - 3^(d-1)), recurse, a_d = A - (a_1 + .. + a_(d-1)) on 2000 random words: {rec_ok}")

# ----------------------------------------------------------------------------------------------------------
# [F] carry reciprocity, cycles
# ----------------------------------------------------------------------------------------------------------
import sympy as sp

u, v = sp.symbols("u v")
rng = random.Random(5)
sym_ok = True
for _ in range(120):
    d = rng.randint(1, 7)
    w = [rng.randint(1, 5) for _ in range(d)]
    A = sum(w)
    lhs = carry(w[::-1], u, v)
    rhs = sp.expand(u ** (d - 1) * v ** A * carry(w, 1 / u, 1 / v, inclusive=True))
    sym_ok &= sp.expand(lhs - rhs) == 0
    # and the exclusive carry is not reciprocal to the reversed exclusive carry (generic words)
print(f"[F1] C_rev(w)(u,v) = u^(d-1) v^A C'_w(1/u,1/v) as a POLYNOMIAL IDENTITY (sympy, 120 random words): {sym_ok}")
w = [1, 2]
print(f"[F2] w=(1,2): C_w={carry(w)}, C'_w={carry(w, inclusive=True)}, C_rev={carry(w[::-1])}, 3*8*C'_w(1/3,1/2)={3 * 8 * carry(w, Fraction(1, 3), Fraction(1, 2), inclusive=True)}")


def syr_word(x: int, k: int, c: int = 1) -> list[int]:
    """valuations of k Syracuse steps x -> (3x + c)/2^v on integers."""
    word = []
    for _ in range(k):
        y = 3 * x + c
        a = (y & -y).bit_length() - 1
        word.append(a)
        x = y >> a
    return word


for name, w, y0 in (("-1", [1], -1), ("{-5,-7}", [1, 2], -5), ("7-cycle of -17", [1, 1, 1, 2, 1, 1, 4], -17)):
    k, A = len(w), sum(w)
    C, Cr = carry(w), carry(w[::-1])
    den = 2 ** A - 3 ** k
    orbit_word = syr_word(y0, k)
    x = y0
    orb = [x]
    for a in orbit_word:
        x = (3 * x + 1) >> a
        orb.append(x)
    rot = any(w[i:] + w[:i] == w[::-1] for i in range(k))
    print(f"[F3] {name}: orbit {orb} has word {orbit_word} (note's word {w}: {orbit_word == w}); C_w={C}, y0=C_w/(2^A-3^k)={Fraction(C, den)}; "
          f"C_rev={Cr}, y0'={Fraction(Cr, den)} ({'integer' if Cr % den == 0 else 'not an integer'}); reversed word is a rotation: {rot}")
# -13801 on the integer cycle of x -> (3x + 139)/2^v with the reversed word
x = -13801
wd = syr_word(x, 7, 139)
xx = x
for a in wd:
    xx = (3 * xx + 139) >> a
print(f"[F4] x -> (3x+139)/2^v from -13801: word {wd} (reversed seven-word (4,1,1,2,1,1,1): {wd == [4, 1, 1, 2, 1, 1, 1]}), returns to {xx} (cycle: {xx == x}); 139 = 3^7 - 2^11 = {3 ** 7 - 2 ** 11}; gcd(13801,139) = {math.gcd(13801, 139)}")
# the 2-adic source class: x = (2^A y - C_w)/3^d for odd y; the note says x = -C_w 3^-d mod 2^(A+1)
rng = random.Random(9)
okA = okA1 = True
for _ in range(2000):
    d = rng.randint(1, 6)
    w = [rng.randint(1, 5) for _ in range(d)]
    A, C = sum(w), carry(w)
    y = 2 * rng.randint(-10 ** 6, 10 ** 6) + 1
    num = 2 ** A * y - C
    if num % 3 ** d:  # need y = C_w 2^-A mod 3^d for an integer source
        y += 2 * (((C * pow(2, -A, 3 ** d) - y) * pow(2, -1, 3 ** d)) % 3 ** d)
        num = 2 ** A * y - C
    x = num // 3 ** d
    okA &= (x - (-C * pow(3, -d, 2 ** A))) % 2 ** A == 0
    okA1 &= (x - (-C * pow(3, -d, 2 ** (A + 1)))) % 2 ** (A + 1) == 0
print(f"[F5] x = -C_w 3^-d mod 2^A on 2000 sources: {okA}; mod 2^(A+1) as the note states: {okA1} (the class mod 2^(A+1) is (2^A - C_w) 3^-d)")

# ----------------------------------------------------------------------------------------------------------
# [G] resonance at h = 10 by cost and by last valuation; LD rate
# ----------------------------------------------------------------------------------------------------------
h, AMAX = 10, 50
joint = np.zeros((1, AMAX + 1))
joint[0, 0] = 1.0
for lev in range(1, h + 1):
    m = 3 ** lev
    new = np.zeros((m, AMAX + 1))
    base = (3 * np.arange(joint.shape[0]) + 1) % m
    for a in range(1, AMAX + 1):
        z = (pow(2, -a, m) * base) % m
        new[z, a:] += (2.0 ** -a) * joint[:, : AMAX + 1 - a]
    joint = new
mu10 = joint.sum(axis=1)
mod = 3 ** h
t_star = pow(2, 12, mod)
ph = np.exp(-2j * math.pi * t_star * np.arange(mod) / mod)
contrib = ph @ joint
tot = contrib.sum()
amp10 = np.abs(np.fft.fft(mu10))
print(f"[G1] h=10: |mu_hat(2^12)| = {abs(tot):.5f} (FFT {amp10[t_star]:.5f}); generic t=7: {amp10[7]:.5f}; own law vs section-A law: {np.abs(mu10 - laws[10]).max():.1e}")
order = np.argsort(-np.abs(contrib))[:10]
print("[G2] by total cost A: " + ", ".join(f"A={A}: {abs(contrib[A]):.4f} (P(A)={joint[:, A].sum():.4f})" for A in order))
print(f"[G3] sum_A |.| = {np.abs(contrib).sum():.4f}; coherent fraction = {abs(tot) / np.abs(contrib).sum():.3f}; |.|-weighted mean A/h = {(np.abs(contrib) * np.arange(AMAX + 1)).sum() / np.abs(contrib).sum() / h:.3f}")
pm = 3 ** (h - 1)
mu9 = laws[9]
base = (3 * np.arange(pm) + 1) % mod
parts = []
for a in range(1, 30):
    z = (pow(2, -a, mod) * base) % mod
    parts.append((2.0 ** -a) * (np.exp(-2j * math.pi * t_star * z / mod) @ mu9))
print("[G4] by last valuation a: " + ", ".join(f"a={a + 1}: {abs(c):.4f}" for a, c in enumerate(parts[:6])) + f"; sum = {abs(sum(parts)):.5f}")
# exact P(A = k) for h geometric(1/2) parts: C(k-1, h-1) 2^-k
print("[G5] check P(A) at h=10: " + ", ".join(f"A={k}: {math.comb(k - 1, 9) / 2 ** k:.4f}" for k in (12, 13, 14, 15, 16, 17)))


def rate(alpha: float) -> float:
    return alpha * math.log(2) + (alpha - 1) * math.log(alpha - 1) - alpha * math.log(alpha)


print("[G6] Cramer rate of A/h for i.i.d. geometric(1/2) valuations, I(a) = a ln2 + (a-1) ln(a-1) - a ln a: "
      + ", ".join(f"I({a:.2f})={rate(a):.4f} (e^-I={math.exp(-rate(a)):.3f})" for a in (1.2, 1.3, 1.35, 1.4, 1.48, 1.5, 1.6, 1.7)))
print("[G7] exact -ln P(A = round(1.35 h))/h at h = 10, 14, 20, 40: "
      + ", ".join(f"h={H}: {-math.log(math.comb(round(1.35 * H) - 1, H - 1) / 2 ** round(1.35 * H)) / H:.4f}" for H in (10, 14, 20, 40)))

# ----------------------------------------------------------------------------------------------------------
# [F6] the S20 addendum's Proposition 6 (trace of a rational cycle is reversal-invariant): quick check
# ----------------------------------------------------------------------------------------------------------
def cycle_sum(word, u, v):
    return sum(carry(word[j:] + word[:j], u, v) for j in range(len(word)))


rng = random.Random(13)
ok6 = all(cycle_sum(w, u, v) == cycle_sum(w[::-1], u, v)
          for w, u, v in (([rng.randint(1, 6) for _ in range(rng.randint(1, 8))], Fraction(rng.randint(1, 7), rng.randint(1, 7)), Fraction(rng.randint(1, 7), rng.randint(1, 7))) for _ in range(300)))
w7 = [1, 1, 1, 2, 1, 1, 4]
tr = sum(Fraction(carry(w7[j:] + w7[:j]), 2 ** 11 - 3 ** 7) for j in range(7))
trr = sum(Fraction(carry(w7[::-1][j:] + w7[::-1][:j]), 2 ** 11 - 3 ** 7) for j in range(7))
print(f"[F6] Proposition 6 (added to the note during the audit): sum_j C_(rot_j w)(u,v) reversal-invariant on 300 random words/rational (u,v): {ok6}; "
      f"seven-cycle trace {tr} = reversed-word trace {trr} (mechanism: the coefficient of u^(d-i) is the sum of v^(block sum) over cyclic blocks of length i-1, a multiset reversal preserves)")

# ----------------------------------------------------------------------------------------------------------
# [H] decay models on M(h), h <= 18 (own values to 14, the deep .outs beyond; they agree to 1e-6)
# ----------------------------------------------------------------------------------------------------------
M = np.array([P14[h][0] for h in range(1, 15)] + note_M[14:])
hs = np.arange(1, 19)
ratios = M[1:] / M[:-1]
print("[H1] ratios M(h)/M(h-1), h=2..18: " + ", ".join(f"{r:.3f}" for r in ratios))
print(f"[H2] window means: h=6..9: {ratios[4:8].mean():.4f}; h=10..13: {ratios[8:12].mean():.4f}; h=14..17: {ratios[12:16].mean():.4f}; h=14..18: {ratios[12:17].mean():.4f}; "
      f"least-squares slope of the ratio over h=11..18: {np.polyfit(np.arange(11, 19), ratios[9:17], 1)[0]:+.4f} per level; the ratio at h=18 ({ratios[16]:.3f}) is the largest so far")
print("[H3] local exponents log2(M(h)/M(2h)): " + ", ".join(f"{h}->{2 * h}: {math.log2(M[h - 1] / M[2 * h - 1]):.3f}" for h in (5, 6, 7, 8, 9)))


def fit_report(name, sel, model_log, params):
    pred = np.exp(model_log(hs[sel], *params))
    rel = pred / M[sel] - 1
    return f"{name}: max rel. residual {np.abs(rel).max() * 100:.1f}%, rms {np.sqrt((rel ** 2).mean()) * 100:.1f}%"


sel = (hs >= 7)
# (i) pure power law C h^-alpha on 7..18
a1, b1 = np.polyfit(np.log(hs[sel]), np.log(M[sel]), 1)
print("[H4] " + fit_report(f"pure power C h^-alpha on 7..18 (alpha = {-a1:.3f})", sel, lambda x, a, b: a * np.log(x) + b, (a1, b1)))
# (ii) shifted power law C (h + c)^-alpha, c on a grid
best = None
for c in np.arange(0, 15.01, 0.25):
    a, b = np.polyfit(np.log(hs[sel] + c), np.log(M[sel]), 1)
    rel = np.exp(a * np.log(hs[sel] + c) + b) / M[sel] - 1
    r = np.sqrt((rel ** 2).mean())
    if best is None or r < best[0]:
        best = (r, c, -a, b)
r, c_s, al_s, b_s = best
print("[H5] " + fit_report(f"shifted power C (h + c)^-alpha on 7..18 (c = {c_s:.2f}, alpha = {al_s:.3f})", sel, lambda x, a, bb: a * np.log(x + c_s) + bb, (-al_s, b_s))
      + f"; its ratios at h=14..18: {', '.join(f'{((hh - 1 + c_s) / (hh + c_s)) ** al_s:.3f}' for hh in (14, 15, 16, 17, 18))}; local exponents 5->10..9->18: "
      + ", ".join(f"{al_s * math.log2((2 * hh + c_s) / (hh + c_s)):.2f}" for hh in (5, 6, 7, 8, 9)))
# the same fits on 7..17 / 11..17 only, and their prediction of the new level 18
best17 = None
s17 = (hs >= 7) & (hs <= 17)
for c in np.arange(0, 15.01, 0.25):
    a, b = np.polyfit(np.log(hs[s17] + c), np.log(M[s17]), 1)
    rel = np.exp(a * np.log(hs[s17] + c) + b) / M[s17] - 1
    r = np.sqrt((rel ** 2).mean())
    if best17 is None or r < best17[0]:
        best17 = (r, c, -a, b)
_, c17, al17, b17 = best17
g17 = (hs >= 11) & (hs <= 17)
a11, b11 = np.polyfit(hs[g17], np.log(M[g17]), 1)
print(f"[H5b] fitted on 7..17 only: shifted power (c = {c17:.2f}, alpha = {al17:.3f}) predicts M(18) = {math.exp(-al17 * math.log(18 + c17) + b17):.5f} (ratio {((17 + c17) / (18 + c17)) ** al17:.3f}); "
      f"geometric on 11..17 (r = {math.exp(a11):.4f}) predicts M(18) = {math.exp(a11 * 18 + b11):.5f} (ratio {math.exp(a11):.3f}); measured M(18) = {M[17]:.5f} (ratio {ratios[16]:.3f})")
# (iii) geometric C r^h on 11..18 and on 7..18
for lo in (11, 7):
    s2 = hs >= lo
    a, b = np.polyfit(hs[s2], np.log(M[s2]), 1)
    print("[H6] " + fit_report(f"geometric C r^h on {lo}..18 (r = {math.exp(a):.4f})", s2, lambda x, aa, bb: aa * x + bb, (a, b)))
# (iv) geometric with power prefactor C h^-beta r^h on 7..18
X = np.column_stack([hs[sel], np.log(hs[sel]), np.ones(sel.sum())])
coef, *_ = np.linalg.lstsq(X, np.log(M[sel]), rcond=None)
print("[H7] " + fit_report(f"C h^-beta r^h on 7..18 (r = {math.exp(coef[0]):.4f}, beta = {-coef[1]:.3f})", sel, lambda x, a, bb, cc: a * x + bb * np.log(x) + cc, tuple(coef)))
# the note's 3.75 h^-2
rel = 3.75 * hs ** -2.0 / M - 1
print("[H8] 3.75 h^-2: relative error h=7..14: " + ", ".join(f"{rel[h - 1] * 100:+.1f}%" for h in range(7, 15)) + "; h=15..18: " + ", ".join(f"{rel[h - 1] * 100:+.1f}%" for h in range(15, 19)))
# extrapolations: what the two readings predict
ag, bg = np.polyfit(hs[hs >= 11], np.log(M[hs >= 11]), 1)
print("[H9] extrapolation of the two readings: h=25, 34, 50: shifted power (7..18) "
      + ", ".join(f"{math.exp(-al_s * math.log(H + c_s) + b_s):.2e}" for H in (25, 34, 50)) + "; geometric (11..18) "
      + ", ".join(f"{math.exp(ag * H + bg):.2e}" for H in (25, 34, 50))
      + "; the two readings first differ by a factor 2 at h = "
      + str(next(H for H in range(19, 300) if math.exp(ag * H + bg) < 0.5 * math.exp(-al_s * math.log(H + c_s) + b_s))))
# ratio bound vs ratio: rescaled ratio to see the trend against 1
print("[H10] 1 - ratio, h=10..18: " + ", ".join(f"{1 - r:.3f}" for r in ratios[8:17]) + "  (a shifted power law has 1 - ratio ~ alpha/(h + c) -> 0; geometric has 1 - r constant)")
print(f"DONE {stamp()}")
