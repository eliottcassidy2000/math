#!/usr/bin/env python3
"""Five structural mirrors for the 3-adic Syracuse law (S20, 2026-09-29).

(1) Fourier side (the circulant / Lambert mirror): with mu_n the level-n law of Tao's Syracuse random
    variable on Z/3^n and mu_hat_n(t) = E e(t Y_n / 3^n), the exact one-step recursion
        mu_hat_n(t) = sum_(a>=1) 2^-a e(t 2^-a / 3^n) mu_hat_(n-1)(t 2^-a mod 3^(n-1))
    is a twisted geometric random walk on the cyclic unit group in the frequency variable: the
    multiplication by 2^-a is the circulant step, the phase e(t 2^-a/3^n) the twist.  We compute the
    conductor profile M_n(h) = max_(cond t = 3^h) |mu_hat_n(t)| (Mazur's "primitive Fourier coefficient
    bound at level h") and the Fourier mass above conductor 3^m against the ell^1 fine-scale distance.
(2) The one-step geometric Gauss sum G_j(t) = c_j sum_(r=1)^(L_j) 2^-r e(t 2^-r / 3^j) (L_j = 2 3^(j-1),
    c_j = 1/(1 - 2^-L_j)), its 2-adic reading e(t m_r / 2^r + t/(2^r 3^j)) with m_r = -3^-j mod 2^r
    (the reciprocal-prime evaluation), and sup_t |G_j(t)| over units.
(3) Same-length collisions (the Moore-bound / tessellation mirror): words of the same length d and
    costs A, A' with A + A' <= (n - d) log_2 3 + 1 land on distinct classes mod 3^n (PROVED); the actual
    first collision of the depth-d layer of the tree of 1.
(4) Carry-polynomial reciprocity (the reciprocal-polynomial mirror): with the exclusive-prefix carry
    C_w(u, v) = sum_j u^(d-j) v^(a_1+...+a_(j-1)) (the affine identity 3^d x + C_w = 2^A y) and the inclusive-prefix
    carry C'_w(u, v) = sum_j u^(d-j) v^(a_1+...+a_j), word reversal is polynomial reciprocity:
    C_rev(w)(u, v) = u^(d-1) v^A C'_w(1/u, 1/v).  The reversed word of the seven-cycle of -17 is a rational, not an
    integer, cycle.
Run: python 04-computation/experiments/collatz_five_mirrors_20260929.py   (about a minute, < 2 GB)
"""
from __future__ import annotations

import cmath
import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from mazur_harmonic_mass_deep_20260928 import level_up  # noqa: E402

LOG23 = math.log(3, 2)


def conductor_level(t: int, n: int) -> int:
    """h with 3^(n-h) || t, i.e. t = 3^(n-h) u, 3 !| u; t = 0 has level 0."""
    if t == 0:
        return 0
    h = n
    while t % 3 == 0:
        t //= 3
        h -= 1
    return h


def carry(word, u=3, v=2, inclusive=False):
    """C_w(u, v) = sum_j u^(d-j) v^(a_1+...+a_(j-1)); inclusive=True uses a_1+...+a_j."""
    d = len(word)
    s, pref = 0, 0
    for j in range(1, d + 1):
        if inclusive:
            pref += word[j - 1]
        s += u ** (d - j) * v ** pref
        if not inclusive:
            pref += word[j - 1]
    return s


if __name__ == "__main__":
    NMAX = 14
    mu = np.array([0.0, 1 / 3, 2 / 3])
    laws = {1: mu}
    for n in range(2, NMAX + 1):
        mu = level_up(mu, n)
        laws[n] = mu
    print("== (1) Fourier profile of the 3-adic Syracuse law: M_n(h) = max_(cond t = 3^h) |mu_hat_n(t)| ==")
    for n in (4, 6, 8, 10, 12, 14):
        mod = 3 ** n
        mh = np.fft.fft(laws[n])  # mu_hat(t) = sum_y mu(y) e^(-2 pi i t y / 3^n); modulus is what we need
        amp = np.abs(mh)
        t = np.arange(mod)
        lev = np.zeros(mod, dtype=np.int64)
        tt = t.copy()
        # level = n - v_3(t) for t > 0
        v3 = np.zeros(mod, dtype=np.int64)
        m = tt.copy()
        for _ in range(n):
            div = (m % 3 == 0) & (m > 0)
            v3[div] += 1
            m[div] //= 3
        lev = np.where(t == 0, 0, n - v3)
        prof = []
        mass_above = []
        for h in range(1, n + 1):
            sel = lev == h
            prof.append((h, float(amp[sel].max()), int(t[sel][np.argmax(amp[sel])]), float((amp[sel] ** 2).sum())))
        # fine-scale ell^1 distance to the lift of mu_m, and the Fourier mass above conductor 3^m (Parseval: sum |mu_hat|^2 over cond > 3^m = 3^n ||mu_n - lift mu_m||_2^2)
        line = "; ".join(f"h={h}: max {mx:.4f} at t={am} (mass {ms:.4f})" for h, mx, am, ms in prof)
        print(f"   n={n}: {line}")
        if n == 14:
            print("   profile M(h) at n=14: " + ", ".join(f"{mx:.5f}" for h, mx, am, ms in prof))
            print("   ratios M(h)/M(h-1): " + ", ".join(f"{prof[i][1]/prof[i-1][1]:.4f}" for i in range(1, len(prof))))
            print("   argmax t divided by 3^(n-h): " + ", ".join(f"h={h}: {am // 3 ** (n - h)}" for h, mx, am, ms in prof))
        if n == 12:
            print("   |mu_hat_12(2^s)| for s = 0..40: " + ", ".join(f"{amp[pow(2, s, mod)]:.4f}" for s in range(0, 41)))
        for mm in (1, 2, n // 2):
            lift = np.tile(laws[mm], 3 ** (n - mm)) / 3 ** (n - mm)
            l1 = float(np.abs(laws[n] - lift).sum())
            fmass = sum(ms for h, mx, am, ms in prof if h > mm)
            l2sq = float(((laws[n] - lift) ** 2).sum()) * mod
            print(f"      m={mm}: ||mu_n - lift mu_m||_1 = {l1:.4f}; Fourier mass above conductor 3^m = {fmass:.4f} = 3^n ||.||_2^2 = {l2sq:.4f}")
    # the exact one-step recursion at n = 6 for a few t
    n = 6
    mod, pm = 3 ** n, 3 ** (n - 1)
    mh6 = np.fft.fft(laws[6])
    mh5 = np.fft.fft(laws[5])
    L = 2 * 3 ** (n - 1)
    ok = True
    for tval in (1, 2, 7, 100, 728):
        s = 0
        for a in range(1, L + 1):
            inv = pow(2, -a, mod)
            s += 2.0 ** (-a) / (1 - 2.0 ** (-L)) * cmath.exp(-2j * math.pi * tval * inv / mod) * mh5[(tval * inv) % pm]
        ok &= abs(s - mh6[tval]) < 1e-9
    print(f"   one-step recursion mu_hat_n(t) = sum_a 2^-a e(t 2^-a/3^n) mu_hat_(n-1)(t 2^-a mod 3^(n-1)) verified at n=6: {ok}")

    print("== (2) one-step geometric Gauss sums G_j(t) = c_j sum_r 2^-r e(t 2^-r/3^j) over units t ==")
    for j in range(1, 13):
        mod = 3 ** j
        L = 2 * 3 ** (j - 1)
        R = min(L, 60)  # terms beyond r = 60 are below 1e-18
        invs = np.array([pow(2, -r, mod) for r in range(1, R + 1)], dtype=np.int64)
        w = 2.0 ** (-np.arange(1, R + 1)) / (1 - 2.0 ** (-L))
        units = np.array([t for t in range(1, mod) if t % 3], dtype=np.int64)
        # G(t) = sum_r w_r e(t inv_r / mod)
        phases = np.exp(2j * math.pi * (np.outer(units, invs) % mod) / mod)
        G = phases @ w
        absG = np.abs(G)
        k = int(np.argmax(absG))
        # 2-adic reading for t = units[k]: e(t m_r/2^r + t/(2^r 3^j)), m_r = -3^-j mod 2^r
        t = int(units[k])
        s2 = 0
        for r in range(1, R + 1):
            m_r = (-pow(3, -j, 2 ** r)) % (2 ** r)
            s2 += w[r - 1] * cmath.exp(2j * math.pi * (t * m_r / 2 ** r + t / (2 ** r * mod)))
        print(f"   j={j:2d}: sup_t |G_j| = {absG.max():.5f} at t = {t} (t/3^j = {t/mod:.4f}; 2-adic reading agrees: {abs(s2 - G[k]) < 1e-9}); mean |G_j| = {absG.mean():.4f}; share of units with |G_j| > 0.9: {(absG > 0.9).mean():.4f}")
    print("   (the supremum tends to 1: for t = 2^s, s large, t 2^-r is a small integer for r <= s and the phases align; no uniform one-step gap)")

    print("== (3) same-length collisions mod 3^n in the tree of 1 ==")
    for n in range(4, 13):
        mod = 3 ** n
        # depth-d layer with valuations <= 24: classes and costs; first depth with two nodes in one class
        layer = [(1, 0)]
        found = None
        for d in range(1, 14):
            nxt = []
            for y, A in layer:
                r = y % 3
                if r == 0:
                    continue
                a = 2 - (1 if r == 2 else 0)
                while a <= 24:
                    nxt.append((((1 << a) * y - 1) // 3, A + a))
                    a += 2
            layer = nxt
            seen = {}
            best = None
            for y, A in layer:
                c = y % mod
                if c in seen:
                    tot = seen[c] + A
                    if best is None or tot < best[0]:
                        best = (tot, d, c)
                else:
                    seen[c] = A
            if best is not None:
                found = best
                break
        bound = lambda d: (n - d) * LOG23 + 1  # noqa: E731
        if found:
            tot, d, c = found
            print(f"   n={n:2d}: first same-length collision at depth d={d}, minimal cost sum A+A' = {tot} (class {c}); proved no-collision bound at that depth: A + A' <= {bound(d):.2f}")
        else:
            print(f"   n={n:2d}: no same-length collision to depth 13 with valuations <= 24")

    print("== (4) carry-polynomial reciprocity and reversed cycles ==")
    import random
    rng = random.Random(1)
    ok = True
    for _ in range(200):
        d = rng.randint(1, 8)
        w = [rng.randint(1, 5) for _ in range(d)]
        A = sum(w)
        u, v = Fraction(3), Fraction(2)
        lhs = Fraction(carry(w[::-1]))
        rhs = u ** (d - 1) * v ** A * Fraction(carry(w, Fraction(1, 3), Fraction(1, 2), inclusive=True))
        ok &= lhs == rhs
        # and the exclusive carry is NOT reciprocal to the reversed one (the shift by one valuation)
    print(f"   C_rev(w)(3,2) = 3^(d-1) 2^A C'_w(1/3, 1/2) (inclusive-prefix carry) on 200 random words: {ok}")
    print(f"   example w = (1,2): C_w = {carry([1, 2])}, C'_w = {carry([1, 2], inclusive=True)}, C_rev = {carry([2, 1])} = 3 * 8 * C'_w(1/3, 1/2) = {3 * 8 * carry([1, 2], Fraction(1, 3), Fraction(1, 2), inclusive=True)}")
    for name, w in (("-1", [1]), ("{-5,-7}", [1, 2]), ("7-cycle of -17", [1, 1, 1, 2, 1, 1, 4])):
        k, A = len(w), sum(w)
        C, Cr = carry(w), carry(w[::-1])
        den = 2 ** A - 3 ** k
        print(f"   {name}: word {w}, C_w = {C}, y_0 = C_w/(2^A - 3^k) = {Fraction(C, den)}; reversed word {w[::-1]}, C_rev = {Cr}, y_0' = {Fraction(Cr, den)} "
              f"({'integer' if Cr % den == 0 else 'rational, not integer'}; rotation of the original: {any(w[i:] + w[:i] == w[::-1] for i in range(k))})")
    print("DONE")
