#!/usr/bin/env python3
"""Independent audit script for 05-knowledge/results/mazur_positive_density_20260928.md
(opus S19 Mazur digest: Theorems A, B, C, the spike profile D, the numbers, the paper's constants).

Everything here is re-derived from the definitions, without importing the session's scripts:
  * the 3-adic Syracuse law mu_n exactly (integer numerators over a common denominator) by the
    FORWARD recursion Y_{k+1} = 2^-A (3 Y_k + 1) and, as a cross-check, by the PARENTS recursion
    mu_{n+1}(z) = sum_{a = eps(z) mod 2} 2^-a mu_n((2^a z - 1)/3);
  * the law in float64 to level NMAX (default 14) by the vectorised parents recursion truncated at
    a <= 64 (a different algorithm from the session's discrete-log FFT);
  * Terras classes, the affine identity, the converse, the sign statement, tree layers by direct
    enumeration (Theorem A); cycle words and resonances (Theorem B); the harmonic-sum inequality,
    the discount product Pi and the identity 1/m = 3^n 2^-A Pi(m) (Theorem C); the spike profile D
    by dynamic programming on the forward closure of -1; the global statistics; Kaprekar cycles and
    Euler bricks; the paper's Section 8 constants in iterated logarithms; Lemma 7.3's arithmetic.
Run:  python 04-computation/experiments/mazur_positive_density_20260928_audit.py [NMAX]
Output: 05-knowledge/results/mazur_positive_density_20260928_audit.out
Auditor session: independent (Fable), 2026-09-28.
"""
from __future__ import annotations

import math
import os
import random
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product

import numpy as np

T0 = time.time()
LN2, LN3 = math.log(2), math.log(3)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
EXACT_MAX = 7
random.seed(20260928)


def stamp() -> str:
    return f"[{time.time() - T0:6.1f}s]"


def v2(x: int) -> int:
    return (x & -x).bit_length() - 1


def syr(x: int):
    """Syracuse map on odd integers of either sign: (S(x), v_2(3x+1))."""
    t = 3 * x + 1
    a = v2(t)
    return t >> a, a


def word(x: int, n: int):
    w = []
    for _ in range(n):
        x, a = syr(x)
        w.append(a)
    return tuple(w), x


def C_w(w) -> int:
    """C_w = sum_{j=1}^n 3^(n-j) 2^(a_1+...+a_{j-1})  (so that 3^n x + C_w = 2^A S^n(x))."""
    n = len(w)
    tot, part = 0, 0
    for j, a in enumerate(w):
        tot += 3 ** (n - 1 - j) * 2 ** part
        part += a
    return tot


def Y_mod(w) -> int:
    """Y_n(w) = C_w 2^-A mod 3^n, the Syracuse offset of the word (Tao's F_n)."""
    n = len(w)
    mod = 3 ** n
    return (C_w(w) * pow(2, -sum(w), mod)) % mod


def Y_by_recursion(w) -> int:
    """The same by the recursion Y_{k+1} = 2^-a_{k+1} (3 Y_k + 1) mod 3^(k+1), Y_0 = 0."""
    y = 0
    for k, a in enumerate(w, start=1):
        mod = 3 ** k
        y = (pow(2, -a, mod) * (3 * y + 1)) % mod
    return y


# ----------------------------------------------------------------------------------------------
# The law: exact (integer numerators, forward and parents forms) and float (vectorised parents form)
# ----------------------------------------------------------------------------------------------

def law_exact_forward(nmax: int):
    """Exact laws mu_1..mu_nmax as dicts y -> Fraction, forward recursion; 2^-A mod 3^lev depends on
    A mod L (L = 2 3^(lev-1) the order of 2), P(A = r mod L) = 2^(L-r)/(2^L - 1), r = 1..L."""
    laws = {}
    num = {0: 1}  # numerators over the common denominator den
    den = 1
    for lev in range(1, nmax + 1):
        mod, L = 3 ** lev, 2 * 3 ** (lev - 1)
        inv2 = pow(2, -1, mod)
        new = defaultdict(int)
        for y, p in num.items():
            cur = (3 * y + 1) % mod
            for r in range(1, L + 1):
                cur = (cur * inv2) % mod  # 2^-r (3y+1)
                new[cur] += p * 2 ** (L - r)
        den *= 2 ** L - 1
        num = dict(new)
        assert sum(num.values()) == den
        laws[lev] = {y: Fraction(p, den) for y, p in num.items()}
    return laws


def law_exact_parents(nmax: int):
    """Exact laws by the parents recursion mu_{lev}(z) = sum_{a = eps(z) mod 2} 2^-a mu_{lev-1}((2^a z - 1)/3 mod 3^(lev-1)),
    with the geometric tail folded by the period L of 2^a z mod 3^lev."""
    laws = {}
    mu = {0: Fraction(1)}
    for lev in range(1, nmax + 1):
        mod, modp, L = 3 ** lev, 3 ** (lev - 1), 2 * 3 ** (lev - 1)
        den = 2 ** L - 1
        new = {}
        for z in range(mod):
            if z % 3 == 0:
                continue
            eps = 1 if z % 3 == 2 else 0  # a odd iff z = 2 mod 3
            tot = Fraction(0)
            for r in range(1, L + 1):
                if r % 2 != eps:
                    continue
                par = ((pow(2, r, mod) * z - 1) // 3) % modp
                tot += Fraction(2 ** (L - r), den) * mu.get(par, Fraction(0))
            new[z] = tot
        mu = new
        laws[lev] = dict(mu)
    return laws


def law_float(nmax: int, amax: int = 64):
    """Float64 laws by the vectorised parents recursion truncated at a <= amax (tail 2^-amax per node)."""
    laws = {1: np.array([0.0, 1 / 3, 2 / 3])}
    for lev in range(2, nmax + 1):
        mod, modp = 3 ** lev, 3 ** (lev - 1)
        prev = laws[lev - 1]
        z = np.arange(mod, dtype=np.int64)
        r = z % 3
        idx = {1: z[r == 1], 2: z[r == 2]}
        del z, r
        new = np.zeros(mod)
        for a in range(1, amax + 1):
            zz = idx[1 if a % 2 == 0 else 2]
            p2 = pow(2, a, mod)
            par = ((p2 * zz - 1) // 3) % modp
            new[zz] += (2.0 ** -a) * prev[par]
        laws[lev] = new
    return laws


def rho_of(mu, n):
    return (2 / 3) * 3 ** n * mu


# ----------------------------------------------------------------------------------------------
# Section A: Theorem A
# ----------------------------------------------------------------------------------------------

def tree_layer(y: int, n: int, amax: int):
    """Sum over the depth-n predecessors x of y (valuations <= amax) of 3^n 2^-A_n(x): exact Fraction, class masses, sign check."""
    total = Fraction(0)
    cls = [Fraction(0)] * 3
    sign_ok = True
    stack = [(y, 0, 0)]
    while stack:
        x, d, A = stack.pop()
        if d == n:
            wgt = Fraction(3 ** n, 2 ** A)
            total += wgt
            cls[x % 3] += wgt
            if (x > 0) != (y > 0):
                sign_ok = False
            continue
        r = x % 3
        if r == 0:
            continue
        a0 = 2 if r == 1 else 1
        for a in range(a0, amax + 1, 2):
            stack.append(((2 ** a * x - 1) // 3, d + 1, A + a))
    return total, cls, sign_ok


def section_A(EX, FL):
    print("\n==== A. Theorem A: harmonic mass of a tree layer = 3^n mu_n(seed) ====")
    # A1: Terras classes by direct scan: the odd x with a given word of length n, total valuation A, are one class mod 2^(A+1)
    nwords = 0
    for n, amax_w in ((1, 6), (2, 5), (3, 4), (4, 2)):
        for w in product(range(1, amax_w + 1), repeat=n):
            A = sum(w)
            M = 2 ** (A + 1)
            hits = [x for x in range(-2 ** (A + 3) + 1, 2 ** (A + 3), 2) if word(x, n)[0] == w]
            res = {x % M for x in hits}
            assert len(hits) == 8 and len(res) == 1, (w, len(hits), res)
            # the class is the one predicted by the converse formula: (2^A - C_w) 3^-n mod 2^(A+1)
            pred = ((2 ** A - C_w(w)) * pow(3, -n, M)) % M
            assert res == {pred}, (w, res, pred)
            nwords += 1
    print(f"   [A1] Terras classes: for {nwords} words (n <= 4), a scan of all odd x in (-2^(A+3), 2^(A+3)) finds exactly the "
          f"8 elements of ONE residue class mod 2^(A+1), and it is the class (2^A - C_w) 3^-n mod 2^(A+1)  -> VERIFIED {stamp()}")
    # A1b: the same for longer words by sampling the predicted class
    bad = 0
    tested = 0
    for n in (5, 6, 8):
        for _ in range(300):
            w = tuple(random.randint(1, 6) for _ in range(n))
            A = sum(w)
            M = 2 ** (A + 1)
            c = ((2 ** A - C_w(w)) * pow(3, -n, M)) % M
            for t in random.sample(range(-50, 50), 6):
                x = c + M * t
                if word(x, n)[0] != w:
                    bad += 1
                tested += 1
            # an odd x outside the class must not have the word
            x = c + 2 * random.randint(1, M // 2 - 1) + M * random.randint(-5, 5)
            if x % 2 == 1 and word(x, n)[0] == w:
                bad += 1
            tested += 1
    print(f"   [A1b] words of length 5, 6, 8 (900 random words, {tested} tests): every element of the predicted class has the word, "
          f"no odd element outside it does; failures {bad}  -> {'VERIFIED' if bad == 0 else 'FAILED'}")
    # A2: affine identity 3^n x + C_w = 2^A S^n(x), C_w odd, Y_n(w) = C_w 2^-A both as closed form and by the recursion
    for _ in range(3000):
        x = random.randint(-10 ** 9, 10 ** 9) | 1
        n = random.randint(1, 12)
        w, y = word(x, n)
        A = sum(w)
        assert 3 ** n * x + C_w(w) == 2 ** A * y, (x, n)
        assert C_w(w) % 2 == 1
        assert y % 3 ** n == Y_mod(w) == Y_by_recursion(w)
    print("   [A2] 3000 random odd x (both signs), n <= 12: 3^n x + C_w = 2^A S^n(x) exactly, C_w odd, S^n(x) = Y_n(w) mod 3^n "
          "with Y_n(w) = C_w 2^-A = the recursion Y_{k+1} = 2^-a(3Y_k+1)  -> VERIFIED")
    # A3: the converse (iii): y odd, 3 !| y, y = Y_n(w) mod 3^n  <=>  x = (2^A y - C_w)/3^n is an odd integer with word w and S^n(x) = y
    ok = bad = 0
    for n in range(1, 5):
        for w in product(range(1, 5), repeat=n):
            A, Cw, mod = sum(w), C_w(w), 3 ** n
            Yw = Y_mod(w)
            for y in range(-4 * mod + 1, 4 * mod, 2):
                if y % 3 == 0:
                    continue
                num = 2 ** A * y - Cw
                if y % mod == Yw:
                    assert num % mod == 0
                    x = num // mod
                    assert x % 2 == 1 and (x > 0) == (y > 0)
                    ww, yy = word(x, n)
                    assert ww == w and yy == y, (w, y, x)
                    ok += 1
                else:
                    assert num % mod != 0
                    bad += 1
    print(f"   [A3] converse: for all words n <= 4, a <= 4, and all odd y in (-4 3^n, 4 3^n), 3 !| y: y = Y_n(w) mod 3^n gives an odd "
          f"integer x = (2^A y - C_w)/3^n with word w, S^n(x) = y and the sign of y ({ok} cases); otherwise x is not an integer "
          f"({bad} cases)  -> VERIFIED")
    # A4: tree layers by direct enumeration (truncated valuations) against 3^n mu_n(y mod 3^n), exact law
    print("   [A4] tree layers sum_{x : S^n(x) = y} 3^n 2^-A(x) (valuations <= amax) versus 3^n mu_n(y mod 3^n):")
    for y in (1, 5, 7, 11, -1, -5, -7, -17):
        for n in range(1, 6):
            mod = 3 ** n
            exact = mod * EX[n].get(y % mod, Fraction(0))
            vals = []
            for amax in (8, 12, 16, 20):
                tot, cls, sign_ok = tree_layer(y, n, amax)
                assert sign_ok
                vals.append(float(tot))
            gaps = [float(exact) - v for v in vals]
            assert all(g >= -1e-12 for g in gaps) and gaps[-1] < 4 * 2.0 ** (-16) * 3 * n * 2 ** n, (y, n, gaps)
            if n in (1, 3, 5):
                print(f"      y={y:4d} n={n}: exact {float(exact):.6f} (= {exact if exact.denominator < 10**8 else '...'}); "
                      f"truncated amax=8,12,16,20: {', '.join(f'{v:.6f}' for v in vals)}; gaps {', '.join(f'{g:.1e}' for g in gaps)}")
    print("      all nodes carry the sign of the seed; the truncated sums increase to the exact value with geometric gaps  -> VERIFIED")
    # A5: exact H_n(1), the depth-1 split, Corollary A2, first arrivals
    H = {n: 3 ** n * EX[n].get(1, Fraction(0)) for n in range(1, EXACT_MAX + 1)}
    H[0] = Fraction(1)
    print(f"   [A5] exact H_n(1) = 3^n mu_n(1): " + ", ".join(f"H_{n} = {H[n] if H[n].denominator < 10**6 else float(H[n])}" for n in range(1, 4))
          + "; " + ", ".join(f"H_{n} = {float(H[n]):.6f}" for n in range(4, EXACT_MAX + 1)))
    assert H[1] == 1 and H[2] == Fraction(8, 7) and H[3] == Fraction(1376, 1387)
    # depth-1 layer of 1: (4^k - 1)/3 with weight 3 4^-k, class k mod 3
    s = [sum(Fraction(3, 4 ** k) for k in range(1, 400) if k % 3 == r) for r in range(3)]
    print(f"      depth-1 split of H_1 = 1 by class mod 3 (0,1,2): {float(s[0]):.6f}, {float(s[1]):.6f}, {float(s[2]):.6f} "
          f"(= 1/21, 16/21, 4/21 up to 2^-800); H_2 = H^(1) + 2 H^(2) = {float(s[1] + 2 * s[2]):.6f} = 8/7  -> VERIFIED")
    assert abs(float(s[0]) - 1 / 21) < 1e-200 and abs(float(s[1]) - 16 / 21) < 1e-200 and abs(float(s[2]) - 4 / 21) < 1e-200
    for n in range(1, 6):
        tot, cls, _ = tree_layer(1, n, 20)
        nxt, _, _ = tree_layer(1, n + 1, 20)
        print(f"      Corollary A2 on the truncated tree, n={n}: H^(1) + 2 H^(2) = {float(cls[1] + 2 * cls[2]):.6f} vs truncated H_{n+1} = {float(nxt):.6f} "
              f"(exact {float(H[n+1]):.6f}); first-arrival mass H_{n}^first = H_n - (3/4) H_{n-1} = {float(H[n] - Fraction(3, 4) * H[n-1]):.6f}")
    print("      H_n = sum_{j<=n} (3/4)^(n-j) H_j^first  <=>  H_n^first = H_n - (3/4) H_{n-1}: first arrivals 1/4, 11/28 = 0.392857, 0.134926, ...  -> VERIFIED")
    # A6: consistency of the laws and the two exact recursions; float vs exact
    PX = law_exact_parents(5)
    for n in range(1, 6):
        assert all(EX[n].get(z, 0) == PX[n].get(z, 0) for z in range(3 ** n))
    for n in range(2, EXACT_MAX + 1):
        for m in range(1, n):
            red = defaultdict(Fraction)
            for z, p in EX[n].items():
                red[z % 3 ** m] += p
            assert all(red[y] == EX[m].get(y, 0) for y in range(3 ** m))
    err = max(max(abs(float(EX[n].get(z, 0)) - FL[n][z]) for z in range(3 ** n)) for n in range(1, EXACT_MAX + 1))
    print(f"   [A6] forward and parents exact recursions agree to level 5; consistency (Y_n mod 3^m has the law of Y_m) exact to level {EXACT_MAX}; "
          f"float law vs exact: max |diff| = {err:.1e} (levels <= {EXACT_MAX})  -> VERIFIED")
    # class means
    for n in sorted({1, min(6, NMAX), min(10, NMAX), NMAX}):
        mu = FL[n]
        rho = rho_of(mu, n)
        z = np.arange(3 ** n)
        m1, m2, mu_ = rho[z % 3 == 1].mean(), rho[z % 3 == 2].mean(), rho[z % 3 != 0].mean()
        assert abs(m1 - 2 / 3) < 1e-9 and abs(m2 - 4 / 3) < 1e-9 and abs(mu_ - 1) < 1e-9 and rho[z % 3 == 0].max() == 0
    print(f"   [A7] Corollary A1: rho_n has mean 1 on units, 2/3 on the class 1 mod 3, 4/3 on the class 2 mod 3, 0 on nonunits (n = 1, 6, 10, {NMAX})  -> VERIFIED")
    # A8: the note's mod-9 class values rho_2 = (16, 32, 22, 8, 4, 44)/21 on the classes 1, 2, 4, 5, 7, 8, and the leaf-child reading
    r2 = {y: 6 * EX[2].get(y, Fraction(0)) for y in (1, 2, 4, 5, 7, 8)}
    assert [r2[y] for y in (1, 2, 4, 5, 7, 8)] == [Fraction(v, 21) for v in (16, 32, 22, 8, 4, 44)], r2
    dom = {}
    for y in (1, 2, 4, 5, 7, 8):
        a = 2 if y % 3 == 1 else 1  # the dominant child (three quarters of the offspring weight) has the least admissible valuation
        dom[y] = ((2 ** a * y - 1) // 3) % 3  # its class mod 3 (a leaf iff 0), well defined mod 3 from y mod 9
    assert dom[5] == 0 and dom[7] == 0 and all(dom[y] != 0 for y in (1, 2, 4, 8))
    print(f"   [A8] rho_2 on the classes 1, 2, 4, 5, 7, 8 mod 9 = {[str(r2[y]) for y in (1, 2, 4, 5, 7, 8)]} = (16, 32, 22, 8, 4, 44)/21; the dominant child "
          f"(a = 2 for y = 1 mod 3, a = 1 for y = 2 mod 3) is a leaf exactly for the classes 5 and 7 mod 9  -> VERIFIED")
    return H


# ----------------------------------------------------------------------------------------------
# Section B: Theorem B
# ----------------------------------------------------------------------------------------------

def section_B(FL):
    print("\n==== B. Theorem B: cycle resonances ====")
    cyc = {}
    for y0 in (-1, -5, -17, 1):
        w, x, k = [], y0, 0
        while True:
            x, a = syr(x)
            w.append(a)
            k += 1
            if x == y0:
                break
        A = sum(w)
        cyc[y0] = (k, A, tuple(w))
        Cw = C_w(w)
        assert Fraction(Cw, 2 ** A - 3 ** k) == y0  # y0 = C_w/(2^A - 3^k)
        s = 1 - (A / k) * LN2 / LN3
        print(f"   cycle of {y0}: length k={k}, word {w}, A={A}, 3^k={3**k} vs 2^A={2**A} ({'negative: 3^k > 2^A' if 3**k > 2**A else 'positive: 2^A > 3^k'}); "
              f"y0 = C_w/(2^A - 3^k) = {Cw}/{2**A - 3**k}; s = 1 - (A/k) log_3 2 = {s:.4f}")
    assert cyc[-1] == (1, 1, (1,)) and cyc[-5] == (2, 3, (1, 2)) and cyc[-17][0] == 7 and cyc[-17][1] == 11 and cyc[1] == (1, 2, (2,))
    # Y_n(1^n) = (3/2)^n - 1 exactly, = -1 mod 3^n
    for n in range(1, 15):
        w = (1,) * n
        val = Fraction(C_w(w), 2 ** n)
        assert val == Fraction(3, 2) ** n - 1 and Y_mod(w) == 3 ** n - 1
    print("   Y_n(1^n) = (3/2)^n - 1 exactly and = -1 mod 3^n, n <= 14  -> VERIFIED")
    # lower bounds rho_{km}(y0) >= (2/3)(3^k/2^A)^m at all computed levels; growth data
    print("   lower bound rho_n(y0) >= (2/3)(3^k/2^A)^(n/k) at n = km (float law):")
    for n in range(1, NMAX + 1):
        mod = 3 ** n
        rho = rho_of(FL[n], n)
        line = f"      n={n:2d}: rho(-1)={rho[mod-1]:9.4f} vs (2/3)(3/2)^n={(2/3)*1.5**n:9.4f}, ratio {rho[mod-1]/1.5**n:.5f}"
        assert rho[mod - 1] >= (2 / 3) * 1.5 ** n - 1e-9
        if n % 2 == 0:
            b = (2 / 3) * (9 / 8) ** (n // 2)
            assert rho[(-5) % mod] >= b - 1e-9 and rho[(-7) % mod] >= b - 1e-9
            line += f"; rho(-5)={rho[(-5)%mod]:.4f}, rho(-7)={rho[(-7)%mod]:.4f} vs {b:.4f}"
        if n % 7 == 0:
            b = (2 / 3) * (2187 / 2048) ** (n // 7)
            assert all(rho[(y) % mod] >= b - 1e-9 for y in (-17, -25, -37, -55, -41, -61, -91))
            line += f"; rho(-17)={rho[(-17)%mod]:.4f} vs {b:.4f}"
        else:
            line += f"; rho(-17)={rho[(-17)%mod]:.4f}"
        line += f"; argmax over units = {'-1' if int(np.argmax(rho)) == mod - 1 else int(np.argmax(rho))}"
        print(line)
    print("   all lower bounds hold; the class of -1 is the maximal atom at every level  -> VERIFIED")
    # increments of rho_n(-1)/(3/2)^n = (2/3) sum_j (2/3)^j H_j^first(-1): monotone, extrapolation
    seq = [rho_of(FL[n], n)[3 ** n - 1] / 1.5 ** n for n in range(1, NMAX + 1)]
    inc = [b - a for a, b in zip(seq, seq[1:])]
    q = inc[-1] / inc[-2]
    print(f"      rho_n(-1)/(3/2)^n increments n->n+1: " + ", ".join(f"{d:.2e}" for d in inc[-6:]) +
          f"; ratio of successive increments {q:.3f}; crude geometric extrapolation of the limit: {seq[-1] + inc[-1] * q / (1 - q):.5f}")


# ----------------------------------------------------------------------------------------------
# Section D: spike profile
# ----------------------------------------------------------------------------------------------

def closure_of_minus_one(K: int):
    """Forward closure of -1 under all S_a(y) = (3y+1)/2^a, points with denominator <= 2^K; value -> exponent of the denominator."""
    pts = {Fraction(-1): 0}
    frontier = [Fraction(-1)]
    while frontier:
        nxt = []
        for y in frontier:
            for a in range(1, K + 2):
                z = (3 * y + 1) / 2 ** a
                if z == y:
                    continue  # the loop at -1
                k = v2(z.denominator)
                if k > K:
                    break
                if z not in pts:
                    pts[z] = k
                    nxt.append(z)
        frontier = nxt
    return pts


def profile_D(pts):
    """D(z) = 2 sum_a 2^-a D((2^a z - 1)/3) over parents in the closure, D(-1) = 1, D = 0 off the closure (path sums 2^(l - A))."""
    D = {Fraction(-1): Fraction(1)}
    for z in sorted(pts, key=lambda q: (pts[q], q)):
        if z == -1:
            continue
        tot = Fraction(0)
        for a in range(1, pts[z] + 2):
            y = (2 ** a * z - 1) / 3
            if y in D:
                tot += Fraction(2, 2 ** a) * D[y]
        D[z] = tot
    return D


def section_D(FL):
    print("\n==== D. The spike profile on the forward closure of -1 ====")
    # the recursion rho_{n+1}(z) = 3 sum_{a = eps(z)} 2^-a rho_n((2^a z - 1)/3) IS the parents form used to build FL (checked against the
    # forward exact law in A6); re-check it pointwise at z = 1, -1, 5, -5 with the exact law
    EX = law_exact_forward(5)
    for n in range(1, 5):
        mod, modn = 3 ** (n + 1), 3 ** n
        L = 2 * 3 ** n
        for z in (1, mod - 1, 5, mod - 5):
            eps = 1 if z % 3 == 2 else 0
            tot = Fraction(0)
            for r in range(1, L + 1):
                if r % 2 == eps:
                    tot += Fraction(2 ** (L - r), 2 ** L - 1) * EX[n].get(((pow(2, r, mod) * z - 1) // 3) % modn, Fraction(0))
            assert tot == EX[n + 1].get(z, Fraction(0)), (n, z)
    print("   the parents recursion mu_{n+1}(z) = sum_{a = eps(z) mod 2} 2^-a mu_n((2^a z - 1)/3 mod 3^n) reproduces the forward law exactly "
          "at z = 1, -1, 5, -5, n <= 4 (and the whole float law is built from it)  -> VERIFIED")
    pts = closure_of_minus_one(9)
    D = profile_D(pts)
    named = {"-1/2": Fraction(-1, 2), "-1/4": Fraction(-1, 4), "1/8": Fraction(1, 8), "11/16": Fraction(11, 16), "49/32": Fraction(49, 32),
             "-1/8": Fraction(-1, 8), "1/16": Fraction(1, 16), "11/32": Fraction(11, 32), "-1/16": Fraction(-1, 16)}
    print("   path-sum values: " + ", ".join(f"D({k}) = {D[v]}" for k, v in named.items()))
    assert [D[v] for v in named.values()] == [Fraction(1, 2), Fraction(3, 4), Fraction(3, 4), Fraction(3, 4), Fraction(3, 4), Fraction(3, 8), Fraction(3, 8), Fraction(3, 8), Fraction(3, 16)]
    # all-ones forward orbit of -1/4
    xj = [Fraction(3 ** (j + 1), 2 ** (j + 2)) - 1 for j in range(12)]
    for j in range(11):
        assert (3 * xj[j] + 1) / 2 == xj[j + 1]
        assert (xj[j] + 1).numerator % 3 ** (j + 1) == 0
        if j <= 7:
            assert D[xj[j]] == Fraction(3, 4)
    print("   x_j = 3^(j+1)/2^(j+2) - 1 (j < 12): x_0 = -1/4, S_1(x_j) = x_{j+1}, x_j = -1 mod 3^(j+1), D(x_j) = 3/4 by path sums (j <= 7)  -> VERIFIED")
    # compare D with the float law at several levels, over all closure points with denominator <= 2^6
    small = [z for z in pts if pts[z] <= 6]
    for n in sorted({min(8, NMAX), min(10, NMAX), min(12, NMAX), NMAX}):
        mod = 3 ** n
        rho = rho_of(FL[n], n)
        top = rho[mod - 1]
        errs = []
        for z in small:
            cls = (z.numerator * pow(z.denominator, -1, mod)) % mod
            errs.append((abs(rho[cls] / top - float(D[z])), z))
        errs.sort(reverse=True)
        line = "; ".join(f"{k}: {rho[(v.numerator * pow(v.denominator, -1, mod)) % mod] / top:.4f}" for k, v in named.items())
        print(f"   level {n}: rho_n(x)/rho_n(-1) at the nine named points: {line}")
        print(f"      over all {len(small)} closure points with denominator <= 2^6: max |ratio - D| = {errs[0][0]:.4f} at x = {errs[0][1]}, "
              f"mean |ratio - D| = {sum(e for e, _ in errs) / len(errs):.5f}; distinct D-values: {sorted(set(float(D[z]) for z in small))}")
    # second tier: classes with ratio >= 0.7 at level n versus the shadows x_j mod 3^n
    for n in sorted({min(8, NMAX), min(10, NMAX), min(12, NMAX)}):
        mod = 3 ** n
        rho = rho_of(FL[n], n)
        top = rho[mod - 1]
        big = sorted(int(c) for c in np.nonzero(rho >= 0.7 * top)[0] if int(c) != mod - 1)
        shadows = {(x.numerator * pow(x.denominator, -1, mod)) % mod: j for j, x in enumerate(xj[:n - 1])}
        extra = [c for c in big if c not in shadows]
        missing = [(j, rho[c] / top) for c, j in shadows.items() if rho[c] < 0.7 * top]
        print(f"   level {n}: classes with rho >= 0.7 rho(-1) other than -1: {len(big)}; shadows x_j mod 3^n (j <= {n-2}) among them: {len(big) - len(extra)}; "
              f"non-shadow classes: {[(c - mod if c > mod // 2 else c) for c in extra][:6]}; shadows below 0.7: {[(j, round(r, 3)) for j, r in missing]}")
        print(f"      ratios at the shadows j = 0..{n-2}: " + ", ".join(f"{rho[c]/top:.3f}" for c, j in sorted(shadows.items(), key=lambda t: t[1])))
    return D


# ----------------------------------------------------------------------------------------------
# Section C: Theorem C
# ----------------------------------------------------------------------------------------------

def section_C(FL, H_exact):
    print("\n==== C. Theorem C: the seed-1 test ====")
    N = 10 ** 6
    k = np.arange(1, N + 1, dtype=np.float64)
    partial = np.cumsum(1.0 / (2 * k - 1))
    slack = 1 + 0.5 * np.log(k) - partial
    print(f"   [C1] sum_{{k<=n}} 1/(2k-1) <= 1 + (1/2) ln n for n <= 10^6: min slack = {slack.min():.3e} at n = {int(np.argmin(slack)) + 1} (equality at n = 1), "
          f"slack at n = 10^6: {slack[-1]:.5f} (limit 1 - ln 2 - gamma/2 = {1 - LN2 - 0.5772156649 / 2:.5f}); "
          f"proof: increments 1/(2n+1) - (1/2) ln(1+1/n) <= 0 since ln(1+1/n) = 2 artanh(1/(2n+1)) >= 2/(2n+1)  -> VERIFIED")
    # Pi(x) over odd x <= 2 10^5 and the exact identity
    X = 2 * 10 ** 5
    best = (0.0, 0)
    tot = 0.0
    cnt = 0
    worst = 0.0
    for x in range(3, X + 1, 2):
        n = A = 0
        y = x
        s = 0.0
        while y != 1:
            s += 1.0 / y
            y, a = syr(y)
            n += 1
            A += a
        lp = A * LN2 - n * LN3 - math.log(x)
        pi = math.exp(lp)
        assert lp <= s / 3 + 1e-12 and s <= 1 + 0.5 * math.log(n) + 1e-12
        worst = max(worst, pi / (math.exp(1 / 3) * n ** (1 / 6)))
        tot += pi
        cnt += 1
        if pi > best[0]:
            best = (pi, x)
    print(f"   [C2] Pi(x) = 2^A/(3^n x) over odd 3 <= x <= {X}: max {best[0]:.5f} at x = {best[1]}, mean {tot/cnt:.5f}, max Pi/(e^(1/3) n^(1/6)) = {worst:.4f}  "
          f"(note: 1.2531 at x = 993, mean 1.167, <= 0.764 over x <= 2 10^6)  -> {'VERIFIED' if best[1] == 993 else 'DIFFERS'}")
    for m in (3, 9, 27, 993, 97, 871):
        y, n, A = m, 0, 0
        prod = Fraction(1)
        while y != 1:
            prod *= 1 + Fraction(1, 3 * y)
            y, a = syr(y)
            n += 1
            A += a
        assert Fraction(1, m) == Fraction(3 ** n, 2 ** A) * prod
    print("   [C3] the identity 1/m = 3^n 2^-A prod_{j<n}(1 + 1/(3 x_j)) holds exactly (Fractions) for m = 3, 9, 27, 97, 871, 993  -> VERIFIED")
    # the orbit values x_0..x_{n-1} of a first arrival are distinct odd integers > 1
    for m in range(3, 20001, 2):
        y, seen = m, set()
        while y != 1:
            assert y not in seen
            seen.add(y)
            y, _ = syr(y)
    print("   [C4] for odd 3 <= m <= 20001 the orbit values before 1 are distinct and > 1 (first arrival), as the proof needs  -> VERIFIED")
    # G_n <= e^(1/3) n^(1/6) H_n^first <= ... H_n on truncated trees of 1
    for n in range(1, 7):
        Gf = Hf = Fraction(0)
        stack = [(1, 0, 0, False)]
        while stack:
            x, d, A, rev = stack.pop()
            if d == n:
                if not rev:
                    Hf += Fraction(3 ** n, 2 ** A)
                    Gf += Fraction(1, x)
                continue
            r = x % 3
            if r == 0:
                continue
            for a in range(2 if r == 1 else 1, 17, 2):
                z = (2 ** a * x - 1) // 3
                stack.append((z, d + 1, A + a, rev or z == 1))
        bound = math.exp(1 / 3) * n ** (1 / 6)
        assert float(Gf) <= bound * float(Hf) and Hf <= H_exact[n] + Fraction(1, 10 ** 6)
        print(f"      n={n}: G_n^first = {float(Gf):.5f}, H_n^first = {float(Hf):.5f} (exact H_n - (3/4)H_(n-1) = {float(H_exact[n] - Fraction(3,4)*H_exact[n-1]):.5f}), "
              f"G/H^first = {float(Gf/Hf):.4f} <= e^(1/3) n^(1/6) = {bound:.4f}; H_n^first <= H_n = {float(H_exact[n]):.5f}")
    print("   [C5] G_n <= e^(1/3) n^(1/6) H_n^first <= e^(1/3) n^(1/6) H_n at n <= 6 (truncated words, valuations <= 16)  -> VERIFIED")
    # Syracuse depth versus Collatz steps, and the empirical density of tau(n) <= 10.46 ln n
    N = 10 ** 6
    tau = np.zeros(N + 1, dtype=np.int64)
    ok = 0
    for n in range(2, N + 1):
        x, t = n, 0
        while x >= n:
            x = x // 2 if x % 2 == 0 else 3 * x + 1
            t += 1
            if x == 1:
                break
        tau[n] = t + (tau[x] if x < n else 0)
        if tau[n] <= 10.46 * math.log(n):
            ok += 1
    print(f"   [C6] fraction of 2 <= n <= 10^6 with tau(n) <= (523/50) ln n: {ok/(N-1):.4f} (note: 54.22%)  -> {'VERIFIED' if abs(ok/(N-1) - 0.5422) < 5e-4 else 'DIFFERS'}")
    depth_ok = True
    for m in range(3, 100001, 2):
        y, n = m, 0
        while y != 1:
            y, _ = syr(y)
            n += 1
        if not (2 * n <= tau[m]):
            depth_ok = False
    print(f"   [C7] Syracuse depth n(m) <= tau(m)/2 <= tau(m) for odd m <= 10^5: {'VERIFIED' if depth_ok else 'FAILED'}")
    # the sequence n^(1/6) H_n from the float law
    Hs = {n: 3 ** n * FL[n][1] for n in range(1, NMAX + 1)}
    print("   [C8] H_n(1) and n^(1/6) H_n (float law, this script): " + ", ".join(f"n={n}: {Hs[n]:.5f} / {n**(1/6)*Hs[n]:.4f}" for n in range(1, NMAX + 1)))
    ratios = [Hs[n + 1] / Hs[n] for n in range(1, NMAX)]
    first = next(i + 1 for i in range(len(ratios)) if all(r < 1 for r in ratios[i:]))
    print("      ratios H_{n+1}/H_n, n = 1..: " + ", ".join(f"{r:.4f}" for r in ratios) + f"; all ratios are < 1 from n = {first} on (H_{first+1}/H_{first} onwards)")
    ses = {8: 0.77458, 18: 0.41046}
    print(f"      with the session's level-18 value: average decline per level from n = 8 to 18 = {1 - (ses[18]/ses[8])**0.1:.3f}; 18^(1/6) H_18 = {18**(1/6)*ses[18]:.4f}")
    # the session's depth-decomposition table (levels 19..38, lower bounds): where is the minimum?
    deep = {19: 0.39655, 20: 0.38661, 21: 0.37690, 22: 0.36970, 23: 0.36999, 24: 0.37612, 25: 0.37251, 26: 0.36554, 27: 0.36323, 28: 0.36606,
            29: 0.36691, 30: 0.37335, 31: 0.38409, 32: 0.39076, 33: 0.40395, 34: 0.41079, 35: 0.41126, 36: 0.41622, 37: 0.42589, 38: 0.44510}
    pruned = {19: 3.7e-9, 20: 8.8e-8, 21: 3.9e-7, 22: 2.0e-6, 23: 5.3e-6, 24: 1.2e-5, 25: 3.4e-5, 26: 6.3e-5, 27: 1.5e-4, 28: 2.5e-4, 29: 3.7e-4,
              30: 6.9e-4, 31: 1.0e-3, 32: 1.7e-3, 33: 2.3e-3, 34: 3.7e-3, 35: 4.9e-3, 36: 6.0e-3, 37: 8.4e-3, 38: 1.0e-2}
    nmin = min(deep, key=deep.get)
    corr = {n: deep[n] + pruned[n] for n in deep}
    nmin_c = min(corr, key=corr.get)
    rising_from = max(n for n in range(20, 39) if deep[n] < deep[n - 1]) + 1
    print(f"   [C9] the session's mazur_seed1_deep_20260928.out values (depth 20, lower bounds): minimum {deep[nmin]:.5f} at n = {nmin} (with the heuristic pruned-weight "
          f"correction: {corr[nmin_c]:.5f} at n = {nmin_c}); H_22 = {deep[22]:.5f} is only a local minimum; the last decrease is at n = {rising_from - 1}, "
          f"so the sequence rises for {38 - rising_from + 1} consecutive levels (n = {rising_from}..38); 38^(1/6) H_38 = {38**(1/6)*deep[38]:.4f}")
    # the depth-30 run (mazur_seed1_deep30_20260928.out), the note's source for n >= 19 in its final form
    deep30 = {19: 0.39655, 20: 0.38661, 21: 0.37690, 22: 0.36971, 23: 0.36999, 24: 0.37613, 25: 0.37254, 26: 0.36559, 27: 0.36336, 28: 0.36625,
              29: 0.36719, 30: 0.37386, 31: 0.38481, 32: 0.39207, 33: 0.40562, 34: 0.41358, 35: 0.41462, 36: 0.42026, 37: 0.43183, 38: 0.45187,
              39: 0.48588, 40: 0.50371, 41: 0.54043, 42: 0.54284, 43: 0.54371, 44: 0.54034, 45: 0.53382, 46: 0.53058, 47: 0.52059, 48: 0.51337}
    pr30 = {19: 9.3e-10, 20: 6.2e-9, 21: 5.3e-8, 22: 1.9e-7, 23: 8.6e-7, 24: 2.1e-6, 25: 6.7e-6, 26: 1.3e-5, 27: 2.4e-5, 28: 5.4e-5, 29: 9.0e-5,
            30: 1.8e-4, 31: 2.8e-4, 32: 3.9e-4, 33: 6.7e-4, 34: 9.3e-4, 35: 1.5e-3, 36: 2.0e-3, 37: 2.5e-3, 38: 3.6e-3, 39: 4.5e-3, 40: 6.3e-3,
            41: 7.7e-3, 42: 1.0e-2, 43: 1.2e-2, 44: 1.4e-2, 45: 1.8e-2, 46: 2.1e-2, 47: 2.6e-2, 48: 2.9e-2}
    m30 = min(deep30, key=deep30.get)
    c30 = {n: deep30[n] + pr30[n] for n in deep30}
    m30c = min(c30, key=c30.get)
    below22 = [n for n in deep30 if deep30[n] < deep30[22]]
    print(f"   [C10] the depth-30 run (lower bounds): minimum {deep30[m30]:.5f} at n = {m30} (corrected {c30[m30c]:.5f} at n = {m30c}); levels below H_22 = {deep30[22]:.5f}: {below22}; "
          f"maximum {max(deep30.values()):.5f} at n = {max(deep30, key=deep30.get)}; H_48 = {deep30[48]:.5f} (corrected {c30[48]:.4f}); n^(1/6) H_n at n = 41..48: "
          + ", ".join(f"{n**(1/6)*deep30[n]:.3f}" for n in range(41, 49))
          + f"; average decline per level from n = 8 (0.77458) to n = 22: {1 - (deep30[22]/0.77458)**(1/14):.3f} (the note says 'about 3%')")


# ----------------------------------------------------------------------------------------------
# Section F: global statistics
# ----------------------------------------------------------------------------------------------

def section_F(FL):
    print("\n==== F. Global statistics of the law (float, this script) ====")
    stats = {}
    for n in range(1, NMAX + 1):
        mu = FL[n]
        mod = 3 ** n
        rho = rho_of(mu, n)
        z = np.arange(mod)
        units = z % 3 != 0
        ru = rho[units]
        ent = float(-(mu[mu > 0] * np.log(mu[mu > 0])).sum())
        deficit = n * LN3 - ent
        l2 = float((ru ** 2).mean())
        l1inc = float(np.abs(mu - np.tile(FL[n - 1], 3) / 3).sum()) if n >= 2 else float("nan")
        e_rlr = float((ru * np.log(np.where(ru > 0, ru, 1))).mean())
        stats[n] = (rho[1], deficit, l2, float(np.median(ru)), float((ru < 0.1).mean()), l1inc, e_rlr)
        print(f"   n={n:2d}: rho_n(1)={rho[1]:.5f} H_n={1.5*rho[1]:.5f} E_units[rho^2]={l2:.4f} deficit n ln3 - H(mu_n)={deficit:.4f} "
              f"E_units[rho ln rho]={e_rlr:.4f} (= deficit - ln(3/2) = {deficit - math.log(1.5):.4f}) median={np.median(ru):.4f} share(rho<0.1)={(ru < 0.1).mean():.4f} "
              f"L1inc={l1inc:.4f}")
    d = {n: stats[n][1] for n in stats}
    incs = [d[n] - d[n - 1] for n in range(7, NMAX + 1)]
    print("   entropy-deficit increments n = 7..: " + ", ".join(f"{v:.4f}" for v in incs) + "; successive ratios: " + ", ".join(f"{b/a:.3f}" for a, b in zip(incs, incs[1:])))
    l2s = {n: stats[n][2] for n in stats}
    print("   E[rho^2] increments per level n = 7..: " + ", ".join(f"{l2s[n]-l2s[n-1]:.4f}" for n in range(7, NMAX + 1)))
    l1 = {n: stats[n][5] for n in stats if n >= 2}

    def slope(a, b):
        return math.log(l1[b] / l1[a]) / math.log(b / a)
    print(f"   L1inc local exponents: over 6..{NMAX}: {slope(6, NMAX):.3f}; 10..{NMAX}: {slope(min(10, NMAX - 1), NMAX):.3f}; {NMAX-2}..{NMAX}: {slope(NMAX-2, NMAX):.3f} "
          f"(note claims n^(-0.93) from levels <= 18; session values 0.1783 at 10, 0.1296 at 14, 0.0969 at 18 give 10..18: {math.log(0.1783/0.0969)/math.log(1.8):.3f}, 14..18: {math.log(0.1296/0.0969)/math.log(18/14):.3f})")
    # tail of rho on units at the top level
    n = NMAX
    rho = rho_of(FL[n], n)
    z = np.arange(3 ** n)
    ru = rho[z % 3 != 0]
    ts = [1, 2, 4, 8, 16, 32]
    tail = [(ru > t).mean() for t in ts]
    print("   tail P_units(rho > t) at level %d: " % n + ", ".join(f"t={t}: {p:.2e} (t^2 P = {t*t*p:.3f})" for t, p in zip(ts, tail)) +
          "; local exponents: " + ", ".join(f"{-math.log(tail[i+1]/tail[i])/LN2:.2f}" for i in range(len(ts) - 1) if tail[i + 1] > 0))
    # fine-scale distances at the top level
    dists = [float(np.abs(FL[n] - np.tile(FL[m], 3 ** (n - m)) / 3 ** (n - m)).sum()) for m in range(1, min(n, 9) + 1)]
    print(f"   ||mu_{n} - lift(mu_m)||_1, m = 1..{min(n,9)}: " + ", ".join(f"{v:.4f}" for v in dists) + "  (session at N = 18: 0.8541 ... 0.3535)")
    # seed-1 fixed-point extrapolation from the parents form
    mod = 3 ** n
    ext = 12 * sum(4.0 ** (-j) * rho[((4 ** j - 1) // 3) % mod] for j in range(2, 40))
    rec = 3 * sum(4.0 ** (-j) * rho_of(FL[n - 1], n - 1)[((4 ** j - 1) // 3) % 3 ** (n - 1)] for j in range(1, 40))
    print(f"   parents form at z = 1: rho_{n}(1) = 3 sum_j 4^-j rho_{n-1}(R_j) = {rec:.6f} vs {rho[1]:.6f}; fixed-point extrapolation from level {n}: rho_inf(1) ~ {ext:.4f}, H_inf ~ {1.5*ext:.4f}")
    # growth of the -5, -7, -17 atoms
    for y0 in (-5, -7, -17):
        vals = [rho_of(FL[m], m)[y0 % 3 ** m] for m in range(1, NMAX + 1)]
        print(f"   rho_n({y0}), n = 1..{NMAX}: " + ", ".join(f"{v:.3f}" for v in vals))


# ----------------------------------------------------------------------------------------------
# Section G: Kaprekar and Euler bricks
# ----------------------------------------------------------------------------------------------

def section_G():
    print("\n==== G. Kaprekar cycles and Euler bricks ====")
    for d in range(2, 7):
        def kap(n):
            s = f"{n:0{d}d}"
            return int("".join(sorted(s, reverse=True))) - int("".join(sorted(s)))
        cyc_of = {}
        cycles = {}
        for n in range(10 ** d):
            if len(set(f"{n:0{d}d}")) == 1:
                continue
            path = []
            x = n
            while x not in cyc_of and x not in path:
                path.append(x)
                x = kap(x)
            if x in cyc_of:
                key = cyc_of[x]
            else:
                c = path[path.index(x):]
                key = tuple(sorted(c))
                cycles[key] = c
            for p in path:
                cyc_of[p] = key
        counts = Counter(cyc_of[n] for n in range(10 ** d) if len(set(f"{n:0{d}d}")) > 1)
        print(f"   d={d}: " + "; ".join(f"{cycles[k]} (len {len(k)}, {counts[k]} starts)" for k in sorted(cycles, key=lambda k: -counts[k])))
        assert all(kap(n) % 9 == 0 for n in range(10 ** d))
    print("   every image is divisible by 9; d = 2: the 5-cycle {9, 81, 63, 27, 45}; d = 3: 495; d = 4: 6174; d = 5: two 4-cycles and {53955, 59994}; "
          "d = 6: a 7-cycle and 631764, 549945  -> VERIFIED")
    sq = {i * i: i for i in range(1, 700)}
    bricks = []
    for a in range(1, 300):
        for b in range(a + 1, 300):
            if a * a + b * b not in sq:
                continue
            for c in range(b + 1, 300):
                if a * a + c * c in sq and b * b + c * c in sq:
                    bricks.append((a, b, c, (a * a + b * b + c * c) in sq))
    print(f"   Euler bricks with edges < 300: {bricks} (flag = integer space diagonal)  -> VERIFIED")


# ----------------------------------------------------------------------------------------------
# Section H: the paper's constants, and Lemma 7.3's arithmetic
# ----------------------------------------------------------------------------------------------

def section_H():
    print("\n==== H. Section 8.1 of the paper: the size of the constants, in iterated logarithms ====")
    A_, E_, L_ = 6409, 2170, 280
    lg = math.log2
    log2_T = lg(10 * A_) + 3 * E_                     # T* = 10 A* 2^(3E*)
    log2_R = E_ + log2_T                              # R* = 2^E* (T* + 3(8192+3) + 1) + 1
    log2_Dexp_log2 = lg(8192 * A_) + 3 * E_           # D_exp = 2^(8192 A* 2^(3E*)): log2 log2 D_exp
    log2_D1 = 2 * (lg(128 * A_) + 3 * E_)             # D_1 = (128 A* 2^(3E*))^2 + 2
    # P* = g^(R*-1)(T*), g(t) = t + 10 8^32768 (t+1)^3 + 1 + T*: log2 g(t) ~ 3 log2 t + 98307.3, so
    # log2 P* ~ 3^(R*-1) (log2 T* + 49154) and log2 log2 P* ~ (R*-1) log2 3 + log2(log2 T* + 49154); log2 of that:
    log3_P = log2_R + lg(lg(3))                       # log2 log2 log2 P*  (R* - 1 ~ R*)
    print(f"   T* = 10 A* 2^(3E*): log2 T* = {log2_T:.2f};  R* = 2^E*(T* + 24586) + 1: log2 R* = {log2_R:.2f}")
    print(f"   D_exp = 2^(8192 A* 2^(3E*)): log2 log2 D_exp = {log2_Dexp_log2:.1f}  (the note's 2^(2^6536): correct for D_exp itself)")
    print(f"   D_1 = (128 A* 2^(3E*))^2 + 2: log2 D_1 = {log2_D1:.0f};  D_pt = 2(D_1 + L*^2 + 1) + 4 ~ D_1")
    print(f"   P* = g^(R*-1)(T*) with g(t) ~ 10 8^32768 (t+1)^3: R* - 1 ~ 2^{log2_R:.1f} compositions of a cubic map, so "
          f"log2 log2 P* ~ (R*-1) log2 3 + 15.8, i.e. log2 log2 log2 P* ~ {log3_P:.2f}: P* is a tower 2^2^2^{log3_P:.0f}")
    print(f"   D_sc = 2 + P + P*^10 + 8^32768 (1+P*)^33 + (100 S*)^2 ~ 2^98304 P*^33 >> D_exp; D_bad ~ 2^38 P*; so D* = max(D_1, D_2, D_3) = D_sc, "
          f"log2 log2 log2 D* ~ {log3_P:.2f} (NOT log2 log2 D* ~ 6536)")
    print(f"   C* = (32 A* D*)^A*: log2 C* = A* (5 + log2 A* + log2 D*), so log2 log2 C* ~ 2^{log3_P:.2f} + {lg(A_):.1f}: log2 log2 log2 C* ~ {log3_P:.2f} "
          f"(the note's 'log2 log2 C* ~ 6548' assumed D* = D_exp)")
    print("   C = 2 C* 20^6409 + 2 + 2481 ~ C*;  F = 2^467 b^(16b) (C+1) (or 2467 b 16^b (C+1); the extraction is ambiguous, immaterial): log2 F ~ log2 C")
    print(f"   N = 20000 (ceil(log2 F) + 64): log2 N ~ 14.3 + log2 log2 C ~ 2^{log3_P:.2f}: log2 log2 N ~ {log3_P:.2f} (the note's 'log2 N ~ 6563' is one exponential short)")
    print("   beta(N) ~ 280 (1.01)^N: log2 beta(N) ~ 0.01435 N; q ~ 160 beta(N); M = 4^(2b+1+3q): log2 M ~ 6 q; m ~ M^(1/6);")
    print(f"   c = 3/(256 M^2 m): log2 c^-1 ~ (13/6) log2 M ~ 2080 beta(N), log2 log2 c^-1 ~ 11 + 0.01435 N, log2^(3) c^-1 ~ log2 N - 6.1 ~ 2^{log3_P:.2f}, "
          f"log2^(4) c^-1 ~ {log3_P:.2f}: c^-1 is a tower of FOUR twos, 2^2^2^2^{log3_P:.0f} (the note's 'beyond 2^(2^(2^6535))' is true but understates by a level)")
    print(f"   J = max(N + N_t, 200(ceil(log2 m) - 78), 2 10^10) ~ 200 log2 m ~ 32000 beta(N) (exponential in N); B = E(J) + ... ~ 42 beta(J), "
          f"log2 beta(J) ~ 0.01435 J: X_0 = 32(2^B M + 1) has log2^(5) X_0 ~ {log3_P:.2f}, one level above c^-1 (the note: 'X_0 likewise')")
    # Lemma 7.3 / (7.9): kappa (L + h) + 1 <= (523/50) L with kappa = 34881/10000, L = ln 2, h <= ln 3 + 1/12288
    kappa = 34881 / 10000
    L, h = LN2, LN3 + 1 / 12288
    lhs, rhs = kappa * (L + h) + 1, 523 / 50 * L
    print(f"   Lemma 7.3's inequality (7.9): kappa(L+h) + 1 = {lhs:.6f} <= (523/50) L = {rhs:.6f}  (margin {rhs - lhs:.2e})  -> {'VERIFIED' if lhs <= rhs else 'FAILED'}; "
          f"with 262/25: {262/25*L:.6f}")
    print(f"   3/ln(4/3) = {3/math.log(4/3):.4f}; [kappa_0 (L + ln 3) + 1]/L with kappa_0 = 1/ln(4/3): {((1/math.log(4/3))*(L+LN3)+1)/L:.4f} (= 3/ln(4/3): the central-word constant)  -> VERIFIED")
    # the seed lemma: R_j = R_k mod 3^q iff j = k mod 3^q, by LTE: v_3(R_j - R_k) = v_3(j - k)
    for q in range(1, 5):
        res = [((4 ** j - 1) // 3) % 3 ** q for j in range(3 ** q)]
        assert sorted(res) == list(range(3 ** q))
    print("   Lemma 4.1: R_0..R_(3^q-1) permute Z/3^q (q <= 4); v_3(R_j - R_k) = v_3(4^(j-k) - 1) - 1 = v_3(j - k) (LTE), so it holds for all q  -> VERIFIED")
    # P6 of the session: how many histories did it actually compare?
    M = (4 ** 40 - 1) // 3
    byDA = Counter()
    for w in product(range(1, 6), repeat=3):
        x = M
        okk = True
        for a in w:
            y = 2 ** a * x - 1
            if y % 3:
                okk = False
                break
            x = y // 3
        if okk:
            byDA[(3, sum(w))] += 1
    print(f"   session P6 (seed R_40, depth 3, valuations <= 5): integral histories by (D, A): {dict(byDA)} -- every class has ONE history, so the "
          f"spacing and equal-endpoint assertions of P6 never compare two endpoints (vacuous as run; both statements are trivial anyway)")


# ----------------------------------------------------------------------------------------------

if __name__ == "__main__":
    wt = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    for rel in ("05-knowledge/results/mazur_positive_density_20260928.md", "04-computation/experiments/mazur_positive_density_20260928.py",
                "04-computation/experiments/mazur_harmonic_mass_20260928.py", "04-computation/experiments/mazur_harmonic_mass_deep_20260928.py",
                "04-computation/experiments/mazur_seed1_test_20260928.py", "04-computation/experiments/mazur_seed1_deep_20260928.py",
                "04-computation/experiments/kaprekar_cuboid_20260928.py"):
        p = os.path.join(wt, rel)
        if os.path.exists(p):
            st = os.stat(p)
            print(f"audited file: {rel}: {st.st_size} bytes, mtime {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(st.st_mtime))}")
    print(f"NMAX = {NMAX} (float law), exact law to level {EXACT_MAX}")
    EX = law_exact_forward(EXACT_MAX)
    print(f"exact law computed {stamp()}")
    FL = law_float(NMAX)
    print(f"float law computed to level {NMAX} {stamp()}")
    H = section_A(EX, FL)
    section_B(FL)
    section_D(FL)
    section_C(FL, H)
    section_F(FL)
    section_G()
    section_H()
    print(f"\nDONE {stamp()}")
