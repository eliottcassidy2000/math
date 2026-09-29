#!/usr/bin/env python3
"""Small-scale reproduction of the elementary components of L. Mazur,
"Explicit positive-density Collatz convergence in logarithmic time" (v2.1, 2026-09-06):
the seed lemma (4^j - 1)/3, inverse histories with prescribed valuations, the affine
identity and the weight 3^d 2^-A, the transfer operator on Z/3^t and Lemma 2.1, the
exact agreement of transfer fibres with actual inverse orbits, the reference density
rho_q (exact rationals), the mean identity (4.1), the residue pullback of a seed, the
deterministic spread of endpoints in residue classes, the time constant 3/log(4/3)
along typical inverse histories, and the empirical density of tau(n) <= 10.46 log n.
None of this tests the explicit mixing constant (Section 8.1 of the paper) or the
Lean development; it tests what is checkable at small scale.

Session: opus, collatz-poset-dag-20260927 (S19), 2026-09-28.
Run: python 04-computation/experiments/mazur_positive_density_20260928.py
"""
from __future__ import annotations

import math
import random
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product

LOG43 = math.log(4 / 3)


def v2(x: int) -> int:
    return (x & -x).bit_length() - 1


def S(x: int) -> tuple[int, int]:
    """Syracuse map on odd x: (S(x), a(x))."""
    y = 3 * x + 1
    a = v2(y)
    return y >> a, a


def tau(n: int) -> int:
    """Ordinary Collatz steps (T(n) = n/2 or 3n+1) to reach 1."""
    t = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        t += 1
    return t


# ----------------------------------------------------------------------------
# P1: the seeds (4^j - 1)/3
# ----------------------------------------------------------------------------

def v3(x: int) -> int:
    c = 0
    while x % 3 == 0:
        x //= 3
        c += 1
    return c


def part1():
    print("== P1: seeds R_j = (4^j - 1)/3 (Lemma 4.1) ==")
    for j in range(1, 200):
        R = (4 ** j - 1) // 3
        assert R % 2 == 1 and S(R) == (1, 2 * j)
        assert v3(4 ** j - 1) == 1 + v3(j)
    for q in range(1, 6):
        res = [((4 ** j - 1) // 3) % 3 ** q for j in range(3 ** q)]
        assert sorted(res) == list(range(3 ** q)), q
    print("   R_j odd, S(R_j) = 1 with valuation 2j, v_3(4^j - 1) = 1 + v_3(j) (j < 200); R_0..R_(3^q - 1) permute Z/3^q for q <= 5")


# ----------------------------------------------------------------------------
# P2: inverse histories, affine identity, weights
# ----------------------------------------------------------------------------

def inverse_history(R: int, word: tuple[int, ...]):
    """Word (a_1, ..., a_d) listed from the R end outward: R_(i+1) = (2^(a_(i+1)) R_i - 1)/3.
    Returns the list of endpoints [R = R_0, R_1, ..., R_d] or None if a step is not integral."""
    hist = [R]
    for a in word:
        y = 2 ** a * hist[-1] - 1
        if y % 3:
            return None
        hist.append(y // 3)
    return hist


def affine_source(R: int, word: tuple[int, ...]) -> Fraction:
    """x = (2^A R - C)/3^d with C = sum_j 3^(j-1) 2^(A - (a_1 + ... + a_j)) (as a rational)."""
    d = len(word)
    A = sum(word)
    C = 0
    part = 0
    for j, a in enumerate(word, start=1):
        part += a
        C += 3 ** (j - 1) * 2 ** (A - part)
    return Fraction(2 ** A * R - C, 3 ** d)


def part2():
    print("\n== P2: inverse histories with prescribed valuations; affine identity; weight x * 3^d 2^-A <= R ==")
    rng = random.Random(1)
    checked = 0
    for _ in range(4000):
        R = rng.randrange(1, 10 ** 6) | 1
        d = rng.randint(1, 8)
        word = tuple(rng.randint(1, 6) for _ in range(d))
        hist = inverse_history(R, word)
        x = affine_source(R, word)
        if hist is None:
            assert x.denominator != 1 or x <= 0 or any(w % 2 == 0 for w in [1]), (R, word)  # not integral (or would be)
            continue
        src = hist[-1]
        assert src == x, (R, word, src, x)
        # forward orbit from the source returns to R with the reversed valuations
        y = src
        vals = []
        for _ in range(d):
            y, a = S(y)
            vals.append(a)
        assert y == R and tuple(vals) == tuple(reversed(word))
        # weight and source charge: x * 3^d / 2^A <= R  (strictly: 3^d x < 2^A R)
        assert 3 ** d * src < 2 ** sum(word) * R
        checked += 1
    print(f"   {checked} integral histories: source = affine formula, forward orbit returns with reversed valuations exactly, 3^d x < 2^A R")
    # existence iff residue condition, step by step: 2^a R_i = 1 mod 3
    for R in range(1, 200, 2):
        for a in range(1, 8):
            y = 2 ** a * R - 1
            assert (y % 3 == 0) == ((2 ** a * R) % 3 == 1)
    print("   one inverse step exists iff 2^a R = 1 (mod 3), i.e. a has the parity fixed by R mod 3; multiples of 3 have no inverse step")


# ----------------------------------------------------------------------------
# P3: transfer operators on Z/3^t and Lemma 2.1; agreement with actual fibres
# ----------------------------------------------------------------------------

def F_w(z: int, word: tuple[int, ...], t: int) -> int:
    """F_w: G_t -> G_(t+d), F_w(z) = 3^d 2^-A z + sum_j 3^(j-1) 2^-(a_1+...+a_j) mod 3^(t+d)."""
    d = len(word)
    A = sum(word)
    mod = 3 ** (t + d)
    inv2A = pow(2, -A, mod)
    val = 3 ** d * inv2A * z
    part = 0
    for j, a in enumerate(word, start=1):
        part += a
        val += 3 ** (j - 1) * pow(2, -part, mod)
    return val % mod


def transfer(g: dict, word: tuple[int, ...], t: int) -> dict:
    """(T_w g)(y) = 3^d 2^-A sum_(z: F_w(z) = y) g(z) on G_(t+d), exact rationals."""
    d = len(word)
    A = sum(word)
    out = defaultdict(Fraction)
    wgt = Fraction(3 ** d, 2 ** A)
    for z in range(3 ** t):
        out[F_w(z, word, t)] += wgt * g[z]
    return out


def mean(g: dict, q: int) -> Fraction:
    return sum(g.values(), Fraction(0)) / 3 ** q


def part3():
    print("\n== P3: transfer operators: injectivity, Lemma 2.1, and exact agreement with inverse-orbit fibres ==")
    rng = random.Random(2)
    for _ in range(200):
        t = rng.randint(1, 3)
        d = rng.randint(1, 3)
        word = tuple(rng.randint(1, 5) for _ in range(d))
        imgs = [F_w(z, word, t) for z in range(3 ** t)]
        assert len(set(imgs)) == 3 ** t  # injective
        g = {z: Fraction(rng.randint(-5, 5)) for z in range(3 ** t)}
        Tg = transfer(g, word, t)
        assert mean(Tg, t + d) == Fraction(1, 2 ** sum(word)) * mean(g, t)
        gabs = {z: abs(v) for z, v in g.items()}
        assert mean({y: abs(v) for y, v in Tg.items()}, t + d) == Fraction(1, 2 ** sum(word)) * mean(gabs, t)
    print("   F_w injective and <T_w g> = 2^-A <g> (with and without absolute values), 200 random cases")
    # agreement: for an integer endpoint R, the word w admits an actual inverse history iff R mod 3^(t+d) lies in Im(F_w)
    # (with t = 0: the image is a single residue class mod 3^d); the source residue is the unique preimage
    for R in range(1, 3000, 2):
        for d in range(1, 4):
            for word in product(range(1, 5), repeat=d):
                hist = inverse_history(R, word)
                y = R % 3 ** d
                img = {F_w(0, word, 0)}  # t = 0: G_0 is a point
                assert (hist is not None) == (y in img), (R, word)
    print("   for every odd R < 3000 and every word of length <= 3 with valuations <= 4: an inverse history exists iff R mod 3^d is in Im(F_w)")


# ----------------------------------------------------------------------------
# P4: the reference density rho_q, exactly
# ----------------------------------------------------------------------------

def rho(q: int) -> dict:
    """rho_q(y) = (2/3) 3^q mu_q(y), mu_q the law of Y_q = 2^-A (3 Y_(q-1) + 1) mod 3^q, A geometric(1/2) on {1,2,...}.
    Exact: 2^-A mod 3^q depends on A mod L, L = 2*3^(q-1); P(A = r mod L) = 2^-r/(1 - 2^-L)."""
    mu = {0: Fraction(1)}  # level 0: Y_0 = 0 on G_0
    for lev in range(1, q + 1):
        mod = 3 ** lev
        L = 2 * 3 ** (lev - 1)
        new = defaultdict(Fraction)
        for y, p in mu.items():
            for r in range(1, L + 1):
                pr = Fraction(2 ** (L - r), 2 ** L - 1)  # P(A = r mod L) for r in 1..L
                y2 = (pow(2, -r, mod) * (3 * y + 1)) % mod
                new[y2] += p * pr
        mu = dict(new)
    return {y: Fraction(2, 3) * 3 ** q * mu.get(y, Fraction(0)) for y in range(3 ** q)}


def part4():
    print("\n== P4: the reference density rho_q (exact) ==")
    rhos = {q: rho(q) for q in range(1, 6)}
    for q, r in rhos.items():
        assert mean(r, q) == Fraction(2, 3)
        assert all(r[y] == 0 for y in range(3 ** q) if y % 3 == 0)
        units = [float(r[y]) for y in range(3 ** q) if y % 3]
        print(f"   q={q}: <rho_q> = 2/3, rho_q = 0 on nonunits; on units min {min(units):.4f} max {max(units):.4f}")
    # l1 distance between rho_q and the lift of rho_m
    for m in range(1, 5):
        for q in range(m + 1, 6):
            dist = sum(abs(rhos[q][y] - rhos[m][y % 3 ** m]) for y in range(3 ** q)) / 3 ** q
            print(f"   <|rho_{q} - rho_{m} o pi|> = {float(dist):.5f}", end=";")
    print()


# ----------------------------------------------------------------------------
# P5: the mean identity (4.1) and the residue pullback of a seed
# ----------------------------------------------------------------------------

def first_crossing_words(b: int, maxlen: int):
    """Small-scale analogue of W(b, u, K): words whose partial valuation sums first reach 2*len at length in
    [b, maxlen] ... we use the simplest prefix-free family: words of length exactly b (all valuations 1..4).
    Prefix-freeness holds trivially for fixed length."""
    return [w for w in product(range(1, 5), repeat=b)]


def part5():
    print("\n== P5: mean identity (4.1) for a prefix-free family, and the exact residue pullback of a seed ==")
    # family V: all words of length 2 with valuations in 1..4 (prefix-free); p(V) = sum 2^-A
    V = list(product(range(1, 5), repeat=2))
    pV = sum(Fraction(1, 2 ** sum(w)) for w in V)
    k = 2
    rk = rho(k)
    Q = k + 2
    comp = defaultdict(Fraction)
    for w in V:
        Tw = transfer(rk, w, k)
        for y, val in Tw.items():
            comp[y] += val
    m = mean(comp, Q)
    assert m == Fraction(2, 3) * pV
    print(f"   family V = words of length 2, valuations 1..4: p(V) = {pV} = {float(pV):.4f}; <sum_w T_w rho_2> = (2/3) p(V) = {float(m):.4f}  (Lemma 2.1 summed)")
    assert all(comp[y] == 0 for y in range(3 ** Q) if y % 3 == 0)
    best = max((y for y in range(3 ** Q) if y % 3), key=lambda y: comp[y])
    print(f"   the composed function vanishes on nonunits; max over units {float(comp[best]):.4f} at y = {best} (mod {3**Q}); mean over units {float(m * Fraction(3, 2)):.4f}")
    # residue pullback: choose the seed M = R_j in the class y = best mod 3^Q with S(M) = 1, and compute the ACTUAL
    # weighted history sum: sum over w in V admitting an inverse history from M of 3^d 2^-A * rho_k(source mod 3^k)
    j = next(j for j in range(1, 3 ** Q + 1) if ((4 ** j - 1) // 3) % 3 ** Q == best)
    M = (4 ** j - 1) // 3
    actual = Fraction(0)
    n_hist = 0
    for w in V:
        hist = inverse_history(M, w)
        if hist is None:
            continue
        n_hist += 1
        src = hist[-1]
        actual += Fraction(3 ** len(w), 2 ** sum(w)) * rk[src % 3 ** k]
    assert actual == comp[best]
    print(f"   seed M = R_{j} = {M} (S(M) = 1) in the best class: actual weighted history sum over {n_hist} admitted words = {float(actual):.4f} = composed function value (exact)")


# ----------------------------------------------------------------------------
# P6: deterministic spread of endpoints (3.9) at small scale
# ----------------------------------------------------------------------------

def part6():
    print("\n== P6: endpoints of histories with fixed (depth, valuation sum) in one residue class mod 3^q are spaced >= 3^q ==")
    M = (4 ** 40 - 1) // 3  # a seed with S(M) = 1
    words = list(product(range(1, 6), repeat=3))  # depth 3, valuations 1..5
    byDA = defaultdict(list)
    for w in words:
        h = inverse_history(M, w)
        if h is None:
            continue
        byDA[(3, sum(w))].append((h[-1], w))
    worst = 0
    for (D, A), lst in byDA.items():
        W = Fraction(3 ** D, 2 ** A)
        ends = sorted(set(x for x, _ in lst))
        # all W*x lie in an interval [M - beta_max, M]
        vals = [W * x for x in ends]
        assert all(v <= M for v in vals)
        spread = float(M - min(vals))
        worst = max(worst, spread)
        for q in (1, 2, 3):
            cls = defaultdict(list)
            for x in ends:
                cls[x % 3 ** q].append(x)
            for c, xs in cls.items():
                xs.sort()
                assert all(b - a >= 3 ** q for a, b in zip(xs, xs[1:]))
        # equal endpoints at fixed depth have the same word
        seen = {}
        for x, w in lst:
            assert seen.setdefault(x, w) == w
    print(f"   seed R_40, depth 3, valuations 1..5: {sum(len(v) for v in byDA.values())} histories in {len(byDA)} (D, A) classes; W x <= M always, max offset M - W x = {worst:.3f}; same-class endpoints spaced >= 3^q (q <= 3); equal endpoints share their word")


# ----------------------------------------------------------------------------
# P7: the time constant along typical inverse histories; empirical density of tau <= 10.46 log n
# ----------------------------------------------------------------------------

def part7():
    print("\n== P7: ordinary-time constant along inverse histories from seeds; empirical density of tau(n) <= (523/50) log n ==")
    rng = random.Random(3)

    def build(seed_j, depth, central):
        """Inverse history from R_j with unit intermediate sources. central=True: valuations 2 at even-parity
        steps and 1 or 3 (equally likely) at odd-parity steps, so the valuation sum is about 2*depth (the paper's
        central words); central=False: geometric valuations on the admitted parity (mean about 2.17)."""
        x = (4 ** seed_j - 1) // 3
        vals = []
        for _ in range(depth):
            par = 1 if (2 * x) % 3 == 1 else 0  # odd valuations admitted iff 2x = 1 mod 3
            cands = [1, 3, 5, 7] if par else [2, 4, 6, 8]
            if central:
                order = ([1, 3] if rng.random() < 0.5 else [3, 1]) + [5, 7] if par else [2, 4, 6, 8]
            else:
                a0 = cands[0]
                k = 0
                while rng.random() < 0.25:
                    k += 1
                order = [a0 + 2 * k] + [c for c in cands if c != a0 + 2 * k]
            for a in order:
                x2 = (2 ** a * x - 1) // 3
                if x2 % 3 != 0:
                    break
            else:
                return None, None
            x = x2
            vals.append(a)
        return x, vals

    # Section 7 of the paper: for a history of d inverse steps with total valuation A from the seed M,
    # tau(x) = d + A + tau(M) exactly (ordinary steps), and A log 2 - d log 3 <= log x - log M <= A log 2 - d log 3 + d log(1 + 1/(3 min x));
    # so tau(x)/log x -> (d + A)/(A log 2 - d log 3) as the seed's share vanishes, = 3/log(4/3) = 10.43 at A = 2d.
    for central in (True, False):
        ratios, limits, Aover = [], [], []
        for trial in range(300):
            j = rng.choice([2, 4, 5, 7, 8, 10, 11, 13])
            M = (4 ** j - 1) // 3
            depth = rng.randint(40, 120)
            x, vals = build(j, depth, central)
            if x is None or len(vals) < 30:
                continue
            d, A = len(vals), sum(vals)
            assert tau(x) == d + A + tau(M)
            # the orbit values x_0 = x, ..., x_(d-1) before M
            y, xs = x, []
            for _ in range(d):
                xs.append(y)
                y, _a = S(y)
            assert y == M
            # M = x 3^d prod(1 + 1/(3 x_j)) / 2^A, so log x - log M = A log 2 - d log 3 - sum log(1 + 1/(3 x_j))
            top = A * math.log(2) - d * math.log(3)
            bot = top - sum(math.log(1 + 1 / (3 * z)) for z in xs)
            assert bot - 1e-9 <= math.log(x) - math.log(M) <= top + 1e-9
            ratios.append(tau(x) / math.log(x))
            limits.append((d + A) / (A * math.log(2) - d * math.log(3)))
            Aover.append(A / d)
        print(f"   {'central-type' if central else 'geometric-parity'} histories, {len(ratios)} samples: mean valuation {sum(Aover)/len(Aover):.3f}; "
              f"tau(x) = d + A + tau(M) and the log bounds hold exactly; tau/log x mean {sum(ratios)/len(ratios):.2f} (seed-diluted), "
              f"seed-free constant (d+A)/(A log 2 - d log 3) mean {sum(limits)/len(limits):.2f}")
    print(f"   at A = 2d the seed-free constant is 3/log(4/3) = {3/LOG43:.3f}; the paper's 523/50 = 10.46 is this constant with rounding losses")
    N = 10 ** 6
    ok = 0
    for n in range(2, N):
        if tau(n) <= 523 / 50 * math.log(n):
            ok += 1
    print(f"   n < {N}: fraction with tau(n) <= 10.46 log n: {ok/(N-2):.4f} (the theorem asserts a positive lower density only beyond an astronomical X_0)")


# ----------------------------------------------------------------------------
# P8: the size of the constants
# ----------------------------------------------------------------------------

def part8():
    print("\n== P8: the paper's constants, in bit lengths ==")
    A_, E_ = 6409, 2170
    # log2 of D_exp = 2^(8192 A* 2^(3E*)): log2 D_exp = 8192 * A* * 2^(3E*)
    log2_Dexp_log2 = math.log2(8192 * A_) + 3 * E_  # log2(log2 D_exp)
    print(f"   D_exp = 2^(8192 A* 2^(3E*)): log2(log2 D_exp) = {log2_Dexp_log2:.1f}, i.e. D_exp ~ 2^(2^{log2_Dexp_log2:.0f})")
    # C* = (32 A* D*)^A*  ->  log2 C* = A* (5 + log2 A* + log2 D*) ~ A* * 2^6523
    print(f"   C* = (32 A* D*)^A*: log2(log2 C*) ~ {math.log2(A_) + log2_Dexp_log2:.1f}; C ~ 2 C* 20^6409")
    # F = 2^467 b^(16b) (C+1) with b = 280; N = 20000 (ceil(log2 F) + 64): log2 N ~ log2(20000) + log2(log2 F)
    print(f"   N = 20000 (ceil(log2 F) + 64) with log2 F ~ log2 C: log2 N ~ {math.log2(20000) + math.log2(A_) + log2_Dexp_log2:.1f}")
    print("   beta(N) = 280 (1.01)^N, q ~ 1.6 sum beta, M = 4^(2b+1+3q), c = 3/(256 M^2 m): c^-1 exceeds 2^(2^(2^6535)); X_0 = 32(2^B M + 1) likewise. No numerical range.")


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
    part7()
    part8()
    print("\nDONE")
