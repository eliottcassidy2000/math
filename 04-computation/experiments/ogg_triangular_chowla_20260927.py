#!/usr/bin/env python3
"""Ogg's fifteen primes by the genus formula and by supersingular j-invariants;
integer recurrences (378 = T_27, 1/30 and 1/42, Hurwitz orders); the nodal cubic
y^2 = x^2 (x+1) against triangular numbers (Pell family, an elliptic curve of
positive rank, the point (3, 6)); the Chowla map s'(n) = sigma(n) - n - 1
(aliquot with 1 and n discounted): parity alternation, termination, cycles;
the Paley heptagon inside PSL(2,7).

Session: opus, collatz-poset-dag-20260927 (S17), 2026-09-27.  Exact integer
arithmetic; numpy sieve for sigma; sympy for factorizations beyond the sieve.

Run: python 04-computation/experiments/ogg_triangular_chowla_20260927.py
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations

import numpy as np
from sympy import divisor_sigma, factorint, isprime, primerange

OGG = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]


# ----------------------------------------------------------------------------
# P1: Ogg's list from the genus formula (class numbers) and from supersingular j
# ----------------------------------------------------------------------------

def legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def class_number(D: int) -> int:
    """h(D) = number of primitive reduced positive definite forms of discriminant D < 0."""
    assert D < 0 and D % 4 in (0, 1)
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if c == a and b < 0:
                continue
            if math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    return h


def genus_X0(p: int) -> Fraction:
    mu = p + 1
    # Kronecker conventions: (-4/2) = 0, (-3/2) = -1, (-3/3) = 0, (-4/3) = -1
    nu2 = 1 if p == 2 else 1 + legendre(-1, p)
    nu3 = 0 if p == 2 else (1 if p == 3 else 1 + legendre(-3, p))
    return 1 + Fraction(mu, 12) - Fraction(nu2, 4) - Fraction(nu3, 3) - 1


def fixed_points_wp(p: int) -> int:
    if p == 2:
        return class_number(-8) + class_number(-4)
    if p % 4 == 1:
        return class_number(-4 * p)
    return class_number(-4 * p) + class_number(-p)


def genus_X0plus(p: int) -> Fraction:
    g = genus_X0(p)
    f = fixed_points_wp(p)
    return (2 * g + 2 - f) / Fraction(4)


def supersingular_j_in_Fp(p: int) -> tuple[int, int]:
    """(number of supersingular j-invariants in char p, number of them in F_p),
    via the Hasse polynomial on the Legendre lambda-line, roots in F_{p^2}."""
    m = (p - 1) // 2
    coeffs = [math.comb(m, k) ** 2 % p for k in range(m + 1)]
    # F_{p^2} = F_p[s]/(s^2 - d), d a non-residue
    d = next(x for x in range(2, p) if legendre(x, p) == -1)

    def mul(u, v):
        return ((u[0] * v[0] + u[1] * v[1] * d) % p, (u[0] * v[1] + u[1] * v[0]) % p)

    def inv(u):
        # (a + b s)^-1 = (a - b s)/(a^2 - b^2 d)
        n = (u[0] * u[0] - u[1] * u[1] * d) % p
        ni = pow(n, p - 2, p)
        return (u[0] * ni % p, (-u[1] * ni) % p)

    def H(lam):
        acc = (0, 0)
        pw = (1, 0)
        for c in coeffs:
            acc = ((acc[0] + c * pw[0]) % p, (acc[1] + c * pw[1]) % p)
            pw = mul(pw, lam)
        return acc

    js = set()
    for a in range(p):
        for b in range(p):
            lam = (a, b)
            if lam in ((0, 0), (1, 0)):
                continue
            if H(lam) != (0, 0):
                continue
            # j = 256 (lam^2 - lam + 1)^3 / (lam^2 (lam - 1)^2)
            l2 = mul(lam, lam)
            num = ((l2[0] - lam[0] + 1) % p, (l2[1] - lam[1]) % p)
            num3 = mul(mul(num, num), num)
            lm1 = ((lam[0] - 1) % p, lam[1])
            den = mul(l2, mul(lm1, lm1))
            j = mul(num3, inv(den))
            j = (256 * j[0] % p, 256 * j[1] % p)
            js.add(j)
    return len(js), sum(1 for j in js if j[1] == 0)


def part1():
    print("== P1: Ogg's fifteen primes ==")
    rows = []
    for p in primerange(2, 120):
        g = genus_X0(p)
        f = fixed_points_wp(p)
        gp = genus_X0plus(p)
        assert g.denominator == 1 and gp.denominator == 1, (p, g, gp)
        rows.append((p, int(g), f, int(gp)))
    genus0 = [p for p, g, f, gp in rows if gp == 0]
    print(f"   primes < 120 with g(X_0(p)^+) = 0 by the genus formula: {genus0}")
    assert genus0 == OGG
    print("   p: g(X_0(p)), fixed points of w_p (class numbers), g(X_0^+(p)) for p <= 100:")
    print("   " + "; ".join(f"{p}:{g},{f},{gp}" for p, g, f, gp in rows if p <= 100))
    # supersingular j-invariants: all in F_p iff p in OGG (p >= 5)
    ss = []
    for p in primerange(5, 100):
        n, nfp = supersingular_j_in_Fp(p)
        ss.append((p, n, nfp))
    all_in = [p for p, n, nfp in ss if n == nfp]
    print(f"   supersingular j-invariants (p from 5 to 97): all in F_p exactly for {all_in}")
    assert all_in == [p for p in OGG if p >= 5]
    print("   counts (p: total, in F_p): " + ", ".join(f"{p}:{n},{nfp}" for p, n, nfp in ss))
    # Paley half-row class number formula for p = 3 mod 4, p > 3
    ok = []
    for p in primerange(7, 200):
        if p % 4 != 3:
            continue
        S = sum(legendre(a, p) for a in range(1, (p + 1) // 2))
        h = Fraction(S, 2 - legendre(2, p))
        assert h == class_number(-p), (p, S, h, class_number(-p))
        ok.append(p)
    print(f"   h(-p) = (sum of (a/p), 0 < a < p/2)/(2 - (2/p)) verified for p = 3 mod 4 up to 199 ({len(ok)} primes)")
    print("   h(-p) for the Ogg primes = 3 mod 4: " + ", ".join(f"{p}:{class_number(-p)}" for p in OGG if p % 4 == 3 and p > 3))


# ----------------------------------------------------------------------------
# P2: integer recurrences
# ----------------------------------------------------------------------------

def tri_index(n: int):
    r = math.isqrt(8 * n + 1)
    return (r - 1) // 2 if r * r == 8 * n + 1 else None


def part2():
    print("\n== P2: integers around the objects ==")
    s = sum(OGG)
    print(f"   sum of Ogg primes = {s} = T_{tri_index(s)} = {factorint(s)}; primes <= 31 sum {sum(p for p in OGG if p <= 31)}, the four sporadic {[p for p in OGG if p > 31]} sum {sum(p for p in OGG if p > 31)}")
    monster = {2: 46, 3: 20, 5: 9, 7: 6, 11: 2, 13: 3, 17: 1, 19: 1, 23: 1, 29: 1, 31: 1, 41: 1, 47: 1, 59: 1, 71: 1}
    sopfr = sum(p * e for p, e in monster.items())
    print(f"   |M| prime factors with multiplicity sum to {sopfr} = {factorint(sopfr)}; distinct-prime sum {s}")
    missing = [p for p in primerange(2, 72) if p not in OGG]
    print(f"   primes < 72 not in the list: {missing}; the list = all primes <= 31 plus {[p for p in OGG if p > 31]}")
    for name, tr in (("(2,3,3)", (2, 3, 3)), ("(2,3,4)", (2, 3, 4)), ("(2,3,5)", (2, 3, 5)), ("(2,3,6)", (2, 3, 6)), ("(2,3,7)", (2, 3, 7))):
        x = sum(Fraction(1, e) for e in tr) - 1
        print(f"   triangle {name}: 1/2+1/3+1/q - 1 = {x}" + (f"; group order 2/x = {2 / x}" if x > 0 else (" (Euclidean)" if x == 0 else f"; Hurwitz bound 1/|x| = {1 / abs(x)}")))
    for n in (42, 168, 504, 1092, 30, 378, 24, 120, 210):
        print(f"   {n} = {factorint(n)}" + (f" = T_{tri_index(n)}" if tri_index(n) else ""))
    for q in (7, 8, 13):
        print(f"   |PSL(2,{q})| = {q * (q * q - 1) // (1 if q % 2 == 0 else 2)}")
    print(f"   Hurwitz genera 3, 7, 14, 17, 118: {[factorint(g) for g in (3, 7, 14, 17, 118)]}; 84(g-1) = {[84 * (g - 1) for g in (3, 7, 14)]}")


# ----------------------------------------------------------------------------
# P3: the nodal cubic and triangular numbers
# ----------------------------------------------------------------------------

def is_square(n: int) -> bool:
    return n >= 0 and math.isqrt(n) ** 2 == n


def ec_add(P, Q, a, b):
    """Addition on y^2 = x^3 + a x + b over Q; None is the identity."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and y1 == -y2:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + a) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - x1 - x2
    y3 = lam * (x1 - x3) - y1
    return (x3, y3)


def ec_mul(n, P, a, b):
    R = None
    Q = P
    if n < 0:
        n, Q = -n, (P[0], -P[1])
    while n:
        if n & 1:
            R = ec_add(R, Q, a, b)
        Q = ec_add(Q, Q, a, b)
        n >>= 1
    return R


def part3():
    print("\n== P3: the nodal cubic y^2 = x^2 (x + 1) and triangular numbers ==")
    print("   rational points: x = t^2 - 1, y = t (t^2 - 1) = 6 C(t+1, 3); group of nonsingular points = G_m via (y+x)/(y-x) = (t+1)/(t-1)")
    # x triangular: t^2 - 1 = T_m  <=>  (2m+1)^2 - 8 t^2 = -7
    xs = [t for t in range(1, 200001) if tri_index(t * t - 1) is not None]
    print(f"   t <= 2*10^5 with x = t^2 - 1 triangular (Pell u^2 - 8t^2 = -7): t = {xs[:12]} ...; x = {[t*t-1 for t in xs[:8]]} = T_m, m = {[tri_index(t*t-1) for t in xs[:8]]}")
    # y triangular: t(t^2-1) = T_k  <=>  8 t^3 - 8 t + 1 = (2k+1)^2  <=>  Y^2 = X^3 - 4X + 1 with X = 2t
    ys = [t for t in range(1, 1000001) if is_square(8 * t ** 3 - 8 * t + 1)]
    print(f"   t <= 10^6 with y = t(t^2-1) triangular: t = {ys}; y = {[t*(t*t-1) for t in ys]} = T_k, k = {[tri_index(t*(t*t-1)) for t in ys]}")
    both = [t for t in ys if tri_index(t * t - 1) is not None]
    print(f"   both coordinates triangular (t <= 10^6): t = {both}, points {[(t*t-1, t*(t*t-1)) for t in both]}")
    # the curve E: Y^2 = X^3 - 4X + 1: integer points |X| <= 10^6 and the group structure of the ones found
    a, b = -4, 1
    pts = []
    for X in range(-2, 10 ** 6 + 1):
        v = X ** 3 - 4 * X + 1
        if is_square(v):
            pts.append((X, math.isqrt(v)))
    print(f"   E: Y^2 = X^3 - 4X + 1, discriminant {-16 * (4 * a ** 3 + 27 * b * b)} = {factorint(-16 * (4 * a ** 3 + 27 * b * b))}; integer points with -2 <= X <= 10^6: {pts}")
    P = (Fraction(0), Fraction(1))
    Q = (Fraction(2), Fraction(1))
    span = {}
    for n in range(-12, 13):
        for m in range(-12, 13):
            R = ec_add(ec_mul(n, P, a, b), ec_mul(m, Q, a, b), a, b)
            if R is not None and R[0].denominator == 1 and R[1].denominator == 1:
                span[(int(R[0]), abs(int(R[1])))] = (n, m)
    expl = {pt: span.get(pt) for pt in pts}
    print(f"   integer points as nP + mQ with P = (0,1), Q = (2,1), |n|,|m| <= 12: {expl}")
    # is Q a multiple of P or P of Q?  check small multiples
    multsP = {(int(R[0]), abs(int(R[1]))) for n in range(1, 40) for R in [ec_mul(n, P, a, b)] if R is not None and R[0].denominator == 1}
    multsQ = {(int(R[0]), abs(int(R[1]))) for n in range(1, 40) for R in [ec_mul(n, Q, a, b)] if R is not None and R[0].denominator == 1}
    print(f"   integer multiples of P (n < 40): {sorted(multsP)}; of Q: {sorted(multsQ)}")
    # square pyramidal: sum of squares a square only at n = 24 (checked to 10^6)
    cann = [n for n in range(1, 10 ** 6) if is_square(n * (n + 1) * (2 * n + 1) // 6)]
    print(f"   n <= 10^6 with 1^2 + ... + n^2 a square: {cann} (Lucas's cannonball problem; 24 -> 70^2, the Leech vector)")


# ----------------------------------------------------------------------------
# P4: the Chowla map s'(n) = sigma(n) - n - 1
# ----------------------------------------------------------------------------

def sigma_sieve(N: int) -> np.ndarray:
    s = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N + 1):
        s[d::d] += d
    return s


def v2(x: int) -> int:
    return (x & -x).bit_length() - 1 if x else 10 ** 9


def part4(N: int = 10 ** 6, SIEVE: int = 4 * 10 ** 6):
    print("\n== P4: the Chowla map s'(n) = sigma(n) - n - 1 (aliquot with 1 and n discounted) ==")
    sig = sigma_sieve(SIEVE)
    cache = {}

    def sp(n: int) -> int:
        if n <= SIEVE:
            return int(sig[n]) - n - 1
        if n not in cache:
            cache[n] = int(divisor_sigma(n)) - n - 1
        return cache[n]

    # parity theorem: for n = 2^a m even, s'(n) odd iff m not a square; for odd n, s'(n) even iff n not a square
    for n in range(2, 300001):
        s = sp(n)
        if n % 2 == 0:
            m = n >> v2(n)
            assert (s % 2 == 1) == (not is_square(m)), n
        else:
            assert (s % 2 == 0) == (not is_square(n)), n
    print("   parity theorem checked for n <= 3*10^5: s'(even) is odd unless the odd part is a square; s'(odd) is even unless n is a square")
    ends = Counter()
    cycles = {}
    maxlen = (0, 0)
    maxpeak = (0, 0)
    lengths = []
    beyond = 0
    for n in range(2, N + 1):
        seen = {}
        x = n
        k = 0
        peak = n
        while True:
            if x == 0:
                ends['0'] += 1
                break
            if x == 1:
                ends['1'] += 1
                break
            if x in seen:
                cyc = []
                y = x
                while True:
                    cyc.append(y)
                    y = sp(y)
                    if y == x:
                        break
                key = tuple(sorted(cyc))
                cycles[key] = cycles.get(key, 0) + 1
                ends[f'cycle{len(cyc)}'] += 1
                break
            seen[x] = k
            x = sp(x)
            if x > SIEVE:
                beyond += 1
            k += 1
            peak = max(peak, x)
            assert k < 10000
        lengths.append(k)
        if k > maxlen[0]:
            maxlen = (k, n)
        if peak / n > maxpeak[0]:
            maxpeak = (peak / n, n)
    print(f"   n <= {N}: endings {dict(ends)}; longest {maxlen[0]} steps (n = {maxlen[1]}); largest peak/n {maxpeak[0]:.1f} (n = {maxpeak[1]}); mean length {sum(lengths)/len(lengths):.2f}; values beyond the sieve {beyond}")
    print(f"   cycles reached (members: number of starts): {[(k, v) for k, v in sorted(cycles.items())]}")
    # drift by parity, and parity persistence, on n <= N
    ev, od = [], []
    for n in range(2, N + 1):
        s = sp(n)
        if s > 0:
            (ev if n % 2 == 0 else od).append(math.log2(s / n))
    print(f"   mean log2(s'(n)/n): even n {sum(ev)/len(ev):+.3f}, odd n {sum(od)/len(od):+.3f} (aliquot s(n): even -0.048, odd -5.34)")
    for M in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
        same = sum(1 for n in range(M, 2 * M) if sp(n) > 0 and sp(n) % 2 == n % 2)
        tot = sum(1 for n in range(M, 2 * M) if sp(n) > 0)
        print(f"   parity persists under s' in [M, 2M), M = {M}: {same/tot:.5f} (squares' share ~ {1/math.sqrt(M):.5f})")
    # odd abundant numbers: the only odd growth
    oa = [n for n in range(1, 100001, 2) if int(sig[n]) > 2 * n + 1]
    print(f"   odd n <= 10^5 with s'(n) > n: {len(oa)}, first {oa[:6]}")


# ----------------------------------------------------------------------------
# P5: the Paley heptagon inside PSL(2,7)
# ----------------------------------------------------------------------------

def part5():
    print("\n== P5: Paley heptagon and PSL(2,7) ==")
    p = 7
    qr = {(x * x) % p for x in range(1, p)}
    arcs = {(i, j) for i in range(p) for j in range(p) if i != j and (j - i) % p in qr}
    auts = 0
    for perm in permutations(range(p)):
        if all((perm[i], perm[j]) in arcs for (i, j) in arcs):
            auts += 1
    print(f"   |Aut(Paley tournament P_7)| = {auts} = F_21; |PSL(2,7)| = 168 = 8 * 21; the Klein quartic's group contains the heptagon's as the Sylow-7 normalizer")


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    print("\nDONE")
