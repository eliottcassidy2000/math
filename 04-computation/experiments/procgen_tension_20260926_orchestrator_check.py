#!/usr/bin/env python3
"""Orchestrator audit of lane `tension` (ranks = tensions, Theorem M, the
numerator bound), written from the note's statements; the lane's scripts were
not read.

  1. Theorem D(ii): min over h of max over edges [w(s) + h(t) - h(s)] equals
     rho_max log 3 - log 2, by an LP (scipy linprog) against exact Karp, for all
     strategies at levels 2-4 and random ones at level 5.
  2. Corollary K: for a random bank of 2-adic centers and a periodic h, the integer
     shadow of the expanding orbit of 1^(L-1)0 makes R(T^L n) - R(n) = log(T^L n/n) > 0.
  3. Theorem M(v): the words 1100 (k = 4) and 11100 (k = 5, 6, 7) are maximal-density
     cycles of G_(sigma_k) that are not Christoffel words.
  4. Q3 refutation and the exact maxima: max rho_max over class-(i) strategies is
     3/5 at level 4 and 5/8 at level 5 (exhaustive), exceeding F_4 = 1/2, F_5 = 3/5;
     and it respects the numerator bound G_(2^(k-2)).
"""
import math, random, itertools
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog

C3 = math.log(2) / math.log(3)


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def targets(k, sig):
    N, H = 1 << k, 1 << (k - 1)
    return [((s // 2) % H) if s % 2 == 0 else (((3 * s + sig[s // 2]) // 2) % H) for s in range(N)], H


def karp(k, sig):
    N = 1 << k
    t, H = targets(k, sig)
    NEG = -10 ** 9
    D = [[0] * N] + [[NEG] * N for _ in range(N)]
    for j in range(1, N + 1):
        Dj, Dp = D[j], D[j - 1]
        for u in range(N):
            if Dp[u] > NEG:
                val = Dp[u] + (u & 1)
                for v in (t[u], t[u] + H):
                    if val > Dj[v]:
                        Dj[v] = val
    best = None
    for v in range(N):
        if D[N][v] <= NEG:
            continue
        worst = min(Fraction(D[N][v] - D[j][v], N - j) for j in range(N) if D[j][v] > NEG)
        if best is None or worst > best:
            best = worst
    return best


def lp_defect(k, sig):
    N = 1 << k
    t, H = targets(k, sig)
    w = [math.log(1.5) if s % 2 else -math.log(2) for s in range(N)]
    # variables: h_0..h_(N-1), tau ; minimize tau s.t. w(s) + h(t) - h(s) - tau <= 0
    A, b = [], []
    for s in range(N):
        for tt in (t[s], t[s] + H):
            row = [0.0] * (N + 1)
            row[tt] += 1.0
            row[s] -= 1.0
            row[N] = -1.0
            A.append(row); b.append(-w[s])
    c = [0.0] * N + [1.0]
    bounds = [(None, None)] * N + [(None, None)]
    bounds[0] = (0, 0)
    res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
    assert res.status == 0
    return res.x[N]


print("1. Theorem D(ii): least rank defect = max cycle mean")
random.seed(7)
cnt = 0
for k in (2, 3, 4):
    H = 1 << (k - 1)
    for sig in itertools.product((1, -1), repeat=H):
        rho = karp(k, list(sig))
        lam = float(rho) * math.log(3) - math.log(2)
        d = lp_defect(k, list(sig))
        assert abs(d - lam) < 1e-9, (k, sig, d, lam)
        cnt += 1
for _ in range(40):
    sig = [random.choice((1, -1)) for _ in range(16)]
    rho = karp(5, sig)
    assert abs(lp_defect(5, sig) - (float(rho) * math.log(3) - math.log(2))) < 1e-9
    cnt += 1
check(True, f"LP minimum of max_edge[w + dh] = rho_max log 3 - log 2 on {cnt} strategies (all of levels 2-4, 40 random at level 5)")
check(abs(lp_defect(6, [1] * 32) - math.log(1.5)) < 1e-9, "Collatz at level 6: least defect log(3/2) (the loop at -1)")

print("2. Corollary K: finite dyadic banks fail")


def v2(x):
    if x == 0:
        return 10 ** 9
    v = 0
    while x % 2 == 0:
        x //= 2; v += 1
    return v


def collatz(n):
    return (3 * n + 1) // 2 if n % 2 else n // 2


random.seed(11)
for L in (5, 7, 9):
    # periodic point x = c/(2^L - 3^(L-1)) of the word 1^(L-1) 0, as a 2-adic integer mod 2^M
    c = 0
    for j in range(L - 1):
        c = 3 * c + 2 ** j
    D = 2 ** L - 3 ** (L - 1)
    M = 200
    x0 = (c * pow(D, -1, 2 ** M)) % 2 ** M
    # a bank of 12 random rational 2-adic centers (numerator/odd denominator), and a random periodic h
    K = 6
    hvals = [random.uniform(-5, 5) for _ in range(2 ** K)]
    bank = []
    for _ in range(12):
        num, den = random.randint(-10 ** 6, -1), random.choice([1, 3, 5, 7, 9, 11])
        bank.append(((num * pow(den, -1, 2 ** M)) % 2 ** M, random.uniform(-3, 3)))
    V = max(v2((x0 - b) % 2 ** M) for b, _ in bank)
    assert V < M - 50
    N = K + L + V + 1
    n = x0 % 2 ** N + 2 ** N * random.randint(10 ** 30, 10 ** 31)
    def R(m):
        return math.log(m) + hvals[m % 2 ** K] + sum(cc * v2((m - b) % 2 ** M) for b, cc in bank)
    y = n
    for _ in range(L):
        y = collatz(y)
    assert y % 2 ** (K + V + 1) == n % 2 ** (K + V + 1)
    inc = R(y) - R(n)
    assert abs(inc - math.log(y / n)) < 1e-9 and inc > 0, (L, inc)
    print(f"   L={L}: V={V}, R(T^L n) - R(n) = {inc:.6f} = log(T^L n/n) -> log(3^(L-1)/2^L) = {math.log(3**(L-1)/2**L):.6f}")
check(True, "for L = 5, 7, 9 the integer shadow of the orbit of 1^(L-1)0 increases every such rank over one period")

print("3. Theorem M(v): non-Christoffel maximizers of sigma_k")


def bad_mask(k):
    M = 1 << k
    bad = []
    for r in range(M):
        x, a, ok = r, 0, True
        for j in range(1, k + 1):
            if x % 2:
                a += 1
            x = (3 * x + 1) // 2 if x % 2 else x // 2
            if not 3 ** a > 2 ** j:
                ok = False
                break
        bad.append(ok)
    return bad


def sigma_k(k):
    b = bad_mask(k)
    return [-1 if b[r] else 1 for r in range(1, 1 << k, 2)]


def word_cycle_in_graph(k, sig, word):
    """periodic Collatz point of the word; check its residues form a closed walk of G_sigma with + signs."""
    p, a = len(word), sum(word)
    c = 0
    for j, bit in enumerate(word):
        if bit:
            c = 3 * c + 2 ** j
    x = Fraction(c, 2 ** p - 3 ** a)
    N = 1 << k
    t, H = targets(k, sig)
    res = []
    y = x
    for j in range(p):
        r = (y.numerator * pow(y.denominator, -1, N)) % N
        res.append(r)
        assert (y.numerator % 2) == word[j]
        y = (3 * y + 1) / 2 if y.numerator % 2 else y / 2
    assert y == x
    for j in range(p):
        s, s2 = res[j], res[(j + 1) % p]
        if s % 2 and sig[s // 2] != 1:
            return False
        if s2 not in (t[s], t[s] + H):
            return False
    return True


cases = [(4, [1, 1, 0, 0], Fraction(1, 2))] + [(k, [1, 1, 1, 0, 0], Fraction(3, 5)) for k in (5, 6, 7)]
for k, word, dens in cases:
    sig = sigma_k(k)
    assert word_cycle_in_graph(k, sig, word), (k, word)
    assert karp(k, sig) == dens == Fraction(sum(word), len(word))
check(True, "1100 is a maximal cycle of sigma_4 (density 1/2), 11100 of sigma_5..sigma_7 (density 3/5); neither is Christoffel")

print("4. max rho_max over class (i): exceeds F_k (Q3 refuted)")


def F_k(k):
    return max(Fraction(a, d) for d in range(1, k + 1) for a in range(d + 1) if 3 ** a < 2 ** d)


def G_num(nmax):
    return max(Fraction(a, d) for a in range(1, nmax + 1) for d in range(a, 3 * a) if 3 ** a < 2 ** d)


for k in (3, 4):
    H = 1 << (k - 1)
    vals = []
    for sig in itertools.product((1, -1), repeat=H):
        rho = karp(k, list(sig))
        if rho < Fraction(C3).limit_denominator(10 ** 9):
            vals.append(rho)
    print(f"   k={k}: {len(vals)} class-(i) strategies, max rho_max = {max(vals)}, F_k = {F_k(k)}, G_(2^(k-2)) = {G_num(2 ** (k - 2))}")
    assert max(vals) <= G_num(2 ** (k - 2))
check(max(vals) == Fraction(3, 5) and F_k(4) == Fraction(1, 2), "level 4: a class-(i) strategy has rho_max = 3/5 > F_4 = 1/2 (= G_2, the numerator bound)")
