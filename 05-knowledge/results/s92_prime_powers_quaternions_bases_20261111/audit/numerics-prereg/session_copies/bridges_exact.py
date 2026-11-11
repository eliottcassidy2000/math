#!/usr/bin/env python3
"""S92 (THM-5135 and the note, sections 4-8): exact checks of the bridges.  Everything is exact (integers or
Fractions) unless marked NUMERICAL.  Exit status 1 on failure.

S1  Multiplication table modulo odd m: the number of (a, b) with a b = t is N_m(t) = sum_{d | gcd(t, m)} d phi(m/d),
    and the additive self-convolution of the table (the F1 seam a1 a2 + u b1 b2, u a unit) has
    (N_m * N_m)(c) = m^3 L_m(c), L_m(c) = sum_{d | gcd(m, c)} d^-1 prod_{p | m/d} (1 - p^-2).
S2  Sums of 2r squares against the hyperbolic form of rank r modulo odd p (and p^2): equal class counts iff r is
    even or p = 1 (mod 4).
S3  The prime 2: four squares and a1 a2 + b1 b2 differ modulo 4 and 8; eight squares and the rank-4 form differ;
    the E8 form x^T C x / 2 (C the Cartan matrix) and the rank-4 form agree modulo 2, 4, 8.
S4  Global counts: r_4 = additive convolution of r_2 with itself = 8 sigma(n) (n odd); bicomplex norm count =
    multiplicative convolution of r_2 with itself, zero for n = 3 (mod 4); Hurwitz 24 sigma(n) (n odd); E8 shells
    240 sigma_3(n) = (8 pi^4/3) n^3 * sigma_3(n)/(zeta(4) n^3); r_8(n) = 16 sigma_3(n) (n odd).  NUMERICAL: Ingham's
    sum D(n) = sum_{m < n} d(m) d(n - m) against sigma(n).
S5  Two-port kernels of THM-5133: u^T G v with uniform states has the rank-2 law at every p not dividing det G.
S6  Central extensions of F_2^2 by F_2: 16 normalized cocycles, 8 classes, groups C2^3, C4 x C2 (3), D4 (3), Q8;
    the class is the set of units of square -1; commutative iff that set is even; antisymmetric sign tables are
    the 8 tournaments on the three points (6 transitive = D4 with the unit of square -1 in the middle, 2 cyclic =
    Q8), symmetric ones are signed triangles (R^4 and the bicomplex numbers C + C).
S7  Signed-integer dictionary: odd m >= 1 <-> n = (chi_4(m) m - 1)/4 in Z; Collatz step, tripling, products;
    labels of the arcs of the five known cycles of 3n + 1 on Z; spectra of the orders 3..6 as sets of integers.
S8  Primes with a balanced divisible state (p | 2^t + 1) = primes of even order of 2; first fan and first
    p-divisible state.
S9  Brahmagupta identity for self-joins and cross-joins.
S10 F1 in base beta, sheet s = 0 (NUMERICAL census, exact structure): triples whose smaller state is p-divisible
    lie in the class 2N + beta = 0 (mod p); remainder of the first-order law.
S11 Tao's digit condition for the packing polynomial 1 + alpha_1 z + alpha_2 z^2 of the joins.
"""
import itertools
import math
import sys
from fractions import Fraction as F

import numpy as np
from numba import njit

OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg)
    OK &= bool(cond)


def primes_of(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def phi(n):
    r = n
    for p in primes_of(n):
        r = r // p * (p - 1)
    return r


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def L(m, e):
    g = math.gcd(m, e % m) if e % m else m
    tot = F(0)
    for d in divisors(g):
        t = F(1, d)
        for p in primes_of(m // d):
            t *= 1 - F(1, p * p)
        tot += t
    return tot


def cconv(a, b):
    q = len(a)
    out = np.zeros(q, dtype=object)
    for i in range(q):
        if a[i]:
            for j in range(q):
                if b[j]:
                    out[(i + j) % q] += a[i] * b[j]
    return out


def dist_prod(q, u=1):
    d = np.zeros(q, dtype=object)
    for a in range(q):
        for b in range(q):
            d[(u * a * b) % q] += 1
    return d


def dist_sq(q):
    d = np.zeros(q, dtype=object)
    for a in range(q):
        d[(a * a) % q] += 1
    return d


def power_conv(d, r):
    out = d
    for _ in range(r - 1):
        out = cconv(out, d)
    return out


# ---------------------------------------------------------------- S1
def s1():
    print("S1 multiplication tables modulo odd m and their additive square")
    for m in (9, 15, 21, 25, 27, 45, 49, 63, 75, 81, 105):
        tab = dist_prod(m)
        good = all(int(tab[t]) == sum(d * phi(m // d) for d in divisors(math.gcd(t, m) if t else m)) for t in range(m))
        sq = cconv(tab, tab)
        good2 = all(F(int(sq[c])) == m ** 3 * L(m, c) for c in range(m))
        tw = cconv(tab, dist_prod(m, 2))
        tw4 = cconv(tab, dist_prod(m, 4))
        good3 = all(int(tw[c]) == int(sq[c]) == int(tw4[c]) for c in range(m))
        check(good and good2 and good3, f"m = {m:3d}: N_m(t) = sum d phi(m/d); (N_m * N_m)(c) = m^3 L_m(c); same with the units 2, 4")


# ---------------------------------------------------------------- S2
def s2():
    print("S2 sums of 2r squares against the hyperbolic form of rank r, class counts modulo q")
    for q in (3, 5, 7, 11, 13, 9, 25):
        p = primes_of(q)[0]
        sq, hy = dist_sq(q), dist_prod(q)
        for r in (1, 2, 3, 4):
            if q > 13 and r > 2:
                continue
            a = power_conv(sq, 2 * r)
            b = power_conv(hy, r)
            same = all(int(a[c]) == int(b[c]) for c in range(q))
            expect = (r % 2 == 0) or (p % 4 == 1)
            check(same == expect, f"q = {q:2d}, r = {r}: equal = {same} (expected {expect}: r even or p = 1 mod 4)")


# ---------------------------------------------------------------- S3
@njit(cache=True)
def e8_classes(C, q):
    cnt = np.zeros(q, np.int64)
    x = np.zeros(8, np.int64)
    total = q ** 8
    for idx in range(total):
        t = idx
        for i in range(8):
            x[i] = t % q
            t //= q
        s = 0
        for i in range(8):
            for j in range(8):
                s += x[i] * C[i, j] * x[j]
        cnt[(s // 2) % q] += 1
    return cnt


def s3():
    print("S3 the prime 2")
    for q in (4, 8):
        four = power_conv(dist_sq(q), 4)
        hyp2 = power_conv(dist_prod(q), 2)
        check(any(int(four[c]) != int(hyp2[c]) for c in range(q)),
              f"modulo {q}: four squares {[int(v) for v in four]} differ from a1 a2 + b1 b2 {[int(v) for v in hyp2]}")
    C = np.zeros((8, 8), np.int64)
    for i in range(8):
        C[i, i] = 2
    for (i, j) in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)):      # E8 Dynkin diagram
        C[i, j] = C[j, i] = -1
    det = round(np.linalg.det(C.astype(float)))
    check(det == 1, f"the E8 Cartan matrix has determinant {det}")
    for q in (2, 4, 8):
        e8 = e8_classes(C, q)
        hyp4 = power_conv(dist_prod(q), 4)
        eight = power_conv(dist_sq(q), 8)
        check(all(int(e8[c]) == int(hyp4[c]) for c in range(q)),
              f"modulo {q}: the E8 form and the hyperbolic form of rank 4 have the same class counts {[int(v) for v in e8]}")
        if q > 2:
            check(any(int(eight[c]) != int(hyp4[c]) for c in range(q)), f"modulo {q}: eight squares differ from both")


# ---------------------------------------------------------------- S4
def s4():
    print("S4 global counts")
    LIM = 400
    r2 = [0] * (LIM + 1)
    B = int(LIM ** 0.5) + 1
    for a in range(-B, B + 1):
        for b in range(-B, B + 1):
            n = a * a + b * b
            if n <= LIM:
                r2[n] += 1
    sig = [0] * (LIM + 1)
    sig3 = [0] * (LIM + 1)
    for d in range(1, LIM + 1):
        for k in range(d, LIM + 1, d):
            sig[k] += d
            sig3[k] += d ** 3
    r4 = [sum(r2[a] * r2[n - a] for a in range(n + 1)) for n in range(LIM + 1)]
    check(all(r4[n] == 8 * sig[n] for n in range(1, LIM + 1, 2)), f"r_4(n) = (r_2 additively convolved with r_2)(n) = 8 sigma(n) for odd n <= {LIM}")
    bic = [sum(r2[a] * r2[n // a] for a in range(1, n + 1) if n % a == 0) for n in range(LIM + 1)]
    check(all(bic[n] == 0 for n in range(3, LIM + 1, 4)) and bic[1] == 16 and bic[5] == 2 * 4 * 8,
          f"bicomplex norm |z1|^2 |z2|^2 = n on Z[i]^2: the multiplicative convolution of r_2 with itself; it is 0 for every n = 3 (mod 4), n <= {LIM}; values at 1, 5, 9, 13, 25: {[bic[n] for n in (1, 5, 9, 13, 25)]}")
    # Hurwitz quaternions of norm n: integer vectors and half-integer vectors
    def hurwitz(n):
        cnt = r4[n]
        B2 = int((4 * n) ** 0.5) + 1
        odd = [v for v in range(-B2, B2 + 1) if v % 2]
        sq = {}
        for a in odd:
            for b in odd:
                sq[a * a + b * b] = sq.get(a * a + b * b, 0) + 1
        cnt += sum(c * sq.get(4 * n - k, 0) for k, c in sq.items())
        return cnt
    check(all(hurwitz(n) == 24 * sig[n] for n in range(1, 60, 2)), "Hurwitz order: 24 sigma(n) elements of norm n for odd n < 60 (one third on the Lipschitz grid)")
    r8 = [sum(r4[a] * r4[n - a] for a in range(n + 1)) for n in range(41)]
    check(all(r8[n] == 16 * sig3[n] for n in range(1, 41, 2)), "r_8(n) = 16 sigma_3(n) for odd n <= 39")
    check(r8[2] != 16 * sig3[2], f"r_8(2) = {r8[2]} is not 16 sigma_3(2) = {16 * sig3[2]} (the prime 2 is special for Z^8)")
    # E8 = D8+: integer vectors with even sum, and half-integer vectors with even sum; shells of norm 2n
    def e8_shell(n):
        tot = 0
        target = 2 * n
        vals = list(range(-int(target ** 0.5) - 1, int(target ** 0.5) + 2))

        def rec(k, rem, par, half):
            if k == 0:
                return 1 if (rem == 0 and par == 0) else 0
            s = 0
            for v in vals:
                w = (2 * v + 1) if half else 2 * v          # twice the coordinate
                if w * w <= rem:
                    s += rec(k - 1, rem - w * w, (par + (v if not half else 0)) % 2, half)
            return s
        tot += rec(8, 4 * target, 0, False)
        # half-integer vectors: coordinates (2v+1)/2, sum even <=> sum of (2v+1) = 0 mod 4
        def rec2(k, rem, ssum):
            if k == 0:
                return 1 if (rem == 0 and ssum % 4 == 0) else 0
            s = 0
            for v in vals:
                w = 2 * v + 1
                if w * w <= rem:
                    s += rec2(k - 1, rem - w * w, ssum + w)
            return s
        tot += rec2(8, 4 * target, 0)
        return tot
    shells = [e8_shell(n) for n in (1, 2, 3)]
    check(shells == [240 * sig3[n] for n in (1, 2, 3)], f"E8 shells of norm 2n, n = 1, 2, 3: {shells} = 240 sigma_3(n)")
    zeta4 = math.pi ** 4 / 90
    check(all(abs((8 * math.pi ** 4 / 3) * n ** 3 * sig3[n] / (zeta4 * n ** 3) - 240 * sig3[n]) < 1e-6 for n in (1, 2, 3, 6, 9)),
          "240 sigma_3(n) = real density (8 pi^4/3) n^3 times the rank-4 series sigma_3(n)/(zeta(4) n^3) at ALL primes")
    check(all(abs(math.pi ** 2 * n * (8 / math.pi ** 2) * sig[n] / n - r4[n]) < 1e-6 for n in range(1, 100, 2)),
          "r_4(n) = real density pi^2 n times (8/pi^2) sigma(n)/n (odd primes), 2-adic density 1, for odd n")
    # Ingham (NUMERICAL)
    N0, W = 20000, 400
    top = N0 + W
    dv = np.zeros(top + 1, np.int64)
    for d in range(1, top + 1):
        dv[d::d] += 1
    sg = np.zeros(top + 1, np.int64)
    for d in range(1, top + 1):
        sg[d::d] += d
    ns = np.arange(N0, top)
    D = np.array([int(np.dot(dv[1:n], dv[n - 1:0:-1])) for n in ns], dtype=float)
    ratio = D / ((6 / math.pi ** 2) * sg[ns] * np.log(ns) ** 2)
    corr = np.corrcoef(D / ns, sg[ns] / ns)[0, 1]
    print(f"  NUMERICAL Ingham: n in [{N0}, {top}): D(n) / ((6/pi^2) sigma(n) ln^2 n) has mean {ratio.mean():.3f}, min {ratio.min():.3f}, "
          f"max {ratio.max():.3f} (lower-order terms in ln n); correlation of D(n)/n with sigma(n)/n: {corr:.4f}")
    check(corr > 0.9, "NUMERICAL illustration of Ingham's theorem (CITED): D(n)/n follows sigma(n)/n (correlation above 0.9 on the window; "
                      "the lower-order terms in ln n depend on n through other divisor sums)")


# ---------------------------------------------------------------- S5
def s5():
    print("S5 two-port kernels u^T G v with uniform states")
    for name, G in (("Gaussian kernel [[5,2],[2,1]] (det 1)", ((5, 2), (2, 1))),
                    ("first hyperbolic kernel [[11,6],[6,3]] (det -3)", ((11, 6), (6, 3))),
                    ("Session260 diag(1, 2) (det 2)", ((1, 0), (0, 2)))):
        for q in (9, 25, 49):
            cnt = [0] * q
            for u1 in range(q):
                for u2 in range(q):
                    w1 = (G[0][0] * u1 + G[1][0] * u2) % q
                    w2 = (G[0][1] * u1 + G[1][1] * u2) % q
                    g = math.gcd(math.gcd(w1, w2), q)
                    for c in range(0, q, g):
                        cnt[c] += q * g
            same = all(F(cnt[c] * q, q ** 4) == L(q, c) for c in range(q))
            det = G[0][0] * G[1][1] - G[0][1] * G[1][0]
            expect = det % primes_of(q)[0] != 0
            check(same == expect, f"{name}, modulo {q}: rank-2 law = {same} (expected {expect})")


# ---------------------------------------------------------------- S6
def s6():
    print("S6 central extensions of F_2^2 by F_2")
    V = range(4)
    nz = [1, 2, 3]
    cocycles = []
    for bits in range(1 << 9):
        f = [[0] * 4 for _ in V]
        k = 0
        for a in nz:
            for b in nz:
                f[a][b] = (bits >> k) & 1
                k += 1
        if all((f[a][b] + f[a ^ b][c] + f[b][c] + f[a][b ^ c]) % 2 == 0 for a in V for b in V for c in V):
            cocycles.append(f)
    check(len(cocycles) == 16, f"normalized 2-cocycles: {len(cocycles)}")
    cob = [[(0 if (a == 0 or b == 0 or a == b) else 1) for b in V] for a in V]

    def key(f):
        return tuple(f[a][b] for a in nz for b in nz)
    classes = {}
    for f in cocycles:
        g2 = [[(f[a][b] + cob[a][b]) % 2 for b in V] for a in V]
        rep = min(key(f), key(g2))
        classes.setdefault(rep, []).append(f)
    check(len(classes) == 8 and all(len(v) == 2 for v in classes.values()), f"classes modulo the coboundary: {len(classes)}, two cocycles each")
    summary = {}
    for rep, fs in sorted(classes.items()):
        f = fs[0]
        # group on V x F_2
        elems = [(a, z) for a in V for z in (0, 1)]

        def mul(x, y):
            return (x[0] ^ y[0], (x[1] + y[1] + f[x[0]][y[0]]) % 2)

        def order(x):
            y, k = x, 1
            while y != (0, 0):
                y = mul(y, x)
                k += 1
            return k
        census = tuple(sorted(order(x) for x in elems))
        abel = all(mul(x, y) == mul(y, x) for x in elems for y in elems)
        squares = tuple(f[a][a] for a in nz)            # 1 = the unit has square -1
        w = sum(squares)
        name = {(True, 0): "C2^3", (True, 2): "C4 x C2", (False, 1): "D4", (False, 3): "Q8"}[(abel, w)]
        n4 = sum(1 for o in census if o == 4)
        check({"C2^3": 0, "C4 x C2": 4, "D4": 2, "Q8": 6}[name] == n4, f"squares -1 at {[nz[i] for i in range(3) if squares[i]]}: {name}, {n4} elements of order 4, commutative = {abel}")
        # real algebra e_a e_b = (-1)^f(a,b) e_(a+b): zero divisors with coordinates in {-1, 0, 1}
        zd = False
        vecs = [v for v in itertools.product((-1, 0, 1), repeat=4) if any(v)]
        for x in vecs:
            for y in vecs:
                z = [0] * 4
                for a in V:
                    if x[a]:
                        for b in V:
                            if y[b]:
                                z[a ^ b] += x[a] * y[b] * (-1) ** f[a][b]
                if not any(z):
                    zd = True
                    break
            if zd:
                break
        check(zd == (name != "Q8"), f"    zero divisors: {zd}")
        # sign tables of the two cocycles of the class on the three points
        kinds = []
        for ff in fs:
            anti = all(ff[a][b] != ff[b][a] for a in nz for b in nz if a < b)
            symm = all(ff[a][b] == ff[b][a] for a in nz for b in nz if a < b)
            if anti:
                out = {a: sum(1 for b in nz if b != a and ff[a][b] == 0) for a in nz}
                if sorted(out.values()) == [1, 1, 1]:
                    kinds.append("cyclic tournament")
                else:
                    mid = [a for a in nz if out[a] == 1][0]
                    kinds.append(f"transitive tournament with middle vertex {mid}")
                    check(squares[nz.index(mid)] == 1 and w == 1, "    the middle vertex is the unit of square -1")
            elif symm:
                kinds.append("signed triangle with " + str(sum(ff[a][b] for a in nz for b in nz if a < b)) + " negative edges")
            else:
                kinds.append("mixed")
        print(f"    sign tables of the class: {kinds}")
        summary[name] = summary.get(name, 0) + 1
        check((w % 2 == 1) == (not abel), "    commutative iff the number of units of square -1 is even")
    check(summary == {"C2^3": 1, "C4 x C2": 3, "D4": 3, "Q8": 1}, f"classes by group: {summary} (1 + 3 + 3 + 1 by the number 0, 2, 1, 3 of units of square -1)")


# ---------------------------------------------------------------- S7
def s7():
    print("S7 signed-integer dictionary")

    def lab(n):
        return abs(4 * n + 1)

    def unlab(m):
        e = 1 if m % 4 == 1 else -1
        return (e * m - 1) // 4
    check(all(unlab(lab(n)) == n for n in range(-10 ** 5, 10 ** 5)) and sorted(lab(n) for n in range(-500, 500)) == list(range(1, 2000, 2)),
          "n -> |4n + 1| is a bijection from Z onto the odd positive integers (inverse m -> (chi_4(m) m - 1)/4)")

    def col(n):
        return 3 * n + 1 if n % 2 else n // 2
    good = True
    for n in range(-3000, 3000):
        ms = 4 * n + 1
        img = 4 * col(n) + 1
        good &= img == (3 * ms + 2 if n % 2 else (ms + 1) // 2)
        good &= lab(-3 * n - 1) == 3 * lab(n)
        for k in range(-30, 30):
            good &= lab(4 * n * k + n + k) == lab(n) * lab(k)
    check(good, "signed label m* = 4n + 1: odd step m* -> 3 m* + 2, even step m* -> (m* + 1)/2; tripling the label is n -> -3n - 1; "
                "products of labels are n * k = 4nk + n + k")
    cycles = []
    for start in (0, 1, -1, -5, -17):
        c, n = [start], col(start)
        while n != start:
            c.append(n)
            n = col(n)
        cycles.append(c)
    labels = [[lab(n) for n in c] for c in cycles]
    for c, l in zip(cycles, labels):
        print(f"    cycle {c}: labels {l}")
    alll = [x for l in labels for x in l]
    check(len(alll) == 29 and len(set(alll)) == 29 and [x for x in alll if x in (7, 21)] == [7],
          "the 29 arcs of the five known cycles of 3n + 1 on Z have 29 distinct labels; exactly one of them (7, the arc -2 -> -1) is one of the two missing path counts")
    check(lab(5) == 21 and col(5) == 16 and lab(-2) == 7 and col(-2) == -1, "21 is the label of the arc 5 -> 16 and 7 of the arc -2 -> -1")
    tr = [(4 ** k - 1) // 3 for k in range(1, 7)]
    check(all(lab(tr[i]) == tr[i + 1] and col(tr[i]) == 4 ** (i + 1) for i in range(5)),
          f"trunk: the arcs (4^k - 1)/3 -> 4^k have labels (4^(k+1) - 1)/3: {tr}")
    # spectra of the orders 3..6
    spectra = {}
    for n in (3, 4, 5, 6):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        hs = set()
        full = (1 << n) - 1
        for bits in range(1 << len(pairs)):
            out = [0] * n
            for k, (i, j) in enumerate(pairs):
                if (bits >> k) & 1:
                    out[i] |= 1 << j
                else:
                    out[j] |= 1 << i
            dp = [[0] * n for _ in range(1 << n)]
            for v in range(n):
                dp[1 << v][v] = 1
            for mask in range(1, 1 << n):
                row = dp[mask]
                for v in range(n):
                    c = row[v]
                    if c:
                        nxt = out[v] & ~mask
                        while nxt:
                            w = (nxt & -nxt).bit_length() - 1
                            dp[mask | (1 << w)][w] += c
                            nxt &= nxt - 1
            hs.add(sum(dp[full]))
        spectra[n] = sorted(hs)
        ints = sorted(unlab(h) for h in hs)
        print(f"    order {n}: H in {spectra[n]}")
        print(f"             integers n = (chi_4(H) H - 1)/4: {ints}")
    check(spectra[4] == [1, 3, 5] and sorted(unlab(h) for h in spectra[4]) == [-1, 0, 1],
          "order 4: the path counts 1, 3, 5 are the arcs leaving 0, -1, 1 (0 -> 0, -1 -> -2, 1 -> 4), one in each of the three cycles through 0, -1, 1")
    lo, hi = -(max(spectra[6]) + 1) // 4, (max(spectra[6]) - 1) // 4
    miss5 = sorted(set(range(-4, 4)) - set(unlab(h) for h in spectra[5]))
    miss6 = sorted(set(range(lo, hi + 1)) - set(unlab(h) for h in spectra[6]))
    print(f"    integers of [-4, 3] missing at order 5: {miss5}; integers of [{lo}, {hi}] missing at order 6: {miss6}")
    check(miss5 == [-2], "order 5 covers the integers of [-4, 3] except -2 (the count 7)")


# ---------------------------------------------------------------- S8
def s8():
    print("S8 primes with balanced divisible states")
    size = 1 << 16
    psi = [0] * size
    for h in range(1, size):
        hh = h >> 1
        psi[h] = psi[hh] + (hh if h % 2 == 0 else (hh // 2 + 1))
    xs = [2 * psi[h] + 1 for h in range(size)]
    check(all((xs[h] == 2 * h + 1) == (h == 0 or h & (h - 1) == 0) for h in range(size)), "x = y exactly at h = 0 and at the powers of 2 (fans, x = y = 2h + 1)")
    rows, agree = [], True
    first_is_fan = []
    for p in range(3, 200, 2):
        if any(p % d == 0 for d in range(3, int(p ** 0.5) + 1, 2)):
            continue
        o, t = 1, 2 % p
        while t != 1:
            t = t * 2 % p
            o += 1
        m = next((k for k in range(1, 2 * p) if (pow(2, k, p) + 1) % p == 0), None)
        agree &= (m is not None) == (o % 2 == 0) and (m is None or m == o // 2)
        fan = 1 << (m - 1) if m else None
        first = next((h for h in range(1, size) if xs[h] % p == 0 and (2 * h + 1) % p == 0), None)
        if m and fan < size:
            agree &= xs[fan] % p == 0 and (2 * fan + 1) % p == 0
            if first == fan:
                first_is_fan.append(p)
        rows.append((p, o, m, fan, first))
    check(agree, "p | 2^t + 1 for some t iff ord_p(2) is even; the least t is ord_p(2)/2 and the first balanced p-divisible state is the fan h = 2^(t-1)")
    print("    (p, ord_p(2), t, first fan, first p-divisible state) for p < 60: " + "; ".join(str(r) for r in rows if r[0] < 60))
    print(f"    primes < 200 whose first p-divisible state is the fan: {first_is_fan}")
    return xs


# ---------------------------------------------------------------- S9
def s9(xs):
    print("S9 Brahmagupta identity")
    good = True
    for hA in range(0, 300, 7):
        for hB in range(1, 300, 11):
            xA, yA, xB, yB = xs[hA], 2 * hA + 1, xs[hB], 2 * hB + 1
            for s in range(0, 5):
                MA, MB, MAB = xA * xA + 2 ** s * yA * yA, xB * xB + 2 ** s * yB * yB, xA * xB + 2 ** s * yA * yB
                good &= MA * MB == MAB ** 2 + 2 ** s * (xA * yB - xB * yA) ** 2
    check(good, "M_s(A, A) M_s(B, B) = M_s(A, B)^2 + 2^s (x_A y_B - x_B y_A)^2 on one-apex states")


# ---------------------------------------------------------------- S10
@njit(cache=True)
def base_head_census(psi, beta, X, primes, cnt, exact):
    H = psi.shape[0]
    tot = 0
    for h1 in range(H):
        y1 = 2 * h1 + 1
        X1 = 1 + beta * psi[h1]
        if beta * (y1 * y1 - 1) // 2 + X1 > X:
            if beta * (y1 * y1 - 1) // 2 > X:
                break
            continue
        for h2 in range(h1, H):
            y2 = 2 * h2 + 1
            if beta * (y1 * y2 - 1) // 2 + X1 > X:
                break
            N = X1 * (1 + beta * psi[h2]) + beta * (y1 * y2 - 1) // 2
            if N > X:
                continue
            w = 1 if h2 == h1 else 2
            tot += w
            for j in range(primes.shape[0]):
                p = primes[j]
                t = 1 if (X1 % p == 0 and y1 % p == 0) else 0
                cnt[j, t, N % p] += w
                if t == 1 and (2 * N + beta) % p != 0:
                    exact[j] += 1
    return tot


def s10():
    print("S10 F1 in base beta, sheet s = 0 (N_beta <= 2^22)")
    X = 1 << 22
    size = 1 << 12
    psi = np.zeros(size, np.int64)
    for h in range(1, size):
        hh = h >> 1
        psi[h] = psi[hh] + (hh if h % 2 == 0 else (hh // 2 + 1))
    for beta in (2, 3, 5, 10):
        primes = np.array([p for p in (3, 5, 7, 11, 13) if beta % p], np.int64)
        cnt = np.zeros((len(primes), 2, 16), np.int64)
        exact = np.zeros(len(primes), np.int64)
        tot = base_head_census(psi, beta, X, primes, cnt, exact)
        cells = []
        good = True
        for j, p in enumerate(primes):
            p = int(p)
            c0 = (-beta * pow(2, -1, p)) % p
            D1, Dp = int(cnt[j, 0].sum()), int(cnt[j, 1].sum())
            rem = max(abs(int(cnt[j, 0, c] + cnt[j, 1, c]) - (D1 / p + (Dp if c == c0 else 0))) for c in range(p)) * p / tot
            fac0 = (cnt[j, 0, c0] + cnt[j, 1, c0]) * p / tot
            good &= exact[j] == 0 and rem < 0.25
            cells.append(f"p={p}: class {c0} factor {fac0:.3f}, F p^2 = {Dp / tot * p * p:.2f}, rem {rem:.3f}")
        check(good, f"beta = {beta:2d}: {tot} triples; p-divisible smaller states give 2N + beta = 0 (mod p) exactly; " + "; ".join(cells))


# ---------------------------------------------------------------- S11
def s11(xs):
    print("S11 Tao's digit condition for the joins")
    vec = {}
    rows = []
    for h1 in range(0, 64):
        for h2 in range(h1, 64):
            for s in range(0, 5):
                p1, p2 = (xs[h1] - 1) // 2, (xs[h2] - 1) // 2
                c3 = (2 ** s * (2 * h1 + 1) * (2 * h2 + 1) - (1 if s == 0 else 0)) // 2
                a1, a2 = p1 + p2 + c3, p1 * p2
                rows.append((a1, a2))
                vec.setdefault(1 + 2 * a1 + 4 * a2, set()).add((a1, a2))
    good = all(divmod(divmod(1 + a1 * b + a2 * b * b, b)[0], b) == (a2, a1) and (1 + a1 * b + a2 * b * b) % b == 1
               for (a1, a2) in rows for b in (max(2, max(a1, a2) + 1), max(a1, a2) + 7))
    check(good, "for every base b >= 2 with b > max(alpha_1, alpha_2) the base-b numeral of N_b is (alpha_2)(alpha_1)(1): the count displays the packing numbers")
    coll = sorted((H, sorted(v)) for H, v in vec.items() if len(v) > 1)
    print(f"    base 2 (carries): least path counts with two packing vectors (alpha_1, alpha_2) in the family: {coll[:4]}")
    check(len(coll) > 0, "in base 2 different packing vectors collide (the count alone does not determine the state)")
    for b in (2, 3, 10, 100):
        share = sum(1 for (a1, a2) in rows if max(a1, a2) < b) / len(rows)
        print(f"    base {b:3d}: share of the {len(rows)} joins (h1 <= h2 < 64, s <= 4) with digits at most b - 1: {share:.4f}")


def main():
    s1()
    s2()
    s3()
    s4()
    s5()
    s6()
    s7()
    xs = s8()
    s9(xs)
    s10()
    s11(xs)
    print("ALL CHECKS PASSED" if OK else "CHECK FAILED")
    sys.exit(0 if OK else 1)


if __name__ == "__main__":
    main()
