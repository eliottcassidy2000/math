#!/usr/bin/env python3
"""Verifier of the corrections to the S91 README: small exact checks of statements quoted in the corrected text.

Own code; nothing is imported from the session or from the auditors.  Every check prints PASS / FAIL.
"""
import math
from fractions import Fraction
from itertools import product

RES = []


def check(name, cond, detail=""):
    RES.append(bool(cond))
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + str(detail) + "]") if detail != "" else ""), flush=True)


# ---------------------------------------------------------------- one-apex states
def psi_table(n):
    ps = [0] * n
    for h in range(1, n):
        k = h >> 1
        ps[h] = ps[k] + k if h % 2 == 0 else ps[k] + k // 2 + 1
    return ps


NMAX = 1 << 19
PS = psi_table(NMAX)


def xy(h):
    return 2 * PS[h] + 1, 2 * h + 1


# (1) x = y exactly at h = 0 and h = 2^j
eq = [h for h in range(NMAX) if PS[h] == h]
check("1  x(h) = y(h) exactly at h = 0 and the powers of two (h < 2^19), with x = y = 2^(j+1) + 1",
      eq == [0] + [1 << j for j in range(19)] and all(xy(1 << j)[0] == (1 << (j + 1)) + 1 for j in range(19)))


def first_div(p, count=1):
    out = []
    for h in range(NMAX):
        x, y = xy(h)
        if x % p == 0 and y % p == 0:
            out.append(h)
            if len(out) == count:
                break
    return out


want = {3: 1, 5: 2, 7: 38, 11: 16, 13: 32, 17: 8, 19: 218, 29: 1435, 31: 8664, 37: 1350, 41: 266, 43: 64, 127: 109283,
        257: 128, 65537: 32768, 1039: 405729}
got = {p: first_div(p)[0] for p in want}
check("2  first p-divisible states quoted in sections 1 and 3", got == want, got)
check("2b 13: divisible states 32, then 721", first_div(13, 2) == [32, 721], first_div(13, 2))
check("2c 61: first divisible state 5886 > 61^2 = 3721; 73: first divisible state printed",
      first_div(61) == [5886] and len(first_div(73)) == 1, (first_div(61), first_div(73)))
bal = {3: 1, 5: 2, 13: 32, 19: 256, 27: 256, 41: 512, 11: 16, 17: 8, 43: 64, 257: 128}
check("2d balanced states h = 2^(m-1) with p | 2^m + 1: 3: 1, 5: 2, 11: 16, 13: 32, 17: 8, 19 and 27: 256, 41: 512, 43: 64, 257: 128",
      all(xy(h)[0] == xy(h)[1] and xy(h)[0] % p == 0 for p, h in bal.items()))
check("2e first balanced divisible states of 3, 13, 19, 27 are 1, 32, 256, 256; 5 has h = 2 before 32; 19 has the "
      "unbalanced h = 218 (x, y) = (361, 437)",
      xy(218) == (361, 437) and xy(2) == (5, 5) and xy(32) == (65, 65) and xy(256) == (513, 513)
      and min(1 << j for j in range(19) if xy(1 << j)[0] % 27 == 0) == 256
      and min(1 << j for j in range(19) if xy(1 << j)[0] % 19 == 0) == 256
      and min(1 << j for j in range(19) if xy(1 << j)[0] % 13 == 0) == 32)

# (3) words with x(h) = n
words = {n: [h for h in range(2 * n) if xy(h)[0] == n] for n in (3, 5, 7, 13, 19, 27, 1039)}
check("3  one-apex words with x = 3, 5, 7, 13, 19, 27, 1039: [1], [2, 3], none, none, [10, 13], [17, 23], none",
      words == {3: [1], 5: [2, 3], 7: [], 13: [], 19: [10, 13], 27: [17, 23], 1039: []}, words)

# (4) orders of 2
def ordm(m):
    k, v = 1, 2 % m
    while v != 1:
        v = v * 2 % m
        k += 1
    return k


check("4  ord_p(2) = 2, 4, 3, 12, 18, 18, 519 for p = 3, 5, 7, 13, 19, 27, 1039",
      [ordm(p) for p in (3, 5, 7, 13, 19, 27, 1039)] == [2, 4, 3, 12, 18, 18, 519])
check("4b ord_m(2) = 3, 10, 5, 14, 4, 18, 18, 7, 8, 2 for m = 7, 93, 31, 381, 15, 219, 171, 127, 51, 3, and 9 at 73",
      [ordm(m) for m in (7, 93, 31, 381, 15, 219, 171, 127, 51, 3)] == [3, 10, 5, 14, 4, 18, 18, 7, 8, 2] and ordm(73) == 9)
check("4c 29, 37, 61 divide 2^14 + 1, 2^18 + 1, 2^30 + 1; 5 * 13 = 2^6 + 1; 27 * 19 = 2^9 + 1; 41 | 2^10 + 1",
      (2 ** 14 + 1) % 29 == 0 and (2 ** 18 + 1) % 37 == 0 and (2 ** 30 + 1) % 61 == 0 and 5 * 13 == 2 ** 6 + 1
      and 27 * 19 == 2 ** 9 + 1 and (2 ** 10 + 1) % 41 == 0)

# (5) p-adic factor for uniform states modulo p^k (exact enumeration)
def padic(p, k, u):
    m = p ** k
    cnt = [0] * m
    for a1, b1 in product(range(m), repeat=2):
        row = [0] * m
        for a2 in range(m):
            base = a1 * a2 % m
            for b2 in range(m):
                row[(base + u * b1 * b2) % m] += 1
        for c in range(m):
            cnt[c] += row[c]
    tot = m ** 4
    return [Fraction(cnt[c] * m, tot) for c in range(m)]


def val(c, p, k):
    if c == 0:
        return k
    v = 0
    while c % p == 0:
        c //= p
        v += 1
    return v


okf = True
detail = {}
for p, k, u in ((3, 2, 1), (3, 2, 2), (5, 2, 1), (5, 2, 3), (3, 3, 1), (7, 2, 3)):
    f = padic(p, k, u)
    m = p ** k
    for c in range(m):
        v = val(c, p, k)
        if v < k:
            okf &= f[c] == (1 + Fraction(1, p)) * (1 - Fraction(1, p ** (v + 1)))
            okf &= f[c] == (1 - Fraction(1, p * p)) * Fraction(sum(p ** i for i in range(v + 1)), p ** v)
    cls = sum(f[c] for c in range(0, m, p)) / (m // p)
    okf &= cls == 1 - Fraction(1, p * p) + Fraction(1, p)
    okf &= sum(f) == m
    detail[(m, u)] = [str(f[1]), str(f[p]), str(f[0])]
check("5  uniform states modulo p^k (9, 25, 27, 49): factor (1 - p^-2) sigma(p^v)/p^v = (1 + 1/p)(1 - p^(-v-1)) at "
      "v < k; mean over the class 0 mod p is 1 - 1/p^2 + 1/p; mean over all classes is 1", okf, detail)
check("5b modulo 9 the three values are 8/9, 32/27, 35/27 = 0.889, 1.185, 1.296",
      detail[(9, 1)] == ["8/9", "32/27", "35/27"])

# (6) 8/pi^2 and the generic limits of the table
e8 = 8 / math.pi ** 2
prod_ = 1.0
for q in range(3, 200000, 2):
    if all(q % d for d in range(3, int(q ** 0.5) + 1, 2)):
        prod_ *= 1 - 1 / q ** 2
check("6  prod over odd p of (1 - 1/p^2) = 8/pi^2 = 0.8106", abs(prod_ - e8) < 1e-5 and round(e8, 4) == 0.8106, (prod_, e8))
gen = {3: "0.375", 5: "0.208", 7: "0.146", 11: "0.092", 13: "0.077", 17: "0.059", 19: "0.053", 29: "0.035", 31: "0.032",
       43: "0.023", 257: "0.004", 41: "0.024"}
check("6b generic limits p/(p^2 - 1) of the section 1 table (29: 0.035) and 41: 0.024",
      all("%.3f" % (p / (p * p - 1)) == s for p, s in gen.items()))
check("6c ratio (1 - 1/9 + 1/3)/(1 - 1/9) = 1.375", Fraction(11, 9) / Fraction(8, 9) == Fraction(11, 8))


# (7) squarefree-level products at 1039 and the sheet weights
def li2(z):
    return sum(z ** k / k ** 2 for k in range(1, 400))


k0 = math.pi ** 2 / 12 + 2 * li2(-0.5) - li2(-0.25)
k1 = li2(-0.25) - li2(-0.5)
w0t, w1t = k0 / (k0 + k1), k1 / (k0 + k1)
w0p = k0 / (math.pi ** 2 / 12)
check("7  kappa_0 = 0.16154, kappa_1 = 0.21251; triple shares (0.43186, 0.56814); word-pair shares (0.19641, 0.80359); "
      "g_bar = 2 (kappa_0 + kappa_1) = 0.7481",
      abs(k0 - 0.1615389173) < 1e-9 and abs(k1 - 0.2125139092) < 1e-9 and abs(w0t - 0.43186124) < 1e-7
      and abs(w0p - 0.19640777) < 1e-7 and abs(2 * (k0 + k1) - 0.7481056530) < 1e-9, (k0, k1, w0t, w0p))


def ssf(n):
    """(8/pi^2) prod over odd p | n of (1 + p/(p^2 - 1))."""
    s = e8
    d = 3
    m = n
    while m % 2 == 0:
        m //= 2
    while d * d <= m:
        if m % d == 0:
            s *= 1 + d / (d * d - 1)
            while m % d == 0:
                m //= d
        d += 2
    if m > 1:
        s *= 1 + m / (m * m - 1)
    return s


S1, S0 = ssf(1039), ssf(1040)
check("7b N = 1039: S_1 = 0.811, S_0 = 1.055; with the TRIPLE shares w_0 S_0 + w_1 S_1 = 0.917 "
      "(with the word-pair shares of section 2 it would be 0.859)",
      "%.3f %.3f %.3f %.3f" % (S1, S0, w0t * S0 + w1t * S1, w0p * S0 + (1 - w0p) * S1) == "0.811 1.055 0.917 0.859",
      (S1, S0, w0t * S0 + w1t * S1, w0p * S0 + (1 - w0p) * S1))
lim = 1.0
for p in (3, 5, 7, 11, 13, 17):
    lim *= 1 - 1 / p ** 2
lim0 = 1.0
for p in (3, 5, 7, 11, 13, 17):
    # N coprime to p: on sheet 0 the class of N + 1 is distinguished with probability 1/(p - 1)
    lim0 *= (1 - 1 / p ** 2) + (1 / p) / (p - 1)
check("7c N coprime to 255255: sheet-1 factor 0.8213, sheet-0 factor 1.0714, word pairs 0.870",
      "%.4f %.4f %.3f" % (lim, lim0, (1 - w0p) * lim + w0p * lim0) == "0.8213 1.0714 0.870",
      (lim, lim0, (1 - w0p) * lim + w0p * lim0))

# (8) the prime 2: signed transfer matrix of q mod 2
import numpy as np

M = np.array([[1, 0, 1, 0], [1, 0, -1, 0], [0, 1, 0, -1], [0, 1, 0, 1]], dtype=object)
M12 = np.linalg.matrix_power(M, 12)
check("8  M M^T = 2I and M^12 = -64 I for the signed transfer matrix of q mod 2",
      (M.dot(M.T) == 2 * np.eye(4, dtype=object)).all() and (M12 == -64 * np.eye(4, dtype=object)).all())
QQ = [h - PS[h] for h in range(1 << 18)]
S = [sum((-1) ** (QQ[h] & 1) for h in range(1 << n)) for n in range(19)]
D = [sum(1 if xy(h)[0] % 4 == 1 else -1 for h in range(1 << n)) for n in range(19)]
check("8b S_(n+12) = -64 S_n and D_(n+12) = -64 D_n for n <= 6 (direct sums over h < 2^18); S_0..S_11 and D_0..D_11 as in THM-5128 (E)",
      all(S[n + 12] == -64 * S[n] and D[n + 12] == -64 * D[n] for n in range(7))
      and S[:12] == [1, 2, 2, 0, 0, 0, -8, -16, -16, -32, -64, -64] and D[:12] == [1, 0, 2, 4, 0, 0, 8, 0, -16, 0, 0, -64],
      (S[:12], D[:12]))


def corank(n):
    rows = []
    for i in range(n):
        r = 0
        for j in range(n):
            if abs(i - j) in (1, 2):
                r |= 1 << j
        rows.append(r)
    rank = 0
    for bit in range(n):
        piv = None
        for idx in range(rank, n):
            if (rows[idx] >> bit) & 1:
                piv = idx
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for idx in range(n):
            if idx != rank and (rows[idx] >> bit) & 1:
                rows[idx] ^= rows[rank]
        rank += 1
    return n - rank


cr = [corank(n) for n in range(1, 37)]
check("8c corank of the band form: 1, 0, 1, 2, 1, 0 with period 6 (n = 1..36); the rank is not periodic",
      cr == [1, 0, 1, 2, 1, 0] * 6, cr[:12])


def two_adic_unitary(j):
    """T_k T_k^* = 2I for odd k, fails for even k != 0, on Z/2^(j+1)."""
    n = 1 << (j + 1)
    good_odd, bad_even = True, True
    for k in range(1, 1 << j):
        T = np.zeros((n, n), dtype=complex)
        for w in range(n):
            c = (w + (w & 1)) // 2
            T[w, (2 * w) % n] += 1
            T[w, (2 * w + 1) % n] += np.exp(2j * np.pi * k * c / (1 << j))
        ok = np.allclose(T.dot(T.conj().T), 2 * np.eye(n))
        if k % 2:
            good_odd &= ok
        else:
            bad_even &= not ok
    return good_odd, bad_even


check("8d modulo 2^j, j = 2..5: T_k T_k^* = 2I for every odd k and for no even k != 0",
      all(two_adic_unitary(j) == (True, True) for j in range(2, 6)))
bnd = True
QA = np.array(QQ[:1 << 16], dtype=np.int64)
for j in range(1, 6):
    for k in range(1, 1 << j, 2):
        cs = np.cumsum(np.exp(2j * np.pi * ((k * QA) % (1 << j)) / (1 << j)))
        for n in range(1, 17):
            bnd &= abs(cs[(1 << n) - 1]) <= 2 ** ((n + j + 1) / 2) + 1e-6
check("8e |sum_{h < 2^n} e(k q(h)/2^j)| <= 2^((n+j+1)/2) for odd k, j <= 5, n <= 16", bnd)

# (9) self-joins mod 8
r8 = [sorted({(x * x + (y * y << s) - (1 if s == 0 else 0)) % 8 for h in range(4096) for (x, y) in [xy(h)]}) for s in range(6)]
check("9  self-join values x^2 + 2^s y^2 - [s = 0] are 1, 3, 5, 1, 1, 1 (mod 8) for s = 0..5; 7 = 1039 = 7 (mod 8); 27 = 3^2 + 2 * 3^2",
      r8 == [[1], [3], [5], [1], [1], [1]] and 1039 % 8 == 7 and 27 == 9 + 18, r8)

# (10) Busch-Moon minimum
def f(n):
    return {0: 3 * 5 ** ((n - 3) // 3), 1: 5 ** ((n - 1) // 3), 2: 9 * 5 ** ((n - 5) // 3)}[n % 3]


check("10 f(3..12) = 3, 5, 9, 15, 25, 45, 75, 125, 225, 375; f(13) = 625, f(14) = 1125, f(15) = 1875, f(16) = 3125; "
      "largest order with f(n) <= 513, 1025, 2049 is 12, 13, 15",
      [f(n) for n in range(3, 13)] == [3, 5, 9, 15, 25, 45, 75, 125, 225, 375] and (f(13), f(14), f(15), f(16)) == (625, 1125, 1875, 3125)
      and [max(n for n in range(3, 30) if f(n) <= D_) for D_ in (513, 1025, 2049)] == [12, 13, 15])

# (11) strip counts
a = {5: 188, 6: 556, 7: 1660}
for n in range(8, 14):
    a[n] = 2 * a[n - 1] + a[n - 2] + 6 * a[n - 3]
V, U = [2], [0]          # ((-1 + sqrt(-7))/2)^n = (V_n + U_n sqrt(-7))/2
for n in range(1, 60):
    v, u = V[-1], U[-1]
    V.append((-v - 7 * u) // 2)
    U.append((v - u) // 2)
check("11 strip counts: recurrence gives a(8) = 5004, a(10) = 44972; 42 a(n) = 32 * 3^n - 9 V_n - 21 U_n for n = 5..13; "
      "V_n^2 + 7 U_n^2 = 2^(n+2); |U_n| = 1 exactly at n = 1, 2, 3, 5, 13 for n < 60",
      a[8] == 5004 and a[10] == 44972 and all(42 * a[n] == 32 * 3 ** n - 9 * V[n] - 21 * U[n] for n in range(5, 14))
      and all(V[n] ** 2 + 7 * U[n] ** 2 == 2 ** (n + 2) for n in range(60))
      and [n for n in range(60) if abs(U[n]) == 1] == [1, 2, 3, 5, 13], [a[n] for n in range(5, 11)])
check("11b x^2 + x + 2 has discriminant -7 and roots of norm 2; 3x^2 + 9xy + 7y^2 = X^2 + XY + Y^2 at (x + y, x + 2y)",
      1 - 8 == -7 and all(3 * x * x + 9 * x * y + 7 * y * y == (x + y) ** 2 + (x + y) * (x + 2 * y) + (x + 2 * y) ** 2
                           for x in range(-5, 6) for y in range(-5, 6)))

# (12) Euler-Sombor arithmetic
check("12 norms a^2 + ab + b^2 at (1,2), (1,3), (2,3), (3,3), (2,2), (1,1) = 7, 13, 19, 27, 12, 3; 76 - 75 = 1; "
      "0 < eps < delta; 2 + 5 + 13 + 33 + 89 + 240 + 657 = 1039; 9 * 2^(s-3) >= 2^s + 1 > 9 * 2^(s-4) for s = 3..12",
      [i * i + i * j + j * j for i, j in ((1, 2), (1, 3), (2, 3), (3, 3), (2, 2), (1, 1))] == [7, 13, 19, 27, 12, 3]
      and 0 < 2 * 19 ** 0.5 - 5 * 3 ** 0.5 < 13 ** 0.5 + 2 * 3 ** 0.5 - 7 ** 0.5 - 19 ** 0.5
      and 2 + 5 + 13 + 33 + 89 + 240 + 657 == 1039
      and all(9 * 2 ** (s - 3) >= 2 ** s + 1 for s in range(3, 13))
      and all(9 * 2 ** (s - 4) < 2 ** s + 1 for s in range(4, 13)))

# (13) registered tolerance against the measured dispersion
nodd = 58720256
se = lambda v, p: (v * p / nodd) ** 0.5
check("13 se at 65537: 0.022 / 0.021 measured (0.446 / 0.390), 0.029 registered (0.75); at 43691: 0.018 / 0.017; "
      "the registered 4 se band is 5.2 / 5.5 measured standard errors",
      ["%.3f" % se(v, 65537) for v in (0.446, 0.390, 0.75)] == ["0.022", "0.021", "0.029"]
      and ["%.3f" % se(v, 43691) for v in (0.446, 0.390)] == ["0.018", "0.017"]
      and ["%.1f" % (4 * (0.75 / v) ** 0.5) for v in (0.446, 0.390)] == ["5.2", "5.5"])

print("ALL CHECKS PASSED (%d checks)" % len(RES) if all(RES) else "SOME CHECKS FAILED")
raise SystemExit(0 if all(RES) else 1)
