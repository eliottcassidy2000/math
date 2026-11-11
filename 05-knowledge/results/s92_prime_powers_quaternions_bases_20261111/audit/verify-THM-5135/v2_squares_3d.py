"""Verifier of THM-5135 (S92 audit follow-up): statement (2), (3)(a), (3)(c), (3)(d) and the step of (3)(d)
(zero counts by rank), the densities, the classical counts.  Own code, exact.
Run: python v2_squares_3d.py > v2_squares_3d.out
"""
import itertools
import math
from fractions import Fraction as F

import numpy as np

FAILS = []


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def table(q, scale=1):
    d = [0] * q
    for a in range(q):
        for b in range(q):
            d[(scale * a * b) % q] += 1
    return d


def squares(q):
    d = [0] * q
    for a in range(q):
        d[(a * a) % q] += 1
    return d


def conv(a, b):
    q = len(a)
    out = [0] * q
    for i in range(q):
        if a[i]:
            for j in range(q):
                if b[j]:
                    out[(i + j) % q] += a[i] * b[j]
    return out


def conv_power(d, r):
    out = d
    for _ in range(r - 1):
        out = conv(out, d)
    return out


# ---------------------------------------------------------------- [1] statement (2)
print("[1] (2): 2r squares against the hyperbolic form of rank r at odd p")
ok = True
okz = True
for q in (3, 5, 7, 11, 13, 9, 25, 27, 49):
    p = min(d for d in range(2, q + 1) if q % d == 0)
    for r in (1, 2, 3, 4):
        if q > 13 and r > 3 and q != 25:
            continue
        a = conv_power(squares(q), 2 * r)
        b = conv_power(table(q), r)
        expect = (r % 2 == 0) or (p % 4 == 1)
        ok &= ((a == b) == expect)
        if q == p:
            eta_sq = 1 if expect else -1
            okz &= b[0] == p ** (2 * r - 1) + (p ** r - p ** (r - 1))
            okz &= a[0] == p ** (2 * r - 1) + eta_sq * (p ** r - p ** (r - 1))
check(ok, "equal class counts iff r even or p = 1 (mod 4): moduli 3, 5, 7, 11, 13, 9, 25 (r <= 4), 27, 49 (r <= 3)")
check(okz, "zeros modulo p: p^(2r-1) + eta (p^r - p^(r-1)), eta = +1 for the hyperbolic form, -1 for the other class (p <= 13, r <= 4)")
check(all(F(-1, 4) ** r == F((-1) ** r, 4 ** r) for r in range(1, 6)), "discriminant of the hyperbolic form of rank r: (-1/4)^r, in the class of (-1)^r")

# ---------------------------------------------------------------- [2] the prime 2
print("[2] (3)(a), (3)(c): the prime 2")
for k in range(1, 7):
    q = 2 ** k
    four = conv_power(squares(q), 4)
    seam2 = conv_power(table(q), 2)
    eight = conv(four, four)
    seam4 = conv(seam2, seam2)
    d4 = sum(1 for c in range(q) if four[c] != seam2[c])
    d8 = sum(1 for c in range(q) if eight[c] != seam4[c])
    dens4 = {F(four[c] * q, q ** 4) for c in range(1, q, 2)}
    check(d4 == q and d8 == q and dens4 == {F(1)},
          f"modulo {q:2d}: four squares differ from the rank-2 seam in {d4} of {q} classes, eight squares from the rank-4 seam in {d8}; "
          f"four-square factor at the odd classes: {sorted(dens4)}")
    if k == 1:
        check(four == [8, 8] and seam2 == [10, 6] and eight[0] == 128 and seam4[0] == 136,
              f"modulo 2: four squares {four} against {seam2}; zeros of eight squares {eight[0]} against {seam4[0]}")
    if k == 2:
        check(four == [32, 64, 96, 64] and seam2 == [88, 48, 72, 48], f"modulo 4: {four} against {seam2}")
    if k == 3:
        check(four == [128, 512, 768, 512, 384, 512, 768, 512] and seam2 == [736, 384, 576, 384, 672, 384, 576, 384], f"modulo 8: {four} against {seam2}")

# ---------------------------------------------------------------- [3] classical counts
print("[3] Jacobi, Hurwitz, eight squares (direct counts)")
LIM = 600
r2 = [0] * (LIM + 1)
B = int(LIM ** 0.5) + 1
for a in range(-B, B + 1):
    for b in range(-B, B + 1):
        n = a * a + b * b
        if n <= LIM:
            r2[n] += 1
sig1 = [0] * (LIM + 1)
sig3 = [0] * (LIM + 1)
sodd = [0] * (LIM + 1)
for d in range(1, LIM + 1):
    for kk in range(d, LIM + 1, d):
        sig1[kk] += d
        sig3[kk] += d ** 3
        if d % 2:
            sodd[kk] += d
r4 = [sum(r2[a] * r2[n - a] for a in range(n + 1)) for n in range(LIM + 1)]
check(all(r4[n] == 8 * sig1[n] for n in range(1, LIM + 1, 2)), f"r_4(n) = sum_(a+b=n) r_2(a) r_2(b) = 8 sigma_1(n), odd n <= {LIM}")
check(all(r4[n] == 8 * sum(d for d in range(1, n + 1) if n % d == 0 and d % 4) for n in range(1, LIM + 1)), f"r_4(n) = 8 sum_(d | n, 4 not | d) d, n <= {LIM}")
R8 = 300
r8 = [sum(r4[a] * r4[n - a] for a in range(n + 1)) for n in range(R8 + 1)]
check(all(r8[n] == 16 * sum((-1) ** (n + d) * d ** 3 for d in range(1, n + 1) if n % d == 0) for n in range(1, R8 + 1)),
      f"r_8(n) = 16 sum_(d | n) (-1)^(n+d) d^3, n <= {R8}")
check(all(r8[n] == 16 * sig3[n] for n in range(1, R8 + 1, 2)) and r8[2] == 112 and 16 * sig3[2] == 144,
      f"r_8(n) = 16 sigma_3(n) for odd n <= {R8}; r_8(2) = {r8[2]}, 16 sigma_3(2) = {16 * sig3[2]}")
check(all(r8[n] == 16 * sig3[n] - 32 * sum(d ** 3 for d in range(1, n + 1, 2) if n % d == 0) and r8[n] < 16 * sig3[n] for n in range(2, R8 + 1, 2)),
      f"even n <= {R8}: r_8(n) = 16 sigma_3(n) - 32 sum_(d | n, d odd) d^3 < 16 sigma_3(n)")


def hurwitz(n):
    cnt = r4[n]
    B2 = int((4 * n) ** 0.5) + 1
    odd = [v for v in range(-B2, B2 + 1) if v % 2]
    sq = {}
    for a in odd:
        for b in odd:
            sq[a * a + b * b] = sq.get(a * a + b * b, 0) + 1
    return cnt + sum(c * sq.get(4 * n - kk, 0) for kk, c in sq.items())


HL = 150
hw = [0] + [hurwitz(n) for n in range(1, HL + 1)]
check(all(hw[n] == 24 * sodd[n] for n in range(1, HL + 1)) and all(hw[n] == 24 * sig1[n] for n in range(1, HL + 1, 2)),
      f"Hurwitz order: 24 * (sum of the odd divisors) for n <= {HL}; = 24 sigma_1(n) at odd n")

# ---------------------------------------------------------------- [4] densities
print("[4] (3)(d): 2-adic densities at the odd classes, Bernoulli ladder")
q = 16
g = np.arange(q, dtype=np.int64)
A, Bb, Cc, D = np.meshgrid(g, g, g, g, indexing="ij")
lip = np.bincount(((A * A + Bb * Bb + Cc * Cc + D * D) % q).ravel(), minlength=q)
hur = np.bincount(((A * A + Bb * Bb + Cc * Cc + D * D + A * D + Bb * D + Cc * D) % q).ravel(), minlength=q)
sm = np.bincount(((A * Bb + Cc * D) % q).ravel(), minlength=q)
dl = {F(int(lip[c]) * q, q ** 4) for c in range(1, q, 2)}
dh = {F(int(hur[c]) * q, q ** 4) for c in range(1, q, 2)}
ds = {F(int(sm[c]) * q, q ** 4) for c in range(1, q, 2)}
check(dl == {F(1)} and dh == {F(3, 2)} and ds == {F(3, 4)}, f"modulo 16, odd classes: Lipschitz {sorted(dl)}, Hurwitz order {sorted(dh)}, seam {sorted(ds)}")
# Hurwitz basis check: a + b i + c j + d (1 + i + j + k)/2 has norm a^2 + b^2 + c^2 + d^2 + ad + bd + cd
okn = True
for (a, b, c, d) in itertools.product(range(-2, 3), repeat=4):
    n4 = (2 * a + d) ** 2 + (2 * b + d) ** 2 + (2 * c + d) ** 2 + d ** 2
    okn &= n4 == 4 * (a * a + b * b + c * c + d * d + a * d + b * d + c * d)
check(okn, "norm form of the Hurwitz order in the basis 1, i, j, (1 + i + j + k)/2: a^2 + b^2 + c^2 + d^2 + ad + bd + cd")
cntH = {}
for (a, b, c, d) in itertools.product(range(-9, 10), repeat=4):
    n = a * a + b * b + c * c + d * d + a * d + b * d + c * d
    if n <= 12:
        cntH[n] = cntH.get(n, 0) + 1
check(all(cntH[n] == hw[n] for n in range(1, 13)), f"the same form counted directly, n <= 12: {[cntH[n] for n in range(1, 13)]}")
bern = {2: F(1, 6), 4: F(-1, 30), 6: F(1, 42), 8: F(-1, 30)}
vals = []
okb = True
for r in (2, 4, 6, 8):
    zeta = sum(j ** -r for j in range(1, 400001)) + (400000.5 ** (1 - r)) / (r - 1)
    lhs = (2 * math.pi) ** r / (math.factorial(r - 1) * zeta)
    rhs = 2 * r / abs(bern[r])
    vals.append(int(rhs))
    okb &= abs(lhs - float(rhs)) < 1e-4
check(okb and vals == [24, 240, 504, 480], f"(2 pi)^r/((r-1)! zeta(r)) = 2r/|B_r| = {vals} for r = 2, 4, 6, 8")
check(abs((8 * math.pi ** 4 / 3) / (math.pi ** 4 / 90) - 240) < 1e-9 and abs((2 * math.pi) ** 4 / 6 - 8 * math.pi ** 4 / 3) < 1e-9,
      "(2 pi)^4/3! = 8 pi^4/3 and (8 pi^4/3)/zeta(4) = 240")
check(abs(2 * math.pi ** 2 * 1.5 * (8 / math.pi ** 2) - 24) < 1e-12 and abs(4 * math.pi ** 2 * (6 / math.pi ** 2) - 24) < 1e-12,
      "Hurwitz order: 2 pi^2 n * (3/2) * (8/pi^2) sigma_1(n)/n = 24 sigma_1(n) = 4 pi^2 n * sigma_1(n)/(zeta(2) n)")

# ---------------------------------------------------------------- [5] the step of (3)(d)
print("[5] step of (3)(d): zero counts of quaternary quadratic forms by rank")


def rank_mod_p(M, p):
    M = [row[:] for row in M]
    n = len(M)
    rk = 0
    for col in range(n):
        piv = next((r for r in range(rk, n) if M[r][col] % p), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][col], -1, p)
        M[rk] = [(x * inv) % p for x in M[rk]]
        for r in range(n):
            if r != rk and M[r][col] % p:
                f = M[r][col]
                M[r] = [(x - f * y) % p for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk


def det_mod_p(M, p):
    n = len(M)
    M = [[x % p for x in row] for row in M]
    det = 1
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col]), None)
        if piv is None:
            return 0
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
            det = -det
        det = (det * M[col][col]) % p
        inv = pow(M[col][col], -1, p)
        for r in range(col + 1, n):
            f = (M[r][col] * inv) % p
            M[r] = [(x - f * y) % p for x, y in zip(M[r], M[col])]
    return det % p


pairs = [(i, j) for i in range(4) for j in range(i, 4)]
for p in (3,):
    g = np.arange(p, dtype=np.int64)
    X = np.meshgrid(g, g, g, g, indexing="ij")
    mono = np.stack([(X[i] * X[j]).ravel() for (i, j) in pairs])           # 10 x p^4
    seam_zero = p ** 3 + p ** 2 - p
    by_rank = {}
    n_seam = 0
    seam_all_nondeg_square = True
    squares_mod_p = {(x * x) % p for x in range(1, p)}
    for coef in itertools.product(range(p), repeat=10):
        vals = (np.array(coef, dtype=np.int64) @ mono) % p
        z = int((vals == 0).sum())
        # polar (Gram) matrix: 2 c_ii on the diagonal, c_ij off it
        M = [[0] * 4 for _ in range(4)]
        for (i, j), c in zip(pairs, coef):
            if i == j:
                M[i][i] = (2 * c) % p
            else:
                M[i][j] = M[j][i] = c
        rk = rank_mod_p(M, p)
        by_rank.setdefault(rk, set()).add(z)
        if z == seam_zero:
            n_seam += 1
            d = det_mod_p(M, p)
            if rk != 4 or d not in squares_mod_p:
                seam_all_nondeg_square = False
    expect = {0: {p ** 4}, 1: {p ** 3}, 2: {2 * p ** 3 - p ** 2, p ** 2}, 3: {p ** 3}, 4: {p ** 3 + p ** 2 - p, p ** 3 - p ** 2 + p}}
    check(by_rank == expect, f"p = {p}: all {p ** 10} forms; zero counts by rank {dict(sorted((k, sorted(v)) for k, v in by_rank.items()))}")
    check(n_seam == 21060 and seam_all_nondeg_square,
          f"p = {p}: forms with the seam's number of zeros ({seam_zero}): {n_seam}, all non-degenerate with polar determinant a square")
import random
random.seed(5135)
for p in (5, 7):
    g = np.arange(p, dtype=np.int64)
    X = np.meshgrid(g, g, g, g, indexing="ij")
    mono = np.stack([(X[i] * X[j]).ravel() for (i, j) in pairs])
    by_rank = {}
    okseam = True
    squares_mod_p = {(x * x) % p for x in range(1, p)}
    forms = [tuple(random.randrange(p) for _ in range(10)) for _ in range(3000)]
    # add degenerate forms on purpose: forms in fewer variables after a random linear change
    for s in (0, 1, 2, 3):
        for _ in range(300):
            T = [[random.randrange(p) for _ in range(4)] for _ in range(4)]
            base = {(i, j): (random.randrange(p) if j < s else 0) for (i, j) in pairs}
            # Q(x) = base(T x): expand
            coef = {pq: 0 for pq in pairs}
            for (i, j), c in base.items():
                if c:
                    for a in range(4):
                        for b in range(4):
                            key = (min(a, b), max(a, b))
                            coef[key] = (coef[key] + c * T[i][a] * T[j][b]) % p
            forms.append(tuple(coef[pq] for pq in pairs))
    for coef in forms:
        vals = (np.array(coef, dtype=np.int64) @ mono) % p
        z = int((vals == 0).sum())
        M = [[0] * 4 for _ in range(4)]
        for (i, j), c in zip(pairs, coef):
            if i == j:
                M[i][i] = (2 * c) % p
            else:
                M[i][j] = M[j][i] = c
        rk = rank_mod_p(M, p)
        by_rank.setdefault(rk, set()).add(z)
        if z == p ** 3 + p ** 2 - p:
            okseam &= (rk == 4 and det_mod_p(M, p) in squares_mod_p)
    expect = {0: {p ** 4}, 1: {p ** 3}, 2: {2 * p ** 3 - p ** 2, p ** 2}, 3: {p ** 3}, 4: {p ** 3 + p ** 2 - p, p ** 3 - p ** 2 + p}}
    check(all(by_rank.get(k, set()) <= expect[k] for k in range(5)) and set(by_rank) == {0, 1, 2, 3, 4} and okseam,
          f"p = {p}: {len(forms)} sampled forms (all ranks): zero counts by rank {dict(sorted((k, sorted(v)) for k, v in by_rank.items()))}; "
          f"p^3 + p^2 - p only at rank 4 with square polar determinant")
check(all(2 * p ** 3 - p ** 2 - (p ** 3 + p ** 2 - p) == p * (p - 1) ** 2 for p in range(2, 200)), "2 p^3 - p^2 - (p^3 + p^2 - p) = p (p - 1)^2")
# F_2
g = np.arange(2, dtype=np.int64)
X = np.meshgrid(g, g, g, g, indexing="ij")
mono = np.stack([(X[i] * X[j]).ravel() for (i, j) in pairs])
n10 = 0
n10_nondeg = 0
zero_hist = {}
for coef in itertools.product(range(2), repeat=10):
    vals = (np.array(coef, dtype=np.int64) @ mono) % 2
    z = int((vals == 0).sum())
    M = [[0] * 4 for _ in range(4)]
    for (i, j), c in zip(pairs, coef):
        if i != j:
            M[i][j] = M[j][i] = c
    nondeg = det_mod_p(M, 2) == 1
    zero_hist[(z, nondeg)] = zero_hist.get((z, nondeg), 0) + 1
    if z == 10:
        n10 += 1
        n10_nondeg += nondeg
check(n10 == 280 and n10_nondeg == 280, f"F_2: 1024 forms; with 10 zeros: {n10}, of which with non-degenerate polar form: {n10_nondeg}")
print(f"       F_2 histogram (zeros, polar form non-degenerate) -> forms: {dict(sorted(zero_hist.items()))}")
check(20160 // 72 == 280 and 20160 % 72 == 0, "280 = |GL_4(2)| / |O+_4(2)| = 20160 / 72")

# ---------------------------------------------------------------- [6] Ingham
print("[6] Ingham: sum_(m < n) d(m) d(n - m) against (6/pi^2) sigma_1(n) ln^2 n")
TOP = 1000000
dv = np.zeros(TOP + 1, dtype=np.int64)
for d in range(1, TOP + 1):
    dv[d::d] += 1
for n, want in ((720720, 0.70), (999983, 1.17)):
    Dn = int(np.dot(dv[1:n], dv[n - 1:0:-1]))
    s1 = sum(d + n // d for d in range(1, int(n ** 0.5) + 1) if n % d == 0 and d * d != n) + (int(n ** 0.5) if int(n ** 0.5) ** 2 == n else 0)
    ratio = Dn / ((6 / math.pi ** 2) * s1 * math.log(n) ** 2)
    check(abs(ratio - want) < 0.005, f"n = {n}: count {Dn}, ratio to the main term {ratio:.4f} (quoted {want})")
small = True
for n in range(2, 41):
    brute = sum(1 for a1 in range(1, n) for a2 in range(1, n) if a1 * a2 < n for b1 in range(1, n) if (n - a1 * a2) % b1 == 0)
    small &= brute == int(np.dot(dv[1:n], dv[n - 1:0:-1]))
check(small, "the sum is the number of solutions of a1 a2 + b1 b2 = n in positive integers (n <= 40)")

print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
