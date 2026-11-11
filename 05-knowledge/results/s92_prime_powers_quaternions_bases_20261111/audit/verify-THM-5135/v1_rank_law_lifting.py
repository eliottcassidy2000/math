"""Verifier of THM-5135 (S92 audit follow-up): statement (1) with the class 0, the table identity for every m,
the lifting lemma of Proofs (3) and its use for E8.  Own code, exact integers / Fractions.
Run: python v1_rank_law_lifting.py > v1_rank_law_lifting.out
"""
import itertools
import math
import random
from fractions import Fraction as F

import numpy as np

FAILS = []


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def vp(n, p):
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def table(q, scale=1):
    d = [0] * q
    for a in range(q):
        for b in range(q):
            d[(scale * a * b) % q] += 1
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


def sigma(j, n):
    return sum(d ** j for d in range(1, n + 1) if n % d == 0)


def phi(n):
    return sum(1 for a in range(1, n + 1) if math.gcd(a, n) == 1)


def primes_of(n):
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out


def L(m, e):
    g = math.gcd(m, e) if e % m else m
    tot = F(0)
    for d in range(1, g + 1):
        if g % d == 0:
            t = F(1, d)
            for p in primes_of(m // d):
                t *= 1 - F(1, p * p)
            tot += t
    return tot


# ---------------------------------------------------------------- [1] statement (1)
print("[1] statement (1): rank-r law modulo p^k, both displayed values (classes of valuation v < k and the class 0)")
for (p, k) in ((2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3), (5, 1), (5, 2), (7, 1), (7, 2), (11, 1)):
    q = p ** k
    tab = table(q)
    for r in range(1, 6):
        if q ** 1 > 30 and r > 4:
            continue
        cnt = conv_power(tab, r)
        ok = True
        for c in range(q):
            fac = F(cnt[c] * q, q ** (2 * r))
            v = vp(c, p) if c else k
            if v < k:
                want = (1 - F(1, p ** r)) * F(sigma(r - 1, p ** v), p ** ((r - 1) * v))
                want2 = (1 - F(1, p ** r)) * sum(F(1, p ** ((r - 1) * i)) for i in range(v + 1))
                ok &= (fac == want == want2)
            else:
                want = (1 - F(1, p ** r)) * F(sigma(r - 1, p ** (k - 1)), p ** ((r - 1) * (k - 1))) + F(1, p ** ((r - 1) * k))
                ok &= (fac == want)
        check(ok, f"p^k = {q:3d}, r = {r}: both formulas hold at every class (class 0 factor {F(cnt[0] * q, q ** (2 * r))})")
# r = 1 closed form and the value quoted in the editor's report (35/27 at 9, r = 2)
check(all(F(table(p ** k)[p ** v % p ** k] * p ** k, p ** (2 * k)) == (1 - F(1, p)) * (v + 1)
          for (p, k) in ((2, 4), (3, 3), (5, 2)) for v in range(k)), "r = 1: factor (1 - 1/p)(v + 1) at a class of valuation v < k")
check(F(conv_power(table(9), 2)[0] * 9, 9 ** 4) == F(35, 27), "class 0 modulo 9, r = 2: factor 35/27")
check(F(conv_power(table(2), 4)[0] * 2, 2 ** 8) == F(17, 16) and F(conv_power(table(2), 4)[0], 2 ** 8) == F(17, 32),
      "class 0 modulo 2, r = 4: factor 17/16 = 15/16 + 1/8, share 17/32 (136 of 256)")
check(F(conv_power(table(3), 4)[0] * 3, 3 ** 8) == F(83, 81), "class 0 modulo 3, r = 4: factor 83/81 = 1 + 3^-3 - 3^-4")

# ---------------------------------------------------------------- [2] table identity for every m
print("[2] N_m(t) = sum_{d | gcd(t, m)} d phi(m/d); (N_m * N_m)(c) = m^3 L_m(c) for every m; unit twists; L_(p^k) closed forms")
bad_t = bad_sq = bad_tw = 0
MMAX = 96
for m in range(1, MMAX + 1):
    tab = table(m)
    for t in range(m):
        g = math.gcd(t, m) if t else m
        if tab[t] != sum(d * phi(m // d) for d in range(1, g + 1) if g % d == 0):
            bad_t += 1
    sq = conv(tab, tab)
    for c in range(m):
        if F(sq[c]) != m ** 3 * L(m, c):
            bad_sq += 1
    for u in (2, 3, 5, 7):
        if math.gcd(u, m) == 1 and m > 1:
            if conv(tab, table(m, u)) != sq:
                bad_tw += 1
check(bad_t == 0, f"N_m(t) formula, every m <= {MMAX}: mismatches {bad_t}")
check(bad_sq == 0, f"(N_m * N_m)(c) = m^3 L_m(c), every m <= {MMAX} (even m included): mismatches {bad_sq}")
check(bad_tw == 0, f"same with one table multiplied by a unit (2, 3, 5, 7 when prime to m): mismatches {bad_tw}")
okL = True
for (p, k) in ((2, 1), (2, 3), (2, 5), (3, 1), (3, 2), (3, 4), (5, 3), (7, 2)):
    q = p ** k
    for e in range(q):
        v = vp(e, p) if e else k
        if v < k:
            okL &= L(q, e) == (1 - F(1, p * p)) * F(sigma(1, p ** v), p ** v) == (1 + F(1, p)) * (1 - F(1, p ** (v + 1)))
        else:
            okL &= L(q, e) == (1 - F(1, p * p)) * F(sigma(1, p ** (k - 1)), p ** (k - 1)) + F(1, p ** k) == 1 + F(1, p) - F(1, p ** (k + 1))
check(okL, "L_(p^k)(e) = (1 - p^-2) sigma_1(p^v)/p^v (v < k), = 1 + 1/p - p^(-k-1) (v >= k)")
check(all(L(2 ** k, c) == F(3, 4) for k in range(1, 7) for c in range(1, 2 ** k, 2)), "at 2^k the factor at the odd classes is 3/4 (k <= 6)")
# Euler product for r >= 2 (local factors of an integer n)
okE = True
for r in (2, 3, 4):
    for n in (1, 2, 12, 45, 360):
        prod = 1.0
        sieve = bytearray([1]) * 20001
        for pp in range(2, 20001):
            if sieve[pp]:
                for j in range(pp * pp, 20001, pp):
                    sieve[j] = 0
                v = vp(n, pp)
                prod *= (1 - pp ** -r) * sigma(r - 1, pp ** v) / pp ** ((r - 1) * v)
        zeta = sum(j ** -r for j in range(1, 200001)) + (200000.5 ** (1 - r)) / (r - 1)
        okE &= abs(prod - sigma(r - 1, n) / (zeta * n ** (r - 1))) < 2e-4
check(okE, "Euler product of the factors = sigma_(r-1)(n)/(zeta(r) n^(r-1)) for r = 2, 3, 4 (primes < 20000, tolerance 2e-4)")

# ---------------------------------------------------------------- [3] lifting lemma
print("[3] lifting lemma: N_k(c) = p^((n-1)(k-1)) (N_1(c mod p) - [p | c]) + [p^2 | c] p^n N_(k-2)(c/p^2), k >= 2")


def class_counts(coef, n, q):
    """coef[(i, j)] (i <= j) integer coefficients of Q; returns the class counts of Q modulo q over (Z/q)^n."""
    grids = np.meshgrid(*[np.arange(q, dtype=np.int64)] * n, indexing="ij")
    val = np.zeros_like(grids[0])
    for (i, j), c in coef.items():
        if c:
            val = (val + c * grids[i] * grids[j]) % q
    return np.bincount(val.ravel(), minlength=q).astype(object)


def polar_det(coef, n):
    M = [[0] * n for _ in range(n)]
    for (i, j), c in coef.items():
        if i == j:
            M[i][i] = 2 * c
        else:
            M[i][j] = M[j][i] = c
    # exact integer determinant (fraction-free via Fractions)
    A = [[F(x) for x in row] for row in M]
    det = F(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if A[r][i] != 0), None)
        if piv is None:
            return 0
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            det = -det
        det *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            for c in range(i, n):
                A[r][c] -= f * A[i][c]
    return int(det)


def recursion(N1, p, n, K):
    """class counts modulo p^k for k = 0..K from those modulo p, by the lemma."""
    out = {0: [1], 1: list(N1)}
    for k in range(2, K + 1):
        q = p ** k
        prev = out[k - 2]
        cur = []
        for c in range(q):
            val = p ** ((n - 1) * (k - 1)) * (N1[c % p] - (1 if c % p == 0 else 0))
            if c % (p * p) == 0:
                val += p ** n * prev[(c // (p * p)) % (p ** (k - 2))]
            cur.append(val)
        out[k] = cur
    return out


random.seed(92)
tested = 0
viol = 0
controls_bad = 0
controls = 0
for (p, n, K, trials) in ((2, 2, 6, 12), (2, 4, 4, 12), (3, 2, 4, 10), (3, 3, 3, 10), (3, 4, 3, 6), (5, 2, 3, 8), (5, 3, 2, 6), (7, 2, 2, 6)):
    done = 0
    while done < trials:
        coef = {(i, j): random.randint(-6, 6) for i in range(n) for j in range(i, n)}
        det = polar_det(coef, n)
        if det % p == 0:
            # negative control: the recursion is not claimed; record whether it happens to fail
            if controls < 40:
                N1 = class_counts(coef, n, p)
                rec = recursion(N1, p, n, 2)
                controls += 1
                if list(class_counts(coef, n, p * p)) != rec[2]:
                    controls_bad += 1
            continue
        N1 = class_counts(coef, n, p)
        rec = recursion(N1, p, n, K)
        for k in range(2, K + 1):
            tested += 1
            if list(class_counts(coef, n, p ** k)) != rec[k]:
                viol += 1
        done += 1
check(viol == 0, f"random integer-valued forms with polar determinant prime to p (p = 2, 3, 5, 7; n = 2, 3, 4): {tested} (form, k) cases, violations {viol}")
print(f"       (control: forms with p | polar determinant: the recursion fails for {controls_bad} of {controls}; it is not claimed there)")

# the two hyperbolic seams as quadratic forms
for (p, r, K) in ((2, 2, 5), (2, 3, 3), (3, 2, 3), (5, 2, 2)):
    coef = {(i, j): 0 for i in range(2 * r) for j in range(i, 2 * r)}
    for i in range(r):
        coef[(2 * i, 2 * i + 1)] = 1
    N1 = class_counts(coef, 2 * r, p)
    rec = recursion(N1, p, 2 * r, K)
    ok = all(rec[k] == conv_power(table(p ** k), r) for k in range(1, K + 1))
    check(ok, f"seam of rank {r} at p = {p}: the recursion from the counts modulo p gives the r-fold table power modulo p^k, k <= {K}")

# ---------------------------------------------------------------- [4] E8
print("[4] E8: x^T C x / 2 against the rank-4 seam")
C = np.zeros((8, 8), dtype=np.int64)
for i in range(8):
    C[i, i] = 2
for (i, j) in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)):
    C[i, j] = C[j, i] = -1
detC = round(float(np.linalg.det(C.astype(float))))
check(detC == 1, f"Cartan matrix of E8: determinant {detC}")
coefE = {(i, j): 0 for i in range(8) for j in range(i, 8)}
for i in range(8):
    coefE[(i, i)] = 1
for (i, j) in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)):
    coefE[(i, j)] = -1
e8 = {}
for q in (2, 4):
    e8[q] = [int(x) for x in class_counts(coefE, 8, q)]
# modulo 8 in chunks over the first two coordinates
q = 8
cnt8 = np.zeros(q, dtype=np.int64)
grids = np.meshgrid(*[np.arange(q, dtype=np.int64)] * 6, indexing="ij")
rest = np.zeros_like(grids[0])
for (i, j), c in coefE.items():
    if c and i >= 2 and j >= 2:
        rest = (rest + c * grids[i - 2] * grids[j - 2]) % q
for x0 in range(q):
    for x1 in range(q):
        val = rest + x0 * x0 + x1 * x1 - x0 * x1
        # cross terms between (x0, x1) and the others: only (1, 2)
        val = (val - x1 * grids[0]) % q
        cnt8 += np.bincount(val.ravel(), minlength=q)
e8[8] = [int(x) for x in cnt8]
seam = {q: conv_power(table(q), 4) for q in (2, 4, 8)}
check(e8[2] == [136, 120] == seam[2], f"modulo 2: E8 {e8[2]} = seam {seam[2]}")
check(e8[4] == [17536, 15360, 17280, 15360] == seam[4], f"modulo 4: E8 {e8[4]} = seam")
check(e8[8] == [2246656, 1966080, 2211840, 1966080, 2242560, 1966080, 2211840, 1966080] == seam[8], f"modulo 8: E8 {e8[8]} = seam")
rec = recursion([136, 120], 2, 8, 12)
check(rec[2] == e8[4] and rec[3] == e8[8], "the recursion started from 136, 120 reproduces the enumerated counts modulo 4 and 8")
okk = True
for k in range(1, 13):
    q = 2 ** k
    for c in range(q):
        v = vp(c, 2) if c else k
        if v < k:
            want = (1 - F(1, 16)) * F(sigma(3, 2 ** v), 2 ** (3 * v))
        else:
            want = (1 - F(1, 16)) * F(sigma(3, 2 ** (k - 1)), 2 ** (3 * (k - 1))) + F(1, 2 ** (3 * k))
        okk &= F(rec[k][c] * q, q ** 8) == want
check(okk, "the recursion from 136, 120 gives the law (1) with r = 4 modulo 2^k for k <= 12, the class 0 included")
check([rec[4][c] for c in (1, 2, 4, 8, 0)] == [251658240, 283115520, 287047680, 287539200, 287604736],
      "modulo 16 by valuation 0..4: 251658240, 283115520, 287047680, 287539200, 287604736 (as in the audit's a3_e8.out)")
# contrasts modulo 2
coefA = {(i, j): 0 for i in range(8) for j in range(i, 8)}
for i in range(3):
    coefA[(2 * i, 2 * i + 1)] = 1
coefA[(6, 6)] = coefA[(7, 7)] = coefA[(6, 7)] = 1
zA = int(class_counts(coefA, 8, 2)[0])
coefZ = {(i, j): (1 if i == j else 0) for i in range(8) for j in range(i, 8)}
zZ = int(class_counts(coefZ, 8, 2)[0])
check(zA == 120 and zZ == 128, f"zeros modulo 2: three hyperbolic planes and x^2 + xy + y^2: {zA}; Z^8: {zZ}")
check(polar_det(coefE, 8) == 1 and polar_det({(i, j): (1 if (j == i + 1 and i % 2 == 0) else 0) for i in range(8) for j in range(i, 8)}, 8) == 1,
      "polar determinants: E8 form 1, rank-4 seam 1")

print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
