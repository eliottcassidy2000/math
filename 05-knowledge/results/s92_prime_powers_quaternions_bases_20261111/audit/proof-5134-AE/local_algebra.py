# local_algebra.py -- auditor proof-5134-AE.  Exact checks of the local algebra used in THM-5134 (A), (C), (D), (E).
#  1. c_d m^2 = J_2(m/d) = #{(alpha, beta) mod m : gcd(alpha, beta, m) = d}, sum_d c_d = 1   (brute force, odd m <= 135, and 225)
#  2. the linear form w -> a w1 + b w2 on (Z/m)^2 takes every multiple of d = gcd(a, b, m) on exactly m d vectors
#     and no other value (brute force over ALL (a, b), odd m <= 45 and m = 49, 63, 75, 81)
#  3. m^3 L_m(c) = #{(a1, a2, b1, b2) mod m : a1 a2 + b1 b2 = c}  for every odd m <= 135 and m = 225 (all classes),
#     and with a unit u: a1 a2 + u b1 b2 (u = 2, 4: the F1 form with s = 1, 2)
#  4. L_m multiplicative; prime-power formulas (v < k and v >= k); mean of L_m over the classes = 1; the mean of
#     L_(p^k) over the classes divisible by p is 1 - 1/p^2 + 1/p
#  5. even moduli (control): the formula fails for m = 2, 4, 6 (so "odd" is needed)
import math, sys
from fractions import Fraction as Fr
import numpy as np


def primes_of(n):
    ps = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        ps.append(n)
    return ps


def c_d(m, d):
    v = Fr(1, d * d)
    for p in primes_of(m // d):
        v *= (1 - Fr(1, p * p))
    return v


def L_m(m, e):
    g = math.gcd(m, e % m) if e % m else m
    return sum(d * c_d(m, d) for d in range(1, m + 1) if g % d == 0)


def vp(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def L_formula(m, e):
    # product over p^k || m of the stated prime-power formulas
    val = Fr(1)
    mm = m
    for p in primes_of(m):
        k = vp(m, p)
        e0 = e % m
        v = k if e0 % (p ** k) == 0 else vp(e0, p)
        if v < k:
            val *= (1 + Fr(1, p)) * (1 - Fr(1, p ** (v + 1)))
            sig = sum(Fr(p ** i) for i in range(v + 1))
            assert (1 + Fr(1, p)) * (1 - Fr(1, p ** (v + 1))) == (1 - Fr(1, p * p)) * sig / p ** v
        else:
            val *= 1 + Fr(1, p) - Fr(1, p ** (k + 1))
            sig = sum(Fr(p ** i) for i in range(k))
            assert 1 + Fr(1, p) - Fr(1, p ** (k + 1)) == (1 - Fr(1, p * p)) * sig / p ** (k - 1) + Fr(1, p ** k)
    return val


def count_form(m, u=1):
    # number of (a1, a2, b1, b2) mod m with a1 a2 + u b1 b2 = c, for every c: additive convolution of the product law
    a = np.arange(m)
    prod = (a[:, None] * a[None, :]) % m
    t = np.bincount(prod.ravel(), minlength=m).astype(object)          # #{(a1, a2): a1 a2 = r}
    tu = [0] * m
    for r in range(m):
        tu[(u * r) % m] += int(t[r])
    out = [0] * m
    for r1 in range(m):
        if t[r1]:
            for r2 in range(m):
                out[(r1 + r2) % m] += int(t[r1]) * tu[r2]
    return out


def main():
    odd = [m for m in range(3, 136, 2)] + [225]
    # 1
    for m in odd:
        cnt = {}
        for al in range(m):
            for be in range(m):
                d = math.gcd(math.gcd(al, be), m)
                cnt[d] = cnt.get(d, 0) + 1
        for d in range(1, m + 1):
            if m % d == 0:
                assert cnt.get(d, 0) == c_d(m, d) * m * m, (m, d)
        assert sum(c_d(m, d) for d in range(1, m + 1) if m % d == 0) == 1
    print("1. c_d m^2 = J_2(m/d) = #{gcd(alpha, beta, m) = d} and sum_d c_d = 1: exact for the %d odd moduli 3..135 and 225" % len(odd))
    # 2
    tested = 0
    for m in [mm for mm in range(3, 46, 2)] + [49, 63, 75, 81]:
        w1 = np.arange(m)[:, None]
        w2 = np.arange(m)[None, :]
        for a in range(m):
            for b in range(m):
                d = math.gcd(math.gcd(a, b), m)
                h = np.bincount(((a * w1 + b * w2) % m).ravel(), minlength=m)
                for c in range(m):
                    assert h[c] == (m * d if c % d == 0 else 0), (m, a, b, c)
                tested += 1
    print("2. fibres of w -> a w1 + b w2 mod m: every multiple of d = gcd(a, b, m) on exactly m d vectors, nothing else: exact, %d forms (all (a, b); odd m <= 45 and 49, 63, 75, 81)" % tested)
    # 3
    n3 = 0
    for m in odd:
        for u in (1, 2, 4):
            if math.gcd(u, m) != 1:
                continue
            cf = count_form(m, u)
            assert sum(cf) == m ** 4
            for c in range(m):
                assert cf[c] == m ** 3 * L_m(m, c), (m, u, c)
            n3 += 1
    print("3. #{a1 a2 + u b1 b2 = c mod m} = m^3 L_m(c) for every class c: exact, %d (m, u) pairs (odd m <= 135 and 225; u = 1, 2, 4)" % n3)
    # 4
    for m in odd + [243, 343, 375, 441, 625, 675, 729, 1125]:
        tot = Fr(0)
        for e in range(m):
            assert L_m(m, e) == L_formula(m, e), (m, e)
            tot += L_m(m, e)
        assert tot == m
        for p in primes_of(m):
            k = vp(m, p)
            pk = p ** k
            # mean of L_(p^k) over the classes divisible by p
            s = sum(L_m(pk, e) for e in range(0, pk, p))
            assert s / (pk // p) == 1 - Fr(1, p * p) + Fr(1, p)
    print("4. L_m = product of the prime-power formulas (v < k: (1+1/p)(1-p^(-v-1)) = (1-p^-2) sigma(p^v)/p^v; v >= k: 1+1/p-p^(-k-1)); mean over classes = 1; mean of L_(p^k) over p | e is 1-1/p^2+1/p: exact")
    print("   values:  L_9 = %s   L_27 = %s   L_25 = %s" % ([str(L_m(9, e)) for e in (1, 3, 0)], [str(L_m(27, e)) for e in (1, 3, 9, 0)], [str(L_m(25, e)) for e in (1, 5, 0)]))
    print("   floats:  L_9 = %s   L_27 = %s   L_25 = %s   L_49 = %s   L_45(g0=1,3,5,9,15,45) = %s" % (
        ["%.4f" % float(L_m(9, e)) for e in (1, 3, 0)], ["%.4f" % float(L_m(27, e)) for e in (1, 3, 9, 0)],
        ["%.4f" % float(L_m(25, e)) for e in (1, 5, 0)], ["%.4f" % float(L_m(49, e)) for e in (1, 7, 0)],
        ["%.4f" % float(L_m(45, e)) for e in (1, 3, 5, 9, 15, 0)]))
    # 5
    for m in (2, 4, 6, 8, 12):
        cf = count_form(m, 1)
        bad = [c for c in range(m) if cf[c] != m ** 3 * L_m(m, c)]
        print("5. control, even m = %d, form a1 a2 + b1 b2 with ALL residues: classes where count != m^3 L_m: %s (the formula is a statement about uniform states; it holds for even m too: %s)" % (m, bad, len(bad) == 0))


if __name__ == "__main__":
    main()
