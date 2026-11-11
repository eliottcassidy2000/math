# Auditor numerics-prereg: independent exact checks of the uniform-state identities quoted in THM-5134 (D), (F)
# and THM-5135 (1), (2) (own code; cyclic convolutions of exact integer distributions).
import math
from fractions import Fraction as Fr
import numpy as np
from aud_common import divisors


def cconv(a, b):
    q = len(a)
    out = np.zeros(q, dtype=object)
    for i in range(q):
        if a[i]:
            out += np.roll(b, i) * a[i]
    return out


def table(q, u=1):
    d = np.zeros(q, dtype=object)
    for a in range(q):
        for b in range(q):
            d[(u * a * b) % q] += 1
    return d


def squares(q):
    d = np.zeros(q, dtype=object)
    for a in range(q):
        d[(a * a) % q] += 1
    return d


def power(d, r):
    out = d
    for _ in range(r - 1):
        out = cconv(out, d)
    return out


def L_exact(m, e):
    g = math.gcd(m, e % m) if e % m else m
    tot = Fr(0)
    for d in divisors(g):
        t = Fr(1, d)
        r = m // d
        p = 3
        while r > 1:
            if r % p == 0:
                t *= 1 - Fr(1, p * p)
                while r % p == 0:
                    r //= p
            p += 2
        tot += t
    return tot


ok = True
print("(i) additive square of the multiplication table modulo odd m equals m^3 L_m(c) (THM-5134 (D), THM-5135 (1)); also with the unit 2")
for m in (9, 15, 21, 25, 27, 33, 35, 45, 49, 63, 75, 81, 105, 121, 125, 135, 147, 175, 225):
    t = table(m)
    sq = cconv(t, t)
    sq2 = cconv(t, table(m, 2))
    good = all(Fr(int(sq[c])) == m ** 3 * L_exact(m, c) for c in range(m)) and all(int(sq[c]) == int(sq2[c]) for c in range(m))
    ok &= good
    print("    m = %3d: %s" % (m, good))

print("(ii) class counts of a1 a2 + 2 b1 b2 and of four squares modulo q (THM-5134 (F))")
for q in (9, 27, 81, 25, 125, 49, 121, 169):
    a = cconv(table(q), table(q, 2))
    b = power(squares(q), 4)
    good = all(int(a[c]) == int(b[c]) for c in range(q))
    ok &= good
    print("    q = %3d: equal for every class: %s" % (q, good))

print("(iii) THM-5135 (2): sum of 2r squares against the hyperbolic form of rank r, class counts modulo q; expected equal iff r even or p = 1 (mod 4)")
for q, p in ((25, 5), (27, 3), (49, 7), (125, 5), (121, 11), (169, 13)):
    cells = []
    for r in (1, 2, 3, 4):
        a = power(squares(q), 2 * r)
        b = power(table(q), r)
        same = all(int(a[c]) == int(b[c]) for c in range(q))
        expect = (r % 2 == 0) or (p % 4 == 1)
        ok &= (same == expect)
        cells.append("r=%d: equal=%s (expected %s)" % (r, same, expect))
    print("    q = %3d: %s" % (q, "; ".join(cells)))

print("(iv) THM-5135 (1): p^k P(u.v = c) for the rank-r seam, uniform states, against (1 - p^-r) sum_{i<=v} p^(-(r-1) i), v < k")
for q, p, k in ((9, 3, 2), (27, 3, 3), (25, 5, 2), (49, 7, 2), (8, 2, 3), (16, 2, 4), (32, 2, 5)):
    t = table(q)
    for r in (1, 2, 3, 4):
        d = power(t, r)
        tot = q ** (2 * r)
        good = True
        for c in range(1, q):
            v = 0
            cc = c
            while cc % p == 0:
                v += 1
                cc //= p
            want = (1 - Fr(1, p ** r)) * sum(Fr(1, p ** ((r - 1) * i)) for i in range(v + 1))
            good &= Fr(int(d[c]) * q, tot) == want
        ok &= good
        print("    p^k = %2d, r = %d: %s" % (q, r, good))
print("ALL MY CHECKS PASSED" if ok else "SOME CHECK FAILED")
