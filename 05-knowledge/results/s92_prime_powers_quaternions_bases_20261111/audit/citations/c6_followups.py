#!/usr/bin/env python3
"""Auditor 'citations' - follow-up checks behind individual findings.

F1  (N_m * N_m)(c) = m^3 L_m(c) for EVERY m <= 150, even m included (THM-5135 (1) states it for odd m only).
F2  The 2-adic classification behind THM-5135 (3)(b): among even unimodular forms of rank 8 over Z_2 the two classes
    H^4 and H^3 + A (A = x^2 + xy + y^2) are told apart modulo 2 already (136 against 120 zeros); E8 has 136.
    Z^8 (odd lattice, same rank and discriminant 1 as E8) has 128: 'rank and discriminant' does not classify at 2.
F3  Bicomplex orders: Z[i]^2 (maximal) against the Z-span of the signed units 1, i, j, ij: counts of norm n, n <= 60;
    for odd n the span has exactly half of the count of Z[i]^2; for n = 2 (mod 4) it has none.
F4  The state behind the early boost of 13 in base 5 (THM-5134 (H), last sentence): first 13-divisible state.
F5  Siegel's formula as a consistency check: r(n) = real density * prod_p local density for Z^4 (odd n) and the
    Hurwitz order (2-adic density 3/2 at odd n), from exact class counts modulo 2^k.
F6  Eisenstein normalisation: (2 pi)^r / ((r-1)! zeta(r)) = 2r/|B_r| for r = 2, 4, 6, 8 (r = 4 gives 240; r = 2 gives 24).
"""
import math
import itertools
from fractions import Fraction


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def table_counts(m):
    t = [0] * m
    for a in range(m):
        for b in range(m):
            t[(a * b) % m] += 1
    return t


def cconv(a, b, m):
    out = [0] * m
    for i in range(m):
        ai = a[i]
        if ai:
            for j in range(m):
                out[(i + j) % m] += ai * b[j]
    return out


def L_m(m, e):
    g = math.gcd(m, e % m) if e % m else m
    tot = Fraction(0)
    for d in divisors(g):
        q = m // d
        f = Fraction(1, d)
        pp = 2
        qq = q
        while qq > 1:
            if qq % pp == 0:
                f *= (1 - Fraction(1, pp * pp))
                while qq % pp == 0:
                    qq //= pp
            pp += 1
        tot += f
    return tot


def main():
    # F1
    bad = 0
    for m in range(1, 151):
        t = table_counts(m)
        sq = cconv(t, t, m)
        for c in range(m):
            if Fraction(sq[c], m ** 3) != L_m(m, c):
                bad += 1
    print("F1 (N_m * N_m)(c) = m^3 L_m(c) for every 1 <= m <= 150 (even m included): failures %d" % bad)

    # F2
    def zeros_mod2(form, nvar):
        z = 0
        for x in itertools.product((0, 1), repeat=nvar):
            if form(x) % 2 == 0:
                z += 1
        return z

    h4 = lambda x: x[0] * x[1] + x[2] * x[3] + x[4] * x[5] + x[6] * x[7]
    h3a = lambda x: x[0] * x[1] + x[2] * x[3] + x[4] * x[5] + x[6] * x[6] + x[6] * x[7] + x[7] * x[7]
    z8 = lambda x: sum(v * v for v in x)
    edges = [(0, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (1, 3)]
    e8 = lambda x: sum(v * v for v in x) - sum(x[a] * x[b] for a, b in edges)
    print("F2 zeros modulo 2 among the 256 vectors: H^4 %d, H^3 + A %d, E8 (x^T C x / 2) %d, Z^8 %d"
          % (zeros_mod2(h4, 8), zeros_mod2(h3a, 8), zeros_mod2(e8, 8), zeros_mod2(z8, 8)))
    print("   determinants: H^4: (-1)^4 = 1; H^3 + A: (-1)^3 * 3 = -3 (not a 2-adic square class of 1); E8: 1; Z^8: 1")

    # F3
    LIM = 60
    r2 = [0] * (LIM + 1)
    B = int(LIM ** 0.5) + 1
    for a in range(-B, B + 1):
        for b in range(-B, B + 1):
            if a * a + b * b <= LIM:
                r2[a * a + b * b] += 1
    full = [0] * (LIM + 1)
    span = [0] * (LIM + 1)
    gauss = [(a, b) for a in range(-B, B + 1) for b in range(-B, B + 1) if 0 < a * a + b * b <= LIM]
    for (a, b) in gauss:
        n1 = a * a + b * b
        for (c, d) in gauss:
            n = n1 * (c * c + d * d)
            if n <= LIM:
                full[n] += 1
                if (a - c) % 2 == 0 and (b - d) % 2 == 0:
                    span[n] += 1
    okc = all(full[n] == sum(r2[a] * r2[n // a] for a in divisors(n)) for n in range(1, LIM + 1))
    half = all(2 * span[n] == full[n] for n in range(1, LIM + 1, 2))
    none2 = all(span[n] == 0 for n in range(2, LIM + 1, 4))
    print("F3 Z[i]^2: count = sum_{ab=n} r_2(a) r_2(b) for n <= %d: %s" % (LIM, okc))
    print("   Z-span of {1, i, j, ij} = {(z1, z2) : z1 = z2 (mod 2)}: half of the Z[i]^2 count at every odd n: %s; "
          "zero at every n = 2 (mod 4): %s" % (half, none2))
    print("   n = 1..16: Z[i]^2 %s" % full[1:17])
    print("              span   %s" % span[1:17])
    print("   units: Z[i]^2 has %d (C_4 x C_4), the span has %d (C_4 x C_2)" % (full[1], span[1]))

    # F4
    def psi(h):
        if h == 0:
            return 0
        if h % 2 == 0:
            return psi(h // 2) + h // 2
        k = h // 2
        return psi(k) + k // 2 + 1

    for b, p in ((5, 13), (3, 7), (10, 7), (2, 7)):
        first = None
        for h in range(0, 5000):
            X = 1 + b * psi(h)
            y = 2 * h + 1
            if X % p == 0 and y % p == 0:
                first = (h, X, y)
                break
        print("F4 base %2d, p = %2d: first p-divisible state (h, X, y) = %s%s"
              % (b, p, first, "  balanced" if first and first[1] == first[2] else "  not balanced"))

    # F5
    def squares_counts(m, nsq):
        s1 = [0] * m
        for a in range(m):
            s1[(a * a) % m] += 1
        out = [0] * m
        out[0] = 1
        for _ in range(nsq):
            out = cconv(out, s1, m)
        return out

    m = 16
    four = squares_counts(m, 4)
    d2_lip = {c: Fraction(four[c] * m, m ** 4) for c in range(1, m, 2)}
    # Hurwitz norm form in the basis (i, j, k, w = (1+i+j+k)/2):  N = b^2 + c^2 + d^2 + e^2 + e (b + c + d)  with
    # element b i + c j + d k + e w -> coordinates (e/2, b + e/2, c + e/2, d + e/2)
    cnt = [0] * m
    for b_, c_, d_, e_ in itertools.product(range(m), repeat=4):
        cnt[(b_ * b_ + c_ * c_ + d_ * d_ + e_ * e_ + e_ * (b_ + c_ + d_)) % m] += 1
    d2_hur = {c: Fraction(cnt[c] * m, m ** 4) for c in range(1, m, 2)}
    print("F5 2-adic densities at odd classes modulo 16: Z^4 (Lipschitz) %s ; Hurwitz order %s"
          % (sorted(set(str(v) for v in d2_lip.values())), sorted(set(str(v) for v in d2_hur.values()))))
    print("   Siegel: Lipschitz pi^2 n * 1 * (8/pi^2) sigma(n)/n = 8 sigma(n); Hurwitz (covolume 1/2) 2 pi^2 n * (3/2) * "
          "(8/pi^2) sigma(n)/n = 24 sigma(n)")

    # F6
    bern = {2: Fraction(1, 6), 4: Fraction(-1, 30), 6: Fraction(1, 42), 8: Fraction(-1, 30)}
    for r in (2, 4, 6, 8):
        zeta_r = abs(bern[r]) * (2 * math.pi) ** r / (2 * math.factorial(r))
        lhs = (2 * math.pi) ** r / (math.factorial(r - 1) * zeta_r)
        print("F6 r = %d: (2 pi)^r / ((r-1)! zeta(r)) = %.6f ; 2r/|B_r| = %s" % (r, lhs, 2 * r / abs(bern[r])))


if __name__ == "__main__":
    main()
