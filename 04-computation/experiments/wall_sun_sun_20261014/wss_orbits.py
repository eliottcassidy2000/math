#!/usr/bin/env python3
"""THM-4656 checks: metallic Wieferich primes and orbit explosion of unit-invariant Cayley tournaments on Z[mu]/(p^2).

mu = mu_a is a root of x^2 - a x - 1 (a = 1 golden, a = 2 silver/Pell, a = 3 bronze, ...), D = a^2 + 4.  For a prime p not
dividing 2D, pi_a(p) = order of mu in (Z[mu]/p)^* (= Pisano-type period of U_n(a,-1)), and pi_a(p^2) in {pi, p pi}.
p is an a-Wall-Sun-Sun prime iff pi_a(p^2) = pi_a(p).  Let pi' be the odd part of pi_a(p), g = mu^(pi/pi').
Cayley tournaments on (R, +), R = Z[mu]/(p^2), with g S = S: exponent E = (#<g>-orbits on R \ {0})/2, and
  inert p ((D/p) = -1): E = (p+1)(p^2-1)/(2 pi')  normally,  (p^2+1)(p^2-1)/(2 pi')  at a-WSS primes;
  split p ((D/p) = +1): on the cyclic group Z_(p^2) with base mu_1 (a root mod p^2), THM-4644's formula
      E = (p-1)/d'  normally,  (p+1)(p-1)/(2 d')  at a-WSS primes   (d = ord_p mu_1, d' odd part).
(1) brute-force orbit counts vs formulas for small p, several a, including a-WSS primes (where the jump must appear);
(2) a-WSS primes below a bound; (3) golden 'Fermat' primes (pi(p) a power of 2) and the Lucas numbers L_(2^j);
(4) golden Paley primes: p = 3 mod 4 split with <mu_1> = QR_p (golden-invariant tournaments on Z_p = Paley and reverse)."""
import sys
from sympy import primerange, factorint, legendre_symbol, isprime

def mat_mult(A, B, m):
    return ((A[0]*B[0] + A[1]*B[2]) % m, (A[0]*B[1] + A[1]*B[3]) % m, (A[2]*B[0] + A[3]*B[2]) % m, (A[2]*B[1] + A[3]*B[3]) % m)

def mat_pow(a, e, m):
    R = (1, 0, 0, 1); B = (a % m, 1, 1, 0)
    while e:
        if e & 1: R = mat_mult(R, B, m)
        B = mat_mult(B, B, m); e >>= 1
    return R

def period(a, m, cands):
    """least k in candidate multiples with Q_a^k = I mod m (Q_a = [[a,1],[1,0]] acts as multiplication by mu)."""
    for k in sorted(cands):
        if mat_pow(a, k, m) == (1, 0, 0, 1): return k
    return None

def divisors_of(n):
    ds = [1]
    for q, e in factorint(n).items():
        ds = [d * q ** k for d in ds for k in range(e + 1)]
    return ds

def pi_p(a, p):
    D = a * a + 4
    ls = legendre_symbol(D % p, p) if D % p else 0
    bound = p - 1 if ls == 1 else 2 * (p + 1)
    return period(a, p, divisors_of(bound)), ls

def is_wss(a, p, pi):
    return mat_pow(a, pi, p * p) == (1, 0, 0, 1)

def odd_part(n):
    while n % 2 == 0: n //= 2
    return n

def orbit_exponent_inert(a, p, pi):
    """brute force: R = Z[mu]/(p^2) as pairs (x, y) = x + y mu, mu^2 = a mu + 1; orbits of <g> on R \\ {0}."""
    m = p * p; po = odd_part(pi); k = pi // po
    # g = mu^k as pair: compute via matrix: mu^k = U_(k) mu + U_(k-1)  (multiplication matrix on (1, mu))
    def mul(u, v):      # (x1 + y1 mu)(x2 + y2 mu) = x1x2 + y1y2 + (x1y2 + y1x2 + a y1 y2) mu
        x1, y1 = u; x2, y2 = v
        return ((x1 * x2 + y1 * y2) % m, (x1 * y2 + y1 * x2 + a * y1 * y2) % m)
    g = (1, 0); base = (0, 1); e = k
    while e:
        if e & 1: g = mul(g, base)
        base = mul(base, base); e >>= 1
    seen = set(); orbits = 0; neg_fixed = 0
    for x in range(m):
        for y in range(m):
            if (x, y) == (0, 0) or (x, y) in seen: continue
            orb = []
            u = (x, y)
            while u not in seen:
                seen.add(u); orb.append(u); u = mul(u, g)
            orbits += 1
            if ((-x) % m, (-y) % m) in set(orb): neg_fixed += 1
    return orbits, neg_fixed

def orbit_exponent_split(a, p):
    """Z_(p^2) with base mu_1 (root of x^2 - a x - 1 mod p^2, Hensel); orbits of <g1> on Z_(p^2) \\ {0}, g1 = mu_1^(d/d')."""
    m = p * p
    r = next(x for x in range(p) if (x * x - a * x - 1) % p == 0)
    for _ in range(3):                                   # Hensel lift
        f = (r * r - a * r - 1) % m; fp = (2 * r - a) % m
        r = (r - f * pow(fp, -1, m)) % m
    assert (r * r - a * r - 1) % m == 0
    d = next(t for t in sorted(divisors_of(p - 1)) if pow(r, t, p) == 1)
    dp = odd_part(d); g = pow(r, d // dp, m)
    seen = set(); orbits = 0; neg_fixed = 0
    for x in range(1, m):
        if x in seen: continue
        orb = []; u = x
        while u not in seen:
            seen.add(u); orb.append(u); u = u * g % m
        orbits += 1
        if (-x) % m in orb: neg_fixed += 1
    wss = pow(r, d, m) == 1
    return orbits, neg_fixed, d, dp, wss

def main():
    out = []
    print("(1) brute-force orbit counts vs formulas")
    nchk = 0
    for a, plist in [(1, [3, 7, 11, 13, 17, 19, 23, 29, 31, 37]), (2, [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]), (3, [5, 7, 11, 17, 19, 23, 29])]:
        for p in plist:
            D = a * a + 4
            if (2 * D) % p == 0: continue
            pi, ls = pi_p(a, p)
            wss = is_wss(a, p, pi)
            po = odd_part(pi)
            if ls == -1:
                orbits, negfix = orbit_exponent_inert(a, p, pi)
                pred = ((p * p + 1) if wss else (p + 1)) * (p * p - 1) // po
                tag = 'inert'
            else:
                orbits, negfix, d, dp, wss2 = orbit_exponent_split(a, p)
                assert wss2 == wss, (a, p)
                pred = ((p + 1) * (p - 1) // dp) if wss else (2 * (p - 1) // dp)
                tag = 'split'
            ok = (orbits == pred and negfix == 0)
            nchk += 1
            print(f"  a={a} p={p:3d} {tag:5s} pi={pi:4d} odd part {po:4d} {'WSS ' if wss else '    '} orbits {orbits:7d} predicted {pred:7d} "
                  f"negation-fixed {negfix}  E = {orbits // 2}  {'ok' if ok else 'MISMATCH'}")
            assert ok
    print(f"  {nchk} cases agree (including the a-WSS primes, where the count jumps)")

    print("(2) a-Wall-Sun-Sun primes below bounds")
    for a, B in [(1, 300000), (2, 300000), (3, 300000), (4, 100000), (5, 100000), (6, 100000), (7, 100000)]:
        D = a * a + 4
        found = [p for p in primerange(3, B) if (2 * D) % p and is_wss(a, p, pi_p(a, p)[0])]
        print(f"  a={a} (D={D}): {found} below {B}")

    print("(3) golden 'Fermat' primes: pi(p) a power of 2, p < 10^6, and Lucas numbers L_(2^j)")
    gf = []
    for p in primerange(3, 10 ** 6):
        if p == 5: continue
        pi, _ = pi_p(1, p)
        if pi & (pi - 1) == 0: gf.append((p, pi))
    def lucas_pow2_mod(j, p):                            # L_(2^j) mod p for j >= 1: L_2 = 3, L_(2n) = L_n^2 - 2 (n even)
        x = 3 % p
        for _ in range(j - 1): x = (x * x - 2) % p
        return x
    # pi(p) = 2^k  <=>  p | L_(2^(k-2))   (audit E, W7: the earlier factor-list check stopped at j = 6 and was wrong)
    bad = [p for p, pi in gf if pi < 8 or lucas_pow2_mod(pi.bit_length() - 3, p) != 0]
    print("  primes with pi(p) = 2^k:", gf)
    print("  every such p divides L_(2^(k-2)) (direct check mod p); mismatches:", bad)

    print("(4) golden Paley primes: p = 3 mod 4, p = +-1 mod 5, <mu_1> = QR_p (index 2, odd order), p < 2000")
    gp = []
    for p in primerange(7, 2000):
        if p % 4 != 3 or p % 5 not in (1, 4): continue
        r = next(x for x in range(p) if (x * x - x - 1) % p == 0)
        d = next(t for t in sorted(divisors_of(p - 1)) if pow(r, t, p) == 1)
        r2 = (1 - r) % p
        d2 = next(t for t in sorted(divisors_of(p - 1)) if pow(r2, t, p) == 1)
        if d == (p - 1) // 2 or d2 == (p - 1) // 2: gp.append(p)
    print("  ", gp)


if __name__ == '__main__':
    main()
