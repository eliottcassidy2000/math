#!/usr/bin/env python3
"""Pisano periods, Wall-Sun-Sun primes, and the Fibonacci analogue of the Wieferich orbit explosion (THM-4644 part 2).

(1) Pisano periods pi(n) (period of Fibonacci mod n = order of Q = [[0,1],[1,1]] in GL_2(Z/n)) for n <= NMAX:
    max of pi(n)/n (classically <= 6, equality iff n = 2*5^k), the split/inert/ramified behaviour at primes
    (pi(p) | p - 1 for p = +-1 mod 5, pi(p) | 2(p + 1) but not p + 1 for p = +-2 mod 5, pi(5) = 20).
(2) Wall-Sun-Sun test p^2 | F_(p - (p/5)) for all primes p <= PMAX (fast doubling mod p^2).
(3) Fibonacci-Wieferich orbit count (split primes p = +-1 mod 5): phi = a root of x^2 - x - 1 in Z/p^2,
    d = ord_p(phi), d' = odd part, g = phi^(d/d').  The circulant tournaments on Z_(p^2) with gS = S number 2^E with
    E = (p - 1)/d' unless p is Wall-Sun-Sun (then E = (p + 1)(p - 1)/(2 d')), by the proof of THM-4644 part 2 with
    2 replaced by phi.  Checked by direct orbit enumeration for p <= ORBMAX (both roots).
(4) 'Fibonacci Paley primes': p = 3 mod 4, p = +-1 mod 5 with ord_p(phi) = (p - 1)/2 for the square root phi: the
    phi-invariant circulant tournaments on Z_p are exactly the Paley tournament and its reverse (p = 11 is the
    eleven-square torus group G_5 = <x + 1, phi x> = Aff(P_11)).
Session opus-2026-10-14-S36.
"""
import sys, math
from sympy import primerange, n_order, legendre_symbol

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
PMAX = int(float(sys.argv[2])) if len(sys.argv) > 2 else 2_000_000
ORBMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 200


def pisano(n):
    if n == 1:
        return 1
    a, b, k = 0, 1, 0
    while True:
        a, b = b, (a + b) % n
        k += 1
        if a == 0 and b == 1:
            return k


def fib_pair(k, m):
    """(F_k, F_{k+1}) mod m by fast doubling"""
    if k == 0:
        return 0, 1
    a, b = fib_pair(k >> 1, m)
    c = a * ((2 * b - a) % m) % m
    d = (a * a + b * b) % m
    return (d, (c + d) % m) if k & 1 else (c, d)


def hensel_phi(p):
    """the two roots of x^2 - x - 1 modulo p^2 (p = +-1 mod 5)"""
    roots = [x for x in range(p) if (x * x - x - 1) % p == 0]
    out = []
    for x in roots:
        # f(x + p t) = f(x) + p t f'(x) mod p^2, f'(x) = 2x - 1
        fx = x * x - x - 1
        t = (-(fx // p) * pow(2 * x - 1, -1, p)) % p
        y = x + p * t
        assert (y * y - y - 1) % (p * p) == 0
        out.append(y)
    return out


def orbit_count(n, g):
    """number of orbits of multiplication by g (a unit) on Z_n minus 0"""
    seen = [False] * n
    cnt = 0
    for x in range(1, n):
        if not seen[x]:
            cnt += 1
            y = x
            while not seen[y]:
                seen[y] = True
                y = y * g % n
    return cnt


def main():
    # (1) Pisano
    best = (0, None)
    viol = []
    for n in range(1, NMAX + 1):
        pn = pisano(n)
        if pn * 1.0 / n > best[0] + 1e-12:
            best = (pn / n, n)
        if pn == 6 * n and not (n % 2 == 0 and all(f in (2, 5) for f in [q for q in range(2, n + 1) if n % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))])):
            viol.append(n)
    six = [n for n in range(1, NMAX + 1) if pisano(n) == 6 * n]
    print(f'(1) n <= {NMAX}: max pi(n)/n = {best[0]} (first at n = {best[1]}); n with pi(n) = 6n: {six}')
    examples = {n: pisano(n) for n in (10, 11, 12, 19, 20, 50, 55, 89, 121, 125)}
    print('    examples pi(n):', examples)
    split_ok = inert_ok = True
    for p in primerange(7, NMAX + 1):
        pp = pisano(p)
        if p % 5 in (1, 4):
            split_ok &= ((p - 1) % pp == 0)
        else:
            inert_ok &= ((2 * (p + 1)) % pp == 0 and (p + 1) % pp != 0)
    print(f'    split p = +-1 mod 5: pi(p) | p - 1 at all p <= {NMAX}: {split_ok};  inert p = +-2 mod 5: pi(p) | 2(p+1), not | p+1: {inert_ok};  pi(5) = {pisano(5)}')
    # (2) Wall-Sun-Sun search
    wss = []
    cnt = 0
    for p in primerange(7, PMAX + 1):
        e = p - legendre_symbol(5, p)
        if fib_pair(e, p * p)[0] == 0:
            wss.append(p)
        cnt += 1
    print(f'(2) Wall-Sun-Sun primes (p^2 | F_(p - (p/5))) among the {cnt} primes 7 <= p <= {PMAX}: {wss}')
    # (3) orbit explosion, split primes
    bad = []
    rows = []
    for p in primerange(11, ORBMAX + 1):
        if p % 5 not in (1, 4):
            continue
        for phi in hensel_phi(p):
            d = n_order(phi % p, p)
            dp = d
            while dp % 2 == 0:
                dp //= 2
            g = pow(phi, d // dp, p * p)
            wss_p = pow(phi, p - 1, p * p) == 1
            N = orbit_count(p * p, g)
            assert N % 2 == 0
            E = N // 2
            pred = (p + 1) * (p - 1) // (2 * dp) if wss_p else (p - 1) // dp
            if E != pred:
                bad.append((p, phi, E, pred))
            rows.append((p, phi, d, dp, E))
    print(f'(3) split primes 11 <= p <= {ORBMAX}: orbit-count formula E = (p-1)/d\' (no Wall-Sun-Sun prime) holds for both roots at all '
          f'{len(rows)} (p, phi) pairs; failures: {bad}')
    print('    sample (p, phi mod p^2, ord_p phi, odd part, E):', rows[:6])
    # (4) Fibonacci Paley primes
    fpp = []
    for p in primerange(7, NMAX + 1):
        if p % 4 == 3 and p % 5 in (1, 4):
            for phi in hensel_phi(p):
                if legendre_symbol(phi % p, p) == 1 and n_order(phi % p, p) == (p - 1) // 2:
                    fpp.append(p)
    print(f'(4) primes p <= {NMAX}, p = 3 mod 4, p = +-1 mod 5, whose square golden root phi generates the squares '
          f'(phi-invariant circulants on Z_p = Paley and reverse): {fpp[:40]} ... ({len(fpp)} primes)')


if __name__ == '__main__':
    main()
