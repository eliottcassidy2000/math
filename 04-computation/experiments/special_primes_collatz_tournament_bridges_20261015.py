#!/usr/bin/env python3
"""Special prime families read through the repository's Collatz / tournament objects (session opus-2026-10-15-S37).

(1) Non-generous primes (OEIS A055578: the least primitive root g of p is not a primitive root of p^2).  Equivalent
    to g^(p-1) = 1 mod p^2: g is a Wieferich base at p.  THM-4644 part 2 / THM-4656 part 1 then say that the
    circulant tournaments on Z_(p^2) invariant under the odd part of g explode exactly at these p.  Check the known
    members 2, 40487, 6692367337 and that no other prime below LIMIT is non-generous.
(2) Wagstaff numbers W_p = (2^p + 1)/3 (p odd) are exactly the odd n with 3n - 1 a power of 2: the first odd
    ancestors of 1 under the 3x-1 map, as (4^k - 1)/3 are under 3x+1.  Check the 3x-1 itinerary.
(3) Woodall W_n = n 2^n - 1 and Cullen C_n = n 2^n + 1: under T(x) = x/2, (3x+1)/2, W_n makes n odd steps and lands at
    n 3^n - 1 (the 2-adic fixed point -1 of the odd branch), while C_n follows the trivial cycle word (OE)^j through
    3^j n 2^(n-2j) + 1 for 2j <= n - 1 (the 2-adic fixed point 1).  Check n <= 60.
(4) Leyland numbers x^y + y^x: for x = 2, 2^y + y^2 can be prime only if y = 1 or y = 3 mod 6.  List the primes with
    y < YMAX (probable-prime test), and the Leyland primes below LEY.
(5) Supersingular primes (Ogg: X_0(p)^+ has genus 0; the primes dividing |Monster|).  For p = 3 mod 4, p > 3, the
    number of F_p-rational supersingular j-invariants is (h(-p) + h(-4p))/2, and Dirichlet's class number formula makes
    h(-p) = I(p)/(2 - (2/p)), where I(p) = sum_(0<a<p/2) (a/p) is the half-interval score imbalance of vertex 0 in the
    Paley tournament P_p.  The total number of supersingular j-invariants is floor(p/12) + eps_p.  Check, for all
    p < PMAX, that 'every supersingular j is F_p-rational' happens exactly for the supersingular primes, and print the
    Paley half-imbalance reading for p = 3 mod 4.
(6) Stern primes (A042978): primes not of the form q + 2a^2 (q prime or q = 1 excluded; a >= 1).  List below SMAX.
"""
import sys, time
from math import gcd, isqrt
from sympy import isprime, primerange, primitive_root, factorint
from sympy.ntheory import legendre_symbol

LIMIT = 10**7
YMAX = 400
LEY = 10**12
PMAX = 2000
SMAX = 10**6


def least_primitive_root(p):
    if p == 2:
        return 1
    phi = p - 1
    qs = list(factorint(phi).keys())
    g = 2
    while True:
        if all(pow(g, phi // q, p) != 1 for q in qs):
            return g
        g += 1


def part1():
    t0 = time.time()
    found = []
    for p in primerange(3, LIMIT):
        g = least_primitive_root(p)
        if pow(g, p - 1, p * p) == 1:
            found.append((p, g))
    big = 6692367337
    gb = least_primitive_root(big)
    print(f'(1) non-generous primes 2 < p < {LIMIT}: {found}  [{time.time() - t0:.0f} s]')
    print(f'    p = 6692367337: least primitive root {gb}, g^(p-1) mod p^2 == 1: {pow(gb, big - 1, big * big) == 1}')
    print('    p = 2: least primitive root 1, which is not a primitive root of 4 (non-generous by convention)')
    return found


def part2():
    ok = True
    for p in range(3, 200, 2):
        W = (2 ** p + 1) // 3
        x = 3 * W - 1
        ok &= (x & (x - 1)) == 0 and x == 2 ** p
    # converse: odd n with 3n - 1 = 2^m forces m odd and n = (2^m + 1)/3
    conv = all(((2 ** m + 1) % 3 == 0) == (m % 2 == 1) for m in range(1, 200))
    wag = [p for p in primerange(3, 1300) if isprime((2 ** p + 1) // 3)]
    print(f'(2) 3*W_p - 1 = 2^p for odd p < 200: {ok}; 3 | 2^m + 1 iff m odd (m < 200): {conv}')
    print(f'    Wagstaff primes (2^p+1)/3, prime p < 1300: p = {wag}')


def T(x):
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2


def part3():
    okW = okC = True
    for n in range(2, 61):
        x = n * 2 ** n - 1
        for _ in range(n):
            okW &= (x % 2 == 1)
            x = T(x)
        okW &= (x == n * 3 ** n - 1)
        y = n * 2 ** n + 1
        j = 0
        while 2 * j <= n - 1:
            okC &= (y == 3 ** j * n * 2 ** (n - 2 * j) + 1)
            if 2 * j + 2 > n - 1:
                break
            okC &= (y % 2 == 1)
            y = T(y)
            okC &= (y % 2 == 0)
            y = T(y)
            j += 1
    print(f'(3) Woodall n 2^n - 1: n odd T-steps then n 3^n - 1 (2 <= n <= 60): {okW}; '
          f'Cullen n 2^n + 1 follows (OE)^j through 3^j n 2^(n-2j) + 1: {okC}')


def part4():
    ys = [y for y in range(1, YMAX) if isprime(2 ** y + y * y)]
    forced = all(y == 1 or y % 6 == 3 for y in ys)
    ley = set()
    x = 2
    while x ** 2 + 2 ** x <= LEY:
        y = 2
        while y <= x:
            v = x ** y + y ** x
            if v > LEY:
                break
            if isprime(v):
                ley.add((v, x, y))
            y += 1
        x += 1
    print(f'(4) 2^y + y^2 prime, y < {YMAX}: y = {ys}; all y = 1 or 3 mod 6: {forced}')
    print(f'    Leyland primes x^y + y^x below {LEY:.0e}: {sorted(ley)}')


def class_number(D):
    # number of reduced primitive positive definite forms of discriminant D < 0
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if b < 0 and a == c:
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    return h


def n_supersingular(p):
    r = {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]
    return p // 12 + r


SUPERSINGULAR = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]


def part5():
    allrat = []
    rows = []
    for p in primerange(5, PMAX):
        if p % 4 == 3:
            ratc = (class_number(-p) + class_number(-4 * p)) // 2
        else:
            ratc = class_number(-4 * p) // 2
        tot = n_supersingular(p)
        if ratc == tot:
            allrat.append(p)
        if p % 4 == 3 and p < 200:
            I = sum(legendre_symbol(a, p) for a in range(1, (p - 1) // 2 + 1))
            chi2 = legendre_symbol(2, p)
            h = class_number(-p)
            rows.append((p, I, h, I == h * (2 - chi2), ratc, tot))
    print(f'(5) primes 5 <= p < {PMAX} with every supersingular j-invariant F_p-rational: {allrat}')
    print(f'    equals the supersingular primes >= 5: {allrat == [q for q in SUPERSINGULAR if q >= 5]}')
    print('    p = 3 mod 4 (p < 200): (p, Paley half-imbalance I(p), h(-p), I = (2-(2/p)) h, #F_p-rational ss j, #ss j)')
    for r in rows:
        print('      ', r)


def part6():
    ps = list(primerange(2, SMAX))
    pset = set(ps)
    stern = []
    for p in ps:
        a = 1
        rep = False
        while 2 * a * a < p:
            if (p - 2 * a * a) in pset:
                rep = True
                break
            a += 1
        if not rep:
            stern.append(p)
    oddcomp = []
    for n in range(9, 20001, 2):
        if n in pset:
            continue
        a = 0
        rep = False
        while 2 * a * a < n:
            if (n - 2 * a * a) in pset:
                rep = True
                break
            a += 1
        if not rep:
            oddcomp.append(n)
    print(f'(6) Stern primes below {SMAX:.0e} (not q + 2a^2, q prime, a >= 1): {stern}')
    print(f'    odd composites below 20001 not of the form q + 2a^2 (a >= 0): {oddcomp}')


if __name__ == '__main__':
    part2()
    part3()
    part4()
    part5()
    part6()
    part1()
