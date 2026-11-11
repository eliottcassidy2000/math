#!/usr/bin/env python3
"""Auditor bridges-note, script A4 (independent code): remaining checks for sections 3-9.

(1) class counts of a^2+b^2+c^2+d^2 against a1 a2 + u b1 b2 modulo 9, 27, 81, 25, 125, 49, 121 (equal) and 4, 8 (differ)
(2) (N_m * N_m)(c) = m^3 L_m(c) for odd m (own formula for L_m)
(3) base-b balanced states X = 1 + b psi(h) = 2h + 1; early p-divisible states in bases 3, 5, 10
(4) Pi_2(C_3)/H(C_3)
(5) every labelled 8-tournament with the fixed-point-free anti-involution (01)(23)(45)(67): H mod 4, H mod 8
(6) head-series values at 7 and 21 (THM-5134 (E), quoted in section 9)
"""
import itertools
from math import gcd, pi
import numpy as np
from numba import njit, prange


def conv(a, b, q):
    out = np.zeros(q, dtype=object)
    for x in range(q):
        if a[x] == 0:
            continue
        for y in range(q):
            out[(x + y) % q] += a[x] * b[y]
    return out


def part1():
    print("=== A4.1 four squares against the seam a1 a2 + u b1 b2 ===")
    for q in (9, 27, 81, 25, 125, 49, 121, 4, 8, 16):
        sq = np.zeros(q, dtype=object)
        for a in range(q):
            sq[(a * a) % q] += 1
        four = conv(conv(sq, sq, q), conv(sq, sq, q), q)
        res = []
        for u in (1, 2, 4, 3, 5):
            if gcd(u, q) != 1 and q % 2 == 1:
                continue
            pr = np.zeros(q, dtype=object)
            pu = np.zeros(q, dtype=object)
            for a in range(q):
                for b in range(q):
                    pr[(a * b) % q] += 1
                    pu[(u * a * b) % q] += 1
            seam = conv(pr, pu, q)
            res.append((u, list(four) == list(seam)))
        print("  q = %3d: four squares = seam for u in %s" % (q, res))


def L(m, e):
    tot = 0.0
    from fractions import Fraction
    tot = Fraction(0)
    g = gcd(m, e % m) if e % m else m
    for d in range(1, m + 1):
        if g % d == 0 and m % d == 0:
            term = Fraction(1, d)
            md = m // d
            p = 2
            t = md
            primes = set()
            while p * p <= t:
                while t % p == 0:
                    primes.add(p)
                    t //= p
                p += 1
            if t > 1:
                primes.add(t)
            for p in primes:
                term *= (1 - Fraction(1, p * p))
            tot += term
    return tot


def part2():
    print()
    print("=== A4.2 additive square of the multiplication table = m^3 L_m ===")
    from fractions import Fraction
    for m in (3, 5, 9, 15, 21, 25, 27, 35, 45, 49, 63, 75, 81, 105):
        tab = np.zeros(m, dtype=object)
        for a in range(m):
            for b in range(m):
                tab[(a * b) % m] += 1
        sq = conv(tab, tab, m)
        ok = all(Fraction(int(sq[c])) == m ** 3 * L(m, c) for c in range(m))
        print("  m = %3d: %s" % (m, ok))


def part3():
    print()
    print("=== A4.3 base-b balanced and early divisible states ===")
    LIM = 1 << 22
    psi = np.zeros(LIM, np.int64)
    for h in range(1, LIM):
        k = h >> 1
        psi[h] = psi[k] + k if h % 2 == 0 else psi[k] + k // 2 + 1
    hs = np.arange(LIM)
    for b in (2, 3, 4, 5, 6, 10):
        bal = hs[(1 + b * psi) == (2 * hs + 1)]
        print("  b = %2d: balanced states below 2^22: %s" % (b, list(map(int, bal[:14]))))
    print("  psi(h) > h/2 for 1 <= h < 2^22:", bool(np.all(2 * psi[1:] > hs[1:])))
    for b, p in ((3, 7), (10, 7), (5, 13), (2, 7), (4, 7), (5, 7), (6, 7)):
        X = 1 + b * psi
        y = 2 * hs + 1
        div = hs[(X % p == 0) & (y % p == 0)]
        h0 = int(div[0])
        print("  base %2d, p = %2d: first p-divisible state h = %d (X = %d, y = %d, balanced: %s); next: %s"
              % (b, p, h0, 1 + b * int(psi[h0]), 2 * h0 + 1, 1 + b * int(psi[h0]) == 2 * h0 + 1, list(map(int, div[1:4]))))


@njit(cache=True)
def hcount(adj, n):
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            c = dp[mask, v]
            if c == 0:
                continue
            out = adj[v] & ~mask
            w = 0
            while out:
                if out & 1:
                    dp[mask | (1 << w), w] += c
                out >>= 1
                w += 1
    s = 0
    for v in range(n):
        s += dp[full, v]
    return s


def part4():
    print()
    print("=== A4.4 Pi_2(C_3)/H(C_3) ===")
    arcs = {(0, 1), (1, 2), (2, 0)}

    def H(vs):
        if len(vs) == 0:
            return 1
        return sum(1 for p in itertools.permutations(vs) if all((p[i], p[i + 1]) in arcs for i in range(len(p) - 1)))
    pi2 = sum(H([v for v in range(3) if (m >> v) & 1]) * H([v for v in range(3) if not (m >> v) & 1]) for m in range(8))
    print("  Pi_2(C_3) = %d, H(C_3) = %d, a = %s ; sigma(3)/3 = 4/3" % (pi2, H([0, 1, 2]), pi2 / H([0, 1, 2])))


@njit(parallel=True, cache=True)
def sc8(orb_a, orb_b, orb_c, orb_d, npair, nfree):
    # free bits: 4 for the pairs {i, sigma i}, 12 for the orbits of other pairs
    total = 1 << (npair + nfree)
    out = np.zeros(total, np.int64)
    for code in prange(total):
        adj = np.zeros(8, np.int64)
        for k in range(npair):
            i = 2 * k
            if (code >> k) & 1:
                adj[i] |= (1 << (i + 1))
            else:
                adj[i + 1] |= (1 << i)
        for k in range(nfree):
            u, v, su, sv = orb_a[k], orb_b[k], orb_c[k], orb_d[k]
            if (code >> (npair + k)) & 1:
                adj[u] |= (1 << v)        # u -> v  forces  sigma(v) -> sigma(u)
                adj[sv] |= (1 << su)
            else:
                adj[v] |= (1 << u)        # v -> u  forces  sigma(u) -> sigma(v)
                adj[su] |= (1 << sv)
        out[code] = hcount(adj, 8)
    return out


def part5():
    print()
    print("=== A4.5 self-converse 8-tournaments (anti-involution (01)(23)(45)(67)) ===")
    sig = [1, 0, 3, 2, 5, 4, 7, 6]
    seen = set()
    oa, ob, oc, od = [], [], [], []
    for u in range(8):
        for v in range(u + 1, 8):
            if sig[u] == v:
                continue
            key = frozenset((frozenset((u, v)), frozenset((sig[u], sig[v]))))
            if key in seen:
                continue
            seen.add(key)
            oa.append(u)
            ob.append(v)
            oc.append(sig[u])
            od.append(sig[v])
    print("  orbits of non-fixed pairs:", len(oa))
    Hs = sc8(np.array(oa), np.array(ob), np.array(oc), np.array(od), 4, len(oa))
    print("  labelled tournaments with this anti-involution: %d ; all H = 1 mod 4: %s ; residues mod 8: %s"
          % (len(Hs), bool(np.all(Hs % 4 == 1)), dict(zip(*map(lambda a: list(map(int, a)), np.unique(Hs % 8, return_counts=True))))))
    print("  distinct H values: %d, min %d, max %d" % (len(set(map(int, Hs))), int(Hs.min()), int(Hs.max())))


def part6():
    print()
    print("=== A4.6 head and tail series at 7 and 21 ===")
    def sig(n):
        return sum(d for d in range(1, n + 1) if n % d == 0)

    def odd(n):
        while n % 2 == 0:
            n //= 2
        return n
    w0 = 854415478 / (854415478 + 1083261684)
    w1 = 1 - w0
    for N in (7, 21, 15, 31, 63):
        n0, n1 = odd(N + 1), odd(N)
        S0 = 8 / pi ** 2 * sig(n0) / n0
        S1 = 8 / pi ** 2 * sig(n1) / n1
        print("  N = %2d: S_0 = %.4f (n = %d), S_1 = %.4f (n = %d), w0 S0 + w1 S1 = %.4f (triple shares at 2^28: %.4f, %.4f)"
              % (N, S0, n0, S1, n1, w0 * S0 + w1 * S1, w0, w1))
    print("  8/pi^2 = %.4f : S_0 is minimal exactly when N + 1 is a power of 2 (N = 1, 3, 7, 15, 31, ...)" % (8 / pi ** 2))


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()
