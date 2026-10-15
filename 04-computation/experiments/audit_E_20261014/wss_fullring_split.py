#!/usr/bin/env python3
"""Audit E (2026-10-14), THM-4656: is the inert formula of part 2 also valid for SPLIT primes when the Cayley tournaments
live on the whole ring R_2 = Z[mu]/(p^2) = (Z/p^2)^2 (as the title says), rather than on one coordinate Z_(p^2)
(as part 1 does)?  Prediction (audit): for every prime p not dividing 2D, split or inert,
    #<g>-orbits on R_2 \\ {0} = (p^2 - 1)(p^(2 - eps) + 1)/pi',   eps = 0 iff p is a-Wall-Sun-Sun,
so E = (p+1)(p^2-1)/(2 pi') normally and (p^2+1)(p^2-1)/(2 pi') at a-WSS primes -- the same as for inert primes.
Method: Burnside on R_2 = Z/p^2 x Z/p^2 via mu -> (mu_1, mu_2) (Hensel lifts of the two roots), g = mu^(pi/pi'),
#Fix(h) = p^min(v(h_1 - 1), 2) * p^min(v(h_2 - 1), 2) - 1; negation-fixed orbits via Burnside on <g> x {+-1}."""
from sympy import primerange

def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n

def vp(x, p, cap=2):
    x %= p ** cap
    if x == 0:
        return cap
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v

def lift(r, a, p):
    m = p * p
    for _ in range(4):
        f = (r * r - a * r - 1) % m
        r = (r - f * pow((2 * r - a) % m, -1, m)) % m
    assert (r * r - a * r - 1) % m == 0
    return r

def order_mod(x, m, bound):
    k, y = 1, x % m
    while y != 1:
        y = y * x % m; k += 1
        if k > bound:
            raise ValueError
    return k

def check(a, p):
    m = p * p
    roots = [x for x in range(p) if (x * x - a * x - 1) % p == 0] if p < 10 ** 5 else None
    if roots is None:
        # sympy-free root finding for large p: Tonelli-Shanks on D = a^2 + 4
        D = (a * a + 4) % p
        from sympy.ntheory.residue_ntheory import sqrt_mod
        s = sqrt_mod(D, p)
        inv2 = pow(2, -1, p)
        roots = [((a + s) * inv2) % p, ((a - s) * inv2) % p]
    assert len(roots) == 2
    m1, m2 = lift(roots[0], a, p), lift(roots[1], a, p)
    d1 = order_mod(m1, p, p); d2 = order_mod(m2, p, p)
    from math import lcm
    pi = lcm(d1, d2)
    pi2 = lcm(order_mod(m1, m, p * p), order_mod(m2, m, p * p))
    wss = (pi2 == pi)
    e = pi // odd_part(pi)
    g1, g2 = pow(m1, e, m), pow(m2, e, m)
    # enumerate <g>
    G = []; h1, h2 = 1, 1
    while True:
        G.append((h1, h2))
        h1, h2 = h1 * g1 % m, h2 * g2 % m
        if (h1, h2) == (1, 1):
            break
    s1 = sum(p ** vp(x - 1, p) * p ** vp(y - 1, p) - 1 for x, y in G)
    s2 = sum(p ** vp(-x - 1, p) * p ** vp(-y - 1, p) - 1 for x, y in G)
    N = s1 // len(G); assert N * len(G) == s1
    Npm = (s1 + s2) // (2 * len(G)); assert Npm * 2 * len(G) == s1 + s2
    po = odd_part(pi)
    eps = 0 if wss else 1
    pred = (p * p - 1) * (p ** (2 - eps) + 1) // po
    return pi, wss, N, 2 * Npm - N, pred

def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)

ok = bad = 0
for a in (1, 2, 3, 6):
    D = a * a + 4
    for p in primerange(3, 200):
        if (2 * D) % p == 0 or legendre(D, p) != 1:
            continue
        pi, wss, N, nf, pred = check(a, p)
        good = (N == pred and nf == 0)
        ok += good; bad += (not good)
        if wss or not good:
            print(f"  a={a} p={p} split pi={pi} {'WSS' if wss else ''} full-ring orbits {N} predicted {pred} negation-fixed {nf} {'ok' if good else 'MISMATCH'}")
print(f"split primes < 200 (a = 1, 2, 3, 6): {ok} agree with the inert formula on the full ring, {bad} mismatches")
for a, p in [(2, 1546463), (6, 643)]:
    pi, wss, N, nf, pred = check(a, p)
    print(f"  a={a} p={p} split pi={pi} WSS={wss}: full-ring orbits {N} predicted {pred} (non-WSS value {(p*p-1)*(p+1)//odd_part(pi)}) negation-fixed {nf}")
