#!/usr/bin/env python3
"""Audit E (2026-10-14), THM-4656 parts 1, 2, 4 and the a-Wall-Sun-Sun lists, by methods different from wss_orbits.py.

(A) pi_a(p), pi_a(p^2) by direct iteration of the pair (U_n, U_(n+1)) of U_(n+1) = a U_n + U_(n-1) until it returns
    to (0, 1) (no divisor search, no matrix powers).
(B) Number of <g>-orbits on R_2 \\ {0} by BURNSIDE (average number of fixed points), not by orbit tracing:
      inert p:  R_2 = Z[mu]/(p^2) as pairs x + y mu;  #Fix(h) = #{r : (h - 1) r = 0} - 1, counted by brute force
                for p <= 13 and by the valuation rule p^(2 v(h-1)) - 1 (cross-checked) for larger p;
      split p:  Z_(p^2) with the Hensel lift mu_1 found by brute-force search over x = r + t p;
                #Fix(h) = p^min(v_p(h-1), 2) - 1.
    The number of negation-fixed orbits is 2 N(<g> x {+-1}) - N(<g>) (Burnside again).
    Formulas (THM-4656):  inert  N = (p+1)(p^2-1)/pi'  or  (p^2+1)(p^2-1)/pi' at a-WSS primes;
                          split  N = 2(p-1)/d'         or  (p+1)(p-1)/d'      at a-WSS primes.
    Includes the Pell-WSS primes 13, 31 AND 1546463 (a = 2) and the bronze 241 (a = 3).
(C) a-WSS primes by the criterion p^2 | U_(p - (D/p))(a, -1) (rank of apparition), a = 1..7.
(D) golden Paley primes p < 2000: p = 3 mod 4, p = +-1 mod 5, the residue root of x^2 - x - 1 has order (p-1)/2.
(E) small facts: Q^5 = 3I mod 5; ladder resistances; triangular Fibonacci numbers; cannonball sums."""
import sys
from fractions import Fraction as F
from math import gcd, isqrt
from sympy import primerange, isprime, factorint

def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def period_iter(a, m):
    """least k > 0 with (U_k, U_(k+1)) = (0, 1) mod m."""
    x, y, k = 0, 1 % m, 0
    while True:
        x, y = y, (a * y + x) % m
        k += 1
        if x == 0 and y == 1 % m:
            return k

def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n

def vp(x, p, cap):
    x %= p ** cap
    if x == 0:
        return cap
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v

# ---------------- (B) Burnside ----------------
def burnside_inert(a, p, pi, brute):
    m = p * p
    def mul(u, v):
        x1, y1 = u; x2, y2 = v
        return ((x1 * x2 + y1 * y2) % m, (x1 * y2 + y1 * x2 + a * y1 * y2) % m)
    # g = mu^(pi/pi')
    k = pi // odd_part(pi)
    g = (1, 0)
    for _ in range(k):
        g = mul(g, (0, 1))
    G = [(1, 0)]
    while True:
        nxt = mul(G[-1], g)
        if nxt == (1, 0):
            break
        G.append(nxt)
    def fix_count(c):            # #{r != 0 : c r = 0}
        if brute:
            cnt = 0
            for x in range(m):
                for y in range(m):
                    if mul(c, (x, y)) == (0, 0):
                        cnt += 1
            return cnt - 1
        v = min(vp(c[0], p, 2), vp(c[1], p, 2))
        return p ** (2 * v) - 1
    s1 = s2 = 0
    for h in G:
        s1 += fix_count(((h[0] - 1) % m, h[1]))
        s2 += fix_count(((-h[0] - 1) % m, (-h[1]) % m))
    N = s1 // len(G); assert N * len(G) == s1
    Npm = (s1 + s2) // (2 * len(G)); assert Npm * 2 * len(G) == s1 + s2
    return N, 2 * Npm - N, len(G)

def burnside_split(a, p):
    m = p * p
    r = next(x for x in range(p) if (x * x - a * x - 1) % p == 0)
    lifts = [x for x in range(r, m, p) if (x * x - a * x - 1) % m == 0]
    assert len(lifts) == 1
    mu1 = lifts[0]
    d = next(t for t in range(1, p) if pow(mu1, t, p) == 1)
    g = pow(mu1, d // odd_part(d), m)
    G = [1]
    while True:
        nxt = G[-1] * g % m
        if nxt == 1:
            break
        G.append(nxt)
    s1 = sum(p ** vp(h - 1, p, 2) - 1 for h in G)
    s2 = sum(p ** vp(-h - 1, p, 2) - 1 for h in G)
    N = s1 // len(G); assert N * len(G) == s1
    Npm = (s1 + s2) // (2 * len(G)); assert Npm * 2 * len(G) == s1 + s2
    return N, 2 * Npm - N, len(G), d, mu1

print("(A)+(B) pi by iteration, WSS by pi(p^2) = pi(p), orbit counts by Burnside vs THM-4656 formulas")
cases = [(1, p) for p in primerange(3, 60)] + [(2, p) for p in primerange(3, 60)] + [(3, p) for p in primerange(3, 45)] \
      + [(a, p) for a in (4, 5, 6, 7) for p in primerange(3, 30)]
nok = nbad = 0
for a, p in cases:
    D = a * a + 4
    if (2 * D) % p == 0:
        continue
    pi = period_iter(a, p)
    pi2 = period_iter(a, p * p)
    assert pi2 in (pi, p * pi)
    wss = (pi2 == pi)
    ls = legendre(D, p)
    if ls == -1:
        N, negfix, Gord = burnside_inert(a, p, pi, brute=(p <= 7))
        if p <= 13 and p > 7:
            Nb, nfb, _ = burnside_inert(a, p, pi, brute=True)
            assert (Nb, nfb) == (N, negfix)
        po = odd_part(pi)
        pred = ((p * p + 1) if wss else (p + 1)) * (p * p - 1) // po
        assert pi % 2 == 0 and (2 * (p + 1)) % pi == 0
    else:
        N, negfix, Gord, d, mu1 = burnside_split(a, p)
        dp = odd_part(d)
        assert odd_part(pi) == dp and (p - 1) % pi == 0
        pred = ((p + 1) * (p - 1) // dp) if wss else (2 * (p - 1) // dp)
    ok = (N == pred and negfix == 0)
    nok += ok; nbad += (not ok)
    if wss or not ok or (a, p) in ((1, 13),):
        print(f"  a={a} p={p} {'inert' if ls == -1 else 'split'} pi={pi} pi(p^2)={pi2} {'WSS' if wss else ''} |<g>|={Gord} "
              f"Burnside orbits {N} formula {pred} negation-fixed {negfix} {'ok' if ok else 'MISMATCH'}")
print(f"  {nok} cases agree, {nbad} mismatches")

# Larger a-WSS primes: Burnside via valuations (no brute force), pi via the order of mu computed by factoring
def mat_mul(A, B, m):
    return ((A[0]*B[0] + A[1]*B[2]) % m, (A[0]*B[1] + A[1]*B[3]) % m, (A[2]*B[0] + A[3]*B[2]) % m, (A[2]*B[1] + A[3]*B[3]) % m)
def mat_pow(a, e, m):
    R = (1, 0, 0, 1); B = (a % m, 1, 1, 0)
    while e:
        if e & 1:
            R = mat_mul(R, B, m)
        B = mat_mul(B, B, m); e >>= 1
    return R
def order_mu(a, p, mod):
    D = a * a + 4
    ls = legendre(D, p)
    M = (p - 1) if ls == 1 else 2 * (p + 1)
    if mod == p * p:
        M *= p
    o = M
    for q in factorint(M):
        while o % q == 0 and mat_pow(a, o // q, mod) == (1, 0, 0, 1):
            o //= q
    return o, ls

for a, p in [(2, 1546463), (3, 241), (6, 191), (6, 643)]:
    pi, ls = order_mu(a, p, p)
    pi2, _ = order_mu(a, p, p * p)
    wss = pi2 == pi
    if ls == -1:
        N, negfix, Gord = burnside_inert(a, p, pi, brute=False)
        pred = ((p * p + 1) if wss else (p + 1)) * (p * p - 1) // odd_part(pi)
        normal = (p + 1) * (p * p - 1) // odd_part(pi)
    else:
        N, negfix, Gord, d, mu1 = burnside_split(a, p)
        pred = ((p + 1) * (p - 1) // odd_part(d)) if wss else (2 * (p - 1) // odd_part(d))
        normal = 2 * (p - 1) // odd_part(d)
    print(f"  a={a} p={p} {'inert' if ls == -1 else 'split'} pi={pi} WSS={wss} |<g>|={Gord}: Burnside orbits {N}, formula {pred}, "
          f"non-WSS value would be {normal}; negation-fixed {negfix}; {'ok' if N == pred and negfix == 0 else 'MISMATCH'}")

# ---------------- (C) a-WSS primes by rank of apparition ----------------
print("(C) a-WSS primes via p^2 | U_(p-(D/p))(a,-1)")
def U_mod(a, n, m):
    # U_n via matrix power: [[a,1],[1,0]]^n = [[U_(n+1), U_n],[U_n, U_(n-1)]]
    return mat_pow(a, n, m)[1]
known = {}
for a, B in [(1, 300000), (2, 300000), (3, 300000), (4, 100000), (5, 100000), (6, 100000), (7, 100000)]:
    D = a * a + 4
    found = []
    for p in primerange(3, B):
        if (2 * D) % p == 0:
            continue
        e = p - legendre(D, p)
        if U_mod(a, e, p * p) == 0:
            found.append(p)
    known[a] = found
    print(f"  a={a} D={D}: {found} below {B}")
for a, p in [(2, 1546463), (5, 21011), (6, 134339)]:
    D = a * a + 4
    e = p - legendre(D, p)
    print(f"  check a={a} p={p}: p^2 | U_(p-(D/p)) ? {U_mod(a, e, p * p) == 0}")

# ---------------- (D) golden Paley primes ----------------
gp_ = []
for p in primerange(3, 2000):
    if p % 4 != 3 or p % 5 not in (1, 4):
        continue
    roots = [x for x in range(p) if (x * x - x - 1) % p == 0]
    qr = [r for r in roots if legendre(r, p) == 1]
    assert len(roots) == 2 and len(qr) == 1
    r = qr[0]
    o = next(t for t in range(1, p) if pow(r, t, p) == 1)
    if o == (p - 1) // 2:
        gp_.append(p)
print(f"(D) golden Paley primes < 2000: {len(gp_)}: {gp_[:12]} ...")
# p = 11 details
r11 = [x for x in range(11) if (x * x - x - 1) % 11 == 0 and legendre(x, 11) == 1][0]
print(f"    p=11: residue root {r11}, <{r11}> = {sorted({pow(r11, k, 11) for k in range(10)})}, QR_11 = {sorted({x*x % 11 for x in range(1, 11)})}")

# ---------------- (E) small facts ----------------
def mm(A, B, m):
    return [[sum(A[i][k] * B[k][j] for k in range(2)) % m for j in range(2)] for i in range(2)]
Q = [[1, 1], [1, 0]]
P = [[1, 0], [0, 1]]
sc = []
for k in range(1, 31):
    P = mm(P, Q, 5)
    if P[0][1] == 0 and P[1][0] == 0 and P[0][0] == P[1][1]:
        sc.append((k, P[0][0]))
print("(E) Q^k mod 5 scalar for (k, scalar):", sc, "; -1 mod 5 =", (-1) % 5, "is a square:", any(x * x % 5 == 4 for x in range(5)))
# ladder: n sections, R_1 = 2, R_(k+1) = 1 + (1 || R_k)
fib = [0, 1]
for _ in range(200):
    fib.append(fib[-1] + fib[-2])
R = F(2); okl = True
for nn in range(1, 40):
    if R != F(fib[2 * nn + 1], fib[2 * nn]):
        okl = False
    R = 1 + R / (1 + R)
print("    ladder R_n = F_(2n+1)/F_(2n) for n = 1..39:", okl)
tri = [f for f in fib[1:2000] if f >= 1 and isqrt(8 * f + 1) ** 2 == 8 * f + 1]
print("    triangular Fibonacci numbers among F_1..F_200:", sorted(set(tri)))
sq = [k for k in range(1, 2 * 10 ** 6) if isqrt(k * (k + 1) * (2 * k + 1) // 6) ** 2 == k * (k + 1) * (2 * k + 1) // 6]
print("    n <= 2*10^6 with 1^2+...+n^2 a square:", sq, "; 1^2+..+24^2 =", 24 * 25 * 49 // 6, "= 70^2 =", 70 ** 2)
