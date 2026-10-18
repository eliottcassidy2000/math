#!/usr/bin/env python3
"""Fibonacci started at (5,1) and (7,1) instead of (1,1): exact facts, checked (mac-mini-2026-10-18).

G_0 = a, G_1 = b, G_(n+1) = G_n + G_(n-1).  Notation: N(a,b) = b^2 - a b - a^2 (Cassini / golden norm form; repo's
O4 convention 'norm b^2 + ab - a^2' differs by the orientation of the pair, the class up to sign is the same).
Checks printed:
  1. first terms; closed forms G_n = a F_(n-1) + b F_n; Binet constants; Cassini G_n^2 - G_(n-1) G_(n+1) = (-1)^(n+1) N.
  2. Lucas/Fibonacci decompositions found by solving 2x2 systems.
  3. primes dividing N: the sequence is an eigen-sequence mod p (G_n = a lam^n, lam a root of x^2 - x - 1 mod p), so
     p never divides any term; periods mod p and mod small m versus the Pisano period pi(m).
  4. golden classes mod 9 (Q-orbits on primitive pairs, classified by N up to sign: {+-1}, {+-2}, {+-4} -- clock_duality
     O4) for Fibonacci, Lucas, (5,1), (7,1) and all (a,1); digital-root cycles; where the Solfeggio digits 3, 6, 9 occur.
  5. the norm values a^2 + a - 1 of the (a,1) family against the prime factors of the nine Solfeggio numbers (NUMEROLOGY
     check, with the split/inert type of each prime in Q(sqrt5)).
"""
from sympy import factorint, isprime, primerange

def seq(a, b, n):
    g = [a, b]
    while len(g) < n: g.append(g[-1] + g[-2])
    return g

def period(a, b, m):
    x, y = a % m, b % m; s = (x, y); k = 0
    while True:
        x, y = y, (x + y) % m; k += 1
        if (x, y) == s: return k

F = seq(0, 1, 200); L = seq(2, 1, 200)
def Fm(n): return F[n] if n >= 0 else (-1) ** (n + 1) * F[-n]
SOLF = [174, 285, 396, 417, 528, 639, 741, 852, 963]

for (a, b) in [(1, 1), (2, 1), (5, 1), (7, 1)]:
    G = seq(a, b, 60); N = b * b - a * b - a * a
    print('== start (%d,%d): N = b^2-ab-a^2 = %d = %s' % (a, b, N, factorint(abs(N)) if N else 0))
    print('   first 24 terms:', G[:24])
    assert all(G[n] == a * Fm(n - 1) + b * F[n] for n in range(1, 60)), 'closed form'
    assert all(G[n] ** 2 - G[n - 1] * G[n + 1] == (-1) ** (n + 1) * N for n in range(1, 58)), 'Cassini'
    # Lucas decomposition G_n = x L_n + y L_(n+1) (solve at n = 0, 1): 2x + y = a, x + 3y = b -> x = (3a-b)/5, y = (2b-a)/5
    if (2 * b - a) % 5 == 0:
        y = (2 * b - a) // 5; x = (3 * a - b) // 5
        assert all(G[n] == x * L[n] + y * L[n + 1] for n in range(58)); print('   Lucas form: G_n = %d L_n + %d L_(n+1)' % (x, y))
    else:
        print('   no integral Lucas form (5 does not divide 2b - a); G_n = %d F_(n-1) + %d F_n' % (a, b))
    primes = [p for p in factorint(abs(N))] if N else []
    for p in primes:
        roots = [x for x in range(p) if (x * x - x - 1) % p == 0]
        lam = (b * pow(a, -1, p)) % p if a % p else None
        ok = lam is not None and all(G[n] % p == a * pow(lam, n, p) % p for n in range(60))
        hit = any(g % p == 0 for g in G)
        print('   mod %d: roots of x^2-x-1 = %s; eigen-sequence G_n = %d * %s^n mod %d: %s; period %d (pi(%d) = %d); any term divisible by %d: %s'
              % (p, roots, a % p, lam, p, ok, period(a, b, p), p, period(0, 1, p), p, hit))
    for m in (2, 3, 4, 5, 8, 9, 10, 11, 27):
        print('   period mod %d: %d (Fibonacci pi = %d)' % (m, period(a, b, m), period(0, 1, m)), end=';')
    print()
    dr = [g % 9 or 9 for g in G[:24]]
    print('   digital roots (24-cycle):', dr, '; 3/6/9 at n mod 4 in', sorted({n % 4 for n in range(24) if dr[n] in (3, 6, 9)}))
    cls = {1: '+-1', 8: '+-1', 2: '+-2', 7: '+-2', 4: '+-4', 5: '+-4'}.get(N % 9, 'non-primitive')
    print('   golden class mod 9 (norm up to sign):', cls)

print('== golden class mod 9 of the (a,1) family, a = 1..12: N = 1 - a - a^2')
for a in range(1, 13):
    N = 1 - a - a * a
    print('   a=%2d N=%4d  class %s' % (a, N, {1: '+-1', 8: '+-1', 2: '+-2', 7: '+-2', 4: '+-4', 5: '+-4'}.get(N % 9, 'non-primitive')))
print('== (5,1) mod 9 is in the Q-orbit of 4*(Fibonacci)?', any(((4 * F[k]) % 9, (4 * F[k + 1]) % 9) == (5, 1) for k in range(24)),
      '; of 2*(Fibonacci)?', any(((2 * F[k]) % 9, (2 * F[k + 1]) % 9) == (5, 1) for k in range(24)),
      '; (7,1) = (F_9, F_10) mod 9?', (F[9] % 9, F[10] % 9) == (7, 1))
norms = {a * a + a - 1: a for a in range(1, 40)}
print('== Solfeggio prime factors (p > 3) vs the norms a^2 + a - 1 of the (a,1) starts (NUMEROLOGY check)')
for s in SOLF:
    fs = {p: e for p, e in factorint(s).items() if p > 3}
    desc = []
    for p in fs:
        t = 'ramified' if p == 5 else ('split' if p % 5 in (1, 4) else 'inert')
        desc.append('%d (%s%s)' % (p, t, ', = N(%d,1)' % norms[p] if p in norms else ''))
    print('   %d = %s ; large primes: %s' % (s, factorint(s), ', '.join(desc) if desc else 'none'))
split = [p for p in primerange(7, 1000)]
frac = sum(1 for p in split if p % 5 in (1, 4)) / len(split)
print('   share of split primes among primes 7..1000: %.3f (Chebotarev: 1/2)' % frac)
