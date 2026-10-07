#!/usr/bin/env python3
"""Is the full base-phi->binary reading F(n) prime more often for prime n?  Mechanism hunt."""
import sys
from math import log
from sympy import isprime
sys.setrecursionlimit(10000)
from phi_binary_primes import strings
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
F = {}
for n in range(2, NMAX + 1):
    ip, fp = strings(n)
    F[n] = (int(ip + fp, 2), len(ip), len(fp))
def stats(ns, label):
    c = sum(isprime(F[n][0]) for n in ns)
    exp = sum(2 / log(F[n][0]) for n in ns)   # F is odd: heuristic 2/ln F
    print(f'  {label:28s} N={len(ns):6d}  P(F prime)={c/len(ns):.4f}  heuristic 2/lnF={exp/len(ns):.4f}')
ps = [n for n in F if isprime(n)]
cs = [n for n in F if not isprime(n)]
stats(ps, 'n prime'); stats(cs, 'n composite')
stats([n for n in cs if n % 2 == 0], 'n even composite')
stats([n for n in cs if n % 2 == 1], 'n odd composite')
for q in (2, 3, 5, 7, 11, 13):
    a = sum(F[n][0] % q == 0 for n in F if n % q == 0) / max(1, sum(1 for n in F if n % q == 0))
    b = sum(F[n][0] % q == 0 for n in F if n % q != 0) / max(1, sum(1 for n in F if n % q != 0))
    print(f'  P({q} | F(n)): given {q}|n: {a:.4f}   given {q}!|n: {b:.4f}')
# does n | F(n)?  and gcd structure
import math
g = sum(1 for n in F if F[n][0] % n == 0)
print('  #n with n | F(n):', g, ' examples', [n for n in list(F)[:3000] if F[n][0] % n == 0][:15])
gc = sum(1 for n in cs if math.gcd(F[n][0], n) > 1) / len(cs)
print(f'  P(gcd(F(n), n) > 1 | n composite) = {gc:.4f};  P(gcd(F(p),p)>1 | p prime) =', sum(1 for n in ps if math.gcd(F[n][0], n) > 1) / len(ps))
