#!/usr/bin/env python3
"""S43b referee: P_N(x)^2 = sum_k C(N+k,2k) C(2k,k)^2 ((x^2-1)/4)^k (Z.-H. Sun's identity (1.7), as cited by KLMSY),
exact polynomial identity for N <= 40 (sympy), and F_m == P_n(x)^2 (mod p^2), x^2 = 1 - 64/m, at the 8 CM points, p <= 300."""
import sympy as sp
from math import comb
x = sp.symbols('x')
ok = all(sp.expand(sp.legendre(N, x)**2 - sum(comb(N+k, 2*k)*comb(2*k, k)**2*((x**2-1)/4)**k for k in range(N+1))) == 0
         for N in range(41))
print('identity for 0 <= N <= 40:', ok)
y = sp.symbols('y')
bad = []
for p in [q for q in range(5, 301) if sp.isprime(q)]:
    n = (p-1)//2; pp = p*p
    P2 = sp.Poly(sp.expand(sp.legendre(n, x)**2).subs(x**2, y), y)   # even polynomial in x -> polynomial in y = x^2
    coeffs = P2.all_coeffs()[::-1]
    for m in [1, -8, 16, -64, 64, 256, -512, 4096]:
        im = pow(m % pp, -1, pp); yv = (1 - 64*im) % pp
        val = 0
        for i, c in enumerate(coeffs):
            c = sp.Rational(c); val = (val + int(c.p)*pow(int(c.q), -1, pp)*pow(yv, i, pp)) % pp
        F = sum(comb(2*k, k)**3*pow(im, k, pp) for k in range(n+1)) % pp
        if (F - val) % pp: bad.append((m, p))
print('F_m == P_n(x)^2 mod p^2 at 8 CM points, 5 <= p <= 300: failures', bad[:6] if bad else 'none')
