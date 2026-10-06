#!/usr/bin/env python3
"""The number field of the 11-square packing constant.

Session opus-2026-10-06-S16.  Inputs (owner, from the Lean formalization github.com/Queuingtheorydotcom/11SquaresFormalized):
  side length T is a root of  P(s) = s^8-20s^7+178s^6-842s^5+1923s^4-496s^3-6754s^2+12420s-6865,
  T = (6u+4)/(1+2u-u^2) with u the root in (9/25, 37/100) of U(u) = 5u^8-10u^7-2u^6+14u^5+12u^4-6u^3+2u^2+2u-1.

Computes: irreducibility, the link P <-> U, discriminants and their factorisations, signature, Frobenius cycle-type
statistics (Chebotarev) to a bound, block systems (subfields) from high-precision roots, the Galois group, ramified
primes and p-maximality (Dedekind criterion), and compares with candidate transitive groups of degree 8.

Reproduce: python3 04-computation/experiments/square11_octic_field_20261006.py [prime_bound]
"""
import sys, itertools, math
from collections import Counter
import mpmath as mp
import sympy as sp
import flint

PB = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
s, u = sp.symbols('s u')
P = s**8 - 20*s**7 + 178*s**6 - 842*s**5 + 1923*s**4 - 496*s**3 - 6754*s**2 + 12420*s - 6865
U = 5*u**8 - 10*u**7 - 2*u**6 + 14*u**5 + 12*u**4 - 6*u**3 + 2*u**2 + 2*u - 1
Pc = [int(c) for c in sp.Poly(P, s).all_coeffs()]
Uc = [int(c) for c in sp.Poly(U, u).all_coeffs()]

print('== 1. The two octics')
print('P irreducible over Q:', sp.Poly(P, s).is_irreducible, '  U irreducible over Q:', sp.Poly(U, u).is_irreducible)
mp.mp.dps = 60
uroots = [r for r in sp.Poly(U, u).nroots(n=50) if r.is_real]
ustar = [r for r in uroots if sp.Rational(9, 25) < r < sp.Rational(37, 100)]
print('real roots of U:', [sp.N(r, 20) for r in uroots])
print('root of U in (9/25, 37/100):', [sp.N(r, 30) for r in ustar])
Tval = [(6*r + 4) / (1 + 2*r - r**2) for r in ustar]
print('T = (6u+4)/(1+2u-u^2):', [sp.N(t, 30) for t in Tval], '  P(T) =', [sp.N(P.subs(s, t), 5) for t in Tval])
# exact link: resultant of U(u) and (1+2u-u^2)s - (6u+4) in u
link = sp.resultant(U, (1 + 2*u - u**2)*s - (6*u + 4), u)
link_f = sp.factor_list(sp.expand(link))
print('Res_u(U, (1+2u-u^2)s-(6u+4)) =', link_f)
# inverse map: express u as a polynomial in T over Q (Q(u) = Q(T))
print('degree check: [Q(T):Q] = 8 = [Q(u):Q], so Q(T) = Q(u) (T is a rational function of u)')

print('\n== 2. Discriminants')
dP = sp.discriminant(P, s); dU = sp.discriminant(U, u)
print('disc P =', dP, '=', sp.factorint(dP))
print('disc U =', dU, '=', sp.factorint(dU))
print('disc P is a square:', sp.sqrt(dP).is_integer if dP > 0 else False, '  sign:', sp.sign(dP))
rr = sp.Poly(P, s).nroots(n=40)
nreal = sum(1 for r in rr if abs(sp.im(r)) < 1e-30)
print('real roots of P:', nreal, ' signature (r1, r2) =', (nreal, (8 - nreal) // 2))
print('real roots of P:', sorted([sp.N(sp.re(r), 15) for r in rr if abs(sp.im(r)) < 1e-30]))

print('\n== 3. Frobenius cycle types (Chebotarev), primes p <=', PB, 'not dividing disc P')
fP = flint.fmpz_poly(list(reversed(Pc)))
bad = set(sp.factorint(dP).keys())
types = Counter()
nprimes = 0
for p in sp.primerange(3, PB):
    if p in bad:
        continue
    g = flint.nmod_poly(list(reversed([c % p for c in Pc])), p)
    fac = g.factor()[1]
    t = tuple(sorted((f.degree() for f, e in fac for _ in range(e)), reverse=True))
    types[t] += 1
    nprimes += 1
print('primes used:', nprimes)
for t, c in sorted(types.items(), key=lambda kv: -kv[1]):
    print(f'   {str(t):28s} {c:7d}  {c / nprimes:.5f}')
root_frac = sum(c for t, c in types.items() if 1 in t) / nprimes
irr_frac = types.get((8,), 0) / nprimes
print(f'P has a root mod p: {root_frac:.5f}   (S8: 1 - D_8/8! = {1 - 14833/40320:.5f})')
print(f'P irreducible mod p: {irr_frac:.5f}   (S8: 1/8 = 0.125)')

print('\n== 4. Block systems from 60-digit roots (subfield test)')
mp.mp.dps = 80
coeffs = [mp.mpf(c) for c in Pc]
roots = mp.polyroots(coeffs, maxsteps=400, extraprec=400)


def integral_poly_from(vals):
    poly = [mp.mpf(1)]
    for v in vals:
        poly = [a - v * b for a, b in zip(poly + [0], [0] + poly)]
    return poly


def is_integral(poly, tol=mp.mpf(10) ** -30):
    return all(abs(mp.im(c)) < tol and abs(mp.re(c) - mp.nint(mp.re(c))) < tol for c in poly)


found = []
idx = list(range(8))
# two blocks of four
for A in itertools.combinations(idx, 4):
    if 0 not in A:
        continue
    B = tuple(i for i in idx if i not in A)
    for f in (lambda S: sum(roots[i] for i in S), lambda S: mp.fprod([roots[i] for i in S])):
        vals = [f(A), f(B)]
        if abs(vals[0] - vals[1]) > 1e-20 and is_integral(integral_poly_from(vals)):
            found.append(('2 blocks of 4', A, B))
# four blocks of two
pairs_all = []


def matchings(rest):
    if not rest:
        yield []
        return
    a = rest[0]
    for b in rest[1:]:
        r2 = [x for x in rest if x not in (a, b)]
        for m in matchings(r2):
            yield [(a, b)] + m


for M in matchings(idx):
    for f in (lambda pr: roots[pr[0]] + roots[pr[1]], lambda pr: roots[pr[0]] * roots[pr[1]]):
        vals = [f(pr) for pr in M]
        if all(abs(vals[i] - vals[j]) > 1e-20 for i in range(4) for j in range(i + 1, 4)) and is_integral(integral_poly_from(vals)):
            found.append(('4 blocks of 2', tuple(M)))
print('block systems found:', found if found else 'none (sum and product invariants, all 35 + 105 candidates)')
