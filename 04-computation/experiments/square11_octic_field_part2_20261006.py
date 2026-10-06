#!/usr/bin/env python3
"""Part 2: certify Gal(P) = S8 with explicit primes (Dedekind + Jordan), compute the field discriminant by
Dedekind's criterion, factor P at the ramified primes, the quadratic resolvent and its class number, and the
exact Artin-type 'entanglement' of root counts with the sign character.

Reproduce: python3 04-computation/experiments/square11_octic_field_part2_20261006.py
"""
import math
from collections import Counter
from fractions import Fraction as Fr
import sympy as sp
import flint

Pc = [1, -20, 178, -842, 1923, -496, -6754, 12420, -6865]          # high -> low
P = flint.fmpz_poly(list(reversed(Pc)))
s = sp.symbols('s')
Psym = sum(c * s ** (8 - i) for i, c in enumerate(Pc))
dP = int(sp.discriminant(Psym, s))
fac = sp.factorint(dP)
print('disc P =', dP, fac)


def pattern(p):
    g = flint.nmod_poly([c % p for c in reversed(Pc)], p)
    f = g.factor()[1]
    return tuple(sorted(((h.degree(), e) for h, e in f), reverse=True))


print('\n== Gal(P) = S8 certificate (Dedekind: unramified Frobenius cycle type = factorisation type mod p)')
want = {'(7,1)': None, '(2,1^6)': None, '(8)': None, '(5,3)': None}
for p in sp.primerange(3, 100000):
    if dP % p == 0:
        continue
    t = pattern(p)
    if any(e > 1 for d, e in t):
        continue
    degs = tuple(d for d, e in t)
    if degs == (7, 1) and want['(7,1)'] is None:
        want['(7,1)'] = p
    if degs == (2, 1, 1, 1, 1, 1, 1) and want['(2,1^6)'] is None:
        want['(2,1^6)'] = p
    if degs == (8,) and want['(8)'] is None:
        want['(8)'] = p
    if degs == (5, 3) and want['(5,3)'] is None:
        want['(5,3)'] = p
    if all(v is not None for v in want.values()):
        break
print('smallest unramified primes with these factorisation types:', want)
print('Proof: P irreducible => G transitive; Frobenius at p =', want['(7,1)'], 'is a 7-cycle => G 2-transitive => primitive;',
      'Frobenius at p =', want['(2,1^6)'], 'is a transposition; a primitive group containing a transposition is S_n (Jordan). So G = S8.')

print('\n== Dedekind criterion: is Z[T] p-maximal?')


def dedekind_maximal(p):
    gp = flint.nmod_poly([c % p for c in reversed(Pc)], p)
    lead, facs = gp.factor()
    G = flint.fmpz_poly([1])
    H = flint.fmpz_poly([1])
    for h, e in facs:
        hz = flint.fmpz_poly([int(c) for c in h.coeffs()])
        G = G * hz
        for _ in range(e - 1):
            H = H * hz
    Fz = (G * H - P)
    coeffs = [int(c) for c in Fz.coeffs()]
    assert all(c % p == 0 for c in coeffs)
    Fbar = flint.nmod_poly([c // p % p for c in coeffs], p)
    Gbar = flint.nmod_poly([int(c) % p for c in G.coeffs()], p)
    Hbar = flint.nmod_poly([int(c) % p for c in H.coeffs()], p)
    g = Fbar.gcd(Gbar).gcd(Hbar)
    return g.degree() == 0, facs


for p in sorted(fac):
    if p < 0:
        continue
    ok, facs = dedekind_maximal(p)
    print(f'   p = {p:>7}: v_p(disc P) = {fac[p]}, factorisation mod p = {[(h.degree(), e) for h, e in facs]}, Z[T] p-maximal: {ok}')

print('\n   Conclusion for d_K (disc P = index^2 * d_K):')
print('   p-maximal at p  =>  v_p(d_K) = v_p(disc P); not p-maximal with v_p(disc P) = 2  =>  v_p(d_K) = 0 (index divisible by p).')

print('\n== Quadratic resolvent (sign character of S8 = the GL(1) piece)')
sq = -1
for p, e in fac.items():
    if p > 0 and e % 2 == 1:
        sq *= p
D0 = sq
D = D0 if D0 % 4 == 1 else 4 * D0
print('squarefree part of disc P:', D0, '  quadratic resolvent field Q(sqrt(', D0, ')), fundamental discriminant', D)


def class_number_negative(Dd):
    """h(D) for a negative fundamental discriminant by counting reduced primitive forms (a,b,c), b^2-4ac = D."""
    h = 0
    a = 1
    while 3 * a * a <= -Dd:
        for b in range(-a + 1, a + 1):
            if (b * b - Dd) % (4 * a):
                continue
            c = (b * b - Dd) // (4 * a)
            if c < a:
                continue
            if b < 0 and (a == c):
                continue
            if math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    return h


hD = class_number_negative(D)
print('class number h(Q(sqrt(D))) =', hD)
# Kronecker symbol check: (D/p) = +1 iff Frobenius at p is an even permutation
mism = 0
cnt = 0
for p in sp.primerange(3, 20000):
    if dP % p == 0:
        continue
    t = pattern(p)
    degs = [d for d, e in t]
    parity = sum(d - 1 for d in degs) % 2          # parity of the permutation
    chi = sp.jacobi_symbol(D % p, p) if p != 2 else None
    cnt += 1
    if (chi == 1) != (parity == 0):
        mism += 1
print(f'sign of Frobenius = Kronecker (D/p) on {cnt} primes: mismatches {mism}')

print('\n== Artin-type entanglement: root count of P mod p, conditioned on the sign character (exact, S8)')
# number of permutations of S8 with k fixed points, split by parity
from itertools import permutations
counts = Counter()
for perm in permutations(range(8)):
    k = sum(1 for i, x in enumerate(perm) if i == x)
    # parity
    seen = [False] * 8
    par = 0
    for i in range(8):
        if not seen[i]:
            j = i
            L = 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                L += 1
            par += L - 1
    counts[(k, par % 2)] += 1
for k in range(9):
    e, o = counts[(k, 0)], counts[(k, 1)]
    if e + o:
        print(f'   {k} roots mod p: density {Fr(e + o, 40320)} = {(e + o) / 40320:.6f};  given (D/p)=+1: {Fr(e, 20160)} = {e / 20160:.6f};  given (D/p)=-1: {Fr(o, 20160)} = {o / 20160:.6f}')
