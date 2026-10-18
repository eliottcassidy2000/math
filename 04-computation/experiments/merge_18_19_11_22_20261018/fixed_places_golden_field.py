#!/usr/bin/env python3
"""Conjugation-fixed primes of the golden contact field F = Q(zeta15)(sqrt-11)(w), w^3 + 2w^2 - 4w - 4 = 0
(mac-mini-2026-10-18-sixcolour; THM-4558 Lemma R screen of the six-colour programme).

Complex conjugation c acts by zeta -> zeta^-1, j -> -j, w -> w (w real). For p not dividing 2*3*5*11*37:
  * Q(zeta15): a prime over p is c-fixed iff -1 lies in <p mod 15>, i.e. iff p = 14 (mod 15); then f = 2 and c acts
    on the residue field F_(p^2) as Frobenius.
  * K = Q(zeta15, j): over such a prime, x^2 + 11 has roots +-s in F_(p^2); c maps the prime (j - s) to (j + s^p),
    so it is fixed iff s^p = -s iff (-11/p) = -1 iff (p/11) = -1.
  * F = K(w): over a c-fixed prime of K, the prime attached to a root r of the cubic in F_(p^2) is fixed iff
    r^p = r, i.e. r in F_p. If the cubic is irreducible mod p, the unique prime above has residue F_(p^6) and c acts
    as z -> z^(p^3).
Lemma R (THM-4558): a c-fixed prime, inert over the real subfield, with residue F_(q^2), gives chi(G(F)) <= kappa(q),
kappa(q) = chi(Cay(F_(q^2), mu_(q+1))). The ramified primes 2, 3, 5, 11, 37 are checked separately (none is c-fixed in
F: 2, 37 are not fixed in Q(zeta15); 3, 5 split in K into conjugate pairs; 11 is not fixed in Q(zeta5))."""
from sympy import primerange, legendre_symbol

def cubic_roots_mod(p):
    return [x for x in range(p) if (x ** 3 + 2 * x * x - 4 * x - 4) % p == 0]

print('p | p mod 15 | (p/11) | cubic roots mod p | c-fixed primes of F above p : residue field, Lemma R bound')
rows = []
for p in primerange(7, 3000):
    if p in (11, 37): continue
    if p % 15 != 14: continue
    leg = legendre_symbol(p % 11, 11)
    roots = cubic_roots_mod(p)
    # K = Q(zeta15, j) has [K:Q]/(e f) = 16/2 = 8 primes above such p (f = 2), each c-fixed when (p/11) = -1; over each,
    # exactly one prime of F is c-fixed: the one attached to the F_p root of the cubic (residue F_(p^2)) if the cubic
    # has a root mod p, or the unique prime above it (residue F_(p^6)) if the cubic is irreducible mod p.
    if leg == 1:
        desc = 'none (j-primes swapped by c)'
    elif roots:
        desc = '8 c-fixed primes with residue F_(%d^2): chi(G(F)) <= kappa(%d)' % (p, p)
        rows.append((p, len(roots)))
    else:
        desc = '8 c-fixed primes with residue F_(%d^6): chi(G(F)) <= kappa(%d^3)' % (p, p)
    print('%5d | %2d | %+d | %s | %s' % (p, p % 15, leg, roots, desc))
print('smallest q with a residue-F_(q^2) fixed place:', rows[:6])
