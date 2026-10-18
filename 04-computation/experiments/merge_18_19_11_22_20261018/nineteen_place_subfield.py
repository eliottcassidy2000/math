#!/usr/bin/env python3
"""Why the 19-adic separation of THM-4690(5) happens (mac-mini-2026-10-18-sixcolour): in the subfield
F' = Q(zeta5, sqrt-11, w) of F = Q(zeta15, sqrt-11, w), the primes above 19 of type h = w + 3 = 7 are fixed by complex
conjugation (zeta5 -> zeta5^-1, j -> -j, w -> w) and have residue field F_361 with conjugation acting as Frobenius.
So Lemma R (THM-4558) applies to G(F'): chi(G(F')) <= kappa(19) = 5, and every unit vector of F' is a unit there.
In F the same primes are NOT fixed (19 = 4 mod 15, <4> = {1, 4} does not contain -1).
Checks, on the 48 homomorphisms Z[zeta15, j, w] -> F_361 grouped into Frobenius orbits:
  * for F: no orbit is conjugation-stable;
  * for F' (restrict zeta15 -> zeta5 = zeta15^3): an orbit is stable iff its w-image is 4 (h = 7).
Also prints the clock-tower denominators of the exact ideal energies (S41 audit a4 / clock_duality root_clock.json)."""
P = 19; NR = (-11) % P
def fm(x, y): return ((x[0] * y[0] + NR * x[1] * y[1]) % P, (x[0] * y[1] + x[1] * y[0]) % P)
def fpow(x, k):
    r = (1, 0)
    for _ in range(k): r = fm(r, x)
    return r
def add(x, y): return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)
els = [(a, b) for a in range(P) for b in range(P)]
z15 = [x for x in els if fpow(x, 15) == (1, 0) and fpow(x, 5) != (1, 0) and fpow(x, 3) != (1, 0)]
jr = [(0, 1), (0, P - 1)]
def cub(x): return add(add(fpow(x, 3), fm((2, 0), fpow(x, 2))), add(fm((P - 4, 0), x), (P - 4, 0)))
wr = [x for x in els if cub(x) == (0, 0)]
frob = lambda t: tuple(fpow(c, P) for c in t)
def orbit(t): return frozenset([t, frob(t)])
neg = lambda x: ((-x[0]) % P, (-x[1]) % P)
conjF = lambda t: (fpow(t[0], 14), neg(t[1]), t[2])                  # zeta -> zeta^-1, j -> -j, w -> w
homsF = [(z, j, w) for z in z15 for j in jr for w in wr]
stableF = [t for t in homsF if orbit(conjF(t)) == orbit(t)]
print('F : primes above 19 = %d, conjugation-stable = %d' % (len({orbit(t) for t in homsF}), len({orbit(t) for t in stableF})))
homsFp = {(fpow(z, 3), j, w) for (z, j, w) in homsF}                 # restrict to zeta5 = zeta15^3
conjFp = lambda t: (fpow(t[0], 4), neg(t[1]), t[2])                  # zeta5 -> zeta5^-1 = zeta5^4
orbsFp = {orbit(t) for t in homsFp}
stab = [o for o in orbsFp if orbit(conjFp(next(iter(o)))) == o]
types = sorted({next(iter(o))[2] for o in stab})
print("F': primes above 19 = %d, conjugation-stable = %d, their w-images %s (w = 4 means h = 7)" % (len(orbsFp), len(stab), types))
for o in stab:
    t = next(iter(o)); c = conjFp(t)
    print('   stable prime', sorted(o), ': conj(hom) = Frobenius(hom)?', c == frob(t) or c == t)
# clock-tower denominators of the exact ideal energies
from fractions import Fraction as Fr
E = [Fr(34, 7), Fr(32898, 9709), Fr(572764898734162, 222399981598543)]
for d, e in enumerate(E):
    den = e.denominator
    print('E_%d = %s : denominator %d ; (2^(2*3^(d+1)) - 1) / den = %s' % (d, e, den, Fr(2 ** (2 * 3 ** (d + 1)) - 1, den)))
print('9709 = (2^18-1)/27:', (2 ** 18 - 1) // 27 == 9709, '; 87381 = (2^18-1)/3:', (2 ** 18 - 1) // 3 == 87381, '; 222399981598543 = (2^54-1)/81:', (2 ** 54 - 1) // 81 == 222399981598543)
