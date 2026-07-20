#!/usr/bin/env python3
"""gmc4_quadratic_decide_kps_S128c117.py -- kind-pasteur-2026-07-20-S128c117

THE DECISIVE TEST: are the surviving QUADRATIC candidates genuine counterexamples,
or do their Q-moments die?

Two candidates from the coefficient search kept E[P^m] = 0 all the way to m = 8, which
is much deeper than the m <= 3 sieve that found them.  But their first few Q-moments
read  [-1, 2, 0]  and  [-1, 4, 0]  -- with a ZERO at m = 3.  That is exactly the shape
of a sequence that is about to vanish permanently, and if E[Q P^m] = 0 for all large m
then the Mathieu property HOLDS for them and they are NOT counterexamples.

Being a counterexample requires:
    (a) E[P^m] = 0 for every m,  AND
    (b) SOME polynomial Q with E[Q P^m] != 0 for arbitrarily large m.
Condition (b) must be tested against a RANGE of Q, not just Q = Z2 -- the Mathieu
condition quantifies over all Q, so failing for one Q proves nothing.

Also hardened here: the 5-TERM counterexample found by dropping a term from the
published 6-term P, which if it survives is a strict improvement on the source.
"""
import sys
from math import factorial
from sympy import symbols, expand, Poly, Rational

MM = int(sys.argv[1]) if len(sys.argv) > 1 else 12
z1, z1c, z2, z2c = symbols('z1 z1c z2 z2c')
V = (z1, z1c, z2, z2c)


def E(expr):
    p = Poly(expand(expr), *V)
    tot = 0
    for (a, b, c, d), co in p.terms():
        if a == b and c == d:
            tot += co * factorial(a) * factorial(c)
    return tot


QS = [("Z2", z2), ("Z1", z1), ("Z2b", z2c), ("Z1b", z1c),
      ("Z1Z2", z1 * z2), ("Z1bZ2b", z1c * z2c), ("Z1Z2b", z1 * z2c),
      ("|Z1|^2", z1 * z1c), ("Z2^2", z2**2), ("1", 1)]

cands = {
    "quad-2": z1 * z1c + z1 * z2c + z1 - z1c * z2 - z1c * z2c - z1c
              - z2 * z2c - z2 - z2c**2 - z2c,
    "quad-3": z1**2 + z1 * z1c + z1 * z2 + z1 - z1c * z2c - z1c - z2 * z2c
              - z2 - z2c**2 - z2c,
}

print("=" * 78)
print("THE SURVIVING QUADRATICS: do their Q-moments die?")
print("=" * 78)
for name, e in cands.items():
    deg = Poly(expand(e), *V).total_degree()
    zs = [E(expand(e**m)) for m in range(1, MM + 1)]
    print()
    print("  %s  (total degree %d)" % (name, deg))
    print("     E[P^m], m=1..%d : %s" % (MM, zs))
    print("     E[P^m] = 0 for all m tested : %s" % all(v == 0 for v in zs))
    powers = [expand(e**m) for m in range(1, MM + 1)]
    alive = []
    for qn, q in QS:
        qs = [E(expand(q * pw)) for pw in powers]
        nz = [m for m, v in enumerate(qs, 1) if v != 0]
        last = max(nz) if nz else None
        print("     Q = %-8s E[QP^m] = %-42s last non-zero m: %s"
              % (qn, str(qs[:7]), last))
        if last is not None and last >= MM - 2:
            alive.append(qn)
    print("     Q's still non-zero near m = %d : %s" % (MM, alive if alive else "NONE"))
    if not alive:
        print("     -> every Q tested has E[Q P^m] = 0 for large m, i.e. the MATHIEU")
        print("        PROPERTY HOLDS here.  NOT a counterexample.  The m <= 3 sieve")
        print("        that produced it was too shallow, and the hit is WITHDRAWN.")
    else:
        print("     -> Q-moments persist: this WOULD be a quadratic counterexample.")

print()
print("=" * 78)
print("THE 5-TERM COUNTEREXAMPLE (drop one term from the published 6)")
print("=" * 78)
P6 = (1 + z2) * (z1c * (1 - z1) + z2c)
terms = Poly(expand(P6), *V).terms()


def build(sub):
    out = 0
    for mono, co in sub:
        t = co
        for v, k in zip(V, mono):
            t *= v**k
        out += t
    return out


for i in (2, 3):
    Pd = build([terms[j] for j in range(len(terms)) if j != i])
    zs = [E(expand(Pd**m)) for m in range(1, MM + 1)]
    print()
    print("  drop term %d -> P = %s" % (i, Pd))
    print("     %d terms, degree %d" % (len(Poly(expand(Pd), *V).terms()),
                                        Poly(expand(Pd), *V).total_degree()))
    print("     E[P^m] = 0 for m=1..%d : %s" % (MM, all(v == 0 for v in zs)))
    powers = [expand(Pd**m) for m in range(1, MM + 1)]
    for qn, q in QS[:6]:
        qs = [E(expand(q * pw)) for pw in powers]
        nz = [m for m, v in enumerate(qs, 1) if v != 0]
        print("     Q = %-8s E[QP^m] = %-40s last non-zero: %s"
              % (qn, str(qs[:6]), max(nz) if nz else None))
