#!/usr/bin/env python3
"""gmc4_extend_kps_S128c117.py -- kind-pasteur-2026-07-20-S128c117

EXTENDING THE VERIFIED GMC(4) COUNTEREXAMPLE.

Confirmed in gmc4_verify: with Z1, Z2 standard complex Gaussians (4 real) and
W_j = conj(Z_j),
    P = (1+Z2)(conj(Z1)(1-Z1) + conj(Z2)),  Q = Z2
gives E[P^m] = 0 and E[Q P^m] = m! for every m.

FOUR THINGS TRIED HERE.

(A) THE GENERATING-FUNCTION FORM.  E[P^m] = 0 for all m >= 1 says exactly
        E[exp(tP)] = 1,
    and E[Q P^m] = m! says
        E[Q exp(tP)] = sum_m t^m = 1/(1-t).
    So the counterexample is: a random variable whose exponential moment generating
    function is IDENTICALLY 1, but which acquires a SIMPLE POLE AT t = 1 after tilting
    by Q.  That is a much sharper way to state it than the moment list, and it says the
    obstruction is a pole at a specific finite t, not a growth rate.  Checked as a
    truncated series.

(B) IS THE "1" IN EACH FACTOR FORCED?  Scan P(a,b) = (1 + a*Z2)(conj(Z1)(1 - b*Z1) +
    conj(Z2)) over small rationals a, b.  If E[P^m] = 0 only at a = b = 1 the
    construction is rigid; if a whole curve works, there is a family.

(C) IS CUBIC / 6 TERMS MINIMAL AT n = 4?  Every 5-term subpolynomial of the expanded P
    is tested (a term is load-bearing iff dropping it breaks E[P^m] = 0), and the
    charge-balanced quadratics are searched directly.

(D) DOES THE CHAIN EXTEND?  Try the 3-complex-variable analogue
        P3 = (1+Z3)( (1+Z2)(conj(Z1)(1-Z1) + conj(Z2)) + conj(Z3) )
    with Q = Z3, to see whether a longer chain gives a different moment sequence --
    i.e. whether m! is the first term of a family or an isolated value.
"""
import sys
from math import factorial
from itertools import combinations
from fractions import Fraction as F
from sympy import symbols, expand, Poly, Rational, simplify, series, exp, nsimplify

MM = int(sys.argv[1]) if len(sys.argv) > 1 else 6

z1, z1c, z2, z2c, z3, z3c = symbols('z1 z1c z2 z2c z3 z3c')
VARS4 = (z1, z1c, z2, z2c)
VARS6 = (z1, z1c, z2, z2c, z3, z3c)


def E(expr, vs=VARS4):
    p = Poly(expand(expr), *vs)
    tot = 0
    for mono, co in p.terms():
        ok = True
        f = 1
        for k in range(0, len(vs), 2):
            if mono[k] != mono[k + 1]:
                ok = False
                break
            f *= factorial(mono[k])
        if ok:
            tot += co * f
    return tot


P = (1 + z2) * (z1c * (1 - z1) + z2c)
Q = z2

print("=" * 78)
print("(A) GENERATING-FUNCTION FORM")
print("=" * 78)
print("  E[P^m] = 0 for all m  <=>  E[exp(tP)] = 1")
print("  E[Q P^m] = m!         <=>  E[Q exp(tP)] = sum t^m = 1/(1-t)")
print()
print("  %-4s %-14s %-14s" % ("m", "E[P^m]/m!", "E[Q P^m]/m!"))
for m in range(0, MM + 1):
    Pm = expand(P**m)
    a = E(Pm)
    b = E(expand(Q * Pm))
    print("  %-4d %-14s %-14s" % (m, Rational(a, factorial(m)), Rational(b, factorial(m))))
print()
print("  So E[exp(tP)] = 1 exactly, and E[Q exp(tP)] = 1 + t + t^2 + ... = 1/(1-t):")
print("  an mgf that is identically 1, acquiring a SIMPLE POLE AT t = 1 when tilted by Q.")
print("  The failure is located at a finite t, not at infinity.")

print()
print("=" * 78)
print("(B) IS THE CONSTRUCTION RIGID?  P(a,b) = (1 + a Z2)(conj(Z1)(1 - b Z1) + conj(Z2))")
print("=" * 78)
cands = [Rational(x) for x in (-2, -1, Rational(-1, 2), 0, Rational(1, 2), 1, 2)]
good = []
for a in cands:
    for b in cands:
        Pab = (1 + a * z2) * (z1c * (1 - b * z1) + z2c)
        if all(E(expand(Pab**m)) == 0 for m in range(1, 5)):
            qs = [E(expand(Q * expand(Pab**m))) for m in range(1, 5)]
            good.append((a, b, qs))
print("  (a,b) with E[P^m] = 0 for m = 1..4:")
for a, b, qs in good:
    print("     a = %-6s b = %-6s   E[Q P^m] = %s" % (a, b, qs))
print("  -> %d parameter pairs work." % len(good))

print()
print("=" * 78)
print("(C) MINIMALITY: is every term of the 6-term P load-bearing?")
print("=" * 78)
pe = Poly(expand(P), *VARS4)
terms = pe.terms()
print("  P has %d terms." % len(terms))
from sympy import Mul, Pow
def build(sub):
    e = 0
    for mono, co in sub:
        t = co
        for v, k in zip(VARS4, mono):
            t *= v**k
        e += t
    return e
drop_ok = []
for i in range(len(terms)):
    sub = [terms[j] for j in range(len(terms)) if j != i]
    Pd = build(sub)
    if all(E(expand(Pd**m)) == 0 for m in range(1, 5)):
        drop_ok.append(i)
print("  terms whose removal PRESERVES E[P^m] = 0 (m<=4) : %s" % (drop_ok if drop_ok else "none"))
print("  -> all six terms load-bearing : %s" % (not drop_ok))

print()
print("  charge-balanced QUADRATIC search (all P with monomials of total degree <= 2,")
print("  coefficients in {-1,0,1}, at least one term):")
mons2 = []
for a in range(3):
    for b in range(3):
        for c in range(3):
            for d in range(3):
                if 1 <= a + b + c + d <= 2:
                    mons2.append((a, b, c, d))
found2 = []
import itertools
for coef in itertools.product((-1, 0, 1), repeat=len(mons2)):
    if all(c == 0 for c in coef):
        continue
    e = 0
    for co, mono in zip(coef, mons2):
        if co:
            t = co
            for v, k in zip(VARS4, mono):
                t *= v**k
            e += t
    if all(E(expand(e**m)) == 0 for m in range(1, 4)):
        q = [E(expand(Q * expand(e**m))) for m in range(1, 4)]
        if any(x != 0 for x in q):
            found2.append((e, q))
            if len(found2) >= 3:
                break
print("  quadratic P with E[P^m]=0 (m<=3) AND some E[Q P^m] != 0 : %d found"
      % len(found2))
for e, q in found2[:3]:
    print("     P = %s   E[Q P^m] = %s" % (e, q))

print()
print("=" * 78)
print("(D) DOES THE CHAIN EXTEND?  three complex variables")
print("=" * 78)
P3 = (1 + z3) * ((1 + z2) * (z1c * (1 - z1) + z2c) + z3c)
Q3 = z3
print("  P3 = (1+Z3)((1+Z2)(conj(Z1)(1-Z1)+conj(Z2)) + conj(Z3))")
rows = []
for m in range(1, 5):
    P3m = expand(P3**m)
    a = E(P3m, VARS6)
    b = E(expand(Q3 * P3m), VARS6)
    rows.append((m, a, b))
    print("     m = %-3d  E[P3^m] = %-10s  E[Q3 P3^m] = %-12s  m! = %d"
          % (m, a, b, factorial(m)))
print("  chain still has E[P^m] = 0 : %s" % all(r[1] == 0 for r in rows))
print("  Q-moments equal m! again   : %s" % all(r[2] == factorial(r[0]) for r in rows))
