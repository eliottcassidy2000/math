#!/usr/bin/env python3
"""gmc4_four_term_mechanism_kps_S128c117.py -- kind-pasteur-2026-07-20-S128c117

THE 4-TERM COUNTEREXAMPLE, AND WHY IT WORKS.

Searching subsets of the published 6-term P found a 4-TERM survivor, and it factors:

    P4 = (1 + Z)(conj(Z) - Y),    Y := |Z1|^2,   Z := Z2

so Z1 enters ONLY through |Z1|^2, which is Exp(1).  Four terms, still cubic, still
4 real Gaussians, still E[P^m] = 0 and E[Z P^m] = m!.  The source's "cubic, 6 terms"
is therefore not minimal.

THE MECHANISM IS AN ALTERNATING BINOMIAL IDENTITY.  With Z a standard complex Gaussian
and Y independent with moments mu_i = E[Y^i],

    E[P^m] = sum_{i=0}^{m} (-1)^i C(m,i)^2 (m-i)! mu_i           (*)

because E[Z^j conj(Z)^k] = delta_{jk} j!.  Now if mu_i = i! then

    C(m,i)^2 (m-i)! i! = m! C(m,i),

so (*) collapses to  m! * sum_i (-1)^i C(m,i) = 0  for every m >= 1.  The vanishing is
the binomial theorem, nothing more -- and it needs mu_i = i! exactly.

AND THAT FORCES Y.  Reading (*) as a triangular linear system in mu determines mu_i
one at a time: mu_1 = 1, mu_2 = 2, mu_3 = 6, ... i.e. mu_i = i! is the UNIQUE moment
sequence making E[P^m] vanish for all m.  So Y must be Exp(1), and |Z1|^2 is its
natural realisation -- which costs exactly TWO real Gaussians.  With two more for Z,
n = 4 is not incidental to this construction; it is forced by it.

Checked here: the factorisation, the moment identity (*), the forcing of mu_i = i!,
and whether a single real Gaussian could supply Y (it cannot).
"""
import sys
from math import factorial
from sympy import symbols, expand, Poly, Rational, binomial, simplify, nsimplify

MM = int(sys.argv[1]) if len(sys.argv) > 1 else 10
z1, z1c, z2, z2c = symbols('z1 z1c z2 z2c')
V = (z1, z1c, z2, z2c)


def E(expr):
    p = Poly(expand(expr), *V)
    t = 0
    for (a, b, c, d), co in p.terms():
        if a == b and c == d:
            t += co * factorial(a) * factorial(c)
    return t


P4 = (1 + z2) * (z2c - z1 * z1c)
print("=" * 78)
print("(1) THE 4-TERM WITNESS  P4 = (1 + Z2)(conj(Z2) - |Z1|^2)")
print("=" * 78)
pe = Poly(expand(P4), *V)
print("  expanded: %s" % expand(P4))
print("  terms = %d, total degree = %d" % (len(pe.terms()), pe.total_degree()))
zs = [E(expand(P4**m)) for m in range(1, MM + 1)]
qs = [E(expand(z2 * expand(P4**m))) for m in range(1, MM + 1)]
print("  E[P^m], m=1..%d : %s" % (MM, zs))
print("  E[Z2 P^m]        : %s" % qs)
print("  E[P^m] = 0 for all m : %s" % all(v == 0 for v in zs))
print("  E[Z2 P^m] = m!       : %s" % all(qs[m - 1] == factorial(m)
                                          for m in range(1, MM + 1)))

print()
print("=" * 78)
print("(2) THE MOMENT IDENTITY  E[P^m] = sum_i (-1)^i C(m,i)^2 (m-i)! mu_i")
print("=" * 78)
ok = True
for m in range(1, 9):
    lhs = E(expand(P4**m))
    rhs = sum((-1)**i * binomial(m, i)**2 * factorial(m - i) * factorial(i)
              for i in range(m + 1))
    if lhs != rhs:
        ok = False
    print("  m = %-3d  direct = %-6s   formula (mu_i = i!) = %-6s   agree: %s"
          % (m, lhs, rhs, lhs == rhs))
print("  identity holds : %s" % ok)
print()
print("  and the collapse:  C(m,i)^2 (m-i)! i! = m! C(m,i), so the sum is")
print("  m! * sum_i (-1)^i C(m,i) = 0.  Checked:")
for m in range(1, 6):
    coll = all(binomial(m, i)**2 * factorial(m - i) * factorial(i)
               == factorial(m) * binomial(m, i) for i in range(m + 1))
    print("     m = %-3d  C(m,i)^2 (m-i)! i! = m! C(m,i) for all i : %s" % (m, coll))

print()
print("=" * 78)
print("(3) THE MOMENTS OF Y ARE FORCED -- solve (*) = 0 for mu one at a time")
print("=" * 78)
mu = [Rational(1)]
for m in range(1, 9):
    s = sum((-1)**i * binomial(m, i)**2 * factorial(m - i) * mu[i]
            for i in range(m))
    coeff = (-1)**m * binomial(m, m)**2 * factorial(0)
    mu_m = Rational(-s, coeff)
    mu.append(mu_m)
    print("  m = %-3d  forces mu_%d = %-8s   (m! = %d)  match: %s"
          % (m, m, mu_m, factorial(m), mu_m == factorial(m)))
print("  -> mu_i = i! is the UNIQUE solution, so Y must be Exp(1).")

print()
print("=" * 78)
print("(4) COULD ONE REAL GAUSSIAN SUPPLY Y?  (would give n = 3)")
print("=" * 78)
print("  Y = c X^2 with X ~ N(0,1) has E[Y^i] = c^i (2i-1)!!:")
for c in (Rational(1), Rational(1, 2), Rational(2)):
    ms = [c**i * (factorial(2 * i) // (2**i * factorial(i))) for i in range(1, 5)]
    print("     c = %-5s moments %s   vs i! = [1, 2, 6, 24]   match: %s"
          % (c, ms, ms == [1, 2, 6, 24]))
print("  A one-Gaussian square has moments (2i-1)!! c^i, and (2i-1)!! != i! beyond")
print("  i = 1 for any c (3 vs 2 already at i = 2).  So Y needs a genuine Exp(1),")
print("  i.e. TWO real Gaussians, and with two more for Z the count n = 4 is FORCED")
print("  by this construction rather than chosen.  Beating n = 4 requires a different")
print("  mechanism, not a cheaper Y.")
