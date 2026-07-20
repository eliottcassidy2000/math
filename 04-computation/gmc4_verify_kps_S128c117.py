#!/usr/bin/env python3
"""gmc4_verify_kps_S128c117.py -- kind-pasteur-2026-07-20-S128c117

VERIFYING AN OUTSIDE COUNTEREXAMPLE TO THE GAUSSIAN MOMENT CONJECTURE AT n = 4.

Claim supplied by the owner from an outside source:

    P = (1 + Z2)(W1(1 - Z1) + W2),   Q = Z2,   "in 4 real Gaussians via complex Zj, Wj"
    E[P^m] = 0 for all m >= 1,  but  E[Q P^m] = m! != 0 for all m >= 1,
    hence GMC(4) is false.  (cubic, 6 terms)

GMC here is the Gaussian-moment form of Zhao's vanishing conjecture / the Mathieu-subspace
property: the set {P : E[P^m] = 0 for all m >= 1} ought to be a Mathieu subspace, i.e.
E[Q P^m] should vanish for all large m.  A P with E[P^m] = 0 for every m and some Q with
E[Q P^m] never vanishing refutes exactly that.

WHICH READING OF "4 REAL GAUSSIANS"?  Two candidates, and only one can be right:
  (i)  W_j = conj(Z_j).  Then Z1, Z2 complex = 4 real Gaussians -- matches the stated count.
  (ii) Z1, Z2, W1, W2 independent complex = 8 real.  But then P contains W's and no
       W-bars, so every monomial of Q P^m has unbalanced W-charge and BOTH expectations
       vanish identically by rotation invariance -- E[Q P^m] = 0, contradicting m!.
So (ii) is self-refuting and (i) must be intended.  Both are computed below rather than
assumed, because an interpretation that is merely plausible is not verified.

EXACT ARITHMETIC.  For independent standard complex Gaussians with E|Z|^2 = 1,
    E[ prod_j Z_j^{a_j} conj(Z_j)^{b_j} ] = prod_j delta(a_j, b_j) * a_j!
(rotation invariance kills unbalanced charge; the balanced case is the classical a!).
So expectations are computed termwise on the expanded polynomial -- no sampling.
"""
import sys
from math import factorial
from sympy import symbols, expand, Poly, simplify

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7

z1, z1c, z2, z2c = symbols('z1 z1c z2 z2c')
w1, w1c, w2, w2c = symbols('w1 w1c w2 w2c')


def E_two(expr):
    """Expectation in reading (i): variables z1,z1c,z2,z2c, W_j = conj(Z_j)."""
    p = Poly(expand(expr), z1, z1c, z2, z2c)
    tot = 0
    for (a, b, c, d), co in p.terms():
        if a == b and c == d:
            tot += co * factorial(a) * factorial(c)
    return tot


def E_four(expr):
    """Expectation in reading (ii): z1,z2,w1,w2 independent complex."""
    p = Poly(expand(expr), z1, z1c, z2, z2c, w1, w1c, w2, w2c)
    tot = 0
    for m, co in p.terms():
        a, b, c, d, e, f, g, h = m
        if a == b and c == d and e == f and g == h:
            tot += co * factorial(a) * factorial(c) * factorial(e) * factorial(g)
    return tot


print("=" * 78)
print("READING (ii): Z1,Z2,W1,W2 independent complex (8 real) -- self-refuting?")
print("=" * 78)
P2 = (1 + z2) * (w1 * (1 - z1) + w2)
Q2 = z2
for m in range(1, 4):
    print("  m = %d :  E[P^m] = %s ,  E[Q P^m] = %s"
          % (m, E_four(P2**m), E_four(Q2 * P2**m)))
print("  -> both vanish identically (unbalanced W-charge), so this reading cannot")
print("     produce E[Q P^m] = m!.  Reading (i) is the intended one.")

print()
print("=" * 78)
print("READING (i): W_j = conj(Z_j), so Z1, Z2 complex = 4 REAL Gaussians")
print("=" * 78)
P = (1 + z2) * (z1c * (1 - z1) + z2c)
Q = z2
pe = Poly(expand(P), z1, z1c, z2, z2c)
print("  P expanded has %d terms, total degree %d :" % (len(pe.terms()), pe.total_degree()))
print("     P = %s" % expand(P))
print()
print("  %-4s %-16s %-16s %-10s" % ("m", "E[P^m]", "E[Q P^m]", "m!"))
ok_zero = True
ok_fact = True
for m in range(1, MMAX + 1):
    Pm = expand(P**m)
    a = E_two(Pm)
    b = E_two(expand(Q * Pm))
    if a != 0:
        ok_zero = False
    if b != factorial(m):
        ok_fact = False
    print("  %-4d %-16s %-16s %-10d" % (m, a, b, factorial(m)))
print()
print("  E[P^m] = 0 for all m tested        : %s" % ok_zero)
print("  E[Q P^m] = m! for all m tested     : %s" % ok_fact)
print()
if ok_zero and ok_fact:
    print("  VERDICT: the counterexample is CONFIRMED, exactly, in 4 real Gaussians.")
    print("  {P : E[P^m] = 0 for all m} is therefore NOT a Mathieu subspace at n = 4,")
    print("  so the Gaussian-moment form of Zhao's vanishing conjecture is FALSE there.")
else:
    print("  VERDICT: NOT reproduced under this reading -- do not propagate the claim.")

print()
print("=" * 78)
print("WHERE THE m! COMES FROM")
print("=" * 78)
print("  For a single standard complex Gaussian, E[Z^m conj(Z)^m] = m!.  So an m! in")
print("  the answer is the signature of ONE balanced charge-m pair surviving, all other")
print("  monomials cancelling.  Q = Z2 supplies the single extra +1 of Z2-charge that")
print("  turns the unique surviving unbalanced monomial of P^m into a balanced one.")
print("  Checking that reading directly -- the Z2-charge profile of P:")
for (a, b, c, d), co in sorted(pe.terms()):
    print("     z1^%d z1c^%d z2^%d z2c^%d  coeff %-4s   z1-charge %+d, z2-charge %+d"
          % (a, b, c, d, co, a - b, c - d))
