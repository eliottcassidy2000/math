#!/usr/bin/env python3
"""gmc4_family_kps_S128c117.py -- kind-pasteur-2026-07-20-S128c117

FOLLOW-UP: three claims from the previous run need hardening, and one of them is
probably false.

(1) THE ONE-PARAMETER FAMILY.  The scan found E[P^m] = 0 not only at the published
    point but along the whole DIAGONAL a = b of
        P_a = (1 + a Z2)(conj(Z1)(1 - a Z1) + conj(Z2)),
    with Q-moments 1, 2a, 6a^2, 24a^3 at a = 1, 2, 1/2, -1, -2 ... i.e. apparently
        E[Q P_a^m] = a^{m-1} m!,   so   E[Q exp(t P_a)] = t/(1 - a t).
    If that holds, the pole is at t = 1/a and can be placed at ANY non-zero real --
    the failure is not tied to t = 1.  Tested here to higher m and more a.
    (Note the previous run's prose said 1/(1-t); the data say the m = 0 term is 0, so
    the correct generating function is t/(1-t).  Corrected.)

(2) TERM MINIMALITY.  Two of the six terms could be dropped while KEEPING E[P^m] = 0.
    But that is only half of what a counterexample needs -- the Q-moments must stay
    non-zero.  Checked here.

(3) THE "QUADRATIC COUNTEREXAMPLES" ARE ALMOST CERTAINLY ARTEFACTS.  The previous
    search filtered on E[P^m] = 0 for m <= 3 only, which is a very weak sieve: it is
    three linear-ish conditions on a large coefficient space, so hits are expected by
    chance.  Re-tested to m = 8.  If they die, the honest conclusion is that the sieve
    was too shallow, NOT that quadratics work -- and the depth at which they die is the
    number to report.
"""
import sys
from math import factorial
from sympy import symbols, expand, Poly, Rational

MM = int(sys.argv[1]) if len(sys.argv) > 1 else 8
z1, z1c, z2, z2c = symbols('z1 z1c z2 z2c')
V = (z1, z1c, z2, z2c)


def E(expr):
    p = Poly(expand(expr), *V)
    tot = 0
    for (a, b, c, d), co in p.terms():
        if a == b and c == d:
            tot += co * factorial(a) * factorial(c)
    return tot


Q = z2
print("=" * 78)
print("(1) THE ONE-PARAMETER FAMILY  P_a = (1 + a Z2)(conj(Z1)(1 - a Z1) + conj(Z2))")
print("=" * 78)
print("  %-8s %-8s %s" % ("a", "E[P^m]", "E[Q P^m] vs a^{m-1} m!"))
for a in (Rational(1), Rational(2), Rational(1, 2), Rational(-1), Rational(3),
          Rational(-3, 2)):
    Pa = (1 + a * z2) * (z1c * (1 - a * z1) + z2c)
    zeros, match = True, True
    obs, pred = [], []
    for m in range(1, MM + 1):
        Pm = expand(Pa**m)
        if E(Pm) != 0:
            zeros = False
        o = E(expand(Q * Pm))
        p_ = a**(m - 1) * factorial(m)
        obs.append(o)
        pred.append(p_)
        if o != p_:
            match = False
    print("  a = %-6s E[P^m]=0 all m<=%d : %-5s   E[QP^m] = a^{m-1} m! : %s"
          % (a, MM, zeros, match))
    print("           observed  %s" % obs[:5])
    print("           predicted %s" % pred[:5])
print()
print("  => E[Q exp(t P_a)] = sum_{m>=1} a^{m-1} t^m = t/(1 - a t):")
print("     a SIMPLE POLE AT t = 1/a, placeable anywhere non-zero.")

print()
print("=" * 78)
print("(2) TERM MINIMALITY -- dropping a term must keep BOTH properties")
print("=" * 78)
P = (1 + z2) * (z1c * (1 - z1) + z2c)
terms = Poly(expand(P), *V).terms()


def build(sub):
    e = 0
    for mono, co in sub:
        t = co
        for v, k in zip(V, mono):
            t *= v**k
        e += t
    return e


for i in range(len(terms)):
    sub = [terms[j] for j in range(len(terms)) if j != i]
    Pd = build(sub)
    zs = [E(expand(Pd**m)) for m in range(1, MM + 1)]
    qs = [E(expand(Q * expand(Pd**m))) for m in range(1, 5)]
    allz = all(v == 0 for v in zs)
    anyq = any(v != 0 for v in qs)
    print("  drop term %d (%s) : E[P^m]=0 all m<=%d : %-5s ; some E[QP^m]!=0 : %-5s%s"
          % (i, build([terms[i]]), MM, allz, anyq,
             "   <-- STILL A COUNTEREXAMPLE" if (allz and anyq) else ""))
print("  -> a 5-term counterexample exists : %s"
      % any(all(E(expand(build([terms[j] for j in range(len(terms)) if j != i])**m)) == 0
                for m in range(1, MM + 1))
            and any(E(expand(Q * expand(build([terms[j] for j in range(len(terms))
                                               if j != i])**m))) != 0
                    for m in range(1, 5))
            for i in range(len(terms))))

print()
print("=" * 78)
print("(3) THE 'QUADRATIC' HITS -- were they real, or a too-shallow sieve?")
print("=" * 78)
cands = [
    z1 * z1c + z1 * z2 + z1 * z2c - z1c**2 - z1c * z2 - z1c * z2c - z1c
    - z2 * z2c - z2 - z2c**2 - z2c,
    z1 * z1c + z1 * z2c + z1 - z1c * z2 - z1c * z2c - z1c - z2 * z2c - z2
    - z2c**2 - z2c,
    z1**2 + z1 * z1c + z1 * z2 + z1 - z1c * z2c - z1c - z2 * z2c - z2
    - z2c**2 - z2c,
]
for idx, e in enumerate(cands, 1):
    first_fail = None
    for m in range(1, MM + 1):
        if E(expand(e**m)) != 0:
            first_fail = m
            break
    print("  candidate %d : E[P^m] = 0 fails first at m = %s"
          % (idx, first_fail if first_fail else "never (up to %d)" % MM))
print()
print("  If each dies at m = 4 or 5, the m <= 3 sieve was simply too shallow -- three")
print("  conditions on a coefficient space that large will produce hits by chance.")
print("  The honest reading is then: NO quadratic counterexample was found, and the")
print("  earlier hits are withdrawn.  Reporting the death depth is the useful part,")
print("  because it says how deep a sieve has to be to mean anything here.")
