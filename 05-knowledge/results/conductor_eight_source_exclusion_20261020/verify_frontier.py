#!/usr/bin/env python3
"""Exact corroboration of the universal conductor-eight arguments.

Only the declared bounded numerical-semigroup/characteristic universes are
enumerated. The symbolic identities are universal. No search for polynomial
Keller maps is used as a proof.
"""
import itertools
import json
import math
import sympy as S

checks = []


def check(label, condition):
    if not condition:
        raise RuntimeError("FAILED: " + label)
    checks.append(label)


def zero(label, expr):
    check(label, S.expand(expr) == 0)


def generated(gens, bound):
    vals = {0}
    for n in range(1, bound + 1):
        if any(n >= g and n - g in vals for g in gens):
            vals.add(n)
    return vals


def primitive(gens):
    vals = generated(gens, max(gens))
    return tuple(n for n in range(1, max(gens) + 1)
                 if n in vals and not any(i in vals and n-i in vals
                                          for i in range(1, n)))


targets = [(2, 9), (3, 5), (4, 5, 6)]
closed = []
for subset in itertools.combinations(range(1, 7), 3):
    vals = {0, *subset}
    if all((n in vals) != (7-n in vals) for n in range(8)) and all(
            a+b >= 8 or a+b in vals for a in vals for b in vals):
        closed.append(tuple(n for n in range(8) if n not in vals))
check("all symmetric genus-four gap sets",
      sorted(closed) == [(1, 2, 3, 7), (1, 2, 4, 7), (1, 3, 5, 7)])

patterns = {}
for gens in targets:
    allowed = sorted(generated(gens, 9) - {0})
    found = []

    def extend(seq):
        d = math.gcd(*seq)
        if d == 1:
            if primitive(seq) == gens:
                found.append(tuple(seq))
            return
        # Classical bound: every later entry <= conductor-1 = 7.
        old = math.gcd(*seq[:-1])
        previous_n = old // d
        for nxt in [a for a in allowed if a <= 7]:
            new_d = math.gcd(d, nxt)
            if new_d == d:
                continue
            n = d // new_d
            if n*nxt not in generated(seq, n*nxt):
                continue
            if nxt >= seq[-1]*previous_n:
                continue
            extend(seq + [nxt])

    for a, b in itertools.permutations(allowed, 2):
        if a % b and b % a:
            extend([a, b])
    patterns[str(gens)] = sorted(found)

expected = {
    str((2, 9)): [(2, 9), (6, 9, 2), (9, 2), (9, 6, 2)],
    str((3, 5)): [(3, 5), (5, 3)],
    str((4, 5, 6)): [(4, 6, 5), (6, 4, 5)],
}
check("all conductor-eight characteristic patterns", patterns == expected)
for gens in targets:
    vals = generated(gens, 40)
    check("target conductor and genus " + str(gens),
          [n for n in range(40) if n not in vals] == list(
              next(g for g in closed if all(n not in vals for n in g))))
    check("target final conductor " + str(gens),
          7 not in vals and all(n in vals for n in range(8, 41)))
    check("target genus four " + str(gens),
          len([n for n in range(40) if n not in vals]) == 4)

t, x = S.symbols("t x")
a, b, c, e, d, f = S.symbols("a b c e d f")
A, B, L = S.symbols("A B L")
ap, bp, cp, ep, dp, fp = S.symbols("ap bp cp ep dp fp")
p0 = t**4+a*t**2+b*t
qgen = t**6+A*t**5+c*t**4+B*t**3+e*t**2+L*t
zero("wrong t5 remainder degree11",
     S.Poly(qgen**2-p0**3, t).nth(11)-2*A)
q1 = qgen.subs(A, 0)
rho = 2*c-3*a
r1 = S.Poly(q1**2-p0**3-rho*p0*q1, t)
zero("wrong t3 remainder degree9", r1.nth(9)-(2*B-3*b))
q2 = q1.subs(B, 3*b/S.Integer(2))
sigma = a*c-c*c+2*e
r2 = S.Poly(q2**2-p0**3-rho*p0*q2-sigma*p0**2, t)
zero("wrong t coefficient remainder degree7",
     r2.nth(7)-(3*a*b/S.Integer(2)-2*b*c+2*L))
q0 = q2.subs(L, b*(c-3*a/S.Integer(4)))
r3 = S.Poly(q0**2-p0**3-rho*p0*q0-sigma*p0**2, t)
kappa = -a**3+a*a*c-a*e-3*b*b/S.Integer(4)
zero("degree6 remainder coefficient", r3.nth(6)-kappa)
r4 = S.Poly(r3.as_expr()-kappa*q0, t)
D = e-a*c+3*a*a/S.Integer(4)
zero("degree5 remainder coefficient", r4.nth(5)+b*D)
r5 = S.Poly(r4.as_expr()-r4.nth(4)*p0, t)
zero("delta-three boundary recovery",
     r5.as_expr().subs(e, a*c-3*a*a/4)-b**3*t*(3*a+2*t*t)/16)
check("b-zero boundary is even", all(
    S.Poly(z.subs(b, 0), t).nth(n) == 0
    for z in [p0, q0] for n in range(1, 7, 2)))

for odd, delta in [(11, 7), (9, 6), (7, 5), (5, 4), (3, 3)]:
    degrees = [0, 6, odd, 6+odd]
    check("basis residues for odd remainder " + str(odd),
          sorted(n % 4 for n in degrees) == [0, 1, 2, 3])
    check("exact codimension for odd remainder " + str(odd),
          sum(n-n % 4 for n in degrees)//4 == delta)
u, v = S.symbols("u v")
basis_transform = S.Matrix([[1, 0, -v, 0], [0, 1, -u, -v],
                            [0, 0, 1, -u], [0, 0, 0, 1]])
check("unimodular module change", basis_transform.det() == 1)

p, q = p0+d, q0+f
vars_, derivatives = [a, b, c, e, d, f], [ap, bp, cp, ep, dp, fp]
Dx = lambda z: sum(S.diff(z, v)*dv for v, dv in zip(vars_, derivatives))
J = S.Poly(Dx(p)*S.diff(q, t)-S.diff(p, t)*Dx(q), t)
zero("full constant Jacobian coefficient",
     J.nth(0)-b*((c-3*a/4)*dp-fp))
equations = {bp: 0, cp: 3*ap/2, dp: a*ap/2, ep: c*ap, fp: e*ap/2}
for n in [7, 5, 3, 2]:
    zero("solved full Jacobian coefficient " + str(n), J.nth(n).subs(equations))
zero("full final Jacobian identity",
     J.as_expr().subs(equations)+ap*b*(3*a*a-4*a*c+6*b*t+4*e)/8)
zero("remaining t coefficient", J.nth(1).subs(equations)+3*ap*b*b/4)

realizations = [
    ("quadratic-nonic", t*t, t**9+t, (2, 9)),
    ("cubic-quintic", t**3, t**5+t, (3, 5)),
    ("quartic-sextic", t**4+t, t**6+S.Rational(3, 2)*t**3+t*t, (4, 5, 6)),
    ("sextic-nonic", t**6+4*t*t, t**9+6*t**5+6*t, (2, 9)),
]
for name, P, Q, gens in realizations:
    check(name + " immersion", S.gcd(S.diff(P, t), S.diff(Q, t)) == 1)
P, Q = realizations[2][1:3]
rem = Q*Q-P**3-2*P*P+S.Rational(3, 4)*Q-P
zero("quartic sextic explicit remainder",
     rem-(-t**5+t**3/8-5*t*t/4-t))
check("quartic sextic remainder degree", S.degree(rem, t) == 5)
P, Q = realizations[3][1:3]
zero("sextic nonic quadratic recovery", Q*Q-P**3-8*P-4*t*t)
z = S.symbols("z")
g = z**4+6*z*z+6
check("sextic nonic four distinct nonzero node squares",
      g.subs(z, 0) != 0 and S.gcd(g, S.diff(g, z)) == 1)

# Explicit adversarial witnesses for each higher-delta odd remainder.
hostiles = [
    ("delta7 wrong t5", t**4, t**6+t**5, 11, 7),
    ("delta6 wrong t3", t**4, t**6+t**3, 9, 6),
    ("delta5 wrong t1", t**4+t, t**6+S.Rational(3, 2)*t**3+t, 7, 5),
]
for name, P, Q, n, delta in hostiles:
    RR = S.Poly(Q*Q-P**3, t)
    if RR.degree() % 4 == 0:
        RR = S.Poly(RR.as_expr()-RR.LC()*P**(RR.degree()//4), t)
    check(name + " exact odd remainder", RR.degree() == n)
    check(name + " primitive fraction field witness", math.gcd(4, n) == 1)
P, Q = t**4+t, t**6+S.Rational(3, 2)*t**3
zero("delta3 adjacent hostile",
     Q*Q-P**3+S.Rational(3, 4)*Q-t**3/8)

# Genuine positive Keller controls for the quadratic lemma; arbitrary
# polynomials F(p) are permitted and removed by constant target shears.
for F in [lambda p: 0, lambda p: p**3+2*p, lambda p: p**5-p*p+7]:
    P = t*t+2*x+3
    Q = 3*t+F(P)
    JJ = S.diff(P, x)*S.diff(Q, t)-S.diff(P, t)*S.diff(Q, x)
    zero("genuine quadratic Keller positive control " + str(S.degree(Q, t)), JJ-6)
    zero("quadratic target reduction " + str(S.degree(Q, t)), Q-F(P)-3*t)

result = {
    "status": "PASS",
    "count": len(checks),
    "symmetric_gap_sets": sorted(closed),
    "characteristic_patterns": patterns,
    "universal_claims": [
        "genus-four semigroup and nonprincipal-pattern classification",
        "all listed patterns have explicit immersed polynomial realizations",
        "generic <4,5,6> monic Keller exclusion",
        "actual monic quadratic Keller coordinate is an automorphism",
    ],
    "not_claimed": [
        "full conductor-eight Keller exclusion",
        "generic quadratic generator is itself a Keller coordinate",
        "polynomial realizability from abstract characteristic admissibility",
    ],
    "checks": checks,
}
print(json.dumps(result, indent=2, sort_keys=True))
