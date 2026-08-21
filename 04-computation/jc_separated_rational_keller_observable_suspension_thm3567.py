#!/usr/bin/env python3
"""Exact companion for proved THM-3567's separated-family obstruction."""

import sympy as sp


failures = []
gates = 0


def gate(condition, label):
    global gates
    gates += 1
    if not bool(condition):
        failures.append(label)


x, y, T, Pvar, Avar, Bvar = sp.symbols("x y T P A B")

# A cubic with rational triple fibre, rational simple critical points, and
# distinct rational critical values.
f = sp.expand(x * (x - 3) * (x - 8))
d = sp.diff(f, x)
q = y / d

gate(f == x**3 - 11 * x**2 + 24 * x, "cubic expansion")
gate(sp.factor(d) == (x - 6) * (3 * x - 4), "critical factorization")
gate(sp.factor(sp.diff(f, x, 2).subs(x, 6)) != 0, "critical point six simple")
gate(
    sp.factor(sp.diff(f, x, 2).subs(x, sp.Rational(4, 3))) != 0,
    "critical point four-thirds simple",
)
gate(f.subs(x, 6) == -36, "first critical value")
gate(
    f.subs(x, sp.Rational(4, 3)) == sp.Rational(400, 27),
    "second critical value",
)

critical_value_polynomial = sp.expand((T + 36) * (27 * T - 400))
discriminant = sp.discriminant(T * (T - 3) * (T - 8) - Pvar, T)
gate(
    sp.expand(discriminant + critical_value_polynomial.subs(T, Pvar)) == 0,
    "cubic discriminant/critical-value polynomial",
)

delta_pullback = sp.factor(critical_value_polynomial.subs(T, f))
H = sp.factor(delta_pullback / d**2)
gate(H == (x + 1) * (3 * x - 25), "critical-value pullback quotient")
gate(sp.factor(delta_pullback - d**2 * H) == 0, "double critical divisor")
gate(
    sp.gcd(sp.Poly(d, x), sp.Poly(H, x)).degree() == 0,
    "critical and companion divisors disjoint",
)

# The separated rational plane map has unit Jacobian on d != 0.
jacobian = sp.factor(
    sp.diff(f, x) * sp.diff(q, y) - sp.diff(f, y) * sp.diff(q, x)
)
gate(jacobian == 1, "rational Keller identity")

triple_points = (0, 3, 8)
for point in triple_points:
    gate(f.subs(x, point) == 0, ("triple fibre P", point))
    gate(d.subs(x, point) != 0, ("triple fibre regularity", point))
    gate(q.subs({x: point, y: 0}) == 0, ("triple fibre Q", point))

# The two full-field polynomial observables of positive Q-degree.
A_pullback = sp.factor(delta_pullback * q)
B_pullback = sp.factor(delta_pullback * q**2)
gate(sp.factor(A_pullback - d * H * y) == 0, "A polynomial observable")
gate(sp.factor(B_pullback - H * y**2) == 0, "B polynomial observable")
gate(
    sp.factor(A_pullback**2 - delta_pullback * B_pullback) == 0,
    "quadratic suspension relation",
)

# Finite controls for the all-n divisibility lemma used in the paper proof.
parity_controls = 0
for n in range(1, 17):
    exponent = (n + 1) // 2
    sharp = sp.cancel(delta_pullback**exponent / d**n)
    _, sharp_denominator = sp.fraction(sharp)
    gate(sp.factor(sharp_denominator) == 1, ("sharp observable exponent", n))

    below = sp.cancel(delta_pullback ** (exponent - 1) / d**n)
    _, below_denominator = sp.fraction(below)
    gate(sp.factor(below_denominator) != 1, ("one-below pole survives", n))

    if n % 2 == 0:
        expected = H ** (n // 2)
    else:
        expected = d * H ** ((n + 1) // 2)
    gate(sp.factor(sharp - expected) == 0, ("parity normal form", n))
    parity_controls += 1

# The full-field intersection surface S: A^2=Delta(P)B has exactly two nodes.
Delta_P = critical_value_polynomial.subs(T, Pvar)
surface = sp.expand(Avar**2 - Delta_P * Bvar)
gradient = [sp.diff(surface, variable) for variable in (Pvar, Avar, Bvar)]
gate(
    sp.gcd(sp.Poly(Delta_P, Pvar), sp.Poly(sp.diff(Delta_P, Pvar), Pvar)).degree()
    == 0,
    "critical-value polynomial squarefree",
)

singular_points = (
    {Pvar: sp.Integer(-36), Avar: sp.Integer(0), Bvar: sp.Integer(0)},
    {Pvar: sp.Rational(400, 27), Avar: sp.Integer(0), Bvar: sp.Integer(0)},
)
for index, point in enumerate(singular_points):
    gate(surface.subs(point) == 0, ("singular point lies on S", index))
    gate(
        all(component.subs(point) == 0 for component in gradient),
        ("Jacobian singularity", index),
    )
    hessian = sp.Matrix(
        [
            [sp.diff(surface, u, v) for v in (Pvar, Avar, Bvar)]
            for u in (Pvar, Avar, Bvar)
        ]
    )
    gate(sp.factor(hessian.subs(point).det()) != 0, ("ordinary node", index))

# Each critical fibre has a regular companion line.  The full line collapses
# to the corresponding node, so the completion is not quasi-finite.
gate(sp.factor(f + 36) == (x - 6) ** 2 * (x + 1), "first critical fibre")
gate(
    sp.factor(27 * f - 400) == (3 * x - 25) * (3 * x - 4) ** 2,
    "second critical fibre",
)

companion_controls = (
    (sp.Integer(-1), sp.Integer(-36)),
    (sp.Rational(25, 3), sp.Rational(400, 27)),
)
for index, (source_x, target_P) in enumerate(companion_controls):
    gate(d.subs(x, source_x) != 0, ("companion root regular", index))
    gate(f.subs(x, source_x) == target_P, ("companion target value", index))
    gate(H.subs(x, source_x) == 0, ("companion H zero", index))
    gate(A_pullback.subs(x, source_x) == 0, ("companion A line collapse", index))
    gate(B_pullback.subs(x, source_x) == 0, ("companion B line collapse", index))

# At a regular target value, the cubic is separable and the open completion
# is the original degree-three finite-etale cover.
gate(
    sp.discriminant(f - Pvar, x).subs(Pvar, 0) == 14400,
    "regular triple-fibre discriminant",
)
gate(sp.degree(f, x) == 3, "generic degree three")

# Degree two is the sharp boundary for the companion-line claim.  With
# f=x^2 and the scaled Delta=4P, the completion is the quotient map
# (x,y)->(x^2,2xy,y^2).  It is etale over points on Delta=0 away from the
# node, so Delta!=0 is a sufficient etale open but not the maximal one.
f2 = x**2
d2 = 2 * x
q2 = y / d2
A2_pullback = sp.cancel(4 * f2 * q2)
B2_pullback = sp.cancel(4 * f2 * q2**2)
surface2 = Avar**2 - 4 * Pvar * Bvar
degree_two_source = {x: 0, y: 1}
degree_two_target = {Pvar: 0, Avar: 0, Bvar: 1}
gate(sp.factor(sp.diff(f2, x) * sp.diff(q2, y)) == 1, "degree-two Keller identity")
gate(sp.factor(A2_pullback**2 - 4 * f2 * B2_pullback) == 0, "degree-two relation")
gate(f2.subs(degree_two_source) == 0, "degree-two boundary lies on Delta zero")
gate(surface2.subs(degree_two_target) == 0, "degree-two boundary target lies on S")
gate(
    any(
        sp.diff(surface2, variable).subs(degree_two_target) != 0
        for variable in (Pvar, Avar, Bvar)
    ),
    "degree-two boundary target smooth",
)
gate(
    sp.det(
        sp.Matrix(
            [
                [sp.diff(A2_pullback, x), sp.diff(A2_pullback, y)],
                [sp.diff(B2_pullback, x), sp.diff(B2_pullback, y)],
            ]
        )
    ).subs(degree_two_source)
    != 0,
    "degree-two boundary map etale",
)
gate(
    all(
        sp.diff(surface2, variable).subs({Pvar: 0, Avar: 0, Bvar: 0}) == 0
        for variable in (Pvar, Avar, Bvar)
    ),
    "degree-two origin remains nodal",
)

if failures:
    raise RuntimeError(("THM-3567 exact companion failures", failures[:8], len(failures)))

print("THM-3567 separated rational Keller full-field companion: PASS")
print("cubic seed: f=x(x-3)(x-8); Q=y/f'(x); Jac(f,Q)=1")
print("triple collision: (0,0),(3,0),(8,0) -> (0,0)")
print("critical values: -36, 400/27; both simple and distinct")
print("full-field intersection surface: A^2=(P+36)(27P-400)B")
print("singular locus: exactly (-36,0,0) and (400/27,0,0)")
print("collapsed companion lines: x=-1 and x=25/3")
print("parity divisibility controls: n=1..16, sharp ceil(n/2)")
print("parity controls:", parity_controls)
print("exact gates:", gates)
print("degree-two boundary: Delta=0 has etale points away from the node")
print("consequence: this separated completion is singular; d>=3 is not quasi-finite")
print("scope: full target field of separated seeds only; no JC(2) conclusion")
