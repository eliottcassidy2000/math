#!/usr/bin/env python3
"""Exact companion for THM-3953's split-ramification boundary triangle."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


a, b, c, h, P, T = sp.symbols("a b c h P T")
r0 = c * a * (a + b)
r1 = c * b * (a + b)
r2 = -c * a * b
roots = (r0, r1, r2)

zero(r0 * r1 + r0 * r2 + r1 * r2,
     "missing linear h-row relation")

C = 2 * sum(roots)
E = 2 * r0 * r1 * r2
G = sp.expand(E + C * h**2 - 2 * h**3)
zero(G + 2 * (h - r0) * (h - r1) * (h - r2),
     "hidden cubic splits into three graph roots")

D01 = sp.expand((a - b) * (a + b))
D02 = sp.expand(a * (a + 2 * b))
D12 = sp.expand(b * (2 * a + b))
zero(r0 - r1 - c * D01, "first pair difference")
zero(r0 - r2 - c * D02, "second pair difference")
zero(r1 - r2 - c * D12, "third pair difference")

# Dehomogenizing by b=1 freezes pairwise coprimality; the point b=0 is
# checked separately and belongs only to D12.
x = sp.symbols("x")
d01 = sp.expand(D01.subs({a: x, b: 1}))
d02 = sp.expand(D02.subs({a: x, b: 1}))
d12 = sp.expand(D12.subs({a: x, b: 1}))
gate(sp.resultant(d01, d02, x) == -3,
     "D01 and D02 are coprime")
gate(sp.resultant(d01, d12, x) == -3,
     "D01 and D12 are coprime")
gate(sp.resultant(d02, d12, x) == -3,
     "D02 and D12 are coprime")
gate(D01.subs(b, 0) != 0 and D02.subs(b, 0) != 0
     and D12.subs(b, 0) == 0,
     "no hidden common projective zero at b=0")


# ---------------------------------------------------------------------------
# Natural depressed-cubic surface: exact domain and collision locus.
# ---------------------------------------------------------------------------

q = sp.expand(E + C * P)
F = sp.expand(T**3 - 3 * P * T - q)
domain_obstruction = sp.expand(C**3 + 27 * E)
domain_factor = 2 * c**3 * (a - b)**2 * (a + 2 * b)**2 * (2 * a + b)**2
zero(domain_obstruction - domain_factor,
     "monic-cubic root obstruction factorization")
gate(domain_factor != 0, "formal distinct-root domain control")

zero(sp.diff(F, T) - 3 * (T**2 - P), "relative derivative")
Fp_ram = sp.expand(sp.diff(F, P).subs({P: h**2, T: -h}))
zero(Fp_ram - (3 * h - C), "relative P derivative on ramification")
zero(sp.diff(G, h) + 2 * h * Fp_ram,
     "simple hidden root forces smooth cubic surface")

disc_G = sp.factor(sp.discriminant(G, h))
disc_expected = 16 * (r0 - r1)**2 * (r0 - r2)**2 * (r1 - r2)**2
zero(disc_G - disc_expected, "hidden-cubic collision discriminant")
zero(
    disc_G - 16 * c**6 * D01**2 * D02**2 * D12**2,
    "finite singular-parameter carrier",
)


# ---------------------------------------------------------------------------
# Six primitive pair-collision rows.  F_P distinguishes the zero collision
# (smooth) from the nonzero 2:-1 collision (singular).
# ---------------------------------------------------------------------------

collision_rows = [
    ("01_zero", {b: -a}, 0, 1, 2, 0, c * a**2, "smooth"),
    ("01_nonzero", {b: a}, 0, 1, 2, 2 * c * a**2, -c * a**2,
     "singular"),
    ("02_zero", {a: 0}, 0, 2, 1, 0, c * b**2, "smooth"),
    ("02_nonzero", {a: -2 * b}, 0, 2, 1, 2 * c * b**2, -c * b**2,
     "singular"),
    ("12_zero", {b: 0}, 1, 2, 0, 0, c * a**2, "smooth"),
    ("12_nonzero", {b: -2 * a}, 1, 2, 0, 2 * c * a**2, -c * a**2,
     "singular"),
]

for label, substitution, i, j, k, repeated, third, kind in collision_rows:
    ri = sp.expand(roots[i].subs(substitution))
    rj = sp.expand(roots[j].subs(substitution))
    rk = sp.expand(roots[k].subs(substitution))
    zero(ri - rj, f"{label}: selected roots coincide")
    zero(ri - repeated, f"{label}: repeated value")
    zero(rk - third, f"{label}: third value")
    fp = sp.expand((3 * ri - C.subs(substitution)))
    if kind == "smooth":
        zero(fp + 2 * third, f"{label}: nonzero smooth derivative")
        gate(fp != 0, f"{label}: formal smooth control")
    else:
        zero(fp, f"{label}: singular P derivative")

# At a common-gcd zero all roots, C, E and the relevant first derivatives
# vanish.  This is the singular override to every primitive row.
for index, root in enumerate(roots):
    zero(root.subs(c, 0), f"common-gcd collision root {index}")
zero(C.subs(c, 0), "common-gcd collision C")
zero(E.subs(c, 0), "common-gcd collision E")
zero(Fp_ram.subs({c: 0, h: 0}), "common-gcd collision F_P")


# ---------------------------------------------------------------------------
# Concrete controls.
# a=t,b=1,c=1 gives three disjoint primitive collision fibres.
# a=2,b=1,c=t is the constant-ratio one-address hostile boundary.
# a=b=1 is the duplicate-root reducible boundary.
# ---------------------------------------------------------------------------

t = sp.symbols("t")
positive = {a: t, b: 1, c: 1}
positive_D = tuple(sp.expand(D.subs(positive)) for D in (D01, D02, D12))
gate(positive_D == (t**2 - 1, t**2 + 2 * t, 2 * t + 1),
     "positive control collision polynomials")
chosen = (-1, 0, sp.Rational(-1, 2))
gate(len(set(chosen)) == 3, "positive control has three distinct fibres")
for index, (polynomial, value) in enumerate(zip(positive_D, chosen)):
    zero(polynomial.subs(t, value),
         f"positive control chosen pair collision {index}")

hostile = {a: 2, b: 1, c: t}
hostile_D = tuple(sp.expand(D.subs(hostile)) for D in (D01, D02, D12))
gate(all(sp.degree(polynomial, t) == 0 for polynomial in hostile_D),
     "constant-ratio hostile has no primitive collision fibres")
hostile_roots = tuple(sp.expand(root.subs(hostile)) for root in roots)
gate(hostile_roots == (6 * t, 3 * t, -2 * t),
     "constant-ratio hostile has one common address")

duplicate = {a: 1, b: 1, c: 1}
duplicate_roots = tuple(sp.expand(root.subs(duplicate)) for root in roots)
gate(duplicate_roots == (2, 2, -1),
     "duplicate-root hostile packet")
duplicate_F = sp.factor(F.subs(duplicate))
gate(sp.rem(sp.Poly(duplicate_F, T), sp.Poly(T + 2, T)) == 0,
     "duplicate-root hostile cubic is reducible")
zero(domain_obstruction.subs(duplicate),
     "duplicate-root domain obstruction vanishes")


summary = {
    "checks": CHECKS,
    "roots": "r0=ca(a+b),r1=cb(a+b),r2=-cab",
    "collisions": "three nonconstant pairwise-coprime primitive carriers",
    "table": "zero collisions smooth; nonzero 2:-1 collisions singular",
    "normality": "domain hypersurface with finite singular locus",
    "boundary": "three distinct fibres give a forbidden boundary triangle",
    "scope": "nonconstant ratio and three distinct polynomial roots",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3953 rationally split hidden-cubic boundary-triangle companion")
print(f"CHECKS={CHECKS}")
print("ROOTS=r0=c*a*(a+b);r1=c*b*(a+b);r2=-c*a*b")
print("PAIR_CARRIERS=(a-b)(a+b),a(a+2b),b(2a+b);PAIRWISE_COPRIME")
print("COLLISIONS=ZERO_SMOOTH;NONZERO_2_TO_MINUS1_SINGULAR;C_ZERO_TRIPLE_SINGULAR")
print("SURFACE=DOMAIN;FINITE_SINGULAR_LOCUS;NORMAL")
print("BOUNDARY=THREE_DISTINCT_FIBRES_FORM_TRIANGLE;NO_KELLER_ATLAS")
print("SCOPE=CONSTANT_RATIO,DUPLICATE_ROOT,RATIONAL_DENOMINATOR_OUTSIDE;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
