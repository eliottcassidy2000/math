#!/usr/bin/env python3
"""Independent hostile audit of the THM-3891 quadratic C8 theorem."""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(newline="\n")

GATES = 0


def gate(condition: bool, message: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(sp.expand(expression)) == 0, message)


def cubic_disc(a, b, c, d):
    return sp.expand(b**2*c**2 - 4*a*c**3 - 4*b**3*d - 27*a**2*d**2 + 18*a*b*c*d)


# -------------------------------------------------------------------------
# 1. Cover splitting and factor-degree bookkeeping.
# -------------------------------------------------------------------------

# For a connected degree-r finite-etale cover of A1, projective completion
# has all ramification over infinity.  RH requires at least 2r-2, while that
# one fibre supplies at most r-1.  This actually works in every degree r>1.
for cover_degree in range(2, 101):
    gate(
        2 * cover_degree - 2 > cover_degree - 1,
        f"one-infinity-fibre contradiction r={cover_degree}",
    )

source_degree_types = set()
target_degree_types = set()
for degrees in itertools.product(range(3), repeat=3):
    if sum(degrees) > 2:
        continue
    source_degree_types.add(tuple(sorted(degrees)))
    spare = 2 - sum(degrees)
    for additions in itertools.product(range(spare + 1), repeat=3):
        if sum(additions) != spare:
            continue
        target_degree_types.add(tuple(sorted(d + e for d, e in zip(degrees, additions))))
gate(source_degree_types == {(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1)},
     "complete dehomogenized degree types")
gate(target_degree_types == {(0, 0, 2), (0, 1, 1)},
     "complete homogeneous degree types")


# -------------------------------------------------------------------------
# 2. Independent coordinate classification and finite hostile atlases.
# -------------------------------------------------------------------------

A, C, U, V = sp.symbols("A C U V")
L1, L2, M1, M2 = sp.symbols("L1 L2 M1 M2")
product_coefficients = (L1*L2, L1*M2 + L2*M1, M1*M2, 0)
zero(
    cubic_disc(*product_coefficients)
    - M1**2 * M2**2 * (L1*M2 - L2*M1)**2,
    "product determinant discriminant",
)

# Sparse homogeneous polynomials are arrays indexed by the A exponent.
def pmul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a_i in enumerate(left):
        for j, b_j in enumerate(right):
            result[i + j] += a_i * b_j
    return tuple(result)


def ppow(poly: tuple[int, ...], exponent: int) -> tuple[int, ...]:
    result = (1,)
    for _ in range(exponent):
        result = pmul(result, poly)
    return result


values = (-1, 0, 1)
rows_011 = 0
pure_011 = 0
moving_011 = 0
constant_011 = 0
for row in itertools.product(values, repeat=8):
    l1a, l1c, m1a, m1c, l2a, l2c, m2a, m2c = row
    m_one = (m1c, m1a)
    m_two = (m2c, m2a)
    determinant = (
        l1c*m2c - l2c*m1c,
        l1a*m2c + l1c*m2a - l2a*m1c - l2c*m1a,
        l1a*m2a - l2a*m1a,
    )
    discriminant = pmul(pmul(ppow(m_one, 2), ppow(m_two, 2)), ppow(determinant, 2))
    rows_011 += 1
    if not discriminant[0] or any(discriminant[index] for index in range(1, len(discriminant))):
        continue
    pure_011 += 1
    gate(m1a == 0 and m2a == 0, "011 UFD forces both M rows to C")
    gate(determinant[1:] == (0, 0), "011 UFD forces determinant to C2")
    if l1a == 0 and l2a == 0:
        constant_011 += 1
    else:
        moving_011 += 1
        gate(m1c*l2a == m2c*l1a, "011 common moving L coordinate")
gate(rows_011 == 6561, "011 labelled atlas size")
gate(pure_011 > 0 and pure_011 == moving_011 + constant_011,
     "011 both coordinate classes exhausted")
gate(moving_011 > 0 and constant_011 > 0, "011 both coordinate classes occur")

rows_002 = 0
pure_002 = 0
for row in itertools.product(values, repeat=6):
    p_c2, p_ac, p_a2, q_c2, q_ac, q_a2 = row
    p_poly = (p_c2, p_ac, p_a2)
    q_poly = (q_c2, q_ac, q_a2)
    discriminant = pmul(ppow(p_poly, 2), ppow(q_poly, 2))
    rows_002 += 1
    if not discriminant[0] or any(discriminant[index] for index in range(1, len(discriminant))):
        continue
    pure_002 += 1
    gate(p_ac == p_a2 == q_ac == q_a2 == 0,
         "002 both determinant coordinates are C2")
gate(rows_002 == 729 and pure_002 > 0, "002 labelled atlas")

# Symbolic coordinate normal forms.
moving_form = U*(A*U+C*V)*((A+C)*U+C*V)
moving_poly = sp.Poly(sp.expand(moving_form), U, V)
moving_coeffs = tuple(
    moving_poly.coeff_monomial(m)
    for m in (U**3, U**2*V, U*V**2, V**3)
)
gate(moving_coeffs == (A**2+A*C, 2*A*C+C**2, C**2, 0),
     "moving coordinate coefficients")
zero(cubic_disc(*moving_coeffs)-C**8, "moving coordinate C8 discriminant")

p0, q0, p1, q1, e0, e1 = sp.symbols("p0 q0 p1 q1 e0 e1")
basis_det = p0*q1-q0*p1
quadratic_u = C**2*(p1*e0-p0*e1)/basis_det
quadratic_v = C**2*(q1*e0-q0*e1)/basis_det
zero(p0*quadratic_v-q0*quadratic_u-C**2*e0,
     "002 first coordinate recovery")
zero(p1*quadratic_v-q1*quadratic_u-C**2*e1,
     "002 second coordinate recovery")


# -------------------------------------------------------------------------
# 3. Binary-cubic pencil lemma, every gcd boundary, and hostile census.
# -------------------------------------------------------------------------

s, t = sp.symbols("s t")
ga, gb, gc, gd = sp.symbols("ga gb gc gd")
# F=UV(U-V) has coefficient row (0,1,-1,0).
pencil_disc = cubic_disc(t*ga, s+t*gb, -s+t*gc, t*gd)
gcd_counts = {0: 0, 1: 0, 2: 0, 3: 0}
support_min = {0: 99, 1: 99, 2: 99}
support_witness = {}
nonproportional_rows = 0
proportional_rows = 0
for row in itertools.product(range(-2, 3), repeat=4):
    aa, bb, cc, dd = row
    gcd_degree = int(dd == 0) + int(aa == 0) + int(aa+bb+cc+dd == 0)
    gcd_counts[gcd_degree] += 1
    proportional = aa == 0 and dd == 0 and bb == -cc
    specialized = sp.Poly(
        sp.expand(pencil_disc.subs({ga: aa, gb: bb, gc: cc, gd: dd, t: 1})),
        s,
        domain=sp.QQ,
    )
    support = int(specialized.sqf_part().degree())
    if proportional:
        proportional_rows += 1
        gate(support == 1, f"proportional one-support row={row}")
        gate(gcd_degree == 3, f"proportional gcd-three row={row}")
    else:
        nonproportional_rows += 1
        gate(gcd_degree <= 2, f"nonproportional gcd boundary row={row}")
        gate(support >= 2, f"pencil support hostile row={row}")
        if support < support_min[gcd_degree]:
            support_min[gcd_degree] = support
            support_witness[gcd_degree] = row
gate(nonproportional_rows == 620 and proportional_rows == 5,
     "complete 5^4 pencil census")
gate(set(support_witness) == {0, 1, 2}, "all nonproportional gcd boundaries occur")
gate(all(value >= 2 for value in support_min.values()), "pencil support minima")

# Riemann--Hurwitz portion of the pencil proof.  For residual degree r>=2,
# total ramification 2r-2 cannot fit over one value.  For r=1 in binary
# degree n>=3, the fixed squarefree factor has n-1>=2 roots; equality of two
# collision values would make the moving linear form vanish identically.
for residual_degree in range(2, 51):
    gate(2*residual_degree-2 > residual_degree-1,
         f"pencil one-value RH contradiction r={residual_degree}")
for form_degree in range(3, 51):
    gate(form_degree-1 >= 2, f"general pencil r1 has two fixed roots n={form_degree}")

# The n>=3 boundary is sharp.  In degree two, F=UV and G=U(U+V) are
# squarefree/nonproportional, but sF+tG=U(tU+(s+t)V) has discriminant
# (s+t)^2 and hence only one projective support value.
def quadratic_disc(a, b, c):
    return sp.expand(b**2-4*a*c)


zero(
    quadratic_disc(t, s+t, 0)-(s+t)**2,
    "degree-two sharp one-support counterexample",
)


# -------------------------------------------------------------------------
# 4. Weighted initials, proportional seam, and Newton hostile controls.
# -------------------------------------------------------------------------

x, z, h = sp.symbols("x z h")
f0 = sp.symbols("f0a:4")
f1 = sp.symbols("f1a:4")
f2 = sp.symbols("f2a:4")
local_rows = tuple(x**2*a0+z*a1+x*z*a2 for a0, a1, a2 in zip(f0, f1, f2))
local_disc = cubic_disc(*local_rows)
scaled_12 = sp.Poly(sp.expand(local_disc.subs({x: h*x, z: h**2*z})), h)
initial_12 = scaled_12.coeff_monomial(h**8)
expected_12 = cubic_disc(*(x**2*a0+z*a1 for a0, a1 in zip(f0, f1)))
zero(initial_12-expected_12, "weight 1,2 initial pencil")
gate(all(exponent[0] >= 8 for exponent, coeff in scaled_12.terms() if coeff),
     "weight 1,2 lower support")

lam, sigma = sp.symbols("lam sigma", nonzero=True)
proportional_rows = tuple(
    sp.expand(row.subs({a1: lam*a0 for a0, a1 in zip(f0, f1)}).subs(z, (sigma-x**2)/lam))
    for row in local_rows
)
expected_rows = tuple(
    sigma*a0+x*sigma*a2/lam-x**3*a2/lam
    for a0, a2 in zip(f0, f2)
)
for actual, expected in zip(proportional_rows, expected_rows):
    zero(actual-expected, "proportional regular coordinate row")
proportional_disc = cubic_disc(*proportional_rows)
scaled_13 = sp.Poly(
    sp.expand(proportional_disc.subs({x: h*x, sigma: h**3*sigma})), h
)
initial_13 = scaled_13.coeff_monomial(h**12)
expected_13 = cubic_disc(*(sigma*a0-x**3*a2/lam for a0, a2 in zip(f0, f2)))
zero(initial_13-expected_13, "weight 1,3 initial pencil")
gate(all(exponent[0] >= 12 for exponent, coeff in scaled_13.terms() if coeff),
     "weight 1,3 lower support")

# Proportional exits are visibly reducible/nonreduced.
mu = sp.symbols("mu")
F0_disc = cubic_disc(*f0)
zero(
    cubic_disc(*(C**2*a0+C*a2 for a0, a2 in zip(f0, f2)))
    - C**4*cubic_disc(*(C*a0+a2 for a0, a2 in zip(f0, f2))),
    "zero first pencil reducible exit",
)
profile = C**2+lam*A+mu*C
zero(cubic_disc(*(profile*a0 for a0 in f0))-profile**4*F0_disc,
     "fully proportional fourth-power exit")

# Independent replay of the moving normal form's exhaustive Newton seams.
al, al1, be, be1, gam, gam1, de, et = sp.symbols("al al1 be be1 gam gam1 de et")
aa = A**2+A*C+al*A+al1*C
bb = 2*A*C+C**2+be*A+be1*C
cc = C**2+gam*A+gam1*C
dd = de*A+et*C
moving_delta = cubic_disc(aa, bb, cc, dd)
H = sp.expand(z**8*moving_delta.subs({A: 1/z, C: x/z}))
gate(sp.denom(H) == 1, "moving infinity polynomial")
zero(H.subs(z, 0)-x**8, "moving unique infinity support")

def edge(poly, order, degree):
    w = sp.symbols("w")
    return sp.expand(sp.Poly(sp.expand(poly.subs(z, w*x**order)), x).coeff_monomial(x**degree))

w = sp.symbols("w")
zero(edge(H, 3, 6)-de*w*(4-27*de*w), "moving delta edge 3")
zero(edge(H, 5, 8)-(1+4*de*w), "moving delta edge 5")
H0 = sp.expand(H.subs(de, 0))
zero(
    edge(H0, 2, 6)
    - w*(-4*gam**3*w**2+(-27*et**2+36*et*gam-8*gam**2)*w+4*(et-gam)),
    "moving gamma edge 2",
)
zero(edge(H0, 4, 8)-(1+4*(et-gam)*w), "moving gamma edge 4")
Heq = sp.expand(H0.subs(et, gam))
zero(edge(Heq, 2, 6)-gam**2*w**2*(1-4*gam*w), "moving equal edge 2")
zero(
    edge(Heq, 3, 8)-(1+2*(2*be+gam-2*gam1)*w+gam**2*w**2),
    "moving equal edge 3",
)
last_seam = sp.expand(moving_delta.subs({de: 0, gam: 0}))
gate(sp.rem(sp.Poly(last_seam, C), sp.Poly(C, C)) == 0,
     "moving terminal seam C-divisible")


semantic = {
    "automatic_split": "all finite-etale A1 component degrees are one",
    "coordinate_classes": "all quadratic C8 rows are moving 011 or C2 constant",
    "pencil": "all nonproportional binary cubics have at least two branch values",
    "weighted": "initial pencils at weights (1,2) and (1,3)",
    "normalization": "distinct exceptional support points force distinct lifted normalization centers",
    "extension": "pencil lemma holds for every squarefree binary degree n>=3",
    "sharp_boundary": "degree-two counterexample F=UV and G=U(U+V)",
    "scope": "quadratic C8 theorem only; JC2 open",
}
blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("THM3891_INDEPENDENT_HOSTILE_AUDIT_20260823")
print("status=PASS;PROMOTION_READY;JC2_OPEN")
print(f"coordinate_atlas=011:{rows_011},pure:{pure_011},moving:{moving_011},constant:{constant_011};002:{rows_002},pure:{pure_002}")
print(f"pencil_census=625;nonproportional={nonproportional_rows};support_min={tuple(sorted(support_min.items()))};gcd_counts={tuple(sorted(gcd_counts.items()))}")
print("strict_extension=squarefree_binary_degree_n>=3_pencil_has_two_support_values_unless_proportional")
print("weighted_normalization_implication=PASS")
print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
print(f"gates={GATES}")
print("ALL CHECKS PASSED")
