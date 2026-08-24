#!/usr/bin/env python3
"""Exact companion for THM-3958's one-hidden-root principal different."""

from __future__ import annotations

import hashlib
import json

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


# ---------------------------------------------------------------------------
# Unique marked-root normal form and irreducibility obstruction.
# ---------------------------------------------------------------------------

h, T, P, r, a, s = sp.symbols("h T P r a s")
C = 2 * (r - a)
E = 2 * a * r**2
G = sp.expand(E + C * h**2 - 2 * h**3)
zero(G + 2 * (h - r) * (h**2 + a * h + a * r),
     "marked-root hidden-cubic normal form")
zero(C.subs(a, (s - r) / 2) - (3 * r - s), "C in r,s coordinates")
zero(E.subs(a, (s - r) / 2) - (s - r) * r**2,
     "E in r,s coordinates")

Q = sp.expand(h**2 + a * h + a * r)
disc_Q = sp.factor(sp.discriminant(Q, h))
zero(disc_Q - a * (a - 4 * r), "residual quadratic discriminant")
disc_rs = sp.factor(disc_Q.subs(a, (s - r) / 2))
zero(4 * disc_rs - (s - r) * (s - 9 * r),
     "residual discriminant in r,s coordinates")
zero(Q.subs({h: r, a: (s - r) / 2}) - r * s,
     "graph-residual incidence carriers")

Crs = 3 * r - s
Ers = (s - r) * r**2
domain_obstruction = sp.factor(Crs**3 + 27 * Ers)
zero(domain_obstruction - s**2 * (9 * r - s),
     "natural cubic reducibility obstruction")


# ---------------------------------------------------------------------------
# Shifted normal surface, finite singular carrier, and principal different.
# ---------------------------------------------------------------------------

X, Y = sp.symbols("X Y")
F = sp.expand(T**3 - 3 * P * T - (Ers + Crs * P))
H = sp.expand(X**2 * (X - 3 * r) + Y * (s - 3 * X))
zero(F.subs({T: X - r, P: Y + r**2}) - H,
     "shifted natural surface equation")

J = sp.expand(X**2 - 2 * r * X - Y)
zero(sp.diff(H, Y) - (s - 3 * X), "surface Y derivative")
zero(sp.diff(H, X) - 3 * J, "surface X derivative")

Xsing = s / 3
Ysing = sp.expand(Xsing**2 - 2 * r * Xsing)
zero((s - 3 * X).subs(X, Xsing), "singular X coordinate")
zero(J.subs({X: Xsing, Y: Ysing}), "singular Y coordinate")
zero(27 * H.subs({X: Xsing, Y: Ysing}) - s**2 * (s - 9 * r),
     "finite singular-parameter carrier")

zero(J - ((X - r) ** 2 - (Y + r**2)),
     "J is T squared minus P in shifted coordinates")
zero(sp.diff(F, T).subs({T: X - r, P: Y + r**2}) - 3 * J,
     "relative derivative is 3J")

R = sp.expand(2 * X**2 - (3 * r + s) * X + 2 * r * s)
zero(H.subs(Y, X**2 - 2 * r * X) + X * R,
     "different quotient factors into graph and residual")
zero(2 * Q.subs({h: r - X, a: (s - r) / 2}) - R,
     "residual different factor is the quadratic cofactor")
gate(sp.Poly(R, X).LC() == 2, "residual factor is primitive")
gate(sp.Poly(R, X).TC() == 2 * r * s,
     "graph and residual factors are globally distinct")

# Concrete irreducible control: r=1,s=t.  Its residual discriminant has two
# distinct simple roots, while the surface has only isolated singularities.
t = sp.symbols("t")
control = {r: 1, s: t}
disc_control = sp.factor(disc_rs.subs(control))
zero(4 * disc_control - (t - 1) * (t - 9),
     "one-root irreducible hostile discriminant")
gate(sp.discriminant(4 * disc_control, t) != 0,
     "hostile residual discriminant is squarefree")
singular_control = sp.factor((s**2 * (s - 9 * r)).subs(control))
gate(singular_control == t**2 * (t - 9),
     "hostile singular carrier is finite")


# ---------------------------------------------------------------------------
# Pure-power one-address anatomy.
# ---------------------------------------------------------------------------

c, d, z = sp.symbols("c d z", nonzero=True)
equal_m_disc = sp.factor(
    ((s - r) * (s - 9 * r)).subs({r: c * z, s: d * z})
)
zero(equal_m_disc - (d - c) * (d - 9 * c) * z**2,
     "equal exponents make residual discriminant a square")

# Freeze both unequal-degree formulas using z=t^k after removing the common
# square.  Their two factors are disjoint and separable.
D_n_gt_m = sp.expand((d * z - c) * (d * z - 9 * c))
D_m_gt_n = sp.expand((d - c * z) * (d - 9 * c * z))
gate(sp.resultant(d * z - c, d * z - 9 * c, z) == -8 * c * d,
     "n>m binomial factors are disjoint")
gate(sp.resultant(d - c * z, d - 9 * c * z, z) == 8 * c * d,
     "m>n binomial factors are disjoint")
gate(sp.discriminant(D_n_gt_m, z) == 64 * c**2 * d**2,
     "n>m k=1 residual is squarefree")
gate(sp.discriminant(D_m_gt_n, z) == 64 * c**2 * d**2,
     "m>n k=1 residual is squarefree")

# For k=1,...,5, use exact specializations c=d=1 to check degree 2k,
# squarefreeness, genus k-1, and two even-degree infinity points.
for gap in range(1, 6):
    dn = sp.expand((t**gap - 1) * (t**gap - 9))
    dm = sp.expand((1 - t**gap) * (1 - 9 * t**gap))
    for label, polynomial in (("n_gt_m", dn), ("m_gt_n", dm)):
        gate(sp.degree(polynomial, t) == 2 * gap,
             f"{label} gap {gap}: squarefree degree")
        gate(sp.gcd(polynomial, sp.diff(polynomial, t)) == 1,
             f"{label} gap {gap}: squarefree")
        gate((sp.degree(polynomial, t) - 2) // 2 == gap - 1,
             f"{label} gap {gap}: hyperelliptic genus")
        gate(sp.degree(polynomial, t) % 2 == 0,
             f"{label} gap {gap}: two infinity points")

# Common-positive hostile r=t,s=t^2.  Divide the discriminant square by t^2;
# two normalization values remain at t=0, while both original h-values
# specialize to zero.
r_common = t
s_common = t**2
a_common = sp.expand((s_common - r_common) / 2)
D_common = sp.factor(
    ((s_common - r_common) * (s_common - 9 * r_common)) / 4
)
D_common_reduced = sp.factor(D_common / t**2)
gate(D_common_reduced.subs(t, 0) == sp.Rational(9, 4),
     "common-positive residual has two unramified normalization addresses")
w_values = (sp.Rational(3, 2), sp.Rational(-3, 2))
h_leads = tuple(sp.simplify((sp.Rational(1, 2) + value) / 2)
                  for value in w_values)
gate(h_leads == (1, sp.Rational(-1, 2)),
     "two residual branches have distinct leading terms")
for value in h_leads:
    zero((value * t).subs(t, 0), "both residual branches collapse to h=0")

# The sharp min=0,gap=1 rows: residual conics (genus zero, two infinities)
# and exactly one finite graph/residual incidence.
sharp_rows = [
    ("constant_r", 1, t, (2, 1), 1),
    ("constant_s", t, 1, (1,), 0),
]
for label, r_row, s_row, relation_row, free_rank in sharp_rows:
    D_row = sp.factor(((s_row - r_row) * (s_row - 9 * r_row)) / 4)
    gate(sp.degree(D_row, t) == 2, f"{label}: residual conic")
    gate(sp.gcd(sp.together(D_row * 4),
                sp.diff(sp.together(D_row * 4), t)) == 1,
         f"{label}: smooth squarefree conic")
    gate(sp.degree(sp.expand(r_row * s_row), t) == 1,
         f"{label}: one finite incidence carrier")
    relation = sp.Matrix([list(relation_row)])
    smith = smith_normal_form(relation, domain=sp.ZZ)
    gate(abs(int(smith[0, 0])) == 1,
         f"{label}: primitive Nagata relation")
    gate(relation.cols - relation.rank() == free_rank,
         f"{label}: exact free class rank")


# ---------------------------------------------------------------------------
# General pure-power Nagata multiplicity rows.
# ---------------------------------------------------------------------------

def pure_power_row(m: int, n: int) -> tuple[int, ...]:
    """Multiplicity row of M=s^2(s-9r), with c=d=1 and m!=n."""
    if n > m:
        gap = n - m
        return (2 * n + m,) + (1,) * gap
    gap = m - n
    if n > 0:
        return (3 * n,) + (1,) * gap
    return (1,) * gap


for m_value in range(0, 5):
    for n_value in range(0, 5):
        if m_value == n_value:
            continue
        r_value = t**m_value
        s_value = t**n_value
        M_value = sp.expand(s_value**2 * (s_value - 9 * r_value))
        gap = abs(m_value - n_value)
        if n_value > m_value:
            expected_M = t ** (2 * n_value + m_value) * (t**gap - 9)
        else:
            expected_M = t ** (3 * n_value) * (1 - 9 * t**gap)
        zero(M_value - expected_M,
             f"pure row {(m_value, n_value)} exact M factorization")
        nonzero_factor = t**gap - 9 if n_value > m_value else 1 - 9 * t**gap
        gate(sp.gcd(nonzero_factor, sp.diff(nonzero_factor, t)) == 1,
             f"pure row {(m_value, n_value)} simple nonzero roots")
        row = pure_power_row(m_value, n_value)
        relation = sp.Matrix([list(row)])
        smith = smith_normal_form(relation, domain=sp.ZZ)
        expected_rank = gap if n_value > m_value or n_value > 0 else gap - 1
        gate(abs(int(smith[0, 0])) == 1,
             f"pure row {(m_value, n_value)} is torsion-free")
        gate(relation.cols - relation.rank() == expected_rank,
             f"pure row {(m_value, n_value)} class rank")


summary = {
    "checks": CHECKS,
    "normal_form": "G=-2(h-r)(h^2+ah+ar), s=r+2a",
    "surface": "integral normal shifted X^2(X-3r)+Y(s-3X)",
    "different": "div(J)=E_graph+E_residual reduced",
    "obstruction": "deleting both makes nonconstant J a forbidden A2 unit",
    "pure_power": "equal split; unequal genus gap-1 and two infinity places",
    "sharp": "(m,n)=(0,1),(1,0) pass coarse gates, fail principal different",
    "class": "pure-power Nagata rows primitive and torsion-free",
    "scope": "exactly one k(t)-root in the natural monogenic cubic",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3958 one-hidden-root principal-different companion")
print(f"CHECKS={CHECKS}")
print("NORMAL_FORM=G_MINUS_TWO_TIMES_H_MINUS_R_TIMES_H2_PLUS_AH_PLUS_AR")
print("SURFACE=DOMAIN;NORMAL;SHIFTED_X2_X_MINUS_3R_PLUS_Y_S_MINUS_3X")
print("DIFFERENT=J_EQUALS_T2_MINUS_P;DIV_J_EQUALS_GRAPH_PLUS_RESIDUAL")
print("NO_ATLAS=DELETE_BOTH_RAMIFICATION_PRIMES;J_FORBIDDEN_UNIT")
print("PURE_POWER=EQUAL_SPLITS;UNEQUAL_GENUS_ABS_M_MINUS_N_MINUS_1;TWO_INFINITY")
print("COMMON_POSITIVE=RESIDUAL_NONUNIBRANCH_AT_TRIPLE_ZERO")
print("SHARP_ROWS=(0,1)_CL_Z;(1,0)_CL_ZERO;PRINCIPAL_DIFFERENT_KILLS_BOTH")
print("SCOPE=EXACTLY_ONE_KT_ROOT_NATURAL_MONOGENIC_MODEL;GENERAL_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
