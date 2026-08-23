#!/usr/bin/env python3
"""Canonical independent exact hostile checker for THM-3834.

This checker reconstructs the quotient bracket from the scalar bivector
formula.  It uses denominator-free branch identities and explicit local
unit/repeated-root controls; it does not import or execute the canonical
companion.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(sp.expand(expression)) == 0, message)


def coeff(expression: sp.Expr, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expression), e).coeff_monomial(e**degree)


r, z, e = sp.symbols("r z e")


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    """Scalar-coordinate form, independent of the canonical matrix code."""
    return sp.expand(
        3 * r**2 * (sp.diff(left, r) * sp.diff(right, z)
                    - sp.diff(left, z) * sp.diff(right, r))
        + 9 * z**2 * (sp.diff(left, r) * sp.diff(right, e)
                      - sp.diff(left, e) * sp.diff(right, r))
        + (3 + 6 * r * e) * (sp.diff(left, z) * sp.diff(right, e)
                              - sp.diff(left, e) * sp.diff(right, z))
    )


surface = r**2 * e - z**3 + r
relation = z**3 - r**2 * e - r
zero(bracket(surface, r), "Casimir against r")
zero(bracket(surface, z), "Casimir against z")
zero(bracket(surface, e), "Casimir against e")

f = sp.Function("f")(e)
g = sp.Function("g")(e)
h = sp.Function("h")(e)
kap = sp.Function("kappa")(e)
p = sp.Function("p")(e)
q = sp.Function("q")(e)
S = sp.Function("S")(e)
T = sp.Function("T")(e)
X = sp.Function("X")(e)
Y = sp.Function("Y")(e)

A0 = e**2 - z / 3
C0 = e**3 - e - e * z / 2
A = A0 + r * g + z**2 * f + r * z * p + r * z**2 * S + r**2 * z**2 * X
C = C0 + r * h + z**2 * kap + r * z * q + r * z**2 * T + r**2 * z**2 * Y

zero(bracket(A, C) + bracket(C, A), "ordered bracket is antisymmetric")
zero(bracket(C, -A) - bracket(A, C), "determinant-one target swap preserves bracket value")
gate(sp.degree(C0, e) == 3 and sp.degree(A0, e) == 2,
     "target swap does not preserve the fixed-charge slice")

normal = sp.rem(bracket(A, C) - 1, relation, z)
poly = sp.Poly(sp.expand(normal), r, z)
gate(max(z_degree for _, z_degree in poly.monoms()) <= 2, "monic quotient normal form")
zero(sp.rem(bracket(A, C) - 1 - normal, relation, z), "no quotient remainder loss")


def bucket(r_degree: int, z_degree: int = 0) -> sp.Expr:
    return poly.coeff_monomial(r**r_degree * z**z_degree)


arm_expected = (3 * e**2 - 1) * f - 2 * e * kap + sp.Rational(1, 12)
zero(bucket(0, 1) / 6 - arm_expected, "arm identity")
zero(bucket(7) - 30 * e**2 * (X * sp.diff(Y, e) - Y * sp.diff(X, e)),
     "top Wronskian")

one_sided = {X: 0}
zero(bucket(7).subs(one_sided).doit(), "X=0 top bucket is identically zero")
f0 = sp.Rational(1, 12)
zero(arm_expected.subs({e: 0, f: f0}), "arm origin value")

r6 = sp.expand(bucket(6).subs(one_sided).doit())
r5 = sp.expand(bucket(5).subs(one_sided).doit())
r4 = sp.expand(bucket(4).subs(one_sided).doit())
zero(r6 + 3 * e * (-7 * e * S * sp.diff(Y, e)
                    + 10 * e * Y * sp.diff(S, e) + 2 * S * Y),
     "one-sided r6 identity")

# S != 0: exact 10/7 tower and complete integral, including zero and
# nonzero values of the integration constant.
alpha, beta = sp.symbols("alpha beta", nonzero=True)
gamma = sp.symbols("gamma")
v = sp.Function("v")(e)
Y10 = alpha * e**6 * v**10
S7 = beta * e**4 * v**7
zero(r6.subs({Y: Y10, S: S7}).doit(), "10/7 tower solves r6")

e_residues = [
    (yy, ss) for yy in range(10) for ss in range(7)
    if 7 * yy - 10 * ss == 2
]
other_residues = [
    (yy, ss) for yy in range(10) for ss in range(7)
    if 7 * yy - 10 * ss == 0
]
gate(e_residues == [(6, 4)], "10/7 e-prime fundamental residue class")
gate(other_residues == [(0, 0)], "10/7 other-prime fundamental residue class")

Tint = (
    sp.Rational(10, 7) * alpha / beta * e**2 * v**3 * f
    + sp.Rational(2, 7) * alpha * e**5 * v**10
    + gamma * e**4 * v**7
)
integrated_r5 = sp.expand(r5.subs({Y: Y10, S: S7, T: Tint}).doit())
zero(integrated_r5, "complete r5 integral for arbitrary gamma")
zero(integrated_r5.subs(gamma, 0), "zero integration constant seam")
zero(integrated_r5.subs(gamma, 11), "nonzero integration constant seam")

W = sp.Function("W")(e)
hom_r5 = sp.expand(
    r5.subs({Y: Y10, S: S7, T: Tint + W}).doit()
    - r5.subs({Y: Y10, S: S7, T: Tint}).doit()
)
zero(hom_r5 - 21 * e**2 * (S7 * sp.diff(W, e) - W * sp.diff(S7, e)),
     "r5 homogeneous kernel")
zero(hom_r5.subs(W, 0).doit(), "zero homogeneous correction seam")
zero(hom_r5.subs(W, 13 * S7).doit(), "nonzero scalar homogeneous correction")
zero(S7 * sp.diff(W, e) - W * sp.diff(S7, e)
     - S7**2 * sp.diff(W / S7, e), "kernel is derivative of W/S")

r4_tower = sp.expand(r4.subs({Y: Y10, S: S7, T: Tint}).doit())
P4 = (
    30 * alpha * beta**2 * e**7 * v**14 * sp.diff(v, e)
    + 10 * alpha * beta**2 * e**6 * v**15
    + 20 * alpha * beta * e**4 * f * v**7 * sp.diff(v, e)
    - 20 * alpha * beta * e**4 * v**8 * sp.diff(f, e)
    + 20 * alpha * beta * e**3 * f * v**8
    + 120 * alpha * e * f**2 * sp.diff(v, e)
    - 30 * alpha * e * f * v * sp.diff(f, e)
    + 60 * alpha * f**2 * v
    - 196 * beta**2 * e**3 * kap * v**4 * sp.diff(v, e)
    + 49 * beta**2 * e**3 * v**5 * sp.diff(kap, e)
    - 98 * beta**2 * e**2 * kap * v**5
    + 196 * beta * gamma * e**3 * f * v**4 * sp.diff(v, e)
    - 49 * beta * gamma * e**3 * v**5 * sp.diff(f, e)
    + 98 * beta * gamma * e**2 * f * v**5
)
# Denominator-free form of r4=(3e^3v^2/(7beta))*P4.
zero(7 * beta * r4_tower - 3 * e**3 * v**2 * P4,
     "denominator-free full r4 factor extraction")

# The unit has repeated nonzero roots; d=0 explicitly tests the origin-unit
# boundary.  Gamma remains symbolic, so both its zero and nonzero seams are
# included in every replay.
unit = (1 + e)**3 * (2 + e)**2
unit0 = unit.subs(e, 0)
f_probe = f0 + 3 * e + 5 * e**2
kap_probe = 2 - e + e**3
for d in range(9):
    probe = sp.expand(P4.subs({
        v: e**d * unit,
        f: f_probe,
        kap: kap_probe,
    }).doit())
    expected_lead = 60 * alpha * f0**2 * unit0 * (2 * d + 1)
    zero(coeff(probe, d) - expected_lead, f"origin lead with repeated-root unit, d={d}")
    gate(expected_lead != 0, f"odd origin multiplier nonzero, d={d}")

# S = 0: shifted 5/2 law.  The deliberately unshifted profile is a hostile
# negative control for the fixed-charge +2Yf term.
r5_s0 = sp.expand(r5.subs(S, 0).doit())
r5_s0_expected = -6 * e * (
    5 * e * Y * sp.diff(f, e) - 2 * e * f * sp.diff(Y, e) + 2 * Y * f
)
zero(r5_s0 - r5_s0_expected, "shifted S=0 r5 equation")

w = sp.Function("w")(e)
f52 = beta * w**2
Y52 = alpha * e * w**5
zero(r5_s0_expected.subs({f: f52, Y: Y52}).doit(), "shifted 5/2 tower")
unshifted = sp.expand(r5_s0_expected.subs({f: f52, Y: alpha * w**5}).doit())
zero(unshifted + 12 * alpha * beta * e * w**7,
     "unshifted hostile control exposes fixed-charge term")
gate(unshifted != 0, "unshifted tower is not a solution")

s0_e_residues = [
    (aa, bb) for aa in range(2) for bb in range(5)
    if 5 * aa - 2 * bb == -2
]
s0_other_residues = [
    (aa, bb) for aa in range(2) for bb in range(5)
    if 5 * aa - 2 * bb == 0
]
gate(s0_e_residues == [(0, 1)], "shifted 5/2 e-prime fundamental residue")
gate(s0_other_residues == [(0, 0)], "shifted 5/2 other-prime fundamental residue")

r4_s0_tower = sp.expand(r4.subs({S: 0, f: f52, Y: Y52}).doit())
E4 = (
    7 * e * T * sp.diff(w, e) - 2 * e * w * sp.diff(T, e) + T * w
    + alpha * w**5 * (3 * e * sp.diff(w, e) - w)
)
zero(r4_s0_tower + 6 * beta * e * w * E4,
     "denominator-free S=0 r4 factor extraction")
zero(E4.subs(T, alpha * w**5).doit(), "exact r4 particular solution")

W0 = sp.Function("W0")(e)
Ehom = 7 * e * W0 * sp.diff(w, e) - 2 * e * w * sp.diff(W0, e) + W0 * w
zero(Ehom.subs(W0, 0).doit(), "zero r4 homogeneous correction seam")
w_probe = 1 + e + e**2
for t in range(7):
    W_probe = e**t * (2 + 3 * e)
    local = sp.expand(Ehom.subs({w: w_probe, W0: W_probe}).doit())
    expected = 2 * (1 - 2 * t)
    zero(coeff(local, t) - expected, f"r4 homogeneous origin coefficient, t={t}")
    gate(expected != 0, f"r4 homogeneous half-order forbidden, t={t}")

# Side profiles, recomputed from the quotient normal form.
r4z2 = sp.expand(bucket(4, 2).subs(one_sided).doit())
r4z1 = sp.expand(bucket(4, 1).subs(one_sided).doit())
p_equation = -2 * Y * sp.diff(p, e) + p * sp.diff(Y, e)
g_equation = (-10 * e * Y * sp.diff(g, e)
              + 3 * e * g * sp.diff(Y, e) + 2 * Y * g)
zero(r4z2 - 15 * e * p_equation, "p Wronskian bucket")
zero(r4z1 - 3 * g_equation, "g shifted Wronskian bucket")

Y_probe = e * w_probe**5
for t in range(7):
    profile = e**t * (2 + 3 * e)
    p_local = sp.expand(p_equation.subs({Y: Y_probe, p: profile}).doit())
    g_local = sp.expand(g_equation.subs({Y: Y_probe, g: profile}).doit())
    p_expected = 2 * (1 - 2 * t)
    g_expected = 10 * (1 - 2 * t)
    zero(coeff(p_local, t) - p_expected, f"p origin coefficient, t={t}")
    zero(coeff(g_local, t + 1) - g_expected, f"g origin coefficient, t={t}")
    gate(p_expected != 0 and g_expected != 0,
         f"p/g half-orders forbidden, t={t}")
zero(p_equation.subs(p, 0).doit(), "p=0 seam")
zero(g_equation.subs(g, 0).doit(), "g=0 seam")

r3z2 = sp.expand(bucket(3, 2).subs({
    X: 0, S: 0, Y: Y52, f: f52, T: alpha * w**5, p: 0, g: 0,
}).doit())
zero(r3z2 + 10 * alpha * e**2 * w**4 * sp.diff(w, e),
     "terminal r3z2 identity")
zero(r3z2.subs(w, 2).doit(), "constant-w positive control")
gate(sp.expand(r3z2.subs(w, 1 + e).doit()) != 0,
     "nonconstant-w hostile control")

zero(arm_expected.subs(f, f0) - e * (e / 4 - 2 * kap),
     "arm uses cancellation-safe factor e")
zero(arm_expected.subs({f: f0, kap: e / 8}), "constant arm solution")

a = sp.symbols("a", nonzero=True)
terminal = sp.expand(bucket(3, 1).subs({
    X: 0, S: 0, Y: a * e, f: f0, T: a, kap: e / 8, p: 0, g: 0,
}).doit())
zero(terminal + 60 * a * e**3, "final r3z monomial")
gate(not terminal.has(h, q), "unused h and q cannot cancel the terminal")
gate(terminal != 0, "Y-nonzero terminal is nonzero")
zero(terminal.subs(a, 0), "forbidden Y=0 seam would remove terminal")

# Ordered support partition: the X!=0 row splits by L=0/L!=0, while X=0
# splits into Y=0/Y!=0.  This records scope, not a symmetry reduction.
partition = {
    "X=Y=0": "THM-3821",
    "X!=0,L=0": "THM-3828",
    "X!=0,L!=0": "THM-3829",
    "X=0,Y!=0": "THM-3834",
}
gate(len(partition) == 4 and set(partition.values()) == {
    "THM-3821", "THM-3828", "THM-3829", "THM-3834"
}, "fixed ordered ansatz partition")

semantic = {
    "proof_commit": "d8e3633b8561a31085163400981208bdac2a38ba",
    "promotion_commit": "f6cd25000d84ad7a0379cd00f52ff47425678a36",
    "theorem_blob": "048bc8341cce4520c3ad2e2c8c30dcf52a577da6",
    "orientation": "ordered X=0,Y!=0 treated directly; target swap exits fixed-charge slice",
    "S_nonzero": "10/7 UFD tower; complete gamma integral; denominator-free r4; odd origin lead",
    "S_zero": "shifted 5/2 tower; unique T; p=g=0; w constant; -60ae3",
    "zero_seams": "gamma=0;W=0;p=0;g=0;d=0 checked;Y=0 explicitly forbidden",
    "scope": "only displayed fixed r2z2 ansatz; higher slots and general planar JC open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "independent checker is assertion-free")

print("audit=THM-3834-independent-hostile-checker")
print("orientation=ordered_direct;swap_not_in_fixed_charge_slice")
print("S_nonzero=UFD_10_7;gamma_zero_nonzero;denominator_free_r4;d0_and_repeated_root_unit")
print("S_zero=shifted_5_2;unshifted_negative_control;W_p_g_zero_nonzero;terminal")
print("partition=THM3821+THM3828+THM3829+THM3834;fixed_ansatz_only")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
