#!/usr/bin/env python3
"""Exact companion for THM-3813's uniform quartic r-repair obstruction."""

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
    gate(sp.cancel(expression) == 0, message)


r, z, e, u = sp.symbols("r z e u")
b0, b1, b2, b3, b4 = sp.symbols("b0 b1 b2 b3 b4")
variables = (r, z, e)
surface = r**2 * e - z**3 + r
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 + 6 * r * e],
        [-9 * z**2, -3 - 6 * r * e, 0],
    ]
)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


g = b0 + b1 * e + b2 * e**2 + b3 * e**3 + b4 * e**4
gp = sp.diff(g, e)
A = e**2 - z / 3 + r * g
K_source = 1 + 2 * r * e
critical = [sp.factor(bracket(A, q)) for q in variables]
zero(
    critical[0] - (r**2 - 9 * z**2 * (2 * e + r * gp)),
    "quartic r Hamiltonian",
)
zero(
    critical[1] - (3 * g * r**2 - 3 * K_source * (2 * e + r * gp)),
    "quartic z Hamiltonian",
)
zero(critical[2] - (9 * g * z**2 - K_source), "quartic e Hamiltonian")
zero(bracket(A, surface), "surface Casimir")


# Universal two-equation compression and residual resultant.
G, D = sp.symbols("G D")
K = 1 + 2 * u
P = sp.expand(G * u**2 - K * (2 * e**3 + u * e * D))
Q = sp.expand(e**2 * K**3 - 729 * G**3 * u**2 * (1 + u) ** 2)
resultant_generic = sp.factor(sp.resultant(P, Q, u))
H_generic = sp.factor(resultant_generic / (G**3 * e**4))
zero(resultant_generic - G**3 * e**4 * H_generic, "universal resultant factor")
zero(
    H_generic.subs(G, 0)
    - D * e * (-729 * e * (D - 4 * e**2) ** 3 - 2),
    "universal nonzero-arm law",
)

H = sp.Poly(sp.expand(H_generic.subs({G: g, D: gp})), e)
gate(H.degree() == 22, "quartic residual degree twenty-two")
gate(H.LC() == 19131876 * b4**5, "quartic residual leading coefficient")
gate(
    H.as_expr().subs({b0: 2, b1: -3, b2: 5, b3: -7, b4: 11}) != 0,
    "positive residual control",
)


# Boundary-only support (with arbitrary repeated roots) forces H|e*g*H'.
# Compute the Euclidean quotient over Q(b0,...,b4), then clear exactly its
# only parameter denominator, 81*b4^3.
quotient, remainder = sp.div(
    sp.Poly(sp.expand(e * g * sp.diff(H.as_expr(), e)), e, domain="EX"),
    sp.Poly(H.as_expr(), e, domain="EX"),
)
expected_quotient = (
    22 * b4 * e**4
    + (53 * b3 + 4) * e**3 / 3
    + (44 * b2 / 3 - 4 * b3**2 / (9 * b4)
       - 4 * b3 / (9 * b4) + 8 / (9 * b4)) * e**2
    + (13 * b1 - 4 * b2 * b3 / (3 * b4)
       + 8 * b3**3 / (27 * b4**2) - 8 * b3 / (9 * b4**2)
       + 16 / (27 * b4**2)) * e
    + 38 * b0 / 3 - 4 * b1 * b3 / (3 * b4)
    + 4 * b1 / (3 * b4) - 8 * b2**2 / (9 * b4)
    + 28 * b2 * b3**2 / (27 * b4**2)
    - 20 * b2 * b3 / (27 * b4**2) - 8 * b2 / (27 * b4**2)
    - 16 * b3**4 / (81 * b4**3) + 16 * b3**3 / (81 * b4**3)
    + 16 * b3**2 / (27 * b4**3) - 80 * b3 / (81 * b4**3)
    + 32 / (81 * b4**3)
)
zero(quotient.as_expr() - expected_quotient, "quartic logarithmic quotient")
gate(remainder.degree() == 21, "quartic logarithmic remainder degree")
q_clear = sp.factor(81 * b4**3 * quotient.as_expr())
gate(sp.Poly(q_clear, e, b0, b1, b2, b3, b4).is_zero is False,
     "cleared quotient is polynomial")
R = sp.Poly(
    sp.expand(81 * b4**3 * e * g * sp.diff(H.as_expr(), e) - q_clear * H.as_expr()),
    e,
)
zero(R.as_expr() - 81 * b4**3 * remainder.as_expr(), "cleared remainder")
gate(R.degree() == 21, "cleared remainder degree")


# The residual system has a hidden Z/7 grading.  Every b4!=0 point maps to
# the weight-zero coordinates below, with T=b4^7!=0.  Pull each top remainder
# coefficient into these coordinates and remove its common Laurent b4-power.
C0, C1, C2, C3, T, s = sp.symbols("C0 C1 C2 C3 T s")
parameter_pullback = {
    b0: C0 / s**3,
    b1: C1 / s**2,
    b2: C2 / s,
    b3: C3,
    b4: s,
}


def invariant_primitive(expression: sp.Expr) -> tuple[int, sp.Expr]:
    pulled = sp.expand(sp.cancel(expression.subs(parameter_pullback)))
    terms = sp.Add.make_args(pulled)
    powers = [int(term.as_powers_dict().get(s, 0)) for term in terms]
    minimum = min(powers)
    normalized = sp.expand(pulled / s**minimum)
    invariant = sp.Integer(0)
    for term in sp.Add.make_args(normalized):
        power = int(term.as_powers_dict().get(s, 0))
        gate(power % 7 == 0, "invariant pullback has one residue class")
        invariant += term / s**power * T ** (power // 7)
    polynomial = sp.Poly(sp.expand(invariant), C0, C1, C2, C3, T)
    content, primitive = polynomial.primitive()
    gate(content != 0, "nonzero invariant primitive content")
    return minimum, primitive.as_expr()


top_degrees = tuple(range(21, 14, -1))
top_packet: list[sp.Expr] = []
laurent_powers: list[int] = []
for degree in top_degrees:
    power, primitive = invariant_primitive(R.coeff_monomial(e**degree))
    laurent_powers.append(power)
    top_packet.append(primitive)
gate(laurent_powers == [4, 3, 2, 1, 0, -1, -2], "top Laurent powers")

# Direct back-substitution verifies that no coefficient or boundary stratum
# was lost in forming each primitive invariant equation.
for degree, power, primitive in zip(top_degrees, laurent_powers, top_packet):
    raw_pullback = sp.expand(
        sp.cancel(R.coeff_monomial(e**degree).subs(parameter_pullback) / s**power)
    )
    primitive_lift = sp.expand(primitive.subs(T, s**7))
    ratio = sp.cancel(raw_pullback / primitive_lift)
    gate(not ratio.has(C0, C1, C2, C3, s), "primitive differs by scalar only")
    gate(ratio != 0, "primitive scalar is nonzero")


# Five top equations define I.  Exact Groebner division of the next two
# equations leaves incompatible affine laws for C2 after localizing at T.
I_basis = sp.groebner(
    top_packet[:5], C0, C1, C2, C3, T, order="grevlex"
)
gate(len(I_basis.polys) == 28, "five-equation Groebner basis size")
normal_16 = sp.factor(I_basis.reduce(top_packet[5])[1])
normal_15 = sp.factor(I_basis.reduce(top_packet[6])[1])
expected_16 = -T * (2688 * C2 - 15308821) / 4
expected_15 = 27 * T * (378784 * C2 - 2308816611) / 8
zero(normal_16 - expected_16, "e16 exact normal form")
zero(normal_15 - expected_15, "e15 exact normal form")
cross_difference = 15308821 * 378784 - 2308816611 * 2688
gate(cross_difference == -407362596704, "incompatible C2 cross-difference")
gate(cross_difference != 0, "two localized affine laws are incompatible")


# Conceptual hostile seam: the only factor that appears when one maximizes
# contact at a nonzero simple root is not a simple-root resonance at all.
# The arm law alpha=-2/(729*w^3) and 531441*w^7+16=0 force g'(alpha)=0.
alpha, w = sp.symbols("alpha w", nonzero=True)
alpha_arm = -sp.Rational(2, 729) / w**3
simple_derivative = 4 * alpha_arm**2 + w
zero(
    sp.expand(simple_derivative * w**6 - (w**7 + sp.Rational(16, 531441))),
    "apparent simple-arm resonance is a repeated root",
)


# Exact reconstruction from any residual root eta with eta*g(eta)!=0.
r_rec = u / e
z_rec = 9 * G * u * (1 + u) / (e * K)
zero(z_rec**2 - K / (9 * G) + Q / (9 * G * e**2 * K**2),
     "Q reconstructs z square")
zero(
    surface.subs({r: r_rec, z: z_rec})
    - u * (1 + u) * Q / (e**3 * K**3),
    "Q reconstructs surface",
)
A_e_symbol = 2 * e + r_rec * D
C_e_rec = sp.factor(9 * G * z_rec**2 - K)
C_z_rec = sp.factor(3 * G * r_rec**2 - 3 * K * A_e_symbol)
zero(C_e_rec + Q / (e**2 * K**2), "Q kills e Hamiltonian")
zero(C_z_rec - 3 * P / e**2, "P kills z Hamiltonian")
C_r_rec = sp.factor(r_rec**2 - 9 * z_rec**2 * A_e_symbol)
zero(
    C_r_rec - P / (G * e**2)
    - Q * (u**2 - P / G) / (e**4 * K**3),
    "P and Q kill r Hamiltonian",
)
zero(Q.subs(u, 0) - e**2, "exclude u=0")
zero(Q.subs(u, -1) + e**2, "exclude u=-1")
zero(Q.subs(u, -sp.Rational(1, 2)) + 729 * G**3 / 16, "exclude K=0")
zero(sp.LC(sp.Poly(Q, u)) + 729 * G**3, "Q excludes infinity")


# Frozen exact packet and hostile controls.  Repeated roots and sparse quartics
# must be rejected by the same seven-equation mechanism, not by squarefreeness.
for profile in (
    (1, -4, 6, -4, 1),       # (e-1)^4
    (0, 0, 1, -2, 1),        # e^2(e-1)^2
    (0, 0, 0, 0, 1),         # e^4
    (2, -3, 5, -7, 11),      # dense hostile control
):
    values = dict(zip((b0, b1, b2, b3, b4), profile))
    evaluated = [sp.expand(R.coeff_monomial(e**degree).subs(values))
                 for degree in top_degrees]
    gate(any(value != 0 for value in evaluated), f"hostile quartic {profile}")

packet_blob = "\n".join(sp.srepr(sp.expand(item)) for item in top_packet).encode()
basis_blob = "\n".join(
    sp.srepr(sp.expand(poly.as_expr())) for poly in I_basis.polys
).encode()
semantic = {
    "carrier": "A=e^2-z/3+r*sum(b_i e^i,0<=i<=4); c=1",
    "new_case": "b4!=0; degree<=3 inherited through THM-3807",
    "resultant": "Res_u(P,Q)=e^4*g^3*H; deg(H)=22; LC=19131876*b4^5",
    "criterion": "boundary-only roots, with repetitions, imply H divides e*g*Hprime",
    "invariants": "C0=b0*b4^3,C1=b1*b4^2,C2=b2*b4,C3=b3,T=b4^7",
    "obstruction": "top seven equations; NF16 and NF15 force incompatible C2 values",
    "cross_difference": cross_difference,
    "recovery": "r=u/e; z=9g*u*(u+1)/(e*(1+2u))",
    "open": "degree>=5 g; mixed z2/r corrections; different arm profiles",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3813-quartic-r-repairs-of-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+r*(b0+b1e+b2e^2+b3e^3+b4e^4)")
print("new_case=b4!=0;degree_le_3_inherited=THM-3807")
print("resultant=Res_u(P,Q)=e^4*g^3*H;deg_H=22;LC_H=19131876*b4^5")
print("boundary_criterion=all_H_roots_in_V(e*g)_implies_H_divides_e*g*Hprime")
print("invariants=C0:b0*b4^3,C1:b1*b4^2,C2:b2*b4,C3:b3,T:b4^7")
print("top_equations=degrees_21_to_15;first_five_form_Groebner_ideal")
print("normal_e16=-T*(2688*C2-15308821)/4")
print("normal_e15=27*T*(378784*C2-2308816611)/8")
print(f"cross_difference={cross_difference}")
print("recovery=r=u/e;z=9g*u*(u+1)/(e*(1+2u))")
print("open=degree_ge_5_g;mixed_z2_r;different_arm_profile;Darboux_pair")
print(f"invariant_packet_sha256={hashlib.sha256(packet_blob).hexdigest()}")
print(f"groebner_basis_sha256={hashlib.sha256(basis_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
