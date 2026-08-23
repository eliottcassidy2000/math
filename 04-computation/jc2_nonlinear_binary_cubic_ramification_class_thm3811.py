#!/usr/bin/env python3
"""Exact companion for THM-3811's nonlinear cubic ramification packet."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: object, rhs: object, label: str) -> None:
    difference = lhs - rhs  # type: ignore[operator]
    if isinstance(difference, sp.MatrixBase):
        reduced = difference.applyfunc(lambda entry: sp.cancel(sp.expand(entry)))
        check(reduced == sp.zeros(*difference.shape), label)
        return
    check(sp.cancel(sp.expand(difference)) == 0, label)


# Universal Delone--Faddeev matrices, index, discriminant, and norm.
a, b, c, d = sp.symbols("a b c d")
x, y, p, q1, q2, T = sp.symbols("x y p q1 q2 T")
I3 = sp.eye(3)
M_omega = sp.Matrix(
    [
        [0, -a * c, -a * d],
        [1, b, 0],
        [0, -a, 0],
    ]
)
M_theta = sp.Matrix(
    [
        [0, -a * d, -b * d],
        [0, 0, d],
        [1, 0, -c],
    ]
)
check(M_omega * M_theta == M_theta * M_omega,
      "universal multiplication commutes")
index_matrix = sp.Matrix.hstack(
    sp.Matrix([1, 0, 0]),
    sp.Matrix([0, x, y]),
    (x * M_omega + y * M_theta) ** 2 * sp.Matrix([1, 0, 0]),
)
binary = a * x**3 + b * x**2 * y + c * x * y**2 + d * y**3
same(index_matrix.det(), -binary, "universal index sign")

Delta = (
    b**2 * c**2 - 4 * a * c**3 - 4 * b**3 * d
    - 27 * a**2 * d**2 + 18 * a * b * c * d
)
basis = [I3, M_omega, M_theta]
trace_pairing = sp.Matrix(
    3, 3, lambda i, j: sp.trace(basis[i] * basis[j])
)
same(trace_pairing.det(), Delta, "universal trace discriminant")

norm_universal = sp.factor((p * I3 + q1 * M_omega + q2 * M_theta).det())
norm_formula = (
    a**2 * d * q1**3 + 2 * a * b * d * q1**2 * q2
    - a * c**2 * q1**2 * q2 - 2 * a * c * d * q1 * q2**2
    + a * c * p * q1**2 - a * d**2 * q2**3
    + 3 * a * d * p * q1 * q2 + b**2 * d * q1 * q2**2
    - b * c * p * q1 * q2 + b * d * p * q2**2
    + b * p**2 * q1 - c * p**2 * q2 + p**3
)
same(norm_universal, norm_formula, "universal norm formula")


# Nonlinear packet.
A, C = sp.symbols("A C")
packet = {a: A, b: C, c: 7 * A, d: C**2 - 3 * A}
Mw = sp.simplify(M_omega.subs(packet))
Mt = sp.simplify(M_theta.subs(packet))
same(Mw**2, -7 * A**2 * I3 + C * Mw - A * Mt,
     "nonlinear omega-square")
same(Mw * Mt, (3 * A**2 - A * C**2) * I3,
     "nonlinear mixed law")
same(Mt**2,
     (3 * A * C - C**3) * I3 + (C**2 - 3 * A) * Mw - 7 * A * Mt,
     "nonlinear theta-square")

index_packet = sp.factor(index_matrix.det().subs(packet))
index_expected = -(A * x**3 + C * x**2 * y + 7 * A * x * y**2
                   + (C**2 - 3 * A) * y**3)
same(index_packet, index_expected, "nonlinear index form")
check(all(coefficient.subs({A: 0, C: 0}) == 0 for coefficient in
          (A, C, 7 * A, C**2 - 3 * A)),
      "square-zero coefficient fibre")

Delta_packet = sp.factor(Delta.subs(packet))
Delta_linear = A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
Delta_expected = (
    Delta_linear + C**2 * (162 * A**3 + 126 * A**2 * C - 4 * C**3)
    - 27 * A**2 * C**4
)
same(Delta_packet, Delta_expected, "nonlinear discriminant packet")
squarefree_gcd = sp.gcd(
    sp.gcd(sp.Poly(Delta_packet, A, C),
           sp.Poly(sp.diff(Delta_packet, A), A, C)),
    sp.Poly(sp.diff(Delta_packet, C), A, C),
)
check(squarefree_gcd == 1, "nonlinear discriminant geometrically reduced")
same(Delta_packet.subs(A, 0), -4 * C**5,
     "no discriminant component lies in A=0")

# The degree-four initial form is the rational four-line Veronese packet.
initial_terms = sum(
    coefficient * A**monomial[0] * C**monomial[1]
    for monomial, coefficient in sp.Poly(Delta_packet, A, C).terms()
    if sum(monomial) == 4
)
same(initial_terms, Delta_linear, "four-line tangent cone")


# Generic cubic field and the specialization rational-root gate.
char_omega = sp.factor(Mw.charpoly(T).as_expr())
f = T**3 - C * T**2 + 7 * A**2 * T + 3 * A**3 - A**2 * C**2
same(char_omega, f, "nonlinear omega characteristic polynomial")
check(sp.Poly(f, T, domain=sp.QQ.frac_field(A, C)).is_irreducible,
      "generic characteristic cubic irreducible over Q(A,C)")
beta = sp.symbols("beta")
special_linear_root = sp.Poly(
    sp.expand(f.subs({A: 1, T: C + beta})), C
)
same(special_linear_root.nth(2), beta - 1,
     "specialized rational-root quadratic coefficient")
same(special_linear_root.nth(1), 2 * beta**2 + 7,
     "specialized rational-root linear coefficient")
check((2 * beta**2 + 7).subs(beta, 1) == 9,
      "only leading root beta=1 fails next coefficient")
same(sp.discriminant(f, T), A**2 * Delta_packet,
     "power-basis index-square discriminant")


# Rational normalization of the irreducible ramification curve.
q = sp.symbols("q")
Rq = (q - 3) * (q + 1) * (q + 2)
Dq = 3 * q**2 + 7
Aq = -2 * q**2 * Rq / Dq**2
Cq = -q * Rq / Dq
wq = Aq * q
thetaq = Aq * (q**2 - 7) / 2
same(Aq * (q**3 + 7 * q + 3), Cq * (Cq + q**2),
     "rational A-nonzero chart")
same(f.subs({A: Aq, C: Cq, T: wq}), 0,
     "ramification root lies on characteristic cubic")
same(sp.diff(f, T).subs({A: Aq, C: Cq, T: wq}), 0,
     "ramification root is double")
same(Delta_packet.subs({A: Aq, C: Cq}), 0,
     "ramification parametrization lies on discriminant")
same(wq**2, -7 * Aq**2 + Cq * wq - Aq * thetaq,
     "ramification omega-square relation")
same(wq * thetaq, 3 * Aq**2 - Aq * Cq**2,
     "ramification mixed relation")
same(thetaq**2,
     3 * Aq * Cq - Cq**3 + (Cq**2 - 3 * Aq) * wq - 7 * Aq * thetaq,
     "ramification theta-square relation")

for root in (-2, -1, 0, 3):
    check(Dq.subs(q, root) != 0, f"normalization denominator at q={root}")
    same(Aq.subs(q, root), 0, f"vertex A at q={root}")
    same(Cq.subs(q, root), 0, f"vertex C at q={root}")
    same(wq.subs(q, root), 0, f"vertex omega at q={root}")
    same(thetaq.subs(q, root), 0, f"vertex theta at q={root}")

companion_root = Cq - 2 * wq
same(
    f.subs({A: Aq, C: Cq}),
    (T - wq) ** 2 * (T - companion_root),
    "double ramification plus simple companion factorization",
)
companion_derivative = sp.factor(
    sp.diff(f, T).subs({A: Aq, C: Cq, T: companion_root})
)
check(companion_derivative != 0,
      "companion derivative is a nonzero rational function")

# The power-basis different has the expected polluted norm A^2*Delta;
# the index square is precisely why it does not principalize ramification.
different_omega = 3 * Mw**2 - 2 * C * Mw + 7 * A**2 * I3
same(different_omega.det(), -A**2 * Delta_packet,
     "omega different norm has index-square pollution")


# Exact affine-profile exponent-one norm-Pell gate.
P0, PA, PC = sp.symbols("P0 PA PC")
Q0, QA, QC = sp.symbols("Q0 QA QC")
R0, RA, RC = sp.symbols("R0 RA RC")
kappa, kappa_inverse = sp.symbols("kappa kappa_inverse")
p_aff = P0 + PA * A + PC * C
q1_aff = Q0 + QA * A + QC * C
q2_aff = R0 + RA * A + RC * C
norm_aff = sp.expand(norm_formula.subs(packet).subs({
    p: p_aff,
    q1: q1_aff,
    q2: q2_aff,
}))
coefficient_equations = sp.Poly(
    norm_aff - kappa * Delta_packet, A, C
).coeffs()
check(len(coefficient_equations) == 38,
      "affine norm-Pell coefficient equation count")
groebner_variables = (
    P0, PA, PC, Q0, QA, QC, R0, RA, RC, kappa, kappa_inverse
)
groebner_gate = sp.groebner(
    coefficient_equations + [kappa * kappa_inverse - 1],
    *groebner_variables,
    order="grevlex",
)
check(len(groebner_gate.polys) == 1
      and groebner_gate.polys[0].as_expr() == 1,
      "affine exponent-one norm-Pell ideal is the unit ideal")

# Total-degree floor: the degree-2 coefficient d makes 3D+5 the largest
# possible norm degree, while Delta has genuine top degree six.
same(sp.Poly(Delta_packet, A, C).total_degree(), 6,
     "nonlinear discriminant total degree six")
for n in range(2, 21):
    lower = 2 * n - 1
    check(3 * (lower - 1) + 5 < 6 * n <= 3 * lower + 5,
          f"torsion norm degree floor n={n}")


source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "criterion": "U*/S*=kernel of the deleted-prime lattice map to Cl(S); constant units iff ramification classes are Z-independent",
    "norm": "for one branch, torsion iff div(u)=nE, equivalently Norm(u)=kappa*Delta^n plus zero companion valuation",
    "packet": "(a,b,c,d)=(A,C,7A,C^2-3A), normal nonmonogenic squarefree S3 cubic",
    "ramification": "one irreducible rational E; q=-2,-1,0,3 are four normalization points over the square-zero vertex; one simple companion",
    "finite_gate": "affine coefficient profiles cannot solve the exponent-one norm equation; Groebner basis [1]",
    "open": "infinite order of E, S*=k*, and affineness of the etale complement remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet")
print("units=U*/S*=kernel(Z^ramification->Cl(S))")
print("packet=(A,C,7A,C^2-3A);normal_nonmonogenic_squarefree_S3")
print("branch=one_irreducible_rational_E;four_vertex_preimages=-2,-1,0,3")
print("companion=one_explicit_simple_root_C-2Aq")
print("norm_gate=affine_exponent1_Groebner_[1];higher_order_degree_floor=2n-1")
print("open=class_E_infinite_order+S_units_constant+etale_complement_affine")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
