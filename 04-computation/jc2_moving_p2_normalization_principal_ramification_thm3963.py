#!/usr/bin/env python3
"""Exact companion for THM-3963's moving P^2 normalization closure.

Reproduction:
  python3 04-computation/jc2_moving_p2_normalization_principal_ramification_thm3963.py
  python3 -O 04-computation/jc2_moving_p2_normalization_principal_ramification_thm3963.py
"""

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


v, w, P, T, c, t = sp.symbols("v w P T c t")


# ---------------------------------------------------------------------------
# Generic normalization parameter and its integral replacement at c=0.
# ---------------------------------------------------------------------------

P_v = v**2 * (3 + c*v)
T_v = v * (3 + c*v)
F = T**3 - 3*P*T - c*P**2

zero(F.subs({P: P_v, T: T_v}), "generic v parametrization lies on F=0")
zero(P_v-v*T_v, "v=P/T in the cubic function field")
gate(sp.degree(P_v, v) == 3 and sp.LC(sp.Poly(P_v, v)) == c,
     "nonzero c gives a separable degree-three P-map")

w_v = c*v
relations = {
    "wT=cP": w*T-c*P,
    "w(w+3)=cT": w*(w+3)-c*T,
    "P(w+3)=T2": P*(w+3)-T**2,
    "c2P=w2(w+3)": c**2*P-w**2*(w+3),
    "w monic": w**3+3*w**2-c**2*P,
}
for label, relation in relations.items():
    zero(relation.subs({P: P_v, T: T_v, w: w_v}),
         f"normalization relation {label}")
zero(sp.Poly(relations["w(w+3)=cT"], w).monic().as_expr()
     - relations["w(w+3)=cT"],
     "w satisfies a monic quadratic over the order")
gate(sp.Poly(relations["w monic"], w).monic().as_expr() ==
     relations["w monic"], "the target-monic cubic relation is compatible")
zero((relations["w(w+3)=cT"]*T - (w+3)*relations["wT=cP"])
     + c*(T**2-P*(w+3)),
     "first chart derives P=T^2/(w+3) by domain cancellation")
zero(w*relations["wT=cP"] + c*relations["P(w+3)=T2"]
     - T*relations["w(w+3)=cT"] + 3*relations["wT=cP"],
     "normalization syzygy consistency")
gate(sp.gcd(w, w+3) == 1, "D(w) and D(w+3) cover the normalization")

# The normalization is R-free on (1,T,w).  The multiplication matrices also
# independently recover its discriminant and the old-order index.
M_T = sp.Matrix([[0, 3*P, c*P], [1, 0, 0], [0, P, 0]])
M_w = sp.Matrix([[0, c*P, 0], [0, 0, c], [1, 0, -3]])
M_basis = [sp.eye(3), M_T, M_w]
trace_pairing = sp.Matrix([
    [sp.trace(left*right) for right in M_basis] for left in M_basis
])
disc_B = sp.factor(trace_pairing.det())
disc_A = sp.factor(sp.discriminant(F, T))
gate(disc_B == -27*P*(c**2*P-4), "normalization discriminant")
gate(disc_A == -27*P**3*(c**2*P-4), "monogenic order discriminant")
transition_A_to_B = sp.Matrix([[1, 0, 3*P], [0, 1, 0], [0, 0, P]])
gate(transition_A_to_B.det() == P and disc_A == P**2*disc_B,
     "old-to-normal basis index is exactly P")
zero((w**3+3*w**2-4) - (w+2)**2*(w-1),
     "branch polynomial has ramified E2 and unramified companion")

# Nagata's vertical relation lattice has two primes over each distinct root
# of c and one primitive sum relation per root.
for root_count in range(1, 6):
    relation = sp.zeros(root_count, 2*root_count)
    for index in range(root_count):
        relation[index, 2*index] = 1
        relation[index, 2*index+1] = 1
    smith = smith_normal_form(relation, domain=sp.ZZ)
    gate(relation.rank() == root_count and
         all(abs(int(smith[i, i])) == 1 for i in range(root_count)),
         f"{root_count}-root vertical class ledger is torsion-free rank {root_count}")


# ---------------------------------------------------------------------------
# The two exact affine charts and their smoothness.
# ---------------------------------------------------------------------------

# On D(w+3): P=T^2/(w+3), with one equation cT=w(w+3).
f0 = c*T-w*(w+3)
P_chart0 = T**2/(w+3)
zero(sp.cancel(relations["wT=cP"].subs(P, P_chart0)*(w+3) + T*f0),
     "D(w+3) presentation recovers wT=cP")
zero(sp.cancel(F.subs(P, P_chart0)*(w+3)**2 + T**3*f0),
     "D(w+3) presentation recovers the cubic equation")

# f0=f0_T=f0_w=0 is impossible already without the t derivative.
f0_T = sp.diff(f0, T)
f0_w = sp.diff(f0, w)
gate(f0_T == c and f0_w == -2*w-3,
     "first-chart decisive Jacobian rows")
G0 = sp.groebner([f0, f0_T, f0_w], T, w, c, order="grevlex")
gate(len(G0.polys) == 1 and G0.polys[0].as_expr() == 1,
     "D(w+3) chart is universally smooth, including zeros of c")

# On D(w): T=cP/w, with one monic equation w^3+3w^2-c^2P=0.
f_inf = c**2*P-w**2*(w+3)
T_chart_inf = c*P/w
zero(sp.cancel(relations["w(w+3)=cT"].subs(T, T_chart_inf)*w
     + f_inf), "D(w) presentation recovers w(w+3)=cT")
zero(sp.cancel(F.subs(T, T_chart_inf)*w**3
     - c*P**2*f_inf), "D(w) presentation recovers the cubic equation")
f_inf_P = sp.diff(f_inf, P)
f_inf_w = sp.diff(f_inf, w)
gate(f_inf_P == c**2, "second-chart P derivative")
zero(f_inf_w + 3*w*(w+2), "second-chart w derivative")
G_inf = sp.groebner([f_inf, f_inf_P, f_inf_w, sp.Symbol("z")*w-1],
                    P, w, c, sp.Symbol("z"), order="grevlex")
gate(len(G_inf.polys) == 1 and G_inf.polys[0].as_expr() == 1,
     "D(w) chart is universally smooth, including zeros of c")


# ---------------------------------------------------------------------------
# Relative ramification and the principal forbidden-unit prime.
# ---------------------------------------------------------------------------

# On the first localization, the relation matrix over k[P,t] comes from
# cT-w(w+3)=0 and (w+3)P-T^2=0.
jacobian_matrix_0 = sp.Matrix([
    [c, -(2*w+3)],
    [-2*T, P],
])
det0 = sp.expand(jacobian_matrix_0.det())
zero(det0 + 3*T*(w+2) + relations["wT=cP"],
     "first-chart ramification determinant identity avoids division by c")

# On D(w), relative different is the derivative of the monic w equation.
different_inf = sp.diff(w**3+3*w**2-c**2*P, w)
zero(different_inf-3*w*(w+2), "second-chart relative different")

E2_relation = sp.factor(f_inf.subs(w, -2))
zero(E2_relation-(c**2*P-4), "principal E2 quotient relation")
gate(sp.Poly(c**2*P-4, P).degree() == 1,
     "E2 quotient is the domain k[t,1/c]")
zero(c*(c*P/4)-1 + (4-c**2*P)/4,
     "c is automatically a unit on E2")
gate(different_inf.subs(w, -2) == 0,
     "the whole prime E2=(w+2) is ramified")
gate(sp.factor(sp.diff(different_inf, w).subs(w, -2)) == -6,
     "E2 is reduced in the ramification divisor")

# The conductor-origin packet is visible but is not needed for the unit kill.
zero(relations["w monic"].subs(P, 0) - w**2*(w+3),
     "the target zero section has normalization addresses w=0,-3")
gate(different_inf.subs(w, -3) == 9,
     "the w=-3 companion address is etale")
gate(sp.factor(det0.subs({w: 0, P: 0})) == -6*T,
     "the w=0 section ramifies only at T=0 in the moving vertical fibre")


# ---------------------------------------------------------------------------
# Hostile c=t fibre: the escaping v=-3/t branch is retained by w=-3.
# ---------------------------------------------------------------------------

F_t = sp.expand(F.subs(c, t))
zero(F_t.subs(t, 0) - T*(T**2-3*P),
     "c=t has a reducible special fibre despite an integral total surface")
gate(f0.subs({c: t, t: 0, w: 0}) == 0,
     "w=0 chart retains the first vertical special-fibre component")
gate(f_inf.subs({c: t, t: 0, w: -3}) == 0,
     "w=-3 chart retains the branch that escapes in the v coordinate")
gate((w+2).subs(w, -3) == -1,
     "the added w=-3 vertical address is not on the principal ramification prime")


summary = {
    "checks": CHECKS,
    "family": "F=T^3-3PT-c(t)P^2;c nonzero",
    "normalization": "B=A[w],w=cv;two smooth charts D(w+3),D(w)",
    "conductor": "B/A=R/(P);zero section splits into E0 and E3",
    "class": "if c has s distinct roots then Cl(B)=Z^s and B*=k*",
    "moving_zero": "w=-3 retains the v=-3/c branch across c=0",
    "ramification": "E2=V(w+2),reduced principal,quotient k[t,1/c]",
    "obstruction": "delete E2 makes nonconstant w+2 a forbidden A2 unit",
    "scope": "moving P^2 conductor-debt family;JC2 open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3963 moving-P2 normalization exact companion")
print(f"CHECKS={CHECKS}")
print("FAMILY=T3_MINUS_3PT_MINUS_C_OF_T_P2;C_NONZERO")
print("NORMALIZATION=B_EQUALS_A_ADJOIN_W;W_EQUALS_CV;TWO_SMOOTH_CHARTS")
print("CONDUCTOR=B_OVER_A_EQUALS_R_OVER_P;ZERO_SECTION_SPLITS_E0_PLUS_E3")
print("CLASS=C_HAS_S_ROOTS_IMPLIES_CL_Z_TO_S;UNITS_SCALAR")
print("MOVING_ZERO=W_MINUS3_RETAINS_ESCAPING_V_MINUS3_OVER_C")
print("RAMIFICATION=E2_EQUALS_V_W_PLUS2;REDUCED_PRINCIPAL;K_T_1_OVER_C")
print("NO_ATLAS=DELETE_E2;W_PLUS2_FORBIDDEN_NONCONSTANT_UNIT")
print("SCOPE=MOVING_P2_CONDUCTOR_DEBT;GENERAL_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
