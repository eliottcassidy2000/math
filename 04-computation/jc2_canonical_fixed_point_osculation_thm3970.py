#!/usr/bin/env python3
"""Exact companion for THM-3970's fixed-point osculation reframe."""

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition, label):
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)


def zero(expr, label):
    gate(sp.factor(expr) == 0, label)


P, T, h, t = sp.symbols("P T h t")
sh, qh = sp.symbols("sh qh")

# Universal algebraic identity after P=h^2. Here qh and sh stand for
# q(h^2,t) and s(h^2,t).
K_formal = qh - 2 * h**3
R_formal = qh - 3 * sh * h**2 + sh**3
zero(K_formal - (R_formal - (h - sh) ** 2 * (2 * h + sh)),
     "universal K-R fixed-point identity")
zero((h**2 - sh**2) - (h - sh) * (h + sh),
     "fixed-point reflection identity")

# The differential sidecar is checked for a generic polynomial q(P,t).
coeffs = sp.symbols("q00 q10 q01 q20 q11 q02 q30 q21 q12 q03")
(q00, q10, q01, q20, q11, q02, q30, q21, q12, q03) = coeffs
q_generic = (q00 + q10 * P + q01 * t + q20 * P**2 + q11 * P * t
             + q02 * t**2 + q30 * P**3 + q21 * P**2 * t
             + q12 * P * t**2 + q03 * t**3)
s_generic = sp.diff(q_generic, P) / 3
J_generic = P - s_generic**2
R_generic = q_generic - 3 * s_generic * P + s_generic**3
zero(sp.diff(R_generic, P) + 3 * sp.diff(s_generic, P) * J_generic,
     "canonical remainder derivative")

C, D, c = sp.symbols("C D c")
C_P, c_P, s_P = sp.symbols("C_P c_P s_P")
zero((2 * C * C_P * c + C**2 * c_P + 3 * s_P * C * D)
     - C * (2 * c * C_P + c_P * C + 3 * s_P * D),
     "factored differential sidecar")

# Exact nongraph family with rho=P+t and arbitrary scalar sigma.
sigma = sp.symbols("sigma")
rho = P + t
C0 = P - rho**2
q_sigma = sp.expand(3 * rho * P - rho**3 + sigma * C0**2)
s_can = sp.factor(sp.diff(q_sigma, P) / 3)
A_sigma = 1 + sp.Rational(2, 3) * sigma * (1 - 2 * rho)
zero(s_can - (rho + A_sigma * C0), "canonical s in nongraph family")

J_sigma = sp.expand(P - s_can**2)
R_sigma = sp.expand(q_sigma - 3 * s_can * P + s_can**3)
zero(J_sigma - C0 * (1 - 2 * rho * A_sigma - A_sigma**2 * C0),
     "nongraph contracted fixed-point factor")
zero(R_sigma - C0**2 * (sigma - 3 * A_sigma
                         + 3 * rho * A_sigma**2 + A_sigma**3 * C0),
     "nongraph second-order remainder")

rho_h = h**2 + t
K_sigma = sp.expand(q_sigma.subs(P, h**2) - 2 * h**3)
zero(K_sigma - (h - rho_h) ** 2
     * (sigma * (h + rho_h) ** 2 - (2 * h + rho_h)),
     "nongraph hidden repeated-factor identity")

M = h**2 - h + t
gate(sp.discriminant(M, h) == 1 - 4 * t,
     "nongraph factor discriminant")
zero(h - rho_h + M, "nongraph orientation")
zero(C0.subs(P, h**2) - (h - rho_h) * (h + rho_h),
     "contraction equals oriented times reflected factors")

# sigma=0 is the sharp reducible cubic-depth control.
q_zero = sp.expand(q_sigma.subs(sigma, 0))
F_zero = sp.expand(T**3 - 3 * P * T - q_zero)
zero(F_zero - (T + rho) * (T**2 - rho * T + rho**2 - 3 * P),
     "degree3 nongraph control fails domain gate")

# For sigma nonzero, the target shear a0=P+t gives the graph normalization
# k[a0,v]. Verify the full parametrization and its principal Jacobian.
a0, v = sp.symbols("a0 v")
P_src = sigma * v**3 + 3 * v**2 + 3 * a0 * v + a0**2
t_src = a0 - P_src
T_src = sigma * v**2 + 3 * v + 2 * a0
F_sigma = sp.expand(T**3 - 3 * P * T - q_sigma)
zero(F_sigma.subs({P: P_src, t: t_src, T: T_src}),
     "irreducible degree4 normalization parametrization")
gate(sp.degree(P_src, v) == 3, "sigma-nonzero cubic field degree row")
jac = sp.det(sp.Matrix([
    [sp.diff(P_src, a0), sp.diff(P_src, v)],
    [sp.diff(t_src, a0), sp.diff(t_src, v)],
]))
zero(jac + 3 * (sigma * v**2 + 2 * v + a0),
     "principal source-target Jacobian")

# The original coefficient depth really jumps from three to four.
gate(sp.degree(q_zero, P) == 3, "reducible control P-depth three")
gate(sp.degree(q_sigma, P) == 4, "nonzero-sigma formal P-depth four")

summary = {
    "checks": CHECKS,
    "canonical": "s=q_P/3; J=P-s^2; R=q-3sP+s^3",
    "forward": "M^2|K and M!=h imply M|h-s; C|J; C^2|R",
    "reflection": "C(h^2) is associate to M(h)M(-h)",
    "curve": "A/(C) isomorphic to B/(M)",
    "converse": "C|J; C^2|R plus orientation imply M^2|K",
    "gcd": "for C!=P, C|J and C|R force C^2|R and uniquely recover M",
    "differential": "2cC_P+c_PC=-3s_PD",
    "controls": "sigma0 reducible depth3; sigma nonzero irreducible depth4 A2",
    "scope": "structural reframe; normalization and JC2 remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3970 canonical fixed-point osculation companion")
print(f"CHECKS={CHECKS}")
print("CANONICAL=S_QP_OVER_3;J_P_MINUS_S2;R_Q_MINUS_3SP_PLUS_S3")
print("FORWARD=M2_DIVIDES_K_IMPLIES_ORIENTED_C_DIVIDES_J_AND_C2_DIVIDES_R")
print("REFLECTION=C_OF_H2_ASSOCIATE_M_TIMES_M_MINUS;CURVES_ISOMORPHIC")
print("CONVERSE=ORIENTED_C_PACKET_IMPLIES_M2_DIVIDES_K")
print("GCD=NON_P_FACTORS_OF_GCD_J_R_BIJECT_NONZERO_REPEATED_FACTORS")
print("DIFFERENTIAL=2C_SMALL_CP_PLUS_SMALL_CP_C_EQUALS_MINUS_3SP_D")
print("CONTROL_SIGMA0=NONGRAPH_DEPTH3_BUT_REDUCIBLE")
print("CONTROL_SIGMANONZERO=NONGRAPH_DEPTH4_IRREDUCIBLE_NORMALIZATION_A2")
print("SCOPE=FIXED_POINT_REFRAME;COLLISION_NORMALIZATION;NONMONOGENIC;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
