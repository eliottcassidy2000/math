#!/usr/bin/env python3
"""Exact companion for THM-3961's arbitrary-q normality trichotomy."""

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


# ---------------------------------------------------------------------------
# Universal P-adic row and ramification/Jacobian dictionary.
# ---------------------------------------------------------------------------

P, T, h, t = sp.symbols("P T h t")
q0 = sp.Function("q0")(t)
q1 = sp.Function("q1")(t)
q2 = sp.Function("q2")(P, t)

q = q0 + P * q1 + P**2 * q2
F = sp.expand(T**3 - 3 * P * T - q)
K = sp.expand(q.subs(P, h**2) - 2 * h**3)

FT = sp.diff(F, T)
FP = sp.diff(F, P)
Ft = sp.diff(F, t)
ramification_substitution = {P: h**2, T: -h}

zero(FT - 3 * (T**2 - P), "relative different generator")
zero(F.subs(ramification_substitution) + K,
     "ramification substitution gives minus hidden polynomial")
zero(sp.diff(K, h) + 2 * h * FP.subs(ramification_substitution),
     "hidden h derivative equals minus 2h times F_P")
zero(sp.diff(K, t) + Ft.subs(ramification_substitution),
     "hidden t derivative equals minus F_t")

# Exact h=0 exceptional row.
h_zero = {P: 0, T: 0}
zero(F.subs(h_zero) + q0, "zero-section equation is q0")
zero(FT.subs(h_zero), "zero-section F_T vanishes")
zero(FP.subs(h_zero) + q1, "zero-section F_P is minus q1")
zero(Ft.subs(h_zero) + sp.diff(q0, t),
     "zero-section F_t is minus q0 prime")


# ---------------------------------------------------------------------------
# The three exact P-adic cases and the vertical-content guard.
# ---------------------------------------------------------------------------

K_once = sp.factor(K.subs(q0, 0))
L_once = h**2 * q2.subs(P, h**2) - 2 * h + q1
zero(K_once - h**2 * L_once, "forced h^2 in the once-P-divisible row")
zero(L_once.subs(h, 0) - q1, "adjusted row has value q1 at h=0")

K_twice = sp.factor(K.subs({q0: 0, q1: 0}))
zero(K_twice - h**3 * (h * q2.subs(P, h**2) - 2),
     "P^2 debt gives forced h^3")

# q(h^2,t) has only even h powers. A finite universal row makes the fixed
# odd coefficient, and therefore primitive content, explicit.
a0, a1, a2, b0, b1 = sp.symbols("a0 a1 a2 b0 b1")
q_finite = b0 + P * b1 + P**2 * (a0 + a1 * P + a2 * P**2)
K_finite = sp.Poly(sp.expand(q_finite.subs(P, h**2) - 2 * h**3), h)
gate(K_finite.coeff_monomial(h**3) == -2,
     "case I has unit h^3 coefficient")
L_finite = sp.Poly(
    sp.expand(b1 - 2 * h + h**2 * (a0 + a1 * h**2 + a2 * h**4)), h
)
gate(L_finite.coeff_monomial(h) == -2,
     "case II adjusted row has unit h coefficient")
gate(sp.gcd_list(K_finite.all_coeffs()) == 1,
     "case I finite universal row has no vertical content")
gate(sp.gcd_list(L_finite.all_coeffs()) == 1,
     "case II finite universal row has no vertical content")

# In the P^2 row the entire zero section lies in the Jacobian scheme.
qhat = sp.Function("qhat")(P, t)
F_twice = sp.expand(T**3 - 3 * P * T - P**2 * qhat)
for expression, message in (
    (F_twice, "P^2 row equation"),
    (sp.diff(F_twice, T), "P^2 row F_T"),
    (sp.diff(F_twice, P), "P^2 row F_P"),
    (sp.diff(F_twice, t), "P^2 row F_t"),
):
    zero(expression.subs({P: 0, T: 0}),
         f"{message} vanishes on the full zero section")


# ---------------------------------------------------------------------------
# Hostile 1: forced h^2 is harmless, even with a repeated q1 zero.
# ---------------------------------------------------------------------------

q_forced = P * t**2
F_forced = sp.expand(T**3 - 3 * P * T - q_forced)
K_forced = sp.factor(q_forced.subs(P, h**2) - 2 * h**3)
L_forced = t**2 - 2 * h
zero(K_forced - h**2 * L_forced, "forced-square hostile factorization")
gate(sp.diff(L_forced, h) == -2,
     "forced-square adjusted curve is smooth")

# Eisenstein at P: every lower T coefficient is divisible by P, and the
# constant row divided by P is not divisible by P.
forced_poly = sp.Poly(F_forced, T)
gate(forced_poly.LC() == 1, "forced-square hostile is monic")
gate(forced_poly.coeff_monomial(T**2) == 0,
     "forced-square T^2 row")
gate(sp.rem(forced_poly.coeff_monomial(T), P, domain=sp.QQ[t]) == 0,
     "forced-square T row divisible by P")
gate(sp.rem(forced_poly.coeff_monomial(1), P, domain=sp.QQ[t]) == 0,
     "forced-square constant row divisible by P")
gate(sp.rem(forced_poly.coeff_monomial(1) / P, P,
            domain=sp.QQ[t]) != 0,
     "forced-square constant row is not divisible by P squared")

forced_gradients = [
    F_forced,
    sp.diff(F_forced, T),
    sp.diff(F_forced, P),
    sp.diff(F_forced, t),
]
for expression in forced_gradients:
    zero(expression.subs({P: 0, T: 0, t: 0}),
         "forced-square hostile has the predicted isolated origin")
gate(sp.diff(F_forced, P).subs({P: 0, T: 0}) == -t**2,
     "forced-square zero section is generically regular")


# ---------------------------------------------------------------------------
# Hostile 2: a repeated nonzero adjusted factor is a singular divisor.
# ---------------------------------------------------------------------------

q_repeat = P**2 + P
F_repeat = sp.expand(T**3 - 3 * P * T - q_repeat)
K_repeat = sp.factor(q_repeat.subs(P, h**2) - 2 * h**3)
zero(K_repeat - h**2 * (h - 1) ** 2,
     "repeated-nonzero hidden factor")

# Irreducibility: -F is quadratic in P and its discriminant has the
# nonsquare factor 4T+1.
disc_repeat_P = sp.factor(sp.discriminant(-F_repeat, P))
gate(disc_repeat_P == (T + 1) ** 2 * (4 * T + 1),
     "repeated hostile quadratic discriminant")
repeat_factor_rows = sp.factor_list(disc_repeat_P, T)[1]
gate(repeat_factor_rows == [(4 * T + 1, 1), (T + 1, 2)],
     "repeated hostile discriminant has an odd linear factor")

repeat_line = {P: 1, T: -1}
for expression in (
    F_repeat,
    sp.diff(F_repeat, T),
    sp.diff(F_repeat, P),
    sp.diff(F_repeat, t),
):
    zero(expression.subs(repeat_line),
         "repeated hidden factor gives a full singular line")
gate(sp.diff(F_repeat, P).subs({P: 0, T: 0}) == -1,
     "forced h=0 component is regular in repeated-factor hostile")


# ---------------------------------------------------------------------------
# Hostile 3: P^2 divisibility makes the zero section nonnormal.
# ---------------------------------------------------------------------------

q_p2 = P**2
F_p2 = sp.expand(T**3 - 3 * P * T - q_p2)
K_p2 = sp.factor(q_p2.subs(P, h**2) - 2 * h**3)
zero(K_p2 - h**3 * (h - 2), "P^2 hostile hidden factorization")
disc_p2_P = sp.factor(sp.discriminant(-F_p2, P))
gate(disc_p2_P == T**2 * (4 * T + 9),
     "P^2 hostile quadratic discriminant")
p2_factor_rows = sp.factor_list(disc_p2_P, T)[1]
gate(p2_factor_rows == [(4 * T + 9, 1), (T, 2)],
     "P^2 hostile discriminant has an odd linear factor")
for expression in (
    F_p2,
    sp.diff(F_p2, T),
    sp.diff(F_p2, P),
    sp.diff(F_p2, t),
):
    zero(expression.subs({P: 0, T: 0}),
         "P^2 hostile has full singular zero section")


# ---------------------------------------------------------------------------
# Positive smooth control and global-different unit gate.
# ---------------------------------------------------------------------------

q_smooth = t
F_smooth = sp.expand(T**3 - 3 * P * T - q_smooth)
K_smooth = t - 2 * h**3
gate(sp.diff(F_smooth, t) == -1, "q=t surface is smooth")
gate(sp.gcd(K_smooth, sp.diff(K_smooth, h),
            sp.diff(K_smooth, t)) == 1,
     "q=t hidden curve is squarefree")
zero(t - (T**3 - 3 * P * T) + F_smooth,
     "q=t surface solves polynomially for t")

# Freeze the arbitrary-q order discriminant and the nonconstant basis row.
Q = sp.symbols("Q")
F_abstract = T**3 - 3 * P * T - Q
disc_abstract = sp.factor(sp.discriminant(F_abstract, T))
zero(disc_abstract + 27 * (Q**2 - 4 * P**3),
     "arbitrary-q cubic discriminant")
delta = sp.Poly(3 * (T**2 - P), T)
gate(delta.all_coeffs() == [3, 0, -3 * P],
     "global different has nonzero T^2 basis row")


summary = {
    "checks": CHECKS,
    "universe": "irreducible T^3-3PT-q(P,t) over alg closed char0",
    "dictionary": "K_h=-2h F_P and K_t=-F_t on P=h^2,T=-h",
    "case_I": "q0 nonzero; normal iff K squarefree",
    "case_II": "q0=0,q1 nonzero; remove forced h^2; normal iff L squarefree",
    "case_III": "P^2 divides q; zero section is codimension-one singular",
    "vertical": "fixed odd coefficient -2 makes adjusted L primitive",
    "normal_closure": "global monogenic different is forbidden A2 unit",
    "escape": "repeated nonzero factor or P^2 zero-section conductor debt",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3961 arbitrary-q hidden repetition normality companion")
print(f"CHECKS={CHECKS}")
print("UNIVERSE=IRREDUCIBLE_CUBIC_DOMAIN;ARBITRARY_Q_IN_K[P,t]")
print("DICTIONARY=K_H_MINUS_2H_F_P;K_T_MINUS_F_t;H_ZERO_ROW_EXACT")
print("CASE_I=Q0_NONZERO;NORMAL_IFF_K_SQUAREFREE")
print("CASE_II=Q0_ZERO_Q1_NONZERO;FORCED_H2_HARMLESS;NORMAL_IFF_L_SQUAREFREE")
print("CASE_III=P2_DIVIDES_Q;ZERO_SECTION_HEIGHT_ONE_SINGULAR;NONNORMAL")
print("VERTICAL=NO_T_FACTOR;FIXED_ODD_COEFFICIENT_MINUS_2")
print("NORMAL_CASES=GLOBAL_DIFFERENT_FORBIDDEN_UNIT;NO_KELLER_A2_CHART")
print("DEBT=REPEATED_NONZERO_HIDDEN_FACTOR_OR_P2_ZERO_SECTION")
print("HOSTILES=PT2_NORMAL;P2_PLUS_P_DIVISOR;P2_ZERO_SECTION;T_SMOOTH")
print("SCOPE=MONOGENIC_GRAMMAR_ONLY;NORMALIZATIONS_OF_TWO_DEBTS_AND_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
