#!/usr/bin/env python3
"""Exact companion for THM-3938's centered degree-four root-map closure.

Reproduction:
  python3 04-computation/jc2_centered_degree_four_root_map_thm3938.py
  python3 -O 04-computation/jc2_centered_degree_four_root_map_thm3938.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


u, A, C, T, p, q, k = sp.symbols("u A C T p q k")
h, r, z1, z2, z3, z4 = sp.symbols("h r z1 z2 z3 z4")


def incidence(A_u: sp.Expr, t: sp.Expr, scale: sp.Expr = sp.S.One):
    numerator, denominator = sp.together(t).as_numer_denom()
    resultant = sp.Poly(
        sp.expand(scale * sp.resultant(A - A_u, denominator * T - numerator, u)), T
    )
    a = sp.factor(resultant.coeff_monomial(T**3))
    trace_term = sp.factor(resultant.coeff_monomial(T**2))
    c = sp.factor(-resultant.coeff_monomial(T))
    d = sp.factor(-resultant.coeff_monomial(1) / 2)
    return resultant, a, trace_term, c, d


def color(A_u: sp.Expr, t: sp.Expr, a: sp.Expr, c: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(-(3 * a.subs(A, A_u) * t**2 + c.subs(A, A_u)) / (2 * t)))


def field_trace(A_u: sp.Expr, t: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(t).as_numer_denom()
    resultant = sp.Poly(sp.resultant(A - A_u, denominator * T - numerator, u), T)
    return sp.factor(-resultant.coeff_monomial(T**2) / resultant.coeff_monomial(T**3))


def retained(exponent: int, e: int) -> bool:
    return exponent % e == 0


# ---------------------------------------------------------------------------
# Complete local trace / Riemann--Hurwitz ledger for deg(t)=4.
# ---------------------------------------------------------------------------

gate(retained(-3, 3), "an order-three pole at polynomial infinity violates trace zero")
gate(not retained(-4, 3), "an order-four pole is compatible with e=3 trace")
gate(retained(-4, 2), "an isolated order-four finite pole cannot lie over e=2")
gate(not retained(-3, 2), "an isolated order-three finite pole can lie over e=2")
gate(retained(-3, 3), "an isolated order-three finite pole cannot lie over e=3")
gate(retained(-2, 2) and not retained(-2, 3),
     "an isolated double finite pole requires e=3")
gate(not retained(-1, 2) and not retained(-1, 3),
     "an isolated simple finite pole requires ramification")

pole_rows = {
    "E": "finite infinity;4@e3",
    "F": "finite infinity;3@e2+1@e2",
    "G": "m_inf=2;2@e3",
    "H": "m_inf=2;1@e2+1@e2",
    "I": "m_inf=1;3@e2",
}
gate(len(pole_rows) == 5, "five non-root-regular collision-free rows remain")
finite_partition_costs = {
    0: {(4,): 2, (3, 1): 2, (2, 2): 4, (2, 1, 1): 4, (1, 1, 1, 1): 4},
    1: {(3,): 1, (2, 1): 3, (1, 1, 1): 3},
    2: {(2,): 2, (1, 1): 2},
}
gate([part for part, cost in finite_partition_costs[0].items() if cost <= 2] ==
     [(4,), (3, 1)], "finite-at-infinity rows are exactly E and F")
gate([part for part, cost in finite_partition_costs[1].items() if cost <= 2] ==
     [(3,)], "m_infinity=1 leaves exactly row I")
gate([part for part, cost in finite_partition_costs[2].items() if cost <= 2] ==
     [(2,), (1, 1)], "m_infinity=2 leaves exactly rows G and H")

# The full partial-fraction traces show that no lower jet was silently lost.
trace_E_general = field_trace(u**3, h + z1/u + z2/u**2 + z3/u**3 + z4/u**4)
gate(sp.expand(trace_E_general - 3*h - 3*z3/A) == 0,
     "row E trace kills exactly the constant and u^-3 jets")
trace_F_general = field_trace(
    u**3 - 3*u,
    h + z1/(u-1) + z2/(u-1)**2 + z3/(u-1)**3 + z4/(u+1),
)
gate(sp.factor(trace_F_general - 3*(A*h + 2*z2 + z3 + 2*h)/(A+2)) == 0,
     "row F trace kills exactly h and 2*z2+z3")
trace_G_general = field_trace(u**3, r*u**2 + p*u + h + q/u + k/u**2)
gate(sp.expand(trace_G_general - 3*h) == 0,
     "row G trace kills exactly the constant jet")
trace_H_general = field_trace(
    u**3 - 3*u, r*u**2 + p*u + h + q/(u-1) + k/(u+1)
)
gate(sp.expand(trace_H_general - 6*r - 3*h) == 0,
     "row H trace forces h=-2*r")
trace_I_general = field_trace(
    u**3 + u**2, r*u + h + z1/u + z2/u**2 + z3/u**3
)
gate(sp.expand(trace_I_general + r - 3*h - (2*z2 + 3*z3)/A) == 0,
     "row I trace forces h=r/3 and 2*z2+3*z3=0")


# ---------------------------------------------------------------------------
# Row E: t(infinity) finite and one order-four pole at an e=3 point.
# ---------------------------------------------------------------------------

A_E = u**3
t_E_raw = p / u + q / u**2 + 1 / u**4
R_E, a_E, trace_E, c_E, d_E = incidence(A_E, t_E_raw)
gate(trace_E == 0 and a_E == A**4, "row E trace and primitive leading coefficient")
C_E_raw = color(A_E, t_E_raw, a_E, c_E)
N_E, D_E = sp.together(C_E_raw).as_numer_denom()
t_E_num = sp.together(t_E_raw).as_numer_denom()[0]
gate(sp.factor(sp.resultant(t_E_num, N_E, u)) == 27 * p**2 * q**6,
     "row E polynomial-color resultant")

E_p0 = sp.cancel(C_E_raw.subs(p, 0))
E_p0_num, E_p0_den = sp.together(E_p0).as_numer_denom()
E_p0_rem = sp.factor(
    sp.rem(sp.Poly(E_p0_num, u), sp.Poly(E_p0_den, u)).as_expr()
)
gate(E_p0_rem == 3 / q**4,
     "row E p=0 seam does not divide unless q=0")

t_E = sp.factor(t_E_raw.subs(q, 0))
C_E = sp.factor(C_E_raw.subs(q, 0))
gate(t_E == (p * u**3 + 1) / u**4, "row E survivor root")
gate(C_E == -sp.Rational(3, 2) * u**8 * (p * A_E + 1),
     "row E survivor color")
gate(c_E.subs(q, 0) == 0 and sp.factor(d_E.subs(q, 0)) == (p * A + 1)**3 / 2,
     "row E survivor incidence")
gate(d_E.subs({q: 0, p: 0}) == sp.Rational(1, 2),
     "row E p=0 is a literal scalar endpoint")

I_E = sp.factor(sp.resultant(A - A_E, C - C_E, u))
Phi_E = A**4 * T**3 + C * T**2 + (p * A + 1)**3 / 2
Disc_E = sp.factor(sp.discriminant(Phi_E, T))
gate(sp.expand(Disc_E + 2 * (p * A + 1)**3 * I_E) == 0,
     "row E discriminant is the implicit component to exponent one")
gate(sp.Poly(I_E, C).degree() == 3, "row E cubic A-projection")
gate(sp.Poly(C_E, u).degree() == 11, "nonzero-p row E has one degree-eleven infinity")
gate(sp.discriminant(u**3 + 1 / p, u) == -27 / p**2,
     "row E has three distinct collapsed addresses for p nonzero")


# ---------------------------------------------------------------------------
# Row F: t(infinity) finite, an order-three pole and a simple pole at e=2.
# ---------------------------------------------------------------------------

A_F = u**3 - 3 * u
t_F_raw = (
    p / (u - 1) - sp.Rational(1, 2) / (u - 1)**2
    + 1 / (u - 1)**3 + q / (u + 1)
)
R_F_raw, a_F_raw, trace_F, c_F_raw, d_F_raw = incidence(A_F, t_F_raw)
gate(trace_F == 0 and sp.expand(a_F_raw - 8 * (A - 2) * (A + 2)**3) == 0,
     "row F trace and primitive pole normalization")
C_F_raw = color(A_F, t_F_raw, a_F_raw, c_F_raw)
N_F, D_F = sp.together(C_F_raw).as_numer_denom()
t_F_num = sp.together(t_F_raw).as_numer_denom()[0]
gate(sp.factor(sp.resultant(t_F_num, N_F, u)) ==
     -927712935936 * q**3 * (p + q)**2 * (2 * p + 6 * q + 1)**6,
     "row F polynomial-color resultant")

F_sum_num = sp.factor(N_F.subs(p, -q))
F_sum_t = sp.factor(t_F_num.subs(p, -q))
gate(sp.factor(sp.prem(F_sum_num, F_sum_t, u)) ==
     -49152 * (4 * q + 1)**7 *
     (4 * q**2 * u - 4 * q**2 + 5 * q * u + q + u + 1),
     "row F p=-q seam has an exact pseudo-remainder")
gate(sp.factor(sp.gcd(F_sum_num.subs(q, -sp.Rational(1, 4)),
                      F_sum_t.subs(q, -sp.Rational(1, 4)))) == 4,
     "row F q=-1/4 is the sole degree-drop overlap of the two seams")

F_sub = {p: -3 * q - sp.Rational(1, 2)}
t_F = sp.factor(t_F_raw.subs(F_sub))
# Divide the resultant row by its harmless scalar eight, including the color.
a_F = sp.factor(a_F_raw / 8)
c_F = sp.factor(c_F_raw.subs(F_sub) / 8)
d_F = sp.factor(d_F_raw.subs(F_sub) / 8)
C_F = sp.factor(C_F_raw.subs(F_sub) / 8)
L_F = (4 * q + 1) * A + 8 * q - 2
gate(a_F == (A - 2) * (A + 2)**3, "row F survivor a")
gate(sp.expand(c_F - 3 * (2 * A - 5) * L_F**2 / 4) == 0, "row F survivor c")
gate(sp.expand(d_F + L_F**3 / 16) == 0, "row F survivor d")
gate(sp.expand(C_F - sp.Rational(3, 4) * (u - 1)**3 * (u + 1) *
     (u**4 + 6 * u**3 - 22 * u - 21) * L_F.subs(A, A_F)) == 0,
     "row F survivor color")
gate(d_F.subs(q, -sp.Rational(1, 4)) == 4,
     "row F q=-1/4 is a literal scalar endpoint")
gate(sp.factor(L_F.subs(A, 2)) == 16 * q and L_F.subs(A, -2) == -4,
     "row F coefficient ideal is a unit for the genuine simple pole q nonzero")

I_F = sp.factor(sp.resultant(A - A_F, C - C_F, u))
Phi_F = a_F * T**3 + C * T**2 + c_F * T + d_F
Disc_F = sp.factor(sp.discriminant(Phi_F, T))
gate(sp.expand(Disc_F - L_F**3 * I_F / 4) == 0,
     "row F discriminant is the implicit component to exponent one")
gate(sp.Poly(I_F, C).degree() == 3, "row F cubic A-projection")
gate(sp.Poly(C_F, u).degree() == 11, "generic row F has one degree-eleven infinity")
A0_F = sp.factor((2 - 8 * q) / (4 * q + 1))
gate(sp.factor(A0_F - 2) == -16 * q / (4 * q + 1),
     "row F collapsed target avoids the first critical value for q nonzero")
gate(sp.factor(A0_F + 2) == 4 / (4 * q + 1),
     "row F collapsed target always avoids the second critical value")


# ---------------------------------------------------------------------------
# Row G: m_infinity=2 and one order-two pole at an e=3 point.
# ---------------------------------------------------------------------------

A_G = u**3
t_G_raw = u**2 + p * u + q / u + k / u**2
R_G, a_G, trace_G, c_G, d_G = incidence(A_G, t_G_raw)
gate(trace_G == 0 and a_G == A**2, "row G trace and leading coefficient")
C_G_raw = color(A_G, t_G_raw, a_G, c_G)
N_G, D_G = sp.together(C_G_raw).as_numer_denom()
t_G_num = sp.together(t_G_raw).as_numer_denom()[0]
gate(sp.factor(sp.resultant(t_G_num, N_G, u)) == 81 * k**6 * (-k + p * q)**6,
     "row G polynomial-color resultant")

G_first = {k: p * q}
G_first_color = sp.cancel(C_G_raw.subs(G_first))
G_num, G_den = sp.together(G_first_color).as_numer_denom()
gate(sp.factor(sp.rem(sp.Poly(G_num, u), sp.Poly(G_den, u)).as_expr()) ==
     -3 * p**6 * (p**3 - q),
     "row G first seam leaves the exact final remainder")

G_sub = {q: p**3, k: p**4}
t_G = sp.factor(t_G_raw.subs(G_sub))
C_G = sp.factor(C_G_raw.subs(G_sub))
L_G = A + p**3
Q_G = u**2 - p * u + p**2
gate(t_G == (u + p)**2 * Q_G / u**2, "row G survivor root")
gate(C_G == -sp.Rational(3, 2) * u**4 * Q_G * (u**2 + 3 * p * u + p**2),
     "row G survivor color")
gate(sp.factor(c_G.subs(G_sub)) == 3 * A * p * L_G**2,
     "row G survivor c")
gate(sp.factor(d_G.subs(G_sub)) == L_G**4 / 2,
     "row G survivor d")
gate(d_G.subs({**G_sub, A: 0}) == p**12 / 2,
     "row G coefficient ideal is a unit for p nonzero")

I_G = sp.factor(sp.resultant(A - A_G, C - C_G, u))
Phi_G = A**2 * T**3 + C * T**2 + 3 * A * p * L_G**2 * T + L_G**4 / 2
Disc_G = sp.factor(sp.discriminant(Phi_G, T))
gate(sp.expand(Disc_G + 2 * L_G**4 * I_G) == 0,
     "row G discriminant is the implicit component to exponent one")
gate(sp.Poly(I_G, A, C).total_degree() == 8 and sp.Poly(I_G, C).degree() == 3,
     "row G is an implicit octic with cubic A-projection")
gate(sp.expand((u + p) * Q_G - (A_G + p**3)) == 0,
     "row G two-address factor lies over A=-p^3")
gate(sp.discriminant(Q_G, u) == -3 * p**2,
     "row G has two distinct collapsed addresses for p nonzero")


# ---------------------------------------------------------------------------
# Row H: m_infinity=2 and two simple poles over the two e=2 points.
# ---------------------------------------------------------------------------

A_H = u**3 - 3 * u
t_H_raw = u**2 + p * u - 2 + q / (u - 1) + k / (u + 1)
R_H_raw, a_H_raw, trace_H, c_H_raw, d_H_raw = incidence(A_H, t_H_raw)
gate(trace_H == 0, "row H exact trace constant is -2")
C_H_raw = color(A_H, t_H_raw, a_H_raw, c_H_raw)
N_H, D_H = sp.together(C_H_raw).as_numer_denom()
t_H_num = sp.together(t_H_raw).as_numer_denom()[0]
E_H = k * p + k + 2 * p**2 + p * q - q - 2
gate(sp.factor(sp.resultant(t_H_num, N_H, u)) == -5184 * k**3 * q**3 * E_H**6,
     "row H polynomial-color resultant")
gate(E_H.subs(p, -1) == -2 * q,
     "row H p=-1 cannot satisfy the live seam with q nonzero")

k_first = sp.factor(-(p - 1) * (2 * p + q + 2) / (p + 1))
H_first_color = sp.cancel(C_H_raw.subs(k, k_first))
H_num, H_den = sp.together(H_first_color).as_numer_denom()
H_rem = sp.factor(sp.rem(sp.Poly(H_num, u), sp.Poly(H_den, u)).as_expr())
gate(H_rem == -3 * (p - 1)**3 * (p + 1)**3 *
     (p**4 + p**3 - 3 * p**2 - 5 * p - 2 * q - 2),
     "row H first seam leaves one exact scalar compatibility")

q_H = sp.factor((p - 2) * (p + 1)**3 / 2)
k_H = sp.factor(-(p - 1)**3 * (p + 2) / 2)
H_sub = {q: q_H, k: k_H}
t_H = sp.factor(t_H_raw.subs(H_sub))
a_H = sp.factor(a_H_raw)
c_H = sp.factor(c_H_raw.subs(H_sub))
d_H = sp.factor(d_H_raw.subs(H_sub))
C_H = sp.factor(C_H_raw.subs(H_sub))
L_H = A + p**3 - 3 * p
Q_H = u**2 - p * u + p**2 - 3
gate(t_H == (u + p)**2 * Q_H / ((u - 1) * (u + 1)),
     "row H survivor root")
gate(sp.expand(a_H - (A**2 - 4)) == 0, "row H survivor a")
gate(sp.expand(c_H + 3 * L_H**2 * (-A * p + p**2 + 1)) == 0, "row H survivor c")
gate(sp.expand(d_H - L_H**4 / 2) == 0, "row H survivor d")
gate(sp.expand(C_H + sp.Rational(3, 2) * (u - 1) * (u + 1) * Q_H *
     (p**2 * u**2 - 5 * p**2 + 3 * p * u**3 - 11 * p * u
      + u**4 - 4 * u**2 - 1)) == 0,
     "row H survivor color")
gate(sp.factor(q_H) == (p - 2) * (p + 1)**3 / 2 and
     sp.factor(k_H) == -(p - 1)**3 * (p + 2) / 2,
     "row H pole cancellations occur exactly at p=2,-1,1,-2")
gate(sp.factor(L_H.subs(A, 2)) == (p - 1)**2 * (p + 2) and
     sp.factor(L_H.subs(A, -2)) == (p + 1)**2 * (p - 2),
     "row H coefficient ideal is a unit away from the pole cancellations")

I_H = sp.factor(sp.resultant(A - A_H, C - C_H, u))
Phi_H = a_H * T**3 + C * T**2 + c_H * T + d_H
Disc_H = sp.factor(sp.discriminant(Phi_H, T))
gate(sp.expand(Disc_H + 2 * L_H**4 * I_H) == 0,
     "row H discriminant is the implicit component to exponent one")
gate(sp.Poly(I_H, A, C).total_degree() == 8 and sp.Poly(I_H, C).degree() == 3,
     "row H is an implicit octic with cubic A-projection")
gate(sp.expand((u + p) * Q_H - (A_H + p**3 - 3 * p)) == 0,
     "row H two-address factor lies over A=-(p^3-3p)")
gate(sp.expand(sp.discriminant(Q_H, u) + 3 * (p - 2) * (p + 2)) == 0,
     "row H has two distinct addresses away from p=plus or minus two")
gate(sp.factor(Q_H.subs(u, 1)) == (p - 2) * (p + 1) and
     sp.factor(Q_H.subs(u, -1)) == (p - 1) * (p + 2),
     "row H addresses avoid both finite pole supports in the genuine range")


# ---------------------------------------------------------------------------
# Row I: m_infinity=1 and one order-three pole at an e=2 point.
# ---------------------------------------------------------------------------

A_I = u**3 + u**2
t_I_raw = u + sp.Rational(1, 3) + p / u - sp.Rational(3, 2) * k / u**2 + k / u**3
R_I, a_I, trace_I, c_I, d_I = incidence(A_I, t_I_raw)
gate(trace_I == 0 and a_I == 216 * A**3, "row I trace and primitive leading coefficient")
C_I_raw = color(A_I, t_I_raw, a_I, c_I)
N_I, D_I = sp.together(C_I_raw).as_numer_denom()
t_I_num = sp.together(t_I_raw).as_numer_denom()[0]
gate(sp.factor(sp.resultant(t_I_num, N_I, u)) ==
     29249267520503808 * k**9 * (3 * p + 2)**6 * (81 * k + 18 * p + 4)**3,
     "row I polynomial-color resultant")

I_first_color = sp.cancel(C_I_raw.subs(p, -sp.Rational(2, 3)))
I_num, I_den = sp.together(I_first_color).as_numer_denom()
gate(sp.expand(sp.rem(sp.Poly(I_num, u), sp.Poly(I_den, u)).as_expr() +
               1024 * (81 * k - 40) / 729) == 0,
     "row I p=-2/3 seam forces k=40/81")

I_second = {k: -(18 * p + 4) / 81}
I_second_color = sp.cancel(C_I_raw.subs(I_second))
I2_num, I2_den = sp.together(I_second_color).as_numer_denom()
I2_rem = sp.factor(sp.rem(sp.Poly(I2_num, u), sp.Poly(I2_den, u)).as_expr())
gate(sp.expand(I2_rem - 4 * (9 * p + 2)**4 *
     (81 * p**2 - 162 * p * u - 72 * p - 12 * u - 20) / 243) == 0,
     "row I second seam has the exact linear remainder")
gate(I_second[k].subs(p, -sp.Rational(2, 9)) == 0,
     "row I only vanishing prefactor on the second seam deletes the triple pole")
gate((81 * p**2 - 72 * p - 20).subs(p, -sp.Rational(2, 27)) != 0,
     "row I linear-coefficient cancellation does not cancel the constant remainder")

I_sub = {p: -sp.Rational(2, 3), k: sp.Rational(40, 81)}
t_I = sp.factor(t_I_raw.subs(I_sub))
C_I = sp.factor(C_I_raw.subs(I_sub))
a_I_s = sp.factor(a_I.subs(I_sub))
c_I_s = sp.factor(c_I.subs(I_sub))
d_I_s = sp.factor(d_I.subs(I_sub))
L_I = 27 * A - 20
Q_I = 9 * u**2 + 15 * u + 10
gate(t_I == (3 * u - 2)**2 * Q_I / (81 * u**3), "row I survivor root")
gate(C_I == -sp.Rational(4, 3) * u**3 * Q_I *
     (27 * u**5 + 45 * u**4 - 30 * u**3 - 60 * u**2 + 16),
     "row I survivor color")
gate(a_I_s == 216 * A**3, "row I survivor a")
gate(sp.expand(c_I_s + 8 * (15 * A - 4) * L_I**2 / 243) == 0, "row I survivor c")
gate(sp.expand(d_I_s - 4 * L_I**4 / 19683) == 0, "row I survivor d")
gate(d_I_s.subs(A, 0) != 0, "row I coefficient ideal is the unit ideal")

I_I = sp.factor(sp.resultant(A - A_I, C - C_I, u))
Phi_I = a_I_s * T**3 + C * T**2 + c_I_s * T + d_I_s
Disc_I = sp.factor(sp.discriminant(Phi_I, T))
gate(sp.expand(Disc_I + sp.Rational(16, 19683) * L_I**4 * I_I) == 0,
     "row I discriminant is the implicit component to exponent one")
gate(sp.Poly(I_I, A, C).total_degree() == 10 and sp.Poly(I_I, C).degree() == 3,
     "row I is an implicit decic with cubic A-projection")
gate(sp.expand((3 * u - 2) * Q_I - (27 * A_I - 20)) == 0,
     "row I two-address factor lies over A=20/27")
gate(sp.discriminant(Q_I, u) == -135,
     "row I has two distinct collapsed addresses")


# ---------------------------------------------------------------------------
# Uniform birationality / genuine-branch controls.
# ---------------------------------------------------------------------------

for label, implicit, pole_order, ramification_index in [
    ("E", I_E, 4, 3),
    ("F", I_F, 3, 2),
    ("G", I_G, 2, 3),
    ("H", I_H, 1, 2),
    ("I", I_I, 3, 2),
]:
    gate(pole_order % ramification_index != 0,
         f"row {label} repeated root is not in k(A) by its pole valuation")
    gate(sp.Poly(implicit, C).degree() == 3,
         f"row {label} implicit relation has no parametrization power")


summary = {
    "checks": CHECKS,
    "scope": "centered trace-zero linear-color;degA=3;degt=4",
    "pole_rows": "E:4e3;F:3e2+1e2;G:m2+2e3;H:m2+1e2+1e2;I:m1+3e2",
    "row_E": "q=0;p=0 scalar, else three addresses",
    "row_F": "p=-3q-1/2;q=-1/4 scalar, else three addresses",
    "row_G": "q=p^3,k=p^4;two addresses",
    "row_H": "q=(p-2)(p+1)^3/2,k=-(p-1)^3(p+2)/2;two addresses",
    "row_I": "p=-2/3,k=40/81;two addresses",
    "branch": "each nonscalar implicit component has discriminant exponent one",
    "conclusion": "scalar-unit or source ramification non-unibranch;no A2 Keller open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3938 centered degree-four root-map exact companion")
print(f"CHECKS={CHECKS}")
print("SCOPE=centered trace-zero linear-color;degA=3;deg_t=4")
print("POLE_ROWS=E:4e3;F:3e2+1e2;G:m2+2e3;H:m2+1e2+1e2;I:m1+3e2")
print("ROW_E=q=0;p=0 scalar,else three addresses")
print("ROW_F=p=-3q-1/2;q=-1/4 scalar,else three addresses")
print("ROW_G=q=p^3,k=p^4;two addresses")
print("ROW_H=q=(p-2)(p+1)^3/2;k=-(p-1)^3(p+2)/2;two addresses")
print("ROW_I=p=-2/3;k=40/81;two addresses")
print("BRANCH=each nonscalar implicit component has discriminant exponent one")
print("CONCLUSION=scalar-unit or source ramification non-unibranch;no A2 Keller open")
print(f"SEMANTIC_SHA256={semantic}")
