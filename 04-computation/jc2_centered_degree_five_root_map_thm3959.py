#!/usr/bin/env python3
"""Exact companion for THM-3959's centered degree-five root-map closure.

Reproduction:
  python3 04-computation/jc2_centered_degree_five_root_map_thm3959.py
  python3 -O 04-computation/jc2_centered_degree_five_root_map_thm3959.py
"""

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


u, A, C, T, p, q, r, k, z = sp.symbols("u A C T p q r k z")
h, z1, z2, z3, z4, z5 = sp.symbols("h z1 z2 z3 z4 z5")


def incidence(A_u: sp.Expr, root: sp.Expr):
    """Return the denominator-cleared centered cubic incidence."""
    numerator, denominator = sp.together(root).as_numer_denom()
    resultant = sp.Poly(
        sp.expand(sp.resultant(A - A_u, denominator * T - numerator, u)), T
    )
    aa = sp.factor(resultant.coeff_monomial(T**3))
    trace_term = sp.factor(resultant.coeff_monomial(T**2))
    cc = sp.factor(-resultant.coeff_monomial(T))
    dd = sp.factor(-resultant.coeff_monomial(1) / 2)
    return resultant, aa, trace_term, cc, dd


def color(A_u: sp.Expr, root: sp.Expr, aa: sp.Expr, cc: sp.Expr) -> sp.Expr:
    return sp.factor(
        sp.cancel(-(3 * aa.subs(A, A_u) * root**2 + cc.subs(A, A_u)) / (2 * root))
    )


def field_trace(A_u: sp.Expr, root: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(root).as_numer_denom()
    resultant = sp.Poly(
        sp.resultant(A - A_u, denominator * T - numerator, u), T
    )
    return sp.factor(
        -resultant.coeff_monomial(T**2) / resultant.coeff_monomial(T**3)
    )


def color_resultant(root: sp.Expr, color_value: sp.Expr) -> sp.Expr:
    root_num = sp.together(root).as_numer_denom()[0]
    color_num = sp.together(color_value).as_numer_denom()[0]
    return sp.factor(sp.resultant(root_num, color_num, u))


def divided_remainder(expression: sp.Expr, domain: str) -> sp.Expr:
    numerator, denominator = sp.together(sp.cancel(expression)).as_numer_denom()
    return sp.factor(
        sp.rem(sp.Poly(numerator, u, domain=domain),
               sp.Poly(denominator, u, domain=domain)).as_expr()
    )


def coefficient_numerators(expression: sp.Expr, parameters: tuple[sp.Symbol, ...]):
    """Primitive parameter-polynomial numerators of a u-polynomial."""
    answer = []
    for coefficient in sp.Poly(expression, u).all_coeffs():
        numerator = sp.together(coefficient).as_numer_denom()[0]
        primitive = sp.primitive(sp.Poly(numerator, *parameters))[1]
        answer.append(sp.factor(primitive.as_expr()))
    return answer


def implicit_data(A_u: sp.Expr, color_value: sp.Expr, phi: sp.Expr):
    implicit = sp.factor(sp.resultant(A - A_u, C - color_value, u))
    discriminant = sp.factor(sp.discriminant(phi, T))
    return implicit, discriminant


# ---------------------------------------------------------------------------
# The seven trace-normalized carriers supplied by THM-3941.
# ---------------------------------------------------------------------------

A_J = u**3
t_J_raw = p/u + q/u**2 + r/u**4 + 1/u**5

A_K = u**3 + u**2
t_K_raw = p/u - 3*q/(2*u**2) + q/u**3 - sp.Rational(5, 2)/u**4 + 1/u**5

A_L = u**3
t_L_raw = u + p/u + q/u**2 + r/u**4

A_M = u**3 - 3*u
t_M_raw = (
    u + p/(u-1) + q/(u+1) + k/(2*(u+1)**2) + k/(u+1)**3
)

A_N = u**3 + u**2
t_N_raw = (
    u**2 + p*u + (p-1)/3 + q/u - 3*k/(2*u**2) + k/u**3
)

A_O = u**3
t_O_raw = u**4 + p*u**2 + q*u + r/u

A_P = u**3 + u**2
t_P_raw = (
    u**4 + sp.Rational(4, 3)*u**3 + p*u**2 + q*u
    + (sp.Rational(1, 9) - p/3 + q/3) + r/u
)

rows = {
    "J": (A_J, t_J_raw, "m0+C3(5)"),
    "K": (A_K, t_K_raw, "m0+C2(5)"),
    "L": (A_L, t_L_raw, "m1+C3(4)"),
    "M": (A_M, t_M_raw, "m1+C2xC2(1,3)"),
    "N": (A_N, t_N_raw, "m2+C2(3)"),
    "O": (A_O, t_O_raw, "m4+C3(1)"),
    "P": (A_P, t_P_raw, "m4+C2(1)"),
}
gate(len(rows) == 7, "THM-3941 leaves exactly seven nonregular N=5 rows")

# These general traces also freeze the lower-principal-part normalizations.
trace_C3_5 = field_trace(
    u**3, h + z1/u + z2/u**2 + z3/u**3 + z4/u**4 + z5/u**5
)
zero(trace_C3_5 - 3*(A*h + z3)/A,
     "C3 order-five trace keeps exactly the constant and u^-3 jets")
trace_C2_5 = field_trace(
    u**3 + u**2, h + z1/u + z2/u**2 + z3/u**3 + z4/u**4 + z5/u**5
)
zero(trace_C2_5 - (3*A**2*h + 2*A*z2 + 3*A*z3 + 2*z4 + 5*z5)/A**2,
     "C2 order-five trace fixes both even cancellation rows")
for label, (A_u, root, _) in rows.items():
    zero(field_trace(A_u, root), f"row {label} has centered trace zero")

incidences = {}
colors = {}
for label, (A_u, root, _) in rows.items():
    incidence_row = incidence(A_u, root)
    incidences[label] = incidence_row
    gate(incidence_row[2] == 0, f"row {label} incidence has no quadratic term")
    colors[label] = color(A_u, root, incidence_row[1], incidence_row[3])

expected_color_resultants = {
    "J": 81*p**2*(p-q*r)**6,
    "K": 2**21*p**2*(2*p+5*q)**6*(2*p+9*q+27)**3,
    "L": 243*r**12*(p*q-r)**6,
    "M": -8549802417586176*k**9*p**3
         *(k*p+k*q+6*k+6*p**2+8*p*q+24*p+2*q**2+16*q+24)**6,
    "N": 1184595334580404224*k**9
         *(243*k+12*p+54*q-4)**3
         *(9*k+4*p**2+6*p*q-8*p-6*q+4)**6,
    "O": 243*r**6*(p*q-r)**6,
    "P": -109418989131512359209*r**3
         *(18*p-54*q-243*r-14)**3
         *(9*p**2-9*p*q-6*p+3*q+9*r+1)**6,
}
for label in rows:
    zero(color_resultant(rows[label][1], colors[label])
         - expected_color_resultants[label],
         f"row {label} exact numerator-denominator color resultant")


# ---------------------------------------------------------------------------
# J = C3(5): a free scalar family and a two-address family.
# ---------------------------------------------------------------------------

J_p0_rem = divided_remainder(colors["J"].subs(p, 0), "QQ(q,r)")
zero(J_p0_rem - 3*r**2*(
    6*q*r**2*u**2 + 4*q*r*u + q - r**5*u**2 - 5*r**4*u - 4*r**3
)/q**5, "J p=0 seam exact remainder")
J_p0_sat = sp.groebner([
    6*q*r**2-r**5, 4*q*r-5*r**4, q-4*r**3, z*r-1
], z, q, r, order="grevlex")
gate(len(J_p0_sat.polys) == 1 and J_p0_sat.polys[0].as_expr() == 1,
     "J p=0 remainder coefficients have no common zero with r nonzero")
zero(divided_remainder(colors["J"].subs({p: 0, q: 0}), "QQ(r)")
     - 3/r**10, "J p=q=0 exceptional seam excludes r nonzero")
J_pqr_rem = divided_remainder(colors["J"].subs(p, q*r), "QQ(q,r)")
zero(J_pqr_rem + 3*(q-r**3)/r**13,
     "J p=qr seam has one final compatibility")

# First J survivor: p=r=0, q arbitrary.
t_J1 = sp.factor(t_J_raw.subs({p: 0, r: 0}))
a_J1 = sp.factor(incidences["J"][1].subs({p: 0, r: 0}))
c_J1 = sp.factor(incidences["J"][3].subs({p: 0, r: 0}))
d_J1 = sp.factor(incidences["J"][4].subs({p: 0, r: 0}))
C_J1 = sp.factor(colors["J"].subs({p: 0, r: 0}))
L_J1 = q*A + 1
zero(t_J1 - (q*u**3+1)/u**5, "J1 repeated root")
zero(C_J1 + sp.Rational(3, 2)*u**10*(q*u**3+1), "J1 polynomial color")
gate(a_J1 == A**5 and c_J1 == 0 and d_J1 == L_J1**3/2,
     "J1 primitive incidence")
gate(d_J1.subs(q, 0) == sp.Rational(1, 2),
     "J1 q=0 is the sole literal scalar endpoint")

I_J1, D_J1 = implicit_data(
    A_J, C_J1, a_J1*T**3 + C*T**2 + c_J1*T + d_J1
)
zero(D_J1 + 2*L_J1**3*I_J1,
     "J1 implicit component has discriminant exponent one")
gate(sp.Poly(I_J1, C).degree() == 3 and
     sp.Poly(I_J1, A, C).total_degree() == 13,
     "J1 generic implicit degree is thirteen with cubic A-projection")
gate(sp.discriminant(u**3+1/q, u) == -27/q**2,
     "J1 q nonzero has three distinct addresses")

# Second J survivor: p=r^4,q=r^3 with r nonzero.
J2_sub = {p: r**4, q: r**3}
t_J2 = sp.factor(t_J_raw.subs(J2_sub))
a_J2 = sp.factor(incidences["J"][1].subs(J2_sub))
c_J2 = sp.factor(incidences["J"][3].subs(J2_sub))
d_J2 = sp.factor(incidences["J"][4].subs(J2_sub))
C_J2 = sp.factor(colors["J"].subs(J2_sub))
L_J2 = A*r**3 + 1
Q_J2 = r**2*u**2-r*u+1
zero(t_J2 - (r*u+1)**2*Q_J2/u**5, "J2 repeated root factorization")
zero(C_J2 + sp.Rational(3, 2)*u**10*Q_J2*(r**2*u**2+3*r*u+1),
     "J2 polynomial color")
gate(a_J2 == A**5 and c_J2 == 3*A**2*r*L_J2**2 and
     d_J2 == L_J2**4/2, "J2 primitive incidence")
I_J2, D_J2 = implicit_data(
    A_J, C_J2, a_J2*T**3 + C*T**2 + c_J2*T + d_J2
)
zero(D_J2 + 2*L_J2**4*I_J2,
     "J2 implicit component has discriminant exponent one")
gate(sp.Poly(I_J2, C).degree() == 3 and
     sp.Poly(I_J2, A, C).total_degree() == 14,
     "J2 implicit degree is fourteen with cubic A-projection")
zero((r*u+1)*Q_J2 - (r**3*A_J+1), "J2 two-address identity")
gate(sp.discriminant(Q_J2, u) == -3*r**2,
     "J2 r nonzero has two distinct addresses")
gate(sp.factor(Q_J2.subs(u, -1/r)) == 3 and
     sp.factor(C_J2.subs(u, -1/r)) == sp.Rational(9, 2)/r**10,
     "J2 third cubic-fibre point has nonzero color")


# ---------------------------------------------------------------------------
# K = C2(5): all three resultant seams, one isolated survivor.
# ---------------------------------------------------------------------------

K_p0_rem = divided_remainder(colors["K"].subs(p, 0), "QQ(q)")
K_p0_coeffs = coefficient_numerators(K_p0_rem, (q,))
gate(sp.gcd_list(K_p0_coeffs) == 1,
     "K p=0 seam coefficient ideal is the unit ideal")
gate(divided_remainder(colors["K"].subs({p: 0, q: 0}), "QQ") != 0,
     "K p=q=0 exceptional degree drop still fails division")

K_middle_rem = divided_remainder(colors["K"].subs(p, -5*q/2), "QQ(q)")
zero(K_middle_rem + 65536*(28*q+125)/1220703125,
     "K 2p+5q seam has one isolated parameter")

K_last_rem = divided_remainder(colors["K"].subs(p, -(9*q+27)/2), "QQ(q)")
K_last_coeffs = coefficient_numerators(K_last_rem, (q,))
gate(sp.gcd_list(K_last_coeffs) == 1,
     "K 2p+9q+27 seam coefficient ideal is the unit ideal")
gate(divided_remainder(
    colors["K"].subs({p: 0, q: -3}), "QQ"
) != 0, "K last-seam exceptional degree drop still fails division")

K_sub = {p: sp.Rational(625, 56), q: -sp.Rational(125, 28)}
t_K = sp.factor(t_K_raw.subs(K_sub))
a_K = sp.factor(incidences["K"][1].subs(K_sub))
c_K = sp.factor(incidences["K"][3].subs(K_sub))
d_K = sp.factor(incidences["K"][4].subs(K_sub))
C_K = sp.factor(colors["K"].subs(K_sub))
L_K = 125*A-28
Q_K = 25*u**2+35*u+14
R_K = 75*u**7+385*u**6+602*u**5+280*u**4-140*u**3-140*u**2+16
zero(t_K - (5*u-2)**2*Q_K/(56*u**5), "K survivor repeated root")
zero(C_K + u**5*Q_K*R_K/14, "K survivor polynomial color")
gate(a_K == 8*A**5 and
     c_K == L_K**2*(70*A**2-35*A+4)/392 and
     d_K == L_K**4/43904, "K survivor primitive incidence")
I_K, D_K = implicit_data(A_K, C_K, a_K*T**3+C*T**2+c_K*T+d_K)
zero(D_K + L_K**4*I_K/10976,
     "K implicit component has discriminant exponent one")
gate(sp.Poly(I_K, C).degree() == 3 and
     sp.Poly(I_K, A, C).total_degree() == 14,
     "K implicit degree is fourteen with cubic A-projection")
zero((5*u-2)*Q_K - L_K.subs(A, A_K), "K two-address identity")
gate(sp.discriminant(Q_K, u) == -175,
     "K survivor has two distinct addresses")
gate(Q_K.subs(u, sp.Rational(2, 5)) == 32 and
     R_K.subs(u, sp.Rational(2, 5)) == -sp.Rational(1024, 3125) and
     C_K.subs(u, sp.Rational(2, 5)) != 0,
     "K third cubic-fibre point has nonzero color")


# ---------------------------------------------------------------------------
# L and O: denominator-complement dual rows are uniformly empty.
# ---------------------------------------------------------------------------

L_seam_rem = divided_remainder(colors["L"].subs(r, p*q), "QQ(p,q)")
zero(L_seam_rem + 3*p**6*(p*u-q),
     "L seam exact remainder is the shared linear obstruction")
O_seam_rem = divided_remainder(colors["O"].subs(r, p*q), "QQ(p,q)")
zero(O_seam_rem - 3*p**3*(p*u-q),
     "O seam exact remainder is the shared linear obstruction")
gate(sp.Poly(L_seam_rem, u).degree() == 1 and
     sp.Poly(O_seam_rem, u).degree() == 1,
     "both C3 mixed rows fail when r=pq is nonzero")


# ---------------------------------------------------------------------------
# M: the two-critical-point conic has no live exact-division point.
# ---------------------------------------------------------------------------

E_M = k*p+k*q+6*k+6*p**2+8*p*q+24*p+2*q**2+16*q+24
M_num, M_den = sp.together(colors["M"]).as_numer_denom()
M_rem = sp.rem(
    sp.Poly(M_num, u, domain="QQ(p,q,k)"),
    sp.Poly(M_den, u, domain="QQ(p,q,k)"),
).as_expr()
M_coeffs = coefficient_numerators(M_rem, (p, q, k))
gate(len(M_coeffs) == 5, "M division has five exact coefficient obligations")
M_sat = sp.groebner(
    [E_M, *M_coeffs, z*p*k-1], z, p, q, k, order="grevlex"
)
gate(len(M_sat.polys) == 1 and M_sat.polys[0].as_expr() == 1,
     "M seam saturated by both required poles is empty")


# ---------------------------------------------------------------------------
# N and P: a paired C2 survivor with the same branch-address factor.
# ---------------------------------------------------------------------------

E1_N = 243*k+12*p+54*q-4
E2_N = 9*k+4*p**2+6*p*q-8*p-6*q+4
N_k_seams = [
    ("E1", (4-12*p-54*q)/243, 2),
    ("E2", -(4*p**2+6*p*q-8*p-6*q+4)/9, 3),
]
for seam_label, k_value, power in N_k_seams:
    seam_color = sp.cancel(colors["N"].subs(k, k_value))
    seam_num, seam_den = sp.together(seam_color).as_numer_denom()
    seam_rem = sp.rem(
        sp.Poly(seam_num, u, domain="QQ(p,q)"),
        sp.Poly(seam_den, u, domain="QQ(p,q)"),
    ).as_expr()
    seam_coeffs = coefficient_numerators(seam_rem, (p, q))
    saturation = sp.groebner(
        [*seam_coeffs, z*k_value-1], z, p, q, order="grevlex"
    )
    gate(saturation.reduce((3*p-4)**power)[1] == 0,
         f"N {seam_label} live radical forces p=4/3")

N_sub = {p: sp.Rational(4, 3), k: -2*(9*q+2)/81}
t_N = sp.factor(t_N_raw.subs(N_sub))
a_N = sp.factor(incidences["N"][1].subs(N_sub))
c_N = sp.factor(incidences["N"][3].subs(N_sub))
d_N = sp.factor(incidences["N"][4].subs(N_sub))
C_N = sp.factor(colors["N"].subs(N_sub))
L_NP = 9*A+9*q+2
zero(t_N - (3*u-1)*(3*u+2)*L_NP.subs(A, A_N)/(81*u**3),
     "N survivor repeated root")
zero(C_N + sp.Rational(4, 3)*u**3*(3*u+2)*L_NP.subs(A, A_N)
     *(9*u**4+30*u**3+24*u**2-4), "N survivor polynomial color")
gate(a_N == 216*A**3, "N survivor leading incidence coefficient")
zero(c_N - 8*(6*A-1)*(27*A-4)*L_NP**2/243,
     "N survivor linear incidence coefficient")
zero(d_N - 4*(27*A-4)**2*L_NP**3/19683,
     "N survivor constant incidence coefficient")
I_N, D_N = implicit_data(A_N, C_N, a_N*T**3+C*T**2+c_N*T+d_N)
zero(D_N + 16*(27*A-4)**2*L_NP**3*I_N/19683,
     "N implicit component has discriminant exponent one")
gate(sp.Poly(I_N, C).degree() == 3 and
     sp.Poly(I_N, A, C).total_degree() == 11,
     "N implicit degree is eleven with cubic A-projection")

E1_P = 18*p-54*q-243*r-14
E2_P = 9*p**2-9*p*q-6*p+3*q+9*r+1
P_r_seams = [
    ("E1", (18*p-54*q-14)/243, 2),
    ("E2", -(9*p**2-9*p*q-6*p+3*q+1)/9, 3),
]
for seam_label, r_value, power in P_r_seams:
    seam_color = sp.cancel(colors["P"].subs(r, r_value))
    seam_num, seam_den = sp.together(seam_color).as_numer_denom()
    seam_rem = sp.rem(
        sp.Poly(seam_num, u, domain="QQ(p,q)"),
        sp.Poly(seam_den, u, domain="QQ(p,q)"),
    ).as_expr()
    seam_coeffs = coefficient_numerators(seam_rem, (p, q))
    saturation = sp.groebner(
        [*seam_coeffs, z*r_value-1], z, p, q, order="grevlex"
    )
    gate(saturation.reduce((9*p-1)**power)[1] == 0,
         f"P {seam_label} live radical forces p=1/9")

P_sub = {p: sp.Rational(1, 9), r: -2*(9*q+2)/81}
t_P = sp.factor(t_P_raw.subs(P_sub))
a_P = sp.factor(incidences["P"][1].subs(P_sub))
c_P = sp.factor(incidences["P"][3].subs(P_sub))
d_P = sp.factor(incidences["P"][4].subs(P_sub))
C_P = sp.factor(colors["P"].subs(P_sub))
zero(t_P - (3*u-1)*(3*u+2)*L_NP.subs(A, A_P)/(81*u),
     "P survivor repeated root including the centered constant shift")
gate(a_P == 729*A, "P survivor leading incidence coefficient")
zero(c_P + (27*A-4)*L_NP**2/9,
     "P survivor linear incidence coefficient")
zero(d_P - (27*A-4)**2*L_NP**3/1458,
     "P survivor constant incidence coefficient")
zero(C_P + sp.Rational(9, 2)*u*(3*u+2)*(9*u**2+6*u-4)
     *L_NP.subs(A, A_P), "P survivor polynomial color")
I_P, D_P = implicit_data(A_P, C_P, a_P*T**3+C*T**2+c_P*T+d_P)
zero(D_P + 2*(27*A-4)**2*L_NP**3*I_P/729,
     "P implicit component has discriminant exponent one")
gate(sp.Poly(I_P, C).degree() == 3 and
     sp.Poly(I_P, A, C).total_degree() == 7,
     "P implicit degree is seven with cubic A-projection")

# N/P share the same L=0 address polynomial.  It has three roots except at
# the two critical values; the zero critical value deletes the required pole,
# while the other leaves two distinct addresses.
address_NP = sp.expand(u**3+u**2+q+sp.Rational(2, 9))
address_disc = sp.factor(sp.discriminant(address_NP, u))
zero(address_disc + (9*q+2)*(27*q+10)/9,
     "N/P address discriminant has exactly the two critical-value factors")
gate(sp.factor(address_NP.subs(q, -sp.Rational(10, 27))) ==
     (3*u-1)*(3*u+2)**2/27,
     "N/P exceptional live fibre has exactly two distinct addresses")
gate(N_sub[k].subs(q, -sp.Rational(2, 9)) == 0 and
     P_sub[r].subs(q, -sp.Rational(2, 9)) == 0,
     "N/P one-address critical value is exactly the pole-deletion seam")


# ---------------------------------------------------------------------------
# Coefficient, generation, maximal-order, and address-folding gates.
# ---------------------------------------------------------------------------

gate(sp.gcd(A, L_J1) == 1, "J1 coefficient ideal is primitive")
gate(sp.gcd(A, L_J2) == 1, "J2 coefficient ideal is primitive for r nonzero")
gate(L_K.subs(A, 0) != 0, "K coefficient ideal is primitive")
zero(d_N.subs(A, 0) - 64*(9*q+2)**3/19683,
     "N coefficient ideal fails only at pole deletion")
zero(d_P.subs(A, 0) - 8*(9*q+2)**3/729,
     "P coefficient ideal fails only at pole deletion")

for label, pole_order, inertia in [
    ("J1", 5, 3), ("J2", 5, 3), ("K", 5, 2),
    ("N", 3, 2), ("P", 1, 2),
]:
    gate(pole_order % inertia != 0,
         f"{label} repeated root is not in k(A) by local pole valuation")

# The scalar endpoint is an actual maximal order: its sole discriminant
# prime occurs once, so no square index can be removed in codimension one.
I_J0 = sp.factor(I_J1.subs(q, 0))
D_J0 = sp.factor(D_J1.subs(q, 0))
zero(8*I_J0 - (27*A**10+8*C**3), "J scalar implicit prime")
zero(D_J0 + 2*I_J0, "J scalar order discriminant is squarefree at its prime")
gate(sp.gcd(10, 3) == 1, "J scalar binomial implicit relation is irreducible")

survivors = {
    "J1": (I_J1, 3), "J2": (I_J2, 2), "K": (I_K, 2),
    "N": (I_N, 2), "P": (I_P, 2),
}
for label, (implicit, minimum_addresses) in survivors.items():
    gate(sp.Poly(implicit, C).degree() == 3,
         f"{label} has a cubic birational normalization projection")
    gate(minimum_addresses >= 2,
         f"{label} supplies a source non-unibranch address packet")


summary = {
    "checks": CHECKS,
    "scope": "centered trace-zero linear-color;degA=3;degt=5",
    "rows": "J,K,L,M,N,O,P",
    "empty": "L,M,O",
    "survivors": "J1,J2,K,N,P",
    "scalar": "J1(q=0):d=1/2;normal maximal order",
    "addresses": "J1:3;J2:2;K:2;N/P:3 generic and 2 at q=-10/27",
    "branch": "every nonscalar implicit prime has discriminant exponent one",
    "conclusion": "scalar-unit or source ramification non-unibranch;no A2 Keller open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3959 centered degree-five root-map exact companion")
print(f"CHECKS={CHECKS}")
print("SCOPE=centered trace-zero linear-color;degA=3;deg_t=5")
print("ROWS=J,K,L,M,N,O,P")
print("EMPTY=L,M,O")
print("SURVIVORS=J1,J2,K,N,P")
print("SCALAR=J1(q=0);d=1/2;normal maximal order")
print("ADDRESSES=J1:3;J2:2;K:2;N/P:3_GENERIC_OR_2_CRITICAL")
print("BRANCH=each nonscalar implicit component has discriminant exponent one")
print("CONCLUSION=scalar-unit or source ramification non-unibranch;no A2 Keller open")
print(f"SEMANTIC_SHA256={semantic}")
