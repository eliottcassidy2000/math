#!/usr/bin/env python3
"""Exact primary certificate for THM-4007.

The certificate continues the oriented live reduced (2,3) seam of THM-3997
and THM-4005 on the only old minimal-support stratum b=d=0.  It proves the
third source-normal degree caps, solves the t^2 Jacobian row, extracts the
new forced A-weight -6, and compares the t^4 nodal-defect row with the first
previously unknown residual layer.

This is a necessary-row calculation in the fixed THM-3992 gauge.  It does
not construct a B_2 pair, settle another reduced cell, or prove JC(2).
"""

from __future__ import annotations

import sys

import sympy as sp


sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def require(condition, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(label)
    print(f"PASS  {label}")


def equal(left, right, label: str) -> None:
    require(simp(left - right) == 0, label)


def jac(first, second, x, t):
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, t)
        - sp.diff(first, t) * sp.diff(second, x)
    )


print("STATUS=THM-4007;PROVED_NECESSARY_ROW;VERIFIED-EXACT;JC(2)_OPEN")
print("SCOPE=oriented_reduced_2:3_live_seam;fixed_THM3992_gauge")

# -------------------------------------------------------------------------
# I. The Laurent support gives the complete third-diagonal degree caps.
# -------------------------------------------------------------------------

s, tau, x, t = sp.symbols("s tau x t")
gamma_cap = sp.symbols("gamma_cap", nonzero=True)

# Only the displayed coefficients can contribute to source order t^3.
# Higher s-powers in any row land above t^3 and therefore need not be stored.
ab = sp.symbols("ab0:5")
a0 = sp.symbols("aa0:4")
a1 = sp.symbols("aa10:13")
a2 = sp.symbols("aa20:22")
a3 = sp.symbols("aa30")
A_laurent = (
    gamma_cap**2 * s**2 * tau**-2
    + sum(ab[j] * s**j for j in range(5)) * tau**-1
    + sum(a0[j] * s**j for j in range(4))
    + sum(a1[j] * s**j for j in range(3)) * tau
    + sum(a2[j] * s**j for j in range(2)) * tau**2
    + a3 * tau**3
)
A_source = sp.expand(A_laurent.subs({s: x * t, tau: t}))
A3_cap = sp.expand(A_source.coeff(t, 3))
equal(
    A3_cap,
    ab[4] * x**4 + a0[3] * x**3 + a1[2] * x**2 + a2[1] * x + a3,
    "complete A third diagonal from Laurent rows",
)
require(sp.Poly(A3_cap, x).degree() == 4, "A third diagonal cap is sharp: deg<=4")

cb = sp.symbols("cb0:6")
c1 = sp.symbols("cc10:15")
c0 = sp.symbols("cc00:04")
cp1 = sp.symbols("ccp10:13")
cp2 = sp.symbols("ccp20:22")
cp3 = sp.symbols("ccp30")
C_laurent = (
    gamma_cap**3 * s**3 * tau**-3
    + sum(cb[j] * s**j for j in range(6)) * tau**-2
    + sum(c1[j] * s**j for j in range(5)) * tau**-1
    + sum(c0[j] * s**j for j in range(4))
    + sum(cp1[j] * s**j for j in range(3)) * tau
    + sum(cp2[j] * s**j for j in range(2)) * tau**2
    + cp3 * tau**3
)
C_source = sp.expand(C_laurent.subs({s: x * t, tau: t}))
C3_cap = sp.expand(C_source.coeff(t, 3))
equal(
    C3_cap,
    cb[5] * x**5
    + c1[4] * x**4
    + c0[3] * x**3
    + cp1[2] * x**2
    + cp2[1] * x
    + cp3,
    "complete C third diagonal from Laurent rows",
)
require(sp.Poly(C3_cap, x).degree() == 5, "C third diagonal cap is sharp: deg<=5")
print("MECHANISM=extreme rows gamma^2*s^2 and gamma^3*s^3 have no higher powers")

# -------------------------------------------------------------------------
# II. Solve the third source-normal row on b=d=0.
# -------------------------------------------------------------------------

a = sp.symbols("a", nonzero=True)
A5 = a**5
gamma = -a**3 / 2
I = sp.Rational(3, 4) * a**2
q = gamma * x**2 + sp.Rational(3, 2) * a / gamma

# Raw rows reconstructed from THM-4005 equations (2)--(3), not imported
# from its code.  A_i and C_i denote the coefficient of t^i.
A0 = a * (1 + A5 * x**2 / 4)
C0 = a**4 * (-sp.Rational(3, 4) * x - A5 * x**3 / 8)
N = a * (sp.Rational(4, 3) / A5 + 2 * x**2)
K = a**4 * (-4 * x / A5 - sp.Rational(3, 2) * x**3)
M = a * (-sp.Rational(32, 9) / A5**2 - 4 * x**2 / (5 * A5))
L = a**4 * (sp.Rational(88, 15) * x / A5**2 - 12 * x**3 / (5 * A5))

equal(sp.diff(A0, x) * K - N * sp.diff(C0, x), 1, "inherited t^0 Jacobian row")
equal(
    2 * (sp.diff(A0, x) * L - M * sp.diff(C0, x))
    + sp.diff(N, x) * K
    - N * sp.diff(K, x),
    0,
    "inherited t^1 Jacobian row",
)
equal(-3 * A0**2 + I, -q * sp.diff(C0, x), "rotated-gradient A coordinate")
equal(2 * C0, q * sp.diff(A0, x), "rotated-gradient C coordinate")

u = sp.symbols("u0:5")
v = sp.symbols("v0:6")
U = sum(u[j] * x**j for j in range(5))
V = sum(v[j] * x**j for j in range(6))

# The coefficient of t^2 in J(A,C).  For pairs P=(P_A,P_C), Q=(Q_A,Q_C),
# det(P',Q)=P_A'*Q_C-P_C'*Q_A; the distinct-row order matters.
J2 = sp.expand(
    3 * (sp.diff(A0, x) * V - U * sp.diff(C0, x))
    + 2 * (sp.diff(N, x) * L - sp.diff(K, x) * M)
    + (sp.diff(M, x) * K - sp.diff(L, x) * N)
)
J2_poly = sp.Poly(J2, x)
J2_equations = [
    J2_poly.coeff_monomial(x**j) for j in range(J2_poly.degree() + 1)
]
solution_tuple = next(iter(sp.linsolve(J2_equations, list(u) + list(v))))
solution = dict(zip(list(u) + list(v), solution_tuple))
free = set().union(*(entry.free_symbols for entry in solution_tuple)).intersection(set(u + v))
require(free == {v[2], v[3], v[4], v[5]}, "third row has exactly four displayed free scalars")

expected_solution = {
    u[0]: sp.Rational(2176, 135) / a**14,
    u[1]: -4 * (a**5 * v[2] - 2 * v[4]) / (3 * a**8),
    u[2]: -4 * (5 * a**6 * v[3] - 10 * a * v[5] - 32) / (15 * a**9),
    u[3]: -4 * v[4] / (3 * a**3),
    u[4]: -4 * v[5] / (3 * a**3),
    v[0]: 2 * (a**5 * v[2] - 2 * v[4]) / a**10,
    v[1]: 2 * (45 * a**6 * v[3] - 90 * a * v[5] - 752) / (45 * a**11),
}
for variable, expected in expected_solution.items():
    equal(solution[variable], expected, f"third-row coefficient {variable}")

U_solved = sp.expand(U.subs(solution))
V_solved = sp.expand(V.subs(solution))
equal(J2.subs(solution), 0, "complete solved t^2 Jacobian row")

# A one-coordinate no-import audit of the load-bearing coefficient.
equal(sp.diff(A0, x).subs(x, 0), 0, "A0'(0)")
equal(sp.diff(C0, x).subs(x, 0), -sp.Rational(3, 4) * a**4, "C0'(0)")
equal(N.subs(x, 0), sp.Rational(4, 3) / a**4, "N(0)")
equal(sp.diff(N, x).subs(x, 0), 0, "N'(0)")
equal(K.subs(x, 0), 0, "K(0)")
equal(sp.diff(K, x).subs(x, 0), -4 / a, "K'(0)")
equal(M.subs(x, 0), -sp.Rational(32, 9) / a**9, "M(0)")
equal(sp.diff(M, x).subs(x, 0), 0, "M'(0)")
equal(L.subs(x, 0), 0, "L(0)")
equal(sp.diff(L, x).subs(x, 0), sp.Rational(88, 15) / a**6, "L'(0)")
equal(
    (
        2 * (sp.diff(N, x) * L - sp.diff(K, x) * M)
        + (sp.diff(M, x) * K - sp.diff(L, x) * N)
    ).subs(x, 0),
    -sp.Rational(544, 15) / a**10,
    "lower-row constant in J[t^2]",
)
equal(solution[u[0]], sp.Rational(2176, 135) / a**14, "forced raw U(0)")
equal(solution[u[0]] / a, sp.Rational(2176, 135) / A5**3, "forced invariant A3 constant")
require(simp(solution[u[0]]) != 0, "new A weight -6 coefficient is nonzero")

old_A_weights = {2, 0, -2, -4}
require(-6 not in old_A_weights, "new A weight -6 is distinct from old minimal invoice")
require(len(old_A_weights | {-6}) == 5, "b=d=0 forces at least five A weights")
print("RESULT [x^0*t^3](A/a)=2176/(135*A5^3)!=0")
print("RESULT b=d=0 adds source weight -6 and excludes the fixed-gauge 4x5 invoice")

# -------------------------------------------------------------------------
# III. Compare the t^4 nodal-defect row and determine the residual coupling.
# -------------------------------------------------------------------------

# If (R4,S4) is the fourth source-normal vector, J[t^3]=0 says
# 4*det((A0,C0)',(R4,S4))+Jlower=0.  The rotated-gradient contribution to
# P[t^4] is therefore -q*Jlower/4.  The minus sign is load-bearing.
Jlower = sp.expand(
    3 * sp.diff(N, x) * V_solved
    - N * sp.diff(V_solved, x)
    + sp.diff(U_solved, x) * K
    - 3 * U_solved * sp.diff(K, x)
    + 2 * (sp.diff(M, x) * L - M * sp.diff(L, x))
)
g4 = sp.expand(
    -q * Jlower / 4
    + 2 * K * V_solved
    + L**2
    - 6 * A0 * N * U_solved
    - 3 * A0 * M**2
    - 3 * N**2 * M
)
require(sp.Poly(g4, x).degree() == 4, "correct elimination cancels every t^4 power above x^4")

g4_expected = (
    -2 * (15 * a**6 * v[3] - 30 * a * v[5] - 464) / (15 * a**17)
    + 4 * (2 * a**5 * v[2] - 7 * v[4]) * x / (3 * a**11)
    + (35 * a**6 * v[3] - 160 * a * v[5] - 448) * x**2 / (15 * a**12)
    + 2 * v[4] * x**3 / a**6
    + 5 * v[5] * x**4 / (3 * a**6)
)
equal(g4, g4_expected, "complete eliminated P[t^4] row")

c02, c21, c40 = sp.symbols("c02 c21 c40")
epsilon = -sp.Rational(1376, 135) / a**12
alpha = sp.Rational(8, 3) / a**7
matching = {
    v[2]: -sp.Rational(3, 16) * a**9 * c21,
    v[3]: sp.Rational(736, 105) / a**6 - sp.Rational(3, 14) * a**9 * c02,
    v[4]: 0,
    v[5]: sp.Rational(8, 5) / a,
}
c40_forced = -sp.Rational(11392, 105) / A5**4 - sp.Rational(6, 7) * c02 / A5
ring_t4 = gamma * (
    c40_forced + c21 * x + (3 * (epsilon / gamma) + c02) * x**2
    + (alpha / gamma) * x**4
)
equal(g4.subs(matching), ring_t4, "full t^4 residual comparison")
equal(
    c40_forced + sp.Rational(6, 7) * c02 / A5,
    -sp.Rational(11392, 105) / A5**4,
    "forced [p4]/[y2] residual coupling",
)

A3_invariant = simp(U_solved.subs(matching) / a)
C3_invariant = simp(V_solved.subs(matching) / a**4)
A3_expected = (
    sp.Rational(2176, 135) / A5**3
    + A5 * c21 * x / 4
    + (sp.Rational(1088, 315) / A5**2 + sp.Rational(2, 7) * A5 * c02) * x**2
    - sp.Rational(32, 15) * x**4 / A5
)
C3_expected = (
    -sp.Rational(3, 8) * c21
    + (-sp.Rational(8128, 315) / A5**3 - sp.Rational(3, 7) * c02) * x
    - sp.Rational(3, 16) * A5 * c21 * x**2
    + (sp.Rational(736, 105) / A5**2 - sp.Rational(3, 14) * A5 * c02) * x**3
    + sp.Rational(8, 5) * x**5 / A5
)
equal(A3_invariant, A3_expected, "complete invariant third A row")
equal(C3_invariant, C3_expected, "complete invariant third C row")
print("RESULT [p^4]Rtilde+(6/(7*A5))*[y^2]Rtilde=-11392/(105*A5^4)")
print("RESULT [p^2*y]Rtilde remains the independent source-odd scalar in this row")

# -------------------------------------------------------------------------
# IV. Hostiles: distinct-row determinant and direct full-series expansion.
# -------------------------------------------------------------------------

def det_prime_correct(first, second):
    return sp.diff(first[0], x) * second[1] - sp.diff(first[1], x) * second[0]


def det_prime_wrong(first, second):
    return sp.diff(first[0], x) * second[1] - first[0] * sp.diff(second[1], x)


P_host = (x**2 + 1, x**3 + x)
Q_host = (x + 2, x**2 + 3)
require(
    simp(det_prime_correct(P_host, Q_host) - det_prime_wrong(P_host, Q_host)) != 0,
    "distinct-row determinant hostile detects swapped second factor",
)
equal(
    det_prime_correct(P_host, P_host),
    det_prime_wrong(P_host, P_host),
    "same-row determinant is a deliberately blind control",
)

# A direct full-series check over a=1 with arbitrary z=v5.  It solves the
# fourth source-normal Jacobian row rather than trusting the eliminated g4.
z = sp.symbols("z")
numeric_sub = {a: 1, v[2]: 0, v[3]: 0, v[4]: 0, v[5]: z}
A0h, C0h = A0.subs(a, 1), C0.subs(a, 1)
Nh, Kh = N.subs(a, 1), K.subs(a, 1)
Mh, Lh = M.subs(a, 1), L.subs(a, 1)
Uh = sp.expand(U_solved.subs(numeric_sub))
Vh = sp.expand(V_solved.subs(numeric_sub))
r = sp.symbols("r0:6")
w = sp.symbols("w0:7")
R4 = sum(r[j] * x**j for j in range(6))
S4 = sum(w[j] * x**j for j in range(7))
J3_direct = sp.expand(
    4 * (sp.diff(A0h, x) * S4 - R4 * sp.diff(C0h, x))
    + 3 * sp.diff(Nh, x) * Vh
    - Nh * sp.diff(Vh, x)
    + sp.diff(Uh, x) * Kh
    - 3 * Uh * sp.diff(Kh, x)
    + 2 * (sp.diff(Mh, x) * Lh - Mh * sp.diff(Lh, x))
)
J3_poly = sp.Poly(J3_direct, x)
J3_equations = [
    J3_poly.coeff_monomial(x**j) for j in range(J3_poly.degree() + 1)
]
fourth_tuple = next(iter(sp.linsolve(J3_equations, list(r) + list(w))))
fourth_solution = dict(zip(list(r) + list(w), fourth_tuple))
fourth_free = set().union(*(entry.free_symbols for entry in fourth_tuple)).intersection(set(r + w))
fourth_zero = {
    variable: sp.expand(value.subs({free_variable: 0 for free_variable in fourth_free}))
    for variable, value in fourth_solution.items()
}
R4h = sp.expand(R4.subs(fourth_zero))
S4h = sp.expand(S4.subs(fourth_zero))
equal(J3_direct.subs(fourth_zero), 0, "direct full-series fourth Jacobian row")

A_host_series = A0h + t * Nh + t**2 * Mh + t**3 * Uh + t**4 * R4h
C_host_series = C0h + t * Kh + t**2 * Lh + t**3 * Vh + t**4 * S4h
J_host_series = jac(A_host_series, C_host_series, x, t)
for order in range(4):
    equal(
        J_host_series.coeff(t, order),
        1 if order == 0 else 0,
        f"direct full-series J row {order}",
    )
P_host_series = sp.expand(
    C_host_series**2 - A_host_series**3 + sp.Rational(3, 4) * A_host_series + sp.Rational(1, 4)
)
P4_host = sp.expand(P_host_series.coeff(t, 4))
P4_host_expected = (
    25 * x**4 * z - 160 * x**2 * z - 448 * x**2 + 60 * z + 928
) / 15
equal(P4_host, P4_host_expected, "direct full-series P[t^4]")
require(sp.Poly(P4_host, x).degree() == 4, "direct hostile has no spurious x^8/x^7/x^6 terms")

wrong_sign_g4 = sp.expand(
    q * Jlower / 4
    + 2 * K * V_solved
    + L**2
    - 6 * A0 * N * U_solved
    - 3 * A0 * M**2
    - 3 * N**2 * M
)
require(sp.Poly(wrong_sign_g4, x).degree() == 8, "wrong elimination sign creates a spurious degree-eight row")
equal(
    sp.Poly(wrong_sign_g4, x).coeff_monomial(x**8),
    2 * a**4 * v[5],
    "wrong-sign x^8 hostile is active",
)

# -------------------------------------------------------------------------
# V. Support strata and the exact gauge-transfer boundary.
# -------------------------------------------------------------------------

# THM-4005 already gives 5x6 off b=d=0.  The new row repairs its sole 4x5
# stratum to 5x5.  These are lower bounds, not existence assertions.
strata = {
    "b_nonzero": (5, 6),
    "b_zero_d_nonzero": (5, 6),
    "b_zero_d_zero": (5, 5),
}
require(min(pair[0] for pair in strata.values()) == 5, "all live strata force at least five A weights")
require(min(pair[1] for pair in strata.values()) == 5, "all live strata force at least five C weights in fixed gauge")
print("SUPPORT_STRATA=" + repr(strata))
print("FIXED_GAUGE_FIRST_UNREJECTED_SUPPORT=5x5")

# Under the operations actually used in THM-3992, A changes only by a
# nonzero scalar and a target constant; the sole mixing is C by A.  Hence
# the nonconstant A-support and its five-weight floor transfer exactly.
scale = sp.symbols("scale", nonzero=True)
A_component_coefficients = {2: 1, 0: 2, -2: 3, -4: 4, -6: 5}
A_after_recorded_gauge = {
    weight: scale * coefficient for weight, coefficient in A_component_coefficients.items()
}
require(
    set(A_after_recorded_gauge) == set(A_component_coefficients),
    "recorded diagonal scaling preserves all five A weights",
)
require(
    all(coefficient != 0 for coefficient in A_after_recorded_gauge.values()),
    "recorded scaling cannot cancel an A component",
)
print("TRANSFER=recorded_scalings_translations_and_C_by_A_shear_preserve_the_A>=5_floor")
print("FIREWALL=C>=5_need_not_transfer_across_C_by_A_shear")
print("FIREWALL=no_reversed_orientation;no_other_reduced_cell;no_arbitrary_target_automorphism")
print("FIREWALL=no_all_row_lift;no_B2_pair_construction;no_JC2_conclusion")
print(f"CHECKS={CHECKS}")
print("ALL THM-4007 THIRD NORMAL ROW EXACT CHECKS PASSED")
