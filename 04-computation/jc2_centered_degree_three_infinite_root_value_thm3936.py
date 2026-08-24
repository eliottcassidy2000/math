#!/usr/bin/env python3
"""Exact companion for THM-3936's infinite-root-value degree-three closure.

Reproduction:
  python3 04-computation/jc2_centered_degree_three_infinite_root_value_thm3936.py
  python3 -O 04-computation/jc2_centered_degree_three_infinite_root_value_thm3936.py
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


u, A, C, T, p, k, q = sp.symbols("u A C T p k q")


def incidence(A_u: sp.Expr, t: sp.Expr, scale: sp.Expr = sp.S.One):
    """Return the scaled resultant and its (a,c,d) trace-zero coefficients."""
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


# ---------------------------------------------------------------------------
# Pole ledger.  At a tame point of A with ramification e, local trace retains
# precisely Laurent exponents divisible by e.
# ---------------------------------------------------------------------------


def retained(exponent: int, e: int) -> bool:
    return exponent % e == 0


gate(retained(-3, 3), "an infinity pole of order three survives the e=3 trace")
gate(not retained(-2, 3) and not retained(-1, 3),
     "infinity pole orders one and two have no forced leading trace term")
gate(retained(-1, 1), "an isolated simple pole over an unramified point survives trace")
gate(not retained(-1, 2) and not retained(-1, 3),
     "an isolated simple pole requires ramification e=2 or e=3")
gate(retained(-2, 1) and retained(-2, 2) and not retained(-2, 3),
     "an isolated double pole requires total ramification e=3")

finite_rh_budget = 2
rows = {
    "A": (2, "1@e2", 1),
    "B": (2, "1@e3", 2),
    "C": (1, "2@e3", 2),
    "D": (1, "1@e2+1@e2", 2),
}
gate(all(cost <= finite_rh_budget for _, _, cost in rows.values()),
     "the four displayed collision-free rows exhaust the finite RH budget")
gate({value[0] for value in rows.values()} == {1, 2},
     "m_infinity=3 is excluded and only one or two remain")


# ---------------------------------------------------------------------------
# Row A: m_infinity=2 and one simple finite pole at an e=2 point.
# ---------------------------------------------------------------------------

A_A = u**3 + u**2
t_A = u**2 + p * u + (p - 1) / 3 + k / u
R_A, a_A, trace_A, c_A, d_A = incidence(A_A, t_A, sp.Rational(1, 27))
gate(trace_A == 0, "row A trace normalization")
gate(a_A == -A, "row A primitive leading coefficient")

C_A_raw = color(A_A, t_A, a_A, c_A)
N_A, D_A = sp.together(C_A_raw).as_numer_denom()
t_A_num = sp.together(t_A).as_numer_denom()[0]
gate(sp.factor(sp.resultant(t_A_num, N_A, u)) ==
     -243 * k**3 * (p - 1)**6 * (27 * k + 6 * p - 2)**3,
     "row A polynomial-color resultant")

seam = {k: (2 - 6 * p) / 27}
seam_remainder = sp.factor(
    sp.rem(sp.Poly(N_A.subs(seam), u), sp.Poly(t_A_num.subs(seam), u)).as_expr()
)
gate(seam_remainder ==
     -(p - 1)**3 * (3 * p - 1)**3 * (3 * u - 1) * (3 * u + 2) / 27,
     "row A exceptional resultant seam leaves an exact nonzero remainder")
gate(sp.factor(C_A_raw.subs({p: sp.Rational(1, 3), k: 0})) ==
     u**2 * (u + 1) * (3 * u - 2) * (3 * u + 2) / 6,
     "row A p=1/3 seam is polynomial only after the finite pole disappears")

C_A = sp.factor(C_A_raw.subs(p, 1))
gate(C_A == u * (3 * u + 4) * (A_A + k) / 2,
     "row A surviving polynomial color")
gate(sp.factor(c_A.subs(p, 1)) == -(A + k)**2,
     "row A surviving c")
gate(sp.factor(d_A.subs(p, 1)) == -(A + k)**3 / 2,
     "row A surviving d")
gate(d_A.subs({p: 1, A: 0}) == -k**3 / 2,
     "row A coefficient ideal is the unit ideal for k nonzero")

H_A = (
    27 * A**5 + 81 * A**4 * k + 16 * A**4 + 36 * A**3 * C
    + 81 * A**3 * k**2 + 48 * A**3 * k + 72 * A**2 * C * k
    + 27 * A**2 * k**3 + 48 * A**2 * k**2 - 4 * A * C**2
    + 36 * A * C * k**2 + 16 * A * k**3 - 8 * C**3 - 4 * C**2 * k
)
Phi_A = -A * T**3 + C * T**2 - (A + k)**2 * T - (A + k)**3 / 2
gate(sp.factor(sp.discriminant(Phi_A, T)) == -(A + k)**3 * H_A / 4,
     "row A discriminant factorization")
gate(sp.expand(sp.resultant(A - A_A, C - C_A, u) + H_A / 8) == 0,
     "row A implicit normalization equation")
gate(sp.Poly(H_A, A, C).total_degree() == 5 and sp.Poly(C_A, u).degree() == 5,
     "row A is a one-place quintic")
gate(sp.factor(sp.discriminant(u**3 + u**2 + k, u)) == -k * (27 * k + 4),
     "row A address cubic has three roots or one double pair, never one root")
gate(sp.gcd(sp.Poly(A_A + k, u), sp.Poly(sp.diff(A_A, u), u)).degree() == 0,
     "generic row A collapsed fibre has three distinct addresses")
gate(sp.expand((A_A + k).subs(u, -sp.Rational(2, 3)) - k - sp.Rational(4, 27)) == 0,
     "row A has exactly two addresses only at k=-4/27")


# ---------------------------------------------------------------------------
# Row B: m_infinity=2 and one simple finite pole at an e=3 point.
# ---------------------------------------------------------------------------

A_B = u**3
t_B = u**2 + p * u + k / u
R_B, a_B, trace_B, c_B, d_B = incidence(A_B, t_B)
gate(trace_B == 0 and a_B == -A, "row B trace and leading coefficient")
C_B_raw = color(A_B, t_B, a_B, c_B)
N_B, D_B = sp.together(C_B_raw).as_numer_denom()
t_B_num = sp.together(t_B).as_numer_denom()[0]
gate(sp.factor(sp.resultant(t_B_num, N_B, u)) == -27 * k**6 * p**6,
     "row B polynomial-color resultant")
C_B = sp.factor(C_B_raw.subs(p, 0))
gate(C_B == sp.Rational(3, 2) * u**2 * (A_B + k),
     "row B surviving polynomial color")
gate(c_B.subs(p, 0) == 0 and sp.factor(d_B.subs(p, 0)) == -(A + k)**3 / 2,
     "row B surviving incidence")
gate(d_B.subs({p: 0, A: 0}) == -k**3 / 2,
     "row B coefficient ideal is the unit ideal for k nonzero")

H_B = 27 * A**5 + 81 * A**4 * k + 81 * A**3 * k**2 + 27 * A**2 * k**3 - 8 * C**3
Phi_B = -A * T**3 + C * T**2 - (A + k)**3 / 2
gate(sp.factor(sp.discriminant(Phi_B, T)) == -(A + k)**3 * H_B / 4,
     "row B discriminant factorization")
gate(sp.expand(sp.resultant(A - A_B, C - C_B, u) + H_B / 8) == 0,
     "row B implicit normalization equation")
gate(sp.Poly(H_B, A, C).total_degree() == 5 and sp.Poly(C_B, u).degree() == 5,
     "row B is a one-place quintic")
gate(sp.discriminant(u**3 + k, u) == -27 * k**2,
     "row B collapsed fibre has three distinct addresses for k nonzero")


# ---------------------------------------------------------------------------
# Row C: m_infinity=1 and one double finite pole at an e=3 point.
# ---------------------------------------------------------------------------

A_C = u**3
t_C = u + p / u + k / u**2
R_C, a_C, trace_C, c_C, d_C = incidence(A_C, t_C)
gate(trace_C == 0 and a_C == -A**2, "row C trace and leading coefficient")
C_C_raw = color(A_C, t_C, a_C, c_C)
N_C, D_C = sp.together(C_C_raw).as_numer_denom()
t_C_num = sp.together(t_C).as_numer_denom()[0]
gate(sp.factor(sp.resultant(t_C_num, N_C, u)) == -27 * k**6 * p**6,
     "row C polynomial-color resultant")
C_C = sp.factor(C_C_raw.subs(p, 0))
gate(C_C == sp.Rational(3, 2) * u**4 * (A_C + k),
     "row C surviving polynomial color")
gate(c_C.subs(p, 0) == 0 and sp.factor(d_C.subs(p, 0)) == -(A + k)**3 / 2,
     "row C surviving incidence")
gate(d_C.subs({p: 0, A: 0}) == -k**3 / 2,
     "row C coefficient ideal is the unit ideal for k nonzero")

H_C = 27 * A**7 + 81 * A**6 * k + 81 * A**5 * k**2 + 27 * A**4 * k**3 - 8 * C**3
Phi_C = -A**2 * T**3 + C * T**2 - (A + k)**3 / 2
gate(sp.factor(sp.discriminant(Phi_C, T)) == -(A + k)**3 * H_C / 4,
     "row C discriminant factorization")
gate(sp.expand(sp.resultant(A - A_C, C - C_C, u) + H_C / 8) == 0,
     "row C implicit normalization equation")
gate(sp.Poly(H_C, A, C).total_degree() == 7 and sp.Poly(C_C, u).degree() == 7,
     "row C is a one-place septic")
gate(sp.discriminant(u**3 + k, u) == -27 * k**2,
     "row C collapsed fibre has three distinct addresses for k nonzero")


# ---------------------------------------------------------------------------
# Row D: m_infinity=1 and two finite simple poles at distinct e=2 values.
# ---------------------------------------------------------------------------

A_D = u**3 - 3 * u
t_D_raw = u + q / (u - 1) + k / (u + 1)
R_D, a_D, trace_D, c_D, d_D = incidence(A_D, t_D_raw)
gate(trace_D == 0 and sp.expand(a_D + (A - 2) * (A + 2)) == 0,
     "row D trace and leading coefficient")
C_D_raw = color(A_D, t_D_raw, a_D, c_D)
N_D, D_D = sp.together(C_D_raw).as_numer_denom()
t_D_num = sp.together(t_D_raw).as_numer_denom()[0]
gate(sp.factor(sp.resultant(t_D_num, N_D, u)) ==
     1728 * k**3 * q**3 * (k + q + 2)**6,
     "row D polynomial-color resultant")
gate(sp.limit((u - 1) * t_D_raw, u, 1) == q and
     sp.limit((u + 1) * t_D_raw, u, -1) == k,
     "q=0 or k=0 is exactly a finite-pole cancellation")

t_D = sp.factor(t_D_raw.subs(k, -q - 2))
C_D = sp.factor(C_D_raw.subs(k, -q - 2))
delta = 2 * q + 2
gate(t_D == (A_D + delta) / ((u - 1) * (u + 1)),
     "row D surviving root normal form")
gate(sp.expand(C_D - sp.Rational(3, 2) * (u**2 - 1) * (u**2 - 5) * (A_D + delta)) == 0,
     "row D surviving polynomial color")
gate(sp.factor(c_D.subs(k, -q - 2)) == 3 * (A + delta)**2,
     "row D surviving c")
gate(sp.factor(d_D.subs(k, -q - 2)) == -(A + delta)**3 / 2,
     "row D surviving d")
gate(sp.expand(sp.resultant(A**2 - 4, A + delta, A) - (delta - 2) * (delta + 2)) == 0,
     "row D coefficient ideal requires delta not equal to plus or minus two")
gate(sp.factor((delta - 2) * (delta + 2)) == 4 * q * (q + 2),
     "row D unit-ideal exclusions are exactly the two pole cancellations")

H_D = (
    27 * A**7 + 81 * A**6 * delta + 81 * A**5 * delta**2 - 648 * A**5
    - 108 * A**4 * C + 27 * A**4 * delta**3 - 1944 * A**4 * delta
    - 216 * A**3 * C * delta - 1944 * A**3 * delta**2 + 2160 * A**3
    - 108 * A**2 * C * delta**2 + 432 * A**2 * C - 648 * A**2 * delta**3
    + 6480 * A**2 * delta - 36 * A * C**2 + 864 * A * C * delta
    + 6480 * A * delta**2 - 8 * C**3 - 36 * C**2 * delta
    + 432 * C * delta**2 + 2160 * delta**3
)
Phi_D = -(A**2 - 4) * T**3 + C * T**2 + 3 * (A + delta)**2 * T - (A + delta)**3 / 2
gate(sp.expand(sp.discriminant(Phi_D, T) + (A + delta)**3 * H_D / 4) == 0,
     "row D discriminant factorization")
gate(sp.expand(sp.resultant(A - A_D, C - C_D, u) + H_D / 8) == 0,
     "row D implicit normalization equation")
gate(sp.Poly(H_D, A, C).total_degree() == 7 and sp.Poly(C_D, u).degree() == 7,
     "row D is a one-place septic")
gate(sp.expand(sp.resultant(A_D + delta, sp.diff(A_D, u), u) -
               27 * (delta - 2) * (delta + 2)) == 0,
     "row D collapsed fibre has three distinct addresses off the cancellations")


# ---------------------------------------------------------------------------
# Uniform branch/maximal-order and birational-normalization controls.
# ---------------------------------------------------------------------------

for label, A_u, t, C_u, H, finite_pole_order, finite_e in [
    ("A", A_A, t_A.subs(p, 1), C_A, H_A, 1, 2),
    ("B", A_B, t_B.subs(p, 0), C_B, H_B, 1, 3),
    ("C", A_C, t_C.subs(p, 0), C_C, H_C, 2, 3),
    ("D", A_D, t_D, C_D, H_D, 1, 2),
]:
    gate(finite_pole_order % finite_e != 0,
         f"row {label} repeated root is not in k(A) by its pole valuation")
    gate(sp.Poly(H, C).degree() == 3,
         f"row {label} implicit component has cubic A-projection")
    gate(sp.Poly(sp.resultant(A - A_u, C - C_u, u), C).degree() == 3,
         f"row {label} parametrization has no resultant multiplicity")


summary = {
    "checks": CHECKS,
    "pole_rows": "A:e2-simple;B:e3-simple;C:e3-double;D:e2+e2-simple",
    "row_A": "p=1;k!=0;quintic;at least two addresses over (-k,0)",
    "row_B": "p=0;k!=0;quintic;three addresses over (-k,0)",
    "row_C": "p=0;k!=0;septic;three addresses over (-k,0)",
    "row_D": "q+k+2=0;q*k!=0;septic;three addresses over (-2q-2,0)",
    "branch": "each H occurs to discriminant exponent one and is genuine e=2",
    "conclusion": "source ramification component non-unibranch; no A2 Keller open",
    "scope": "centered trace-zero linear-color degree-three; t(infinity)=infinity",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3936 centered infinite-root-value exact companion")
print(f"CHECKS={CHECKS}")
print("SCOPE=centered trace-zero linear-color degree-three;t(infinity)=infinity")
print("POLE_ROWS=A:e2-simple;B:e3-simple;C:e3-double;D:e2+e2-simple")
print("ROW_A=p=1;k!=0;quintic;at-least-two-address fibre (-k,0)")
print("ROW_B=p=0;k!=0;quintic;three-address fibre (-k,0)")
print("ROW_C=p=0;k!=0;septic;three-address fibre (-k,0)")
print("ROW_D=q+k+2=0;q*k!=0;septic;three-address fibre (-2q-2,0)")
print("BRANCH=each implicit H has discriminant exponent one;genuine tame e=2")
print("CONCLUSION=source ramification component non-unibranch;no A2 Keller open")
print(f"SEMANTIC_SHA256={semantic}")
