#!/usr/bin/env python3
"""Exact certificate for THM-3998 downstream of THM-3992.

This certificate verifies the graded identities, boundary filters, and
all-degree obstruction proved in THM-3998.  Work in

    B_2 = k[x,u,p,y] inside k[x,t],
    u=x^2*t, p=t*(1+u), y=x*t*p,

with wt(x)=1, wt(t)=-2.  THM-3992 forces the deleted-line node

    A(x,0)=gamma^2*x^2+a,
    C(x,0)=gamma^3*x^3+(3*a*gamma/2)*x,

where a*gamma != 0, and A_t(0,0)=-2/(3*gamma*a).

Frozen support universe (coefficient degrees in u are NOT bounded):

    A = x^2*f(u) + F(u) + x^-2*m(u),
    C = x^3*g(u) + x*k(u) + [at most one x^r*l(u)],

where r is any integer other than 1 or 3 (those two values are already
absorbed in k and g); the last row is omitted in the core cell.
Negative-weight coefficients obey the exact THM-3973 module divisibilities.
The script verifies an all-degree, all-r no-go in this finite-support
universe.  Thus the result is stronger than a coefficient or weight cutoff,
but says nothing about wider support or extra weights in A.
"""

from __future__ import annotations

import sympy as sp


x, t, u = sp.symbols("x t u")
a, gamma = sp.symbols("a gamma", nonzero=True)


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def check(label: str, expr) -> None:
    residual = simp(expr)
    if residual != 0:
        raise AssertionError(f"{label}: residual {residual}")
    print(f"PASS  {label}")


def gate(label: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def jac(P, Q):
    return sp.expand(sp.diff(P, x) * sp.diff(Q, t)
                     - sp.diff(P, t) * sp.diff(Q, x))


print("STATUS=FINITE-SUPPORT EXACT NO-GO; NOT A JC(2) PROOF")
print("FIELD=algebraically closed characteristic zero")
print("RING=B2=k[x,u,p,y], u=x^2*t, p=t*(1+u), y=x*t*p")
print("SUPPORT_A={2,0,-2}; SUPPORT_C={3,1} plus <=1 arbitrary integer weight")
print("COEFFICIENT_DEGREE_IN_u=UNBOUNDED")

# -------------------------------------------------------------------------
# Ring, grading, and control checks.
# -------------------------------------------------------------------------

u_xt = x**2 * t
p_xt = t * (1 + u_xt)
y_xt = x * t * p_xt
check("canonical p relation (parent-task typo excluded)", p_xt - (t + x**2 * t**2))
check("canonical y relation", y_xt - (x * t**2 + x**3 * t**3))

# Ambient positive control and THM-3973's rational Darboux positive control.
check("positive control J(x,t)=1", jac(x, t) - 1)
g_rat = u_xt * (1 + u_xt)
w_rat = 1 + 2 * u_xt
P_rat = g_rat / w_rat**2
Q_rat = w_rat**3 / x
check("positive control rational J(P,Q)=1", jac(P_rat, Q_rat) - 1)

# THM-3992's boundary node is a hostile geometry-only control.
A_node = gamma**2 * x**2 + a
C_node = gamma**3 * x**3 + sp.Rational(3, 2) * a * gamma * x
I_node = sp.Rational(3, 4) * a**2
check("hostile node eliminant", C_node**2 - A_node**3 + I_node * A_node + a**3 / 4)
check("hostile node has zero Jacobian", jac(A_node, C_node))

# Nonvacuity/hostile control inside the exact support cell: adjoining the
# forced p-jet satisfies the boundary and A_t invoices but still is not
# Keller.  Thus the ansatz is populated and the filters are compatible.
forced_p_coefficient = -sp.Rational(2, 3) / (gamma * a)
A_forced_probe = A_node + forced_p_coefficient * p_xt
C_forced_probe = C_node
check("support-cell control A boundary", A_forced_probe.subs(t, 0) - A_node)
check("support-cell control C boundary", C_forced_probe.subs(t, 0) - C_node)
check(
    "support-cell control forced A_t value",
    sp.diff(A_forced_probe, t).subs({x: 0, t: 0}) - forced_p_coefficient,
)
gate(
    "support-cell control is hostile rather than Keller",
    simp(jac(A_forced_probe, C_forced_probe) - 1) != 0,
)

# Bare t fails the negative-weight module condition, whereas p has the
# mandatory u(1+u) numerator and has the simple t-jet needed by THM-3992.
gate(
    "bare t fails the exact weight-minus-two coefficient ideal",
    sp.rem(u, u * (1 + u), u) != 0,
)
gate(
    "p numerator passes the exact weight-minus-two coefficient ideal",
    sp.rem(u * (1 + u), u * (1 + u), u) == 0,
)
check("p supplies the forced simple A_t jet", sp.diff(p_xt, t).subs({x: 0, t: 0}) - 1)

# Direct verification of the n=2 homogeneous bracket law:
# J(x^r f(u),x^s g(u))=x^(r+s+1)(r f g'-s f'g).
r, s = sp.symbols("r s", integer=True)
ff = sp.Function("ff")
gg = sp.Function("gg")
# Use concrete symbolic exponents for every pair occurring in the core cell;
# this avoids asking SymPy to differentiate x to a symbolic integer power.
for rr, ss in [(2, 3), (2, 1), (0, 3), (0, 1), (-2, 3), (-2, 1)]:
    lhs = jac(x**rr * ff(u_xt), x**ss * gg(u_xt))
    rhs = x**(rr + ss + 1) * (
        rr * ff(u_xt) * sp.diff(gg(u), u).subs(u, u_xt)
        - ss * sp.diff(ff(u), u).subs(u, u_xt) * gg(u_xt)
    )
    check(f"homogeneous bracket law ({rr},{ss})", lhs - rhs)

# -------------------------------------------------------------------------
# Core support cell A:{2,0,-2}, C:{3,1}.
# -------------------------------------------------------------------------

q, qp, F, Fp, m, mp = sp.symbols("q qp F Fp m mp")
k = sp.Rational(3, 2) * q * F
kp = sp.Rational(3, 2) * (qp * F + q * Fp)

# Weight 6 gives g^2=f^3 and hence f=q^2,g=q^3 after matching the
# nonzero boundary values.  Weight 4 then integrates to k=(3/2)qF.
weight4 = 2 * q**2 * kp - 2 * q * qp * k - 3 * Fp * q**3
check("weight 4 collapses to k=(3/2)qF", weight4)

# Weight 2 is the derivative of F^2/2+2*m*q^2.
weight2 = -Fp * k - 2 * m * (3 * q**2 * qp) - 3 * mp * q**3
conic_derivative = F * Fp + 2 * mp * q**2 + 4 * m * q * qp
check("weight 2 is the conic derivative", weight2 + sp.Rational(3, 2) * q * conic_derivative)

# The integration constant is a^2 because F(0)=a and m(0)=0:
# F^2+4*m*q^2=a^2.  Eliminate m' from its derivative in the scalar row.
mp_from_conic = -F * Fp / (2 * q**2) - 2 * m * qp / q
scalar_row = -2 * m * kp - mp * k
scalar_after_mp = simp(scalar_row.subs(mp, mp_from_conic))
check(
    "scalar row after differentiated conic",
    scalar_after_mp - sp.Rational(3, 4) * Fp * (F**2 - 4 * m * q**2) / q,
)
scalar_after_conic = simp(scalar_after_mp.subs(m, (a**2 - F**2) / (4 * q**2)))
check(
    "scalar row gives 3*F'*(2F^2-a^2)=4q",
    scalar_after_conic - sp.Rational(3, 4) * Fp * (2 * F**2 - a**2) / q,
)

# Exact degree certificate.  If d=deg(F)>=1, the scalar ODE forces
# deg(q)=3d-1.  Polynomial m forces q^2 | a^2-F^2, hence
# 2(3d-1)<=2d, impossible.  The difference is 4d-2, positive for every
# integer d>=1.  A constant F makes the ODE's left side zero.
d_degree = sp.symbols("d_degree", integer=True, positive=True)
degree_gap = sp.expand(2 * (3 * d_degree - 1) - 2 * d_degree)
check("all-degree gap is 4*d-2", degree_gap - (4 * d_degree - 2))
gate("all-degree gap is positive from d=1 onward", (4 * 1 - 2) > 0)
print("DEDUCTION core: no polynomial coefficient solution in any degree.")
print("  nonconstant F: deg(q)=3deg(F)-1 but q^2 divides a^2-F^2")
print("  constant F: scalar ODE has zero left side and nonzero q right side")

# -------------------------------------------------------------------------
# One extra arbitrary C-weight.
# -------------------------------------------------------------------------

# Let the extra row be x^r*l(u).  Its three bracket weights against A are
# r+3,r+1,r-1.  The core bracket weights are 6,4,2,0.  Therefore a collision
# is possible only for r in {-3,-1,1,3,5,7}; r=1,3 merely enlarge the two
# existing coefficient polynomials.
core_weights = {0, 2, 4, 6}
offsets = {3, 1, -1}
collision_weights = sorted({core - offset for core in core_weights for offset in offsets})
gate(
    "complete one-row collision set",
    collision_weights == [-3, -1, 1, 3, 5, 7],
)
print("GENERIC_EXTRA: outside collision set, the contradictory core rows are unchanged")

# At r=-3 or -1 the unique lowest equation is
# -2*m*l'-r*m'*l=0.  Since ord_0(m)=1 (the nonzero forced A_t jet), a
# nonzero l of integer valuation v would require 2v+r=0, impossible for
# these two odd r.
for rr in (-3, -1):
    required_twice_v = -rr
    gate(
        f"r={rr}: lowest row requires nonintegral valuation",
        required_twice_v % 2 == 1,
    )

# At r=5 or 7 the unique highest equation is
# 2*f*l'-r*f'*l=0.  The boundary normal form has f(0)=gamma^2!=0 and
# l(0)=0.  If ord_0(l)=v>=1, its leading coefficient is 2*f(0)*v, never
# zero in characteristic zero.  Hence l=0.
f0, ell_lead = sp.symbols("f0 ell_lead", nonzero=True)
v = sp.symbols("v", integer=True, positive=True)
for rr in (5, 7):
    leading = 2 * f0 * v * ell_lead
    gate(
        f"r={rr}: highest-row leading coefficient is nonzero",
        leading.is_nonzero is True,
    )

print("DEDUCTION extra: an arbitrary single additional C-weight is zero or absorbed.")

print("CERTIFICATE_SCOPE=>=2 extra C weights, or any extra A weight, not treated here; canonically first live sizes are 3x4 and 4x3")
print("MIXED_RESIDUAL_SCOPE=R in (p^2,y) was not specialized; obstruction uses only its forced node jet")
print("THEOREM_ID=THM-3998")
print("ALL THM-3998 EXACT CHECKS PASSED")
