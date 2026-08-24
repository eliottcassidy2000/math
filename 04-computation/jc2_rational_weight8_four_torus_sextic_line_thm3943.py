#!/usr/bin/env python3
"""Exact companion for THM-3943's rational weight-eight line obstruction.

Reproduction:
  python3 04-computation/jc2_rational_weight8_four_torus_sextic_line_thm3943.py
  python3 -O 04-computation/jc2_rational_weight8_four_torus_sextic_line_thm3943.py
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


u, x, ybar = sp.symbols("u x ybar")


# Degtyarev's four-cuspidal trigonal quotient and its normalization.
f = 4*ybar**3 - (24*x**3 + 3)*ybar + 8*x**6 + 20*x**3 - 1
x_u = 3*u/(u**3 + 2)
ybar_u = -(u**6 - 20*u**3 - 8)/(2*(u**3 + 2)**2)
gate(sp.factor(f.subs({x: x_u, ybar: ybar_u})) == 0,
     "four-cuspidal trigonal parametrization lies on f")


# The simple weight-eight list and the genus-zero sublist.  The classification
# and four-torus assertion are cited in the theorem; only the delta ledger is
# repo-derived here.
simple_weight8_delta = {
    "E6+A5+4A2": 3 + 3 + 4,
    "E6+6A2": 3 + 6,
    "2A5+4A2": 2*3 + 4,
    "A5+6A2+A1": 3 + 6 + 1,
    "A5+6A2": 3 + 6,
    "8A2+A1": 8 + 1,
    "8A2": 8,
}
gate(len(simple_weight8_delta) == 7, "seven cited simple weight-eight configurations")
gate([name for name, delta in simple_weight8_delta.items() if delta == 10] ==
     ["E6+A5+4A2", "2A5+4A2", "A5+6A2+A1"],
     "exactly three simple weight-eight configurations are rational")


# ---------------------------------------------------------------------------
# Rigid E6+A5+4A2: sparse normalization coordinates.
# ---------------------------------------------------------------------------

s = sp.symbols("s")
section_rigid = 2*x**2 - sp.Rational(1, 2)
rigid_rad = sp.factor(ybar_u - section_rigid.subs(x, x_u))
gate(rigid_rad == 6*(u - 1)**2*(2*u + 1)/(u**3 + 2)**2,
     "rigid double-cover radicand")

u_s = (s**2 - 1)/2
x_s = sp.factor(x_u.subs(u, u_s))
y_s = sp.factor(sp.sqrt(6)*s*(u_s - 1)/(u_s**3 + 2))
gate(sp.factor(y_s**2 - rigid_rad.subs(u, u_s)) == 0,
     "rigid normalization square root")

D_s = s**6 - 3*s**4 + 3*s**2 + 15
gate(sp.factor(x_s - 12*(s**2 - 1)/D_s) == 0 and
     sp.factor(y_s - 4*sp.sqrt(6)*s*(s**2 - 3)/D_s) == 0,
     "rigid affine normalization coordinates")
gate(sp.factor(f.subs({x: x_s, ybar: y_s**2 + section_rigid.subs(x, x_s)})) == 0,
     "rigid normalization lies on the sextic")

S, T = sp.symbols("S T")
Xr = 12*(S**2 - T**2)*T**4
Yr = 4*sp.sqrt(6)*S*(S**2 - 3*T**2)*T**3
Zr = S**6 - 3*S**4*T**2 + 3*S**2*T**4 + 15*T**6
gate(sp.factor(sp.gcd(sp.gcd(Xr, Yr), Zr)) == 1,
     "rigid projective normalization is basepoint-free")
gate(all(sp.expand(expr).coeff(S, 5).coeff(T, 1) == 0 for expr in (Xr, Yr, Zr)),
     "rigid line pullbacks have no S5T coefficient")


def binary_coeff_vector(expr: sp.Expr) -> sp.Matrix:
    poly = sp.Poly(sp.expand(expr), S, T)
    return sp.Matrix([poly.coeff_monomial(S**j*T**(6-j)) for j in range(7)])


rigid_span = sp.Matrix.hstack(*(binary_coeff_vector(expr) for expr in (Xr, Yr, Zr)))
gate(rigid_span.rank() == 3, "rigid line-pullback plane has dimension three")
gate(sp.Matrix.hstack(rigid_span, binary_coeff_vector(S**6)).rank() == 4,
     "pure S6 is not a rigid line pullback")
gate(sp.Matrix.hstack(rigid_span, binary_coeff_vector(T**6)).rank() == 4,
     "pure T6 is not a rigid line pullback")
# A mixed sixth power has nonzero S5T coefficient 6*p^5*q.  Therefore the
# preceding two endpoint tests exhaust all pure sixth powers.
p_root, q_root = sp.symbols("p_root q_root")
gate(sp.Poly((p_root*S + q_root*T)**6, S, T).coeff_monomial(S**5*T) ==
     6*p_root**5*q_root, "mixed sixth power exposes the missing S5T slot")


# ---------------------------------------------------------------------------
# 2A5+4A2: a quadratic normalization and a norm contradiction.
# ---------------------------------------------------------------------------

a, alpha, B, root = sp.symbols("a alpha B root")
section_two_a5 = a*x**2 + (2-a)*x - sp.Rational(1, 2)
Q_a = (a-2)*u**2 + 2*a*u + 2
two_a5_rad = sp.factor(ybar_u - section_two_a5.subs(x, x_u))
gate(sp.factor(two_a5_rad - 3*(u-1)**2*Q_a/(u**3+2)**2) == 0,
     "2A5+4A2 double-cover radicand")
gate(sp.discriminant(Q_a, u) == 4*(a**2 - 2*a + 4),
     "two-A5 split-conic degeneration polynomial")
gate(Q_a.subs(u, 1) == 3*a and sp.LC(sp.Poly(Q_a, u)) == a-2,
     "a=0 and a=2 are the two ordered rigid cusp-tangency endpoints")

Z = u**3 + 2
norm_two_a5 = sp.expand((Z + 3*alpha*u)**2 - B*(u-1)**2*Q_a)
P_two = sp.Poly(norm_two_a5, u)
gate(P_two.coeff_monomial(u**6) == 1 and P_two.coeff_monomial(u**5) == 0,
     "normalized two-A5 line norm is monic with zero quintic term")
gate(P_two.coeff_monomial(u**3) == 4*(1-B),
     "two-A5 cubic coefficient forces B=1 after root=0")
gate(P_two.coeff_monomial(1) == 4-2*B,
     "two-A5 constant coefficient contradicts B=1")
gate(sp.Poly((u-root)**6, u).coeff_monomial(u**5) == -6*root,
     "a finite one-place norm first forces root zero")
gate(P_two.coeff_monomial(u**3).subs(B, 1) == 0 and
     P_two.coeff_monomial(1).subs(B, 1) == 2,
     "two-A5 finite one-place coefficient contradiction is uniform in a")

# For a!=2 the conic v^2=3Q_a has two points at W=0.  In degree-three
# homogeneous coordinates both map to [0:0:1], so every gamma=0 line contains
# both normalization addresses.  At a=2 the family is the rigid row above.
U, V, W = sp.symbols("U V W")
X2h = 3*U*W**2
Y2h = (U-W)*V*W
Z2h = U**3 + 2*W**3
gate(X2h.subs(W, 0) == 0 and Y2h.subs(W, 0) == 0 and
     Z2h.subs(W, 0) == U**3,
     "both two-A5 conic infinity addresses map to the same plane point")
gate(section_two_a5.subs(a, 2) == section_rigid,
     "a=2 is the displayed rigid endpoint")


# ---------------------------------------------------------------------------
# A5+6A2+A1: tangent parameter tau, finite norm and infinity seams.
# ---------------------------------------------------------------------------

tau, lam = sp.symbols("tau lam")
den_tau = (tau-1)*(tau+2)**2
a_tau = 2*(tau**3 - 3*tau - 1)/den_tau
b_tau = 6*tau*(tau+1)/den_tau
c_tau = -(tau**3 + 3*tau**2 + 8)/(2*den_tau)
section_a5_a1 = a_tau*x**2 + b_tau*x + c_tau
Q_tau = u**2 + 2*(tau+1)*u + tau + 3
a5_a1_rad = sp.factor(ybar_u - section_a5_a1.subs(x, x_u))
expected_rad = (6*(tau-u)**2*(u-1)**2*Q_tau /
                ((tau-1)*(tau+2)**2*(u**3+2)**2))
gate(sp.factor(a5_a1_rad - expected_rad) == 0,
     "A5+6A2+A1 double-cover radicand")
gate(sp.factor(sp.discriminant(Q_tau, u) - 4*(tau-1)*(tau+2)) == 0,
     "A5+A1 conic is smooth at every advertised parameter")

R_tau = (u-tau)**2*(u-1)**2*Q_tau
norm_a5_a1 = sp.expand((Z + 3*alpha*u)**2 - B*R_tau)
P_tau = sp.Poly(norm_a5_a1, u)
gate(P_tau.coeff_monomial(u**6) == 1-B and
     P_tau.coeff_monomial(u**5) == 0,
     "A5+A1 normalized line-norm leading packet")
tau_cubic = P_tau.coeff_monomial(u**3)
tau_constant = P_tau.coeff_monomial(1)
gate(sp.factor(tau_cubic/2 - tau_constant - 2*(B-1)) == 0,
     "finite one-place root zero forces B=1")
gate((1-B).subs(B, 1) == 0,
     "B=1 contradicts the nonzero leading norm of a finite one-place line")

# A line with zero Z coefficient and beta!=0 has norm below.  A finite pure
# sixth power forces root=0; its u4, u3, and constant equations have no common
# tau.  beta=0 is X=0 and meets both infinity points.
norm_gamma_zero = sp.expand(9*alpha**2*u**2 - B*R_tau)
P_zero = sp.Poly(norm_gamma_zero, u)
gate(P_zero.coeff_monomial(u**6) == -B and
     P_zero.coeff_monomial(u**5) == 0,
     "gamma-zero finite norm has nonzero leading coefficient when beta is nonzero")
e4 = sp.factor(P_zero.coeff_monomial(u**4)/B)
e3 = sp.factor(P_zero.coeff_monomial(u**3)/(2*B))
e0 = sp.factor(P_zero.coeff_monomial(1)/B)
gate(sp.factor(e4 - 3*tau*(tau+1)) == 0 and
     sp.factor(e3 - (-tau**3-3*tau**2+2)) == 0 and
     sp.factor(e0 + tau**2*(tau+3)) == 0,
     "gamma-zero finite one-place compatibility equations")
gate(sp.gcd(sp.gcd(e4, e3), e0) == 1,
     "gamma-zero finite compatibility ideal is empty")

# If a nonzero-Z line has all zeros at one infinity address, its norm is a
# nonzero constant.  Vanishing of the degree-six, degree-four, and degree-three
# terms gives B=1, 2alpha+tau(tau+1)=0, and the excluded parameter polynomial.
gate(sp.factor(P_tau.coeff_monomial(u**4).subs(B, 1)/3) ==
     2*alpha + tau**2 + tau,
     "constant norm fixes the affine line slope")
gate(sp.factor(P_tau.coeff_monomial(u**3).subs(B, 1) +
               2*(tau-1)*(tau+2)**2) == 0,
     "constant norm exists only at excluded tangent parameters")

X3h = 3*U*W**2
Y3h = (U-tau*W)*(U-W)*V
Z3h = U**3 + 2*W**3
gate(X3h.subs(W, 0) == 0 and Y3h.subs(W, 0) == U**2*V and
     Z3h.subs(W, 0) == U**3,
     "A5+A1 two conic infinity addresses map to distinct plane points")
Q_tau_h = sp.expand(W**2*Q_tau.subs(u, U/W))
gate(Q_tau_h.subs(W, 0) == U**2,
     "A5+A1 conic has infinity signs V=plus-or-minus U")


summary = {
    "checks": CHECKS,
    "cited_scope": "simple weight8;four torus structures;seven configurations",
    "rational_simple": "E6+A5+4A2;2A5+4A2;A5+6A2+A1",
    "rigid": "missing S5T reduces sixth powers to two impossible endpoints",
    "two_A5": "finite norm coefficient 2;gamma0 has two addresses",
    "A5_A1": "finite norms empty;constant norm only tau=1,-2",
    "boundary": "non-simple J2,3+3A2 and J2,0+4A2 not classified",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3943 rational weight-eight four-torus line obstruction")
print(f"CHECKS={CHECKS}")
print("CITED_SCOPE=simple weight8;four torus structures;seven configurations")
print("RATIONAL_SIMPLE=E6+A5+4A2;2A5+4A2;A5+6A2+A1")
print("RIGID=missing S5T reduces sixth powers to two impossible endpoints")
print("TWO_A5=finite norm leaves constant 2;gamma0 has two addresses")
print("A5_A1=finite norms empty;constant norm only tau=1,-2")
print("BOUNDARY=non-simple J2,3+3A2 and J2,0+4A2 not classified")
print(f"SEMANTIC_SHA256={semantic}")
