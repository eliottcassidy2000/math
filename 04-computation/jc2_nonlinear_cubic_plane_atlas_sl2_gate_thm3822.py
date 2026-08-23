#!/usr/bin/env python3
"""Exact intrinsic SL2 and cheapest elementary-cell probe for THM-3822."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


GATES = 0


def require(condition: object, label: str) -> None:
    global GATES
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    GATES += 1


def equal(left: object, right: object, label: str) -> None:
    difference = left-right  # type: ignore[operator]
    require(sp.cancel(sp.expand(difference)) == 0, label)


# Intrinsic coordinates.  On THM-3811's U they are
# c=C, ell=3*theta+14*A, h=A/D, k=omega/D.
c, ell, h, k, D = sp.symbols("c ell h k D")
det = c*k-h*ell-1
Q = 3*k**2+7*h**2
N = 3+2*h*ell
R = h*(14*k+9*h)
S = k*ell+3*h*c**2
Psi = sp.expand(N*R-Q*S)

# D is the unique simultaneous solution QD=N, RD=S.  Compatibility is Psi.
E1 = Q*D-N
E2 = R*D-S
equal(sp.resultant(E1, E2, D), Psi,
      "SL2 hypersurface is the D-compatibility resultant")

# Q and R generate the unit ideal after imposing det=0.  The two small
# homogeneous certificates put h^3 and k^3 in (Q,R); (ck-h ell)^5=1 then
# puts 1 in that ideal.  A direct Groebner gate guards the gluing argument.
equal(196*h*Q-3*(14*k-9*h)*R, 1615*h**3,
      "h-cube Bezout certificate")
equal((1615*k+882*h)*Q+(-189*k-686*h)*R, 4845*k**3,
      "k-cube Bezout certificate")
unit_basis = sp.groebner([det, Q, R], D, ell, c, k, h, order="grevlex")
require(len(unit_basis.polys) == 1 and unit_basis.polys[0].as_expr() == 1,
        "Q and R cover the intrinsic SL2 surface")

# The two unimodular charts h!=0 and k!=0 are domains.  After eliminating
# ell or c, each is quadratic with the same nonconstant squarefree
# discriminant (up to a square k^2).  Hence the hypersurface has no hidden
# component on which the reconstruction could fail.
phi_h = sp.together(Psi.subs(ell, (c*k-1)/h)).as_numer_denom()[0]
disc_h = sp.factor(sp.discriminant(phi_h, c))
disc_core = sp.factor(disc_h/9)
require(sp.Poly(phi_h, c).degree() == 2, "h-chart quadratic degree")
require(sp.gcd(sp.Poly(disc_core, h, k),
               sp.Poly(sp.diff(disc_core, h), h, k)) == 1,
        "h-chart discriminant squarefree in h")
require(sp.gcd(sp.Poly(disc_core, h, k),
               sp.Poly(sp.diff(disc_core, k), h, k)) == 1,
        "h-chart discriminant squarefree in k")
phi_k = sp.together(Psi.subs(c, (1+h*ell)/k)).as_numer_denom()[0]
equal(sp.discriminant(phi_k, ell), 9*k**2*disc_core,
      "k-chart has same nonsquare discriminant")

# Reconstruct the original nonlinear cubic packet and check all three laws.
A = h*D
W = k*D
Theta = (ell-14*h*D)/3
F = (
    W**2+7*A**2-c*W+A*Theta,
    W*Theta-3*A**2+A*c**2,
    Theta**2-3*A*c+c**3-(c**2-3*A)*W+7*A*Theta,
)
reconstruction_basis = sp.groebner(
    [det, E1, E2], D, ell, c, k, h, order="grevlex"
)
for index, relation in enumerate(F, 1):
    remainder = reconstruction_basis.reduce(sp.expand(9*relation))[1]
    equal(remainder, 0, f"reconstructed cubic law {index}")

# On h=0 the intrinsic equations force ck=1 and ell=0.  Also QD=N gives
# D=k^-2=c^2.  Therefore an etale plane atlas with base pair (A,c) and
# Jacobian kappa must satisfy c^2 J(h,c)=kappa modulo h.
equal(det.subs(h, 0), c*k-1, "companion determinant")
equal(Psi.subs(h, 0), -3*k**3*ell, "companion ell vanishing")
equal(E1.subs(h, 0), 3*(k**2*D-1), "companion D=c^2")

# Cheapest standard elementary big cell:
# E_+(p)E_-(s)=[[1+ps,p],[s,1]].  The intrinsic equation is quadratic in p.
p, s = sp.symbols("p s")
cell = {c: 1+p*s, ell: p, h: s, k: 1}
equal(det.subs(cell), 0, "two-shear determinant")
cell_equation = sp.factor(-Psi.subs(cell)/3)
cell_expected = (
    s**3*(7*s**2+3)*p**2
    +(14*s**4-6*s**3-s**2+1)*p
    +s*(7*s**2-9*s-11)
)
equal(cell_equation, cell_expected, "two-shear intrinsic quadratic")
cell_discriminant = sp.factor(sp.discriminant(cell_expected, p))
disc_expected = (
    84*s**7+232*s**6+120*s**5+161*s**4-12*s**3-2*s**2+1
)
equal(cell_discriminant, disc_expected, "two-shear square discriminant")
require(sp.gcd(sp.Poly(disc_expected, s),
               sp.Poly(sp.diff(disc_expected, s), s)) == 1,
        "discriminant polynomial squarefree")
require(sp.Poly(disc_expected, s).degree() == 7,
        "discriminant has at least two distinct roots")
require(sp.gcd_list([
    s**3*(7*s**2+3),
    14*s**4-6*s**3-s**2+1,
    s*(7*s**2-9*s-11),
]) == 1, "constant-s coefficient packet is nonzero")

# Hostile positive control: determinant + square lift + punctured-arm
# congruence alone are nonempty.  This matrix deliberately fails Psi, proving
# the intrinsic hypersurface is load-bearing.
x, y = sp.symbols("x y")
h0 = x**4*y-1
c0 = x**3*y
k0 = x*(1-h0)
ell0 = -h0
J_hc = sp.diff(h0, x)*sp.diff(c0, y)-sp.diff(h0, y)*sp.diff(c0, x)
equal(c0*k0-h0*ell0, 1, "punctured-arm SL2 positive control")
equal(J_hc, x**6*y, "punctured-arm Jacobian")
require(sp.rem(sp.Poly(c0**2*J_hc-1, y), sp.Poly(h0, y)) == 0,
        "c^2 J(h,c)=1 modulo h positive control")
hostile_Psi = sp.factor(Psi.subs({c: c0, ell: ell0, h: h0, k: k0}))
require(hostile_Psi != 0 and sp.rem(sp.Poly(hostile_Psi, y),
                                   sp.Poly(h0, y)) == 0,
        "positive control fails only the full intrinsic equation")

source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert)
                for node in ast.walk(ast.parse(source))),
        "no inactive Python asserts")

semantic = {
    "intrinsic": "U is the SL2 hypersurface ck-h*ell=1 and Psi=0; D is uniquely reconstructed from QD=N,RD=S",
    "atlas": "reconstruct A=hD; a polynomial solution is etale exactly when J(A,c) is a nonzero constant",
    "arm": "h=0 gives ck=1,ell=0,D=c^2 and c^2*J(h,c)=kappa mod h",
    "cell": "E_+(p)E_-(s) forces a quadratic whose squarefree degree-seven discriminant makes s, then p, constant",
    "scope": "determinant and arm congruence alone do not obstruct; arbitrary SL2 need not lie in E2",
}
semantic_hash = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()
print("object=THM-3822-intrinsic-SL2-and-punctured-arm-gate")
print("intrinsic=ck-h*ell=1;Psi=(3+2h*ell)h(14k+9h)-(3k^2+7h^2)(k*ell+3hc^2)=0")
print("reconstruction=QD=3+2h*ell;h(14k+9h)D=k*ell+3hc^2;A=hD")
print("arm=h=0=>ck=1,ell=0,D=c^2;c^2J(h,c)=kappa_mod_h")
print(f"two_shear_discriminant={disc_expected}")
print("two_shear_verdict=NO_DOMINANT_ATLAS;s_and_p_forced_constant")
print("hostile=arm_and_square_lift_nonempty_but_Psi_load_bearing")
print("scope=standard_Eplus_Eminus_cell_only;SL2_not_equal_E2")
print(f"GATES={GATES}")
print(f"semantic_sha256={semantic_hash}")
