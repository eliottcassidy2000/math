#!/usr/bin/env python3
"""Universal source-factor/unit identities for THM-3574."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


x, y, z, H = sp.symbols("x y z H")
a, b, X = sp.symbols("a b X")

u = 1 + x * y
F1 = sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y))
F2 = sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
F3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)

# Universal target packet.  H is first an independent indeterminate; any
# polynomial H(a,b) may subsequently be substituted, then a=F1,b=F2.
phi = 4 * H * (1 + b * H + 4 * a * H**2)
U = 1 + b * H + 4 * a * H**2
R = 1 + 2 * b * H + 12 * a * H**2
M = 48 * a**2 * H**2 + 8 * a * b * H + 16 * a - b**2
V = 1 - b * H + 12 * a * H**2

require(sp.expand(H**2 * M + 1 - U * V) == 0, "target inverse-unit identity")

L_phi = sp.expand(27 * a**2 * phi**2 + 18 * a * b * phi + 16 * a - b**3 * phi - b**2)
E_phi = sp.expand(L_phi * X**3 + (4 + 3 * b * phi) * X + 2 * phi)
core_packet = sp.expand((R * X + 2 * H) * (R * M * X**2 - 2 * H * M * X + 4 * U))
require(sp.expand(L_phi - R**2 * M) == 0, "target Jelonek factor packet")
require(sp.expand(E_phi - core_packet) == 0, "target core factor packet")

# Substitute the fixed source map while retaining H as an indeterminate.
phi_source = sp.expand(phi.subs({a: F1, b: F2}))
P = sp.expand(F3 + phi_source)
S = sp.expand(x + 2 * u * H)
quotient, remainder = sp.div(P, S, x)
Q = sp.expand(quotient)
require(sp.expand(remainder) == 0 and sp.expand(P - S * Q) == 0, "universal source factorization")

UF = sp.expand(U.subs({a: F1, b: F2}))
RF = sp.expand(R.subs({a: F1, b: F2}))
MF = sp.expand(M.subs({a: F1, b: F2}))
VF = sp.expand(V.subs({a: F1, b: F2}))

# R is a unit modulo every irreducible factor of S, with inverse u.
linear_unit_quotient, linear_unit_remainder = sp.div(sp.expand(RF * u - 1), S, x)
require(sp.expand(linear_unit_remainder) == 0, "R*u=1 modulo S")

# The target cofactor relation vanishes modulo Q.  Combining it with
# H^2*M+1=U*V yields the displayed polynomial inverse K of M modulo Q.
target_quadratic = sp.expand(MF * (RF * x**2 - 2 * H * x) + 4 * UF)
quadratic_relation_quotient, quadratic_relation_remainder = sp.div(target_quadratic, Q, z)
require(sp.expand(quadratic_relation_remainder) == 0, "target quadratic relation modulo Q")
require(
    len(sp.Poly(quadratic_relation_quotient, x, y, z, H).terms()) == 54,
    "target quadratic quotient term count",
)

K = sp.expand(-H**2 - sp.Rational(1, 4) * VF * (RF * x**2 - 2 * H * x))
unit_quotient, unit_remainder = sp.div(sp.expand(1 - MF * K), Q, z)
require(sp.expand(unit_remainder) == 0, "M*K=1 modulo Q")
require(sp.Poly(unit_quotient, x, y, z, H, domain=sp.QQ).total_degree() == 29, "inverse certificate degree")

# Degree controls instantiate the general leading-term proof.
degree_rows = []
for m in range(6):
    Hm = (a**m if m else 1) + b + 1
    Rm = sp.expand(R.subs(H, Hm))
    Mm = sp.expand(M.subs(H, Hm))
    row = (m, sp.degree(Rm, a), sp.degree(Mm, a))
    require(row == (m, 2 * m + 1, 2 * m + 2), f"target degree row m={m}")
    degree_rows.append(row)

# The old h(b) family is the specialization H=-h/2.
h = sp.symbols("h")
D_old = 3 * a * h**2 - b * h + 1
C_old = 12 * a**2 * h**2 - 4 * a * b * h + 16 * a - b**2
require(sp.expand(R.subs(H, -h / 2) - D_old) == 0, "THM-3568 R regression")
require(sp.expand(M.subs(H, -h / 2) - C_old) == 0, "THM-3568 M regression")
require(sp.expand(S.subs(H, -h / 2) - (x - u * h)) == 0, "THM-3565 source-factor regression")

# H=0 is deliberately excluded: the linear factor is the coordinate x and
# its purported target unit R collapses to the constant one.
require(sp.expand(P.subs(H, 0) - x * (2 - 3 * x * y - x**2 * z)) == 0, "zero-row hostile factorization")
require(R.subs(H, 0) == 1, "zero-row hostile constant unit")

print("THM-3574 universal reducible target-graph component unit audit")
print("source factorization: F3+phi_H(F1,F2)=S_H*Q_H")
print("S_H=x+2*(1+xy)*H(F1,F2)")
print("R(F1,F2)*(1+xy)=1 mod S_H")
print("H^2*M+1=U*V on the target")
print("M(F1,F2)*K_H=1 mod Q_H")
print(
    "relation/inverse certificate term counts:",
    len(sp.Poly(linear_unit_quotient, x, y, z, H).terms()),
    len(sp.Poly(quadratic_relation_quotient, x, y, z, H).terms()),
    len(sp.Poly(unit_quotient, x, y, z, H).terms()),
)
print("degree rows (deg_a H,deg_a R,deg_a M):", degree_rows)
print("H=0 hostile: S=x and R=1, so nonzero H is load-bearing")
print("verdict: every irreducible component of every nonzero reducible graph pullback has a nonconstant unit")
print("all active truth gates passed")
