#!/usr/bin/env python3
"""Exact identities for THM-3573's polynomial Pell-parameter descent."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


a, b, X, H = sp.symbols("a b X H")

# The normalized polynomial compiler q=1/H.
U = 1 + b * H + 4 * a * H**2
phi = 4 * H * U
q = 1 / H
compiler = 4 * (q**2 + b * q + 4 * a) / q**3
require(sp.cancel(compiler - phi) == 0, "normalized polynomial compiler")
require(
    sp.cancel(phi * q**3 - 4 * q**2 - 4 * b * q - 16 * a) == 0,
    "THM-3570 compiler equation",
)

# The compiler does more than supply a rational root: the cubic has a compact
# polynomial factorization over C[a,b,H].
R = 1 + 2 * b * H + 12 * a * H**2
M = 48 * a**2 * H**2 + 8 * a * b * H + 16 * a - b**2
L_phi = sp.expand(27 * a**2 * phi**2 + 18 * a * b * phi + 16 * a - b**3 * phi - b**2)
E_phi = sp.expand(L_phi * X**3 + (4 + 3 * b * phi) * X + 2 * phi)
factor_packet = (R * X + 2 * H) * (R * M * X**2 - 2 * H * M * X + 4 * U)
require(sp.expand(L_phi - R**2 * M) == 0, "leading coefficient square packet")
require(sp.expand(E_phi - factor_packet) == 0, "full polynomial factorization")
require(sp.cancel(E_phi.subs(X, -2 * H / R)) == 0, "compiled rational root")

# The quadratic cofactor remembers a nonsquare sidecar.
cofactor = R * M * X**2 - 2 * H * M * X + 4 * U
cofactor_disc = sp.factor(sp.discriminant(cofactor, X))
cofactor_square = 2 + 3 * b * H + 12 * a * H**2
require(
    sp.expand(cofactor_disc + 4 * M * cofactor_square**2) == 0,
    "quadratic cofactor discriminant",
)

# The only apparent nonunit numerator branch is q=(c*a)/D.  Expanding the
# divisibility obstruction through a^1 pins the exact contradiction used in
# the UFD proof.
c, D0, D1, D2 = sp.symbols("c D0 D1 D2", nonzero=True)
D = D0 + a * D1 + a**2 * D2
obstruction = sp.Poly(c**2 * a + b * c * D + 4 * D**2, a)
constant_row = obstruction.coeff_monomial(1)
linear_row = obstruction.coeff_monomial(a)
require(sp.expand(constant_row - D0 * (b * c + 4 * D0)) == 0, "N=c*a constant row")
require(
    sp.expand(linear_row - (c**2 + (b * c + 8 * D0) * D1)) == 0,
    "N=c*a linear row",
)
require(
    sp.expand(linear_row.subs(D0, -b * c / 4) - c * (c - b * D1)) == 0,
    "N=c*a impossible second coefficient",
)

# Exact degree law: deg_a(phi)=3*deg_a(H)+1.  These controls exercise every
# residue class representative m=0,...,12 without assuming the theorem.
for m in range(13):
    Hm = b + 1 if m == 0 else a**m + b + 1
    phim = sp.Poly(sp.expand(phi.subs(H, Hm)), a)
    require(phim.degree() == 3 * m + 1, f"degree law m={m}")

# The complete first new resonance, H=u(b)*a+v(b), and a hostile with exactly
# one rational linear factor.
u, v = sp.symbols("u v", nonzero=True)
H4 = u * a + v
phi4 = sp.Poly(sp.expand(phi.subs(H, H4)), a)
expected_phi4 = (
    16 * u**3 * a**4
    + 48 * u**2 * v * a**3
    + 4 * u * (b * u + 12 * v**2) * a**2
    + 4 * (2 * b * u * v + u + 4 * v**3) * a
    + 4 * v * (b * v + 1)
)
require(sp.expand(phi4.as_expr() - expected_phi4) == 0, "complete degree-four family")

phi_a = sp.expand(phi.subs(H, a))
require(sp.expand(phi_a - (4 * a + 4 * b * a**2 + 16 * a**4)) == 0, "H=a control")
M_a = sp.expand(M.subs(H, a))
minus_M_a = sp.Poly(-M_a, b)
require(
    sp.expand(sp.discriminant(minus_M_a.as_expr(), b) - 64 * a * (4 * a**3 + 1)) == 0,
    "H=a nonsquare cofactor control",
)

# Regression to THM-3565's degree-one family.
h = sp.symbols("h", nonzero=True)
phi_h = sp.factor(phi.subs(H, -h / 2))
require(
    sp.expand(phi_h - (-2 * h**3 * a + b * h**2 - 2 * h)) == 0,
    "THM-3565 specialization",
)

print("THM-3573 polynomial Pell-parameter descent audit")
print("classification: q=1/H, phi=4*H*(1+b*H+4*a*H^2)")
print("degree law: deg_a(phi)=3*deg_a(H)+1")
print("L_phi=R^2*M and E_phi=(R*X+2*H)*(R*M*X^2-2*H*M*X+4*U)")
print("degree-four family: H=u(b)*a+v(b), u!=0")
print("hostile H=a: phi=4*a+4*b*a^2+16*a^4")
print("N=c*a obstruction: D0=-b*c/4 then D1=c/b, not polynomial")
print("all active truth gates passed")
