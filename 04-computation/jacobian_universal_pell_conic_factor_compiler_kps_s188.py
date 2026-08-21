#!/usr/bin/env python3
"""Exact identities for THM-3570's universal Pell-conic compiler."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


a, b, x, phi, q = sp.symbols("a b x phi q", nonzero=True)
t, W = sp.symbols("t W")

# Regard the core cubic equation as a quadratic in phi.
A_phi = 27 * a**2 * x**3
B_phi = (18 * a * b - b**3) * x**3 + 3 * b * x + 2
C_phi = (16 * a - b**2) * x**3 + 4 * x
disc_phi = sp.factor(B_phi**2 - 4 * A_phi * C_phi)
R = 12 * a * x**2 - b**2 * x**2 - 2 * b * x - 1
S = 12 * a * x**2 - b**2 * x**2 + b * x + 2
require(sp.expand(disc_phi + R * S**2) == 0, "quadratic-in-phi discriminant factorization")

# The exceptional divisor S=0 has nonsquare x-discriminant, so it has no
# rational point over C(a,b).  The symbolic row pins the odd divisor.
disc_x_S = sp.factor(sp.discriminant(S, x))
require(sp.expand(disc_x_S - 3 * (3 * b**2 - 32 * a)) == 0, "S=0 nonsquare discriminant")

# Parametrize W^2=-R by the line W=1+t*x through the rational point (0,1).
conic = sp.expand(W**2 + R)
x_t = sp.factor(2 * (b - t) / (t**2 - b**2 + 12 * a))
W_t = sp.factor(1 + t * x_t)
require(sp.factor(conic.subs({x: x_t, W: W_t})) == 0, "Pell conic parametrization")

# The two quadratic-formula values of phi in the t-chart.
phi_plus = sp.factor(
    sp.cancel(((-B_phi + S * W) / (2 * A_phi)).subs({x: x_t, W: W_t}))
)
phi_minus = sp.factor(
    sp.cancel(((-B_phi - S * W) / (2 * A_phi)).subs({x: x_t, W: W_t}))
)
q_sub = t - b
expected_plus = q_sub * (q_sub**2 + 3 * b * q_sub + 36 * a) / (108 * a**2)
expected_minus = 4 * (q_sub**2 + b * q_sub + 4 * a) / q_sub**3
require(sp.cancel(phi_plus - expected_plus) == 0, "first conic branch")
require(sp.cancel(phi_minus - expected_minus) == 0, "second conic branch")

# The involution q -> 12a/q identifies the two branches.
compiler_phi = 4 * (q**2 + b * q + 4 * a) / q**3
compiler_x = -2 * q / (q**2 + 2 * b * q + 12 * a)
require(
    sp.expand(sp.discriminant(q**2 + 2 * b * q + 12 * a, q) - 4 * (b**2 - 12 * a)) == 0,
    "compiler-root denominator nonsquare gate",
)
involuted_phi = sp.factor(compiler_phi.subs(q, 12 * a / q))
first_q = q * (q**2 + 3 * b * q + 36 * a) / (108 * a**2)
require(sp.cancel(involuted_phi - first_q) == 0, "branch involution")

# Direct converse: the compiler always supplies a core root.
L_compiler = sp.expand(
    27 * a**2 * compiler_phi**2
    + 18 * a * b * compiler_phi
    + 16 * a
    - b**3 * compiler_phi
    - b**2
)
E_compiler = sp.expand(
    L_compiler * x**3 + (4 + 3 * b * compiler_phi) * x + 2 * compiler_phi
)
require(sp.factor(E_compiler.subs(x, compiler_x)) == 0, "compiler root identity")

# The equivalent cubic divisibility equation is often the cheapest search
# form for polynomial target graphs.
compiler_equation = sp.factor(phi * q**3 - 4 * q**2 - 4 * b * q - 16 * a)
require(
    sp.cancel(compiler_equation.subs(phi, compiler_phi)) == 0,
    "compiler cubic equation",
)

# Regression to THM-3565's complete first-resonance family.
h = sp.symbols("h", nonzero=True)
h_phi = -2 * h**3 * a + b * h**2 - 2 * h
h_q = -2 / h
h_root = h / (3 * a * h**2 - b * h + 1)
require(sp.cancel(compiler_phi.subs(q, h_q) - h_phi) == 0, "h-family specialization")
require(sp.cancel(compiler_x.subs(q, h_q) - h_root) == 0, "h-family root specialization")

# Independent rational controls, including q depending on a and b.
controls = [sp.Integer(1), b + 1, a + b, (a + 1) / (b + 2), 12 * a / (b - 3)]
for index, q_value in enumerate(controls):
    phi_value = sp.cancel(compiler_phi.subs(q, q_value))
    x_value = sp.cancel(compiler_x.subs(q, q_value))
    L_value = sp.expand(
        27 * a**2 * phi_value**2
        + 18 * a * b * phi_value
        + 16 * a
        - b**3 * phi_value
        - b**2
    )
    E_value = sp.cancel(L_value * x_value**3 + (4 + 3 * b * phi_value) * x_value + 2 * phi_value)
    require(E_value == 0, f"rational compiler control {index}")

print("THM-3570 universal Pell-conic target-graph factor audit")
print("disc_phi=-(12*a*x^2-b^2*x^2-2*b*x-1)*(12*a*x^2-b^2*x^2+b*x+2)^2")
print("Pell conic: W^2=(1+b*x)^2-12*a*x^2")
print("parameter: q=t-b, x=-2*q/(q^2+2*b*q+12*a)")
print("compiler: phi=4*(q^2+b*q+4*a)/q^3")
print("other branch is the same compiler under q -> 12*a/q")
print("equivalent equation: phi*q^3-4*q^2-4*b*q-16*a=0")
print("h-family regression: q=-2/h")
print("all active truth gates passed")
