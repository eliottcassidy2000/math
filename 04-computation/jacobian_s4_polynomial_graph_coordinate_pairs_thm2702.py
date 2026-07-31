#!/usr/bin/env python3
"""Exact companion for THM-2702 polynomial graph Keller slices."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


x, y, T = sp.symbols("x y T")
kappa, c0 = sp.symbols("kappa c0", nonzero=True)
f = sp.Function("f")(x, y)

A = x**2 - 2 * y
B = y**2 - 2 * x * f
d = f

J_Ad = sp.factor(sp.det(sp.Matrix([A, d]).jacobian([x, y])))
J_Bd = sp.factor(sp.det(sp.Matrix([B, d]).jacobian([x, y])))
J_AB = sp.factor(sp.det(sp.Matrix([A, B]).jacobian([x, y])))

fx = sp.diff(f, x)
fy = sp.diff(f, y)
require(sp.expand(J_Ad - 2 * (fx + x * fy)) == 0, "(A,d) Jacobian")
require(sp.expand(J_Bd + 2 * (y * fx + f * fy)) == 0, "(B,d) Jacobian")
require(
    sp.expand(J_AB + 4 * (x * fx + x**2 * fy + f - x * y)) == 0,
    "(A,B) Jacobian",
)

# The triangular coordinate t conjugates partial_x+x partial_y to partial_x.
t = sp.symbols("t")
H = t**5 - 3 * t**2 + 2 * t + 7
t_xy = y - x**2 / 2
f_Ad = sp.expand(kappa * x / 2 + H.subs(t, t_xy))
A_Ad = A
d_Ad = f_Ad
require(
    sp.factor(sp.det(sp.Matrix([A_Ad, d_Ad]).jacobian([x, y])) - kappa) == 0,
    "nonlinear (A,d) control",
)
t_inv = -sp.symbols("AA") / 2
AA, dd = sp.symbols("AA dd")
t_inv = -AA / 2
x_inv = sp.factor(2 * (dd - H.subs(t, t_inv)) / kappa)
y_inv = sp.factor(t_inv + x_inv**2 / 2)
require(sp.factor((x_inv**2 - 2 * y_inv) - AA) == 0, "(A,d) inverse A")
require(
    sp.factor(
        (kappa * x_inv / 2 + H.subs(t, y_inv - x_inv**2 / 2)) - dd
    )
    == 0,
    "(A,d) inverse d",
)

# The (A,B) equation is diagonal in C[x,t] under x*d/dx+1.
f_AB = sp.expand(x * y / 2 - x**3 / 8 - kappa / 4)
B_AB = sp.expand(y**2 - 2 * x * f_AB)
require(
    sp.factor(sp.det(sp.Matrix([A, B_AB]).jacobian([x, y])) - kappa) == 0,
    "(A,B) positive family",
)
require(
    sp.expand(B_AB.subs(y, t + x**2 / 2) - (t**2 + kappa * x / 2)) == 0,
    "(A,B) triangular form",
)
x_AB_inv = sp.factor(2 * (dd - t_inv**2) / kappa)
y_AB_inv = sp.factor(t_inv + x_AB_inv**2 / 2)
require(sp.factor((x_AB_inv**2 - 2 * y_AB_inv) - AA) == 0, "(A,B) inverse A")
require(
    sp.factor(
        (y_AB_inv**2 - 2 * x_AB_inv * f_AB.subs({x: x_AB_inv, y: y_AB_inv}))
        - dd
    )
    == 0,
    "(A,B) inverse B",
)

# Finite coefficient control: on degree <=6, d/dx has exactly the predicted
# affine solution space, while x*d/dx+1 has no kernel.
N = 6
coeffs = {
    (i, j): sp.symbols(f"c_{i}_{j}")
    for i in range(N + 1)
    for j in range(N + 1 - i)
}
generic_xt = sum(coeffs[i, j] * x**i * t**j for i, j in coeffs)
D_generic = sp.Poly(sp.diff(generic_xt, x), x, t)
free_Ad = [(0, j) for j in range(N + 1)]
forced_Ad = [(i, j) for i, j in coeffs if i > 0 and (i, j) != (1, 0)]
require(len(free_Ad) == 7 and len(forced_Ad) == 20, "degree-six (A,d) dimensions")
require(D_generic.coeff_monomial(1) == coeffs[1, 0], "degree-six slice coefficient")
for i, j in forced_Ad:
    require(
        D_generic.coeff_monomial(x ** (i - 1) * t**j) == i * coeffs[i, j],
        "degree-six derivative diagonal",
    )

eigenvalues = sorted({i + 1 for i, j in coeffs})
require(eigenvalues == list(range(1, N + 2)), "Euler-plus-one spectrum")

# If J_(B,d)=kappa, E=y*d_x+f*d_y has E(f)=-kappa/2.  Its
# action on x,y is nilpotent of lengths four and three.  The coefficient of
# T in the Jacobian of exp(T E) is exactly div(E)=f_y.
cc = sp.symbols("cc")
X_flow = x + T * y + T**2 * f / 2 + T**3 * cc / 6
Y_flow = y + T * f + T**2 * cc / 2
flow_jac = sp.expand(sp.det(sp.Matrix([X_flow, Y_flow]).jacobian([x, y])))
require(sp.expand(flow_jac).coeff(T, 1) == fy, "flow divergence coefficient")

# The THM-2696 constant-different graph is outside both positive families.
f_J0 = x * y - c0
controls = []
for pair in ((A, f_J0), (y**2 - 2 * x * f_J0, f_J0), (A, y**2 - 2 * x * f_J0)):
    controls.append(sp.factor(sp.det(sp.Matrix(pair).jacobian([x, y]))))
require(all(sp.Poly(value, x, y).total_degree() > 0 for value in controls), "J0 hostile")

print("THM-2702 polynomial graph coordinate-pair Keller classification")
print("jacobians=(A,d):2(f_x+x f_y);(B,d):-2(y f_x+f f_y);(A,B):-4(x f_x+x^2 f_y+f-xy)")
print("Ad_family=f=(kappa/2)x+H(y-x^2/2);triangular=(A,d)=(-2t,(kappa/2)x+H(t))")
print("Bd_family=EMPTY_for_kappa_nonzero;LND_ladders=x:4,y:3;divergence=f_y")
print("AB_family=f=xy/2-x^3/8-kappa/4;triangular=(A,B)=(-2t,t^2+kappa*x/2)")
print(f"degree_{N}_control=Ad_free_H:{len(free_Ad)}:Ad_forced:{len(forced_Ad)}:Euler_spectrum:{eigenvalues}")
print(f"constant_J0_graph_degrees={[sp.Poly(v, x, y).total_degree() for v in controls]}")
print("ALL CHECKS PASSED")
