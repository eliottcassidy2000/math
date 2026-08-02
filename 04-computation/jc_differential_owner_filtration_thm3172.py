#!/usr/bin/env python3
"""Exact controls for THM-3172's differential-owner filtration.

The proofs are algebraic; these checks exercise the implicit-derivative
recurrence and the hostile and positive controls.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def jacobian(a: sp.Expr, b: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.diff(a, x) * sp.diff(b, y) - sp.diff(a, y) * sp.diff(b, x))


def transverse_numerators(
    f: sp.Expr,
    T: sp.Symbol,
    v: sp.Symbol,
    depth: int,
) -> list[sp.Expr]:
    """Return A_n with root derivative d_v^n xi=A_n/f_T^(2n-1)."""
    f_T = sp.diff(f, T)
    f_v = sp.diff(f, v)

    def lam(h: sp.Expr) -> sp.Expr:
        return sp.expand(f_T * sp.diff(h, v) - f_v * sp.diff(h, T))

    values = [-f_v]
    for n in range(1, depth):
        A_n = values[-1]
        values.append(sp.expand(f_T * lam(A_n) - (2 * n - 1) * A_n * lam(f_T)))
    return values


u, v, T = sp.symbols("u v T")

print("DIFFERENTIAL-OWNER / TRANSVERSE-COVARIANT EXACT AUDIT")

# Target-shear covariance of the base derivation span.
U, V = sp.symbols("U V")
H = U**2 + U
g = u**2 * v + v**3 + u * v**2
g_shear = sp.expand(g.subs({u: U, v: V - H}))
old_du_minus = sp.expand((sp.diff(g, u) - (2 * u + 1) * sp.diff(g, v)).subs({u: U, v: V - H}))
old_dv = sp.expand(sp.diff(g, v).subs({u: U, v: V - H}))
require(sp.expand(sp.diff(g_shear, U) - old_du_minus) == 0, "u-shear chain rule")
require(sp.expand(sp.diff(g_shear, V) - old_dv) == 0, "v-shear chain rule")
print("[PASS] polynomial target shear preserves the generated base-derivation module")

# Family of genuine nonlinear polynomial automorphisms.  Here
# P=x+(y+x^2)^m, Q=y+x^2 and the inverse x=u-v^m is the transverse test.
x, y = sp.symbols("x y")
automorphism_cases = 0
recurrence_checks = 0
for m in range(2, 9):
    w = y + x**2
    P = x + w**m
    Q = w
    require(jacobian(P, Q, x, y) == 1, ("automorphism Jacobian", m))
    root = u - v**m
    f = T - root
    values = transverse_numerators(f, T, v, m + 2)
    f_T = sp.diff(f, T)
    for n, A_n in enumerate(values, start=1):
        implicit = sp.cancel(A_n.subs(T, root) / f_T.subs(T, root) ** (2 * n - 1))
        direct = sp.diff(root, v, n)
        require(sp.expand(implicit - direct) == 0, ("automorphism recurrence", m, n))
        recurrence_checks += 1
    require(values[m - 1] != 0 and values[m] == 0, ("sharp termination", m))
    automorphism_cases += 1
print(
    f"[PASS] {automorphism_cases} nonlinear automorphisms: "
    f"{recurrence_checks} implicit derivatives, sharp A_m!=0 and A_(m+1)=0"
)

# One explicit target shear of the nonlinear m=4 automorphism retains the
# unit Jacobian and hence a unimodular gradient for the sheared mate.
m = 4
w = y + x**2
P = x + w**m
Q = w
H_source = P**2 + 3 * P
S = sp.expand(Q + H_source)
require(jacobian(P, S, x, y) == 1, "sheared automorphism")
print("[PASS] nonlinear automorphism plus quadratic target shear keeps Jac(P,S)=1")

# Cubic punctured branch of the THM-3167 hostile.  The map is Keller on G_m,
# but the first target derivative adjoins X^{-1}; this cannot happen in an
# A^2 polynomial owner, whose only units are constants.
X, Y, r = sp.symbols("X Y r", nonzero=True)
t = sp.symbols("t")
u_cubic = -r * Y / (3 * X**2)
v_cubic = X**3
require(sp.simplify(jacobian(u_cubic, v_cubic, X, Y) - r) == 0, "cubic Laurent Jacobian")
f_split = (T - u) * (T**3 - t)
dx_dv_from_f = -sp.diff(f_split, t) / sp.diff(f_split, T)
dx_dv = sp.cancel(dx_dv_from_f.subs({T: X, t: X**3}))
require(sp.cancel(dx_dv - 1 / (3 * X**2)) == 0, "inverse-spectral cubic derivative")
inverse_x = sp.cancel(3 * X * dx_dv)
require(sp.cancel(X * inverse_x) == 1, "first derivative exposes inverse X")
H_cubic = u_cubic**2 + u_cubic
require(
    sp.simplify(jacobian(u_cubic, v_cubic + H_cubic, X, Y) - r) == 0,
    "cubic target shear",
)
print("[PASS] THM-3167 cubic branch: first differential hull contains X^{-1}")
print("[PASS] the same cubic unit obstruction survives a nonlinear target shear")

# The cubic root derivatives do not terminate.  Check six levels by an
# independent branch derivation (d/dt)=(3X^2)^(-1)d/dX.
f_cubic = T**3 - t
cubic_values = transverse_numerators(f_cubic, T, t, 6)
direct = X
for n, A_n in enumerate(cubic_values, start=1):
    direct = sp.cancel(sp.diff(direct, X) / (3 * X**2))
    implicit = sp.cancel(A_n.subs(T, X) / (3 * X**2) ** (2 * n - 1))
    require(sp.cancel(direct - implicit) == 0, ("cubic recurrence", n))
    require(implicit != 0, ("cubic nontermination", n))
print("[PASS] cubic punctured branch: six independent nonterminating transverse derivatives")

# Falling-factorial closed form for every power branch t=X^m.  This is the
# exact sequence-side contrast with the terminating automorphism tail above:
# positive-integer exponent m terminates, reciprocal exponent 1/m does not.
power_branch_checks = 0
for power in range(2, 8):
    f_power = T**power - t
    values = transverse_numerators(f_power, T, t, 7)
    for n, A_n in enumerate(values, start=1):
        coefficient = sp.prod(sp.Rational(1, power) - j for j in range(n))
        expected = coefficient * X ** (1 - power * n)
        implicit = sp.cancel(
            A_n.subs(T, X) / (power * X ** (power - 1)) ** (2 * n - 1)
        )
        require(sp.cancel(implicit - expected) == 0, ("power branch closed form", power, n))
        require(implicit != 0, ("power branch nontermination", power, n))
        power_branch_checks += 1
print(
    f"[PASS] {power_branch_checks} power-branch values: "
    "partial_t^n(X)=(1/m)_(falling n) X^(1-mn), never zero for m>1"
)

# Polynomial Gate-II hostile u=X^4, v=XY.  Its x-coordinate is constant in
# the v direction (so the transverse A_n test alone is blind), while the u
# derivative exposes X^{-1} and the physical Jacobian is nonconstant.
u_gate2 = X**4
v_gate2 = X * Y
require(sp.expand(jacobian(u_gate2, v_gate2, X, Y) - 4 * X**4) == 0, "Gate-II Jacobian")
dx_du = 1 / (4 * X**3)
require(sp.cancel(X * (4 * X**2 * dx_du)) == 1, "Gate-II differential unit")
gate2_f = T**4 - u
require(transverse_numerators(gate2_f, T, v, 3) == [0, 0, 0], "Gate-II transverse blindness")
print("[PASS] polynomial Gate-II hostile: v-tail is blind but the full differential hull gains X^{-1}")

# THM-2241's sharp ownership warning: transverse termination can occur on a
# Laurent slice even when no polynomial A^2 mate/properness package exists.
h = x * y**2
s = -1 / y
require(sp.cancel(jacobian(h, s, x, y)) == 1, "Laurent slice Jacobian")


def D_h(a: sp.Expr) -> sp.Expr:
    return sp.expand(jacobian(h, a, x, y))


d1 = D_h(x)
d2 = D_h(d1)
d3 = D_h(d2)
require(d1 != 0 and d2 != 0 and d3 == 0, "Laurent finite-tail near miss")
print("[PASS] Laurent near miss h=xy^2, s=-1/y: D_h^3(x)=0 without an A^2 owner")

# A compositional target representative forces a derivative unit.  In the
# cubic branch Phi(X)=X^3, Phi'(X)=3X^2 is invertible only because X is a
# Laurent unit.  This is the exact one-jet mechanism of the owner gate.
phi_prime = 3 * X**2
require(sp.cancel(phi_prime * dx_dv) == 1, "composition derivative unit")
print("[PASS] compositional gate: Phi'(X) * partial_v(X)=1 for Phi(X)=X^3")

print("SUMMARY: 0 failed controls; all arithmetic is exact over QQ(symbols)")
