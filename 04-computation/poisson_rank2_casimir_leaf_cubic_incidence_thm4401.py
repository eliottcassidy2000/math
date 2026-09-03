#!/usr/bin/env python3
"""Primary exact certificate for the reserved THM-4401 candidate.

This independently starts from the three-variable core in arXiv:2608.23777
and checks its Casimir-leaf parameterizations, cubic incidence, discriminant,
affine completion, zero-level split, and a complementary fixed-(T,S) slice.
All checks use exact SymPy arithmetic.  No assertion statement is used, so
``python -O`` executes the same mathematical gates.
"""

from __future__ import annotations

import sympy as sp


checks = 0


def require_zero(name: str, expression: sp.Expr) -> None:
    global checks
    checks += 1
    value = sp.factor(sp.cancel(expression))
    if value != 0:
        raise RuntimeError(f"FAIL {name}: {value}")


# ---------------------------------------------------------------------------
# 1. Reconstruct the core directly from the preprint formulas.
# ---------------------------------------------------------------------------
x, y, beta = sp.symbols("x y beta")
paper_u = x * y
w = 1 + paper_u
R_core = 2 * x - 3 * x**2 * y - x**3 * beta
S_core = y + 3 * x * w**2 * beta + 3 * x * y**2 * (4 + 3 * paper_u)
T_core = -sp.Rational(1, 2) * (
    w**3 * beta + y**2 * w * (4 + 3 * paper_u)
)
alpha = 2 - 3 * x * y - x**2 * beta

require_zero("R=x*alpha", R_core - x * alpha)
require_zero(
    "zero-level CRT Bezout",
    1 - alpha / 2 - x * (3 * y + x * beta) / 2,
)


# ---------------------------------------------------------------------------
# 2. Every nonzero Casimir leaf R=rho is G_m x A^1.
#
# New leaf coordinates (not the preprint's u=xy):
#       t=1/x,  u=t(1+xy)=t+y.
# The inverse is x=1/t, y=u-t, beta=5t^2-3tu-rho*t^3.
# ---------------------------------------------------------------------------
rho, t, u = sp.symbols("rho t u")
leaf_inverse = {
    x: 1 / t,
    y: u - t,
    beta: 5 * t**2 - 3 * t * u - rho * t**3,
}
R_leaf = sp.factor(R_core.subs(leaf_inverse))
S_leaf = sp.factor(S_core.subs(leaf_inverse))
T_leaf = sp.factor(T_core.subs(leaf_inverse))

S_expected = 2 * t + 4 * u - 3 * rho * u**2
T_expected = (rho * u**3 - u**2 - u * t) / 2
require_zero("nonzero-leaf R", R_leaf - rho)
require_zero("nonzero-leaf S", S_leaf - S_expected)
require_zero("nonzero-leaf T", T_leaf - T_expected)
require_zero("forward coordinate u=t+y", u - (t + leaf_inverse[y]))

J_leaf = sp.det(
    sp.Matrix(
        [
            [sp.diff(S_expected, t), sp.diff(S_expected, u)],
            [sp.diff(T_expected, t), sp.diff(T_expected, u)],
        ]
    )
)
require_zero("leaf Jacobian", J_leaf + t)


# ---------------------------------------------------------------------------
# 3. The natural affine completion is the marked-root cubic incidence.
# ---------------------------------------------------------------------------
U, S, T = sp.symbols("U S T")
P = rho * U**3 - 2 * U**2 + S * U + 4 * T
P_at_leaf = P.subs({U: u, S: S_expected, T: T_expected})
Pprime_at_leaf = sp.diff(P, U).subs({U: u, S: S_expected})
require_zero("marked-root equation", P_at_leaf)
require_zero("derivative coordinate P'(u)=2t", Pprime_at_leaf - 2 * t)
require_zero(
    "recover t from target and root",
    t - (S_expected - 4 * u + 3 * rho * u**2) / 2,
)

disc = sp.factor(sp.discriminant(P, U))
disc_expected = 4 * (
    S**2
    - rho * S**3
    + 32 * T
    - 36 * rho * S * T
    - 108 * rho**2 * T**2
)
require_zero("cubic discriminant", disc - disc_expected)

disc_pullback = sp.factor(disc.subs({S: S_expected, T: T_expected}))
disc_pullback_expected = 4 * t**2 * (
    9 * rho**2 * u**2 - 8 * rho * t - 12 * rho * u + 4
)
require_zero("discriminant pullback", disc_pullback - disc_pullback_expected)


# ---------------------------------------------------------------------------
# 4. For rho != 0, affine equivalence to the universal depressed cubic.
# ---------------------------------------------------------------------------
v = u - sp.Rational(2, 3) / rho
p_target = S / rho - sp.Rational(4, 3) / rho**2
q_target = (
    4 * T / rho
    + 2 * S / (3 * rho**2)
    - sp.Rational(16, 27) / rho**3
)
p_source = sp.factor(p_target.subs(S, S_expected))
q_source = sp.factor(q_target.subs({S: S_expected, T: T_expected}))
require_zero("depressed-cubic source p", p_source - (2 * t / rho - 3 * v**2))
require_zero("universal marked-root equation", q_source + v**3 + p_source * v)
require_zero("deleted universal ramification divisor", p_source + 3 * v**2 - 2 * t / rho)

triple_S = sp.Rational(4, 3) / rho
triple_T = -sp.Rational(2, 27) / rho**2
require_zero("triple-root target p=0", p_target.subs(S, triple_S))
require_zero(
    "triple-root target q=0",
    q_target.subs({S: triple_S, T: triple_T}),
)


# ---------------------------------------------------------------------------
# 5. The rho=0 fiber is A^2 disjoint union (G_m x A^1).
# The curved component is exactly THM-3554 after T=-F1/2 and S=F2.
# ---------------------------------------------------------------------------
S_zero = sp.expand(S_expected.subs(rho, 0))
T_zero = sp.expand(T_expected.subs(rho, 0))
require_zero("rho=0 Kummer identity", S_zero**2 + 32 * T_zero - 4 * t**2)

S_plane = sp.factor(S_core.subs(x, 0))
T_plane = sp.factor(T_core.subs(x, 0))
require_zero("x=0 plane S", S_plane - y)
require_zero("x=0 plane T", T_plane + beta / 2 + 2 * y**2)
J_plane = sp.det(
    sp.Matrix(
        [
            [sp.diff(S_plane, y), sp.diff(S_plane, beta)],
            [sp.diff(T_plane, y), sp.diff(T_plane, beta)],
        ]
    )
)
require_zero("x=0 plane Jacobian", J_plane + sp.Rational(1, 2))
require_zero(
    "x=0 plane inverse",
    T_plane.subs({y: S, beta: -2 * T - 4 * S**2}) - T,
)


# ---------------------------------------------------------------------------
# 6. Complementary hostile: fix (T,S)=(tau,sigma).
# For tau != 0, w is a unit and h=x/w gives
#   C_(sigma,tau) = Spec k[h,(1-sigma*h-6*tau*h^2)^(-1)].
# The full two-dimensional slice adds the free D-coordinate.
# ---------------------------------------------------------------------------
h, sigma, tau = sp.symbols("h sigma tau")
d = 1 - sigma * h - 6 * tau * h**2
y_h = sigma + 6 * tau * h
x_h = h / d
beta_h = -2 * tau * d**3 - 4 * y_h**2 * d**2 - 3 * h * y_h**3 * d
fixed_inverse = {x: x_h, y: y_h, beta: beta_h}

require_zero("fixed-(T,S) T", T_core.subs(fixed_inverse) - tau)
require_zero("fixed-(T,S) S", S_core.subs(fixed_inverse) - sigma)
R_h = 2 * h - sigma * h**2 - 4 * tau * h**3
require_zero("fixed-(T,S) R", R_core.subs(fixed_inverse) - R_h)
require_zero("fixed-(T,S) w inverse", w.subs(fixed_inverse) * d - 1)
require_zero("fixed-(T,S) etale derivative", sp.diff(R_h, h) - 2 * d)

# The rational identities producing h and d directly from the core.
h_core = x / w
require_zero("global rational y identity", y - (S_core + 6 * T_core * h_core))
require_zero(
    "global rational d identity",
    w * (1 - S_core * h_core - 6 * T_core * h_core**2) - 1,
)
require_zero(
    "global rational cubic inverse",
    R_core - (2 * h_core - S_core * h_core**2 - 4 * T_core * h_core**3),
)

# Distinguished connected punctured triple slice.
R_collision_slice = sp.factor(R_h.subs({sigma: 0, tau: sp.Rational(1, 8)}))
d_collision_slice = sp.factor(d.subs({sigma: 0, tau: sp.Rational(1, 8)}))
require_zero("collision root h=0", R_collision_slice.subs(h, 0))
require_zero("collision root h=2", R_collision_slice.subs(h, 2))
require_zero("collision root h=-2", R_collision_slice.subs(h, -2))
require_zero("collision roots avoid deleted divisor +2", d_collision_slice.subs(h, 2) + 2)
require_zero("collision roots avoid deleted divisor -2", d_collision_slice.subs(h, -2) + 2)


# A literal realization of all four canonical outputs in a planar Poisson
# algebra is impossible before any geometry: the required bracket matrix has
# rank four, whereas J*Omega_2*J^T has rank at most two.
required_brackets = sp.Matrix(
    [
        [0, 0, -1, 0],
        [0, 0, 0, -1],
        [1, 0, 0, 0],
        [0, 1, 0, 0],
    ]
)
require_zero("required bracket determinant", required_brackets.det() - 1)


print(f"PASS {checks} exact gates")
print("R=rho leaf ring (rho!=0): k[t,t^-1,u]")
print("S =", S_expected)
print("T =", T_expected)
print("Jac_(t,u)(S,T) =", sp.factor(J_leaf))
print("marked-root polynomial =", P)
print("P'(u) on incidence =", sp.factor(Pprime_at_leaf))
print("discriminant =", disc)
print("discriminant pullback =", disc_pullback)
print("rho!=0 completion: finite flat degree 3; critical divisor t=0")
print("rho=0 curved completion: finite flat degree 2; S^2+32*T=4*t^2")
print("rho=0 plane component: (y,beta)->(S,T)=(y,-beta/2-2*y^2), an automorphism")
print("fixed (T,S)=(1/8,0): R(h)=", R_collision_slice)
print("fixed-slice deleted factor =", d_collision_slice)
print("literal planar four-output Poisson bracket matrix determinant =", required_brackets.det())
