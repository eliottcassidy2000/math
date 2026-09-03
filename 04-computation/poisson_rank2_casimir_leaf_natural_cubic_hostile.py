#!/usr/bin/env python3
"""Exact hostile probe for the natural Casimir-leaf cubic coordinate.

On the nonzero core leaf R=a, Long's localized variables give a natural
depressed cubic for w=1+xy over the target (T,S).  This script proves that
identity and exhibits a rational etale source point at which w is a double
root.  It therefore tests only this primitive coordinate, not every possible
cubic presentation of the leaf.
"""

from __future__ import annotations

import sympy as sp


def check(condition: bool, label: str) -> None:
    if not bool(condition):
        raise RuntimeError(f"FAIL: {label}")


x, y, beta = sp.symbols("x y beta")
w = 1 + x * y
Pi = sp.expand(w**3 * beta + y**2 * w * (4 + 3 * x * y))
Sigma = sp.expand(y + 3 * x * w**2 * beta + 3 * x * y**2 * (4 + 3 * x * y))
R = sp.expand(2 * x - 3 * x**2 * y - x**3 * beta)
T = sp.expand(-Pi / 2)
S = Sigma

# Reconstruct the localized formulas in independent variables (x,w,a).
X, W, a = sp.symbols("X W a")
alpha = a / X
b = 2 + 4 * W - 3 * alpha * W**2
c = sp.expand((alpha * W**3 - W**2 - W) / 2)
S_local = sp.cancel(b / X)
T_local = sp.cancel(c / X**2)

K = sp.expand(a * S_local**3 - S_local**2 + 36 * a * S_local * T_local
              + 108 * a**2 * T_local**2 - 32 * T_local)
H = sp.expand(S_local**2 + 24 * T_local)
P_on_source = sp.factor(K * W**3 + H * W + 8 * T_local)
check(P_on_source == 0, "natural depressed cubic identity")

# The same identity reconstructed directly from the polynomial core.
A, B, aa, Z = sp.symbols("A B aa Z")
K_target = aa * B**3 - B**2 + 36 * aa * B * A + 108 * aa**2 * A**2 - 32 * A
P = sp.expand(K_target * Z**3 + (B**2 + 24 * A) * Z + 8 * A)
P_derivative = sp.diff(P, Z)
direct_residual = sp.factor(P.subs({aa: R, A: T, B: S, Z: w}))
check(direct_residual == 0, "direct core substitution into cubic")

# Rational hostile: a valid point of the nonzero R=4 leaf is a double root
# of the natural w-cubic, even though the core and leaf maps remain etale.
point = {x: sp.Integer(1), y: -sp.Rational(4, 3), beta: sp.Integer(2)}
target = tuple(sp.factor(f.subs(point)) for f in (R, T, S))
root = sp.factor(w.subs(point))
check(target == (4, sp.Rational(1, 27), -sp.Rational(2, 3)),
      "rational hostile target")
check(root == -sp.Rational(1, 3), "rational hostile root")

P_hostile = sp.factor(P.subs({aa: target[0], A: target[1], B: target[2]}))
expected_hostile = -sp.Rational(4, 27) * (3 * Z - 2) * (3 * Z + 1)**2
check(sp.expand(P_hostile - expected_hostile) == 0,
      "hostile cubic factorization")
check(P_hostile.subs(Z, root) == 0, "hostile root")
check(sp.diff(P_hostile, Z).subs(Z, root) == 0, "hostile double root")
check(sp.diff(P_hostile, Z, 2).subs(Z, root) != 0, "hostile multiplicity exactly two")

core_jacobian = sp.factor(sp.Matrix([R, T, S]).jacobian([x, y, beta]).det())
check(core_jacobian == 1, "normalized three-dimensional core determinant")

# On R=a!=0, put s=1/x and v=xy.  The leaf chart is G_m x A^1 and
# det d(T,S)/d(s,v)=s^2, so the hostile point s=1 is etale on the leaf.
s, v = sp.symbols("s v")
leaf_substitution = {
    x: 1 / s,
    y: v * s,
    beta: s**2 * (2 - 3 * v - aa * s),
}
T_leaf = sp.factor(T.subs(leaf_substitution))
S_leaf = sp.factor(S.subs(leaf_substitution))
leaf_jacobian = sp.factor(sp.Matrix([T_leaf, S_leaf]).jacobian([s, v]).det())
check(leaf_jacobian == s**2, "leaf target Jacobian")
check(leaf_jacobian.subs({s: 1, aa: 4, v: -sp.Rational(4, 3)}) == 1,
      "hostile point is etale")

print("RANK-TWO POISSON CASIMIR-LEAF NATURAL CUBIC HOSTILE")
print("P_a,T,S(W)=K W^3+(S^2+24T)W+8T")
print("K=a*S^3-S^2+36*a*S*T+108*a^2*T^2-32*T")
print("universal substitution P(R,T,S;1+xy)=0: PASS")
print("hostile source (x,y,beta)=(1,-4/3,2)")
print("hostile target (R,T,S)=(4,1/27,-2/3)")
print("hostile root W=-1/3")
print("hostile factorization:", P_hostile)
print("P_W=0 and P_WW!=0 at hostile root: PASS")
print("core determinant=1; leaf determinant at hostile=1: PASS")
print("verdict=natural W-incidence is not globally its simple-root locus")
print("scope=does_not_exclude_another_primitive_coordinate_or_normalized_incidence")
print("RESULT=PASS")
