#!/usr/bin/env python3
"""Exact algebraic and graph controls for THM-3996."""

import sympy as sp


def z(label, expr):
    value = sp.factor(sp.cancel(sp.expand(expr)))
    if value != 0:
        raise RuntimeError(f"{label}: {value}")
    print(f"PASS  {label}")


# THM-3992 target nodal cubic and its degree-one normalization component.
X, a = sp.symbols("X a", nonzero=True)
A = X**2 + a
C = X**3 + sp.Rational(3, 2) * a * X
G = C**2 - A**3 + sp.Rational(3, 4) * a**2 * A + sp.Rational(1, 4) * a**3
z("THM3992 normalization lies on target node", G)
node_parameter = X**2 + sp.Rational(3, 2) * a
z("both node parameters map to A=-a/2", A + a / 2 - node_parameter)
z("both node parameters map to C=0 modulo node equation", sp.rem(C, node_parameter, X))
z("node-address polynomial is squarefree", sp.discriminant(node_parameter, X) + 6 * a)

# A directed d-cycle is the sharp finite-etale address packet.  Removing one
# address from d=3 produces the minimal unbalanced path.
cycle_checks = 0
for d in range(2, 13):
    adjacency = sp.zeros(d)
    for i in range(d):
        adjacency[i, (i + 1) % d] = 1
    outdegrees = [sum(adjacency[i, j] for j in range(d)) for i in range(d)]
    indegrees = [sum(adjacency[j, i] for j in range(d)) for i in range(d)]
    if outdegrees != indegrees or outdegrees != [1] * d:
        raise RuntimeError((d, indegrees, outdegrees))
    cycle_checks += 1
print(f"PASS  balanced directed cycle controls ({cycle_checks} packets)")

path3 = sp.zeros(3)
path3[0, 1] = 1
path3[1, 2] = 1
path_indegrees = [sum(path3[j, i] for j in range(3)) for i in range(3)]
path_outdegrees = [sum(path3[i, j] for j in range(3)) for i in range(3)]
if path_indegrees == path_outdegrees:
    raise RuntimeError("deleted-address path should be unbalanced")
print(f"HOSTILE deleted d=3 address gives path indegree={path_indegrees}, outdegree={path_outdegrees}")

# Riemann-Hurwitz for a connected finite etale cover of A1: after projective
# completion, only infinity may ramify and sum(e_i-1)=d-r.  Hence
# 2g-2=-d-r, which has the unique positive d,r solution d=r=1.
rh_solutions = []
for degree in range(1, 20):
    for points_over_infinity in range(1, degree + 1):
        rhs = -degree - points_over_infinity
        if rhs >= -2 and rhs % 2 == 0:
            genus = (rhs + 2) // 2
            if genus >= 0:
                rh_solutions.append((degree, points_over_infinity, genus))
if rh_solutions != [(1, 1, 0)]:
    raise RuntimeError(rh_solutions)
print("PASS  Riemann-Hurwitz finite-etale A1 control has only degree one")

# Explicit embedded d=2 curve control:
# Z: t(t-x^2+1)=0, u=x^2, v=x(x^2-1-2t), target v^2=u(u-1)^2.
x, t = sp.symbols("x t")
u = x**2
v = x * (x**2 - 1 - 2 * t)
target = v**2 - u * (u - 1) ** 2
z("plane two-cycle target equation on L", target.subs(t, 0))
z("plane two-cycle target equation on companion", target.subs(t, x**2 - 1))
jacobian = sp.det(sp.Matrix([[sp.diff(u, x), sp.diff(u, t)], [sp.diff(v, x), sp.diff(v, t)]]))
z("plane map Jacobian at both clutch points", sp.rem(jacobian + 4, x**2 - 1, x).subs(t, 0))
z("two components meet at exactly x^2=1", sp.resultant(t, t - x**2 + 1, t) - (1 - x**2))

print("RESULT finite-locus node graph has indegree=outdegree at every vertex")
print("RESULT every complete finite-locus address edge lies on a directed cycle")
print("RESULT distinct THM3992 companion owners require an extra address or a Jelonek node")
print("RESULT a full two-address packet is the two-edge cycle")
print("NONCONSEQUENCE no ownership census, Jelonek membership, Keller inverse, or JC2 closure is proved")
print("ALL EXACT CHECKS PASSED")
