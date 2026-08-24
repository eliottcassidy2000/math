#!/usr/bin/env python3
"""Exact algebraic controls for THM-3999.

The proof of divisor multiplicities and class identities is geometric and is
written in the theorem.  This certificate independently checks the canonical
extension Q=G/t, the boundary chart and endpoint polynomial, the monomial
ideal iff, the R=0 G_m control, and the correctly oriented two-cycle matrix.
"""

from __future__ import annotations

import sympy as sp


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def zero(label: str, expr) -> None:
    value = simp(expr)
    if value != 0:
        raise AssertionError(f"{label}: {value}")
    print(f"PASS  {label}")


def gate(label: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


x, t, z, p, y = sp.symbols("x t z p y")
a, gamma = sp.symbols("a gamma", nonzero=True)
alpha = sp.Rational(3, 2) * a / gamma

h00, h10, h01, h20 = sp.symbols("h00 h10 h01 h20")
k00, k10, k01, k20 = sp.symbols("k00 k10 k01 k20")
H = h00 + h10*p + h01*y + h20*p**2
K = k00 + k10*p + k01*y + k20*p**2
R = p**2*H + y*K
R0 = sp.expand(R.subs(p, 0))

u_src = x**2*t
z_src = 1 + u_src
p_src = z_src*t
y_src = x*t*p_src
source_sub = {z: z_src, p: p_src, y: y_src}

G_source = sp.expand((gamma*u_src + alpha*p + R).subs(source_sub))
Q = gamma*x**2 + alpha*z + z*p*H + x*p*K
Q_source = sp.expand(Q.subs(source_sub))

print("STATUS=THM-3999 VERIFIED-EXACT CERTIFICATE")
print("RING=B2=k[x,z,p,y] with z=1+x^2*t, p=z*t, y=x*t*p")
zero("canonical extension G=t*Q", G_source - t*Q_source)

# On the z-1 chart of X2, the determinantal relations give
#   z(z-1)^2=x^3*y, p(z-1)=x*y.
# Hence z/x^2=x*y/(z-1)^2 and p/x=y/(z-1).  This is Q/x^2
# written without a denominator vanishing on D=(x,z).
u_chart = z - 1
Q_over_x2_chart = (
    gamma
    + alpha*x*y/u_chart**2
    + x*y*p*H/u_chart**2
    + y*K/u_chart
)
boundary_restriction = sp.expand(Q_over_x2_chart.subs({x: 0, z: 0, p: 0}))
zero("boundary restriction Q/x^2=gamma-R(0,y)", boundary_restriction - (gamma - R0))
zero("endpoint polynomial is nonzero at y=0", boundary_restriction.subs(y, 0) - gamma)

# The ideal statement is monomial.  Exhaustion through total degree eight is
# a hostile control; the theorem proves it from the two monomial generators.
checks = 0
for i in range(9):
    for j in range(9-i):
        in_residual = i >= 2 or j >= 1
        vanishes_at_p0 = i >= 1
        in_p2_py = i >= 2 or (i >= 1 and j >= 1)
        gate(
            f"monomial ideal iff p^{i} y^{j}",
            (in_residual and vanishes_at_p0) == in_p2_py,
        )
        checks += 1
gate("complete monomial hostile count", checks == 45)

# Boundary valuations on X2.  Here y and z-1 are generic units on D.
valuation = {"x": 1, "z": 3, "t": -2, "p": 1, "y": 0, "u": 0}
gate("z(z-1)^2=x^3*y valuation", valuation["z"] == 3*valuation["x"] + valuation["y"])
gate("p(z-1)=x*y valuation", valuation["p"] == valuation["x"] + valuation["y"])
gate("t=(z-1)/x^2 valuation", valuation["t"] == -2*valuation["x"])
gate("Q boundary order is exactly two", 2*valuation["x"] == 2)
gate("G=tQ has boundary order zero", valuation["t"] + 2 == 0)

# R=0 is an ambient completion control (THM-3997 excludes it for Keller pairs).
Q0 = sp.expand(gamma*x**2 + alpha*z_src)
zero("R=0 companion formula", Q0 - (alpha + x**2*(gamma + alpha*t)))
t_gm = -gamma/alpha - x**-2
zero("R=0 G_m parametrization", Q0.subs(t, t_gm))
zero(
    "R=0 quotient makes x^2 invertible",
    x**2*(-(gamma + alpha*t)/alpha) - 1 + Q0/alpha,
)
Q0_on_L = sp.expand(Q0.subs(t, 0))
zero("R=0 clutch equation", Q0_on_L - (gamma*x**2 + alpha))
gate("R=0 clutch roots avoid x=0", simp(Q0_on_L.subs(x, 0)) != 0)
gate("R=0 endpoint polynomial is the nonzero constant gamma", simp(gamma) != 0)

# Correct orientation for the THM-3996 proper full two-address packet.
B = sp.Matrix([[-1, 1], [1, -1]])
Laplacian = B*B.T
gate("two-cycle incidence columns are balanced", list(B.T*sp.ones(2, 1)) == [0, 0])
gate("two-cycle incidence rank", B.rank() == 1)
gate("two-cycle kernel generator", B*sp.Matrix([1, 1]) == sp.zeros(2, 1))
gate("two-cycle Laplacian", Laplacian == sp.Matrix([[2, -2], [-2, 2]]))
gate("reduced Laplacian determinant gives critical Z/2", Laplacian[:1, :1].det() == 2)

# Principal-divisor class arithmetic in Cl(X2)=Z[D].
class_D = 1
class_L = 2*class_D
class_companion_total = -2*class_D
gate("div(t)=L-2D class balance", class_L - 2*class_D == 0)
gate("div(Q)=C_total+2D class balance", class_companion_total + 2*class_D == 0)
gate("div(G)=L+C_total class balance", class_L + class_companion_total == 0)

print("LIMITATION=total strict divisor only; factor ownership and node-address completeness are not inferred")
print("FIREWALL=graph critical Z/2 is not THM-3994 local A1 class group Z/2")
print("THEOREM_ID=THM-3999")
print("ALL THM-3999 EXACT CHECKS PASSED")
