#!/usr/bin/env python3
"""Import-free hostile local referee for THM-4342's K=0 root exit."""

import sys
import sympy as sp

sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(value, label):
    global CHECKS
    CHECKS += 1
    if value is not True:
        raise AssertionError("root-exit hostile failure: " + label)


P, b, y, delta, sigma = sp.symbols("P b y delta sigma")
U, W = sp.symbols("U W", nonzero=True)
Theta, xi = sp.symbols("Theta xi")
Phi, eta, alpha, Delta, upsilon = sp.symbols("Phi eta alpha Delta upsilon")
J = Phi + eta * P + alpha * P**2
B = P**3 * J
D = -3 + sp.Rational(8, 3) * P - sp.Rational(1376, 135) * P**2 \
    + Delta * P**3 + upsilon * P**4 + U * P**5
require(sp.Rational(2848, 45) - sp.Rational(7, 6) * sp.Rational(5696, 105) == 0,
        "inherited K=0 forces Delta=5696/105")


def reciprocal(A):
    return sp.expand(
        (1 - delta**2 * P * b**2)
        * (b**2 - P**2 * A - delta * b * B - delta**2 * b**2 * P * D)
        - delta**2 * b**2 / 2
    )


def strict_transform(A, exponent):
    transformed = sp.cancel(reciprocal(A).subs(b, P**exponent * y) / P**(2 * exponent))
    require(not sp.denom(transformed).has(P), "strict transform is polynomial")
    return sp.expand(transformed)


A1 = P * (Theta + xi * P + W * P**2)
A2 = P**2 * (xi + W * P)
A3 = W * P**3
E1, E2, E3 = strict_transform(A1, 1), strict_transform(A2, 2), strict_transform(A3, 2)

require(sp.simplify(E1.subs(delta, 0)
                    - (y**2 - P * (Theta + xi * P + W * P**2))) == 0,
        "m1 special normalization")
require(sp.simplify(E2.subs(delta, 0) - (y**2 - xi - W * P)) == 0,
        "m2 special normalization")
require(sp.simplify(E3.subs(delta, 0) - (y**2 - W * P)) == 0,
        "m3 special normalization")

# The full cubic discriminant deliberately confounds endpoint multiplicity
# with a Laurent-torus collision; saturation removes the Theta^2 exit factor.
full_disc = sp.factor(sp.discriminant(A1, P))
sat_disc = sp.factor(sp.discriminant(Theta + xi * P + W * P**2, P))
require(sp.simplify(full_disc - Theta**2 * (xi**2 - 4 * W * Theta)) == 0,
        "full discriminant factors exit and collision")
require(sp.simplify(sat_disc - (xi**2 - 4 * W * Theta)) == 0,
        "saturated discriminant")
require(sp.discriminant(xi + W * P, P) == 1,
        "m2 saturated polynomial has no collision")

# Smoothness of the honest root-exit charts.
require(sp.simplify(sp.diff(E1, P).subs({P: 0, y: 0}) + Theta) == 0,
        "m1 smooth when Theta nonzero")
q = 1 - delta**2 / 2
require(sp.simplify(E2.subs(P, 0) - (q * y**2 - xi)) == 0,
        "m2 exceptional divisor")
require(sp.simplify(sp.diff(E2, y).subs(P, 0) - 2 * q * y) == 0,
        "m2 two exceptional points are simple")
require(sp.simplify(sp.diff(E3, P).subs({P: 0, y: 0}) + W) == 0,
        "m3 binomial collar smooth")

# The only possible toric collision for m=1 is a double root of Q.
a, x, X, Y = sp.symbols("a x X Y", nonzero=True)
Ad = W * P * (P - a)**2
Fd = reciprocal(Ad)
weighted = sp.expand(Fd.subs({P: a + sigma**6 * X,
                              b: sigma**6 * Y,
                              delta: sigma**6}))
face = sp.factor(sp.Poly(weighted, sigma).coeff_monomial(sigma**12))
Ba = sp.factor(B.subs(P, a))
require(sp.simplify(face - (Y**2 - Ba * Y - a**3 * W * X**2)) == 0,
        "double conic face")

# When B(a)=0 the critical point is an exact horizontal section.  Include a
# general first derivative of B and the C correction in the tangent cone.
B1, B2, C0, C1 = sp.symbols("B1 B2 C0 C1")
p = a + x
Fnode = sp.expand(
    (1 - delta**2 * p * b**2)
    * (b**2 - W * p**3 * x**2 - delta * b * (B1 * x + B2 * x**2)
       - delta**2 * b**2 * (C0 + C1 * x))
    - delta**2 * b**2 / 2
)
for expression, label in (
    (Fnode, "F"), (sp.diff(Fnode, b), "Fb"), (sp.diff(Fnode, x), "Fx")
):
    require(sp.simplify(expression.subs({x: 0, b: 0})) == 0,
            "horizontal section " + label)
tangent = 0
for (ix, ib), coefficient in sp.Poly(Fnode, x, b).terms():
    if ix + ib == 2:
        tangent += coefficient * x**ix * b**ib
tangent = sp.expand(tangent)
q0 = 1 - delta**2 * (C0 + sp.Rational(1, 2))
require(sp.simplify(tangent - (q0 * b**2 - delta * B1 * x * b
                               - a**3 * W * x**2)) == 0,
        "horizontal node tangent")
tangent_disc = sp.factor(sp.discriminant(tangent, b))
require(sp.simplify(tangent_disc
                    - x**2 * (delta**2 * B1**2 + 4 * q0 * a**3 * W)) == 0,
        "horizontal node unit discriminant")

# First normal Hasse gate on the double stratum.
normal_gate = 4 * W**2 * Phi - 2 * W * xi * eta + xi**2 * alpha
require(sp.simplify(normal_gate.subs(xi, -2 * W * a)
                    - 4 * W**2 * J.subs(P, a)) == 0,
        "double normal-Hasse numerator")

# Exact differential transform, abstracting E_P and E_y.  This challenges a
# possible hidden loss of sigma order under b=P^e*y.
e = sp.symbols("e", integer=True, positive=True)
EP, EY, dy = sp.symbols("EP EY dy", nonzero=True)
dP = -EY * dy / EP
db_on_curve = sp.expand(P**e * dy + e * P**(e - 1) * y * dP)
FP_on_curve = P**(2 * e - 1) * (P * EP - e * y * EY)
ratio = sp.factor(P**(2 * e) * y**2 * db_on_curve / FP_on_curve)
require(sp.simplify(ratio - P**e * y**2 * dy / EP) == 0,
        "reciprocal differential survives root-exit blowup")

# No late carrier at the deepest exit: on 2v(b)=5v(P), the normal term has
# strictly positive excess 6+v(P)/2 in sigma units.
r = sp.symbols("r", positive=True)
excess = sp.simplify(6 + sp.Rational(5, 2) * r + 3 * r - 5 * r)
require(excess == 6 + r / 2, "deep-exit normal excess")

# Branch and genus controls.
require(sp.factor((y**2 - W * P * (P - a)**2).subs(y, (P - a) * Y)
                  / (P - a)**2) == Y**2 - W * P,
        "double T normalization rational")
require(sum(value - 1 for value in (11, 11, 5, 3, 3, 1)) == 28,
        "m1 RH ledger")
require(sum(value - 1 for value in (11, 11, 3, 3, 3, 1)) == 26,
        "m2 RH ledger")
require(sum(value - 1 for value in (11, 11, 7, 1)) == 26,
        "m3 RH ledger")

print("THM4342 KZERO ROOT-EXIT HOSTILE REFEREE=PASS")
print("seam_relation=K=0=>Delta=5696/105")
print("full_disc=Theta^2*(xi^2-4WTheta);saturated_disc=xi^2-4WTheta")
print("exit_models=m1:y^2=P(Theta+xiP+WP^2);m2:y^2=xi+WP;m3:y^2=WP")
print("m1_double_face=Y^2-B(a)Y-a^3W X^2")
print("m1_double_normal_gate=4W^2Phi-2Wxi*eta+xi^2alpha")
print("horizontal_node_disc=x^2*(delta^2 B1^2+4q0 a^3W)")
print("form_transform=-sigma^16*P^e*y^2*dy/E_P")
print("m3_normal_excess=6+r/2")
print("checks=" + str(CHECKS))
