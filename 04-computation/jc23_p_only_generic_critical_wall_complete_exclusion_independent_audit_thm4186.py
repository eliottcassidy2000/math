#!/usr/bin/env python3
"""Independent normalized-coordinate audit for THM-4186.

This implementation reconstructs the (X,T) source directly.  It imports no
code or polynomial from the source-chart certificate.  It checks the exact
M=9 filtration, the normalized Hessian bridge and universal fibres, and then
recomputes three exact resultants: one generic control, one genuine projected-
discriminant hostile, and one genuine Q20(-1/6)=0 hostile.
"""

from __future__ import annotations

from hashlib import sha256

import sympy as sp


CHECKS = 0


def need(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuation(poly: sp.Expr, variable: sp.Symbol) -> int:
    return min(monomial[0] for monomial, _ in sp.Poly(poly, variable).terms())


X, T = sp.symbols("X T")
Delta, Phi, Theta, eta = sp.symbols("Delta Phi Theta eta")
K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
P = T + X**2 * T**2
Y = X * T * P
G = sp.expand(
    -X**2 * T / 2
    - 3 * P
    + sp.Rational(8, 3) * P**2
    - sp.Rational(1376, 135) * P**3
    + K * Y**2
    + Phi * P**2 * Y
    + Delta * P**4
    + Theta * P * Y**2
    + eta * P**3 * Y
)
f = sp.cancel(sp.diff(G, X) / T)
h = sp.diff(G, T)
need(sp.denom(f) == 1, "G_X/T is not polynomial")
need(sp.factor(sp.diff(G, X) - T * f) == 0, "G_X=T*f changed")

# Source completeness: the only residual-weight-nine monomials are p^3*y
# and y^3.  Therefore eta=zeta=0 is an M<=8 filtration exit, not a Y-only
# exact-M=9 wall.
admissible = {
    (i, j)
    for i in range(5)
    for j in range(4)
    if 0 < 2 * i + 3 * j <= 9 and (i, j) not in {(0, 1), (1, 1)}
}
need(admissible == {
    (1, 0), (2, 0), (3, 0), (0, 2), (2, 1),
    (4, 0), (1, 2), (3, 1), (0, 3),
}, "complete residual monomial universe changed")
need({pair for pair in admissible if 2 * pair[0] + 3 * pair[1] == 9}
     == {(3, 1), (0, 3)}, "weight-nine top row changed")

# No finite T!=0 common point can escape at X=infinity in the P chamber.
f_poly = sp.Poly(f, X)
h_poly = sp.Poly(h, X)
need((f_poly.degree(), h_poly.degree()) == (8, 9),
     "normalized X-degrees changed")
need(sp.factor(f_poly.LC() - 9 * eta * T**8) == 0,
     "f leading row changed")
need(sp.factor(h_poly.LC() - 9 * eta * T**8) == 0,
     "h leading row changed")

# Exact row-independent bridge.  At a common zero with T!=0 it identifies
# det D(f,h) with det Hess(G)/T, so Keller--Morse reduces the full source
# scheme even when its T-eliminant has repeated roots.
normalized_jacobian = sp.det(sp.Matrix((
    (sp.diff(f, X), sp.diff(f, T)),
    (sp.diff(h, X), sp.diff(h, T)),
)))
normalized_hessian = sp.det(sp.hessian(G, (X, T)))
need(sp.factor(T * normalized_jacobian
               - normalized_hessian - f * sp.diff(G, X, T)) == 0,
     "normalized Hessian bridge changed")

# The two universal fibres, including their values and Morse determinants.
need(sp.factor(f.subs(T, 0) + X) == 0, "T=0 f row changed")
need(sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0,
     "T=0 h row changed")
need(sp.rem(sp.Poly(normalized_hessian.subs(T, 0) - 6, X),
            sp.Poly(X**2 + 6, X)).is_zero,
     "T=0 universal pair ceased to be Morse")

t_universal = -sp.Rational(1, 6)
universal_modulus = sp.Poly(X**2 - 6, X)
for expression, expected, label in (
    (f, 0, "f"),
    (h, 0, "h"),
    (G, sp.Rational(1, 2), "value"),
    (normalized_hessian, -6, "Hessian"),
):
    remainder = sp.rem(
        sp.Poly(sp.expand(expression.subs(T, t_universal) - expected), X),
        universal_modulus,
    )
    need(remainder.is_zero, f"T=-1/6 universal {label} changed")

# The collapsed p=P coordinate has no further critical point: on P=0 with
# T!=0, substitution T=-1/X^2 forces X^2=6 in both critical equations.
need(sp.factor(f.subs(T, -1 / X**2) + (X**2 - 6) / X) == 0,
     "P=0 nonzero-T f row changed")
need(sp.factor(h.subs(T, -1 / X**2) + (X**2 - 6) / 2) == 0,
     "P=0 nonzero-T h row changed")


def normalized_residual(
    name: str,
    values: dict[sp.Symbol, sp.Expr],
) -> tuple[sp.Poly, sp.Expr, sp.Expr, sp.Expr]:
    specialized_G = sp.expand(G.subs(values))
    specialized_f = sp.cancel(sp.diff(specialized_G, X) / T)
    specialized_h = sp.diff(specialized_G, T)
    result = sp.resultant(specialized_f, specialized_h, X)
    need(valuation(result, T) == 56, f"{name}: T-artifact changed")
    quotient = sp.cancel(result / (T**56 * (6 * T + 1)**2))
    need(sp.denom(quotient) == 1, f"{name}: residual is not polynomial")
    residual = sp.Poly(quotient, T)
    specialized_K = sp.factor(K.subs(values))
    need(residual.degree() == 20, f"{name}: residual degree changed")
    need(sp.factor(
        residual.TC() + sp.Rational(3**15, 2**7) * values[eta]**7
    ) == 0, f"{name}: residual constant changed")
    need(sp.factor(
        residual.LC()
        - 72900 * values[eta]**5 * specialized_K**4 * values[Theta]**4
    ) == 0, f"{name}: residual leading row changed")
    need(sp.Poly(specialized_f, X).LC() == 9 * values[eta] * T**8,
         f"{name}: f infinity row changed")
    need(sp.Poly(specialized_h, X).LC() == 9 * values[eta] * T**8,
         f"{name}: h infinity row changed")
    return residual, specialized_G, specialized_f, specialized_h


generic_values = {
    Delta: sp.Integer(2),
    Phi: sp.Integer(5),
    Theta: sp.Integer(7),
    eta: sp.Integer(11),
}
generic_Q, _, _, _ = normalized_residual("generic", generic_values)
need(sp.gcd(generic_Q, generic_Q.diff()).degree() == 0,
     "generic control residual is not squarefree")
need(generic_Q.eval(t_universal) != 0,
     "generic control meets universal fibre")

# Genuine discriminant wall: two distinct Morse points X=0,1 share T=1.
collision_values = {
    Delta: sp.Rational(1271, 180),
    Phi: sp.Rational(1733, 7560),
    Theta: -sp.Rational(206281, 7560),
    eta: -sp.Rational(1733, 7560),
}
collision_Q, collision_G, collision_f, collision_h = normalized_residual(
    "collision", collision_values
)
need(sp.factor(K.subs(collision_values)) == sp.Rational(11891, 216),
     "collision K changed")
need(sp.factor(sp.gcd(collision_Q, collision_Q.diff()).monic().as_expr())
     == T - 1, "collision residual gcd changed")
collision_fibre_gcd = sp.gcd(
    sp.Poly(collision_f.subs(T, 1), X),
    sp.Poly(collision_h.subs(T, 1), X),
).monic()
need(sp.factor(collision_fibre_gcd.as_expr()) == X * (X - 1),
     "collision T=1 fibre changed")
collision_hessian = sp.det(sp.hessian(collision_G, (X, T)))
for x_value in (sp.Integer(0), sp.Integer(1)):
    need(collision_hessian.subs({X: x_value, T: 1}) != 0,
         "collision point ceased to be Morse")
    need(collision_G.subs({X: x_value, T: 1})
         not in {sp.Integer(0), sp.Rational(1, 2)},
         "collision hostile unexpectedly became a target-node inverse")

# Genuine universal-fibre wall: the universal pair is joined by X=1, and the
# extra point is itself Morse.  Q20 has a simple root at -1/6, so the total
# resultant valuation there is three.
universal_values = {
    Delta: sp.Integer(1),
    Phi: -sp.Rational(1176023, 2700),
    Theta: sp.Rational(32981, 450),
    eta: sp.Integer(1),
}
universal_Q, universal_G, universal_f, universal_h = normalized_residual(
    "universal-wall", universal_values
)
need(sp.factor(K.subs(universal_values)) == sp.Rational(5591, 90),
     "universal-wall K changed")
need(universal_Q.eval(t_universal) == 0,
     "universal-wall Q20(-1/6) changed")
need(universal_Q.diff().eval(t_universal) != 0,
     "universal-wall Q20 root ceased to be simple")
need(sp.gcd(universal_Q, universal_Q.diff()).degree() == 0,
     "universal-wall residual ceased to be squarefree")
universal_fibre_gcd = sp.gcd(
    sp.Poly(universal_f.subs(T, t_universal), X),
    sp.Poly(universal_h.subs(T, t_universal), X),
).monic()
need(sp.factor(universal_fibre_gcd.as_expr())
     == (X - 1) * (X**2 - 6),
     "universal-wall fibre gcd changed")
universal_hessian = sp.det(sp.hessian(universal_G, (X, T)))
need(universal_hessian.subs({X: 1, T: t_universal}) != 0,
     "universal-wall extra point ceased to be Morse")
need(universal_G.subs({X: 1, T: t_universal})
     not in {sp.Integer(0), sp.Rational(1, 2)},
     "universal-wall hostile unexpectedly became a target-node inverse")

packet = (8, 5, 4, 3, 2, 2, 1)
critical_length = 24
defect = 18
need(sum(packet) == 25 and sum(index - 1 for index in packet) == defect,
     "packet ledger changed")
need(2 * 21 - critical_length - 1 + 2 < 20,
     "finite response inequality changed")
need(2 * (25 - critical_length) < defect,
     "full response inequality changed")

digests = []
for residual in (generic_Q, collision_Q, universal_Q):
    primitive = sp.primitive(residual)[1]
    coefficient_text = ",".join(str(value) for value in primitive.all_coeffs())
    digests.append(sha256(coefficient_text.encode("ascii")).hexdigest())
semantic = (
    "filtration_top={(3,1),(0,3)};normalized=T^56*(6T+1)^2*Q20;"
    "collision_gcd=T-1;universal_fibre=(X-1)*(X^2-6);"
    "L=24;packet=(8,5,4,3,2,2,1);finite=(21,2);full=25"
)

print("THM4186_NORMALIZED_INDEPENDENT_ACCEPT")
print(f"checks={CHECKS}")
print("weight9_top=(p^3*y,y^3);eta=zeta=0_is_M_le_8")
print("normalized_resultant=T^56*(6T+1)^2*Q20")
print(f"generic_Q20_sha256={digests[0]}")
print(f"collision_Q20_sha256={digests[1]}")
print(f"universal_Q20_sha256={digests[2]}")
print("collision_hostile=gcd(Q20,Q20prime)=T-1;fibre=X*(X-1)")
print("universal_hostile=Q20(-1/6)=0;fibre=(X-1)*(X^2-6)")
print("normalized_hessian_bridge=T*detD(f,h)=detHess(G)+f*G_XT")
print(f"critical_length={critical_length};packet={packet};defect={defect}")
print("finite=(n=21,beta=2);full=n=25")
print(f"semantic_sha256={sha256(semantic.encode('ascii')).hexdigest()}")
