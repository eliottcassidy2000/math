#!/usr/bin/env python3
"""Primary normalized-coordinate certificate for THM-4189.

This path works in the normalized (X,T) chart.  It computes the complete
symbolic Theta-zero critical resultant, reconstructs the Newton polygon and
every boundary face, and checks two exact hostile projection walls: a
reduced algebraic fibre with two distinct Morse points sharing T, and an
extra Morse point in the universal T=-1/6 fibre.
"""

from __future__ import annotations

from hashlib import sha256
from math import gcd

import sympy as sp


CHECKS = 0


def need(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly: sp.Expr, variable: sp.Symbol) -> int:
    terms = sp.Poly(poly, variable).terms()
    need(bool(terms), "zero polynomial has no valuation")
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def exact_quotient(
    numerator: sp.Expr,
    denominator: sp.Expr,
    variable: sp.Symbol,
    label: str,
) -> sp.Poly:
    quotient = sp.cancel(numerator / denominator)
    need(sp.denom(quotient) == 1, f"{label}: quotient is not polynomial")
    return sp.Poly(quotient, variable)


def convex_hull(points: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    ordered = sorted(set(points))

    def cross(origin, left, right):
        return ((left[0] - origin[0]) * (right[1] - origin[1])
                - (left[1] - origin[1]) * (right[0] - origin[0]))

    lower: list[tuple[int, int]] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(points: tuple[tuple[int, int], ...]):
    vertices = convex_hull(points)
    area2 = abs(sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]))
        for index in range(len(vertices))
    )
    need((area2 - boundary + 2) % 2 == 0, "Pick parity changed")
    return vertices, area2, boundary, (area2 - boundary + 2) // 2


def selected_face(poly, x, y, u, v, level):
    return sp.factor(sum(
        coefficient * x**i * y**j
        for (i, j), coefficient in sp.Poly(poly, x, y).terms()
        if u * i + v * j == level
    ))


def reduced_mod(expr, modulus, variable):
    numerator, denominator = sp.together(sp.cancel(expr)).as_numer_denom()
    need(sp.resultant(denominator, modulus, variable) != 0,
         "denominator meets reduction modulus")
    return sp.rem(sp.Poly(numerator, variable), sp.Poly(modulus, variable))


X, T = sp.symbols("X T")
Delta, Phi, eta = sp.symbols("Delta Phi eta")
K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
P = T + X**2 * T**2
Y = X * T * P

# Complete P-only exact-M=9 source on zeta=Theta=0.
G = sp.expand(
    -X**2 * T / 2
    - 3 * P
    + sp.Rational(8, 3) * P**2
    - sp.Rational(1376, 135) * P**3
    + K * Y**2
    + Phi * P**2 * Y
    + Delta * P**4
    + eta * P**3 * Y
)
f = sp.cancel(sp.diff(G, X) / T)
h = sp.diff(G, T)
need(sp.denom(f) == 1, "G_X/T is not polynomial")
need(sp.factor(sp.diff(G, X) - T * f) == 0, "G_X=T*f changed")

# Exact-M=9 completeness and normalized infinity firewall.
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
need((sp.degree(f, X), sp.degree(h, X)) == (8, 9),
     "normalized critical degrees changed")
need(sp.factor(sp.Poly(f, X).LC() - 9 * eta * T**8) == 0,
     "f leading row changed")
need(sp.factor(sp.Poly(h, X).LC() - 9 * eta * T**8) == 0,
     "h leading row changed")

normalized_hessian = sp.det(sp.hessian(G, (X, T)))
normalized_jacobian = sp.det(sp.Matrix((
    (sp.diff(f, X), sp.diff(f, T)),
    (sp.diff(h, X), sp.diff(h, T)),
)))
need(sp.factor(
    T * normalized_jacobian - normalized_hessian - f * sp.diff(G, X, T)
) == 0, "normalized Morse-resultant bridge changed")

# Universal coordinate fibres and their exact values/Morse determinants.
need(sp.factor(f.subs(T, 0) + X) == 0, "T=0 f row changed")
need(sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0,
     "T=0 h row changed")
need(reduced_mod(normalized_hessian.subs(T, 0) - 6, X**2 + 6, X).is_zero,
     "T=0 universal pair ceased to be Morse")

t_universal = -sp.Rational(1, 6)
for expression, expected, label in (
    (f, 0, "f"),
    (h, 0, "h"),
    (G, sp.Rational(1, 2), "value"),
    (normalized_hessian, -6, "Hessian"),
):
    need(reduced_mod(
        expression.subs(T, t_universal) - expected, X**2 - 6, X
    ).is_zero, f"T=-1/6 universal {label} changed")

# Complete symbolic normalized resultant.  Both endpoints are units exactly
# under eta*K!=0, so no coefficient or critical subwall survives here.
resultant = sp.resultant(f, h, X)
need(valuation(resultant, T) == 56, "normalized T-artifact changed")
Q19 = exact_quotient(
    resultant, T**56 * (6 * T + 1)**2, T, "normalized resultant"
)
need(Q19.degree() == 19, "normalized residual degree changed")
need(sp.factor(Q19.TC() + sp.Rational(3**15, 2**7) * eta**7) == 0,
     "normalized residual constant changed")
need(sp.factor(Q19.LC() - 944784 * K**6 * eta**7) == 0,
     "normalized residual leading row changed")

# Newton polygon and exact face audit in the independent (s,p) boundary
# variables.  Phi is absent from every boundary face.
s, p, Q = sp.symbols("s p Q")
H = sp.expand(
    -3 * p
    + sp.Rational(8, 3) * p**2
    - sp.Rational(1376, 135) * p**3
    + K * s**2 * p**2
    + Phi * s * p**3
    + Delta * p**4
    + eta * s * p**4
)
F_Q = sp.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)
control = {Delta: 2, Phi: 5, eta: 7, Q: 11}
support = tuple(
    powers for powers, coefficient
    in sp.Poly(F_Q.subs(control), s, p).terms()
    if coefficient != 0
)
vertices = ((0, 1), (2, 0), (4, 2), (3, 4), (1, 5), (0, 5))
need(polygon_ledger(support) == (vertices, 28, 10, 10),
     "Theta-zero Newton polygon changed")

faces = (
    ((1, 2, 2), s**2 * (1 - Q / 2) - p, "low"),
    ((-1, 1, -2),
     s**2 * ((1 - Q / 2) - K * Q * (s * p)**2), "carrier"),
    ((-2, -1, -10),
     -Q * p**2 * s**3 * (K * s + eta * p**2), "merged"),
    ((-1, -2, -11), Q * eta * p**4 * s * (p - s**2), "top"),
    ((0, -1, -5), Q * p**5 * (Delta + eta * s), "horizontal"),
    ((1, 0, 0),
     p * (-1 + Q * (-3 * p + sp.Rational(8, 3) * p**2
                     - sp.Rational(1376, 135) * p**3 + Delta * p**4)),
     "vertical affine"),
)
for normal, expected, label in faces:
    need(sp.factor(selected_face(F_Q, s, p, *normal) - expected) == 0,
         f"{label} face changed")
need(sp.Poly(K * s + eta * p**2, s).degree() == 1,
     "merged face is not primitive-linear")
need(sp.diff(Delta + eta * s, s) == eta,
     "horizontal face derivative changed")
need((Delta + eta * s).subs(s, -Delta / eta) == 0,
     "horizontal torus root changed")

packet = (8, 7, 4, 2, 2, 1)
defect = sum(index - 1 for index in packet)
critical_length = 19 + 2 + 2
full_n, finite_n, beta = 24, 20, 2
need((sum(packet), defect, critical_length) == (24, 18, 23),
     "packet/length ledger changed")
need(2 * (full_n - critical_length) < defect,
     "full commutator contradiction changed")
need(2 * finite_n - critical_length - 1 + beta < finite_n - 1,
     "finite merger contradiction changed")

# Hostile A: an exact reduced fibre with two distinct Morse source points
# sharing T.  The points separate in p=T+(XT)^2, so this is precisely a
# projection collision, not nonreduced source geometry.
tau = sp.symbols("tau")
P6_expr = (
    16512 * tau**6 - 174912 * tau**5 - 235020 * tau**4
    - 35764 * tau**3 + 10607 * tau**2 + 4530 * tau + 1485
)
P6 = sp.Poly(P6_expr, tau)
need(sp.factor(P6_expr) == P6_expr, "collision polynomial became reducible")
left, right = sp.Rational(315, 1000), sp.Rational(316, 1000)
need(P6.eval(left) > 0 and P6.eval(right) < 0,
     "collision root interval signs changed")
need(P6.count_roots(left, right) == 1,
     "collision interval no longer isolates one real root")

N_tau = (
    33024 * tau**6 + 265344 * tau**5 + 263496 * tau**4
    + 33064 * tau**3 + 3973 * tau**2 - 2100 * tau - 1485
)
collision_values = {
    Delta: (1376 * tau**2 - 240 * tau + 135) / (180 * tau**3),
    Phi: N_tau / (1620 * tau**4 * (tau + 1)**2 * (3 * tau + 1)),
    eta: -N_tau / (1620 * tau**5 * (tau + 1)**2 * (3 * tau + 1)),
}
collision_K = sp.factor(K.subs(collision_values))
for label, expression in (
    ("Delta", collision_values[Delta]),
    ("eta", collision_values[eta]),
    ("K", collision_K),
):
    numerator = sp.together(expression).as_numer_denom()[0]
    need(sp.gcd(P6, sp.Poly(numerator, tau)).degree() == 0,
         f"collision {label} can vanish at tau")

collision_G = sp.cancel(G.subs(collision_values))
collision_f = sp.cancel(f.subs(collision_values))
collision_h = sp.cancel(h.subs(collision_values))
collision_hessian = sp.cancel(normalized_hessian.subs(collision_values))

def tau_remainder(expr):
    numerator = sp.together(expr).as_numer_denom()[0]
    return sp.rem(sp.Poly(numerator, tau), P6)


for x_value in (0, 1):
    need(tau_remainder(collision_f.subs({X: x_value, T: tau})).is_zero,
         "collision f point changed")
    need(tau_remainder(collision_h.subs({X: x_value, T: tau})).is_zero,
         "collision h point changed")
    hessian_numerator = sp.together(
        collision_hessian.subs({X: x_value, T: tau})
    ).as_numer_denom()[0]
    need(sp.gcd(P6, sp.Poly(hessian_numerator, tau)).degree() == 0,
         "collision point ceased to be Morse")
    for value in (0, sp.Rational(1, 2)):
        value_numerator = sp.together(
            collision_G.subs({X: x_value, T: tau}) - value
        ).as_numer_denom()[0]
        need(sp.gcd(P6, sp.Poly(value_numerator, tau)).degree() == 0,
             "collision point became a target-node inverse")

fibre_q = sp.Poly(X * (X - 1), X)
quotients = []
for name, fibre_poly in (
    ("f", collision_f.subs(T, tau)),
    ("h", collision_h.subs(T, tau)),
):
    numerator = sp.together(fibre_poly).as_numer_denom()[0]
    quotient, remainder = sp.div(sp.Poly(numerator, X), fibre_q)
    for coefficient in remainder.all_coeffs():
        need(sp.rem(sp.Poly(coefficient, tau), P6).is_zero,
             f"collision {name} fibre lost X(X-1)")
    quotients.append(quotient)
quotient_resultant = sp.resultant(
    quotients[0].as_expr(), quotients[1].as_expr(), X
)
quotient_numerator = sp.together(quotient_resultant).as_numer_denom()[0]
need(sp.gcd(P6, sp.Poly(quotient_numerator, tau)).degree() == 0,
     "collision fibre acquired an additional common root")
need(sp.rem(sp.Poly(tau**2, tau), P6).as_expr() != 0,
     "collision p-coordinates ceased to be distinct")

# Hostile B: an extra reduced point joins the two universal points at
# T=-1/6, while the residual Q19 itself stays squarefree.
universal_values = {
    Delta: sp.Integer(1),
    Phi: -sp.Rational(12416, 25),
    eta: -sp.Rational(32856, 125),
}
need(sp.factor(K.subs(universal_values)) == sp.Rational(5591, 90),
     "universal hostile K changed")
universal_G = sp.expand(G.subs(universal_values))
universal_f = sp.cancel(sp.diff(universal_G, X) / T)
universal_h = sp.diff(universal_G, T)
universal_resultant = sp.resultant(universal_f, universal_h, X)
universal_Q = exact_quotient(
    universal_resultant, T**56 * (6 * T + 1)**2, T,
    "universal hostile resultant",
)
need(universal_Q.degree() == 19, "universal hostile degree changed")
need(universal_Q.eval(t_universal) == 0,
     "universal hostile Q19(-1/6) changed")
need(universal_Q.diff().eval(t_universal) != 0,
     "universal hostile residual root ceased to be simple")
need(sp.gcd(universal_Q, universal_Q.diff()).degree() == 0,
     "universal hostile residual ceased to be squarefree")
universal_fibre_gcd = sp.gcd(
    sp.Poly(universal_f.subs(T, t_universal), X),
    sp.Poly(universal_h.subs(T, t_universal), X),
).monic()
need(sp.factor(universal_fibre_gcd.as_expr()) == (X - 1) * (X**2 - 6),
     "universal hostile fibre changed")
universal_hessian = sp.det(sp.hessian(universal_G, (X, T)))
need(universal_hessian.subs({X: 1, T: t_universal})
     == -sp.Rational(1375271635, 68024448),
     "universal hostile extra Hessian changed")
need(universal_G.subs({X: 1, T: t_universal})
     == sp.Rational(2050529, 5038848),
     "universal hostile extra value changed")

coefficient_text = "|".join(sp.sstr(value) for value in Q19.all_coeffs())
symbolic_digest = sha256(coefficient_text.encode("ascii")).hexdigest()
semantic = (
    "wall=zeta=Theta=0;normalized=T^56*(6T+1)^2*Q19;"
    "L=23;packet=(8,7,4,2,2,1);finite=(20,2);full=24;"
    "collision=fibre_X*(X-1);universal=(X-1)*(X^2-6)"
)

print("THM4189_NORMALIZED_PRIMARY_EXACT_ACCEPT")
print(f"checks={CHECKS}")
print("wall=zeta=Theta=0;eta*Delta*K!=0;Phi_arbitrary")
print("normalized_resultant=T^56*(6T+1)^2*Q19")
print(f"Q19_constant={sp.factor(Q19.TC())}")
print(f"Q19_leading={sp.factor(Q19.LC())}")
print(f"symbolic_coefficient_sha256={symbolic_digest}")
print("polygon=((0,1),(2,0),(4,2),(3,4),(1,5),(0,5));ledger=(28,10,10)")
print("packet=(8,7,4,2,2,1);defect=18;critical_length=23")
print("carrier=K*W^2=q-1/2;finite=(20,2);full=24")
print("collision_hostile=T=tau;fibre=X*(X-1);source_p=tau,tau+tau^2")
print("universal_hostile=T=-1/6;fibre=(X-1)*(X^2-6);Q19_squarefree")
print(f"semantic_sha256={sha256(semantic.encode('ascii')).hexdigest()}")
