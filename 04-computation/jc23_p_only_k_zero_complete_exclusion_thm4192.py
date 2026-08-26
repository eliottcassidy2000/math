#!/usr/bin/env python3
"""Primary normalized-coordinate exact certificate for THM-4192.

This path is deliberately independent of the source-coordinate critical pair
used by the companion audit.
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
    return sp.rem(sp.Poly(sp.expand(expr), variable), sp.Poly(modulus, variable))


X, T = sp.symbols("X T")
Phi, Theta, eta = sp.symbols("Phi Theta eta")
Delta0 = sp.Rational(5696, 105)
need(sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta0 == 0,
     "K=0 coefficient value changed")

P = T + X**2 * T**2
Y = X * T * P
G = sp.expand(
    -X**2 * T / 2
    - 3 * P
    + sp.Rational(8, 3) * P**2
    - sp.Rational(1376, 135) * P**3
    + Phi * P**2 * Y
    + Delta0 * P**4
    + Theta * P * Y**2
    + eta * P**3 * Y
)

# Exact-M=9 completeness: zeta=0 and eta!=0 leaves the P-only row.
admissible = {
    (i, j)
    for i in range(5)
    for j in range(4)
    if 0 < 2 * i + 3 * j <= 9 and (i, j) not in {(0, 1), (1, 1)}
}
need(admissible == {
    (1, 0), (2, 0), (3, 0), (0, 2), (2, 1),
    (4, 0), (1, 2), (3, 1), (0, 3),
}, "complete source monomial universe changed")
need({pair for pair in admissible if 2 * pair[0] + 3 * pair[1] == 9}
     == {(3, 1), (0, 3)}, "weight-nine top row changed")

f = sp.cancel(sp.diff(G, X) / T)
h = sp.diff(G, T)
need(sp.denom(f) == 1 and sp.factor(sp.diff(G, X) - T * f) == 0,
     "normalized critical pair changed")
need((sp.degree(f, X), sp.degree(h, X)) == (8, 9),
     "normalized critical degrees changed")
need(sp.factor(sp.Poly(f, X).LC() - 9 * eta * T**8) == 0,
     "f leading row changed")
need(sp.factor(sp.Poly(h, X).LC() - 9 * eta * T**8) == 0,
     "h leading row changed")

hessian = sp.det(sp.hessian(G, (X, T)))
critical_jacobian = sp.det(sp.Matrix((
    (sp.diff(f, X), sp.diff(f, T)),
    (sp.diff(h, X), sp.diff(h, T)),
)))
need(sp.factor(T * critical_jacobian - hessian - f * sp.diff(G, X, T)) == 0,
     "normalized Hessian bridge changed")

# Both universal pairs persist in all three coefficient strata.
need(sp.factor(f.subs(T, 0) + X) == 0,
     "T=0 f row changed")
need(sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0,
     "T=0 h row changed")
need(reduced_mod(hessian.subs(T, 0) - 6, X**2 + 6, X).is_zero,
     "T=0 universal pair ceased to be Morse")
for expression, expected, label in (
    (f, 0, "f"),
    (h, 0, "h"),
    (G, sp.Rational(1, 2), "value"),
    (hessian, -6, "Hessian"),
):
    need(reduced_mod(
        expression.subs(T, -sp.Rational(1, 6)) - expected,
        X**2 - 6,
        X,
    ).is_zero, f"T=-1/6 universal {label} changed")

# One symbolic resultant contains the three exhaustive coefficient strata.
resultant = sp.resultant(f, h, X)
base_factor = T**56 * (6 * T + 1)**2

need(valuation(resultant, T) == 56, "Theta-nonzero T artifact changed")
Q18 = exact_quotient(resultant, base_factor, T, "Theta-nonzero resultant")
need(Q18.degree() == 18, "Theta-nonzero residual degree changed")
need(sp.factor(Q18.TC() + sp.Rational(3**15, 2**7) * eta**7) == 0,
     "Q18 constant changed")
need(sp.factor(Q18.LC() - 656100 * Theta**6 * eta**5) == 0,
     "Q18 leading row changed")

resultant0 = sp.expand(resultant.subs(Theta, 0))
need(valuation(resultant0, T) == 56, "Theta-zero T artifact changed")
Q15 = exact_quotient(resultant0, base_factor, T, "Theta-zero resultant")
need(Q15.degree() == 15, "Theta-zero residual degree changed")
need(sp.factor(Q15.TC() + sp.Rational(3**15, 2**7) * eta**7) == 0,
     "Q15 constant changed")
need(sp.factor(Q15.LC() - sp.Rational(6561, 4) * Phi**6 * eta**5) == 0,
     "Q15 leading row changed")

resultant00 = sp.expand(resultant0.subs(Phi, 0))
need(valuation(resultant00, T) == 56, "terminal T artifact changed")
Q14 = exact_quotient(resultant00, base_factor, T, "terminal resultant")
need(Q14.degree() == 14, "terminal residual degree changed")
need(sp.factor(Q14.TC() + sp.Rational(3**15, 2**7) * eta**7) == 0,
     "Q14 constant changed")
need(sp.factor(Q14.LC() - sp.Rational(59049, 4) * eta**9) == 0,
     "Q14 leading row changed")

# Generic squarefree controls are not hypotheses; they prove nonvacuity and
# attack accidental degree/collision artifacts in each stratum.
controls = (
    (sp.Poly(Q18.as_expr().subs({Theta: 1, Phi: 0, eta: 1}), T), "Q18"),
    (sp.Poly(Q15.as_expr().subs({Phi: 1, eta: 1}), T), "Q15"),
    (sp.Poly(Q14.as_expr().subs(eta, 1), T), "Q14"),
)
for control_poly, label in controls:
    need(sp.gcd(control_poly, control_poly.diff()).degree() == 0,
         f"{label} squarefree control failed")
    need(control_poly.eval(-sp.Rational(1, 6)) != 0,
         f"{label} universal-fibre control failed")

# Newton polygons, exact faces, and residue indices.  The final vertical face
# has s=0,p!=0, hence X=0,T=p and is affine rather than a puncture.
s, p, Q = sp.symbols("s p Q")
H = sp.expand(
    -3 * p
    + sp.Rational(8, 3) * p**2
    - sp.Rational(1376, 135) * p**3
    + Phi * s * p**3
    + Delta0 * p**4
    + Theta * s**2 * p**3
    + eta * s * p**4
)
F_Q = sp.expand((s**2 - p) * (1 - Q * H) - Q * s**2 / 2)
vertical_face = p * (
    -1 + Q * (-3 * p + sp.Rational(8, 3) * p**2
              - sp.Rational(1376, 135) * p**3 + Delta0 * p**4)
)

strata = (
    {
        "name": "Theta_nonzero",
        "poly": F_Q,
        "control": {Theta: 7, Phi: 5, eta: 11, Q: 13},
        "ledger": (
            ((0, 1), (2, 0), (4, 3), (3, 4), (1, 5), (0, 5)),
            27, 9, 10,
        ),
        "faces": (
            ((1, 2, 2), s**2 * (1 - Q / 2) - p),
            ((-3, 2, -6),
             s**2 * ((1 - Q / 2) - Q * Theta * s**2 * p**3)),
            ((-1, -1, -7), -Q * p**3 * s**3 * (Theta * s + eta * p)),
            ((-1, -2, -11), Q * eta * p**4 * s * (p - s**2)),
            ((0, -1, -5), Q * p**5 * (Delta0 + eta * s)),
        ),
        "packet": (8, 5, 5, 4, 1),
        "length": 22,
        "defect": 18,
        "gap": 1,
    },
    {
        "name": "Theta_zero_Phi_nonzero",
        "poly": sp.expand(F_Q.subs(Theta, 0)),
        "control": {Phi: 5, eta: 11, Q: 13},
        "ledger": (
            ((0, 1), (2, 0), (3, 3), (3, 4), (1, 5), (0, 5)),
            23, 9, 8,
        ),
        "faces": (
            ((1, 2, 2), s**2 * (1 - Q / 2) - p),
            ((-3, 1, -6),
             s**2 * ((1 - Q / 2) - Q * Phi * s * p**3)),
            ((-1, 0, -3), -Q * s**3 * p**3 * (Phi + eta * p)),
            ((-1, -2, -11), Q * eta * p**4 * s * (p - s**2)),
            ((0, -1, -5), Q * p**5 * (Delta0 + eta * s)),
        ),
        "packet": (8, 4, 4, 2, 1),
        "length": 19,
        "defect": 14,
        "gap": 0,
    },
    {
        "name": "Theta_Phi_zero",
        "poly": sp.expand(F_Q.subs({Theta: 0, Phi: 0})),
        "control": {eta: 11, Q: 13},
        "ledger": (
            ((0, 1), (2, 0), (3, 4), (1, 5), (0, 5)),
            22, 8, 8,
        ),
        "faces": (
            ((1, 2, 2), s**2 * (1 - Q / 2) - p),
            ((-4, 1, -8),
             s**2 * ((1 - Q / 2) - Q * eta * s * p**4)),
            ((-1, -2, -11), Q * eta * p**4 * s * (p - s**2)),
            ((0, -1, -5), Q * p**5 * (Delta0 + eta * s)),
        ),
        "packet": (8, 5, 4, 1),
        "length": 18,
        "defect": 14,
        "gap": 0,
    },
)

summary_rows = []
for row in strata:
    support = tuple(
        powers for powers, coefficient
        in sp.Poly(row["poly"].subs(row["control"]), s, p).terms()
        if coefficient
    )
    need(polygon_ledger(support) == row["ledger"],
         f"{row['name']} polygon ledger changed")
    indices = []
    for (u, v, level), expected in row["faces"]:
        need(sp.factor(selected_face(row["poly"], s, p, u, v, level)
                       - expected) == 0,
             f"{row['name']} face {(u, v, level)} changed")
        indices.append(u + v - level)
    need(tuple(sorted(indices, reverse=True)) == row["packet"],
         f"{row['name']} residue packet changed")
    need(sp.factor(selected_face(row["poly"], s, p, 1, 0, 0)
                   - vertical_face) == 0,
         f"{row['name']} affine vertical face changed")
    need(all(gcd(abs(b[0] - a[0]), abs(b[1] - a[1])) == 1
             for a, b in zip(row["ledger"][0][:-1], row["ledger"][0][1:])),
         f"{row['name']} nonvertical edge ceased to be primitive")
    need(sum(row["packet"]) - row["length"] == row["gap"],
         f"{row['name']} n-L gap changed")
    need(sum(index - 1 for index in row["packet"]) == row["defect"],
         f"{row['name']} defect changed")
    need(2 * (sum(row["packet"]) - row["length"]) < row["defect"],
         f"{row['name']} commutator contradiction changed")
    summary_rows.append(
        f"{row['name']}:{row['packet']}:L{row['length']}:g{row['ledger'][3]}"
    )

# Direct algebraic checks behind rationality of every nonaffine place.  Each
# torus face equation is linear in the primitive edge monomial.
w = sp.symbols("w")
linear_face_polynomials = (
    1 - Q / 2 - Q * Theta * w,
    Theta * w + eta,
    1 - Q / 2 - Q * Phi * w,
    Phi + eta * w,
    1 - Q / 2 - Q * eta * w,
    w - 1,
    Delta0 + eta * w,
)
for face_poly in linear_face_polynomials:
    need(sp.degree(face_poly, w) == 1,
         "a claimed rational face ceased to be primitive-linear")

# Projection hostile A: two distinct Morse critical points share one T-value.
# The real algebraic parameter tau is isolated exactly, and the quotient-fibre
# resultant proves that X=0,1 are the complete common fibre.
tau = sp.symbols("tau")
C3 = sp.Poly(68352 * tau**3 - 9632 * tau**2 + 1680 * tau - 945, tau)
left, right = sp.Rational(51, 200), sp.Rational(32, 125)
need(C3.eval(left) < 0 and C3.eval(right) > 0,
     "collision interval signs changed")
need(C3.count_roots(left, right) == 1 and sp.discriminant(C3) < 0,
     "collision interval no longer isolates the unique real root")
B_tau = (
    68352 * tau**6 + 205056 * tau**5 + 195424 * tau**4
    + 49088 * tau**3 - 7952 * tau**2 + 1680 * tau - 315
)
eta_tau = sp.factor(
    -2 * B_tau / (315 * tau**4 * (tau + 1)**2 * (5 * tau + 2))
)
Phi_tau = sp.factor(-eta_tau * tau)
Theta_tau = sp.factor((
    136704 * tau**7 + 410112 * tau**6 + 390848 * tau**5
    + 98176 * tau**4 - 15904 * tau**3 + 3360 * tau**2
    + 945 * tau + 630
) / (630 * tau**4 * (tau + 1)**2 * (5 * tau + 2)))
collision_values = {Phi: Phi_tau, Theta: Theta_tau, eta: eta_tau}
for label, expression in collision_values.items():
    numerator, denominator = sp.together(expression).as_numer_denom()
    need(sp.gcd(C3, sp.Poly(numerator, tau)).degree() == 0,
         f"collision {label} can vanish")
    need(sp.gcd(C3, sp.Poly(denominator, tau)).degree() == 0,
         f"collision {label} denominator can vanish")


def tau_remainder(expr):
    numerator, denominator = sp.together(expr).as_numer_denom()
    need(sp.gcd(C3, sp.Poly(denominator, tau)).degree() == 0,
         "collision expression denominator meets cubic")
    return sp.rem(sp.Poly(numerator, tau), C3)


collision_G = sp.cancel(G.subs(collision_values))
collision_f = sp.cancel(f.subs(collision_values))
collision_h = sp.cancel(h.subs(collision_values))
collision_hessian = sp.cancel(hessian.subs(collision_values))
for x_value in (0, 1):
    need(tau_remainder(collision_f.subs({X: x_value, T: tau})).is_zero,
         "collision f point changed")
    need(tau_remainder(collision_h.subs({X: x_value, T: tau})).is_zero,
         "collision h point changed")
    hessian_numerator = sp.together(
        collision_hessian.subs({X: x_value, T: tau})
    ).as_numer_denom()[0]
    need(sp.gcd(C3, sp.Poly(hessian_numerator, tau)).degree() == 0,
         "collision point ceased to be Morse")
    for value in (0, sp.Rational(1, 2)):
        value_numerator = sp.together(
            collision_G.subs({X: x_value, T: tau}) - value
        ).as_numer_denom()[0]
        need(sp.gcd(C3, sp.Poly(value_numerator, tau)).degree() == 0,
             "collision point became a target-node inverse")

fibre = sp.Poly(X * (X - 1), X)
quotients = []
for name, expression in (("f", collision_f), ("h", collision_h)):
    numerator = sp.together(expression.subs(T, tau)).as_numer_denom()[0]
    quotient, remainder = sp.div(sp.Poly(numerator, X), fibre)
    for coefficient in remainder.all_coeffs():
        need(sp.rem(sp.Poly(coefficient, tau), C3).is_zero,
             f"collision {name} fibre lost X(X-1)")
    quotients.append(quotient)
quotient_resultant = sp.resultant(
    quotients[0].as_expr(), quotients[1].as_expr(), X
)
quotient_numerator = sp.together(quotient_resultant).as_numer_denom()[0]
need(sp.gcd(C3, sp.Poly(quotient_numerator, tau)).degree() == 0,
     "collision fibre acquired an additional common root")
need(not tau_remainder(tau**2).is_zero,
     "collision source p-coordinates ceased to be distinct")

# Projection hostile B: a third Morse point joins the universal T=-1/6
# fibre in row II.  The residual Q15 is nevertheless squarefree.
t_universal = -sp.Rational(1, 6)
universal_values = {
    Phi: sp.Rational(201808, 1575),
    Theta: sp.Integer(0),
    eta: sp.Rational(2766784, 875),
}
universal_Q15 = sp.Poly(Q15.as_expr().subs(universal_values), T)
need(universal_Q15.eval(t_universal) == 0,
     "universal hostile Q15(-1/6) changed")
need(universal_Q15.diff().eval(t_universal) != 0,
     "universal hostile residual root ceased to be simple")
need(sp.gcd(universal_Q15, universal_Q15.diff()).degree() == 0,
     "universal hostile residual ceased to be squarefree")
universal_G = sp.expand(G.subs(universal_values))
universal_f = sp.cancel(sp.diff(universal_G, X) / T)
universal_h = sp.diff(universal_G, T)
universal_fibre_gcd = sp.gcd(
    sp.Poly(universal_f.subs(T, t_universal), X),
    sp.Poly(universal_h.subs(T, t_universal), X),
).monic()
need(sp.factor(universal_fibre_gcd.as_expr()) == (X - 1) * (X**2 - 6),
     "universal hostile fibre changed")
universal_hessian = sp.det(sp.hessian(universal_G, (X, T)))
need(universal_G.subs({X: 1, T: t_universal})
     == sp.Rational(42257, 91854),
     "universal hostile extra value changed")
need(universal_hessian.subs({X: 1, T: t_universal})
     == -sp.Rational(27296015083, 937461924),
     "universal hostile extra Hessian changed")

coefficient_digest = sha256("|".join(
    sp.sstr(coefficient)
    for polynomial in (Q18, Q15, Q14)
    for coefficient in polynomial.all_coeffs()
).encode("ascii")).hexdigest()
semantic = (
    "K=0;Delta=5696/105;normalized=(Q18,Q15,Q14);"
    "lengths=(22,19,18);packets=((8,5,5,4,1),(8,4,4,2,1),(8,5,4,1));"
    "all_boundary_places=rational;n-L=(1,0,0);defects=(18,14,14);"
    "projection_hostiles=repeated_T,universal_fibre"
)

print("THM4192_K_ZERO_NORMALIZED_EXACT_ACCEPT")
print(f"checks={CHECKS}")
print("strata=Theta!=0;Theta=0,Phi!=0;Theta=Phi=0")
print("normalized_residual_degrees=18,15,14")
print("normalized_constants=-3^15*eta^7/2^7")
print("normalized_leading=656100*Theta^6*eta^5;6561*Phi^6*eta^5/4;59049*eta^9/4")
print("polygons=(27,9,10);(23,9,8);(22,8,8)")
print("packets=(8,5,5,4,1);(8,4,4,2,1);(8,5,4,1)")
print("critical_lengths=22,19,18;degrees=23,19,18;defects=18,14,14")
print("boundary_fields=all_rational;commutator_upper=2,0,0")
print("projection_collision_hostile=T=tau;fibre=X*(X-1);two_Morse_points")
print("projection_universal_hostile=row_II;T=-1/6;fibre=(X-1)*(X^2-6);Q15_squarefree")
print(f"symbolic_coefficient_sha256={coefficient_digest}")
print(f"semantic_sha256={sha256(semantic.encode('ascii')).hexdigest()}")
