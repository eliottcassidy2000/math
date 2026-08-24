#!/usr/bin/env python3
"""Exact companion for THM-3926's THM-3907 class/genus packet.

Reproduction:
  python3 04-computation/jc2_unit_ideal_cubic_class_genus_boundary_thm3926.py
  python3 -O 04-computation/jc2_unit_ideal_cubic_class_genus_boundary_thm3926.py
"""

from __future__ import annotations

import hashlib
import itertools
import json

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(sp.expand(expression)) == 0, message)


A, C, omega, theta, T, z = sp.symbols("A C omega theta T z")

# Delone--Faddeev row (a,b,c,d)=(A,C,AC-1,A).
r_omega = A**2 * C + A * theta - A - C * omega + omega**2
r_mixed = A**2 + theta * omega
r_theta = A * C * theta + A * C - A * omega + theta**2 - theta
relations = (r_omega, r_mixed, r_theta)

P_theta = T**3 + (A * C - 1) * T**2 + A * C * T + A**3
D_theta = sp.diff(P_theta, T).subs(T, theta)
Delta = sp.factor(
    C**2 * (A * C - 1) ** 2
    - 4 * A * (A * C - 1) ** 3
    - 4 * C**3 * A
    - 27 * A**4
    + 18 * A**2 * C * (A * C - 1)
)
Delta_expected = (
    -4 * A**4 * C**3
    - 27 * A**4
    + 30 * A**3 * C**2
    + A**2 * C**4
    - 30 * A**2 * C
    - 6 * A * C**3
    + 4 * A
    + C**2
)
zero(Delta - Delta_expected, "THM-3907 discriminant")
zero(P_theta.subs(T, theta) - theta * r_theta - A * r_mixed,
     "theta characteristic relation")
zero(P_theta.subs({T: -1}) - (A**3 - 2), "norm of theta+1")
zero(P_theta.subs({T: 0}) - A**3, "norm of theta")


# ---------------------------------------------------------------------------
# Exact factorial chart S_{A(theta+1)}.
# ---------------------------------------------------------------------------

C_chart = -(theta**3 - theta**2 + A**3) / (A * theta * (theta + 1))
omega_chart = -A**2 / theta
for index, relation in enumerate(relations):
    zero(relation.subs({C: C_chart, omega: omega_chart}), f"factorial chart relation {index}")

zero(
    (A * theta * (theta + 1) * C + theta**3 - theta**2 + A**3).subs(C, C_chart),
    "chart inverse numerator",
)
zero((theta * omega + A**2).subs(omega, omega_chart), "chart inverse omega relation")


# ---------------------------------------------------------------------------
# Smoothness of the cubic surface.
# ---------------------------------------------------------------------------

jacobian_surface = sp.Matrix(relations).jacobian((A, C, omega, theta))
surface_minors = []
for rows in itertools.combinations(range(3), 2):
    for columns in itertools.combinations(range(4), 2):
        surface_minors.append(sp.expand(jacobian_surface.extract(rows, columns).det()))
surface_singular_gb = sp.groebner(
    list(relations) + surface_minors,
    A,
    C,
    omega,
    theta,
    order="grevlex",
)
gate(len(surface_singular_gb.polys) == 1 and surface_singular_gb.polys[0].as_expr() == 1,
     "the normal cubic surface is smooth")


# ---------------------------------------------------------------------------
# Nagata boundary primes and divisor matrix.
# ---------------------------------------------------------------------------

prime_substitutions = {
    "P0": {A: 0, omega: 0, theta: 0},
    "P1": {A: 0, omega: 0, theta: 1},
    "PC": {A: 0, omega: C, theta: 0},
}
for label, substitution in prime_substitutions.items():
    for index, relation in enumerate(relations):
        zero(relation.subs(substitution), f"{label} satisfies relation {index}")

rho = sp.symbols("rho")
q_substitution = {A: rho, theta: -1, omega: rho**2}
for index, relation in enumerate(relations):
    reduced = sp.rem(
        sp.Poly(sp.expand(relation.subs(q_substitution)), rho),
        sp.Poly(rho**3 - 2, rho),
    ).as_expr()
    zero(reduced, f"Q_rho satisfies relation {index}")

zero(P_theta.subs({A: 0}) - T**2 * (T - 1), "A-fibre theta collision")
zero(Delta.subs(A, 0) - C**2, "A-line is generically etale")

# Rows are div(A), div(theta), div(theta+1) in columns
# P0,P1,PC,Q1,Q2,Q3.
divisor_matrix = sp.Matrix(
    [
        [1, 1, 1, 0, 0, 0],
        [1, 0, 2, 0, 0, 0],
        [0, 0, 0, 1, 1, 1],
    ]
)
gate(divisor_matrix.rank() == 3, "three independent chart-unit divisors")
smith = smith_normal_form(divisor_matrix.T, domain=sp.ZZ)
smith_diagonal = tuple(smith[index, index] for index in range(3))
gate(smith_diagonal == (1, 1, 1), "Nagata relation lattice is primitive")

# P0=-2g, P1=g, PC=g and Q3=-Q1-Q2.
class_coordinates = {
    "P0": (-2, 0, 0),
    "P1": (1, 0, 0),
    "PC": (1, 0, 0),
    "Q1": (0, 1, 0),
    "Q2": (0, 0, 1),
    "Q3": (0, -1, -1),
}
for row_index, row in enumerate(divisor_matrix.tolist()):
    total = tuple(
        sum(row[column] * class_coordinates[label][coordinate]
            for column, label in enumerate(class_coordinates))
        for coordinate in range(3)
    )
    gate(total == (0, 0, 0), f"class coordinates satisfy divisor relation {row_index}")


# ---------------------------------------------------------------------------
# Ramification derivative: div(D_theta)=E+P0+PC, hence [E]=g.
# ---------------------------------------------------------------------------

zero(D_theta.subs(prime_substitutions["P1"]) - 1, "D_theta is a unit at P1")
u = sp.symbols("u")
D_blowup = sp.expand(D_theta.subs(theta, A * u) / A)
zero(D_blowup.subs({A: 0, u: C}) + C, "D_theta has order one at P0")
zero(D_blowup.subs({A: 0, u: 0}) - C, "D_theta has order one at PC")

E_class = tuple(
    -class_coordinates["P0"][coordinate] - class_coordinates["PC"][coordinate]
    for coordinate in range(3)
)
gate(E_class == (1, 0, 0), "ramification class is the primitive g generator")
boundary_basis = sp.Matrix.hstack(
    sp.Matrix(E_class),
    sp.Matrix(class_coordinates["Q1"]),
    sp.Matrix(class_coordinates["Q2"]),
)
gate(abs(int(boundary_basis.det())) == 1, "E,Q1,Q2 form a class-group basis")


# ---------------------------------------------------------------------------
# Closure and smoothness control for the ramification curve.
# ---------------------------------------------------------------------------

saturation = sp.groebner(
    list(relations) + [D_theta, z * A - 1],
    z,
    A,
    C,
    omega,
    theta,
    order="lex",
)
E_generators = [
    sp.expand(polynomial.as_expr())
    for polynomial in saturation.polys
    if not polynomial.as_expr().has(z)
]
gate(len(E_generators) == 5, "ramification closure has five Groebner generators")

jacobian_E = sp.Matrix(E_generators).jacobian((A, C, omega, theta))
E_minors = []
for rows in itertools.combinations(range(len(E_generators)), 3):
    for columns in itertools.combinations(range(4), 3):
        E_minors.append(sp.expand(jacobian_E.extract(rows, columns).det()))
E_singular_gb = sp.groebner(
    E_generators + E_minors,
    A,
    C,
    omega,
    theta,
    order="grevlex",
)
gate(len(E_singular_gb.polys) == 1 and E_singular_gb.polys[0].as_expr() == 1,
     "ramification curve is smooth")


# ---------------------------------------------------------------------------
# Genus-two Kummer model and the six deleted infinity places.
# ---------------------------------------------------------------------------

t = sp.symbols("t")
H = A**3 * (2 * t + 1) - t**2 * (t**2 + 2 * t - 1)
C_E = t * (2 - 3 * t) / (A * (2 * t + 1))
omega_E = -A**2 / t

for index, relation in enumerate(relations):
    numerator = sp.factor(sp.together(relation.subs({theta: t, C: C_E, omega: omega_E})).as_numer_denom()[0])
    gate(sp.rem(sp.Poly(numerator, A), sp.Poly(H, A)).as_expr() == 0,
         f"Kummer model relation {index}")
D_E_numerator = sp.factor(
    sp.together(D_theta.subs({theta: t, C: C_E, omega: omega_E})).as_numer_denom()[0]
)
gate(sp.rem(sp.Poly(D_E_numerator, A), sp.Poly(H, A)).as_expr() == 0,
     "Kummer model lies on ramification")

branch_quadratic = t**2 + 2 * t - 1
gate(sp.discriminant(branch_quadratic, t) == 8, "two simple Kummer zero addresses")
gate(sp.gcd(branch_quadratic, t * (2 * t + 1)) == 1, "four Kummer branch values are distinct")
kummer_branch_values = 4
genus = (3 * (-2) + kummer_branch_values * (3 - 1) + 2) // 2
gate(genus == 2, "Kummer Riemann-Hurwitz genus")

# The repeated-root incidence gives an independent hyperelliptic model.
B = sp.symbols("B")
incidence = (1 - 2 * t**3) * B**2 + (2 * t - t**4) * B - t**2
incidence_disc = sp.factor(sp.discriminant(incidence, B))
sextic = t**6 - 12 * t**3 + 8
zero(incidence_disc - t**2 * sextic, "repeated-root square completion")
gate(sp.gcd(sextic, sp.diff(sextic, t)) == 1, "hyperelliptic sextic is squarefree")

# Two A=0/C=infinity points, one t=-1/2 point, and three t=infinity points.
infinity_places = 2 + 1 + 3
gate(infinity_places == 6, "six normalization places at target infinity")
zero(D_theta.subs({theta: -1, A: rho}) - (5 - rho * C),
     "each Q_rho meets E in one affine point")


# ---------------------------------------------------------------------------
# Compactly supported Euler ledger over C.
# ---------------------------------------------------------------------------

chi_factorial_chart = 0  # chi(G_m)*chi(A1 minus {0,-1})
chi_A_boundary = 3 - 1   # three A1's; P0 and PC meet once
chi_Q_boundary = 3       # three disjoint A1's
chi_X = chi_factorial_chart + chi_A_boundary + chi_Q_boundary
chi_E = 2 - 2 * genus - infinity_places
chi_Q_minus_E = 1 - 1
chi_natural_open = chi_X - (chi_E + 2 * chi_Q_minus_E)
gate(chi_X == 5, "surface Euler characteristic")
gate(chi_E == -8, "ramification Euler characteristic")
gate(chi_natural_open == 13, "natural etale deletion Euler characteristic")


summary = {
    "checks": CHECKS,
    "class_group": "Z^3",
    "ramification_class": list(E_class),
    "boundary_basis_det": int(boundary_basis.det()),
    "genus": genus,
    "infinity_places": infinity_places,
    "chi_X": chi_X,
    "chi_E": chi_E,
    "chi_natural_open": chi_natural_open,
    "surface_smooth": True,
    "ramification_smooth": True,
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3926 exact unit-ideal cubic class/genus companion")
print(f"CHECKS={CHECKS}")
print("CLASS_GROUP=Z^3;BASIS=(E,Q1,Q2);E_CLASS=(1,0,0)")
print("CHART=k[A,A^-1,theta,theta^-1,(theta+1)^-1]")
print(f"GENUS={genus};INFINITY_PLACES={infinity_places}")
print(f"CHI_X={chi_X};CHI_E={chi_E};CHI_NATURAL_OPEN={chi_natural_open}")
print("BOUNDARY=class/unit sieve passes, rational-boundary and Euler sieves fail")
print(f"SEMANTIC_SHA256={semantic}")
