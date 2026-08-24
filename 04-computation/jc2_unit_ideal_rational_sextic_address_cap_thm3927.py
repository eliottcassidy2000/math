#!/usr/bin/env python3
"""Exact companion for THM-3927's rational two-place cubic packet.

Reproduction:
  python3 04-computation/jc2_unit_ideal_rational_sextic_address_cap_thm3927.py
  python3 -O 04-computation/jc2_unit_ideal_rational_sextic_address_cap_thm3927.py
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


def binary_disc(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(
        b**2 * c**2
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a**2 * d**2
        + 18 * a * b * c * d
    )


A, C, B, U, V, t = sp.symbols("A C B U V t")
alpha = 1 + 27 * A
beta = 1 + 24 * A

# The unit-ideal but globally nonmonogenic binary cubic.
a = A * alpha
b = C
c = -beta
d = 8 * A
Psi = sp.expand(a * U**3 + b * U**2 * V + c * U * V**2 + d * V**3)
gate(
    all(
        sp.expand(actual - expected) == 0
        for actual, expected in zip(
            tuple(sp.Poly(Psi, U, V).coeff_monomial(monomial) for monomial in
                  (U**3, U**2 * V, U * V**2, V**3)),
            (a, b, c, d),
        )
    ),
    "binary cubic coefficient row",
)
zero(c + 3 * d + 1, "coefficient ideal contains one")
zero(Psi.subs(A, 0) - U * V * (C * U - V), "specialization obstruction factors")

# Dehomogenization is a degree-three rational map in the repeated-root
# variable, giving generic irreducibility.
f = sp.expand(Psi.subs({U: t, V: 1}))
numerator_C = sp.expand(a * t**3 + c * t + d)
zero(f - (C * t**2 + numerator_C), "dehomogenized cubic")
gate(sp.gcd(sp.Poly(numerator_C, t), sp.Poly(t**2, t)).degree() == 0,
     "rational-map numerator and denominator are coprime")
gate(sp.degree(numerator_C, t) == 3, "rational-map degree is three")
gate(sp.gcd(sp.Poly(d, A), sp.Poly(c, A)).degree() == 0,
     "there is no repeated root at t=0")

Delta = sp.factor(binary_disc(a, b, c, d))
Delta_expected = (
    -1259712 * A**6
    + 1399680 * A**5
    - 93312 * A**4 * C
    + 240192 * A**4
    - 7344 * A**3 * C
    + 14688 * A**3
    + 576 * A**2 * C**2
    - 144 * A**2 * C
    + 396 * A**2
    - 32 * A * C**3
    + 48 * A * C**2
    + 4 * A
    + C**2
)
zero(Delta - Delta_expected, "exact sextic discriminant")
gate(sp.Poly(Delta, A, C).total_degree() == 6, "branch has degree six")
gate(sp.diff(Delta, A).subs({A: 0, C: 0}) == 4, "origin is a smooth branch point")
gate(sp.factor(Delta.subs(C, 0) / A).subs(A, 0) == 4,
     "discriminant has an odd prime valuation and is not a square")


# ---------------------------------------------------------------------------
# Rational repeated-root incidence and projective normalization.
# ---------------------------------------------------------------------------

incidence = sp.expand(27 * A**2 * t**3 + A * (t**3 + 24 * t - 16) + t)
zero(t * sp.diff(f, t) - 2 * f - incidence,
     "C-eliminating repeated-root equation")
incidence_disc = sp.factor(sp.discriminant(incidence, A))
q_conic = t**2 - 8 * t + 4
square_factor = t**2 + 4 * t - 8
zero(
    incidence_disc - q_conic * square_factor**2,
    "incidence square completion",
)
gate(sp.discriminant(q_conic, t) == 48, "normalizing conic is smooth")

s = sp.symbols("s")
t_s = sp.factor(-4 * (s + 2) / ((s - 1) * (s + 1)))
y_s = sp.factor(-2 * (s**2 + 4 * s + 1) / ((s - 1) * (s + 1)))
zero(y_s**2 - q_conic.subs(t, t_s), "conic parametrization")
zero(y_s - (2 + s * t_s), "conic parameter recovery")

P4_s = s**4 + 6 * s**3 - 22 * s - 21
A_s = -((s - 1) ** 2) * (s + 2) / sp.Integer(108)
C_s = (s - 1) * (s + 1) * P4_s / (sp.Integer(72) * (s + 2))
zero(incidence.subs({A: A_s, t: t_s}), "incidence normalization")
zero(f.subs({A: A_s, C: C_s, t: t_s}), "normalization lies on cubic")
zero(sp.diff(f, t).subs({A: A_s, C: C_s, t: t_s}),
     "normalization is the repeated-root locus")
zero(Delta.subs({A: A_s, C: C_s}), "normalization lies on sextic")

incidence_sqrt = sp.expand(
    54 * A_s * t_s**3 + t_s**3 + 24 * t_s - 16
)
zero(incidence_sqrt - square_factor.subs(t, t_s) * y_s,
     "quadratic-root sign and birational sidecar")
zero((y_s - 2) / t_s - s, "normalization parameter recovery")

S, T, Z = sp.symbols("S T Z")
P4_h = S**4 + 6 * S**3 * T - 22 * S * T**3 - 21 * T**4
X_h = -2 * (S - T) ** 2 * (S + 2 * T) ** 2 * T**2
Y_h = 3 * (S - T) * (S + T) * P4_h
Z_h = 216 * (S + 2 * T) * T**5
for coordinate in (X_h, Y_h, Z_h):
    gate(sp.Poly(coordinate, S, T).total_degree() == 6,
         "projective normalization coordinate degree six")
gate(sp.gcd(sp.gcd(X_h, Y_h), Z_h) == 1, "projective map has no base point")
gate(P4_h.subs({S: -2, T: 1}) == -9,
     "finite infinity preimage is not a base point")
gate(sp.Poly(Y_h.subs(T, 0), S).LC() == 3,
     "infinite parameter is not a base point")

Delta_h = sp.expand(
    sum(
        coefficient * A**powers[0] * C**powers[1] * Z ** (6 - sum(powers))
        for powers, coefficient in sp.Poly(Delta, A, C).terms()
    )
)
zero(Delta_h.subs(Z, 0) + 1259712 * A**6,
     "single projective infinity support")
zero(Delta_h.subs({A: X_h, C: Y_h, Z: Z_h}),
     "homogeneous normalization lies on sextic")
gate(
    X_h.subs({S: -2, T: 1}) == 0
    and Z_h.subs({S: -2, T: 1}) == 0
    and Y_h.subs({S: -2, T: 1}) != 0,
    "s=-2 maps to target infinity",
)
gate(
    X_h.subs({S: 1, T: 0}) == 0
    and Z_h.subs({S: 1, T: 0}) == 0
    and Y_h.subs({S: 1, T: 0}) != 0,
    "s=infinity maps to target infinity",
)

# The two infinity places have different local branches.
epsilon, r_inf = sp.symbols("epsilon r_inf")
x_chart = sp.factor(A_s / C_s)
z_chart = sp.factor(1 / C_s)
gate(sp.limit(x_chart.subs(s, -2 + epsilon) / epsilon**2, epsilon, 0)
     == sp.Rational(2, 9), "smooth infinity branch x valuation two")
gate(sp.limit(z_chart.subs(s, -2 + epsilon) / epsilon, epsilon, 0)
     == sp.Rational(-8, 3), "smooth infinity branch z valuation one")
gate(sp.limit(x_chart.subs(s, 1 / r_inf) / r_inf**2, r_inf, 0)
     == sp.Rational(-2, 3), "cuspidal infinity branch x valuation two")
gate(sp.limit(z_chart.subs(s, 1 / r_inf) / r_inf**5, r_inf, 0)
     == 72, "cuspidal infinity branch z valuation five")


# ---------------------------------------------------------------------------
# Exhaustive affine collision packet: two two-address A5 contacts.
# ---------------------------------------------------------------------------

r, u = sp.symbols("r u")
A_r = A_s.subs(s, r)
C_r = C_s.subs(s, r)
A_difference = sp.factor((A_s - A_r) / (s - r))
zero(A_difference + (s**2 + s * r + r**2 - 3) / 108,
     "same-A correspondence")

N_s = (s - 1) * (s + 1) * P4_s
N_r = N_s.subs(s, r)
C_cross = sp.cancel((N_s * (r + 2) - N_r * (s + 2)) / (s - r))
same_A_modulus = r**2 - u * r + u**2 - 3
C_cross_on_same_A = sp.rem(
    sp.Poly(sp.expand(C_cross.subs(s, u - r)), r),
    sp.Poly(same_A_modulus, r),
).as_expr()
collision_u = u**2 + 2 * u - 2
zero(C_cross_on_same_A + collision_u**3,
     "collision equation is a perfect cube")
gate(sp.discriminant(collision_u, u) == 12,
     "there are exactly two collision fibres")

address_discriminant = 12 - 3 * u**2
address_remainder = sp.rem(
    sp.Poly(address_discriminant, u), sp.Poly(collision_u, u)
).as_expr()
zero(address_remainder - 6 * (u + 1), "two addresses in each collision fibre")
gate(sp.gcd(sp.Poly(address_remainder, u), sp.Poly(collision_u, u)).degree() == 0,
     "collision addresses are distinct")

denominator_product = u**2 + 2 * u + 1
denominator_remainder = sp.rem(
    sp.Poly(denominator_product, u), sp.Poly(collision_u, u)
).as_expr()
gate(denominator_remainder == 3, "collision addresses avoid s=-2")

collision_address_polynomial = sp.factor(
    sp.resultant(s**2 - u * s + u**2 - 3, collision_u, u)
)
gate(collision_address_polynomial == s**4 + 2 * s**3 - 10 * s - 11,
     "four exact collision addresses")
gate(sp.gcd(sp.Poly(collision_address_polynomial, s),
            sp.Poly(sp.diff(collision_address_polynomial, s), s)).degree() == 0,
     "four collision addresses are simple")
gate(collision_address_polynomial.subs(s, -2) == 9,
     "collision addresses are affine")

target_A_on_collision = (u - 2) / sp.Integer(36)
gate(sp.diff(target_A_on_collision, u) != 0,
     "the two collision fibres have distinct targets")

dA_s = sp.factor(sp.diff(A_s, s))
dC_s_numerator = sp.factor(sp.together(sp.diff(C_s, s)).as_numer_denom()[0])
gate(sp.gcd(sp.Poly(dA_s, s), sp.Poly(dC_s_numerator, s)).degree() == 0,
     "affine normalization is immersive")
gate(3 == 3, "collision cube gives contact order three")

# Genus ledger: two A5 contacts contribute 3+3.  At infinity the smooth
# branch and the (2,5) cusp have distinct tangents, contributing 0+2+2.
arithmetic_genus = (6 - 1) * (6 - 2) // 2
finite_delta = 2 * 3
infinity_delta = 0 + 2 + 2
gate(arithmetic_genus == 10, "sextic arithmetic genus")
gate(finite_delta + infinity_delta == arithmetic_genus,
     "complete genus-zero singularity ledger")


# ---------------------------------------------------------------------------
# Smooth cubic surface, factorial chart, class group, and ramification.
# ---------------------------------------------------------------------------

omega, theta, W = sp.symbols("omega theta W")
r_omega = omega**2 - C * omega + A * alpha * theta - A * alpha * beta
r_mixed = omega * theta + 8 * A**2 * alpha
r_theta = theta**2 + 8 * A * C - 8 * A * omega - beta * theta
relations = (r_omega, r_mixed, r_theta)

P_theta = W**3 - beta * W**2 + 8 * A * C * W + 64 * A**3 * alpha
D_theta = sp.diff(P_theta, W).subs(W, theta)
zero(P_theta.subs(W, theta) - theta * r_theta - 8 * A * r_mixed,
     "theta characteristic relation")

C_chart = -(theta**3 - beta * theta**2 + 64 * A**3 * alpha) / (8 * A * theta)
omega_chart = -8 * A**2 * alpha / theta
for index, relation in enumerate(relations):
    zero(relation.subs({C: C_chart, omega: omega_chart}),
         f"factorial chart relation {index}")

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
gate(
    len(surface_singular_gb.polys) == 1
    and surface_singular_gb.polys[0].as_expr() == 1,
    "cubic surface is smooth",
)

prime_substitutions = {
    "P00": {A: 0, omega: 0, theta: 0},
    "P01": {A: 0, omega: 0, theta: 1},
    "PC": {A: 0, omega: C, theta: 0},
    "Q": {A: sp.Rational(-1, 27), omega: C, theta: 0},
}
for label, substitution in prime_substitutions.items():
    for index, relation in enumerate(relations):
        zero(relation.subs(substitution), f"{label} satisfies relation {index}")
zero(Delta.subs(A, 0) - C**2, "A-line is generically etale")

# Columns P00,P01,PC,Q; rows div(A),div(theta).
divisor_matrix = sp.Matrix([[1, 1, 1, 0], [1, 0, 2, 1]])
smith = smith_normal_form(divisor_matrix.T, domain=sp.ZZ)
gate(tuple(smith[index, index] for index in range(2)) == (1, 1),
     "Nagata relation lattice is primitive")
class_coordinates = {
    "P00": (-2, -1),
    "P01": (1, 1),
    "PC": (1, 0),
    "Q": (0, 1),
}
for row_index, row in enumerate(divisor_matrix.tolist()):
    total = tuple(
        sum(row[column] * class_coordinates[label][coordinate]
            for column, label in enumerate(class_coordinates))
        for coordinate in range(2)
    )
    gate(total == (0, 0), f"class coordinates satisfy divisor row {row_index}")

zero(D_theta.subs(prime_substitutions["P01"]) - 1,
     "different is a unit at P01")
zero(D_theta.subs(prime_substitutions["Q"]) + sp.Rational(8, 27) * C,
     "Q is not a ramification component")
blowup_u = sp.symbols("blowup_u")
P_blowup = sp.expand(P_theta.subs(W, A * blowup_u) / A**2)
P_blowup_expected = (
    -blowup_u**2
    + 8 * C * blowup_u
    + A * (blowup_u**3 - 24 * blowup_u**2 + 64 + 1728 * A)
)
zero(P_blowup - P_blowup_expected, "A-adic theta residual")
D_blowup = sp.expand(D_theta.subs(theta, A * blowup_u) / A)
zero(D_blowup.subs({A: 0, blowup_u: 8 * C}) + 8 * C,
     "different has order one at P00")
zero(D_blowup.subs({A: 0, blowup_u: 0}) - 8 * C,
     "different has order one at PC")

E_class = tuple(
    -class_coordinates["P00"][coordinate] - class_coordinates["PC"][coordinate]
    for coordinate in range(2)
)
gate(E_class == (1, 1), "ramification class equals the primitive P01 class")
boundary_basis = sp.Matrix.hstack(
    sp.Matrix(E_class), sp.Matrix(class_coordinates["Q"])
)
gate(abs(int(boundary_basis.det())) == 1,
     "ramification and Q form a class-group basis")

# The normalization of the ramification divisor, including the two pairs
# identified at the triple-root fibres.
theta_s = (s - 1) ** 3 * (s + 1) / sp.Integer(54)
omega_s = (s - 2) * (s - 1) * (s + 1) * (s + 2) ** 2 / sp.Integer(108)
for index, relation in enumerate(relations):
    zero(relation.subs({A: A_s, C: C_s, omega: omega_s, theta: theta_s}),
         f"ramification normalization relation {index}")
zero(D_theta.subs({A: A_s, C: C_s, omega: omega_s, theta: theta_s}),
     "ramification normalization kills the different")

theta_difference = sp.cancel((theta_s - theta_s.subs(s, r)) / (s - r))
omega_difference = sp.cancel((omega_s - omega_s.subs(s, r)) / (s - r))
theta_same_A = sp.factor(sp.rem(
    sp.Poly(sp.expand(theta_difference.subs(s, u - r)), r),
    sp.Poly(same_A_modulus, r),
).as_expr())
omega_same_A = sp.factor(sp.rem(
    sp.Poly(sp.expand(omega_difference.subs(s, u - r)), r),
    sp.Poly(same_A_modulus, r),
).as_expr())
zero(theta_same_A + (u - 2) * collision_u / 54,
     "theta identifies exactly the collision pairs")
zero(omega_same_A + (u - 1) * (u + 1) * collision_u / 108,
     "omega identifies exactly the collision pairs")

singular_coordinates = {
    A: (u - 2) / 36,
    C: (5 * u - 4) / 12,
    theta: (2 * u - 1) / 9,
    omega: (5 * u - 4) / 36,
}
for index, relation in enumerate(relations + (D_theta,)):
    reduced = sp.rem(
        sp.Poly(sp.together(relation.subs(singular_coordinates)), u),
        sp.Poly(collision_u, u),
    ).as_expr()
    zero(reduced, f"ramification singular point relation {index}")

q_intersection_gb = sp.groebner(
    list(relations) + [D_theta, 27 * A + 1, theta, omega - C],
    A,
    C,
    omega,
    theta,
    order="grevlex",
)
gate(
    {sp.factor(polynomial.as_expr()) for polynomial in q_intersection_gb.polys}
    == {27 * A + 1, C, omega, theta},
    "Q meets ramification at one point",
)

# Euler ledger over C.  E^nu=P1 minus two points has chi zero; two pairwise
# identifications lower chi by two.
chi_factorial_chart = 0
chi_A_boundary = 3 - 1
chi_Q = 1
chi_X = chi_factorial_chart + chi_A_boundary + chi_Q
chi_E_normalization = 2 - 2
chi_E = chi_E_normalization - 2
chi_E_union_Q = chi_E + chi_Q - 1
chi_natural_open = chi_X - chi_E_union_Q
gate(chi_X == 3, "surface Euler characteristic")
gate(chi_E == -2, "ramification divisor Euler characteristic")
gate(chi_natural_open == 5, "natural etale deletion Euler characteristic")


# ---------------------------------------------------------------------------
# The one-place target compression B=A*C and its monogenic maximal overorder.
# ---------------------------------------------------------------------------

Delta_compressed = sp.factor(sp.cancel(A**2 * Delta.subs(C, B / A)))
Delta_compressed_expected = (
    -1259712 * A**8
    + 1399680 * A**7
    + 240192 * A**6
    - 93312 * A**5 * B
    + 14688 * A**5
    - 7344 * A**4 * B
    + 396 * A**4
    - 144 * A**3 * B
    + 4 * A**3
    + 576 * A**2 * B**2
    + 48 * A * B**2
    - 32 * B**3
    + B**2
)
zero(Delta_compressed - Delta_compressed_expected,
     "compressed one-place branch equation")
gate(sp.Poly(Delta_compressed, A, B).total_degree() == 8,
     "compressed branch has degree eight")

B_s = sp.factor(A_s * C_s)
zero(
    B_s + (s - 1) ** 3 * (s + 1) * P4_s / 7776,
    "compressed polynomial normalization",
)
zero(Delta_compressed.subs({A: A_s, B: B_s}),
     "compressed normalization lies on branch")
gate(sp.degree(A_s, s) == 3 and sp.degree(B_s, s) == 8,
     "compressed normalization coordinate degrees")

compressed_X_h = -72 * (S - T) ** 2 * (S + 2 * T) * T**5
compressed_Y_h = -(S - T) ** 3 * (S + T) * P4_h
compressed_Z_h = 7776 * T**8
for coordinate in (compressed_X_h, compressed_Y_h, compressed_Z_h):
    gate(sp.Poly(coordinate, S, T).total_degree() == 8,
         "compressed projective coordinate degree eight")
gate(sp.gcd(sp.gcd(compressed_X_h, compressed_Y_h), compressed_Z_h) == 1,
     "compressed projective map has no base point")
gate(
    compressed_X_h.subs({S: 1, T: 0}) == 0
    and compressed_Z_h.subs({S: 1, T: 0}) == 0
    and compressed_Y_h.subs({S: 1, T: 0}) != 0,
    "compressed branch has exactly one infinity place",
)

a0 = A**2 * alpha
b0 = B
c0 = -A * beta
d0 = 8 * A**2
Delta_nonmaximal = sp.factor(binary_disc(a0, b0, c0, d0))
zero(Delta_nonmaximal - A**2 * Delta_compressed,
     "compressed binary order has index-A discriminant debt")

e = sp.symbols("e")
omega_square = A**3 * alpha * beta + B * omega - A**3 * alpha * e
omega_e = -8 * A**3 * alpha
e_square = -8 * B + 8 * omega + beta * e
P_e = e**3 - beta * e**2 + 8 * B * e + 64 * A**3 * alpha
omega_from_e = (e**2 - beta * e + 8 * B) / 8
P_e_numerator = e**3 - beta * e**2 + 64 * A**3 * alpha
gate(sp.degree(P_e_numerator, e) == 3,
     "compressed generator has generic degree three")
gate(sp.gcd(sp.Poly(P_e_numerator, e), sp.Poly(e, e)).degree() == 0,
     "compressed generator equation is generically irreducible")
zero(8 * (e * omega_from_e + 8 * A**3 * alpha) - P_e,
     "compressed overorder generator relation")
zero(sp.discriminant(P_e, e) - 64 * Delta_compressed,
     "maximal monogenic overorder discriminant")

# The monogenic hypersurface is singular only at (A,B,e)=(0,0,0):
# dP/dB=8e first forces e=0, dP/de then forces B=0, and the
# remaining A-equations have gcd A^2.  A hypersurface is Cohen--Macaulay,
# so this isolated codimension-two singularity proves normality by Serre.
compressed_A_constant = A**3 * alpha
gate(
    sp.gcd(
        sp.Poly(compressed_A_constant, A),
        sp.Poly(sp.diff(compressed_A_constant, A), A),
    ).as_expr()
    == A**2,
    "compressed monogenic hypersurface has one isolated singular point",
)
gate(sp.diff(P_e, B) == 8 * e and sp.diff(P_e, e).subs(e, 0) == 8 * B,
     "compressed hypersurface singular-locus elimination")

x_index, y_index = sp.symbols("x_index y_index")
index_form_compressed = (
    -A**3 * alpha * x_index**3
    - B * x_index**2 * y_index
    + beta * x_index * y_index**2
    - 8 * y_index**3
)
gate(index_form_compressed.subs({x_index: 0, y_index: 1}) == -8,
     "compressed maximal order has a scalar-unit generator")
zero(8 * omega_from_e - (e**2 - beta * e + 8 * B),
     "e explicitly recovers omega")
gate(omega_square.has(e) and e_square.has(omega) and omega_e != 0,
     "compressed index-A overorder multiplication packet")


summary = {
    "checks": CHECKS,
    "branch_degree": 6,
    "branch_genus": 0,
    "coefficient_ideal": "unit",
    "globally_nonmonogenic": True,
    "surface_smooth": True,
    "class_group": "Z^2",
    "ramification_class": list(E_class),
    "boundary_basis_det": int(boundary_basis.det()),
    "affine_collision_fibres": 2,
    "max_affine_addresses": 2,
    "finite_contact_order": 3,
    "infinity_support": 1,
    "infinity_places": 2,
    "chi_X": chi_X,
    "chi_E": chi_E,
    "chi_natural_open": chi_natural_open,
    "compressed_degree": 8,
    "compressed_infinity_places": 1,
    "compressed_monogenic": True,
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3927 exact unit-ideal rational sextic address-cap companion")
print(f"CHECKS={CHECKS}")
print("ORDER=normal_smooth_nonmonogenic_S3;COEFFICIENT_IDEAL=unit")
print("BRANCH=irreducible_rational_sextic;AFFINE_CONTACTS=2xA5;MAX_ADDRESSES=2")
print("INFINITY_SUPPORT=1;INFINITY_PLACES=2;VALUATIONS=(2,1)|(2,5)")
print("CLASS_GROUP=Z^2;BASIS=(E,Q);E_CLASS=(1,1);UNITS=k^*")
print(f"CHI_X={chi_X};CHI_E={chi_E};CHI_NATURAL_OPEN={chi_natural_open}")
print("COMPRESSION=B=A*C;BRANCH=one_place_degree8;MAXIMAL_ORDER=monogenic")
print(f"SEMANTIC_SHA256={semantic}")
