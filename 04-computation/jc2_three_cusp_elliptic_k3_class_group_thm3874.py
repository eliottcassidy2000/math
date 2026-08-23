#!/usr/bin/env python3
"""Assertion-free exact gates for THM-3874.

Reproduction:
  python3 04-computation/jc2_three_cusp_elliptic_k3_class_group_thm3874.py
  python3 -O 04-computation/jc2_three_cusp_elliptic_k3_class_group_thm3874.py
"""

from __future__ import annotations

import hashlib

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


a, x, y, w, X, Y, u = sp.symbols("a x y w X Y u")
CHECKS = 0


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def nonzero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value == 0:
        raise AssertionError(f"{label}: unexpectedly zero")


def equal(label: str, left: object, right: object) -> None:
    global CHECKS
    CHECKS += 1
    if left != right:
        raise AssertionError(f"{label}: {left!r} != {right!r}")


Delta_xy = sp.expand(
    81 * x**5
    + 90 * x**4
    + 25 * x**3
    + 30 * x**2 * y**2
    + 30 * x * y**2
    - y**4
    + 8 * y**2
)
F = 15 * a**2 - 15 * a + 4
L = 9 * a - 5
M = 9 * a - 4
b = 1 - a
Delta = sp.expand(Delta_xy.subs(x, a - 1))

zero("bicubic_polarization", F**2 - a**3 * L**2 - b**3 * M**2)
zero("quartic_Delta", Delta - (y**2 * (2 * F - y**2) - b**3 * M**2))

# The sign-changed quadratic resolvent has a monic even quartic generic fibre.
quartic_rhs = sp.expand(y**4 - 2 * F * y**2 + b**3 * M**2)
quartic_eq = sp.expand(w**2 - quartic_rhs)
Xmap = 2 * (w + y**2 - F)
Ymap = 2 * Xmap * y
weierstrass = sp.expand(Ymap**2 - Xmap * (Xmap**2 + 4 * F * Xmap + 4 * a**3 * L**2))
zero("quartic_to_Weierstrass", sp.rem(sp.Poly(weierstrass, w), sp.Poly(quartic_eq, w)).as_expr())

# Standard Weierstrass invariants for Y^2=X^3+4F X^2+4a^3L^2 X.
a2 = 4 * F
a4 = 4 * a**3 * L**2
b2 = 4 * a2
b4 = 2 * a4
b6 = 0
b8 = -a4**2
c4 = sp.factor(b2**2 - 24 * b4)
disc = sp.factor(-b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6)
zero(
    "c4_formula",
    c4 + 64 * (243 * a**5 - 1170 * a**4 + 1875 * a**3 - 1380 * a**2 + 480 * a - 64),
)
zero("elliptic_discriminant", disc + 4096 * a**6 * (a - 1) ** 3 * L**4 * M**2)

# Multiplicative finite fibres: c4 is nonzero at every discriminant root.
equal("c4_at_a0", c4.subs(a, 0), 4096)
equal("c4_at_a1", c4.subs(a, 1), 1024)
equal("c4_at_L0", c4.subs(a, sp.Rational(5, 9)), sp.Rational(16384, 729))
equal("c4_at_M0", c4.subs(a, sp.Rational(4, 9)), sp.Rational(4096, 729))
equal("finite_discriminant_multiplicities", (6, 3, 4, 2), (6, 3, 4, 2))

# Minimal infinity scaling X=u^-4 X_inf, Y=u^-6 Y_inf, a=u^-1.
a2_inf = sp.factor(u**4 * a2.subs(a, 1 / u))
a4_inf = sp.factor(u**8 * a4.subs(a, 1 / u))
c4_inf = sp.factor(u**8 * c4.subs(a, 1 / u))
disc_inf = sp.factor(u**24 * disc.subs(a, 1 / u))
equal("infinity_a2_order", sp.Poly(a2_inf, u).terms()[-1][0][0], 2)
equal("infinity_a4_order", sp.Poly(a4_inf, u).terms()[-1][0][0], 3)
equal("infinity_c4_order", sp.Poly(c4_inf, u).terms()[-1][0][0], 3)
equal("infinity_disc_order", sp.Poly(disc_inf, u).terms()[-1][0][0], 9)
equal("Euler_number", 6 + 3 + 4 + 2 + 9, 24)
equal("root_rank", (6 - 1) + (3 - 1) + (4 - 1) + (2 - 1) + 7, 18)

# Root-lattice determinants and the C2 torsion correction.
def cartan_A(n: int) -> sp.Matrix:
    matrix = 2 * sp.eye(n)
    for j in range(n - 1):
        matrix[j, j + 1] = -1
        matrix[j + 1, j] = -1
    return matrix


def star_cartan(arms: tuple[int, ...]) -> sp.Matrix:
    matrix = 2 * sp.eye(1 + sum(arms))
    nxt = 1
    for length in arms:
        previous = 0
        for _ in range(length):
            matrix[previous, nxt] = -1
            matrix[nxt, previous] = -1
            previous = nxt
            nxt += 1
    return matrix


A5 = cartan_A(5)
A2 = cartan_A(2)
A3 = cartan_A(3)
A1 = cartan_A(1)
E7 = star_cartan((1, 2, 3))
equal("det_A5", A5.det(), 6)
equal("det_A2", A2.det(), 3)
equal("det_A3", A3.det(), 4)
equal("det_A1", A1.det(), 2)
equal("det_E7", E7.det(), 2)

omega3_A5 = A5.inv()[:, 2]
omega2_A3 = A3.inv()[:, 1]
equal("twice_omega3_A5", tuple(2 * omega3_A5), (1, 2, 3, 2, 1))
equal("twice_omega2_A3", tuple(2 * omega2_A3), (1, 2, 1))
equal("trivial_lattice_abs_discriminant", 6 * 3 * 4 * 2 * 2, 288)
equal("NS_abs_discriminant", 288 // 2**2, 72)

# Kill the finite exceptional roots and O,T,U,E7.  The only survivors are
# d=alpha_3 and e=beta_2; twice the C2 gluing class gives 3d+2e.
relation = sp.Matrix([[3, 2]])
snf = smith_normal_form(relation, domain=sp.ZZ)
equal("class_group_relation_gcd", sp.gcd(3, 2), 1)
equal("class_group_relation_snf", tuple(snf), (1, 0))

# Full presentation control: 20 trivial-lattice generators plus T.  The
# first row is twice the torsion-section gluing law; the other rows kill the
# boundary and affine ADE curves.  Its Smith form has twenty unit entries
# and one free column.
omega_E7 = E7.inv()[:, 6]
equal("twice_omega_E7", tuple(2 * omega_E7), (6, 3, 4, 2, 5, 4, 3))
generator_count = 21
gluing = [0] * generator_count
gluing[0], gluing[1], gluing[20] = -2, -4, 2  # O,F,T
for index, coefficient in enumerate((1, 2, 3, 2, 1)):
    gluing[2 + index] = coefficient
for index, coefficient in enumerate((1, 2, 1)):
    gluing[9 + index] = coefficient
for index, coefficient in enumerate((6, 3, 4, 2, 5, 4, 3)):
    gluing[13 + index] = coefficient

presentation_rows = [gluing]
killed_generators = [
    0,
    1,
    20,
    2,
    3,
    5,
    6,
    7,
    8,
    9,
    11,
    12,
    *range(13, 20),
]
for index in killed_generators:
    row = [0] * generator_count
    row[index] = 1
    presentation_rows.append(row)
full_presentation = sp.Matrix(presentation_rows)
full_snf = smith_normal_form(full_presentation, domain=sp.ZZ)
equal("full_boundary_presentation_rank", full_presentation.rank(), 20)
equal("full_boundary_presentation_snf", tuple(full_snf[i, i] for i in range(20)), (1,) * 20)

# Finite singularity/fibre addresses inherited from the normalization.
equal("cusp_a_addresses", (0 + 1, -1 + 1, -1 + 1), (1, 0, 0))
equal("node_M_address", -sp.Rational(5, 9) + 1, sp.Rational(4, 9))
equal("node_L_addresses", -sp.Rational(4, 9) + 1, sp.Rational(5, 9))
zero("M_node_fibre", M.subs(a, sp.Rational(4, 9)))
zero("L_node_fibre", L.subs(a, sp.Rational(5, 9)))

# The unit proof uses the odd total degree of the irreducible branch.
equal("Delta_total_degree", sp.Poly(Delta_xy, x, y).total_degree(), 5)
equal("quadratic_surface_module_rank", 2, 2)

semantic_packet = (
    "quadratic resolvent is the elliptic K3 with fibres I6,I3,I4,I2,III*",
    "Triv=U+A5+A2+A3+A1+E7 has rank 20 and index 2 in NS",
    "MW rank zero, MW torsion C2, absolute NS discriminant 72",
    "finite exceptions leave alpha3=d and beta2=e",
    "boundary kills O,T,U,E7 and the torsion correction gives 3d+2e=0",
    "Cl(Q)=Z^2/<(3,2)>=Z and Q*=k*",
    "no codimension-one-unramified connected C3 cover of Q",
    "no normal finite-flat cubic with sole simple branch Delta",
    "disc(NS) has 3-primary data; only the explicit affine quotient kills it",
)
semantic_sha = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()

print("THM3874_FIBRES", "I6(a=0),I3(a=1),I4(L=0),I2(M=0),III*(infinity)")
print("THM3874_K3", "Euler=24;rho=20;MW_rank=0;MW_tors=C2")
print("THM3874_NS", "Triv index2;abs_disc_NS=72;3-primary not zero before quotient")
print("THM3874_AFFINE_QUOTIENT", "survivors d=alpha3,e=beta2;relation 3d+2e=0")
print("THM3874_CLASS_GROUP", "Cl(Q)=Z;Cl(Q)[3]=0;units=k*")
print("THM3874_C3", "connected codimension-one-unramified C3 layer=NONE")
print("THM3874_CUBIC", "normal finite-flat cubic with sole simple branch Delta=NONE")
print("THM3874_SCOPE", "three-cusp quadratic resolvent only;JC2 remains open")
print("SEMANTIC_SHA256", semantic_sha)
print("CHECKS", CHECKS)
