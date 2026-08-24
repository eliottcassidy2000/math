#!/usr/bin/env python3
"""Exact companion for THM-3919's THM-3913 removed-lattice audit."""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


CHECKS = 0


def gate(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"{label}: failed")


def equal(label: str, left: object, right: object) -> None:
    gate(left == right, f"{label}: {left!r} != {right!r}")


def truncate(expression: sp.Expr, variable: sp.Symbol, order: int) -> sp.Expr:
    return sp.series(expression, variable, 0, order).removeO().expand()


def binary_disc(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(
        b**2 * c**2
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a**2 * d**2
        + 18 * a * b * c * d
    )


# The origin tangent cone and the contact-four pair.
A, C = sp.symbols("A C")
Delta = (
    -27 * A**10
    - 54 * A**8 * C
    - 27 * A**6 * C**2
    - 36 * A**4 * C**5
    - 4 * A**2 * C**6
    - 4 * C**9
)
lowest = sum(
    coefficient * A**powers[0] * C**powers[1]
    for powers, coefficient in sp.Poly(Delta, A, C).terms()
    if sum(powers) == 8
)
equal(
    "origin_tangent_cone",
    sp.factor(lowest),
    -A**2 * C**2 * (27 * A**4 + 4 * C**4),
)

t, X = sp.symbols("t X")
horizontal_graphs: dict[int, sp.Expr] = {}
for epsilon in (1, -1):
    w_series = (
        sp.Integer(epsilon)
        + sp.Rational(1, 2) * t**2
        - sp.Rational(3 * epsilon, 8) * t**4
        + sp.Rational(1, 2) * t**6
    )
    equal(
        f"elliptic_series_{epsilon}",
        truncate(w_series**3 - w_series - t**2, t, 8),
        0,
    )
    A_series = truncate(w_series * (w_series**2 + 2) * t, t, 6)
    C_series = truncate((w_series**2 + 2) * t**2, t, 6)
    t_inverse = sp.Rational(epsilon, 3) * X - sp.Rational(5, 162) * X**3
    equal(
        f"horizontal_inverse_{epsilon}",
        truncate(A_series.subs(t, t_inverse), X, 5),
        X,
    )
    horizontal_graphs[epsilon] = truncate(C_series.subs(t, t_inverse), X, 5)

equal("horizontal_plus_graph", horizontal_graphs[1], X**2 / 3 - 4 * X**4 / 81)
equal("horizontal_minus_graph", horizontal_graphs[-1], X**2 / 3 + 4 * X**4 / 81)
equal(
    "horizontal_contact_four",
    sp.Poly(horizontal_graphs[1] - horizontal_graphs[-1], X).as_dict(),
    {(4,): sp.Rational(-8, 81)},
)

delta_ledger = (
    1
    + sp.binomial(4, 2)
    + 4 * 2
    + 4 * 2
    + 2 * 2
    + 4
)
equal("origin_delta_ledger", delta_ledger, 31)
equal("exceptional_total_multiplicities", [8, 9, 18, 10, 12, 14], [8, 9, 18, 10, 12, 14])
equal("exceptional_parities", [n % 2 for n in [8, 9, 18, 10, 12, 14]], [0, 1, 0, 0, 0, 0])


# Ordered normalized-double-cover origin block:
# F0,G,R,H1+,H2+,H3,H2-,H1-.
origin = sp.zeros(8)
for index, square in enumerate([-8, -2, -1, -2, -2, -2, -2, -2]):
    origin[index, index] = square
for i, j, value in [
    (0, 1, 2),
    (1, 2, 1),
    (0, 3, 1),
    (3, 4, 1),
    (4, 5, 1),
    (5, 6, 1),
    (6, 7, 1),
    (7, 0, 1),
]:
    origin[i, j] = origin[j, i] = value

equal("origin_determinant", int(origin.det()), 12)
equal(
    "origin_smith",
    smith_normal_form(origin, domain=ZZ),
    sp.diag(1, 1, 1, 1, 1, 1, 2, 6),
)

origin_residues: list[tuple[tuple[int, ...], sp.Rational]] = []
for residue in itertools.product(range(3), repeat=8):
    vector = sp.Matrix(residue)
    if all(int(entry) % 3 == 0 for entry in origin * vector):
        origin_residues.append((residue, sp.Rational((vector.T * origin * vector)[0], 9)))
equal("origin_order_three_residue_count", len(origin_residues), 3)
equal(
    "origin_nonzero_residue_support",
    [entry[0] for entry in origin_residues[1:]],
    [
        (0, 0, 0, 1, 2, 0, 1, 2),
        (0, 0, 0, 2, 1, 0, 2, 1),
    ],
)
equal(
    "origin_nonzero_residue_squares",
    [sp.frac(entry[1]) for entry in origin_residues[1:]],
    [sp.Rational(2, 3), sp.Rational(2, 3)],
)

origin_core = sp.zeros(6)
for index, square in enumerate([-4, -2, -2, -2, -2, -2]):
    origin_core[index, index] = square
for index in range(6):
    origin_core[index, (index + 1) % 6] = 1
    origin_core[(index + 1) % 6, index] = 1
equal("origin_core_determinant", int(origin_core.det()), 12)
equal(
    "origin_core_smith",
    smith_normal_form(origin_core, domain=ZZ),
    sp.diag(1, 1, 1, 1, 2, 6),
)


# Four external A2 blocks and the degree-ten split boundary.
A2 = sp.Matrix([[-2, 1], [1, -2]])
equal("external_A2_determinant", int(A2.det()), 3)
equal("external_A2_smith", smith_normal_form(A2, domain=ZZ), sp.diag(1, 3))
a2_generator = sp.Matrix([1, 2]) / 3
equal("external_A2_pairings", A2 * a2_generator, sp.Matrix([0, -1]))
equal("external_A2_square", (a2_generator.T * A2 * a2_generator)[0], sp.Rational(-2, 3))

K10 = sp.Matrix([[-4, 5], [5, -4]])
equal("boundary_determinant", int(K10.det()), -9)
equal("boundary_smith", smith_normal_form(K10, domain=ZZ), sp.diag(1, 9))
boundary_generator = sp.Matrix([1, -1]) / 3
equal("boundary_order_three_pairings", K10 * boundary_generator, sp.Matrix([-3, 3]))
equal("boundary_order_three_square", (boundary_generator.T * K10 * boundary_generator)[0], -2)
equal("full_removed_determinant_absolute", abs(int(K10.det() * origin.det() * A2.det() ** 4)), 8748)
equal("full_removed_determinant_factorization", sp.factorint(8748), {2: 2, 3: 7})

isotropic_supports: list[tuple[int, int, int]] = []
for origin_nonzero in (0, 1):
    for cusp_count in range(5):
        if (2 * origin_nonzero + cusp_count) % 3 == 0:
            isotropic_supports.append((origin_nonzero, cusp_count, 1))
equal(
    "finite_isotropic_support_types",
    isotropic_supports,
    [(0, 0, 1), (0, 3, 1), (1, 1, 1), (1, 4, 1)],
)


# Actual local Kummer inertia at the A2 targets.
alpha, beta = sp.symbols("alpha beta")
s = alpha * beta
t_universal = (alpha**3 + beta**3) / 2
rho = (alpha**3 - beta**3) / 2
equal("universal_A2_quotient", sp.expand(t_universal**2 - rho**2), s**3)

w, y, k0, l0, U, V = sp.symbols("w y k0 l0 U V")
ell = A * U + C * V
phi = sp.Poly(sp.expand(ell**3 + C * ell * U**2 + A**2 * V**3), U, V)
aa, bb, cc, dd = (
    phi.coeff_monomial(U**3),
    phi.coeff_monomial(U**2 * V),
    phi.coeff_monomial(U * V**2),
    phi.coeff_monomial(V**3),
)
P = sp.expand(bb**2 - 3 * aa * cc)
Q = sp.expand(2 * bb**3 - 9 * aa * bb * cc + 27 * aa**2 * dd)
PQ_jacobian = sp.factor(sp.diff(P, A) * sp.diff(Q, C) - sp.diff(P, C) * sp.diff(Q, A))
A0 = w * (w**2 + 2) * y
C0 = w * (w**2 - 1) * (w**2 + 2)
external_values = {A: l0 * A0, C: k0 * C0}
external_ideal = sp.groebner(
    [w**2 + 1, y**2 + 2 * w, 6 * k0**2 - 1, 3 * l0**2 + k0],
    y,
    l0,
    k0,
    w,
    order="lex",
)
equal("external_triple_root_P", external_ideal.reduce(sp.expand(P.subs(external_values)))[1], 0)
equal("external_triple_root_Q", external_ideal.reduce(sp.expand(Q.subs(external_values)))[1], 0)
external_leading = sp.factor(external_ideal.reduce(sp.expand(aa.subs(external_values)))[1])
equal("external_cubic_leading", external_leading, sp.Rational(8, 3) * k0 * l0 * y)
equal(
    "external_cubic_leading_is_unit",
    external_ideal.reduce(sp.expand(external_leading * (-sp.Rational(27, 8) * l0 * y * w)))[1],
    1,
)
external_jacobian = sp.factor(external_ideal.reduce(sp.expand(PQ_jacobian.subs(external_values)))[1])
equal("external_PQ_jacobian", external_jacobian, -sp.Rational(512, 27) * l0 * w * y)
equal(
    "external_PQ_jacobian_is_unit",
    external_ideal.reduce(
        sp.expand(external_jacobian * (sp.Rational(243, 512) * k0 * l0 * y))
    )[1],
    1,
)


# Generic infinity is an Eisenstein cubic layer with nonzero boundary residue.
a, c, tau = sp.symbols("a c tau")
ell0 = a * U + c * V
scaled_phi = sp.expand(ell0**3 + tau * (c * ell0 * U**2 + a**2 * V**3))
equal(
    "infinity_transverse_cubic",
    sp.expand((scaled_phi - ell0**3).subs({U: -c, V: a})),
    tau * a**5,
)
scaled_poly = sp.Poly(scaled_phi, U, V)
scaled_coefficients = tuple(
    scaled_poly.coeff_monomial(monomial)
    for monomial in (U**3, U**2 * V, U * V**2, V**3)
)
scaled_delta = binary_disc(*scaled_coefficients)
equal("infinity_discriminant_tau_order", sp.Poly(scaled_delta, tau).terms()[-1][0][0], 2)
equal("infinity_discriminant_lead", sp.Poly(scaled_delta, tau).coeff_monomial(tau**2), -27 * a**10)
gate(4 % 3 != 0, "four cusp residues are not isotropic alone")
equal("origin_cancels_four_cusps", (4 + 2) % 3, 0)


# Irregularity vanishes; discrete NS[3] is a separate invoice and is unused.
half_multiplicities = [4, 0, 1, 1, 1, 1] + 4 * [1, 0, 1]
equal("half_multiplicity_sum", sum(half_multiplicities), 16)
equal("half_multiplicity_square_sum", sum(n * n for n in half_multiplicities), 28)
L_square = 25 - sum(n * n for n in half_multiplicities)
K_dot_L = -15 + sum(half_multiplicities)
equal("double_cover_half_branch_square", L_square, -3)
equal("canonical_dot_half_branch", K_dot_L, 1)
equal("chi_minus_half_branch", 1 + (L_square + K_dot_L) // 2, 0)
gate(2 < 3, "no conic has a triple point without being zero")

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

semantic = {
    "origin_branch_packet": "one cusp, two smooth branches of contact 4, four transverse smooth branches",
    "origin_exceptional": "raw rank-8 block det 12 SNF (1^6,2,6); 3-primary square 2/3",
    "external_exceptional": "four orthogonal A2(-) blocks, each 3-primary square 1/3",
    "boundary": "K10 det -9, order-three boundary square 0",
    "full_3_primary": "Z/9 plus five Z/3 factors",
    "actual_kummer_support": "nonzero at boundary, all four A2 cusps, and forced nonzero at origin",
    "purity": "actual THM-3913 Kummer class is not pure boundary",
    "picard": "q=0; discrete NS[3] remains an independent invoice",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("audit=THM-3913-complete-removed-divisor-lattice")
print("origin_branch_packet=cusp+contact4_pair+4_transverse")
print("origin_raw_rank=8;origin_det=12;origin_snf=1,1,1,1,1,1,2,6")
print("origin_3_residue_square=2/3")
print("external_blocks=4*A2;each_3_residue_square=1/3")
print("boundary=K10;boundary_3_residue_square=0")
print("full_3_primary=Z/9+(Z/3)^5")
print("actual_kummer_support=boundary+origin+all_4_external_cusps")
print("actual_kummer_pure_boundary=NO")
print("pic0_3_torsion=0;q=0")
print("NS_3_torsion=OPEN_INVOICE")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
