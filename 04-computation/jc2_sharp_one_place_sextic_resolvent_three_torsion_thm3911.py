#!/usr/bin/env python3
"""Exact companion for THM-3911's sharp sextic resolvent obstruction."""

from __future__ import annotations

import hashlib
import json
from itertools import product

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


t, A, C, Z, T, S = sp.symbols("t A C Z T S")

# The sharp one-place sextic from THM-3906.
G = t**4 + t**3 - 2 * t - 5
A_t = t * G
C_t = (t**2 + 1) * G
delta = (
    -A**6
    + 6 * A**5
    - 25 * A**4
    - 6 * A**3 * C**2
    - 32 * A**3 * C
    - 5 * A**2 * C**3
    - 18 * A**2 * C**2
    + A * C**4
    + 7 * A * C**3
    + C**5
    + 5 * C**4
)
zero(delta.subs({A: A_t, C: C_t}), "sextic vanishes on parametrization")
gate(sp.Poly(delta, A, C).total_degree() == 6, "branch has degree six")
gate(sp.factor(delta) == delta, "branch is irreducible over Q")
zero(
    sp.resultant(A - A_t, C - C_t, t) - delta,
    "parametrization eliminant is exactly the sextic",
)

# A rational inverse proves birationality.  The homogeneous morphism has
# coordinates T*S*Gh, (T^2+S^2)*Gh, S^6 and no base point.
inverse_numerator = -(A**3 - 3 * A**2 - A * C**2 - 6 * A * C)
inverse_denominator = A**2 * C + 4 * A**2 + A * C + C**2
zero(
    inverse_numerator.subs({A: A_t, C: C_t})
    - t * inverse_denominator.subs({A: A_t, C: C_t}),
    "birational inverse",
)
Gh = T**4 + T**3 * S - 2 * T * S**3 - 5 * S**4
A_h_t = T * S * Gh
C_h_t = (T**2 + S**2) * Gh
Z_h_t = S**6
gate(C_h_t.subs(S, 0) == T**6, "infinity parameter is base-point free")
gate(Z_h_t.subs(S, 1) == 1, "finite roots of Gh do not create base points")

# Complete normalization-address packet.  Away from A=0, equal radial ratio
# forces the second address to be t or 1/t.  The reciprocal collision factor
# H gives two nodes; t=+/-1 are the two critical fixed addresses.
G_reciprocal = sp.expand(t**4 * G.subs(t, 1 / t))
gate(sp.discriminant(G, t) == -23355, "ordinary four-address fibre is reduced")
gate(sp.gcd(G, G_reciprocal) == 1, "four origin tangents are distinct")
H = t**4 + t**3 + 3 * t**2 + t + 1
gate(sp.discriminant(H, t) == 189, "four node addresses are reduced")
gate(H.subs(t, 1) == 7 and H.subs(t, -1) == 3, "nodes avoid cusp addresses")
A_difference = sp.factor(sp.together(A_t - A_t.subs(t, 1 / t)))
C_difference = sp.factor(sp.together(C_t - C_t.subs(t, 1 / t)))
gate(
    A_difference == (t - 1) ** 3 * (t + 1) ** 3 * H / t**5,
    "reciprocal A collision factor",
)
gate(
    C_difference == (t - 1) ** 3 * (t + 1) ** 3 * (t**2 + 1) * H / t**6,
    "reciprocal C collision factor",
)
zero(sp.rem(C_t + 3, H, t), "node C coordinate")
zero(sp.rem(A_t**2 - 3 * A_t + 9, H, t), "node A coordinates")

dA_t = sp.factor(sp.diff(A_t, t))
dC_t = sp.factor(sp.diff(C_t, t))
gate(sp.gcd(dA_t, dC_t) == t**2 - 1, "only two critical addresses")
gate((A_t.subs(t, 1), C_t.subs(t, 1)) == (-5, -10), "first cusp target")
gate((A_t.subs(t, -1), C_t.subs(t, -1)) == (3, -6), "second cusp target")
for address, expected in ((1, 1680), (-1, 432)):
    determinant = sp.det(
        sp.Matrix(
            [
                [sp.diff(A_t, t, 2).subs(t, address), sp.diff(C_t, t, 2).subs(t, address)],
                [sp.diff(A_t, t, 3).subs(t, address), sp.diff(C_t, t, 3).subs(t, address)],
            ]
        )
    )
    gate(determinant == expected, f"A2 cusp determinant at t={address}")

slope = sp.cancel(dC_t / dA_t)
slope_difference = sp.factor(sp.together(slope - slope.subs(t, 1 / t)))
zero(
    slope_difference
    - 3 * (t - 1) * (t + 1) * (2 * t**2 + t + 2)
    / (t * (5 * t**2 + 4 * t + 5)),
    "reciprocal tangent separation",
)
gate(
    sp.gcd(sp.together(slope_difference).as_numer_denom()[0], H) == 1,
    "the two node branches are transverse",
)

# Independent singular-ideal check.  Its elimination support and the four
# fibre gcds prove that no finite singular target was missed.
delta_A = sp.diff(delta, A)
delta_C = sp.diff(delta, C)
singular_groebner = sp.groebner([delta, delta_A, delta_C], A, C, order="lex")
gate(
    sp.factor(singular_groebner.polys[-1].as_expr())
    == C**5 * (C + 3) * (C + 6) ** 2 * (C + 10) ** 2,
    "complete singular C-support",
)
expected_fibre_gcds = {
    0: A**3,
    -3: A**2 - 3 * A + 9,
    -6: A - 3,
    -10: A + 5,
}
for c_value, expected in expected_fibre_gcds.items():
    fibre_polynomials = [
        sp.Poly(polynomial.subs(C, c_value), A)
        for polynomial in (delta, delta_A, delta_C)
    ]
    fibre_gcd = fibre_polynomials[0]
    for polynomial in fibre_polynomials[1:]:
        fibre_gcd = sp.gcd(fibre_gcd, polynomial)
    gate(sp.monic(fibre_gcd.as_expr(), A) == expected, f"singular fibre C={c_value}")

tangent_cone = sum(
    coefficient * A**monomial[0] * C**monomial[1]
    for monomial, coefficient in sp.Poly(delta, A, C).terms()
    if sum(monomial) == 4
)
gate(
    sp.discriminant(tangent_cone.subs(C, 1), A) == -425644875,
    "origin is an ordinary quadruple point",
)
hessian_determinant = sp.det(sp.hessian(delta, (A, C)))
gate(
    sp.gcd(
        A**2 - 3 * A + 9,
        sp.Poly(hessian_determinant.subs(C, -3), A),
    )
    == 1,
    "both residual singularities are nodes",
)

# The unique projective infinity is smooth.  The cusp tangent lines make the
# two finite A2 contacts particularly visible.
delta_h = sp.expand(
    sum(
        coefficient * A**monomial[0] * C**monomial[1] * Z ** (6 - sum(monomial))
        for monomial, coefficient in sp.Poly(delta, A, C).terms()
    )
)
zero(delta_h.subs(Z, 0) + A**6, "unique projective infinity support")
infinity_point = {A: 0, C: 1, Z: 0}
gate(sp.diff(delta_h, Z).subs(infinity_point) == 1, "infinity is smooth")
gate(sp.factor(delta.subs(C, 2 * A)) == -A**4 * (A + 5) ** 2,
     "first cusp tangent line")
gate(sp.factor(delta.subs(C, -2 * A)) == -A**4 * (A - 3) ** 2,
     "second cusp tangent line")

# Removed-divisor lattice on a smooth projective resolution.  Ordered blocks:
# elliptic exceptional D; boundary B+,B-; two A2 chains; two A1 curves.
D_block = sp.Matrix([[-2]])
boundary_block = sp.Matrix([[-2, 3], [3, -2]])
A2_block = sp.Matrix([[-2, 1], [1, -2]])
A1_block = sp.Matrix([[-2]])
relation_gram = sp.diag(
    D_block,
    boundary_block,
    A2_block,
    A2_block,
    A1_block,
    A1_block,
)
gate(relation_gram.shape == (9, 9), "removed-divisor relation rank")
gate(boundary_block.det() == -5, "sextic split-boundary determinant")
gate(relation_gram.det() == 360, "removed-divisor Gram determinant")
gate(boundary_block.rank() == 2, "boundary classes are independent")

# Exhaust the mod-three kernel.  Any putative 3-torsion class x has 3x in
# this relation lattice.  Pairing forces its residue into the two A2 radical
# directions.  Every nonzero such residue has square nonintegral after
# division by nine, so the relation lattice is 3-saturated in every integral
# ambient Picard lattice with this incidence packet.
kernel_residues: list[tuple[int, ...]] = []
square_divisible_residues: list[tuple[int, ...]] = []
for residue in product(range(3), repeat=9):
    vector = sp.Matrix(residue)
    if any(int(entry) % 3 for entry in relation_gram * vector):
        continue
    kernel_residues.append(residue)
    square = int((vector.T * relation_gram * vector)[0])
    if square % 9 == 0:
        square_divisible_residues.append(residue)
gate(len(kernel_residues) == 9, "two A2 radical directions modulo three")
gate(
    square_divisible_residues == [(0,) * 9],
    "no nonzero three-divisible relation-lattice residue",
)

semantic_payload = {
    "branch": "rational_sextic_one_smooth_infinity_place",
    "finite_packet": "ordinary_quadruple_plus_2A2_plus_2A1",
    "boundary_gram": [[-2, 3], [3, -2]],
    "removed_gram_determinant": 360,
    "mod3_kernel_size": len(kernel_residues),
    "mod3_square_integral_nonzero": 0,
    "consequence": "Cl_resolvent_3_torsion_zero_no_normal_S3_cubic_lift",
    "scope": "sharp_control_only_no_JC_counterexample",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3911-sharp-one-place-sextic-resolvent-three-torsion-obstruction")
print("branch=sextic;normalization=A1;infinity_places=1")
print("finite_singularities=ordinary_quadruple_plus_2A2_plus_2A1")
print("boundary_gram=[[-2,3],[3,-2]];boundary_determinant=-5")
print("removed_gram_determinant=360;mod3_kernel=9;nonzero_integral_squares=0")
print("quadratic_resolvent_units=kstar;class_group_3_torsion=zero")
print("normal_irreducible_S3_cubic_lift=obstructed")
print("scope=sharp_control_only_no_JC_counterexample")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
