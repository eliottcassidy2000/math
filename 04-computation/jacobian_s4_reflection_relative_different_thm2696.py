#!/usr/bin/env python3
"""Exact checks for THM-2696.

The script verifies the reflection-completed S4 quotient, its graph quartic
and squared-pair resolvent, the intrinsic relative Jacobian/different, a
nontrivial polynomial-generator change, a dominant sandwich control, and the
constant-different planar normalization boundary.  Every truth-bearing check
uses ``require`` so optimized mode retains the evidence.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mmul(left: tuple[tuple[int, int], tuple[int, int]],
         right: tuple[tuple[int, int], tuple[int, int]]) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(2)) % 2 for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def mpow(matrix: tuple[tuple[int, int], tuple[int, int]], exponent: int) -> tuple[tuple[int, int], tuple[int, int]]:
    answer = ((1, 0), (0, 1))
    for _ in range(exponent):
        answer = mmul(answer, matrix)
    return answer


print("THM-2696 REFLECTION-COMPLETED S4 RELATIVE DIFFERENT -- exact checks")

# The C3 rotation and C2 reflection generate GL_2(F_2)=S3 on the three
# nonzero V4 characters.
c3 = ((0, 1), (1, 1))
c2 = ((0, 1), (1, 0))
identity = ((1, 0), (0, 1))
require(mpow(c3, 3) == identity and mpow(c3, 1) != identity, "C3 order failure")
require(mpow(c2, 2) == identity and mpow(c2, 1) != identity, "C2 order failure")
require(mmul(mmul(c2, c3), c2) == mpow(c3, 2), "C2 does not invert C3")
generated = {identity}
changed = True
while changed:
    changed = False
    for old in tuple(generated):
        for generator in (c2, c3):
            new = mmul(old, generator)
            if new not in generated:
                generated.add(new)
                changed = True
require(len(generated) == 6, "C2,C3 do not generate GL2(F2)")
nonzero = ((1, 0), (0, 1), (1, 1))
orbit = {
    ((matrix[0][0] * 1 + matrix[0][1] * 0) % 2,
     (matrix[1][0] * 1 + matrix[1][1] * 0) % 2)
    for matrix in generated
}
require(orbit == set(nonzero), "standard-plane orbit failure")
print("C2_C3_generate_GL2F2_order=6 nonzero_orbit=3: PASS")

# Source roots, quotient coordinates, graph quartic, and squared-pair cubic.
x, y, z, T, U, W = sp.symbols("x y z T U W")
s1_xyz = x + y + z
s2_xyz = x * y + x * z + y * z
s3_xyz = x * y * z
A_xyz = x**2 + y**2 + z**2
B_xyz = x**2 * y**2 + x**2 * z**2 + y**2 * z**2
d_xyz = x * y * z

quartic_roots = (x + y + z, x - y - z, -x + y - z, -x - y + z)
P = sp.expand(sp.prod(T - root for root in quartic_roots))
P_expected = T**4 - 2 * A_xyz * T**2 - 8 * d_xyz * T + A_xyz**2 - 4 * B_xyz
require(sp.expand(P - P_expected) == 0, "graph quartic coefficient failure")

resolvent_roots = (4 * x**2, 4 * y**2, 4 * z**2)
S = sp.expand(sp.prod(U - root for root in resolvent_roots))
S_expected = U**3 - 4 * A_xyz * U**2 + 16 * B_xyz * U - 64 * d_xyz**2
require(sp.expand(S - S_expected) == 0, "squared-pair resolvent failure")

disc_P = sp.factor(sp.discriminant(P, T))
disc_S = sp.factor(sp.discriminant(S, U))
require(sp.expand(disc_P - disc_S) == 0, "quartic/resolvent discriminants differ")
print("graph_quartic_and_squared_pair_resolvent: PASS")
print("disc_graph_quartic_equals_disc_scaled_resolvent: PASS")

# Abstract source quotient coordinates and the relative Jacobian.
s1, s2, s3 = sp.symbols("s1 s2 s3")
A = s1**2 - 2 * s2
B = s2**2 - 2 * s1 * s3
d = s3
F = sp.Matrix((A, B, d))
jac_F = sp.factor(F.jacobian((s1, s2, s3)).det())
J0 = s1 * s2 - s3
require(sp.expand(jac_F - 4 * J0) == 0, "relative Jacobian failure")
require(sp.expand(J0.subs({s1: s1_xyz, s2: s2_xyz, s3: s3_xyz}) -
                          (x + y) * (x + z) * (y + z)) == 0,
        "plus-hyperplane factor failure")

disc_H = sp.discriminant(T**3 - s1 * T**2 + s2 * T - s3, T)
Q = W**3 - A * W**2 + B * W - d**2
disc_Q_pullback = sp.discriminant(Q, W)
require(sp.expand(disc_Q_pullback - disc_H * J0**2) == 0,
        "normalized resolvent discriminant pullback failure")
AA, BB, dd = sp.symbols("AA BB dd")
disc_P_target = 4096 * sp.discriminant(W**3 - AA * W**2 + BB * W - dd**2, W)
require(sp.expand(disc_P_target.subs({AA: A, BB: B, dd: d}) -
                          256 * disc_H * jac_F**2) == 0,
        "scaled discriminant/Jacobian identity failure")
print("relative_Jacobian=4*(s1*s2-s3): PASS")
print("DiscQ_pullback=DiscH*(Jac/4)^2: PASS")
print("DiscP_pullback=256*DiscH*Jac^2: PASS")

# A nontrivial triangular source and target generator change.  Both
# automorphisms have determinant one, so the transported relative different
# must be exactly the pullback of 4*J0.
u1, u2, u3 = sp.symbols("u1 u2 u3")
source_inverse = {s1: u1 - (u2 - u3) ** 2, s2: u2 - u3, s3: u3}
A_u, B_u, d_u = (sp.expand(entry.subs(source_inverse, simultaneous=True)) for entry in F)
F_changed = sp.Matrix((A_u + B_u**2, B_u + d_u, d_u))
jac_changed = sp.factor(F_changed.jacobian((u1, u2, u3)).det())
jac_pullback = sp.factor(jac_F.subs(source_inverse, simultaneous=True))
require(sp.expand(jac_changed - jac_pullback) == 0,
        "polynomial-generator invariance failure")
require(sp.Poly(jac_changed, u1, u2, u3).total_degree() > 0,
        "generator change spuriously made Jacobian constant")
print("triangular_source_target_generator_control_preserves_different: PASS")

# A dominant non-automorphic source sandwich cannot remove the factor either.
r1, r2, r3 = sp.symbols("r1 r2 r3")
source_sandwich = {s1: r1**2 + r2, s2: r2 + r3, s3: r3}
A_r, B_r, d_r = (sp.expand(entry.subs(source_sandwich, simultaneous=True)) for entry in F)
F_sandwich = sp.Matrix((A_r + B_r**2, B_r + d_r, d_r))
jac_sandwich = sp.factor(F_sandwich.jacobian((r1, r2, r3)).det())
expected_sandwich = sp.factor(8 * r1 * J0.subs(source_sandwich, simultaneous=True))
require(sp.expand(jac_sandwich - expected_sandwich) == 0,
        "dominant sandwich factor failure")
print("dominant_polynomial_sandwich_retains_nonconstant_factor: PASS")

# Constant-relative-different planar slice J0=c.  It is A2 on the source,
# but its coefficient image is a singular nonnormal hypersurface whose
# normalization is that A2.
c = sp.symbols("c", nonzero=True)
u, v, At, Bt, dt = sp.symbols("u v A B d")
slice_A = u**2 - 2 * v
slice_B = v**2 - 2 * u**2 * v + 2 * c * u
slice_d = u * v - c
R_c = (
    4 * At**3 * dt**2 - At**2 * Bt**2 - 18 * At * Bt * dt**2
    + 2 * At * Bt * c**2 + 4 * Bt**3 + 27 * dt**4
    - 18 * dt**2 * c**2 + 8 * dt * c**3 - c**4
)
require(sp.expand(R_c.subs({At: slice_A, Bt: slice_B, dt: slice_d})) == 0,
        "constant-different image relation failure")
slice_groebner = sp.groebner(
    (At - slice_A, Bt - slice_B, dt - slice_d),
    u, v, At, Bt, dt,
    order="lex",
    domain=sp.QQ.frac_field(c),
)
elimination_rows = [
    row.as_expr()
    for row in slice_groebner.polys
    if not row.as_expr().has(u) and not row.as_expr().has(v)
]
require(len(elimination_rows) == 1 and
        sp.expand(4 * elimination_rows[0] - R_c) == 0,
        "constant-different contraction kernel failure")

w = dt + c
slice_cubic = u**3 - At * u - 2 * w
slice_quadratic = Bt * u**2 + 2 * At * dt * u + w * (3 * dt - c)
resultant = sp.factor(sp.resultant(slice_cubic, slice_quadratic, u))
require(sp.expand(resultant - w**2 * R_c) == 0,
        "constant-different resultant failure")
linear_subresultant = sp.subresultants(slice_cubic, slice_quadratic, u)[-2]
linear_coefficient = sp.factor(sp.Poly(linear_subresultant, u).coeff_monomial(u))
expected_linear = 4 * At**2 * dt**2 - At * Bt**2 + Bt * c**2 - 2 * Bt * c * dt - 3 * Bt * dt**2
require(sp.expand(linear_coefficient - expected_linear) == 0,
        "birational linear subresultant failure")
coefficient_domain = sp.QQ.frac_field(c)
relation_poly = sp.Poly(R_c, At, Bt, dt, domain=coefficient_domain)
exception_poly = sp.Poly((dt + c) * expected_linear, At, Bt, dt,
                         domain=coefficient_domain)
require(sp.polys.polytools.gcd(relation_poly, exception_poly).total_degree() == 0,
        "an image component is trapped in the exceptional resultant locus")

t = sp.symbols("t", nonzero=True)
singular_target = {At: -2 * t, Bt: t**2, dt: -c}
gradient = [sp.factor(sp.diff(R_c, variable).subs(singular_target)) for variable in (At, Bt, dt)]
expected_gradient = (
    4 * t**2 * (8 * c**2 + t**3),
    4 * t * (8 * c**2 + t**3),
    -8 * c * (8 * c**2 + t**3),
)
require(tuple(gradient) == expected_gradient, "singular-gradient formula failure")

first_preimage = {u: 0, v: t}
second_preimage = {u: t**2 / (2 * c), v: 0}
for label, preimage in (("first", first_preimage), ("second", second_preimage)):
    values = tuple(sp.factor(expr.subs(preimage)) for expr in (slice_A, slice_B, slice_d))
    reduced = tuple(sp.rem(sp.together(value - target).as_numer_denom()[0],
                           t**3 + 8 * c**2, t)
                    for value, target in zip(values, (-2 * t, t**2, -c)))
    require(reduced == (0, 0, 0), f"{label} singular preimage failure")
print("constant_J0_slice_image_relation_resultant_and_birational_recovery: PASS")
print("constant_J0_slice_lex_elimination_kernel=(Rc): PASS")
print("constant_J0_slice_no_extraneous_component_gcd: PASS")
print("constant_J0_slice_singular_double_preimage_t^3=-8c^2: PASS")
print("scope=fixed_S4_quotient_and_polynomial_coordinate_family_only; general_S4_JC2_DC2_open")
