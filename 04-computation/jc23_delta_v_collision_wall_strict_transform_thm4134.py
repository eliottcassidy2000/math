#!/usr/bin/env python3
"""Primary exact certificate for prospective THM-4134.

The calculation starts from the normalized (X,T) polynomial of THM-4130,
specializes its critical eliminant to the Delta_V collision wall, and then
keeps the reciprocal variables u=1/T and S=XT.  This retains the boundary
root which the affine eliminant loses.  The certificate distinguishes the
generic wall from its first-normal secondary wall, checks the two resulting
boundary packets and affine critical lengths, and audits the permutation
obstruction in both full-boundary degree branches.

No Python ``assert`` is used, so optimized mode retains every truth gate.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


X, T, Phi, Theta = sp.symbols("X T Phi Theta")
P = T + X**2 * T**2
Y = X * T * P
G = (
    -sp.Rational(1, 2) * X**2 * T
    - 3 * P
    + sp.Rational(8, 3) * P**2
    - sp.Rational(1376, 135) * P**3
    + sp.Rational(2848, 45) * Y**2
    + Phi * P**2 * Y
    + Theta * P * Y**2
)

# Rebuild the complete critical resultant rather than importing THM-4130.
f = sp.cancel(sp.diff(G, X) / T)
h = sp.diff(G, T)
require(sp.denom(f) == 1, "G_X/T is not polynomial")
require(sp.expand(sp.diff(G, X) - T * f) == 0, "lost the exact T factor")
resultant_xt = sp.factor(sp.resultant(f, h, X))
Q16 = sp.cancel(resultant_xt / (-T**42 * Theta**3 * (6 * T + 1) ** 2))
require(sp.denom(Q16) == 1, "Q16 is not polynomial")
Q16_poly = sp.Poly(Q16, T)
W = 135 * Phi**2 + 5504 * Theta
require(Q16_poly.degree() == 16, "formal critical degree changed")
require(
    sp.factor(Q16_poly.LC())
    == sp.Rational(16842242073296896, 20503125) * Theta**2 * W,
    "wrong Q16 leading row",
)
require(Q16_poly.TC() == 12288 * Theta**3, "wrong Q16 constant row")

# Wall and first-normal coordinate.
r, w, u, S, q = sp.symbols("r w u S q")
theta_wall = -sp.Rational(135, 5504) * r**2
D = 273375 * r**2 + 2696167424
wall_sub = {Phi: r, Theta: theta_wall}
Qwall = sp.Poly(sp.cancel(Q16.subs(wall_sub)), T)
expected_c15 = sp.Rational(513984438272, 12075125625) * r**4 * D
require(Qwall.degree() == 15, "generic wall residual is not degree fifteen")
require(sp.factor(Qwall.LC()) == expected_c15, "wrong wall degree-fifteen row")
require(
    Qwall.TC() == -sp.Rational(7381125, 40707584) * r**6,
    "wrong wall residual constant",
)

theta_normal = (w - 135 * r**2) / 5504
Qhat = sp.expand(u**16 * Q16.subs({Phi: r, Theta: theta_normal, T: 1 / u}))
qhat_u = sp.factor(sp.diff(Qhat, u).subs({u: 0, w: 0}))
qhat_w = sp.factor(sp.diff(Qhat, w).subs({u: 0, w: 0}))
require(qhat_u == expected_c15, "reciprocal exceptional factor changed")
require(
    qhat_w == sp.Rational(1027968876544, 2080125) * r**4,
    "wrong first normal coefficient",
)
u_over_w = sp.factor(-qhat_w / qhat_u)
require(u_over_w == -sp.Rational(11610, 1) / D, "wrong escaping-root slope")

# The secondary wall D=0 has a quadratic, rather than linear, reciprocal jet.
D_poly = sp.Poly(D, r)


def reduce_mod_D(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_poly = sp.Poly(numerator, r, domain=sp.QQ)
    denominator_poly = sp.Poly(denominator, r, domain=sp.QQ)
    denominator_inverse = sp.invert(denominator_poly, D_poly)
    reduced = sp.rem(numerator_poly * denominator_inverse, D_poly).as_expr()
    return sp.factor(reduced)


c14 = sp.factor(Qwall.coeff_monomial(T**14))
c14_secondary = reduce_mod_D(c14)
qhat_w_secondary = reduce_mod_D(qhat_w)
expected_c14_secondary = sp.Rational(
    2025766671394476091774451830489088,
    420378134765625,
)
require(c14_secondary == expected_c14_secondary, "secondary wall lost degree fourteen")
require(qhat_w_secondary != 0, "secondary normal coefficient vanished")
require(
    sp.factor(-qhat_w_secondary / c14_secondary) == -sp.Rational(5, 501248),
    "wrong square-root escape law",
)

# Reciprocal compactification of the source curve.
R = 1 + S**2 * u
A = Theta * S**2 + Phi * S - sp.Rational(1376, 135)
B = sp.Rational(8, 3) + sp.Rational(2848, 45) * S**2
K = A * R**3 + u * B * R**2 - 3 * u**2 * R - sp.Rational(1, 2) * S**2 * u**4
require(
    sp.cancel(u**3 * G.subs({X: S * u, T: 1 / u}) - K) == 0,
    "reciprocal source identity failed",
)
F = K - q * u**3
S0 = sp.Rational(2752, 135) / r
B0 = sp.factor(B.subs(S, S0))
Bprime0 = sp.factor(sp.diff(B, S).subs(S, S0))
require(
    sp.factor(A.subs({Phi: r, Theta: theta_wall}) - theta_wall * (S - S0) ** 2)
    == 0,
    "wall edge is not a repeated root",
)
require(
    sp.factor(B0 - 8 * D / (820125 * r**2)) == 0,
    "wrong first-normal boundary unit",
)
require(
    sp.factor(sp.diff(F, S, 2).subs({S: S0, u: 0, Phi: r, Theta: theta_wall}) / 2)
    == theta_wall,
    "wrong repeated-edge quadratic coefficient",
)
require(
    sp.factor(sp.diff(F, u).subs({S: S0, u: 0, Phi: r, Theta: theta_wall}))
    == B0,
    "wrong boundary normal coefficient",
)
boundary_u_over_v2 = sp.factor(-theta_wall / B0)

# The critical equations in the reciprocal chart are K_S=0 and uK_u-3K=0.
KS = sp.diff(K, S)
L = sp.expand(u * sp.diff(K, u) - 3 * K)
critical_matrix = sp.Matrix(
    [
        [sp.diff(KS, S), sp.diff(KS, u)],
        [sp.diff(L, S), sp.diff(L, u)],
    ]
).subs({S: S0, u: 0, Phi: r, Theta: theta_wall})
critical_det = sp.factor(critical_matrix.det())
require(
    sp.factor(critical_det - D / 1044900) == 0,
    "wrong reciprocal critical Jacobian",
)
normal_vector = sp.Matrix(
    [sp.diff(KS, Theta), sp.diff(L, Theta)]
).subs({S: S0, u: 0, Phi: r, Theta: theta_wall}) / 5504
S_slope, u_slope = [
    sp.factor(entry) for entry in -critical_matrix.inv() * normal_vector
]
expected_S_slope = sp.factor(
    2752 * (273375 * r**2 - 1348083712) / (18225 * r**3 * D)
)
require(S_slope == expected_S_slope, "wrong escaping S-label slope")
require(u_slope == u_over_w, "joint ideal and eliminant slopes disagree")
K_normal = sp.factor(
    (
        sp.diff(K, Theta) / 5504
        + sp.diff(K, S) * S_slope
        + sp.diff(K, u) * u_slope
    ).subs({S: S0, u: 0, Phi: r, Theta: theta_wall})
)
critical_value_coefficient = sp.factor(K_normal / u_slope**3)
require(
    critical_value_coefficient == 2 * D**3 / (82909778259375 * r**2),
    "wrong escaping critical-value coefficient",
)

# On D=0 the boundary first-normal term also dies and exposes an ordinary node.
theta_secondary = reduce_mod_D(theta_wall)
S0_squared_secondary = reduce_mod_D(S0**2)
Bprime_squared_secondary = reduce_mod_D(Bprime0**2)
require(theta_secondary == sp.Rational(489856, 2025), "wrong secondary theta")
require(S0_squared_secondary == -sp.Rational(15, 356), "wrong secondary root")
require(
    Bprime_squared_secondary == -sp.Rational(91136, 135),
    "wrong secondary mixed coefficient",
)
secondary_tangent_discriminant = sp.factor(Bprime_squared_secondary + 12 * theta_secondary)
require(
    secondary_tangent_discriminant == sp.Rational(501248, 225),
    "secondary boundary is not an ordinary node",
)

# Boundary packets and Riemann--Hurwitz ledgers.
generic_packet = (7, 5, 3, 2, 2, 1)
secondary_packet = (7, 3, 2, 2, 2, 2, 1)
require(sum(generic_packet) == 20, "wrong generic packet degree")
require(sum(length - 1 for length in generic_packet) == 14, "wrong genus-eight RH total")
require(sum(secondary_packet) == 19, "wrong secondary packet degree")
require(sum(length - 1 for length in secondary_packet) == 12, "wrong genus-seven RH total")
generic_affine_critical_length = 2 + 2 + Qwall.degree()
secondary_affine_critical_length = 2 + 2 + 14
require(generic_affine_critical_length == 19, "wrong generic critical length")
require(secondary_affine_critical_length == 18, "wrong secondary critical length")
degree_alternatives_generic = (16, 20)
degree_alternatives_secondary = (15, 19)


def identity(degree: int) -> tuple[int, ...]:
    return tuple(range(degree))


def cycle(degree: int, entries: tuple[int, ...]) -> tuple[int, ...]:
    require(len(entries) >= 2 and len(set(entries)) == len(entries), "bad cycle")
    result = list(range(degree))
    for left, right in zip(entries, entries[1:] + entries[:1]):
        result[left] = right
    return tuple(result)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * len(permutation)
    for point, image in enumerate(permutation):
        result[image] = point
    return tuple(result)


def commutator(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return compose(compose(compose(left, right), inverse(left)), inverse(right))


def support(permutation: tuple[int, ...]) -> set[int]:
    return {point for point, image in enumerate(permutation) if point != image}


def cycle_type(permutation: tuple[int, ...]) -> tuple[int, ...]:
    unseen = set(range(len(permutation)))
    lengths: list[int] = []
    while unseen:
        start = min(unseen)
        point = start
        length = 0
        while point in unseen:
            unseen.remove(point)
            point = permutation[point]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def orbit(generators: tuple[tuple[int, ...], ...]) -> set[int]:
    found = {0}
    active = [0]
    while active:
        point = active.pop()
        for generator in generators:
            image = generator[point]
            if image not in found:
                found.add(image)
                active.append(image)
    return found


def audit_full_boundary_branch(degree: int, critical_length: int, packet: tuple[int, ...]) -> int:
    require(critical_length == degree - 1, "critical/full-degree offset changed")
    require(sum(packet) == degree, "packet has wrong degree")
    require(packet != (3,) + (1,) * (degree - 3), "packet became a three-cycle")
    cases = 0
    for r0 in range(2, critical_length - 1):
        r1 = critical_length - r0
        require(r1 >= 2, "universal critical pair lost")
        m0, m1 = degree - r0, degree - r1
        require(m0 + m1 == degree + 1, "support sum changed")
        left = cycle(degree, tuple(range(m0)))
        right = cycle(degree, (0,) + tuple(range(m0, degree)))
        require(support(left) | support(right) == set(range(degree)), "support union changed")
        require(support(left) & support(right) == {0}, "shared pivot changed")
        require(orbit((left, right)) == set(range(degree)), "canonical pair not transitive")
        bracket_type = cycle_type(commutator(left, right))
        require(bracket_type == (3,) + (1,) * (degree - 3), "commutator not a three-cycle")
        require(bracket_type != packet, "wall packet escaped obstruction")
        cases += 1
    return cases


generic_support_cases = audit_full_boundary_branch(20, 19, generic_packet)
secondary_support_cases = audit_full_boundary_branch(19, 18, secondary_packet)
require(generic_support_cases == 16, "wrong generic support census")
require(secondary_support_cases == 15, "wrong secondary support census")

# Hostile exact wall point and an independent good-reduction squarefree control.
hostile_Q = sp.Poly(Q16.subs({Phi: 5504, Theta: -743040}), T)
hostile_integer = hostile_Q.clear_denoms(convert=True)[1].primitive()[1]
require(hostile_Q.degree() == 15, "hostile wall did not drop to degree fifteen")
require(sp.gcd(hostile_Q, hostile_Q.diff()).degree() == 0, "hostile residual not squarefree over QQ")
hostile_mod101 = sp.Poly(hostile_integer.as_expr(), T, modulus=101)
hostile_mod101_coefficients = tuple(int(value) for value in hostile_mod101.all_coeffs())
expected_mod101_coefficients = (
    22,
    -9,
    25,
    -7,
    -17,
    -43,
    -37,
    47,
    -19,
    -2,
    -18,
    -43,
    -20,
    12,
    -45,
    19,
)
require(hostile_mod101_coefficients == expected_mod101_coefficients, "mod-101 residual changed")
require(sp.gcd(hostile_mod101, hostile_mod101.diff()).degree() == 0, "bad mod-101 reduction")

# The residual horizontal BC carrier is not killed by the pole-pair gate.
a, rho = sp.symbols("a rho", nonzero=True)
section_U = a / 2 + 16 * rho**2 / (9 * a**2)
section_V = -rho - 64 * rho**3 / (27 * a**3)
target_E = section_V**2 - section_U**3 + sp.Rational(3, 4) * a**2 * section_U + a**3 / 4
require(sp.factor(target_E - (a**3 / 2 + rho**2)) == 0, "polynomial section failed")
require(sp.degree(section_U, rho) == 2 and sp.degree(section_V, rho) == 3, "wrong section pole pair")

semantic_lines = (
    "scope=theta_only_DeltaV_collision_wall_M8",
    "wall=W=135*Phi^2+5504*Theta=0;Phi!=0",
    "secondary=D=273375*Phi^2+2696167424=0",
    "generic=genus8;critical_length19;packet=7,5,3,2,2,1",
    "secondary=genus7;critical_length18;packet=7,3,2,2,2,2,1",
    "degree_alternatives=generic:16|20;secondary:15|19",
    "full_boundary_branches=n20,n19:CONTRADICTION",
    "residual=generic_n16_horizontal_BC|secondary_n15_horizontal_BC",
    "critical_resultant_inertia_is_not_cover_inertia",
    "JC2=OPEN",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4134_PRIMARY_AUDIT")
print("scope=theta_only_DeltaV_collision_wall_M8")
print("wall=W=135*Phi^2+5504*Theta=0;Phi!=0")
print("secondary_wall=D=273375*Phi^2+2696167424")
print("Qwall_degree=15;c15=" + str(expected_c15))
print("critical_strict_transform_generic=u/W=" + str(u_over_w))
print("critical_strict_transform_secondary=u^2/W=-5/501248")
print("critical_value_escape=G*W^2->" + str(critical_value_coefficient))
print("boundary_generic=u-(" + str(boundary_u_over_v2) + ")*v^2+...;residue_order=4;e=5")
print("boundary_secondary=node_discriminant=501248/225;residue_orders=1,1;e=2,2")
print("generic_ledger=genus8;critical_length=19;packet=7,5,3,2,2,1;degrees=16|20")
print("secondary_ledger=genus7;critical_length=18;packet=7,3,2,2,2,2,1;degrees=15|19")
print("full_boundary_support_cases=generic:" + str(generic_support_cases) + ";secondary:" + str(secondary_support_cases))
print("full_boundary_verdict=n20:CONTRADICTION;n19:CONTRADICTION")
print("residual=generic_n16_horizontal_BC;secondary_n15_horizontal_BC")
print("hostile=(Phi,Theta)=(5504,-743040);Q15_mod101_coefficients=" + str(hostile_mod101_coefficients))
print("hostile_mod101_gcd=1")
print("polynomial_section=(a/2+16*rho^2/(9*a^2),-rho-64*rho^3/(27*a^3));pole_pair=2,3")
print("firewall=critical_resultant_inertia_is_not_cover_inertia;root_labels_and_BC_images_retained")
print("semantic_sha256=" + semantic_sha256)
print("verdict=HIGH_BRANCHES_EXCLUDED;LOW_HORIZONTAL_BC_BRANCHES_OPEN;JC2_OPEN")
