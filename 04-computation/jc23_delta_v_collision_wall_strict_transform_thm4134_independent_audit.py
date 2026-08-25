#!/usr/bin/env python3
"""Independent exact audit for prospective THM-4134.

This script imports no primary implementation.  It works in the rational
coordinates s=XT, p=T+s^2, eliminates s from a different critical pair, and
uses z=1/p for the boundary normal calculation.  A dictionary permutation
implementation independently checks the two full-boundary contradictions.

No Python ``assert`` is used, so optimized mode retains every truth gate.
"""

from __future__ import annotations

import hashlib

import sympy as sy


def check(condition: bool, label: str) -> None:
    if not condition:
        raise ArithmeticError(label)


s, p, Phi, Theta = sy.symbols("s p Phi Theta")
t = p - s**2
kappa = sy.Rational(2848, 45)
H = (
    -3 * p
    + sy.Rational(8, 3) * p**2
    - sy.Rational(1376, 135) * p**3
    + kappa * s**2 * p**2
    + Phi * s * p**3
    + Theta * s**2 * p**3
)
Gsp = -sy.Rational(1, 2) * s**2 / t + H

# Critical equations on t*p != 0, derived directly in (s,p).
Acrit = sy.expand(
    -s + t**2 * p * (2 * kappa * s + Phi * p + 2 * Theta * s * p)
)
Bcrit = sy.expand(
    -6
    + sy.Rational(32, 3) * p
    - sy.Rational(2752, 45) * p**2
    + 6 * kappa * s**2 * p
    + 7 * Phi * s * p**2
    + 8 * Theta * s**2 * p**2
)
Gs = sy.factor(sy.diff(Gsp, s))
Gp = sy.factor(sy.diff(Gsp, p))
check(sy.cancel(Gs - p * Acrit / t**2) == 0, "wrong s-critical equation")
check(
    sy.cancel(2 * t**2 * Gp + s * Acrit - t**2 * Bcrit) == 0,
    "wrong p-critical reduction",
)
resultant_sp = sy.factor(sy.resultant(Acrit, Bcrit, s))
R16 = sy.cancel(resultant_sp / (-sy.Rational(2, 373669453125) * p**2))
check(sy.denom(R16) == 1, "R16 is not polynomial")
R16_poly = sy.Poly(R16, p)
W = 135 * Phi**2 + 5504 * Theta
check(R16_poly.degree() == 16, "formal independent degree changed")
check(
    sy.factor(R16_poly.LC()) == 34012224000000 * Theta**5 * W,
    "wrong independent leading row",
)
check(R16_poly.TC() == 23277095392665600000, "wrong independent constant row")
check(Acrit.subs(p, 0) == -s and Bcrit.subs(p, 0) == -6, "p^2 artifact became finite")

r, w, z, q = sy.symbols("r w z q")
theta_wall = -sy.Rational(135, 5504) * r**2
D = 273375 * r**2 + 2696167424
Rwall = sy.Poly(sy.cancel(R16.subs({Phi: r, Theta: theta_wall})), p)
expected_c15 = -sy.Rational(2758100942197265625, 106054873287491584) * r**10 * D
check(Rwall.degree() == 15, "independent wall residual is not degree fifteen")
check(sy.factor(Rwall.LC()) == expected_c15, "wrong independent wall leading row")
check(Rwall.TC() == 23277095392665600000, "wrong independent wall constant")

theta_normal = (w - 135 * r**2) / 5504
Rhat = sy.expand(z**16 * R16.subs({Phi: r, Theta: theta_normal, p: 1 / z}))
rhat_z = sy.factor(sy.diff(Rhat, z).subs({z: 0, w: 0}))
rhat_w = sy.factor(sy.diff(Rhat, w).subs({z: 0, w: 0}))
check(rhat_z == expected_c15, "independent exceptional factor changed")
check(
    rhat_w == -sy.Rational(372343627196630859375, 1233196201017344) * r**10,
    "wrong independent normal coefficient",
)
z_over_w = sy.factor(-rhat_w / rhat_z)
check(z_over_w == -sy.Rational(11610, 1) / D, "independent escape slope changed")

D_poly = sy.Poly(D, r, domain=sy.QQ)


def reduce_mod_D(expression: sy.Expr) -> sy.Expr:
    numerator, denominator = sy.cancel(expression).as_numer_denom()
    numerator_poly = sy.Poly(numerator, r, domain=sy.QQ)
    denominator_poly = sy.Poly(denominator, r, domain=sy.QQ)
    denominator_inverse = sy.invert(denominator_poly, D_poly)
    return sy.factor(sy.rem(numerator_poly * denominator_inverse, D_poly).as_expr())


c14_secondary = reduce_mod_D(Rwall.coeff_monomial(p**14))
rhat_w_secondary = reduce_mod_D(rhat_w)
check(
    c14_secondary
    == sy.Rational(57910159368111680226053516112070443008, 20503125),
    "independent secondary degree-fourteen row changed",
)
check(rhat_w_secondary != 0, "independent secondary normal coefficient vanished")
check(
    sy.factor(-rhat_w_secondary / c14_secondary) == -sy.Rational(5, 501248),
    "independent square-root escape law changed",
)

# Restore the two strata omitted by the rational critical equations.
Gp_p0 = sy.factor((2 * Gp).subs(p, 0))
check(
    sy.cancel(Gp_p0 - (1 - 6 * s**2) / s**2) == 0,
    "wrong p=0 critical pair",
)
check(sy.cancel((Gsp - sy.Rational(1, 2)).subs(p, 0)) == 0, "wrong p=0 value")
X, T = sy.symbols("X T")
Pxt = T + X**2 * T**2
Yxt = X * T * Pxt
Gxt = (
    -sy.Rational(1, 2) * X**2 * T
    - 3 * Pxt
    + sy.Rational(8, 3) * Pxt**2
    - sy.Rational(1376, 135) * Pxt**3
    + kappa * Yxt**2
    + Phi * Pxt**2 * Yxt
    + Theta * Pxt * Yxt**2
)
check(
    sy.expand(sy.diff(Gxt, T).subs(T, 0) + (X**2 + 6) / 2) == 0,
    "wrong t=0 critical pair",
)
check(Gxt.subs(T, 0) == 0, "wrong t=0 critical value")
generic_affine_critical_length = 2 + 2 + Rwall.degree()
secondary_affine_critical_length = 2 + 2 + 14
check(generic_affine_critical_length == 19, "wrong independent generic length")
check(secondary_affine_critical_length == 18, "wrong independent secondary length")

# Independent boundary chart z=1/p.  Multiplication by the displayed unit
# clears the rational t denominator without changing the local curve.
J = sy.cancel((1 - s**2 * z) * z**3 * (Gsp.subs(p, 1 / z) - q))
check(sy.denom(J) == 1, "boundary normal equation is not polynomial")
s0 = sy.Rational(2752, 135) / r
wall_point = {s: s0, z: 0, Phi: r, Theta: theta_wall}
check(J.subs(wall_point) == 0, "wall point missed the boundary curve")
check(sy.factor(sy.diff(J, s).subs(wall_point)) == 0, "edge root is not repeated")
boundary_z = sy.factor(sy.diff(J, z).subs(wall_point))
boundary_ss = sy.factor(sy.diff(J, s, 2).subs(wall_point) / 2)
boundary_sz = sy.factor(sy.diff(J, s, z).subs(wall_point))
boundary_zz = sy.factor(sy.diff(J, z, 2).subs(wall_point) / 2)
check(
    sy.factor(boundary_z - 8 * D / (820125 * r**2)) == 0,
    "wrong independent first-normal boundary unit",
)
check(boundary_ss == theta_wall, "wrong independent repeated-edge coefficient")
check(
    boundary_sz == sy.Rational(15675392, 6075) / r,
    "wrong independent mixed boundary coefficient",
)
check(reduce_mod_D(boundary_zz) == -3, "wrong secondary z^2 coefficient")
theta_secondary = reduce_mod_D(theta_wall)
boundary_sz_squared_secondary = reduce_mod_D(boundary_sz**2)
secondary_tangent_discriminant = sy.factor(
    boundary_sz_squared_secondary - 4 * theta_secondary * reduce_mod_D(boundary_zz)
)
check(
    secondary_tangent_discriminant == sy.Rational(501248, 225),
    "independent secondary boundary is not nodal",
)

# ds wedge dp/[t(G-q)] becomes -z^2 ds wedge dz/J.  These powers give
# residue orders four on the smooth tangency and one on each node branch.
residue_factor = sy.cancel(
    (-z**-2) * (z / (1 - s**2 * z)) * ((1 - s**2 * z) * z**3)
)
check(residue_factor == -z**2, "wrong independent residue numerator")

generic_packet = (7, 5, 3, 2, 2, 1)
secondary_packet = (7, 3, 2, 2, 2, 2, 1)
check(sum(generic_packet) == 20, "wrong independent generic packet")
check(sum(length - 1 for length in generic_packet) == 14, "wrong independent genus-eight RH")
check(sum(secondary_packet) == 19, "wrong independent secondary packet")
check(sum(length - 1 for length in secondary_packet) == 12, "wrong independent genus-seven RH")


# Dictionary permutations share no implementation with the primary audit.
def make_cycle(degree: int, word: tuple[int, ...]) -> dict[int, int]:
    check(len(word) >= 2 and len(set(word)) == len(word), "bad dictionary cycle")
    result = {point: point for point in range(degree)}
    for index, point in enumerate(word):
        result[point] = word[(index + 1) % len(word)]
    return result


def invert(mapping: dict[int, int]) -> dict[int, int]:
    return {image: point for point, image in mapping.items()}


def multiply(left: dict[int, int], right: dict[int, int]) -> dict[int, int]:
    return {point: left[right[point]] for point in right}


def bracket(left: dict[int, int], right: dict[int, int]) -> dict[int, int]:
    return multiply(multiply(multiply(left, right), invert(left)), invert(right))


def moved(mapping: dict[int, int]) -> frozenset[int]:
    return frozenset(point for point, image in mapping.items() if point != image)


def generated_orbit(generators: tuple[dict[int, int], ...]) -> frozenset[int]:
    found = {0}
    active = [0]
    while active:
        point = active.pop()
        for generator in generators:
            image = generator[point]
            if image not in found:
                found.add(image)
                active.append(image)
    return frozenset(found)


def cycle_lengths(mapping: dict[int, int]) -> tuple[int, ...]:
    unseen = set(mapping)
    lengths: list[int] = []
    while unseen:
        start = min(unseen)
        point = start
        length = 0
        while point in unseen:
            unseen.remove(point)
            point = mapping[point]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def audit_full_boundary_branch(degree: int, critical_length: int, packet: tuple[int, ...]) -> int:
    check(critical_length == degree - 1, "independent critical/full-degree offset changed")
    cases = 0
    for r0 in range(2, critical_length - 1):
        r1 = critical_length - r0
        check(r1 >= 2, "independent universal pair lost")
        m0, m1 = degree - r0, degree - r1
        check(m0 + m1 == degree + 1, "independent support sum changed")
        left = make_cycle(degree, tuple(range(m0)))
        right = make_cycle(degree, (0,) + tuple(range(m0, degree)))
        check(moved(left) & moved(right) == frozenset({0}), "independent pivot changed")
        check(moved(left) | moved(right) == frozenset(range(degree)), "independent union changed")
        check(generated_orbit((left, right)) == frozenset(range(degree)), "independent nontransitive pair")
        bracket_type = cycle_lengths(bracket(left, right))
        check(bracket_type == (3,) + (1,) * (degree - 3), "independent bracket changed")
        check(bracket_type != packet, "independent wall packet escaped")
        cases += 1
    return cases


generic_support_cases = audit_full_boundary_branch(20, 19, generic_packet)
secondary_support_cases = audit_full_boundary_branch(19, 18, secondary_packet)
check(generic_support_cases == 16, "wrong independent generic support census")
check(secondary_support_cases == 15, "wrong independent secondary support census")

# Distinct hostile residual in the p-projection.
hostile_R = sy.Poly(R16.subs({Phi: 5504, Theta: -743040}), p)
hostile_integer = hostile_R.clear_denoms(convert=True)[1].primitive()[1]
check(hostile_R.degree() == 15, "independent hostile degree changed")
check(sy.gcd(hostile_R, hostile_R.diff()).degree() == 0, "independent hostile not squarefree")
hostile_mod101 = sy.Poly(hostile_integer.as_expr(), p, modulus=101)
hostile_mod101_coefficients = tuple(int(value) for value in hostile_mod101.all_coeffs())
expected_mod101_coefficients = (
    27,
    13,
    19,
    -32,
    -45,
    -48,
    46,
    8,
    -5,
    11,
    -27,
    43,
    46,
    -17,
    -40,
    38,
)
check(hostile_mod101_coefficients == expected_mod101_coefficients, "independent mod-101 residual changed")
check(sy.gcd(hostile_mod101, hostile_mod101.diff()).degree() == 0, "independent bad reduction")

# Reconstruct the degree-(2,3) horizontal positive control by addition on the
# quadratic target base change, rather than importing the primary formula.
a, rho = sy.symbols("a rho", nonzero=True)
double_U, double_V = -a, -rho
addition_slope = sy.factor((rho - double_V) / (a / 2 - double_U))
section_U = sy.factor(addition_slope**2 - a / 2 - double_U)
section_V = sy.factor(addition_slope * (a / 2 - section_U) - rho)
check(addition_slope == 4 * rho / (3 * a), "wrong section addition slope")
target_E = section_V**2 - section_U**3 + sy.Rational(3, 4) * a**2 * section_U + a**3 / 4
check(sy.factor(target_E - (a**3 / 2 + rho**2)) == 0, "independent polynomial section failed")
check(sy.degree(section_U, rho) == 2 and sy.degree(section_V, rho) == 3, "wrong independent pole pair")

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

print("THM4134_INDEPENDENT_AUDIT")
print("scope=theta_only_DeltaV_collision_wall_M8")
print("coordinate_route=(s,p);critical_projection=p;boundary_coordinate=z=1/p")
print("wall=W=135*Phi^2+5504*Theta=0;Phi!=0")
print("secondary_wall=D=273375*Phi^2+2696167424")
print("Rwall_degree=15;c15=" + str(expected_c15))
print("critical_strict_transform_generic=z/W=" + str(z_over_w))
print("critical_strict_transform_secondary=z^2/W=-5/501248")
print("boundary_generic=smooth_tangency;residue_numerator=-z^2;residue_order=4;e=5")
print("boundary_secondary=node_discriminant=501248/225;residue_orders=1,1;e=2,2")
print("generic_ledger=genus8;critical_length=19;packet=7,5,3,2,2,1;degrees=16|20")
print("secondary_ledger=genus7;critical_length=18;packet=7,3,2,2,2,2,1;degrees=15|19")
print("full_boundary_support_cases=generic:" + str(generic_support_cases) + ";secondary:" + str(secondary_support_cases))
print("full_boundary_verdict=n20:CONTRADICTION;n19:CONTRADICTION")
print("residual=generic_n16_horizontal_BC;secondary_n15_horizontal_BC")
print("hostile=(Phi,Theta)=(5504,-743040);R15_mod101_coefficients=" + str(hostile_mod101_coefficients))
print("hostile_mod101_gcd=1")
print("polynomial_section=" + str((section_U, section_V)) + ";pole_pair=2,3")
print("firewall=critical_resultant_inertia_is_not_cover_inertia;finite_BC_images_unresolved")
print("semantic_sha256=" + semantic_sha256)
print("verdict=ACCEPT;HIGH_BRANCHES_EXCLUDED;LOW_HORIZONTAL_BC_BRANCHES_OPEN;JC2_OPEN")
