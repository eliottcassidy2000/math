#!/usr/bin/env python3
"""Exact primary audit for the proposed THM-4138 Delta_V wall closure.

This script checks the algebraic/Mordell--Weil carrier, the forced nodal
horizontal image, the marked vanishing-cycle geometry, and the sharp
orbit-merger arithmetic.  The Shioda--Tate and punctured-surface arguments
are printed as typed proof ledgers; they are not replaced by computation.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(poly: sp.Expr, variable: sp.Symbol) -> int:
    expanded = sp.Poly(sp.expand(poly), variable)
    return min(monomial[0] for monomial, coefficient in expanded.terms() if coefficient)


# ---------------------------------------------------------------------------
# 1. The BC divisor and its forced target section.
# ---------------------------------------------------------------------------

a, r, z, U, V = sp.symbols("a r z U V", nonzero=True)
target_E = V**2 - U**3 + sp.Rational(3, 4) * a**2 * U + a**3 / 4

section_U = a / 2 + 16 * r**2 / (9 * a**2)
section_V = -r - 64 * r**3 / (27 * a**3)
section_q = a**3 / 2 + r**2
require(
    sp.factor(target_E.subs({U: section_U, V: section_V}) - section_q) == 0,
    "the displayed 3P section left the target pencil",
)
require(sp.degree(section_U, r) == 2, "wrong U pole order")
require(sp.degree(section_V, r) == 3, "wrong V pole order")

curve_equation = V**2 - (U - a / 2) * (U + a / 4) ** 2
require(
    sp.factor(curve_equation.subs({U: section_U, V: section_V})) == 0,
    "the 3P image left the forced nodal cubic",
)

param_U = a / 2 + z**2
param_V = z * (z**2 + 3 * a / 4)
q_on_curve = sp.factor(target_E.subs({U: param_U, V: param_V}))
require(q_on_curve == a**2 * (8 * a + 9 * z**2) / 16, "wrong q-map on S")
require(sp.degree(q_on_curve, z) == 2, "horizontal component lost quadratic q-degree")

node_U = -a / 4
node_q = sp.factor(q_on_curve.subs(z**2, -3 * a / 4))
require(node_q == 5 * a**3 / 64, "wrong horizontal-node pencil value")
require(
    sp.factor(curve_equation.subs({U: node_U, V: 0})) == 0,
    "horizontal node left S",
)
require(
    sp.factor(sp.diff(curve_equation, U).subs({U: node_U, V: 0})) == 0
    and sp.factor(sp.diff(curve_equation, V).subs({U: node_U, V: 0})) == 0,
    "horizontal singular point is not singular",
)
node_hessian_det = sp.factor(
    sp.det(
        sp.Matrix(
            [
                [sp.diff(curve_equation, U, U), sp.diff(curve_equation, U, V)],
                [sp.diff(curve_equation, V, U), sp.diff(curve_equation, V, V)],
            ]
        ).subs({U: node_U, V: 0})
    )
)
require(node_hessian_det == 3 * a, "horizontal singularity is not an ordinary node")

o1_U = a / 2
o1_gradient = sp.factor(sp.diff(curve_equation, U).subs({U: o1_U, V: 0}))
require(o1_gradient == -9 * a**2 / 16, "S is not smooth at the second target node")
require(
    sp.factor(q_on_curve.subs(z, 0) - a**3 / 2) == 0,
    "S does not pass through o1 at q=a^3/2",
)
o0_intersection_U = sp.factor(param_U.subs(z**2, -8 * a / 9))
require(o0_intersection_U == -7 * a / 18, "wrong S intersection with q=0")
require(o0_intersection_U != -a / 2, "S unexpectedly passes through the o0 node")

# The degree tower B -> normalization(S) -> q has product two.  The e=1
# alternative would give an affine k(q)-point on E_q, contrary to the exact
# target Mordell--Weil computation in THM-4120.
degree_tower_options = ((1, 2), (2, 1))
require(all(d * e == 2 for d, e in degree_tower_options), "bad degree tower")
surviving_degree_tower = (1, 2)


# ---------------------------------------------------------------------------
# 2. Exact Mordell--Weil and polynomial-section audit after normalization.
# ---------------------------------------------------------------------------

t, u, x, y = sp.symbols("t u x y")
weierstrass_A = sp.Integer(-3)
weierstrass_B = t**2 + 2
weierstrass_delta = sp.factor(-16 * (4 * weierstrass_A**3 + 27 * weierstrass_B**2))
require(weierstrass_delta == -432 * t**2 * (t**2 + 4), "wrong normalized discriminant")

c4 = -48 * weierstrass_A
c6 = -864 * weierstrass_B
require(c4 == 144, "wrong finite c4")
require(c6 == -864 * (t**2 + 2), "wrong finite c6")

A_infinity = -3 * u**4
B_infinity = u**4 + 2 * u**6
c4_infinity = -48 * A_infinity
c6_infinity = -864 * B_infinity
delta_infinity = sp.factor(-16 * (4 * A_infinity**3 + 27 * B_infinity**2))
infinity_valuations = (
    valuation(c4_infinity, u),
    valuation(c6_infinity, u),
    valuation(delta_infinity, u),
)
require(infinity_valuations == (4, 4, 8), "infinity is not the IV* row")

fiber_root_ranks = {"I2": 1, "IV*": 6}
mw_rank = 10 - 2 - sum(fiber_root_ranks.values())
require(mw_rank == 1, "wrong Shioda--Tate rank")

height_P = sp.factor(2 - sp.Rational(1, 2) - sp.Rational(4, 3))
require(height_P == sp.Rational(1, 6), "wrong generator height")
trivial_lattice_discriminant = 2 * 3
primitive_pairs = [
    (index, torsion_order)
    for index in range(1, 9)
    for torsion_order in range(1, 9)
    if sp.Rational(index**2 * torsion_order**2, trivial_lattice_discriminant) == height_P
]
require(primitive_pairs == [(1, 1)], "height/determinant did not force primitive torsion-free MW")


Point = tuple[sp.Expr, sp.Expr] | None


def elliptic_add(left: Point, right: Point) -> Point:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if sp.cancel(x1 - x2) == 0:
        if sp.cancel(y1 + y2) == 0:
            return None
        slope = sp.cancel((3 * x1**2 + weierstrass_A) / (2 * y1))
    else:
        slope = sp.cancel((y2 - y1) / (x2 - x1))
    x3 = sp.cancel(slope**2 - x1 - x2)
    y3 = sp.cancel(slope * (x1 - x3) - y1)
    return x3, y3


def on_surface(point: Point) -> bool:
    if point is None:
        return True
    x0, y0 = point
    return sp.factor(y0**2 - x0**3 + 3 * x0 - 2 - t**2) == 0


P: Point = (sp.Integer(1), t)
multiples: dict[int, Point] = {}
current: Point = None
denominator_degree_rows: list[tuple[int, int, int]] = []
for n in range(1, 13):
    current = elliptic_add(current, P)
    require(current is not None and on_surface(current), f"{n}P left the elliptic surface")
    multiples[n] = current
    x_n = sp.cancel(current[0])
    numerator, denominator = sp.fraction(x_n)
    require(
        sp.gcd(sp.Poly(numerator, t), sp.Poly(denominator, t)).degree() == 0,
        f"{n}P x-coordinate was not reduced",
    )
    denominator_degree = int(sp.degree(denominator, t))
    epsilon_2 = int(n % 2 != 0)
    epsilon_3 = int(n % 3 != 0)
    intersection_numerator = n * n - 12 + 3 * epsilon_2 + 8 * epsilon_3
    require(intersection_numerator % 12 == 0, f"nonintegral height intersection for {n}P")
    finite_O_intersection = intersection_numerator // 12
    require(finite_O_intersection >= 0, f"negative height intersection for {n}P")
    require(
        denominator_degree == 2 * finite_O_intersection,
        f"denominator/height mismatch for {n}P",
    )
    denominator_degree_rows.append((n, denominator_degree, finite_O_intersection))

expected_P2 = (-2, -t)
expected_P3 = (
    (4 * t**2 + 9) / 9,
    -t * (8 * t**2 + 27) / 27,
)
expected_P4 = (
    (16 * t**2 + 81) / (4 * t**2),
    (8 * t**4 + 216 * t**2 + 729) / (8 * t**3),
)
for actual, expected, label in (
    (multiples[2], expected_P2, "2P"),
    (multiples[3], expected_P3, "3P"),
    (multiples[4], expected_P4, "4P"),
):
    require(
        actual is not None
        and all(sp.factor(sp.cancel(a0 - e0)) == 0 for a0, e0 in zip(actual, expected)),
        f"wrong coordinates for {label}",
    )

# At IV*, 3P reduces to a nonzero point of the additive identity component.
# Thus no nonzero multiple meets O at infinity in characteristic zero; the
# height intersection above is entirely at finite t.
x3, y3 = multiples[3]  # type: ignore[misc]
X3_infinity = sp.limit(u**2 * x3.subs(t, 1 / u), u, 0)
Y3_infinity = sp.limit(u**3 * y3.subs(t, 1 / u), u, 0)
identity_component_parameter = sp.factor(X3_infinity / Y3_infinity)
require((X3_infinity, Y3_infinity) == (sp.Rational(4, 9), -sp.Rational(8, 27)), "wrong 3P reduction")
require(identity_component_parameter == -sp.Rational(3, 2), "3P reduced to the identity")

polynomial_multiples = tuple(n for n, _, intersection in denominator_degree_rows if intersection == 0)
require(polynomial_multiples == (1, 2, 3), "unexpected polynomial positive multiple")
require(denominator_degree_rows[3] == (4, 2, 1), "4P is not the first denominator hostile")

# P and 2P have pole pair (0,1); only +/-3P matches THM-4122's positive
# reduced (2 rho,3 rho) pole pair, and it has rho=1.
pole_pairs = {
    1: (sp.degree(multiples[1][0], t), sp.degree(multiples[1][1], t)),  # type: ignore[index]
    2: (sp.degree(multiples[2][0], t), sp.degree(multiples[2][1], t)),  # type: ignore[index]
    3: (sp.degree(multiples[3][0], t), sp.degree(multiples[3][1], t)),  # type: ignore[index]
}
require(pole_pairs == {1: (0, 1), 2: (0, 1), 3: (2, 3)}, "wrong integral-section pole pairs")


# ---------------------------------------------------------------------------
# 3. Marked vanishing cycles and the complete meridian set.
# ---------------------------------------------------------------------------

# Scale a=1 and use the same q*=1/4 reference as THM-4130.  The forced BC
# points lie on the unpushed second real vanishing cycle, so a parallel
# push-off is load-bearing rather than cosmetic.
reference_q = sp.Rational(1, 4)
reference_r_squared = reference_q - sp.Rational(1, 2)
reference_U = sp.factor(section_U.subs({a: 1, r**2: reference_r_squared}))
reference_V_squared = sp.factor(section_V.subs(a, 1) ** 2).subs(r**2, reference_r_squared)
reference_V_squared = sp.factor(reference_V_squared)
require(reference_r_squared == -sp.Rational(1, 4), "wrong reference base change")
require(reference_U == sp.Rational(1, 18), "wrong marked U-coordinate")
require(reference_V_squared == -sp.Rational(121, 2916), "wrong marked V-coordinate")
reference_target_rhs = sp.factor(reference_U**3 - sp.Rational(3, 4) * reference_U)
require(reference_target_rhs == reference_V_squared, "marked BC point left the reference fibre")

# Exact topology ledger: three punctures on a genus-one curve have free rank
# four.  The pushed vanishing pair meets once; its regular neighborhood is a
# once-punctured torus and the complementary disk contains O,Q+,Q-.  Hence
# mu_O=[delta_0,delta_1]^{-1} mu_-^{-1} mu_+^{-1} (up to orientation), so
# delta_0,delta_1,mu_+,mu_- are a complete generating set.
punctured_surface_free_rank = 2 * 1 + 3 - 1
require(punctured_surface_free_rank == 4, "wrong punctured-torus rank")
monodromy_generator_supports = ("X", "Y", "T+", "T-")
require(len(monodromy_generator_supports) == punctured_surface_free_rank, "incomplete meridian basis")


# ---------------------------------------------------------------------------
# 4. Orbit-merger obstruction, including identity-generator hostiles.
# ---------------------------------------------------------------------------


def merger_capacity(support_size: int) -> int:
    require(support_size == 0 or support_size >= 2, "a permutation cannot have support one")
    return 0 if support_size == 0 else support_size - 1


def maximum_capacity(degree: int, critical_length: int) -> tuple[int, int, int]:
    support_sum_bound = 2 * degree - critical_length
    best = -1
    best_pair = (-1, -1)
    allowed = (0,) + tuple(range(2, degree + 1))
    for support_x in allowed:
        for support_y in allowed:
            if support_x + support_y > support_sum_bound:
                continue
            capacity = merger_capacity(support_x) + merger_capacity(support_y) + 1 + 1
            if capacity > best:
                best = capacity
                best_pair = (support_x, support_y)
    return support_sum_bound, best, degree - 1


generic_capacity = maximum_capacity(16, 19)
secondary_capacity = maximum_capacity(15, 18)
require(generic_capacity == (13, 14, 15), "generic identity edge case escaped")
require(secondary_capacity == (12, 13, 14), "secondary identity edge case escaped")


def cycle_permutation(degree: int, entries: tuple[int, ...]) -> tuple[int, ...]:
    permutation = list(range(degree))
    for left, right in zip(entries, entries[1:] + entries[:1]):
        permutation[left] = right
    return tuple(permutation)


def transposition(degree: int, left: int, right: int) -> tuple[int, ...]:
    permutation = list(range(degree))
    permutation[left], permutation[right] = permutation[right], permutation[left]
    return tuple(permutation)


def orbit_size(generators: tuple[tuple[int, ...], ...]) -> int:
    found = {0}
    active = [0]
    while active:
        point = active.pop()
        for generator in generators:
            image = generator[point]
            if image not in found:
                found.add(image)
                active.append(image)
    return len(found)


# Sharp hostile at the missing one unit of capacity: one handle generator may
# be the identity.  A cycle on n-2 letters plus two adjacent transpositions
# is transitive and attains capacity n-1.  THM-4138 wins only because its
# support bound is one lower in each residual.
for hostile_degree in (15, 16):
    hostile_cycle = cycle_permutation(hostile_degree, tuple(range(hostile_degree - 2)))
    hostile_t1 = transposition(hostile_degree, hostile_degree - 3, hostile_degree - 2)
    hostile_t2 = transposition(hostile_degree, hostile_degree - 2, hostile_degree - 1)
    require(
        orbit_size((hostile_cycle, hostile_t1, hostile_t2)) == hostile_degree,
        "capacity-threshold hostile is not transitive",
    )
    hostile_capacity = (hostile_degree - 3) + 1 + 1
    require(hostile_capacity == hostile_degree - 1, "capacity hostile missed equality")

generic_packet = (7, 5, 3, 2, 2, 1)
secondary_packet = (7, 3, 2, 2, 2, 2, 1)
require(sum(generic_packet) - 4 == 16, "wrong generic finite-BC deletion")
require(sum(secondary_packet) - 4 == 15, "wrong secondary finite-BC deletion")
require(sum(index - 1 for index in generic_packet) == 14, "wrong generic RH total")
require(sum(index - 1 for index in secondary_packet) == 12, "wrong secondary RH total")


semantic_lines = (
    "scope=theta_only_DeltaV_collision_wall_M8_horizontal_residuals",
    "degree_tower=B_to_S:1;S_to_q:2;e1_excluded_by_EqK_equals_O",
    "quadratic_surface=fibres:I2+2I1+IVstar;MW=ZP;heightP=1/6;torsion=0",
    "polynomial_sections=plusminusP,plusminus2P,plusminus3P;pole23=plusminus3P_only",
    "forced_S=V2-(U-a/2)(U+a/4)2;node_q=5a3/64;BC_values=plusminus3P",
    "punctured_fixed_sheet=parallel_push_off_required;generators=X,Y,Tplus,Tminus",
    "orbit_merger=generic_cap14_lt15;secondary_cap13_lt14",
    "verdict=degree16_and_degree15_horizontal_BC_branches_contradict_transitivity",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4138_PRIMARY_AUDIT")
print("scope=theta_only_DeltaV_collision_wall_M8_horizontal_residuals")
print("degree_tower_options=B_to_S,S_to_q:" + str(degree_tower_options))
print("degree_tower_survivor=B_to_S,S_to_q:" + str(surviving_degree_tower))
print("e1_hostile=CONTRADICTION_WITH_THM4120_EqK_equals_O")
print("normalized_surface=y^2=x^3-3x+2+t^2")
print("discriminant=" + str(weierstrass_delta))
print("fibres=I2@0;I1@t^2=-4;IVstar@infinity;root_ranks=A1+E6")
print("infinity_valuations_c4_c6_Delta=" + str(infinity_valuations))
print("MW_rank=" + str(mw_rank) + ";height_P=" + str(height_P) + ";index,torsion=" + str(primitive_pairs[0]))
print("P=(1,t);2P=" + str(expected_P2) + ";3P=" + str(expected_P3))
print("4P_first_denominator_hostile=" + str(expected_P4))
print("denominator_degree_height_rows=" + str(tuple(denominator_degree_rows)))
print("IVstar_3P_reduction=" + str((X3_infinity, Y3_infinity)) + ";additive_parameter=" + str(identity_component_parameter))
print("polynomial_positive_multiples=" + str(polynomial_multiples) + ";pole_pairs=" + str(pole_pairs))
print("forced_horizontal_curve=V^2=(U-a/2)(U+a/4)^2")
print("q_on_S=" + str(q_on_curve) + ";S_node_q=" + str(node_q) + ";node_hessian_det=" + str(node_hessian_det))
print("target_node_contacts=o1:transverse_to_both_branches;o0_node:missed")
print("reference_q=1/4;BC_U=1/18;BC_V^2=-121/2916;standard_delta1_contains_BC_points")
print("fixed_sheet_repair=parallel_push_off_inside_Milnor_annulus")
print("punctured_torus_relation=[delta0,delta1]*muO*muPlus*muMinus=1")
print("monodromy_generators=X,Y,TPlus,TMinus;muO=surface_relation_word;T_supports=2,2")
print("generic_capacity=support_sum,max_capacity,required=" + str(generic_capacity))
print("secondary_capacity=support_sum,max_capacity,required=" + str(secondary_capacity))
print("capacity_hostile=cycle_on_n_minus_2_plus_two_transpositions_attains_n_minus_1")
print("semantic_sha256=" + semantic_sha256)
print("verdict=GENERIC_N16_CONTRADICTION;SECONDARY_N15_CONTRADICTION;DeltaV_WALL_EMPTY_RELATIVE_TO_DEPENDENCIES")
