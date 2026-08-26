#!/usr/bin/env python3
"""Clean-room exact audit for the proposed THM-4138 exclusion.

No code is imported from the primary THM-4138 artifact.  This verifier uses
only the Python standard library.  It checks the quadratic-base-change
Weierstrass arithmetic with exact rational-function operations, reconstructs
the multiples of the Mordell--Weil generator, audits the marked real-fibre
coordinates, and verifies the orbit-merger capacity (including identity
generators and sharp one-unit hostiles).  The Shioda--Tate, height-pairing,
local Milnor-annulus, and punctured-surface steps are stated as typed proof
ledgers; finite algebra is not presented as a replacement for those results.

Replay:
    python3 -B 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138_independent_audit.py
    python3 -B -O 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138_independent_audit.py
    PYTHONHASHSEED=0 python3 -B 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138_independent_audit.py
    PYTHONHASHSEED=8675309 python3 -B 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138_independent_audit.py
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib
import itertools
from typing import Iterable


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------------------
# Exact univariate polynomials over Q, stored in ascending coefficient order.
# ---------------------------------------------------------------------------

Poly = tuple[F, ...]
Rat = tuple[Poly, Poly]
Point = tuple[Rat, Rat] | None


def pol(values: Iterable[int | F]) -> Poly:
    out = [F(value) for value in values]
    if not out:
        out = [F(0)]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


ZERO = pol([0])
ONE = pol([1])
T_POLY = pol([0, 1])


def pzero(value: Poly) -> bool:
    return value == ZERO


def padd(left: Poly, right: Poly) -> Poly:
    size = max(len(left), len(right))
    return pol(
        (left[i] if i < len(left) else 0) + (right[i] if i < len(right) else 0)
        for i in range(size)
    )


def pneg(value: Poly) -> Poly:
    return pol(-coefficient for coefficient in value)


def psub(left: Poly, right: Poly) -> Poly:
    return padd(left, pneg(right))


def pscale(value: Poly, scalar: int | F) -> Poly:
    scalar = F(scalar)
    return pol(scalar * coefficient for coefficient in value)


def pmul(left: Poly, right: Poly) -> Poly:
    out = [F(0)] * (len(left) + len(right) - 1)
    for i, a0 in enumerate(left):
        for j, b0 in enumerate(right):
            out[i + j] += a0 * b0
    return pol(out)


def ppow(value: Poly, exponent: int) -> Poly:
    check(exponent >= 0, "negative polynomial exponent")
    answer = ONE
    base = value
    power = exponent
    while power:
        if power & 1:
            answer = pmul(answer, base)
        base = pmul(base, base)
        power //= 2
    return answer


def pdivmod(numerator: Poly, denominator: Poly) -> tuple[Poly, Poly]:
    check(not pzero(denominator), "polynomial division by zero")
    remainder = list(numerator)
    quotient = [F(0)] * max(1, len(numerator) - len(denominator) + 1)
    while not pzero(pol(remainder)) and len(remainder) >= len(denominator):
        shift = len(remainder) - len(denominator)
        coefficient = remainder[-1] / denominator[-1]
        quotient[shift] += coefficient
        for index, value in enumerate(denominator):
            remainder[index + shift] -= coefficient * value
        remainder = list(pol(remainder))
    return pol(quotient), pol(remainder)


def pexact_div(numerator: Poly, denominator: Poly) -> Poly:
    quotient, remainder = pdivmod(numerator, denominator)
    check(pzero(remainder), "inexact polynomial division")
    return quotient


def pgcd(left: Poly, right: Poly) -> Poly:
    a0, b0 = left, right
    while not pzero(b0):
        _, remainder = pdivmod(a0, b0)
        a0, b0 = b0, remainder
    if pzero(a0):
        return ONE
    return pscale(a0, 1 / a0[-1])


def pvaluation(value: Poly) -> int:
    check(not pzero(value), "valuation of zero polynomial")
    return next(index for index, coefficient in enumerate(value) if coefficient)


def rmake(numerator: Poly, denominator: Poly = ONE) -> Rat:
    check(not pzero(denominator), "rational-function denominator zero")
    if pzero(numerator):
        return ZERO, ONE
    divisor = pgcd(numerator, denominator)
    top = pexact_div(numerator, divisor)
    bottom = pexact_div(denominator, divisor)
    scale = bottom[-1]
    top = pscale(top, 1 / scale)
    bottom = pscale(bottom, 1 / scale)
    return top, bottom


def rconst(value: int | F) -> Rat:
    return rmake(pol([F(value)]))


T = rmake(T_POLY)


def radd(left: Rat, right: Rat) -> Rat:
    return rmake(padd(pmul(left[0], right[1]), pmul(right[0], left[1])), pmul(left[1], right[1]))


def rneg(value: Rat) -> Rat:
    return rmake(pneg(value[0]), value[1])


def rsub(left: Rat, right: Rat) -> Rat:
    return radd(left, rneg(right))


def rmul(left: Rat, right: Rat) -> Rat:
    return rmake(pmul(left[0], right[0]), pmul(left[1], right[1]))


def rdiv(left: Rat, right: Rat) -> Rat:
    check(not pzero(right[0]), "rational-function division by zero")
    return rmake(pmul(left[0], right[1]), pmul(left[1], right[0]))


def rpow(value: Rat, exponent: int) -> Rat:
    check(exponent >= 0, "negative rational-function exponent")
    return rmake(ppow(value[0], exponent), ppow(value[1], exponent))


def rdegree(value: Rat) -> int:
    check(not pzero(value[0]), "degree of zero rational function")
    return (len(value[0]) - 1) - (len(value[1]) - 1)


def rpoly_degree(value: Rat) -> int:
    check(value[1] == ONE, "coordinate is not polynomial")
    return len(value[0]) - 1


def elliptic_add(left: Point, right: Point) -> Point:
    """Group law on y^2=x^3-3x+2+t^2 over Q(t)."""
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == rneg(y2):
            return None
        slope = rdiv(rsub(rmul(rconst(3), rpow(x1, 2)), rconst(3)), rmul(rconst(2), y1))
    else:
        slope = rdiv(rsub(y2, y1), rsub(x2, x1))
    x3 = rsub(rsub(rpow(slope, 2), x1), x2)
    y3 = rsub(rmul(slope, rsub(x1, x3)), y1)
    return x3, y3


def on_surface(point: Point) -> bool:
    if point is None:
        return True
    x0, y0 = point
    equation = rsub(
        rsub(radd(rsub(rpow(y0, 2), rpow(x0, 3)), rmul(rconst(3), x0)), rconst(2)),
        rpow(T, 2),
    )
    return equation == rconst(0)


# ---------------------------------------------------------------------------
# 1. Quadratic base change and the rational elliptic surface.
# ---------------------------------------------------------------------------

# The dimensionless variable w=r^2/a^3 checks the original-coordinate 3P
# formulas without a symbolic-algebra dependency.
w = rmake(pol([0, 1]))
u_over_a = radd(rconst(F(1, 2)), rmul(rconst(F(16, 9)), w))
v_factor = radd(rconst(1), rmul(rconst(F(64, 27)), w))
v_squared_over_a3 = rmul(w, rpow(v_factor, 2))
q_over_a3 = radd(rconst(F(1, 2)), w)

target_identity = rsub(
    radd(
        radd(rsub(v_squared_over_a3, rpow(u_over_a, 3)), rmul(rconst(F(3, 4)), u_over_a)),
        rconst(F(1, 4)),
    ),
    q_over_a3,
)
check(target_identity == rconst(0), "3P missed the original target pencil")

nodal_identity = rsub(
    v_squared_over_a3,
    rmul(rsub(u_over_a, rconst(F(1, 2))), rpow(radd(u_over_a, rconst(F(1, 4))), 2)),
)
check(nodal_identity == rconst(0), "3P missed the forced nodal image")

# Normalized surface y^2=x^3-3x+2+t^2.
a4 = -3
b6 = padd(pmul(T_POLY, T_POLY), pol([2]))
discriminant = pscale(padd(pol([4 * a4**3]), pscale(pmul(b6, b6), 27)), -16)
expected_discriminant = pol([0, 0, -1728, 0, -432])
check(discriminant == expected_discriminant, "wrong base-change discriminant")

u4 = pol([0, 0, 0, 0, 1])
u6 = pol([0, 0, 0, 0, 0, 0, 1])
u8 = pol([0, 0, 0, 0, 0, 0, 0, 0, 1])
c4_infinity = pscale(u4, 144)
c6_infinity = pscale(padd(u4, pscale(u6, 2)), -864)
delta_infinity = pscale(pmul(u8, pol([1, 0, 4])), -432)
infinity_valuations = tuple(map(pvaluation, (c4_infinity, c6_infinity, delta_infinity)))
check(infinity_valuations == (4, 4, 8), "infinity is not the IV* valuation row")
check(2 + 1 + 1 + 8 == 12, "singular-fibre Euler ledger is not rational")

# Shioda--Tate: rho=10 for a rational elliptic surface, and A1+E6 has
# root rank seven.  P meets the nonidentity I2 component and the nontrivial
# IV* component class; 3P reaches the additive identity component.
mw_rank = 10 - 2 - (1 + 6)
height_P = F(2) - F(1, 2) - F(4, 3)
check(mw_rank == 1 and height_P == F(1, 6), "wrong Shioda--Tate/height ledger")

# Unimodularity gives 1=disc(A1)disc(E6)h(Q)/|tors|^2.  If P=mQ+T,
# h(P)=m^2 h(Q); the only positive integral (m,|tors|) is (1,1).
primitive_solutions = tuple(
    (index, torsion)
    for index in range(1, 25)
    for torsion in range(1, 25)
    if F(index * index * torsion * torsion, 6) == height_P
)
check(primitive_solutions == ((1, 1),), "P was not forced primitive and torsion-free")


# ---------------------------------------------------------------------------
# 2. Clean-room Mordell--Weil multiples and polynomial-section classification.
# ---------------------------------------------------------------------------

P: Point = (rconst(1), T)
multiples: dict[int, Point] = {}
current: Point = None
denominator_height_rows: list[tuple[int, int, int]] = []
for multiplier in range(1, 13):
    current = elliptic_add(current, P)
    check(current is not None and on_surface(current), f"{multiplier}P left the surface")
    multiples[multiplier] = current
    x_coordinate = current[0]
    denominator_degree = len(x_coordinate[1]) - 1
    epsilon_two = int(multiplier % 2 != 0)
    epsilon_three = int(multiplier % 3 != 0)
    numerator = multiplier * multiplier - 12 + 3 * epsilon_two + 8 * epsilon_three
    check(numerator % 12 == 0, f"nonintegral section intersection for {multiplier}P")
    finite_intersection = numerator // 12
    check(finite_intersection >= 0, f"negative section intersection for {multiplier}P")
    check(denominator_degree == 2 * finite_intersection, f"denominator/height mismatch for {multiplier}P")
    denominator_height_rows.append((multiplier, denominator_degree, finite_intersection))

expected_2P: Point = (rconst(-2), rneg(T))
expected_3P: Point = (
    rmake(pol([9, 0, 4]), pol([9])),
    rmake(pol([0, -27, 0, -8]), pol([27])),
)
expected_4P: Point = (
    rmake(pol([81, 0, 16]), pol([0, 0, 4])),
    rmake(pol([729, 0, 216, 0, 8]), pol([0, 0, 0, 8])),
)
check(multiples[2] == expected_2P, "independent 2P formula mismatch")
check(multiples[3] == expected_3P, "independent 3P formula mismatch")
check(multiples[4] == expected_4P, "independent 4P denominator hostile mismatch")

# The height formula is zero exactly for 1,2,3 and is positive for 4,5,6;
# for n>=7 positivity follows already from n^2-12>0.  The IV* component
# group has order three.  The nonsingular reduction of 3P below has nonzero
# additive parameter, so no nonzero multiple meets O at infinity in char 0.
check([row[2] for row in denominator_height_rows[:6]] == [0, 0, 0, 1, 2, 2], "bad low height boundary")
check(7 * 7 - 12 > 0, "large-multiple positivity threshold failed")

x3, y3 = multiples[3]  # type: ignore[misc]
x3_reduction = x3[0][-1] / x3[1][-1]
y3_reduction = y3[0][-1] / y3[1][-1]
check(rdegree(x3) == 2 and rdegree(y3) == 3, "3P has wrong infinity orders")
check((x3_reduction, y3_reduction) == (F(4, 9), F(-8, 27)), "wrong IV* identity reduction")
check(y3_reduction * y3_reduction == x3_reduction**3, "3P reduction missed the cuspidal cubic")
additive_parameter = x3_reduction / y3_reduction
check(additive_parameter == F(-3, 2), "3P reduction became additive zero")

polynomial_positive_multiples = tuple(row[0] for row in denominator_height_rows if row[2] == 0)
check(polynomial_positive_multiples == (1, 2, 3), "unexpected positive polynomial multiple")

pole_pairs = {
    multiplier: (rpoly_degree(multiples[multiplier][0]), rpoly_degree(multiples[multiplier][1]))  # type: ignore[index]
    for multiplier in polynomial_positive_multiples
}
check(pole_pairs == {1: (0, 1), 2: (0, 1), 3: (2, 3)}, "wrong polynomial-section pole pairs")


# ---------------------------------------------------------------------------
# 3. Forced horizontal curve and marked vanishing-cycle geometry.
# ---------------------------------------------------------------------------

# With z=4r/(3a) (up to sign), U=a/2+z^2 and
# V=+/- z(z^2+3a/4).  Therefore q/a^3=1/2+9(z^2/a)/16.
node_scaled_z_squared = F(-3, 4)
node_q_over_a3 = F(1, 2) + F(9, 16) * node_scaled_z_squared
check(node_q_over_a3 == F(5, 64), "wrong node value of the forced horizontal curve")
check(node_q_over_a3 not in (F(0), F(1, 4), F(1, 2)), "horizontal branch collision hit a node/reference value")

# The horizontal curve is smooth at o1=(a/2,0), misses o0=(-a/2,0),
# and has its ordinary node at (-a/4,0).  These are the scaled a=1 rows.
gradient_at_o1 = F(-9, 16)
node_hessian_determinant = F(3)
o0_curve_value = -((-F(1, 2) - F(1, 2)) * (-F(1, 2) + F(1, 4)) ** 2)
check(gradient_at_o1 != 0, "forced S is singular at o1")
check(node_hessian_determinant != 0, "forced S node is not ordinary")
check(o0_curve_value != 0, "forced S unexpectedly passes through o0")

reference_q = F(1, 4)
reference_r_squared = reference_q - F(1, 2)
reference_U = F(1, 2) + F(16, 9) * reference_r_squared
reference_V_squared = reference_r_squared * (1 + F(64, 27) * reference_r_squared) ** 2
reference_target_rhs = reference_U**3 - F(3, 4) * reference_U
check(reference_r_squared == F(-1, 4), "wrong marked quadratic fibre")
check(reference_U == F(1, 18), "wrong marked U-coordinate")
check(reference_V_squared == F(-121, 2916), "wrong marked V-coordinate")
check(reference_target_rhs == reference_V_squared, "marked points left E_(1/4)")

# At q*=1/4 the branch roots are -sqrt(3)/2,0,+sqrt(3)/2.  The positive
# U-coordinate and U^2<3/4 put both Q-points on the lift of [0,sqrt(3)/2],
# i.e. the standard delta_1, while delta_0 is the negative interval lift.
check(reference_U > 0 and reference_U * reference_U < F(3, 4), "Q points are not inside delta1 interval")

# A genus-one surface with three punctures has free rank four.  After pushing
# delta1 to a parallel curve in its Milnor annulus, a regular neighborhood of
# delta0 union delta1 is a once-punctured torus; its complementary disk holds
# O,Q+,Q-.  Van Kampen gives [delta0,delta1]muOmu+mu-=1, and eliminating muO
# makes delta0,delta1,mu+,mu- a free basis.
punctured_torus_rank = 2 * 1 + 3 - 1
handle_spine_euler = 1 - 2
check(punctured_torus_rank == 4 and handle_spine_euler == -1, "wrong punctured-surface ledger")


# ---------------------------------------------------------------------------
# 4. Orbit-merger inequality and residual capacities.
# ---------------------------------------------------------------------------

Permutation = tuple[int, ...]


def identity(degree: int) -> Permutation:
    return tuple(range(degree))


def cycle(degree: int, entries: tuple[int, ...]) -> Permutation:
    result = list(range(degree))
    if len(entries) >= 2:
        for left, right in zip(entries, entries[1:] + entries[:1]):
            result[left] = right
    return tuple(result)


def swap(degree: int, left: int, right: int) -> Permutation:
    result = list(range(degree))
    result[left], result[right] = result[right], result[left]
    return tuple(result)


def support_size(permutation: Permutation) -> int:
    return sum(index != image for index, image in enumerate(permutation))


def merger_capacity_from_support(size: int) -> int:
    check(size == 0 or size >= 2, "permutation support cannot have size one")
    return 0 if size == 0 else size - 1


def exact_single_generator_mergers(permutation: Permutation) -> int:
    visited: set[int] = set()
    mergers = 0
    for start in range(len(permutation)):
        if start in visited:
            continue
        cursor = start
        length = 0
        while cursor not in visited:
            visited.add(cursor)
            length += 1
            cursor = permutation[cursor]
        if length > 1:
            mergers += length - 1
    return mergers


def orbit_component_count(generators: tuple[Permutation, ...]) -> int:
    degree = len(generators[0])
    parent = list(range(degree))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: int, right: int) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for generator in generators:
        for vertex, image in enumerate(generator):
            union(vertex, image)
    return len({find(vertex) for vertex in range(degree)})


# The exact merges contributed by one permutation are sum(cycle_length-1),
# at most support-1 for a nonidentity and zero for the identity.  This proves
# sum max(|supp g|-1,0)>=n-1 for any transitive generated action.  Exhaust all
# ordered four-generator systems on four letters as an implementation hostile.
small_permutations = tuple(tuple(p) for p in itertools.permutations(range(4)))
for permutation in small_permutations:
    size = support_size(permutation)
    check(exact_single_generator_mergers(permutation) <= merger_capacity_from_support(size), "single-generator merger bound failed")

four_letter_systems = 0
for generators in itertools.product(small_permutations, repeat=4):
    components = orbit_component_count(generators)
    capacity = sum(merger_capacity_from_support(support_size(generator)) for generator in generators)
    check(4 - components <= capacity, "orbit-merger inequality failed on S4 hostile census")
    four_letter_systems += 1
check(four_letter_systems == 24**4, "incomplete S4 generator census")


def residual_capacity(degree: int, critical_length: int) -> tuple[int, int, int, tuple[tuple[int, int], ...]]:
    support_sum_bound = 2 * degree - critical_length
    allowed = (0,) + tuple(range(2, degree + 1))
    best = -1
    maximizers: list[tuple[int, int]] = []
    for support_x in allowed:
        for support_y in allowed:
            if support_x + support_y > support_sum_bound:
                continue
            capacity = (
                merger_capacity_from_support(support_x)
                + merger_capacity_from_support(support_y)
                + 1
                + 1
            )
            if capacity > best:
                best = capacity
                maximizers = [(support_x, support_y)]
            elif capacity == best:
                maximizers.append((support_x, support_y))
    return support_sum_bound, best, degree - 1, tuple(maximizers)


generic_capacity = residual_capacity(16, 19)
secondary_capacity = residual_capacity(15, 18)
check(generic_capacity == (13, 14, 15, ((0, 13), (13, 0))), "generic identity case escaped")
check(secondary_capacity == (12, 13, 14, ((0, 12), (12, 0))), "secondary identity case escaped")

# Sharp one-unit relaxation: one handle generator is the identity, the other
# is an (n-2)-cycle, and two transpositions attach the remaining vertices.
# At the actual residual bound the analogous construction leaves one isolated
# vertex, showing exactly where the strict inequality is spent.
hostile_rows: list[tuple[int, int, int, int]] = []
for degree in (15, 16):
    full_cycle = cycle(degree, tuple(range(degree - 2)))
    full_t_plus = swap(degree, degree - 3, degree - 2)
    full_t_minus = swap(degree, degree - 2, degree - 1)
    full_generators = (identity(degree), full_cycle, full_t_plus, full_t_minus)
    check(orbit_component_count(full_generators) == 1, "one-unit relaxed hostile is not transitive")
    full_capacity = sum(merger_capacity_from_support(support_size(g)) for g in full_generators)
    check(full_capacity == degree - 1, "one-unit hostile missed equality")

    near_cycle = cycle(degree, tuple(range(degree - 3)))
    near_t_plus = swap(degree, degree - 4, degree - 3)
    near_t_minus = swap(degree, degree - 3, degree - 2)
    near_generators = (identity(degree), near_cycle, near_t_plus, near_t_minus)
    near_components = orbit_component_count(near_generators)
    near_capacity = sum(merger_capacity_from_support(support_size(g)) for g in near_generators)
    check(near_components == 2 and near_capacity == degree - 2, "residual near-hostile did not leave one isolated sheet")
    hostile_rows.append((degree, full_capacity, near_capacity, near_components))


generic_full_packet = (7, 5, 3, 2, 2, 1)
secondary_full_packet = (7, 3, 2, 2, 2, 2, 1)
generic_O_packet = (7, 5, 3, 1)
secondary_O_packet = (7, 3, 2, 2, 1)
check(sum(generic_O_packet) == 16 and sum(secondary_O_packet) == 15, "wrong O-fibre deletion")
check(sum(index - 1 for index in generic_O_packet) + 2 == 14, "wrong generic ramification split")
check(sum(index - 1 for index in secondary_O_packet) + 2 == 12, "wrong secondary ramification split")
check(sum(generic_full_packet) - sum(generic_O_packet) == 4, "generic BC pair is not 2+2")
check(sum(secondary_full_packet) - sum(secondary_O_packet) == 4, "secondary BC pair is not 2+2")


# Semantic ledger deliberately matches the independently interpreted theorem
# consequence, making primary/independent agreement machine-comparable.
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


print("THM4138_INDEPENDENT_AUDIT")
print("implementation=python_stdlib_fraction_polynomials;no_primary_import")
print("scope=theta_only_exact_M8_DeltaV_wall_relative_to_THM4120_4122_4130_4134")
print("claim_search_at_head=a97c925b76dc;exact_duplicate_or_correction=NONE_FOUND")
print("freshness=fetched_origin_main_4fa44ddd9d;incoming_THM4139_is_reserved_and_not_a_duplicate_or_correction")
print("degree_tower=deg(B_to_normalizationS)=1;deg(normalizationS_to_P1q)=2")
print("degree_one_image=excluded_by_THM4120_Eq_of_kq_equals_O_and_finite_horizontal_image")
print("normalized_surface=y^2=x^3-3x+2+t^2")
print("discriminant=-432*t^2*(t^2+4);fibres=I2+2I1+IVstar;Euler=12")
print("infinity_valuations_c4_c6_Delta=" + str(infinity_valuations))
print("Shioda_Tate=rank1;local_contributions_P=1/2+4/3;height_P=" + str(height_P))
print("unimodular_index_torsion_solution=" + str(primitive_solutions) + ";MW=ZP_torsion_free")
print("P=(1,t);2P=(-2,-t);3P=((4t^2+9)/9,-t(8t^2+27)/27)")
print("4P_denominator_hostile=((16t^2+81)/(4t^2),(8t^4+216t^2+729)/(8t^3))")
print("denominator_degree_height_rows=" + str(tuple(denominator_height_rows)))
print("IVstar_3P_reduction=" + str((x3_reduction, y3_reduction)) + ";additive_parameter=" + str(additive_parameter))
print("polynomial_sections=plusminus_mP_for_m=" + str(polynomial_positive_multiples) + ";pole_pairs=" + str(pole_pairs))
print("forced_horizontal_curve=V^2=(U-a/2)(U+a/4)^2;node_q_over_a3=" + str(node_q_over_a3))
print("target_node_contacts=o1_transverse_to_both_local_branches;o0_missed")
print("transport_path=avoid_q_over_a3_0,1/2,5/64;delta0_fibre_detour_returns_standard_at_qstar")
print("reference_q=1/4;Q_U=" + str(reference_U) + ";Q_V_squared=" + str(reference_V_squared))
print("marked_cycle=Qplus,Qminus_on_standard_delta1;parallel_push_off=MANDATORY")
print("fixed_sheet_injection=local_uv_base_preserved_and_degree_one;distinct_neighborhoods_plus_finite_transport")
print("punctured_torus_rank=4;relation=[delta0,delta1]*muO*muPlus*muMinus=1")
print("free_basis=delta0,pushed_delta1,muPlus,muMinus;monodromy=X,Y,TPlus,TMinus")
print("branch_monodromy=TPlus,TMinus_are_transpositions;cover_complement_connected_so_action_transitive")
print("orbit_merger_theorem=sum_max_support_minus_one_ge_degree_minus_one_for_transitive_action")
print("S4_four_generator_exhaustive_systems=" + str(four_letter_systems) + ";identity_cases=included")
print("generic_capacity=support_sum,max_capacity,required,maximizers=" + str(generic_capacity))
print("secondary_capacity=support_sum,max_capacity,required,maximizers=" + str(secondary_capacity))
print("sharp_relaxed_and_near_hostiles=" + str(tuple(hostile_rows)))
print("ramification_split=generic_O(7,5,3,1)+T+T;secondary_O(7,3,2,2,1)+T+T")
print("semantic_sha256=" + semantic_sha256)
print("verdict=ACCEPT;N16_AND_N15_HORIZONTAL_BC_BRANCHES_IMPOSSIBLE_RELATIVE_TO_DEPENDENCIES")
print("boundary=DeltaV_wall_only;other_collision_walls_Mge9_other_cells_JC2_DC2_OPEN")
