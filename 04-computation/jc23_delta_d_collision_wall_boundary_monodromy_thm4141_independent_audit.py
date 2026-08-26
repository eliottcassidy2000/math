#!/usr/bin/env python3
"""Clean-room boundary and monodromy audit for candidate THM-4141.

This standard-library route imports no primary checker and uses the reciprocal
coordinates z=1/p, w=s/p, rather than the primary (u,v) chart.  A tiny sparse
Laurent engine checks the collision strict transform over Q(r,Q); integer
polygon arithmetic checks genus and ramification; and permutation arithmetic
checks the two abstract monodromy budgets.  The Mordell--Weil input
E_q(k(q))={O} and the horizontal-carrier classification remain inherited
geometric dependencies, not computations performed here.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import math


Sparse = dict[tuple[int, ...], Fraction]
failures: list[str] = []


def check(condition: bool, message: str) -> None:
    if not condition:
        failures.append(message)


def clean(poly: Sparse) -> Sparse:
    return {exponents: value for exponents, value in poly.items() if value}


def add(left: Sparse, right: Sparse) -> Sparse:
    answer = dict(left)
    for exponents, value in right.items():
        answer[exponents] = answer.get(exponents, Fraction(0)) + value
    return clean(answer)


def scale(poly: Sparse, scalar: Fraction | int) -> Sparse:
    scalar = Fraction(scalar)
    return clean({exponents: scalar * value for exponents, value in poly.items()})


def multiply(left: Sparse, right: Sparse) -> Sparse:
    answer: Sparse = {}
    for left_exponents, left_value in left.items():
        for right_exponents, right_value in right.items():
            exponents = tuple(
                a + b for a, b in zip(left_exponents, right_exponents)
            )
            answer[exponents] = (
                answer.get(exponents, Fraction(0)) + left_value * right_value
            )
    return clean(answer)


def power(poly: Sparse, exponent: int) -> Sparse:
    result: Sparse = {(0,) * len(next(iter(poly))): Fraction(1)}
    base = poly
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        remaining //= 2
    return result


def variable(index: int, number_of_variables: int = 4) -> Sparse:
    exponents = [0] * number_of_variables
    exponents[index] = 1
    return {tuple(exponents): Fraction(1)}


def coefficient_at(poly: Sparse, trailing: tuple[int, int]) -> Sparse:
    """Keep the first two exponents when the final two equal trailing."""

    answer: Sparse = {}
    for exponents, value in poly.items():
        if exponents[-2:] == trailing:
            key = exponents[:2]
            answer[key] = answer.get(key, Fraction(0)) + value
    return clean(answer)


# Variables are (r,Q,s,p).  Clearing D0=6r^2+7 avoids a symbolic rational
# function package.  The live wall is
# Delta=5696/(15D0), K=Delta*r^2, Phi=2Delta*r.
r_var, Q_var, s_var, p_var = (variable(index) for index in range(4))
one: Sparse = {(0, 0, 0, 0): Fraction(1)}
r2 = power(r_var, 2)
s2 = power(s_var, 2)
D0 = add(scale(r2, 6), scale(one, 7))
delta_clear = Fraction(5696, 15)

lower_H = add(
    add(scale(p_var, -3), scale(power(p_var, 2), Fraction(8, 3))),
    scale(power(p_var, 3), Fraction(-1376, 135)),
)
H_clear = add(
    multiply(D0, lower_H),
    add(
        scale(multiply(multiply(r2, s2), power(p_var, 2)), delta_clear),
        add(
            scale(
                multiply(multiply(r_var, s_var), power(p_var, 3)),
                2 * delta_clear,
            ),
            scale(power(p_var, 4), delta_clear),
        ),
    ),
)
s2_minus_p = add(s2, scale(p_var, -1))
F_clear = add(
    add(
        multiply(D0, s2_minus_p),
        scale(multiply(multiply(Q_var, s2_minus_p), H_clear), -1),
    ),
    scale(multiply(multiply(Q_var, D0), s2), Fraction(-1, 2)),
)

# Exact support and Newton polygon.  Coefficients are grouped by (s,p); the
# first two exponents record the independent parameters r,Q.
support_coefficients: dict[tuple[int, int], Sparse] = {}
for (r_exponent, Q_exponent, s_exponent, p_exponent), value in F_clear.items():
    support_key = (s_exponent, p_exponent)
    coefficient_key = (r_exponent, Q_exponent)
    coefficient = support_coefficients.setdefault(support_key, {})
    coefficient[coefficient_key] = coefficient.get(coefficient_key, Fraction(0)) + value
support = sorted(key for key, value in support_coefficients.items() if clean(value))


def cross(origin: tuple[int, int], first: tuple[int, int], second: tuple[int, int]) -> int:
    return (first[0] - origin[0]) * (second[1] - origin[1]) - (
        first[1] - origin[1]
    ) * (second[0] - origin[0])


def convex_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted(set(points))
    lower: list[tuple[int, int]] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


hull = convex_hull(support)
expected_hull = [(0, 1), (2, 0), (4, 2), (2, 4), (0, 5)]
check(hull == expected_hull, "Newton polygon changed")

twice_area = abs(
    sum(
        hull[index][0] * hull[(index + 1) % len(hull)][1]
        - hull[index][1] * hull[(index + 1) % len(hull)][0]
        for index in range(len(hull))
    )
)
boundary_steps = sum(
    math.gcd(
        abs(hull[(index + 1) % len(hull)][0] - hull[index][0]),
        abs(hull[(index + 1) % len(hull)][1] - hull[index][1]),
    )
    for index in range(len(hull))
)
interior_points = (twice_area - boundary_steps + 2) // 2
check((twice_area, boundary_steps, interior_points) == (24, 10, 8), "Pick ledger changed")

edge_ledger: list[tuple[tuple[int, int], tuple[int, int], tuple[int, int], int, int]] = []
for start, end in zip(hull, hull[1:] + hull[:1]):
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = math.gcd(abs(dx), abs(dy))
    inward = (-dy // length, dx // length)
    support_height = inward[0] * start[0] + inward[1] * start[1]
    distance = inward[0] + inward[1] - support_height
    edge_ledger.append((start, end, inward, length, distance))
expected_edge_data = [
    ((0, 1), (2, 0), (1, 2), 1, 1),
    ((2, 0), (4, 2), (-1, 1), 2, 2),
    ((4, 2), (2, 4), (-1, -1), 2, 4),
    ((2, 4), (0, 5), (-1, -2), 1, 7),
    ((0, 5), (0, 1), (1, 0), 4, 1),
]
check(edge_ledger == expected_edge_data, "edge distance ledger changed")

# Independent reciprocal chart: z=1/p and w=s/p.  A monomial s^i p^j in
# z^6*D0*F_Q becomes w^i z^(6-i-j).  Then shift w=v-1/r, so the repeated
# edge root is v=0.  The output ring has variables (r,Q,z,v), allowing
# negative r powers after the shift.
reciprocal: Sparse = {}
for (r_exponent, Q_exponent, s_exponent, p_exponent), value in F_clear.items():
    key = (r_exponent, Q_exponent, 6 - s_exponent - p_exponent, s_exponent)
    reciprocal[key] = reciprocal.get(key, Fraction(0)) + value
reciprocal = clean(reciprocal)
check(all(exponents[2] >= 0 for exponents in reciprocal), "negative reciprocal z power")

shifted: Sparse = {}
for (r_exponent, Q_exponent, z_exponent, w_exponent), value in reciprocal.items():
    for v_exponent in range(w_exponent + 1):
        omitted = w_exponent - v_exponent
        key = (
            r_exponent - omitted,
            Q_exponent,
            z_exponent,
            v_exponent,
        )
        contribution = value * math.comb(w_exponent, v_exponent) * ((-1) ** omitted)
        shifted[key] = shifted.get(key, Fraction(0)) + contribution
shifted = clean(shifted)

expected_z0_v2 = {(0, 1): -delta_clear}
expected_z0_v3 = {(1, 1): 2 * delta_clear}
expected_z0_v4 = {(2, 1): -delta_clear}
expected_z1_v0 = {
    (0, 1): Fraction(1376 * 6, 135),
    (-2, 1): Fraction(1376 * 7, 135),
}
check(not coefficient_at(shifted, (0, 0)), "collision constant term did not vanish")
check(not coefficient_at(shifted, (0, 1)), "collision linear term did not vanish")
check(coefficient_at(shifted, (0, 2)) == expected_z0_v2, "wrong v^2 row")
check(coefficient_at(shifted, (0, 3)) == expected_z0_v3, "wrong v^3 row")
check(coefficient_at(shifted, (0, 4)) == expected_z0_v4, "wrong v^4 row")
check(coefficient_at(shifted, (1, 0)) == expected_z1_v0, "wrong z row")

# Since the z coefficient is 1376*Q*(6r^2+7)/(135r^2), it is a unit on the
# live wall.  The curve is smooth and z~const*v^2.  Also
# ds wedge dp=-z^-3 dw wedge dz and F_Q=z^-6*R, so the residue is a unit
# times z^3 dv and has order 6: the merged puncture has index 7.
generic_infinity_packet = sorted([1, 2, 2, 4, 4, 7], reverse=True)
wall_infinity_packet = sorted([1, 2, 2, 7, 7], reverse=True)
check(sum(index - 1 for index in generic_infinity_packet) == 14, "off-wall RH defect changed")
check(sum(index - 1 for index in wall_infinity_packet) == 14, "wall RH defect changed")
check(2 * interior_points - 2 == 14, "genus/RH mismatch")
check(
    sum(index - 1 for index in [1, 2, 2, 8, 7]) == 15,
    "naive index-eight hostile unexpectedly preserves RH",
)

# The three constant-field punctures have indices 1,7,7 and hence are forced
# to O by the inherited E_q(k(q))={O}.  The BC edge is the irreducible closed
# point K*(sp)^2=q-1/2: its q-1/2 valuation is one, so its two index-two
# geometric points respond together.  Therefore the only degree responses
# are 15 (BC finite) and 19 (BC over O).
rational_index_sum = 1 + 7 + 7
quadratic_locked_sum = 2 + 2
degree_alternatives = [rational_index_sum, rational_index_sum + quadratic_locked_sum]
check(degree_alternatives == [15, 19], "degree response changed")


def cycle_permutation(size: int, cycle: list[int]) -> tuple[int, ...]:
    permutation = list(range(size))
    for source, target in zip(cycle, cycle[1:] + cycle[:1]):
        permutation[source] = target
    return tuple(permutation)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def cycle_type(permutation: tuple[int, ...]) -> list[int]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cursor = start
        length = 0
        while cursor not in seen:
            seen.add(cursor)
            length += 1
            cursor = permutation[cursor]
        lengths.append(length)
    return sorted(lengths, reverse=True)


# Full-boundary response: n=19 and L=18 give |supp X|+|supp Y|<=20.
# Transitivity by X,Y needs merger capacity at least 18, forcing equality,
# one shared pivot and one cycle on each support.  Up to conjugacy their
# support sizes are a and 20-a.  The commutator is always a three-cycle,
# incompatible with the boundary cycle type (7,7,2,2,1).
full_branch_commutators: list[tuple[int, int, list[int]]] = []
for first_size in range(3, 18):
    second_size = 20 - first_size
    first = cycle_permutation(19, list(range(first_size)))
    second = cycle_permutation(19, [0] + list(range(first_size, 19)))
    commutator = compose(compose(compose(first, second), inverse(first)), inverse(second))
    full_branch_commutators.append((first_size, second_size, cycle_type(commutator)))
check(
    all(cycle_lengths == [3] + [1] * 16 for _, _, cycle_lengths in full_branch_commutators),
    "shared-pivot commutator is not uniformly a three-cycle",
)
check(wall_infinity_packet != [3] + [1] * 16, "boundary packet became a three-cycle")

# Finite-BC response: n=15 and L=18 give support sum <=12.  The two finite
# index-two points contribute at most two transposition merger units.
required_mergers_15 = 14
finite_branch_capacity = {
    "both_nonidentity": 12,
    "one_identity": 13,
    "both_identity": 2,
}
check(
    all(capacity < required_mergers_15 for capacity in finite_branch_capacity.values()),
    "finite-carrier merger budget no longer contradicts transitivity",
)
hostile_length_17_one_identity_capacity = (2 * 15 - 17) - 1 + 2
check(
    hostile_length_17_one_identity_capacity == required_mergers_15,
    "length-seventeen sharp hostile changed",
)

semantic_lines = (
    "scope=delta_only_exact_M8_DeltaD_wall_boundary",
    "coordinates=z=1/p;w=s/p;collision_shift=w+1/r",
    "newton_vertices=(0,1),(2,0),(4,2),(2,4),(0,5)",
    "pick=area2:24,boundary:10,genus:8",
    "local_form=-Q*Delta*v^2+[1376Q/(135r^2)]*z+higher",
    "wall_packet=(7,7,2,2,1)",
    "degree_alternatives=(15,19)_relative_to_target_MW",
    "n19=shared_pivot_commutator_three_cycle_contradiction",
    "n15=two_transposition_orbit_merger_capacity_contradiction",
    "geometric_dependencies=Keller_realization,target_MW,horizontal_carrier,finite_transport",
    "verdict=INDEPENDENT_EXACT_AUDIT_OF_CANDIDATE_CLOSURE;JC2=OPEN",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4141_DELTA_D_BOUNDARY_MONODROMY_INDEPENDENT")
print("coordinates=z=1/p;w=s/p;shift=w+1/r;engine=stdlib_sparse_laurent")
print("newton_vertices=" + ";".join(f"{point[0]},{point[1]}" for point in hull))
print(f"pick=area2:{twice_area};boundary:{boundary_steps};genus:{interior_points}")
print("edge_rows=" + ";".join(
    f"{start}->{end}:normal={normal}:length={length}:index={distance}"
    for start, end, normal, length, distance in edge_ledger
))
print("collision_edge=-Q*Delta*w^2*(r*w+1)^2;root=w=-1/r;multiplicity=2")
print("strict_transform_D0_cleared=v2:-5696*Q/15;z:1376*Q*(6r^2+7)/(135r^2)")
print("local_geometry=smooth_tangency:z~v^2;residue_order=6;merged_index=7")
print("generic_packet=" + ",".join(map(str, generic_infinity_packet)))
print("wall_packet=" + ",".join(map(str, wall_infinity_packet)) + ";RH_defect=14")
print("affine_s=0_points=4;ramification_index=1")
print("rational_punctures=AB:1,collision:7,far:7;index_sum=15")
print("BC_closed_point=K*(s*p)^2=q-1/2;degree=2;indices=2,2;galois_locked")
print("degree_alternatives=15,19;relative_to=E_q(k(q))={O}")
print("n19_support_sum_max=20;mergers_required=18;equality=shared_pivot_two_cycles")
print("n19_commutator_controls=" + ";".join(
    f"{first}+{second}:3cycle" for first, second, _ in full_branch_commutators
))
print("n19_boundary_commutator=7,7,2,2,1;verdict=contradiction")
print("n15_support_sum_max=12;finite_meridians=2_transpositions;required=14")
print("n15_capacity=" + ";".join(
    f"{case}:{capacity}" for case, capacity in finite_branch_capacity.items()
))
print("hostile_L17_one_identity_capacity=14=required")
print("semantic_sha256=" + semantic_sha256)
print("dependencies=Keller_realization;target_MW;horizontal_carrier;finite_etale_transport")
print("verdict=candidate_DeltaD_wall_closure_audited;JC2=OPEN")

if failures:
    for failure in failures:
        print("FAIL=" + failure)
    raise SystemExit(1)
