#!/usr/bin/env python3
"""Exact finite and symbolic audit for THM-4138.

The certificate checks the two Delta_V wall ledgers, the carrier
ramification-index budget, the sharp orbit-merger bound, and minimal hostile
changes to either numerical margin. It also independently verifies that the
polynomial section displayed in THM-4134 is 3P on the quadratic base change
and that it passes through the nodal target fibre.

No Python ``assert`` is used, so optimized mode retains every truth gate.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def identity(degree: int) -> tuple[int, ...]:
    return tuple(range(degree))


def cycle(degree: int, entries: tuple[int, ...]) -> tuple[int, ...]:
    require(len(entries) >= 2, "cycle is too short")
    require(len(entries) == len(set(entries)), "cycle repeats a letter")
    require(all(0 <= entry < degree for entry in entries), "cycle letter out of range")
    result = list(range(degree))
    for left, right in zip(entries, entries[1:] + entries[:1]):
        result[left] = right
    return tuple(result)


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    require(len(left) == len(right), "degree mismatch")
    return tuple(left[right[index]] for index in range(len(left)))


def support(permutation: tuple[int, ...]) -> set[int]:
    return {point for point, image in enumerate(permutation) if point != image}


def cycles(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(len(permutation)))
    answer: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit: list[int] = []
        point = start
        while point in unseen:
            unseen.remove(point)
            orbit.append(point)
            point = permutation[point]
        answer.append(tuple(orbit))
    return tuple(answer)


def permutation_index(permutation: tuple[int, ...]) -> int:
    return len(permutation) - len(cycles(permutation))


def generated_orbits(generators: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    degree = len(generators[0])
    require(all(len(generator) == degree for generator in generators), "degree mismatch")
    unseen = set(range(degree))
    answer: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        found = {start}
        active = [start]
        while active:
            point = active.pop()
            for generator in generators:
                image = generator[point]
                if image not in found:
                    found.add(image)
                    active.append(image)
        unseen -= found
        answer.append(tuple(sorted(found)))
    return tuple(sorted(answer, key=lambda orbit: orbit[0]))


def transposition(degree: int, left: int, right: int) -> tuple[int, ...]:
    return cycle(degree, (left, right))


def audit_wall(
    *,
    label: str,
    degree: int,
    critical_length: int,
    genus: int,
    packet: tuple[int, ...],
    origin_packet: tuple[int, ...],
) -> str:
    carrier_packet = (2, 2)
    require(critical_length == degree + 3, label + ": critical offset changed")
    require(sum(origin_packet) == degree, label + ": origin degree changed")
    require(sorted(packet) == sorted(origin_packet + carrier_packet), label + ": packet split changed")
    origin_defect = sum(entry - 1 for entry in origin_packet)
    carrier_defect = sum(entry - 1 for entry in carrier_packet)
    require(carrier_defect == 2, label + ": carrier defect changed")
    require(origin_defect + carrier_defect == 2 * genus - 2, label + ": RH not saturated")

    # Sharp positive control: total index n-2 and exactly two orbits.
    x = cycle(degree, tuple(range(degree - 3)))
    y = identity(degree)
    z1 = transposition(degree, degree - 4, degree - 3)
    z2 = transposition(degree, degree - 3, degree - 2)
    fix_sum = (degree - len(support(x))) + (degree - len(support(y)))
    total_index = sum(permutation_index(g) for g in (x, y, z1, z2))
    orbit_packet = generated_orbits((x, y, z1, z2))
    require(fix_sum == critical_length, label + ": fixed-sheet sharp control changed")
    require(total_index == degree - 2, label + ": sharp total index changed")
    require(len(orbit_packet) == 2, label + ": sharp control did not leave two orbits")

    # Same carrier defect when both simple branch points share one target value.
    collided = compose(
        transposition(degree, degree - 4, degree - 3),
        transposition(degree, degree - 2, degree - 1),
    )
    require(permutation_index(collided) == 2, label + ": collided carrier index changed")
    require(len(support(collided)) == 4, label + ": collided carrier is not (2,2)")
    require(len(generated_orbits((x, y, collided))) >= 2, label + ": collided carrier became transitive")

    # Hostile 1: losing one fixed sheet makes transitivity exactly possible.
    hostile_x = cycle(degree, tuple(range(degree - 2)))
    hostile_z1 = transposition(degree, degree - 3, degree - 2)
    hostile_z2 = transposition(degree, degree - 2, degree - 1)
    hostile_fixed = (degree - len(support(hostile_x))) + degree
    hostile_generators = (hostile_x, y, hostile_z1, hostile_z2)
    require(hostile_fixed == degree + 2, label + ": n+2 hostile fixed sum changed")
    require(
        sum(permutation_index(g) for g in hostile_generators) == degree - 1,
        label + ": n+2 hostile index changed",
    )
    require(len(generated_orbits(hostile_generators)) == 1, label + ": n+2 hostile not transitive")

    # Hostile 2: one extra carrier transposition also reaches transitivity.
    hostile_z3 = transposition(degree, degree - 2, degree - 1)
    defect3_generators = (x, y, z1, z2, hostile_z3)
    require(
        sum(permutation_index(g) for g in defect3_generators) == degree - 1,
        label + ": defect-three hostile index changed",
    )
    require(len(generated_orbits(defect3_generators)) == 1, label + ": defect-three hostile not transitive")

    # Scalar bound for every feasible node split.
    feasible_splits = 0
    for r0 in range(2, degree + 1):
        r1 = critical_length - r0
        if not 2 <= r1 <= degree:
            continue
        max_support_sum = 2 * degree - r0 - r1
        require(max_support_sum == degree - 3, label + ": split support budget changed")
        require(max_support_sum - 1 + carrier_defect <= degree - 2, label + ": split escaped index gate")
        feasible_splits += 1
    require(feasible_splits > 0, label + ": no feasible node split")

    return (
        label
        + ":n="
        + str(degree)
        + ";L="
        + str(critical_length)
        + ";origin_defect="
        + str(origin_defect)
        + ";carrier_defect=2;max_total_index="
        + str(degree - 2)
        + ";min_orbits=2;splits="
        + str(feasible_splits)
    )


generic_line = audit_wall(
    label="generic",
    degree=16,
    critical_length=19,
    genus=8,
    packet=(7, 5, 3, 2, 2, 1),
    origin_packet=(7, 5, 3, 1),
)
secondary_line = audit_wall(
    label="secondary",
    degree=15,
    critical_length=18,
    genus=7,
    packet=(7, 3, 2, 2, 2, 2, 1),
    origin_packet=(7, 3, 2, 2, 1),
)

# THM-4134's explicit section is exactly 3P.
a, rho, U, V = sp.symbols("a rho U V", nonzero=True)
q = a**3 / 2 + rho**2
weierstrass_a4 = -sp.Rational(3, 4) * a**2


def elliptic_double(point: tuple[sp.Expr, sp.Expr]) -> tuple[sp.Expr, sp.Expr]:
    x_coord, y_coord = point
    slope = sp.cancel((3 * x_coord**2 + weierstrass_a4) / (2 * y_coord))
    x_double = sp.factor(slope**2 - 2 * x_coord)
    y_double = sp.factor(slope * (x_coord - x_double) - y_coord)
    return x_double, y_double


def elliptic_add(
    left: tuple[sp.Expr, sp.Expr], right: tuple[sp.Expr, sp.Expr]
) -> tuple[sp.Expr, sp.Expr]:
    x_left, y_left = left
    x_right, y_right = right
    slope = sp.cancel((y_right - y_left) / (x_right - x_left))
    x_sum = sp.factor(slope**2 - x_left - x_right)
    y_sum = sp.factor(slope * (x_left - x_sum) - y_left)
    return x_sum, y_sum


P = (a / 2, rho)
twice_P = elliptic_double(P)
thrice_P = elliptic_add(twice_P, P)
expected_twice = (-a, -rho)
expected_thrice = (
    a / 2 + 16 * rho**2 / (9 * a**2),
    -rho - 64 * rho**3 / (27 * a**3),
)
require(all(sp.factor(left - right) == 0 for left, right in zip(twice_P, expected_twice)), "2P identity failed")
require(all(sp.factor(left - right) == 0 for left, right in zip(thrice_P, expected_thrice)), "3P identity failed")


def target_E(point: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    x_coord, y_coord = point
    return sp.factor(y_coord**2 - x_coord**3 + sp.Rational(3, 4) * a**2 * x_coord + a**3 / 4)


require(sp.factor(target_E(P) - q) == 0, "P misses the base-changed surface")
require(sp.factor(target_E(thrice_P) - q) == 0, "3P misses the base-changed surface")
require(sp.factor(thrice_P[0].subs(rho, 0) - a / 2) == 0, "3P misses nodal U")
require(sp.factor(thrice_P[1].subs(rho, 0)) == 0, "3P misses nodal V")

xi = sp.symbols("xi")
local_node = sp.factor(
    (V**2 - (a / 2 + xi) ** 3 + sp.Rational(3, 4) * a**2 * (a / 2 + xi) + a**3 / 4)
    - a**3 / 2
)
require(
    sp.expand(local_node - (V**2 - sp.Rational(3, 2) * a * xi**2 - xi**3)) == 0,
    "wrong nodal local form",
)
section_xi = sp.factor(thrice_P[0] - a / 2)
require(section_xi == 16 * rho**2 / (9 * a**2), "wrong 3P nodal U order")
require(sp.diff(thrice_P[1], rho).subs(rho, 0) == -1, "wrong 3P nodal V tangent")

semantic_lines = (
    "scope=theta_only_exact_M8_DeltaV_wall",
    generic_line,
    secondary_line,
    "carrier_images=distinct_transpositions_or_collided_disjoint_transpositions",
    "punctured_generators=two_vanishing_loops_plus_carrier_meridians",
    "fixed_sheet_offset=L-n=3",
    "transitivity_requires_total_index>=n-1",
    "section_4134=3P;2P=(-a,-rho);3P_hits_node_at_rho=0",
    "firewall=3P_is_hostile_not_Mordell-Weil_classification",
    "verdict=n16,n15:CONTRADICTION;DeltaV_wall=EMPTY;JC2=OPEN",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4138_HORIZONTAL_CARRIER_AUDIT")
print(generic_line)
print(secondary_line)
print("carrier_collision=two_distinct_transpositions_or_one_(2,2);total_index=2")
print("sharp_positive=actual_ledgers_attain_total_index_n-2_and_exactly_two_orbits")
print("hostile_fixed_offset=n+2:total_index_n-1;TRANSITIVE")
print("hostile_carrier_defect=3:total_index_n-1;TRANSITIVE")
print("section_identity=2P=(-a,-rho);THM4134_section=3P")
print("section_collision=rho=0->(U,V)=(a/2,0);local=q-a^3/2=V^2-(3a/2)xi^2-xi^3")
print("firewall=section_identity_is_hostile_control_not_Mordell-Weil_completeness")
print("semantic_sha256=" + semantic_sha256)
print("verdict=n16:CONTRADICTION;n15:CONTRADICTION;DeltaV_wall=EMPTY;JC2=OPEN")
