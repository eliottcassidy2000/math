#!/usr/bin/env python3
"""Exact monodromy, block-lattice, and decomposition controls for THM-2821.

The calculation is deliberately elementary.  SymPy is used only for exact
polynomial identities.  Permutation groups, subgroup intervals, and uniform
block systems are reconstructed from finite tuples with no group-theory CAS.
"""

import ast
from itertools import combinations
from math import gcd
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    return any(
        isinstance(node, ast.Assert)
        for node in ast.walk(
            ast.parse(Path(path).read_text(encoding="utf-8"))
        )
    )


def identity(degree):
    return tuple(range(degree))


def compose(left, right):
    """Return left after right."""
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    answer = [0] * len(permutation)
    for source, target in enumerate(permutation):
        answer[target] = source
    return tuple(answer)


def permutation_from_cycles(degree, cycle_list):
    answer = list(range(degree))
    for cycle in cycle_list:
        zero_based = tuple(entry - 1 for entry in cycle)
        for source, target in zip(
            zero_based, zero_based[1:] + zero_based[:1]
        ):
            answer[source] = target
    return tuple(answer)


def nontrivial_cycles(permutation):
    seen = set()
    answer = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        orbit = []
        point = start
        while point not in seen:
            seen.add(point)
            orbit.append(point + 1)
            point = permutation[point]
        if len(orbit) > 1:
            answer.append(tuple(orbit))
    return tuple(answer)


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        length = 0
        point = start
        while point not in seen:
            seen.add(point)
            length += 1
            point = permutation[point]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def lcm(left, right):
    return left * right // gcd(left, right)


def permutation_order(permutation):
    order = 1
    for length in cycle_type(permutation):
        order = lcm(order, length)
    return order


def generated_group(generators):
    degree = len(generators[0])
    steps = tuple(generators) + tuple(inverse(item) for item in generators)
    unit = identity(degree)
    group = {unit}
    frontier = [unit]
    while frontier:
        current = frontier.pop()
        for step in steps:
            product = compose(current, step)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return frozenset(group)


def conjugate(left, right):
    return compose(compose(left, right), inverse(left))


def is_transitive(group):
    return {element[0] for element in group} == set(range(len(next(iter(group)))))


def all_subgroups(group):
    """Enumerate all subgroups by adjoining every possible group element."""
    degree = len(next(iter(group)))
    unit_subgroup = frozenset({identity(degree)})
    known = {unit_subgroup}
    frontier = [unit_subgroup]
    while frontier:
        subgroup = frontier.pop()
        for element in group.difference(subgroup):
            generated = generated_group(tuple(subgroup) + (element,))
            require(generated.issubset(group), "subgroup escaped ambient group")
            if generated not in known:
                known.add(generated)
                frontier.append(generated)
    return frozenset(known)


def uniform_partitions(points, block_size):
    points = tuple(sorted(points))
    if not points:
        yield ()
        return
    first = points[0]
    for companions in combinations(points[1:], block_size - 1):
        block = frozenset((first,) + companions)
        remainder = tuple(point for point in points if point not in block)
        for tail in uniform_partitions(remainder, block_size):
            yield tuple(sorted((block,) + tail, key=lambda item: min(item)))


def noncrossing_matchings(points):
    """Generate every noncrossing perfect matching in cyclic order."""
    points = tuple(points)
    if not points:
        yield ()
        return
    first = points[0]
    for partner_index in range(1, len(points), 2):
        partner = points[partner_index]
        inside = points[1:partner_index]
        outside = points[partner_index + 1 :]
        for left in noncrossing_matchings(inside):
            for right in noncrossing_matchings(outside):
                yield ((first, partner),) + left + right


def matching_permutation(matching, degree):
    permutation = list(range(degree))
    for left, right in matching:
        permutation[left] = right
        permutation[right] = left
    return tuple(permutation)


def rotate_matching(matching, shift, degree):
    return tuple(
        sorted(
            tuple(sorted(((left + shift) % degree, (right + shift) % degree)))
            for left, right in matching
        )
    )


def image_of_block(permutation, block):
    return frozenset(permutation[point] for point in block)


def invariant_uniform_partitions(group):
    degree = len(next(iter(group)))
    answer = {}
    for block_size in range(2, degree):
        if degree % block_size:
            continue
        systems = []
        for partition in uniform_partitions(range(degree), block_size):
            block_set = frozenset(partition)
            if all(
                image_of_block(element, block) in block_set
                for element in group
                for block in partition
            ):
                systems.append(partition)
        answer[block_size] = tuple(systems)
    return answer


def orbit_block_system(group, subgroup, basepoint=0):
    block = frozenset(element[basepoint] for element in subgroup)
    translates = frozenset(
        image_of_block(element, block) for element in group
    )
    require(
        set().union(*translates) == set(range(len(next(iter(group))))),
        "subgroup translates do not cover the fibre",
    )
    require(
        sum(len(item) for item in translates)
        == len(next(iter(group))),
        "subgroup translates are not disjoint",
    )
    return tuple(sorted(translates, key=lambda item: min(item)))


def one_based_partition(partition):
    return tuple(
        tuple(sorted(point + 1 for point in block)) for block in partition
    )


def subgroup_interval(group):
    stabilizer = frozenset(element for element in group if element[0] == 0)
    subgroups = all_subgroups(group)
    interval = tuple(
        sorted(
            (
                subgroup
                for subgroup in subgroups
                if stabilizer.issubset(subgroup)
            ),
            key=lambda subgroup: (
                len(subgroup),
                sorted(subgroup),
            ),
        )
    )
    return stabilizer, interval


def check_symbolic_carriers():
    x, y, u, z = sp.symbols("x y u z")

    power_cubic = 2 * x**3 - 1
    power_response = sp.cancel(
        power_cubic**2 / (power_cubic**2 - 1)
    )
    power_D = x**3 * (x**3 - 1)
    power_E = x**3 - sp.Rational(1, 2)
    require(
        sp.expand(power_E**2 - power_D) == sp.Rational(1, 4),
        "power accessory identity",
    )
    require(
        sp.cancel(power_E**2 / power_D - power_response) == 0,
        "power response reconstruction",
    )
    require(
        sp.expand(
            power_cubic**2 - 1 - 4 * x**3 * (x**3 - 1)
        )
        == 0,
        "power infinity fibre",
    )
    require(
        sp.cancel(
            power_response
            - (2 * u - 1) ** 2 / ((2 * u - 1) ** 2 - 1)
        ).subs(u, x**3)
        == 0,
        "power degree-three-then-two decomposition",
    )
    require(
        sp.cancel(
            sp.together(
                1 / (1 - power_response) - (1 - power_cubic**2)
            )
        )
        == 0,
        "power target polynomialization",
    )

    cheb3 = 4 * y**3 - 3 * y
    cheb2 = 2 * z**2 - 1
    cheb6 = 32 * y**6 - 48 * y**4 + 18 * y**2 - 1
    cheb_response = sp.cancel(cheb3**2 / (cheb3**2 - 1))
    cheb_D = (y**2 - sp.Rational(1, 4)) ** 2 * (y**2 - 1)
    cheb_E = cheb3 / 4
    require(
        sp.expand(cheb_E**2 - cheb_D) == sp.Rational(1, 16),
        "Chebyshev accessory identity",
    )
    require(
        sp.cancel(cheb_E**2 / cheb_D - cheb_response) == 0,
        "Chebyshev response reconstruction",
    )
    require(
        sp.expand(
            cheb3**2
            - 1
            - (y - 1)
            * (y + 1)
            * (2 * y - 1) ** 2
            * (2 * y + 1) ** 2
        )
        == 0,
        "Chebyshev infinity fibre",
    )
    require(
        sp.expand(cheb6 - (2 * cheb3**2 - 1)) == 0,
        "T6=T2 after T3",
    )
    require(
        sp.expand(
            cheb6
            - (4 * cheb2**3 - 3 * cheb2).subs(z, y)
        )
        == 0,
        "T6=T3 after T2",
    )
    require(
        sp.cancel(cheb_response - (cheb6 + 1) / (cheb6 - 1)) == 0,
        "Chebyshev T6 target identity",
    )
    reverse_outer = sp.cancel(
        u * (4 * u - 3) ** 2
        / (u * (4 * u - 3) ** 2 - 1)
    )
    require(
        sp.cancel(cheb_response - reverse_outer.subs(u, y**2)) == 0,
        "Chebyshev degree-two-then-three decomposition",
    )
    require(
        sp.expand(
            u * (4 * u - 3) ** 2
            - 1
            - (u - 1) * (4 * u - 1) ** 2
        )
        == 0,
        "reverse cubic outer denominator",
    )
    require(
        sp.cancel(
            sp.together(
                1 / (1 - cheb_response) - (1 - cheb3**2)
            )
        )
        == 0,
        "Chebyshev target polynomialization",
    )
    return {
        "power_infinity": sp.factor(power_cubic**2 - 1),
        "cheb_infinity": sp.factor(cheb3**2 - 1),
        "power_polynomial": sp.factor(1 - power_cubic**2),
        "cheb_polynomial": sp.factor(1 - cheb3**2),
        "reverse_outer": sp.factor(reverse_outer),
    }


def check_power_group():
    degree = 6
    sigma_zero = permutation_from_cycles(
        degree, ((1, 4), (2, 5), (3, 6))
    )
    sigma_infinity = permutation_from_cycles(degree, ((1, 2, 3),))
    sigma_one = inverse(compose(sigma_zero, sigma_infinity))
    require(
        compose(compose(sigma_zero, sigma_infinity), sigma_one)
        == identity(degree),
        "power branch product",
    )
    require(cycle_type(sigma_zero) == (2, 2, 2), "power zero inertia")
    require(
        cycle_type(sigma_infinity) == (3, 1, 1, 1),
        "power infinity inertia",
    )
    require(cycle_type(sigma_one) == (6,), "power one inertia")

    group = generated_group((sigma_zero, sigma_infinity))
    require(len(group) == 18, "power monodromy order")
    require(is_transitive(group), "power monodromy transitivity")
    second_rotation = conjugate(sigma_zero, sigma_infinity)
    base = generated_group((sigma_infinity, second_rotation))
    require(len(base) == 9, "power C3-square base")
    require(
        all(permutation_order(element) in (1, 3) for element in base),
        "power base exponent",
    )
    require(
        compose(sigma_infinity, second_rotation)
        == compose(second_rotation, sigma_infinity),
        "power base commutativity",
    )
    require(
        conjugate(sigma_zero, sigma_infinity) == second_rotation
        and conjugate(sigma_zero, second_rotation) == sigma_infinity,
        "power factor swap",
    )
    require(
        all(
            conjugate(element, base_element) in base
            for element in group
            for base_element in base
        ),
        "power base normality",
    )

    systems = invariant_uniform_partitions(group)
    expected = (
        (
            frozenset((0, 1, 2)),
            frozenset((3, 4, 5)),
        ),
    )
    require(systems[2] == (), "power acquired a size-two block system")
    require(systems[3] == expected, "power size-three block system changed")

    stabilizer, interval = subgroup_interval(group)
    require(len(stabilizer) == 3, "power point stabilizer")
    require(
        stabilizer == generated_group((second_rotation,)),
        "power stabilizer is not the second C3 factor",
    )
    require(
        tuple(len(subgroup) for subgroup in interval) == (3, 9, 18),
        "power intermediate subgroup interval",
    )
    require(
        interval[1] == base and interval[2] == group,
        "power named subgroup interval",
    )
    interval_systems = tuple(
        orbit_block_system(group, subgroup) for subgroup in interval
    )
    require(
        tuple(len(system[0]) for system in interval_systems) == (1, 3, 6),
        "power subgroup/block correspondence",
    )
    return {
        "tuple": (sigma_zero, sigma_infinity, sigma_one),
        "group": group,
        "systems": systems,
        "stabilizer": stabilizer,
        "interval": interval,
        "interval_systems": interval_systems,
    }


def check_chebyshev_group():
    degree = 6
    sigma_zero = permutation_from_cycles(
        degree, ((1, 6), (2, 5), (3, 4))
    )
    sigma_infinity = permutation_from_cycles(
        degree, ((2, 6), (3, 5))
    )
    sigma_one = inverse(compose(sigma_zero, sigma_infinity))
    require(
        compose(compose(sigma_zero, sigma_infinity), sigma_one)
        == identity(degree),
        "Chebyshev branch product",
    )
    require(
        cycle_type(sigma_zero) == (2, 2, 2),
        "Chebyshev zero inertia",
    )
    require(
        cycle_type(sigma_infinity) == (2, 2, 1, 1),
        "Chebyshev infinity inertia",
    )
    require(cycle_type(sigma_one) == (6,), "Chebyshev one inertia")

    rotation = sigma_one
    reflection = sigma_infinity
    group = generated_group((rotation, reflection))
    require(len(group) == 12, "Chebyshev monodromy order")
    require(is_transitive(group), "Chebyshev monodromy transitivity")
    require(permutation_order(rotation) == 6, "Chebyshev rotation order")
    require(permutation_order(reflection) == 2, "Chebyshev reflection order")
    require(
        conjugate(reflection, rotation) == inverse(rotation),
        "Chebyshev dihedral relation",
    )

    systems = invariant_uniform_partitions(group)
    expected_size_two = (
        (
            frozenset((0, 3)),
            frozenset((1, 4)),
            frozenset((2, 5)),
        ),
    )
    expected_size_three = (
        (
            frozenset((0, 2, 4)),
            frozenset((1, 3, 5)),
        ),
    )
    require(
        systems[2] == expected_size_two,
        "Chebyshev size-two block system changed",
    )
    require(
        systems[3] == expected_size_three,
        "Chebyshev size-three block system changed",
    )

    stabilizer, interval = subgroup_interval(group)
    require(len(stabilizer) == 2, "Chebyshev point stabilizer")
    require(
        stabilizer == generated_group((reflection,)),
        "Chebyshev stabilizer is not the vertex reflection",
    )
    require(
        tuple(len(subgroup) for subgroup in interval) == (2, 4, 6, 12),
        "Chebyshev intermediate subgroup interval",
    )
    require(
        interval[1]
        == generated_group(
            (reflection, compose(compose(rotation, rotation), rotation))
        ),
        "Chebyshev order-four middle subgroup",
    )
    require(
        interval[2]
        == generated_group((reflection, compose(rotation, rotation))),
        "Chebyshev order-six middle subgroup",
    )
    require(interval[3] == group, "Chebyshev full interval endpoint")
    interval_systems = tuple(
        orbit_block_system(group, subgroup) for subgroup in interval
    )
    require(
        tuple(len(system[0]) for system in interval_systems)
        == (1, 2, 3, 6),
        "Chebyshev subgroup/block correspondence",
    )
    return {
        "tuple": (sigma_zero, sigma_infinity, sigma_one),
        "group": group,
        "systems": systems,
        "stabilizer": stabilizer,
        "interval": interval,
        "interval_systems": interval_systems,
    }


def check_nielsen_convention():
    """Cross-check THM-2817's rho/tau convention on all five matchings.

    THM-2817 takes rho to be the inverse of inertia over 1 and tau to be
    inertia over 0.  Its pole inertia is tau*rho, since
    tau*(tau*rho)*rho^{-1}=1.
    """
    degree = 6
    rho = tuple((index + 1) % degree for index in range(degree))
    matchings = tuple(noncrossing_matchings(range(degree)))
    require(len(matchings) == 5, "hexagon Catalan matching count")
    orbit_rows = {}
    for matching in matchings:
        tau = matching_permutation(matching, degree)
        pole = compose(tau, rho)
        one = inverse(rho)
        require(
            compose(compose(tau, pole), one) == identity(degree),
            "Nielsen branch-product convention",
        )
        group = generated_group((rho, tau))
        systems = invariant_uniform_partitions(group)
        block_sizes = tuple(
            block_size
            for block_size, block_systems in sorted(systems.items())
            if block_systems
        )
        orbit = frozenset(
            rotate_matching(matching, shift, degree)
            for shift in range(degree)
        )
        orbit_rows.setdefault(orbit, set()).add(
            (cycle_type(pole), len(group), block_sizes)
        )
    require(len(orbit_rows) == 2, "hexagon matching rotation-orbit count")
    require(
        all(len(rows) == 1 for rows in orbit_rows.values()),
        "signature changed within a matching rotation orbit",
    )
    signatures = {
        next(iter(rows)) for rows in orbit_rows.values()
    }
    require(
        signatures
        == {
            ((3, 1, 1, 1), 18, (3,)),
            ((2, 2, 1, 1), 12, (2, 3)),
        },
        "Nielsen carrier signatures",
    )
    return signatures


def tuple_text(branch_tuple):
    names = ("sigma_0", "sigma_infinity", "sigma_1")
    return ", ".join(
        f"{name}={nontrivial_cycles(permutation)}"
        for name, permutation in zip(names, branch_tuple)
    )


def main():
    require(not has_asserts(Path(__file__)), "truth-bearing Python assert found")
    symbolic = check_symbolic_carriers()
    power = check_power_group()
    chebyshev = check_chebyshev_group()
    nielsen_signatures = check_nielsen_convention()

    print("THM-2821 SEXTIC MONODROMY AND BLOCK-LATTICE AUDIT")
    print(f"power branch tuple: {tuple_text(power['tuple'])}")
    print("power cycle types: (2^3), (3,1^3), (6)")
    print("power group: (C3 x C3) semidirect C2 (factor swap), order=18")
    print(
        "power nontrivial block systems:",
        one_based_partition(power["systems"][3][0]),
    )
    print(
        "power point-stabilizer interval orders:",
        tuple(len(subgroup) for subgroup in power["interval"]),
    )
    print("power nontrivial decomposition degrees: 3 then 2 only")
    print(f"power infinity fibre: {symbolic['power_infinity']}")
    print(f"power polynomialized cover: {symbolic['power_polynomial']}")

    print(f"Chebyshev branch tuple: {tuple_text(chebyshev['tuple'])}")
    print("Chebyshev cycle types: (2^3), (2^2,1^2), (6)")
    print("Chebyshev group: Dih(C6), order=12")
    print(
        "Chebyshev size-two block system:",
        one_based_partition(chebyshev["systems"][2][0]),
    )
    print(
        "Chebyshev size-three block system:",
        one_based_partition(chebyshev["systems"][3][0]),
    )
    print(
        "Chebyshev point-stabilizer interval orders:",
        tuple(len(subgroup) for subgroup in chebyshev["interval"]),
    )
    print("Chebyshev nontrivial decomposition degrees: 3 then 2; 2 then 3")
    print(f"Chebyshev infinity fibre: {symbolic['cheb_infinity']}")
    print(f"Chebyshev polynomialized cover: {symbolic['cheb_polynomial']}")
    print(f"reverse y^2 outer map: {symbolic['reverse_outer']}")
    print(
        "five noncrossing matchings / two rotation orbits:",
        sorted(nielsen_signatures),
    )
    print("all uniform block systems and intermediate subgroups: PASS")
    print("accessory identities and both functional decompositions: PASS")
    print("assert_nodes=0")


if __name__ == "__main__":
    main()
