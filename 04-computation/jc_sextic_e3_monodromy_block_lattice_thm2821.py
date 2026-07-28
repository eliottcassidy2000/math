#!/usr/bin/env python3
"""Exact referee for THM-2821's two sextic monodromy block lattices."""

from itertools import combinations

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(left, right):
    """Permutation product left after right."""
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    result = [0] * len(permutation)
    for source, target in enumerate(permutation):
        result[target] = source
    return tuple(result)


def power(permutation, exponent):
    result = tuple(range(len(permutation)))
    for _ in range(exponent):
        result = compose(permutation, result)
    return result


def order(permutation):
    result = tuple(range(len(permutation)))
    for exponent in range(1, 10_000):
        result = compose(permutation, result)
        if result == tuple(range(len(permutation))):
            return exponent
    raise RuntimeError("permutation order search escaped")


def cycles(permutation, include_fixed=True):
    seen = set()
    result = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle = []
        point = start
        while point not in seen:
            seen.add(point)
            cycle.append(point)
            point = permutation[point]
        if include_fixed or len(cycle) > 1:
            result.append(tuple(cycle))
    return tuple(result)


def cycle_type(permutation):
    return tuple(sorted((len(cycle) for cycle in cycles(permutation)), reverse=True))


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = compose(current, generator)
            if candidate not in group:
                group.add(candidate)
                frontier.append(candidate)
    return frozenset(group)


def subgroup(generators):
    return generated_group(generators)


def blocks_containing_zero(group):
    """All nontrivial blocks containing the distinguished sheet zero."""
    degree = len(next(iter(group)))
    blocks = []
    for size in range(2, degree):
        for tail in combinations(range(1, degree), size - 1):
            block = frozenset((0,) + tail)
            good = True
            for permutation in group:
                image = frozenset(permutation[index] for index in block)
                if image != block and image & block:
                    good = False
                    break
            if good:
                blocks.append(tuple(sorted(block)))
    return tuple(blocks)


def noncrossing_matchings(points):
    points = tuple(points)
    if not points:
        return ((),)
    first = points[0]
    result = []
    for partner_index in range(1, len(points), 2):
        partner = points[partner_index]
        inside = points[1:partner_index]
        outside = points[partner_index + 1 :]
        for left in noncrossing_matchings(inside):
            for right in noncrossing_matchings(outside):
                result.append(
                    tuple(sorted(((first, partner),) + left + right))
                )
    return tuple(result)


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


def rotation_orbit(matching, degree):
    return frozenset(rotate_matching(matching, shift, degree) for shift in range(degree))


def main():
    degree = 6
    identity = tuple(range(degree))
    rho = (1, 2, 3, 4, 5, 0)
    tau_power = (1, 0, 3, 2, 5, 4)
    tau_chebyshev = (1, 0, 5, 4, 3, 2)

    require(order(rho) == 6, "full-cycle order changed")
    require(cycle_type(compose(tau_power, rho)) == (3, 1, 1, 1),
            "power pole passport changed")
    require(cycle_type(compose(tau_chebyshev, rho)) == (2, 2, 1, 1),
            "Chebyshev pole passport changed")

    group_power = generated_group((rho, tau_power))
    group_chebyshev = generated_group((rho, tau_chebyshev))
    require(len(group_power) == 18, "power monodromy order changed")
    require(len(group_chebyshev) == 12, "Chebyshev monodromy order changed")

    even_cycle = compose(rho, tau_power)
    odd_cycle = compose(tau_power, rho)
    require(cycles(even_cycle, include_fixed=False) == ((0, 2, 4),),
            "power even three-cycle changed")
    require(cycles(odd_cycle, include_fixed=False) == ((1, 3, 5),),
            "power odd three-cycle changed")
    require(compose(even_cycle, odd_cycle) == compose(odd_cycle, even_cycle),
            "power base cycles stopped commuting")
    base = subgroup((even_cycle, odd_cycle))
    require(len(base) == 9, "power C3-square kernel changed")
    require(all(order(element) in (1, 3) for element in base),
            "power kernel is not elementary C3-square")
    require(compose(tau_power, compose(even_cycle, tau_power)) == odd_cycle,
            "power involution stopped swapping the C3 factors")
    require(group_power == frozenset(
        set(base) | {compose(element, tau_power) for element in base}
    ), "power semidirect-product decomposition changed")

    require(compose(tau_chebyshev, compose(rho, tau_chebyshev)) == inverse(rho),
            "Chebyshev dihedral relation changed")
    dihedral_normal_form = {
        power(rho, exponent)
        for exponent in range(6)
    } | {
        compose(power(rho, exponent), tau_chebyshev)
        for exponent in range(6)
    }
    require(group_chebyshev == frozenset(dihedral_normal_form),
            "Chebyshev dihedral normal form changed")

    power_blocks = blocks_containing_zero(group_power)
    chebyshev_blocks = blocks_containing_zero(group_chebyshev)
    require(power_blocks == ((0, 2, 4),),
            f"power block lattice changed: {power_blocks}")
    require(chebyshev_blocks == ((0, 3), (0, 2, 4)),
            f"Chebyshev block lattice changed: {chebyshev_blocks}")

    matchings = noncrossing_matchings(range(degree))
    require(len(matchings) == 5, "Catalan matching count changed")
    orbit_rows = {}
    for matching in matchings:
        tau = matching_permutation(matching, degree)
        passport = cycle_type(compose(tau, rho))
        group_order = len(generated_group((rho, tau)))
        blocks = blocks_containing_zero(generated_group((rho, tau)))
        key = tuple(sorted(rotation_orbit(matching, degree)))
        orbit_rows.setdefault(key, set()).add((passport, group_order, blocks))
    require(len(orbit_rows) == 2, "unmarked rotation-orbit count changed")
    signatures = {
        next(iter(rows))
        for rows in orbit_rows.values()
    }
    require(signatures == {
        ((3, 1, 1, 1), 18, ((0, 2, 4),)),
        ((2, 2, 1, 1), 12, ((0, 3), (0, 2, 4))),
    }, f"unmarked signatures changed: {signatures}")

    x, y, z, u = sp.symbols("x y z u")
    outer = z**2 / (z**2 - 1)
    cubic_power = 2 * x**3 - 1
    cubic_chebyshev = 4 * y**3 - 3 * y
    response_power = sp.cancel(outer.subs(z, cubic_power))
    response_chebyshev = sp.cancel(outer.subs(z, cubic_chebyshev))
    chebyshev_outer_three = sp.cancel(
        u * (4 * u - 3) ** 2 / (u * (4 * u - 3) ** 2 - 1)
    )
    require(sp.cancel(response_chebyshev - chebyshev_outer_three.subs(u, y**2)) == 0,
            "Chebyshev degree-two-first decomposition failed")
    require(sp.degree(sp.fraction(response_power)[0], x) == 6,
            "power response degree changed")
    require(sp.degree(sp.fraction(response_chebyshev)[0], y) == 6,
            "Chebyshev response degree changed")
    require(sp.degree(cubic_power, x) == 3 and sp.degree(cubic_chebyshev, y) == 3,
            "inner cubic degree changed")
    require(sp.degree(sp.fraction(outer)[0], z) == 2,
            "common outer degree changed")
    require(sp.degree(sp.fraction(chebyshev_outer_three)[0], u) == 3,
            "Chebyshev alternate outer degree changed")
    require(identity in group_power and identity in group_chebyshev,
            "identity missing from a monodromy group")

    print("THM-2821 SEXTIC E3 MONODROMY / BLOCK-LATTICE AUDIT")
    print("power_passport=(3,1,1,1); group_order=18; "
          "group=(C3xC3)semidirectC2")
    print(f"power_blocks_containing_0={power_blocks}")
    print("power_decomposition=degree3_then_degree2; no_degree2_inner_block")
    print("chebyshev_passport=(2,2,1,1); group_order=12; "
          "group=C6_semidirect_C2_dihedral")
    print(f"chebyshev_blocks_containing_0={chebyshev_blocks}")
    print("chebyshev_decompositions=degree3_then2_and_degree2_then3")
    print("noncrossing_matchings=5; unmarked_rotation_orbits=2")
    print("scope=response block lattice only; no Keller-map decomposition, "
          "JC2, or DC2")


if __name__ == "__main__":
    main()
