#!/usr/bin/env python3
"""Clean-room exact audit of the two 72/539 profiles.

This implementation deliberately does not import the proposed bridge audit or
either canonical producer.  Coefficient sections are reconstructed by cube-edge
intersection.  Circle sets are reconstructed by a separate exact wall-cell
partition with direct midpoint and boundary evaluation.  Every gate remains
active under ``python -O``.
"""

from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
from math import gcd


R = F(3, 14)
EPS = F(1, 14)
TARGET = F(72, 539)
CHECKS = 0


def require(condition: bool, detail: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(detail)


def dot(a: tuple, b: tuple) -> F:
    return sum((F(x) * F(y) for x, y in zip(a, b)), F())


def section_vertices(
    normal: tuple[int, int, int], level: int | F, lo: int | F, hi: int | F
) -> tuple[tuple[F, F, F], ...]:
    """Vertices of a plane/cube section, obtained only from cube edges."""
    level, lo, hi = F(level), F(lo), F(hi)
    found: set[tuple[F, F, F]] = set()

    # Plane-containing cube vertices cover edges parallel to a zero coefficient.
    for corner in product((lo, hi), repeat=3):
        if dot(normal, corner) == level:
            found.add(corner)

    # Every remaining section vertex is the intersection with a cube edge.
    for free in range(3):
        if normal[free] == 0:
            continue
        fixed = [i for i in range(3) if i != free]
        for values in product((lo, hi), repeat=2):
            point = [F()] * 3
            point[fixed[0]], point[fixed[1]] = values
            residue = level - sum(
                (F(normal[i]) * point[i] for i in fixed), F()
            )
            point[free] = residue / normal[free]
            if lo <= point[free] <= hi:
                found.add(tuple(point))

    answer = tuple(sorted(found))
    require(bool(answer), ("empty section", normal, level, lo, hi))
    require(
        all(dot(normal, point) == level for point in answer),
        ("section equation", normal, level),
    )
    require(
        all(lo <= coordinate <= hi for point in answer for coordinate in point),
        ("section cube containment", normal, level),
    )
    return answer


def cross_coordinate(
    normal: tuple[int, int, int],
    speed: tuple[F, F, F],
    error: tuple[F, F, F],
    axis: int | None = None,
) -> F:
    """A projection coordinate; its section width is axis-independent."""
    i = axis if axis is not None else next(index for index, value in enumerate(normal) if value)
    require(normal[i] != 0, ("zero projection axis", normal, i))
    j, k = (i + 1) % 3, (i + 2) % 3
    return (speed[j] * error[k] - speed[k] * error[j]) / normal[i]


def signed_sectors(pattern: tuple[int, int, int]) -> tuple[tuple[int, int, int], ...]:
    """All mixed-sign sectors, modulo global sign, independently generated."""
    sectors: set[tuple[int, int, int]] = set()
    for order in set(permutations(pattern)):
        for signs in product((-1, 1), repeat=3):
            vector = tuple(sign * value for sign, value in zip(signs, order))
            if not (any(x < 0 for x in vector) and any(x > 0 for x in vector)):
                continue
            first = next(x for x in vector if x)
            if first < 0:
                vector = tuple(-x for x in vector)
            sectors.add(vector)
    return tuple(sorted(sectors))


def defect_data(pattern: tuple[int, int, int]) -> tuple[tuple[int, ...], F]:
    all_units = all(value % 3 for value in pattern)
    bound = (3 * sum(pattern) - 1) // 14
    defects = tuple(
        delta
        for delta in range(-bound, bound + 1)
        if (delta % 3 == 0) == all_units
    )
    return defects, (F(2, 3) if all_units else F(1, 3))


def coefficient_ceiling(
    pattern: tuple[int, int, int]
) -> tuple[F, tuple[tuple[int, int, int], tuple[F, F, F]], tuple[int, ...]]:
    defects, weight = defect_data(pattern)
    best = F(-1)
    witnesses: list[tuple[tuple[int, int, int], tuple[F, F, F]]] = []
    for normal in signed_sectors(pattern):
        speed_vertices = section_vertices(normal, 0, 0, 1)
        error_sections = [section_vertices(normal, delta, -R, R) for delta in defects]
        for speed in speed_vertices:
            value = F()
            for section in error_sections:
                coordinates = [cross_coordinate(normal, speed, error) for error in section]
                value += weight * (max(coordinates) - min(coordinates))
            if value > best:
                best = value
                witnesses = [(normal, speed)]
            elif value == best:
                witnesses.append((normal, speed))
    require(best >= 0 and bool(witnesses), ("coefficient ceiling", pattern))
    return best, min(witnesses), defects


def fixed_coefficient_profile() -> tuple[tuple[int, F, F, F], ...]:
    normal = (-11, 4, 7)
    speed = (F(1), F(1), F(1))
    rows = []
    for delta in (-3, 0, 3):
        section = section_vertices(normal, delta, -R, R)
        # Axis one is the report's declared coordinate (e_0-e_2)/4.
        values = [cross_coordinate(normal, speed, error, 1) for error in section]
        rows.append((delta, min(values), max(values), F(2, 3) * (max(values) - min(values))))
    return tuple(rows)


def circle_distance(value: F) -> F:
    value %= 1
    return min(value, 1 - value)


def is_bad(speed: int, phase: F) -> bool:
    return circle_distance(F(speed) * phase) < EPS


def physical_owners(tails: tuple[int, ...], x: F, sheet: int) -> frozenset[int]:
    return frozenset(speed for speed in tails if is_bad(speed, x + F(sheet, 2)))


def quotient_pair_active(pair: tuple[int, int], y: F) -> bool:
    return all(
        any(is_bad(speed, (y + sheet) / 2) for speed in pair)
        for sheet in (0, 1)
    )


def centered_seventh(value: int) -> int:
    """The representative in {-3,...,3}; inputs used here are integral."""
    return (value + 3) % 7 - 3


def seventh_pair_formula(pair: tuple[int, int]) -> tuple[int, int, F]:
    p, q = pair
    require(p % 2 == q % 2 == 1 and gcd(p, q) == 1, ("primitive odd pair", pair))
    alpha, beta = (p + q) // 2, (q - p) // 2
    d, e = centered_seventh(alpha), centered_seventh(beta)
    mass = F(2, 49) * (1 + F(e * e - d * d, p * q))
    return d, e, mass


def walls_for_physical(tails: tuple[int, ...]) -> tuple[F, ...]:
    walls = {F(0), F(1)}
    for speed in tails:
        for sheet in (0, 1):
            for integer in range(-speed, 2 * speed + 1):
                for sign in (-1, 1):
                    x = (F(integer) + sign * EPS) / speed - F(sheet, 2)
                    if 0 <= x <= 1:
                        walls.add(x)
    return tuple(sorted(walls))


def walls_for_pairs(pairs: tuple[tuple[int, int], ...]) -> tuple[F, ...]:
    walls = {F(0), F(1)}
    for pair in pairs:
        for speed in pair:
            for sheet in (0, 1):
                for integer in range(-speed, 2 * speed + 1):
                    for sign in (-1, 1):
                        y = 2 * (F(integer) + sign * EPS) / speed - sheet
                        if 0 <= y <= 1:
                            walls.add(y)
    return tuple(sorted(walls))


def merge_active_cells(
    walls: tuple[F, ...], active, boundary_active
) -> tuple[tuple[F, F], ...]:
    cells = []
    for left, right in zip(walls, walls[1:]):
        if active((left + right) / 2):
            if cells and cells[-1][1] == left and boundary_active(left):
                cells[-1] = (cells[-1][0], right)
            else:
                cells.append((left, right))
    return tuple(cells)


def pair_components(pair: tuple[int, int]) -> tuple[tuple[F, F], ...]:
    walls = walls_for_pairs((pair,))
    return merge_active_cells(
        walls,
        lambda y: quotient_pair_active(pair, y),
        lambda y: quotient_pair_active(pair, y),
    )


def physical_mass(tails: tuple[int, ...]) -> tuple[F, tuple[tuple[F, F], ...]]:
    walls = walls_for_physical(tails)

    def active(x: F) -> bool:
        return bool(physical_owners(tails, x, 0) and physical_owners(tails, x, 1))

    components = merge_active_cells(walls, active, active)
    mass = sum((right - left for left, right in components), F())
    return mass, components


def same_owner_mass(speed: int) -> F:
    walls = walls_for_physical((speed,))

    def active(x: F) -> bool:
        return is_bad(speed, x) and is_bad(speed, x + F(1, 2))

    return sum(
        (
            right - left
            for left, right in zip(walls, walls[1:])
            if active((left + right) / 2)
        ),
        F(),
    )


def affine_maps(source: tuple[int, ...], target: tuple[int, ...]) -> tuple:
    found = []
    for image in permutations(target):
        slope = F(image[1] - image[0], source[1] - source[0])
        intercept = F(image[0]) - slope * source[0]
        if tuple(slope * x + intercept for x in source) == image:
            found.append((slope, intercept, image))
    return tuple(found)


def common_scale_maps(source: tuple[int, ...], target: tuple[int, ...]) -> tuple:
    found = []
    for image in permutations(target):
        scale = F(image[0], source[0])
        if tuple(scale * x for x in source) == image:
            found.append((scale, image))
    return tuple(found)


def eligible(pattern: tuple[int, int, int]) -> bool:
    return (
        sum(value != 0 for value in pattern) >= 2
        and gcd(gcd(*pattern[:2]), pattern[2]) == 1
        and sum(value % 3 == 0 for value in pattern) <= 1
        and pattern != (0, 1, 1)
    )


def main() -> None:
    profile = fixed_coefficient_profile()
    expected_profile = (
        (-3, F(9, 196), F(3, 28), F(2, 49)),
        (0, F(-3, 77), F(3, 77), F(4, 77)),
        (3, F(-3, 28), F(-9, 196), F(2, 49)),
    )
    require(profile == expected_profile, ("coefficient profile", profile))
    coefficient_value, witness, defects = coefficient_ceiling((4, 7, 11))
    require(defects == (-3, 0, 3), ("defects", defects))
    require(coefficient_value == TARGET, ("coefficient value", coefficient_value))
    require(
        witness[1] == (F(1), F(1), F(1))
        and sum(witness[0]) == 0
        and tuple(sorted(abs(value) for value in witness[0])) == (4, 7, 11),
        ("equal-speed boundary witness", witness),
    )

    pairs = ((1, 7), (1, 11), (7, 11))
    components = {pair: pair_components(pair) for pair in pairs}
    expected_components = {
        (1, 7): ((F(6, 49), F(1, 7)), (F(6, 7), F(43, 49))),
        (1, 11): ((F(6, 77), F(8, 77)), (F(69, 77), F(71, 77))),
        (7, 11): ((F(13, 49), F(2, 7)), (F(5, 7), F(36, 49))),
    }
    require(components == expected_components, ("pair components", components))
    pair_masses = tuple(
        sum((right - left for left, right in components[pair]), F()) for pair in pairs
    )
    require(pair_masses == (F(2, 49), F(4, 77), F(2, 49)), pair_masses)
    require(
        tuple(seventh_pair_formula(pair) for pair in pairs)
        == ((-3, 3, F(2, 49)), (-1, -2, F(4, 77)), (2, 2, F(2, 49))),
        "seventh-rounding coordinates",
    )

    all_cells = [cell for pair in pairs for cell in components[pair]]
    for first, second in zip(sorted(all_cells), sorted(all_cells)[1:]):
        require(first[1] <= second[0], ("pair overlap", first, second))

    quotient_union_mass = sum(pair_masses, F())
    dyadic_mass, dyadic_components = physical_mass((1, 7, 11))
    require(quotient_union_mass == dyadic_mass == TARGET, (quotient_union_mass, dyadic_mass))
    require(len(dyadic_components) == 12, ("physical components", dyadic_components))

    # Parity/owner hostile: only the reinterpreted even coefficient can own both lifts.
    require(tuple(same_owner_mass(s) for s in (1, 7, 11)) == (F(), F(), F()), "odd owner")
    require(same_owner_mass(4) == F(1, 7), ("even owner", same_owner_mass(4)))
    hostile_mass, _ = physical_mass((4, 7, 11))
    require(hostile_mass == F(9, 49) and hostile_mass != TARGET, hostile_mass)

    require(affine_maps((1, 7, 11), (4, 7, 11)) == (), "affine maps")
    require(common_scale_maps((1, 7, 11), (4, 7, 11)) == (), "common scales")

    nearby_coefficient, _, _ = coefficient_ceiling((6, 7, 13))
    nearby_physical, _ = physical_mass((1, 7, 13))
    require(nearby_coefficient == F(233, 1911), nearby_coefficient)
    require(nearby_physical == F(80, 637) == F(240, 1911), nearby_physical)
    require(nearby_coefficient != nearby_physical, "nearby family hostile")

    hits = []
    universe = 0
    for pattern in combinations_with_replacement(range(19), 3):
        if not eligible(pattern):
            continue
        universe += 1
        value, _, _ = coefficient_ceiling(pattern)
        if value == TARGET:
            hits.append(pattern)
    require(universe == 750 and hits == [(4, 7, 11)], (universe, hits))

    print("CLEANROOM_72_539_BRIDGE_AUDIT")
    print("coefficient_profile=-3:[9/196,3/28]->2/49;0:[-3/77,3/77]->4/77;+3:[-3/28,-9/196]->2/49")
    print("coefficient_ceiling=(4,7,11):72/539;defects=-3,0,3;box_hits=1/750")
    print("dyadic_profile=(1,7):2/49;(1,11):4/77;(7,11):2/49;pair_components=6")
    print("physical_(1,7,11)=72/539;physical_components=12;owner_overlap_correction=0")
    print("typed_profile_map=central_atom_unique;outer_atoms_unordered;addresses_and_owners_lost")
    print("affine_maps=0;common_scale_maps=0")
    print("owner_hostile=same_owner_mass_odd:0;speed4:1/7;physical_(4,7,11)=9/49")
    print("nearby_hostile=coefficient_(6,7,13):233/1911;physical_(1,7,13):240/1911")
    print(f"PASS checks={CHECKS}")


if __name__ == "__main__":
    main()
