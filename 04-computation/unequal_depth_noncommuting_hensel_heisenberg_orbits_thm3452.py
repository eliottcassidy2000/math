#!/usr/bin/env python3
"""Exact controls for provisional THM-3452's weighted Heisenberg law."""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import product
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "c0c00380198dca87af2940d7986199885633b9a042aeae25b557656abadb6d4b"
ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / "01-canon/theorems/THM-3442-smooth-hensel-fibre-vector-field-orbit-law.md":
        "5ca39bf0693f7d314fd2907d113eeb8aaca613a5f6e02171cd5a1a4700b3a11a",
    ROOT / "01-canon/theorems/THM-3446-weighted-depth-commuting-hensel-lattice-action.md":
        "b64e56f02accabd0efd63b16172bb96316184ebeaf3738f5302e71083a62f658",
    ROOT / "01-canon/theorems/THM-3449-noncommuting-smooth-hensel-heisenberg-orbit-law.md":
        "fe7e35c96c500b57825d67ee89609af52fd6584e10b06b8a3a3ae004dbc2762c",
}


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def mixed_mul(left, right, prime: int, first_depth: int,
              second_depth: int, central_depth: int):
    """Multiply x^i y^j z^r with z=xyx^-1y^-1."""
    first_modulus = prime ** first_depth
    second_modulus = prime ** second_depth
    central_modulus = prime ** central_depth
    i, j, r = left
    ii, jj, rr = right
    return (
        (i + ii) % first_modulus,
        (j + jj) % second_modulus,
        (r + rr - j * ii) % central_modulus,
    )


def mixed_inverse(element, prime: int, first_depth: int,
                  second_depth: int, central_depth: int):
    first_modulus = prime ** first_depth
    second_modulus = prime ** second_depth
    central_modulus = prime ** central_depth
    i, j, r = element
    return (
        (-i) % first_modulus,
        (-j) % second_modulus,
        (-r - i * j) % central_modulus,
    )


def mixed_power(element, exponent: int, prime: int, first_depth: int,
                second_depth: int, central_depth: int):
    answer = (0, 0, 0)
    base = element
    while exponent:
        if exponent & 1:
            answer = mixed_mul(
                answer, base, prime, first_depth, second_depth, central_depth
            )
        base = mixed_mul(
            base, base, prime, first_depth, second_depth, central_depth
        )
        exponent //= 2
    return answer


def presentation_controls():
    rows = []
    for prime, first_depth, second_depth, central_depth in (
        (2, 2, 1, 1),
        (3, 2, 1, 1),
        (2, 3, 2, 1),
    ):
        first_modulus = prime ** first_depth
        second_modulus = prime ** second_depth
        central_modulus = prime ** central_depth
        elements = tuple(product(
            range(first_modulus),
            range(second_modulus),
            range(central_modulus),
        ))
        identity = (0, 0, 0)
        x_generator = (1, 0, 0)
        y_generator = (0, 1, 0)
        z_generator = (0, 0, 1)
        require(
            len(elements)
            == prime ** (first_depth + second_depth + central_depth),
            ("mixed presentation order", prime, first_depth, second_depth),
        )
        for element in elements:
            inverse = mixed_inverse(
                element, prime, first_depth, second_depth, central_depth
            )
            require(
                mixed_mul(
                    element, inverse, prime, first_depth, second_depth,
                    central_depth,
                ) == identity,
                ("right inverse", element),
            )
            require(
                mixed_mul(
                    inverse, element, prime, first_depth, second_depth,
                    central_depth,
                ) == identity,
                ("left inverse", element),
            )

        x_inverse = mixed_inverse(
            x_generator, prime, first_depth, second_depth, central_depth
        )
        y_inverse = mixed_inverse(
            y_generator, prime, first_depth, second_depth, central_depth
        )
        commutator = mixed_mul(
            mixed_mul(
                mixed_mul(
                    x_generator, y_generator, prime, first_depth,
                    second_depth, central_depth,
                ),
                x_inverse, prime, first_depth, second_depth, central_depth,
            ),
            y_inverse, prime, first_depth, second_depth, central_depth,
        )
        require(commutator == z_generator, ("commutator", commutator))
        require(
            mixed_power(
                x_generator, first_modulus, prime, first_depth,
                second_depth, central_depth,
            ) == identity,
            "first horizontal order bound",
        )
        require(
            mixed_power(
                x_generator, first_modulus // prime, prime, first_depth,
                second_depth, central_depth,
            ) != identity,
            "first horizontal exact order",
        )
        require(
            mixed_power(
                y_generator, second_modulus, prime, first_depth,
                second_depth, central_depth,
            ) == identity,
            "second horizontal order bound",
        )
        require(
            mixed_power(
                y_generator, second_modulus // prime, prime, first_depth,
                second_depth, central_depth,
            ) != identity,
            "second horizontal exact order",
        )
        require(
            mixed_power(
                z_generator, central_modulus, prime, first_depth,
                second_depth, central_depth,
            ) == identity,
            "central order bound",
        )

        horizontal = tuple(product(
            range(first_modulus), range(second_modulus)
        ))
        for first, second, third in product(horizontal, repeat=3):
            first_second = (
                (first[0] + second[0]) % first_modulus,
                (first[1] + second[1]) % second_modulus,
            )
            second_third = (
                (second[0] + third[0]) % first_modulus,
                (second[1] + third[1]) % second_modulus,
            )
            left_cocycle = (
                -first[1] * second[0]
                - first_second[1] * third[0]
            ) % central_modulus
            right_cocycle = (
                -second[1] * third[0]
                - first[1] * second_third[0]
            ) % central_modulus
            require(
                left_cocycle == right_cocycle,
                ("mixed cocycle", first, second, third),
            )

        center_size = sum(
            1
            for i, j, _ in elements
            if i % central_modulus == 0 and j % central_modulus == 0
        )
        expected_center = prime ** (
            first_depth + second_depth - central_depth
        )
        require(center_size == expected_center, ("center", center_size))
        rows.append((
            prime,
            first_depth,
            second_depth,
            central_depth,
            len(elements),
            center_size,
            commutator,
        ))
    return tuple(rows)


def fibre_points(prime: int, dimension: int, carry_depth: int, level: int,
                 base):
    modulus = prime ** level
    scale = prime ** carry_depth
    width = prime ** (level - carry_depth)
    return tuple(
        tuple((base[index] + scale * offset[index]) % modulus
              for index in range(dimension))
        for offset in product(range(width), repeat=dimension)
    )


def generated_orbits(points, generators):
    unseen = set(points)
    orbits = []
    while unseen:
        start = min(unseen)
        orbit = {start}
        frontier = [start]
        while frontier:
            point = frontier.pop()
            for mapping in generators:
                image = mapping(point)
                require(image in unseen or image in orbit,
                        ("generator left fibre", point, image))
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        unseen.difference_update(orbit)
        orbits.append(tuple(sorted(orbit)))
    return tuple(sorted(orbits))


def apply_power(mapping, point, exponent: int):
    image = point
    for _ in range(exponent):
        image = mapping(image)
    return image


def affine_three_maps(prime: int, shallow_depth: int, deep_depth: int,
                      level: int):
    modulus = prime ** level
    shallow_scale = prime ** shallow_depth
    deep_scale = prime ** deep_depth

    def shallow(point):
        x_coord, y_coord, z_coord = point
        return (
            x_coord,
            (y_coord + shallow_scale) % modulus,
            (z_coord + shallow_scale * x_coord) % modulus,
        )

    def shallow_inverse(point):
        x_coord, y_coord, z_coord = point
        return (
            x_coord,
            (y_coord - shallow_scale) % modulus,
            (z_coord - shallow_scale * x_coord) % modulus,
        )

    def deep(point):
        x_coord, y_coord, z_coord = point
        return (
            (x_coord + deep_scale) % modulus,
            y_coord,
            z_coord,
        )

    def deep_inverse(point):
        x_coord, y_coord, z_coord = point
        return (
            (x_coord - deep_scale) % modulus,
            y_coord,
            z_coord,
        )

    def commutator(point):
        return deep_inverse(shallow_inverse(deep(shallow(point))))

    return shallow, deep, commutator


def affine_three_controls():
    rows = []
    right_action_checks = 0
    cases = (
        (3, 1, 2, (2, 3, 4)),
        (2, 2, 3, (3, 4, 5, 6)),
    )
    for prime, shallow_depth, deep_depth, levels in cases:
        total_depth = shallow_depth + deep_depth
        for level in levels:
            modulus = prime ** level
            shallow_scale = prime ** shallow_depth
            deep_scale = prime ** deep_depth
            points = fibre_points(
                prime, 3, shallow_depth, level, (0, 0, 0)
            )
            shallow, deep, commutator = affine_three_maps(
                prime, shallow_depth, deep_depth, level
            )
            expected_commutator = lambda point: (
                point[0],
                point[1],
                (point[2] - shallow_scale * deep_scale) % modulus,
            )
            require(
                all(commutator(point) == expected_commutator(point)
                    for point in points),
                ("unequal A3 commutator", prime, shallow_depth,
                 deep_depth, level),
            )
            commute = all(
                shallow(deep(point)) == deep(shallow(point))
                for point in points
            )
            require(
                commute == (level <= total_depth),
                ("unequal A3 abelian cutoff", prime, shallow_depth,
                 deep_depth, level),
            )
            orbits = generated_orbits(points, (shallow, deep))
            if level <= total_depth:
                group_order = prime ** (2 * level - total_depth)
            else:
                group_order = prime ** (3 * level - 2 * total_depth)
            require(
                {len(orbit) for orbit in orbits} == {group_order},
                ("unequal A3 orbit size", prime, shallow_depth,
                 deep_depth, level),
            )
            require(
                len(orbits) == len(points) // group_order,
                ("unequal A3 bank count", prime, shallow_depth,
                 deep_depth, level),
            )

            if level > total_depth:
                first_depth = level - shallow_depth
                second_depth = level - deep_depth
                central_depth = level - total_depth
                normal_images = set()
                for i, j, r in product(
                    range(prime ** first_depth),
                    range(prime ** second_depth),
                    range(prime ** central_depth),
                ):
                    normal_images.add((
                        deep_scale * j % modulus,
                        shallow_scale * i % modulus,
                        -shallow_scale * deep_scale * r % modulus,
                    ))
                require(
                    len(normal_images) == group_order,
                    ("unequal A3 normal uniqueness", prime, level),
                )
            rows.append((
                prime,
                shallow_depth,
                deep_depth,
                level,
                len(points),
                group_order,
                len(orbits),
                len(orbits[0]),
                commute,
            ))

    # Freeze the right-point convention and the sign of the mixed cocycle on
    # the first unequal dyadic mixed-modulus presentation.
    prime = 2
    shallow_depth = 1
    deep_depth = 2
    level = 4
    modulus = prime ** level
    shallow_scale = prime ** shallow_depth
    deep_scale = prime ** deep_depth
    first_depth = level - shallow_depth
    second_depth = level - deep_depth
    central_depth = level - shallow_depth - deep_depth
    elements = tuple(product(
        range(prime ** first_depth),
        range(prime ** second_depth),
        range(prime ** central_depth),
    ))
    probes = fibre_points(
        prime, 3, shallow_depth, level, (0, 0, 0)
    )[:16]

    def right_normal(point, element):
        i, j, r = element
        x_coord, y_coord, z_coord = point
        return (
            (x_coord + deep_scale * j) % modulus,
            (y_coord + shallow_scale * i) % modulus,
            (
                z_coord
                + shallow_scale * i * x_coord
                - shallow_scale * deep_scale * r
            ) % modulus,
        )

    for point, left, right in product(probes, elements, elements):
        product_element = mixed_mul(
            left, right, prime, first_depth, second_depth, central_depth
        )
        require(
            right_normal(right_normal(point, left), right)
            == right_normal(point, product_element),
            ("mixed right action", point, left, right),
        )
        right_action_checks += 1
    require(right_action_checks == 65536,
            ("right-action check count", right_action_checks))
    return tuple(rows), right_action_checks


def affine_four_boundary_controls():
    rows = []
    for prime, shallow_depth, deep_depth in (
        (2, 1, 2),
        (3, 1, 2),
        (2, 2, 3),
        (5, 1, 3),
    ):
        shallow_scale = prime ** shallow_depth
        deep_scale = prime ** deep_depth
        total_depth = shallow_depth + deep_depth
        boundary = total_depth + shallow_depth
        for level in (boundary, boundary + 1):
            modulus = prime ** level
            probes = tuple(product(range(3), repeat=4))

            def shallow(point):
                x_coord, y_coord, z_coord, w_coord = point
                return (
                    x_coord,
                    (y_coord + shallow_scale) % modulus,
                    (z_coord + shallow_scale * x_coord) % modulus,
                    (w_coord + shallow_scale * z_coord) % modulus,
                )

            def shallow_inverse(point):
                x_coord, y_coord, z_coord, w_coord = point
                return (
                    x_coord,
                    (y_coord - shallow_scale) % modulus,
                    (z_coord - shallow_scale * x_coord) % modulus,
                    (
                        w_coord
                        - shallow_scale * z_coord
                        + shallow_scale * shallow_scale * x_coord
                    ) % modulus,
                )

            def deep(point):
                x_coord, y_coord, z_coord, w_coord = point
                return (
                    (x_coord + deep_scale) % modulus,
                    y_coord,
                    z_coord,
                    w_coord,
                )

            def deep_inverse(point):
                x_coord, y_coord, z_coord, w_coord = point
                return (
                    (x_coord - deep_scale) % modulus,
                    y_coord,
                    z_coord,
                    w_coord,
                )

            def commutator(point):
                return deep_inverse(
                    shallow_inverse(deep(shallow(point)))
                )

            expected = lambda point: (
                point[0],
                point[1],
                (point[2] - shallow_scale * deep_scale) % modulus,
                (
                    point[3]
                    + shallow_scale * shallow_scale * deep_scale
                ) % modulus,
            )
            require(
                all(commutator(point) == expected(point)
                    for point in probes),
                ("A4 exact commutator", prime, shallow_depth,
                 deep_depth, level),
            )
            require(
                all(commutator(deep(point)) == deep(commutator(point))
                    for point in probes),
                ("A4 deep centrality", prime, shallow_depth,
                 deep_depth, level),
            )
            defects = {
                (
                    commutator(shallow(point))[3]
                    - shallow(commutator(point))[3]
                ) % modulus
                for point in probes
            }
            expected_defect = (
                shallow_scale * shallow_scale * deep_scale
            ) % modulus
            require(
                defects == {expected_defect},
                ("A4 shallow triple carry", prime, shallow_depth,
                 deep_depth, level, defects),
            )
            require(
                (defects == {0}) == (level <= boundary),
                ("A4 sharp boundary", prime, shallow_depth,
                 deep_depth, level),
            )
            if level == boundary + 1:
                rows.append((
                    prime,
                    shallow_depth,
                    deep_depth,
                    level,
                    len(probes),
                    shallow_scale * shallow_scale * deep_scale,
                    tuple(sorted(defects)),
                ))
    return tuple(rows)


def rank_mod_two(vectors):
    rows = [list(vector) for vector in vectors]
    rank = 0
    width = len(rows[0]) if rows else 0
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(rows))
             if rows[index][column] % 2),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for index in range(len(rows)):
            if index != rank and rows[index][column] % 2:
                rows[index] = [
                    (left + right) % 2
                    for left, right in zip(rows[index], rows[rank])
                ]
        rank += 1
    return rank


def dyadic_positive_controls():
    field_rows = []
    for y_residue in (0, 1):
        delta_shallow = (0, 1, y_residue)
        sigma_shallow = (0, 1, (y_residue + 1) % 2)
        delta_deep = (y_residue, 1, 0)
        bracket = (1, 0, 1)
        ordinary_rank = rank_mod_two((
            delta_shallow, delta_deep, bracket
        ))
        pair_rank = rank_mod_two((sigma_shallow, delta_deep))
        repaired_rank = rank_mod_two((
            sigma_shallow, delta_deep, bracket
        ))
        require(delta_shallow != (0, 0, 0), "zero shallow carry")
        require(ordinary_rank == 2, ("ordinary triple", y_residue))
        require(pair_rank == 2, ("repaired pair", y_residue))
        require(repaired_rank == 3, ("repaired triple", y_residue))
        field_rows.append((
            y_residue,
            delta_shallow,
            sigma_shallow,
            delta_deep,
            bracket,
            ordinary_rank,
            pair_rank,
            repaired_rank,
        ))

    orbit_rows = []
    for deep_depth in (2, 3):
        deep_scale = 2 ** deep_depth
        total_depth = deep_depth + 1
        for level in range(2, deep_depth + 3):
            modulus = 2 ** level

            def shallow(point):
                x_coord, y_coord, z_coord = point
                return (
                    x_coord,
                    (y_coord + 2) % modulus,
                    (z_coord + 2 * y_coord) % modulus,
                )

            def deep(point):
                x_coord, y_coord, z_coord = point
                return (
                    (x_coord + deep_scale * y_coord) % modulus,
                    (y_coord + deep_scale) % modulus,
                    z_coord,
                )

            if level <= deep_depth:
                group_order = 2 ** (level - 1)
            elif level <= total_depth:
                group_order = 2 ** (2 * level - total_depth)
            else:
                group_order = 2 ** (3 * level - 2 * total_depth)
            bank_counts = set()
            orbit_sizes = set()
            fibre_size = None
            for base in product(range(2), repeat=3):
                points = fibre_points(2, 3, 1, level, base)
                orbits = generated_orbits(points, (shallow, deep))
                fibre_size = len(points)
                bank_counts.add(len(orbits))
                orbit_sizes.update(len(orbit) for orbit in orbits)
            expected_banks = fibre_size // group_order
            require(
                orbit_sizes == {group_order},
                ("dyadic positive orbit size", deep_depth, level,
                 orbit_sizes),
            )
            require(
                bank_counts == {expected_banks},
                ("dyadic positive banks", deep_depth, level, bank_counts),
            )
            orbit_rows.append((
                deep_depth,
                level,
                8,
                fibre_size,
                group_order,
                expected_banks,
                tuple(sorted(orbit_sizes)),
            ))
    return tuple(field_rows), tuple(orbit_rows)


def dyadic_involution_hostile():
    rows = []
    for deep_depth in (2, 3):
        level = deep_depth + 2
        modulus = 2 ** level
        deep_scale = 2 ** deep_depth
        total_depth = deep_depth + 1
        points = tuple(
            (1 + 2 * i, 2 * j, 2 * k)
            for i, j, k in product(
                range(2 ** (level - 1)), repeat=3
            )
        )

        def shallow(point):
            x_coord, y_coord, z_coord = point
            return ((-x_coord) % modulus, y_coord, z_coord)

        def deep(point):
            x_coord, y_coord, z_coord = point
            return (
                x_coord,
                (y_coord + deep_scale) % modulus,
                (z_coord + deep_scale * x_coord) % modulus,
            )

        def deep_inverse(point):
            x_coord, y_coord, z_coord = point
            return (
                x_coord,
                (y_coord - deep_scale) % modulus,
                (z_coord - deep_scale * x_coord) % modulus,
            )

        def commutator(point):
            return deep_inverse(shallow(deep(shallow(point))))

        require(
            all(shallow(shallow(point)) == point for point in points),
            ("involution square", deep_depth),
        )
        require(
            all(
                commutator(point)
                == (
                    point[0],
                    point[1],
                    (point[2] - 2 * deep_scale * point[0]) % modulus,
                )
                for point in points
            ),
            ("involution commutator", deep_depth),
        )
        orbits = generated_orbits(points, (shallow, deep))
        require(
            {len(orbit) for orbit in orbits} == {16},
            ("involution image orbit", deep_depth),
        )
        nominal_order = 2 ** (3 * level - 2 * total_depth)

        first_modulus = 2 ** (level - 1)
        second_modulus = 2 ** (level - deep_depth)
        central_modulus = 2 ** (level - total_depth)
        base = (1, 0, 0)
        normal_images = []
        for i, j, r in product(
            range(first_modulus),
            range(second_modulus),
            range(central_modulus),
        ):
            image = apply_power(shallow, base, i)
            image = apply_power(deep, image, j)
            image = apply_power(commutator, image, r)
            normal_images.append(image)
        distinct_images = len(set(normal_images))
        fixed_words = sum(image == base for image in normal_images)
        require(distinct_images == 16,
                ("involution normal image", deep_depth, distinct_images))
        require(fixed_words == 2 ** deep_depth,
                ("involution kernel", deep_depth, fixed_words))
        require(nominal_order == len(normal_images),
                ("involution nominal order", deep_depth))
        rows.append((
            deep_depth,
            level,
            len(points),
            nominal_order,
            distinct_images,
            fixed_words,
            len(orbits),
            len(orbits[0]),
        ))

    carries = (
        ("delta_shallow", (1, 0, 0)),
        ("sigma_shallow", (0, 0, 0)),
        ("delta_deep", (0, 1, 1)),
        ("bracket", (0, 0, 1)),
        ("ordinary_rank", 3),
        ("repaired_rank", 2),
    )
    return tuple(rows), carries


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "assert node found",
    )
    forbidden_calls = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        and node.func.id in {"eval", "exec", "compile"}
    }
    require(not forbidden_calls, ("forbidden calls", forbidden_calls))
    for path, expected_hash in DEPENDENCIES.items():
        actual_hash = lf_sha256(path)
        require(actual_hash == expected_hash,
                ("dependency changed", path.name, actual_hash))

    presentation = presentation_controls()
    affine_three, right_action_checks = affine_three_controls()
    affine_four = affine_four_boundary_controls()
    dyadic_fields, dyadic_positive = dyadic_positive_controls()
    dyadic_hostile = dyadic_involution_hostile()
    semantic_payload = (
        presentation,
        affine_three,
        right_action_checks,
        affine_four,
        dyadic_fields,
        dyadic_positive,
        dyadic_hostile,
    )
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic_sha256),
        )

    print("unequal-depth noncommuting Hensel--Heisenberg controls")
    print(f"mixed_presentation_rows={presentation}")
    print(f"unequal_A3_rows={affine_three}")
    print(f"right_action_checks={right_action_checks}")
    print(f"swapped_A4_boundary_rows={affine_four}")
    print(f"dyadic_field_rows={dyadic_fields}")
    print(f"dyadic_positive_rows={dyadic_positive}")
    print(f"dyadic_involution_hostile={dyadic_hostile}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
