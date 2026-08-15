#!/usr/bin/env python3
"""Exact controls for THM-3449's noncommuting Hensel--Heisenberg law."""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import product
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "fa4dbd16a9f56fb350388fde53c8b3ef39cdfbc6f3092cd31aed83c506bc8154"
ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / "01-canon/theorems/THM-3442-smooth-hensel-fibre-vector-field-orbit-law.md":
        "5ca39bf0693f7d314fd2907d113eeb8aaca613a5f6e02171cd5a1a4700b3a11a",
    ROOT / "01-canon/theorems/THM-3444-commuting-smooth-hensel-vector-field-lattice-action.md":
        "d6128768b4591273577c351ddad1977506667204ea39008ff24f741d77a4fb32",
    ROOT / "01-canon/theorems/THM-3446-weighted-depth-commuting-hensel-lattice-action.md":
        "b64e56f02accabd0efd63b16172bb96316184ebeaf3738f5302e71083a62f658",
}


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def hp_mul(left, right, prime: int, horizontal_depth: int,
           central_depth: int):
    """Multiplication for x^i y^j z^k, z=xyx^-1y^-1."""
    horizontal_modulus = prime ** horizontal_depth
    central_modulus = prime ** central_depth
    i, j, k = left
    ii, jj, kk = right
    return (
        (i + ii) % horizontal_modulus,
        (j + jj) % horizontal_modulus,
        (k + kk - j * ii) % central_modulus,
    )


def hp_inverse(element, prime: int, horizontal_depth: int,
               central_depth: int):
    horizontal_modulus = prime ** horizontal_depth
    central_modulus = prime ** central_depth
    i, j, k = element
    return (
        (-i) % horizontal_modulus,
        (-j) % horizontal_modulus,
        (-k - i * j) % central_modulus,
    )


def hp_power(element, exponent: int, prime: int, horizontal_depth: int,
             central_depth: int):
    answer = (0, 0, 0)
    base = element
    while exponent:
        if exponent & 1:
            answer = hp_mul(
                answer, base, prime, horizontal_depth, central_depth
            )
        base = hp_mul(base, base, prime, horizontal_depth, central_depth)
        exponent //= 2
    return answer


def presentation_controls():
    rows = []
    for prime, horizontal_depth, central_depth in (
        (2, 2, 1),
        (3, 2, 1),
        (2, 3, 2),
    ):
        horizontal_modulus = prime ** horizontal_depth
        central_modulus = prime ** central_depth
        elements = tuple(product(
            range(horizontal_modulus),
            range(horizontal_modulus),
            range(central_modulus),
        ))
        identity = (0, 0, 0)
        x = (1, 0, 0)
        y = (0, 1, 0)
        z = (0, 0, 1)
        require(
            len(elements)
            == prime ** (2 * horizontal_depth + central_depth),
            ("presentation order", prime, horizontal_depth, central_depth),
        )
        for element in elements:
            inverse = hp_inverse(
                element, prime, horizontal_depth, central_depth
            )
            require(
                hp_mul(element, inverse, prime, horizontal_depth,
                       central_depth) == identity,
                ("right inverse", element),
            )
            require(
                hp_mul(inverse, element, prime, horizontal_depth,
                       central_depth) == identity,
                ("left inverse", element),
            )
        x_inverse = hp_inverse(x, prime, horizontal_depth, central_depth)
        y_inverse = hp_inverse(y, prime, horizontal_depth, central_depth)
        commutator = hp_mul(
            hp_mul(
                hp_mul(x, y, prime, horizontal_depth, central_depth),
                x_inverse,
                prime,
                horizontal_depth,
                central_depth,
            ),
            y_inverse,
            prime,
            horizontal_depth,
            central_depth,
        )
        require(commutator == z, ("commutator", commutator, z))
        require(
            hp_power(x, horizontal_modulus, prime, horizontal_depth,
                     central_depth) == identity,
            "x order bound",
        )
        require(
            hp_power(x, horizontal_modulus // prime, prime, horizontal_depth,
                     central_depth) != identity,
            "x exact order",
        )
        require(
            hp_power(y, horizontal_modulus, prime, horizontal_depth,
                     central_depth) == identity,
            "y order bound",
        )
        require(
            hp_power(z, central_modulus, prime, horizontal_depth,
                     central_depth) == identity,
            "z order bound",
        )
        require(
            hp_power(z, central_modulus // prime, prime, horizontal_depth,
                     central_depth) != identity,
            "z exact order",
        )

        # The bilinear cocycle identity is the associativity proof on the
        # complete horizontal universe; central coordinates then cancel.
        horizontal = tuple(product(range(horizontal_modulus), repeat=2))
        for first, second, third in product(horizontal, repeat=3):
            first_second = (
                (first[0] + second[0]) % horizontal_modulus,
                (first[1] + second[1]) % horizontal_modulus,
            )
            second_third = (
                (second[0] + third[0]) % horizontal_modulus,
                (second[1] + third[1]) % horizontal_modulus,
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
                ("cocycle", first, second, third),
            )

        center_size = sum(
            1
            for i, j, _ in elements
            if i % central_modulus == 0 and j % central_modulus == 0
        )
        expected_center = prime ** (
            2 * (horizontal_depth - central_depth) + central_depth
        )
        require(center_size == expected_center, ("center", center_size))
        rows.append((
            prime,
            horizontal_depth,
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
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        unseen.difference_update(orbit)
        orbits.append(tuple(sorted(orbit)))
    return tuple(sorted(orbits))


def standard_heisenberg_controls():
    rows = []
    right_action_checks = 0
    for prime, carry_depth, level in (
        (2, 1, 2),
        (2, 1, 3),
        (2, 1, 4),
        (3, 1, 2),
        (3, 1, 3),
        (2, 2, 4),
        (2, 2, 5),
    ):
        modulus = prime ** level
        scale = prime ** carry_depth
        width = prime ** (level - carry_depth)
        points = fibre_points(prime, 3, carry_depth, level, (0, 0, 0))

        def g(point):
            x, y, z_coord = point
            return ((x + scale) % modulus, y, z_coord)

        def g_inverse(point):
            x, y, z_coord = point
            return ((x - scale) % modulus, y, z_coord)

        def h(point):
            x, y, z_coord = point
            return (x, (y + scale) % modulus,
                    (z_coord + scale * x) % modulus)

        def h_inverse(point):
            x, y, z_coord = point
            return (x, (y - scale) % modulus,
                    (z_coord - scale * x) % modulus)

        def commutator(point):
            return h_inverse(g_inverse(h(g(point))))

        def expected_commutator(point):
            x, y, z_coord = point
            return (x, y, (z_coord + scale * scale) % modulus)

        require(
            all(commutator(point) == expected_commutator(point)
                for point in points),
            ("A3 commutator", prime, carry_depth, level),
        )
        commute = all(g(h(point)) == h(g(point)) for point in points)
        require(
            commute == (level <= 2 * carry_depth),
            ("A3 abelian cutoff", prime, carry_depth, level),
        )
        orbits = generated_orbits(points, (g, h))
        horizontal_depth = level - carry_depth
        central_depth = max(0, level - 2 * carry_depth)
        group_order = prime ** (2 * horizontal_depth + central_depth)
        require(
            {len(orbit) for orbit in orbits} == {group_order},
            ("A3 orbit size", prime, carry_depth, level),
        )
        require(
            len(orbits) == len(points) // group_order,
            ("A3 orbit banks", prime, carry_depth, level),
        )

        # Pullback multiplication induces a right action on points.  Starting
        # from zero, the normal word G^i H^j K^k has these pairwise-distinct
        # images; the ij term is the right-action word-area correction.
        normal_images = {
            (
                scale * i % modulus,
                scale * j % modulus,
                scale * scale * (i * j + k) % modulus,
            )
            for i, j, k in product(
                range(prime ** horizontal_depth),
                range(prime ** horizontal_depth),
                range(prime ** central_depth),
            )
        }
        require(
            len(normal_images) == group_order,
            ("A3 normal uniqueness", prime, carry_depth, level),
        )

        if (prime, carry_depth, level) == (2, 1, 3):
            elements = tuple(product(range(4), range(4), range(2)))

            def right_normal(point, element):
                i, j, r = element
                x, y, z_coord = point
                moved_x = (x + scale * i) % modulus
                return (
                    moved_x,
                    (y + scale * j) % modulus,
                    (
                        z_coord
                        + scale * j * moved_x
                        + scale * scale * r
                    ) % modulus,
                )

            for point, left, right in product(points, elements, elements):
                product_element = hp_mul(left, right, 2, 2, 1)
                require(
                    right_normal(right_normal(point, left), right)
                    == right_normal(point, product_element),
                    ("right action", point, left, right),
                )
                right_action_checks += 1
        rows.append((
            prime,
            carry_depth,
            level,
            len(points),
            group_order,
            len(orbits),
            len(orbits[0]),
        ))
    require(right_action_checks == 64 * 32 * 32,
            ("right-action check count", right_action_checks))
    return tuple(rows), right_action_checks


def step_three_controls():
    orbit_rows = []
    for prime, carry_depth, level in (
        (2, 1, 2),
        (2, 1, 3),
        (3, 1, 2),
        (3, 1, 3),
    ):
        modulus = prime ** level
        scale = prime ** carry_depth
        points = fibre_points(prime, 4, carry_depth, level, (0, 0, 0, 0))

        def g(point):
            x, y, z_coord, w = point
            return ((x + scale) % modulus, y, z_coord, w)

        def g_inverse(point):
            x, y, z_coord, w = point
            return ((x - scale) % modulus, y, z_coord, w)

        def h(point):
            x, y, z_coord, w = point
            return (
                x,
                (y + scale) % modulus,
                (z_coord + scale * x) % modulus,
                (w + scale * z_coord) % modulus,
            )

        def h_inverse(point):
            x, y, z_coord, w = point
            return (
                x,
                (y - scale) % modulus,
                (z_coord - scale * x) % modulus,
                (w - scale * z_coord + scale * scale * x) % modulus,
            )

        def commutator(point):
            return h_inverse(g_inverse(h(g(point))))

        def expected_commutator(point):
            x, y, z_coord, w = point
            return (
                x,
                y,
                (z_coord + scale * scale) % modulus,
                (w - scale ** 3) % modulus,
            )

        require(
            all(commutator(point) == expected_commutator(point)
                for point in points),
            ("A4 commutator", prime, carry_depth, level),
        )
        require(
            all(commutator(g(point)) == g(commutator(point))
                for point in points),
            ("A4 g centrality", prime, carry_depth, level),
        )
        require(
            all(commutator(h(point)) == h(commutator(point))
                for point in points),
            ("A4 h centrality", prime, carry_depth, level),
        )
        orbits = generated_orbits(points, (g, h))
        horizontal_depth = level - carry_depth
        central_depth = max(0, level - 2 * carry_depth)
        group_order = prime ** (2 * horizontal_depth + central_depth)
        require(
            {len(orbit) for orbit in orbits} == {group_order},
            ("A4 orbit size", prime, carry_depth, level),
        )
        orbit_rows.append((
            prime,
            carry_depth,
            level,
            len(points),
            group_order,
            len(orbits),
            len(orbits[0]),
        ))

    boundary_rows = []
    for prime, carry_depth in ((2, 1), (3, 1), (5, 1), (2, 2)):
        level = 3 * carry_depth + 1
        modulus = prime ** level
        scale = prime ** carry_depth
        width = prime ** (level - carry_depth)
        offsets = range(min(width, 3))
        probes = tuple(
            tuple(scale * value % modulus for value in offset)
            for offset in product(offsets, repeat=4)
        )

        def g(point):
            x, y, z_coord, w = point
            return ((x + scale) % modulus, y, z_coord, w)

        def g_inverse(point):
            x, y, z_coord, w = point
            return ((x - scale) % modulus, y, z_coord, w)

        def h(point):
            x, y, z_coord, w = point
            return (
                x,
                (y + scale) % modulus,
                (z_coord + scale * x) % modulus,
                (w + scale * z_coord) % modulus,
            )

        def h_inverse(point):
            x, y, z_coord, w = point
            return (
                x,
                (y - scale) % modulus,
                (z_coord - scale * x) % modulus,
                (w - scale * z_coord + scale * scale * x) % modulus,
            )

        def commutator(point):
            return h_inverse(g_inverse(h(g(point))))

        require(
            all(
                commutator(point)
                == (
                    point[0],
                    point[1],
                    (point[2] + scale * scale) % modulus,
                    (point[3] - scale ** 3) % modulus,
                )
                for point in probes
            ),
            ("A4 boundary commutator", prime, carry_depth),
        )
        require(
            all(commutator(g(point)) == g(commutator(point))
                for point in probes),
            ("A4 boundary g central", prime, carry_depth),
        )
        h_defects = {
            (h(commutator(point))[3] - commutator(h(point))[3]) % modulus
            for point in probes
        }
        require(
            h_defects == {scale ** 3 % modulus},
            ("A4 first triple carry", prime, carry_depth, h_defects),
        )
        boundary_rows.append((
            prime,
            carry_depth,
            level,
            len(probes),
            scale ** 3,
            tuple(sorted(h_defects)),
        ))
    return tuple(orbit_rows), tuple(boundary_rows)


def dyadic_square_hostile():
    rows = []
    for level in (2, 3):
        modulus = 2 ** level
        scale = 2
        width = 2 ** (level - 1)
        points = tuple(
            (1 + scale * i, scale * j, scale * k)
            for i, j, k in product(range(width), repeat=3)
        )

        def g(point):
            x, y, z_coord = point
            return ((-x) % modulus, y, z_coord)

        def h(point):
            x, y, z_coord = point
            return (x, (y + scale) % modulus,
                    (z_coord + scale * x) % modulus)

        def h_inverse(point):
            x, y, z_coord = point
            return (x, (y - scale) % modulus,
                    (z_coord - scale * x) % modulus)

        def commutator(point):
            return h_inverse(g(h(g(point))))

        require(
            all(
                commutator(point)
                == (point[0], point[1],
                    (point[2] - 4 * point[0]) % modulus)
                for point in points
            ),
            ("dyadic hostile commutator", level),
        )
        require(
            all(g(g(point)) == point for point in points),
            ("dyadic hostile square", level),
        )
        orbits = generated_orbits(points, (g, h))
        expected_orbit = 4 if level == 2 else 16
        require(
            {len(orbit) for orbit in orbits} == {expected_orbit},
            ("dyadic hostile orbit", level),
        )
        rows.append((
            level,
            len(points),
            expected_orbit,
            len(orbits),
        ))
    carries = (
        ("delta_g", (1, 0, 0)),
        ("delta_h", (0, 1, 1)),
        ("bracket", (0, 0, 1)),
        ("sigma_g", (0, 0, 0)),
        ("sigma_h", (0, 1, 1)),
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
        require(lf_sha256(path) == expected_hash,
                ("dependency changed", path.name, lf_sha256(path)))

    presentation = presentation_controls()
    standard, right_action_checks = standard_heisenberg_controls()
    step_three_orbits, step_three_boundary = step_three_controls()
    dyadic = dyadic_square_hostile()
    semantic_payload = (
        presentation,
        standard,
        right_action_checks,
        step_three_orbits,
        step_three_boundary,
        dyadic,
    )
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic_sha256),
        )

    print("noncommuting smooth Hensel--Heisenberg controls")
    print(f"presentation_rows={presentation}")
    print(f"standard_A3_rows={standard}")
    print(f"right_action_checks={right_action_checks}")
    print(f"step3_A4_orbit_rows={step_three_orbits}")
    print(f"step3_A4_boundary_rows={step_three_boundary}")
    print(f"dyadic_square_hostile={dyadic}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
