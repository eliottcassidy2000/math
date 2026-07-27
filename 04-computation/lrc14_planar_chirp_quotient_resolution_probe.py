#!/usr/bin/env python3
"""Exact finite-affine-plane probe for the planar-chirp resolution ladder."""

from __future__ import annotations

P = 13
Point = tuple[int, int]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dot(a: Point, x: Point) -> int:
    return (a[0] * x[0] + a[1] * x[1]) % P


def line_sum(function: dict[Point, int], normal: Point, level: int) -> int:
    return sum(value for x, value in function.items() if dot(normal, x) == level)


def main() -> None:
    points = [(x, y) for x in range(P) for y in range(P)]
    directions = [(1, slope) for slope in range(P)] + [(0, 1)]

    # Each projective dual direction partitions F_p^2 into p affine lines
    # of p points, and any two points lie together in exactly one direction.
    for normal in directions:
        fibres = [
            [x for x in points if dot(normal, x) == level] for level in range(P)
        ]
        require(all(len(fibre) == P for fibre in fibres), "wrong affine fibre size")
        require(
            {x for fibre in fibres for x in fibre} == set(points),
            "affine fibres do not partition the plane",
        )

    pair_direction_checks = 0
    for x in points:
        for y in points:
            if x == y:
                continue
            common = sum(dot(a, x) == dot(a, y) for a in directions)
            require(common == 1, f"nonunique affine direction for {x},{y}")
            pair_direction_checks += 1

    # Exact affine-plane reconstruction:
    # sum over the p+1 lines through x = p*f(x)+sum_y f(y).
    function = {
        x: ((17 * x[0] + 23 * x[1] + 5 * x[0] * x[1]) % 31) - 15
        for x in points
    }
    total = sum(function.values())
    reconstructed = {}
    for x in points:
        through_x = sum(
            line_sum(function, normal, dot(normal, x)) for normal in directions
        )
        numerator = through_x - total
        require(numerator % P == 0, "affine reconstruction lost integrality")
        reconstructed[x] = numerator // P
    require(reconstructed == function, "all-direction affine reconstruction failed")

    # If one dual direction is missing, a nonzero zero-mean function of its
    # coordinate is invisible to every retained line-sum family.
    missing = (0, 1)
    omitted_hostile = {
        x: (1 if dot(missing, x) == 0 else -1 if dot(missing, x) == 1 else 0)
        for x in points
    }
    require(any(omitted_hostile.values()), "omitted-direction hostile vanished")
    retained = [a for a in directions if a != missing]
    for normal in retained:
        require(
            all(
                line_sum(omitted_hostile, normal, level) == 0
                for level in range(P)
            ),
            f"omitted-direction hostile visible to {normal}",
        )
    require(
        any(
            line_sum(omitted_hostile, missing, level) != 0
            for level in range(P)
        ),
        "omitted-direction hostile invisible even in its own direction",
    )

    # Two transverse jet directions are much weaker: the checkerboard has
    # zero row and column sums but is nonzero.
    checkerboard = {x: 0 for x in points}
    checkerboard[(0, 0)] = 1
    checkerboard[(0, 1)] = -1
    checkerboard[(1, 0)] = -1
    checkerboard[(1, 1)] = 1
    for normal in ((1, 0), (0, 1)):
        require(
            all(
                line_sum(checkerboard, normal, level) == 0
                for level in range(P)
            ),
            "two-direction checkerboard became visible",
        )
    require(any(checkerboard.values()), "checkerboard hostile vanished")

    # This is not merely an arbitrary edge field: at one nonzero lag it is
    # realized by a genuine scalar signal on eight disjoint sites.
    h = (3, 3)
    signal: dict[Point, int] = {}
    for start, value in checkerboard.items():
        if value == 0:
            continue
        endpoint = ((start[0] + h[0]) % P, (start[1] + h[1]) % P)
        require(start not in signal and endpoint not in signal, "support collision")
        signal[start] = 1
        signal[endpoint] = value
    edge_field = {
        u: signal.get(((u[0] + h[0]) % P, (u[1] + h[1]) % P), 0)
        * signal.get(u, 0)
        for u in points
    }
    require(edge_field == checkerboard, "checkerboard is not a genuine lag field")

    print("LRC14 planar-chirp quotient-resolution exact probe")
    print(f"affine plane: F_{P}^2, points={len(points)}, directions={len(directions)}")
    print(f"ordered distinct point/direction checks: {pair_direction_checks}")
    print("all 14 directions: exact point reconstruction PASS")
    print("13 directions: omitted-direction zero-mean hostile PASS")
    print("2 transverse directions: genuine-signal checkerboard hostile PASS")
    print("VERDICT: one jet line yields 13-point derivative-coset sums, not atoms")


if __name__ == "__main__":
    main()
