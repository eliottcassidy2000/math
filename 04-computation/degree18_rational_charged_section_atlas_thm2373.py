#!/usr/bin/env python3
"""Exact companion for THM-2373's rational charged-section atlas."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


ZERO = Fraction(0)
WEIGHTS = {
    "y": 1,
    "u": 2,
    "B": 2,
    "C": 3,
    "D": 4,
    "W": 5,
}
ADJACENT = (("B", "C"), ("C", "D"), ("D", "W"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def scale_point(
    point: dict[str, Fraction],
    scalar: Fraction,
) -> dict[str, Fraction]:
    return {
        coordinate: value * scalar ** WEIGHTS[coordinate]
        for coordinate, value in point.items()
    }


def chart_scalar(
    point: dict[str, Fraction],
    chart: tuple[str, str],
) -> Fraction:
    left, right = chart
    require(point[left] != 0 and point[right] != 0, "chart is not live")
    return point[left] / point[right]


def normalized_point(
    point: dict[str, Fraction],
    chart: tuple[str, str],
) -> dict[str, Fraction]:
    return scale_point(point, chart_scalar(point, chart))


def g0(point: dict[str, Fraction]) -> Fraction:
    y = point["y"]
    u = point["u"]
    B = point["B"]
    C = point["C"]
    D = point["D"]
    W = point["W"]
    return (
        -5878656 * W * y
        - 26040609 * u**3
        + 49601160 * B * u**2
        + 1607445 * u**2 * y**2
        - 20995200 * B**2 * u
        - 2857680 * B * u * y**2
        - 52907904 * D * u
        - 138915 * u * y**4
        + 777600 * B**2 * y**2
        + 33592320 * B * D
        - 5598720 * B * C * y
        + 78120 * B * y**4
        + 1959552 * D * y**2
        - 435456 * C * y**3
        + 1127 * y**6
    )


def main() -> None:
    coordinates = ("B", "C", "D", "W")

    # Every support of size at least three meets one adjacent chart.
    support_checks = 0
    for support_size in (3, 4):
        for support in combinations(coordinates, support_size):
            live_charts = [
                chart
                for chart in ADJACENT
                if chart[0] in support and chart[1] in support
            ]
            require(live_charts, "higher support escaped every chart")
            support_checks += 1

    # Exact rational samples verify normalizations, covariance, and transitions.
    samples = 73
    normalization_checks = 0
    covariance_checks = 0
    transition_checks = 0
    invariant_checks = 0
    for seed in range(1, samples + 1):
        point = {
            "y": Fraction((3 * seed) % 17 + 1, 19),
            "u": Fraction((5 * seed) % 23 + 1, 29),
            "B": Fraction((7 * seed) % 31 + 1, 37),
            "C": Fraction((11 * seed) % 41 + 1, 43),
            "D": Fraction((13 * seed) % 47 + 1, 53),
            "W": Fraction((17 * seed) % 59 + 1, 61),
        }
        scalar = Fraction((19 * seed) % 67 + 1, 71)
        scaled = scale_point(point, scalar)
        require(g0(scaled) == scalar**6 * g0(point), "G0 covariance changed")
        covariance_checks += 1

        normalized: dict[tuple[str, str], dict[str, Fraction]] = {}
        scalars: dict[tuple[str, str], Fraction] = {}
        for chart in ADJACENT:
            chart_point = normalized_point(point, chart)
            left, right = chart
            require(
                chart_point[left] == chart_point[right],
                "charged section did not reach its equality slice",
            )
            normalized[chart] = chart_point
            scalars[chart] = chart_scalar(point, chart)
            normalization_checks += 1

        for first_index in range(len(ADJACENT)):
            for second_index in range(first_index + 1, len(ADJACENT)):
                first = ADJACENT[first_index]
                second = ADJACENT[second_index]
                transition = scalars[second] / scalars[first]
                require(
                    scale_point(normalized[first], transition)
                    == normalized[second],
                    "chart transition left the weighted orbit",
                )
                transition_checks += 1

        transition_invariants = (
            point["B"] * point["D"] / point["C"] ** 2,
            point["C"] * point["W"] / point["D"] ** 2,
            point["B"] * point["W"] / (point["C"] * point["D"]),
        )
        scaled_invariants = (
            scaled["B"] * scaled["D"] / scaled["C"] ** 2,
            scaled["C"] * scaled["W"] / scaled["D"] ** 2,
            scaled["B"] * scaled["W"] / (scaled["C"] * scaled["D"]),
        )
        require(
            transition_invariants == scaled_invariants
            and transition_invariants[0] * transition_invariants[1]
            == transition_invariants[2],
            "weight-zero transition invariant changed",
        )
        invariant_checks += 3

    # Directly verify the displayed closed formulas on one symbolic-rational point.
    point = {
        "y": Fraction(2, 3),
        "u": Fraction(3, 5),
        "B": Fraction(5, 7),
        "C": Fraction(7, 11),
        "D": Fraction(11, 13),
        "W": Fraction(13, 17),
    }
    bc = normalized_point(point, ("B", "C"))
    require(
        bc
        == {
            "y": point["B"] * point["y"] / point["C"],
            "u": point["B"] ** 2 * point["u"] / point["C"] ** 2,
            "B": point["B"] ** 3 / point["C"] ** 2,
            "C": point["B"] ** 3 / point["C"] ** 2,
            "D": point["B"] ** 4 * point["D"] / point["C"] ** 4,
            "W": point["B"] ** 5 * point["W"] / point["C"] ** 5,
        },
        "BC closed section formula changed",
    )
    cd = normalized_point(point, ("C", "D"))
    require(
        cd
        == {
            "y": point["C"] * point["y"] / point["D"],
            "u": point["C"] ** 2 * point["u"] / point["D"] ** 2,
            "B": point["C"] ** 2 * point["B"] / point["D"] ** 2,
            "C": point["C"] ** 4 / point["D"] ** 3,
            "D": point["C"] ** 4 / point["D"] ** 3,
            "W": point["C"] ** 5 * point["W"] / point["D"] ** 5,
        },
        "CD closed section formula changed",
    )
    dw = normalized_point(point, ("D", "W"))
    require(
        dw
        == {
            "y": point["D"] * point["y"] / point["W"],
            "u": point["D"] ** 2 * point["u"] / point["W"] ** 2,
            "B": point["D"] ** 2 * point["B"] / point["W"] ** 2,
            "C": point["D"] ** 3 * point["C"] / point["W"] ** 3,
            "D": point["D"] ** 5 / point["W"] ** 4,
            "W": point["D"] ** 5 / point["W"] ** 4,
        },
        "DW closed section formula changed",
    )

    # The B,D-only boundary has a mu_2 stabilizer and no weight-one Laurent monomial.
    hostile_solutions = [
        (left_power, right_power)
        for left_power in range(-100, 101)
        for right_power in range(-100, 101)
        if 2 * left_power + 4 * right_power == 1
    ]
    hostile = {
        "y": ZERO,
        "u": ZERO,
        "B": Fraction(2),
        "C": ZERO,
        "D": Fraction(3),
        "W": ZERO,
    }
    require(
        not hostile_solutions
        and scale_point(hostile, Fraction(-1)) == hostile
        and all(
            hostile[left] == 0 or hostile[right] == 0
            for left, right in ADJACENT
        ),
        "two-sparse stabilizer hostile changed",
    )

    print("THM-2373 degree-eighteen charged-section atlas exact referee")
    print(f"higher-support chart-cover subsets: {support_checks}")
    print(f"exact rational samples: {samples}")
    print(f"G0 weighted-covariance checks: {covariance_checks}")
    print(f"chart normalizations: {normalization_checks}")
    print(f"chart transitions: {transition_checks}")
    print(f"weight-zero invariant checks: {invariant_checks}")
    print("closed BC/CD/DW section formulas: PASS / PASS / PASS")
    print("B,D-only mu2/no-weight-one hostile: PASS")
    print(
        "VERDICT: every higher-support spectral orbit meets a "
        "root-free rational charged chart section"
    )


if __name__ == "__main__":
    main()
