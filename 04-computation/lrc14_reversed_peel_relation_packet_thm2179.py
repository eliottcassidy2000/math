#!/usr/bin/env python3
"""Exact referee for THM-2179.

The all-row argument in THM-2179 is Fourier/Jackson analysis.  This companion
checks every rational constant in that argument, verifies the adjacent
Jackson cutoff, and replays the named reversed-peel hostile row by two exact
interval algorithms.  Assertions remain active under ``python -O``.
"""

from fractions import Fraction


RADIUS = Fraction(3, 41)
DANGER_MASS = 2 * RADIUS
SAFE_MASS = 1 - DANGER_MASS
JACKSON_N = 71
RELATION_HEIGHT = 2 * JACKSON_N - 2
PI_UPPER_NUMERATOR = 355
PI_UPPER_DENOMINATOR = 113
ETA_SIMPLE_CAP = Fraction(297, 50_000)

CORE_FLOORS = {
    0: Fraction(1),
    1: Fraction(35, 41),
    2: Fraction(59, 82),
    3: Fraction(1615, 2706),
    4: Fraction(239, 492),
    5: Fraction(2729, 7380),
    6: Fraction(153101, 568260),
}

HOSTILE_SMALL = (1, 2, 3, 4, 6, 8)
HOSTILE_BODY = (95, 163, 187, 206, 208, 214, 332)
HOSTILE_ROW = HOSTILE_SMALL + HOSTILE_BODY

EXPECTED_BODY_MASS = Fraction(7521335361151, 22863470734060)
EXPECTED_SMALL_MASS = Fraction(20, 41)
EXPECTED_FULL_MASS = Fraction(470973614624713, 3388366362787692)
EXPECTED_EPSILON_SUM = Fraction(
    85546520069055739, 1389230208742953720
)
EXPECTED_LEVEL_ONE_BOUND = -Fraction(
    727156708364069, 33883663627876920
)
EXPECTED_HIGHER_OVERLAP = Fraction(
    1812297618203733, 11294554542625640
)
EXPECTED_WHOLE_COVARIANCE = -Fraction(
    2983319810838331, 138923020874295372
)
EXPECTED_DEPTHS = (
    Fraction(470973614624713, 3388366362787692),
    Fraction(342782540043043, 3080333057079720),
    Fraction(102915338091641, 2258910908525128),
    Fraction(49500703, 7125847970),
    Fraction(1499638563, 103924866973),
    Fraction(666, 634885),
    Fraction(60351927, 5638413685),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def jackson_coefficient(n: int, k: int) -> int:
    """The integer coefficient C_k from THM-2145."""
    require(0 <= k <= 2 * n - 2, "Jackson coefficient index")
    if k <= n:
        return (
            4 * n**3
            - 6 * n * k**2
            + 2 * n
            + 3 * k**3
            - 3 * k
        ) // 6
    return ((2 * n - k) ** 3 - (2 * n - k)) // 6


def jackson_coefficient_by_convolution(n: int, k: int) -> int:
    """Independent convolution of the triangular Fejer coefficients."""
    return sum(
        (n - abs(a)) * (n - abs(k - a))
        for a in range(-(n - 1), n)
        if abs(k - a) < n
    )


def jackson_eta_pi_cap(n: int) -> Fraction:
    """Strict rational eta upper bound obtained from pi < 355/113."""
    c_zero = n * (2 * n * n + 1) // 3
    odd_sum = sum(
        (
            Fraction(jackson_coefficient(n, k), k * k)
            for k in range(1, 2 * n - 2, 2)
        ),
        Fraction(0),
    )
    return Fraction(1, 2) - Fraction(
        4 * PI_UPPER_DENOMINATOR**2,
        PI_UPPER_NUMERATOR**2 * c_zero,
    ) * odd_sum


def circle_fraction(x: Fraction) -> Fraction:
    return x % 1


def is_dangerous(speed: int, time: Fraction) -> bool:
    phase = circle_fraction(speed * time)
    return min(phase, 1 - phase) < RADIUS


def boundary_points(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    points = {Fraction(0), Fraction(1)}
    for speed in speeds:
        for index in range(speed):
            center = Fraction(index, speed)
            half_width = RADIUS / speed
            points.add(circle_fraction(center - half_width))
            points.add(circle_fraction(center + half_width))
    return tuple(sorted(points))


def safe_measure_by_cells(speeds: tuple[int, ...]) -> Fraction:
    points = boundary_points(speeds)
    return sum(
        (
            right - left
            for left, right in zip(points, points[1:])
            if not any(
                is_dangerous(speed, (left + right) / 2)
                for speed in speeds
            )
        ),
        Fraction(0),
    )


def danger_arcs(
    speeds: tuple[int, ...],
) -> list[tuple[Fraction, Fraction]]:
    arcs: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        half_width = RADIUS / speed
        for index in range(speed):
            left = circle_fraction(Fraction(index, speed) - half_width)
            right = left + 2 * half_width
            if right <= 1:
                arcs.append((left, right))
            else:
                arcs.append((left, Fraction(1)))
                arcs.append((Fraction(0), right - 1))
    return arcs


def safe_measure_by_union(speeds: tuple[int, ...]) -> Fraction:
    arcs = sorted(danger_arcs(speeds))
    danger_measure = Fraction(0)
    current_left: Fraction | None = None
    current_right: Fraction | None = None
    for left, right in arcs:
        if current_left is None:
            current_left, current_right = left, right
        elif left <= current_right:
            current_right = max(current_right, right)
        else:
            danger_measure += current_right - current_left
            current_left, current_right = left, right
    if current_left is not None:
        danger_measure += current_right - current_left
    return 1 - danger_measure


def hostile_depth_ledger() -> tuple[
    tuple[Fraction, ...], tuple[Fraction, ...]
]:
    """Mass by small-danger depth on the body-safe set, and per-small masses."""
    points = boundary_points(HOSTILE_ROW)
    depths = [Fraction(0) for _ in range(len(HOSTILE_SMALL) + 1)]
    individual = [Fraction(0) for _ in HOSTILE_SMALL]
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if any(is_dangerous(speed, midpoint) for speed in HOSTILE_BODY):
            continue
        small_state = tuple(
            is_dangerous(speed, midpoint) for speed in HOSTILE_SMALL
        )
        length = right - left
        depths[sum(small_state)] += length
        for index, active in enumerate(small_state):
            if active:
                individual[index] += length
    return tuple(depths), tuple(individual)


def main() -> None:
    require(RELATION_HEIGHT == 140, "Jackson degree changed")

    # Independent closed form versus direct convolution for the full support.
    for k in range(0, RELATION_HEIGHT + 1):
        require(
            jackson_coefficient(JACKSON_N, k)
            == jackson_coefficient_by_convolution(JACKSON_N, k),
            f"Jackson coefficient mismatch at {k}",
        )
        require(
            jackson_coefficient(JACKSON_N, k) > 0,
            f"Jackson multiplier vanished at {k}",
        )

    c_zero = JACKSON_N * (2 * JACKSON_N**2 + 1) // 3
    odd_sum = sum(
        (
            Fraction(jackson_coefficient(JACKSON_N, k), k * k)
            for k in range(1, RELATION_HEIGHT, 2)
        ),
        Fraction(0),
    )
    eta_cap = jackson_eta_pi_cap(JACKSON_N)
    advertised_eta_cap = Fraction(357148519, 60146943550)
    require(c_zero == 238631, "Jackson C_0 changed")
    require(odd_sum > 290903, "odd Jackson sum lost its simple floor")
    require(
        eta_cap < advertised_eta_cap < ETA_SIMPLE_CAP,
        "Jackson eta cap chain failed",
    )

    # N=70 is the hostile adjacent control for this exact global ledger.
    eta_cap_70 = jackson_eta_pi_cap(70)

    plateau_rows = []
    for defect in range(7, 14):
        core_size = 13 - defect
        plateau = CORE_FLOORS[core_size] * SAFE_MASS**defect
        error = (
            13 + core_size * SAFE_MASS**defect
        ) * ETA_SIMPLE_CAP
        margin = plateau - error
        require(margin > 0, f"defect {defect} margin did not close")
        plateau_rows.append((defect, core_size, plateau, margin))

    weakest_plateau = CORE_FLOORS[6] * SAFE_MASS**7
    weakest_margin = weakest_plateau - (
        13 + 6 * SAFE_MASS**7
    ) * ETA_SIMPLE_CAP
    expected_margin = Fraction(
        478970390236831, 39525379884148950000
    )
    require(
        min(row[2] for row in plateau_rows) == weakest_plateau,
        "the weakest lifted plateau moved",
    )
    require(
        min(row[3] for row in plateau_rows) == weakest_margin,
        "the weakest Jackson margin moved",
    )
    require(weakest_margin == expected_margin, "displayed margin changed")
    require(
        weakest_plateau
        - (13 + 6 * SAFE_MASS**7) * eta_cap_70
        < 0,
        "the N=70 adjacent hostile control unexpectedly closes",
    )

    # Exact support: the safe interval transform vanishes precisely at 41Z.
    supported_positive_modes = tuple(
        k for k in range(1, RELATION_HEIGHT + 1) if k % 41
    )
    require(
        len(supported_positive_modes) == 137,
        "height-140 41-unit support census changed",
    )

    # The named failure of all scalar one-peel sign rearrangements.
    depths, danger_masses = hostile_depth_ledger()
    require(depths == EXPECTED_DEPTHS, "hostile depth distribution changed")
    body_mass = sum(depths, Fraction(0))
    full_mass = depths[0]
    small_mass_cells = safe_measure_by_cells(HOSTILE_SMALL)
    body_mass_cells = safe_measure_by_cells(HOSTILE_BODY)
    full_mass_cells = safe_measure_by_cells(HOSTILE_ROW)
    require(body_mass == body_mass_cells == EXPECTED_BODY_MASS, "body mass")
    require(small_mass_cells == EXPECTED_SMALL_MASS, "small mass")
    require(full_mass == full_mass_cells == EXPECTED_FULL_MASS, "full mass")
    require(
        safe_measure_by_union(HOSTILE_SMALL) == EXPECTED_SMALL_MASS,
        "independent small interval union",
    )
    require(
        safe_measure_by_union(HOSTILE_BODY) == EXPECTED_BODY_MASS,
        "independent body interval union",
    )
    require(
        safe_measure_by_union(HOSTILE_ROW) == EXPECTED_FULL_MASS,
        "independent full interval union",
    )

    epsilons = tuple(
        mass - DANGER_MASS * body_mass for mass in danger_masses
    )
    require(all(epsilon > 0 for epsilon in epsilons), "epsilon sign changed")
    epsilon_sum = sum(epsilons, Fraction(0))
    require(epsilon_sum == EXPECTED_EPSILON_SUM, "epsilon sum changed")
    separate_base = (1 - len(HOSTILE_SMALL) * DANGER_MASS) * body_mass
    level_one_bound = separate_base - epsilon_sum
    require(
        level_one_bound == EXPECTED_LEVEL_ONE_BOUND,
        "level-one signed bound changed",
    )

    higher_overlap = sum(
        ((depth - 1) * mass for depth, mass in enumerate(depths) if depth >= 2),
        Fraction(0),
    )
    require(
        higher_overlap == EXPECTED_HIGHER_OVERLAP,
        "higher-overlap correction changed",
    )
    require(
        full_mass == level_one_bound + higher_overlap,
        "depth inclusion-exclusion identity failed",
    )
    whole_covariance = full_mass - body_mass * small_mass_cells
    require(
        whole_covariance == EXPECTED_WHOLE_COVARIANCE,
        "whole-mask covariance changed",
    )

    print("THM-2179 exact reversed-peel relation-packet referee")
    print(
        f"Jackson N={JACKSON_N}, relation_height={RELATION_HEIGHT}, "
        f"C_0={c_zero}"
    )
    print("odd_Jackson_sum>290903")
    print(
        f"eta_pi_cap<{advertised_eta_cap}"
        f"<{ETA_SIMPLE_CAP}"
    )
    print("N=70 fails the same exact global margin ledger")
    print(
        f"weakest_lifted_plateau={weakest_plateau} "
        f"(defect=7, core_size=6)"
    )
    print(f"universal_positive_margin>{weakest_margin}")
    print(
        "nonzero_relation_modes="
        f"{2 * len(supported_positive_modes)} signed 41-units"
    )
    print(f"hostile_small={HOSTILE_SMALL}")
    print(f"hostile_body={HOSTILE_BODY}")
    print(f"hostile_body_mass={body_mass}")
    print(f"hostile_small_mass={small_mass_cells}")
    print(f"hostile_full_mass={full_mass}")
    for speed, epsilon in zip(HOSTILE_SMALL, epsilons):
        print(f"epsilon_{speed}={epsilon}")
    print(f"sum_epsilon=sum_abs_epsilon={epsilon_sum}")
    print(f"separate_charge_base={separate_base}")
    print(f"level_one_signed_bound={level_one_bound}")
    print(f"higher_overlap_correction={higher_overlap}")
    print(f"whole_mask_covariance={whole_covariance}")
    print("independent_exact_interval_evaluators=boundary_cells,interval_union")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
