#!/usr/bin/env python3
"""Exact referee for THM-2179.

The all-row argument in THM-2179 is Fourier/Jackson analysis.  This companion
checks every rational constant in the defect-at-least-seven and defect-six
packets, exhausts the seven-core floor by two exact interval algorithms,
verifies the adjacent Jackson cutoffs, and replays the named reversed-peel
hostile row.  Assertions remain active under ``python -O``.
"""

from fractions import Fraction
from itertools import combinations


RADIUS = Fraction(3, 41)
DANGER_MASS = 2 * RADIUS
SAFE_MASS = 1 - DANGER_MASS
D7_JACKSON_N = 71
D7_RELATION_HEIGHT = 2 * D7_JACKSON_N - 2
D6_JACKSON_N = 91
D6_RELATION_HEIGHT = 2 * D6_JACKSON_N - 2
PI_UPPER_NUMERATOR = 355
PI_UPPER_DENOMINATOR = 113
D7_ETA_SIMPLE_CAP = Fraction(297, 50_000)
D6_ETA_SIMPLE_CAP = Fraction(93, 20_000)

CORE_FLOORS = {
    0: Fraction(1),
    1: Fraction(35, 41),
    2: Fraction(59, 82),
    3: Fraction(1615, 2706),
    4: Fraction(239, 492),
    5: Fraction(2729, 7380),
    6: Fraction(153101, 568260),
    7: Fraction(39965, 211068),
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


def seven_core_masses_by_global_cells() -> tuple[
    dict[tuple[int, ...], Fraction], tuple[Fraction, ...]
]:
    """Evaluate all seven-cores simultaneously on one global arrangement."""
    points = boundary_points(tuple(range(1, 14)))
    masses = {
        core: Fraction(0) for core in combinations(range(1, 14), 7)
    }
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        safe_speeds = tuple(
            speed
            for speed in range(1, 14)
            if not is_dangerous(speed, midpoint)
        )
        width = right - left
        for core in combinations(safe_speeds, 7):
            masses[core] += width
    return masses, points


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
    require(D7_RELATION_HEIGHT == 140, "defect-seven Jackson degree changed")
    require(D6_RELATION_HEIGHT == 180, "defect-six Jackson degree changed")

    # Independent closed form versus direct convolution for both supports.
    for jackson_n, relation_height in (
        (D7_JACKSON_N, D7_RELATION_HEIGHT),
        (D6_JACKSON_N, D6_RELATION_HEIGHT),
    ):
        for k in range(0, relation_height + 1):
            require(
                jackson_coefficient(jackson_n, k)
                == jackson_coefficient_by_convolution(jackson_n, k),
                f"Jackson coefficient mismatch at N={jackson_n}, k={k}",
            )
            require(
                jackson_coefficient(jackson_n, k) > 0,
                f"Jackson multiplier vanished at N={jackson_n}, k={k}",
            )

    c_zero_71 = D7_JACKSON_N * (2 * D7_JACKSON_N**2 + 1) // 3
    odd_sum_71 = sum(
        (
            Fraction(jackson_coefficient(D7_JACKSON_N, k), k * k)
            for k in range(1, D7_RELATION_HEIGHT, 2)
        ),
        Fraction(0),
    )
    eta_cap_71 = jackson_eta_pi_cap(D7_JACKSON_N)
    advertised_eta_cap_71 = Fraction(357148519, 60146943550)
    require(c_zero_71 == 238631, "N=71 Jackson C_0 changed")
    require(
        odd_sum_71 > 290903,
        "N=71 odd Jackson sum lost its simple floor",
    )
    require(
        eta_cap_71 < advertised_eta_cap_71 < D7_ETA_SIMPLE_CAP,
        "N=71 Jackson eta cap chain failed",
    )

    # N=70 is the hostile adjacent control for this exact global ledger.
    eta_cap_70 = jackson_eta_pi_cap(70)

    plateau_rows = []
    for defect in range(7, 14):
        core_size = 13 - defect
        plateau = CORE_FLOORS[core_size] * SAFE_MASS**defect
        error = (
            13 + core_size * SAFE_MASS**defect
        ) * D7_ETA_SIMPLE_CAP
        margin = plateau - error
        require(margin > 0, f"defect {defect} margin did not close")
        plateau_rows.append((defect, core_size, plateau, margin))

    weakest_plateau = CORE_FLOORS[6] * SAFE_MASS**7
    weakest_margin = weakest_plateau - (
        13 + 6 * SAFE_MASS**7
    ) * D7_ETA_SIMPLE_CAP
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

    # Exhaust the 1,716 seven-cores by independent exact interval algorithms.
    seven_core_masses, seven_core_points = (
        seven_core_masses_by_global_cells()
    )
    require(
        len(seven_core_points) == 184,
        "seven-core global boundary count changed",
    )
    seven_core_minimum: Fraction | None = None
    seven_core_minimizers: list[tuple[int, ...]] = []
    seven_core_count = 0
    for core in combinations(range(1, 14), 7):
        core_tuple = tuple(core)
        mass_by_cells = seven_core_masses[core_tuple]
        mass_by_union = safe_measure_by_union(core_tuple)
        require(
            mass_by_cells == mass_by_union,
            f"seven-core evaluators disagree at {core_tuple}",
        )
        seven_core_count += 1
        if seven_core_minimum is None or mass_by_cells < seven_core_minimum:
            seven_core_minimum = mass_by_cells
            seven_core_minimizers = [core_tuple]
        elif mass_by_cells == seven_core_minimum:
            seven_core_minimizers.append(core_tuple)
    expected_seven_core_minimizer = (1, 5, 7, 8, 9, 11, 13)
    require(seven_core_count == 1716, "seven-core census size changed")
    require(
        seven_core_minimum == CORE_FLOORS[7],
        "seven-core floor changed",
    )
    require(
        seven_core_minimizers == [expected_seven_core_minimizer],
        "seven-core minimizer is no longer unique",
    )

    # The seven-core floor closes the d=6 ledger at N=91.
    c_zero_91 = D6_JACKSON_N * (2 * D6_JACKSON_N**2 + 1) // 3
    odd_sum_91 = sum(
        (
            Fraction(jackson_coefficient(D6_JACKSON_N, k), k * k)
            for k in range(1, D6_RELATION_HEIGHT, 2)
        ),
        Fraction(0),
    )
    eta_cap_91 = jackson_eta_pi_cap(D6_JACKSON_N)
    advertised_eta_cap_91 = Fraction(586539659, 126632692550)
    require(c_zero_91 == 502411, "N=91 Jackson C_0 changed")
    require(
        odd_sum_91 > 614083,
        "N=91 odd Jackson sum lost its simple floor",
    )
    require(
        eta_cap_91 < advertised_eta_cap_91 < D6_ETA_SIMPLE_CAP,
        "N=91 Jackson eta cap chain failed",
    )

    defect_six_plateau = CORE_FLOORS[7] * SAFE_MASS**6
    defect_six_error_coefficient = 13 + 7 * SAFE_MASS**6
    defect_six_margin = defect_six_plateau - (
        defect_six_error_coefficient * D6_ETA_SIMPLE_CAP
    )
    expected_defect_six_plateau = Fraction(
        73466285703125, 1002595001939388
    )
    expected_defect_six_error_coefficient = Fraction(
        74619214508, 4750104241
    )
    expected_defect_six_margin = Fraction(
        287560991216713, 1253243752424235000
    )
    require(
        defect_six_plateau == expected_defect_six_plateau,
        "defect-six plateau changed",
    )
    require(
        defect_six_error_coefficient
        == expected_defect_six_error_coefficient,
        "defect-six error coefficient changed",
    )
    require(
        defect_six_margin == expected_defect_six_margin > 0,
        "defect-six margin failed",
    )
    eta_cap_90 = jackson_eta_pi_cap(90)
    require(
        defect_six_plateau
        - defect_six_error_coefficient * eta_cap_90
        < 0,
        "the N=90 defect-six adjacent control unexpectedly closes",
    )

    # Exact support: the safe interval transform vanishes precisely at 41Z.
    supported_positive_modes_140 = tuple(
        k for k in range(1, D7_RELATION_HEIGHT + 1) if k % 41
    )
    require(
        len(supported_positive_modes_140) == 137,
        "height-140 41-unit support census changed",
    )
    supported_positive_modes_180 = tuple(
        k for k in range(1, D6_RELATION_HEIGHT + 1) if k % 41
    )
    require(
        len(supported_positive_modes_180) == 176,
        "height-180 41-unit support census changed",
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
        f"defect>=7: Jackson N={D7_JACKSON_N}, "
        f"relation_height={D7_RELATION_HEIGHT}, C_0={c_zero_71}"
    )
    print("defect>=7: odd_Jackson_sum>290903")
    print(
        f"defect>=7: eta_pi_cap<{advertised_eta_cap_71}"
        f"<{D7_ETA_SIMPLE_CAP}"
    )
    print("defect>=7: N=70 fails the same exact global margin ledger")
    print(
        f"weakest_lifted_plateau={weakest_plateau} "
        f"(defect=7, core_size=6)"
    )
    print(f"universal_positive_margin>{weakest_margin}")
    print(
        "height_140_nonzero_relation_modes="
        f"{2 * len(supported_positive_modes_140)} signed 41-units"
    )
    print(
        f"seven_core_census={seven_core_count}, "
        f"floor={seven_core_minimum}, "
        f"unique_minimizer={seven_core_minimizers[0]}"
    )
    print(
        f"defect=6: Jackson N={D6_JACKSON_N}, "
        f"relation_height={D6_RELATION_HEIGHT}, C_0={c_zero_91}"
    )
    print("defect=6: odd_Jackson_sum>614083")
    print(
        f"defect=6: eta_pi_cap<{advertised_eta_cap_91}"
        f"<{D6_ETA_SIMPLE_CAP}"
    )
    print(f"defect_six_lifted_plateau={defect_six_plateau}")
    print(
        "defect_six_error_coefficient="
        f"{defect_six_error_coefficient}"
    )
    print(f"defect_six_positive_margin>{defect_six_margin}")
    print("defect=6: N=90 fails the same exact global margin ledger")
    print(
        "height_180_nonzero_relation_modes="
        f"{2 * len(supported_positive_modes_180)} signed 41-units"
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
