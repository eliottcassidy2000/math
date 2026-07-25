#!/usr/bin/env python3
"""Exact referee for THM-2166 (hybrid whole-core smoothing).

The six arbitrary speeds are smoothed factorwise with the N=150 Jackson
polynomial (factor height 298).  A retained seven-core E subset {1,...,13}
is first intersected exactly and then smoothed as the single one-dimensional
set G_E with the N=355 Jackson kernel (scalar height 708).

This companion has independent checks for both load-bearing finite objects:

1. all 1,716 core masses and positive-length component counts are computed
   both by iterative rational interval intersection and by a global
   arrangement-cell sweep;
2. Jackson coefficients are computed both from the closed cubic formula and
   by direct convolution of the integer Fejer coefficients.

The support-two conclusion is deliberately *not* attributed to Fourier
sparsity.  Its existence starts from the elementary post-extraction fact
that every seven-subset of {1,...,13} either contains 1 or contains both
members of one of the six pairs (2,3),(4,5),...,(12,13); an independent
set-valued sweep sharpens the coefficient height from 708 to 57 and checks
the exact height-56 hostile profiles.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations


RADIUS = F(1, 14)
ALPHA6 = F(61, 273)
PI2_UPPER = F(355 * 355, 113 * 113)
FAR_N = 150
CORE_N = 355
FAR_ETA_CAP = F(439, 156_250)
CORE_ETA_CAP = F(11_872, 10_000_000)
EXPECTED_BINDING_CORE = (1, 5, 7, 8, 9, 11, 13)
EXPECTED_BINDING_MASS = F(45_107, 229_320)
EXPECTED_BINDING_COMPONENTS = 20
EXPECTED_MARGIN = F(41_050_267, 1_222_741_406_250)
EXPECTED_K_DISTRIBUTION = {
    12: 4,
    14: 31,
    16: 171,
    17: 1,
    18: 452,
    20: 584,
    22: 262,
    24: 152,
    26: 54,
    28: 5,
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def frac(x: F) -> F:
    return x % 1


def is_safe(speed: int, t: F) -> bool:
    x = frac(speed * t)
    return min(x, 1 - x) >= RADIUS


def safe_intervals(speed: int) -> list[tuple[F, F]]:
    """Closed weak-safe intervals in [0,1]."""
    return [
        (F(14 * j + 1, 14 * speed), F(14 * j + 13, 14 * speed))
        for j in range(speed)
    ]


def intersect_intervals(
    left: list[tuple[F, F]], right: list[tuple[F, F]]
) -> list[tuple[F, F]]:
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo <= hi:
            if out and lo <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], hi))
            else:
                out.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return out


def core_by_intersection(core: tuple[int, ...]) -> tuple[F, int, int]:
    intervals = [(F(0), F(1))]
    for speed in core:
        intervals = intersect_intervals(intervals, safe_intervals(speed))
    mass = sum((right - left for left, right in intervals), F(0))
    positive = sum(left < right for left, right in intervals)
    isolated = len(intervals) - positive
    return mass, positive, isolated


def global_boundary_cells() -> tuple[F, ...]:
    points = {F(0), F(1)}
    for speed in range(1, 14):
        for j in range(speed):
            points.add(frac(F(14 * j - 1, 14 * speed)))
            points.add(frac(F(14 * j + 1, 14 * speed)))
    return tuple(sorted(points))


def cores_by_arrangement() -> dict[tuple[int, ...], tuple[F, int]]:
    """Independent mass/component sweep on the 177 global open cells."""
    points = global_boundary_cells()
    cells = tuple(zip(points, points[1:]))
    safe_masks: list[frozenset[int]] = []
    for left, right in cells:
        midpoint = (left + right) / 2
        safe_masks.append(
            frozenset(speed for speed in range(1, 14) if is_safe(speed, midpoint))
        )

    out: dict[tuple[int, ...], tuple[F, int]] = {}
    for core in combinations(range(1, 14), 7):
        flags = [set(core) <= mask for mask in safe_masks]
        mass = sum(
            (right - left for (left, right), safe in zip(cells, flags) if safe),
            F(0),
        )
        # The cell list is cyclic on the circle.  Every nonempty core makes
        # the cells adjacent to 0 unsafe, but use the cyclic formula anyway.
        positive_components = sum(
            safe and not flags[(i - 1) % len(flags)]
            for i, safe in enumerate(flags)
        )
        out[core] = (mass, positive_components)
    return out


def jackson_coefficient_formula(n: int, k: int) -> int:
    require(n >= 2 and 0 <= k <= 2 * n - 2, ("Jackson index", n, k))
    if k <= n:
        numerator = 4 * n**3 - 6 * n * k**2 + 2 * n + 3 * k**3 - 3 * k
    else:
        d = 2 * n - k
        numerator = d**3 - d
    require(numerator % 6 == 0, ("Jackson integrality", n, k))
    value = numerator // 6
    require(value > 0, ("Jackson positivity", n, k))
    return value


def jackson_coefficient_convolution(n: int, k: int) -> int:
    """Independent convolution of A_j=n-|j| for |j|<n."""
    total = 0
    for j in range(-(n - 1), n):
        other = k - j
        if abs(other) < n:
            total += (n - abs(j)) * (n - abs(other))
    return total


def jackson_data(n: int, crosscheck: bool = False) -> tuple[int, F, F]:
    if crosscheck:
        for k in range(2 * n - 1):
            require(
                jackson_coefficient_formula(n, k)
                == jackson_coefficient_convolution(n, k),
                ("Jackson convolution mismatch", n, k),
            )
    c0 = n * (2 * n * n + 1) // 3
    require(c0 == jackson_coefficient_formula(n, 0), ("C0", n))
    odd_sum = sum(
        (
            F(jackson_coefficient_formula(n, k), k * k)
            for k in range(1, 2 * n - 2, 2)
        ),
        F(0),
    )
    # eta_n = 2 int ||x|| J_n(x) dx
    #       = 1/2 - 4/(pi^2 C0) sum_{odd k} C_k/k^2.
    # pi^2 < (355/113)^2 makes the following a strict upper bound.
    eta_upper = F(1, 2) - 4 * odd_sum / (PI2_UPPER * c0)
    require(eta_upper > 0, ("eta positivity", n))
    return c0, odd_sum, eta_upper


def hybrid_margin(mass: F, components: int, far_eta: F, core_eta: F) -> F:
    return (
        (ALPHA6 - 6 * far_eta) * mass
        - 6 * far_eta
        - components * core_eta
    )


def support_two_carrier(core: tuple[int, ...]) -> tuple[int, int | None]:
    """Return 1, or a consecutive pair, giving a {-1,0,1} representation."""
    if 1 in core:
        return 1, None
    for left, right in ((2, 3), (4, 5), (6, 7), (8, 9), (10, 11), (12, 13)):
        if left in core and right in core:
            return left, right
    raise RuntimeError(("seven-core lacks the pigeonhole carrier", core))


def bounded_pair_values(left: int, right: int, height: int) -> frozenset[int]:
    """Independent set-valued implementation of the sparse carry bank."""
    return frozenset(
        left_coefficient * left + right_coefficient * right
        for left_coefficient in range(-height, height + 1)
        for right_coefficient in range(-height, height + 1)
        if abs(left_coefficient * left + right_coefficient * right) <= 708
    )


def sparse_carry_failures(
    cores: tuple[tuple[int, ...], ...], height: int
) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    pair_values = {
        pair: bounded_pair_values(*pair, height)
        for pair in combinations(range(1, 14), 2)
    }
    target = frozenset(range(-708, 709))
    failures = []
    for core in cores:
        represented: set[int] = set()
        for pair in combinations(core, 2):
            represented.update(pair_values[pair])
        if represented != target:
            failures.append((core, tuple(sorted(target - represented))))
    return tuple(failures)


def main() -> None:
    arrangement = cores_by_arrangement()
    require(len(global_boundary_cells()) == 178, "global boundary count")
    require(len(arrangement) == 1716, "seven-core count")

    cores = tuple(combinations(range(1, 14), 7))
    rows: list[tuple[F, tuple[int, ...], F, int, int]] = []
    component_distribution: Counter[int] = Counter()
    one_carriers = pair_carriers = 0
    for core in cores:
        mass, components, isolated = core_by_intersection(core)
        require(
            (mass, components) == arrangement[core],
            ("independent core mismatch", core, mass, components, arrangement[core]),
        )
        carrier = support_two_carrier(core)
        if carrier[1] is None:
            one_carriers += 1
        else:
            pair_carriers += 1
            require(carrier[1] - carrier[0] == 1, ("not consecutive", core, carrier))
        component_distribution[components] += 1
        margin = hybrid_margin(
            mass, components, FAR_ETA_CAP, CORE_ETA_CAP
        )
        rows.append((margin, core, mass, components, isolated))

    require(
        dict(sorted(component_distribution.items())) == EXPECTED_K_DISTRIBUTION,
        ("component distribution", component_distribution),
    )
    require(one_carriers + pair_carriers == 1716, "carrier total")
    require(
        sparse_carry_failures(cores, 57) == (),
        "height-57 support-two carry coverage",
    )
    expected_height_56_failures = (
        ((1, 2, 3, 4, 5, 6, 7), (-706, -705, -699, 699, 705, 706)),
        ((1, 2, 3, 4, 5, 6, 8), (-701, 701)),
        ((1, 2, 4, 5, 6, 8, 10), (-701, 701)),
    )
    require(
        sparse_carry_failures(cores, 56) == expected_height_56_failures,
        "height-56 hostile profiles",
    )

    far_c0, far_sum, far_eta = jackson_data(FAR_N, crosscheck=True)
    core_c0, core_sum, core_eta = jackson_data(CORE_N, crosscheck=True)
    require(far_eta < FAR_ETA_CAP, ("far eta cap", far_eta, FAR_ETA_CAP))
    require(core_eta < CORE_ETA_CAP, ("core eta cap", core_eta, CORE_ETA_CAP))
    require(2 * FAR_N - 2 == 298, "far height")
    require(2 * CORE_N - 2 == 708, "core height")

    rows.sort()
    worst = rows[0]
    require(worst[0] == EXPECTED_MARGIN, ("binding margin", worst))
    require(worst[1] == EXPECTED_BINDING_CORE, ("binding core", worst))
    require(worst[2] == EXPECTED_BINDING_MASS, ("binding mass", worst))
    require(worst[3] == EXPECTED_BINDING_COMPONENTS, ("binding K", worst))
    require(all(row[0] > 0 for row in rows), "nonpositive hybrid margin")
    require(
        sum(row[0] == EXPECTED_MARGIN for row in rows) == 1,
        "hybrid-margin minimizer is not unique",
    )
    adjacent_margin = hybrid_margin(
        EXPECTED_BINDING_MASS,
        EXPECTED_BINDING_COMPONENTS,
        far_eta,
        jackson_data(354, crosscheck=True)[2],
    )
    require(adjacent_margin < 0, ("adjacent N=354 ledger", adjacent_margin))

    # Hostile type check.  The natural seven-torus product does NOT have
    # coefficient-vector support <=2: at vector (1,1,1,0,0,0,0) its
    # coefficient is hat(1_G,1)^3 hat(1_G,0)^4, nonzero because 7 does not
    # divide 1.  The theorem only asserts a post-extraction arithmetic
    # encoding of a scalar |nu|<=708 mode.
    hostile_vector = (1, 1, 1, 0, 0, 0, 0)
    require(sum(value != 0 for value in hostile_vector) == 3, "hostile support")
    require(1 % 7 != 0, "hostile one-mode should be nonzero")

    print("THM-2166 HYBRID WHOLE-CORE SMOOTHING -- EXACT REFEREE")
    print(
        "universe: all C(13,7)=1716 retained cores; "
        "178 boundary points / 177 open cells"
    )
    print(
        "independent core paths: iterative closed-interval intersection "
        "and global arrangement-cell runs -- exact agreement"
    )
    print(f"positive-component_distribution={dict(sorted(component_distribution.items()))}")
    print(
        f"support carriers: contains_1={one_carriers}, "
        f"consecutive_pair={pair_carriers}"
    )
    print(
        "independent set-valued sparse bank: height 57 covers all "
        "1417 carries for every core; height 56 has exactly the 3 "
        "recorded hostile cores"
    )
    print()
    print(
        f"far Jackson: N={FAR_N}, H={2*FAR_N-2}, C0={far_c0}, "
        f"floor(odd_sum)={far_sum.numerator // far_sum.denominator}, "
        f"eta_upper<{FAR_ETA_CAP}"
    )
    print(
        f"core Jackson: N={CORE_N}, H={2*CORE_N-2}, C0={core_c0}, "
        f"floor(odd_sum)={core_sum.numerator // core_sum.denominator}, "
        f"eta_upper<{CORE_ETA_CAP}"
    )
    print(
        "independent Jackson paths: cubic coefficient formula and direct "
        "integer Fejer convolution -- exact agreement"
    )
    print()
    print(
        f"binding_core={worst[1]} mass={worst[2]} "
        f"positive_components={worst[3]} isolated={worst[4]}"
    )
    print(
        "minimum certified margin="
        f"{worst[0]} (~{float(worst[0]):.12g}) > 0"
    )
    print(
        "consequence: common scalar frequency 0<|nu|<=708; "
        "far coefficient height<=298; core arithmetic support<=2 and "
        "coefficient height<=57"
    )
    print(
        "HOSTILE TYPE CHECK: support<=2 is NOT Fourier-vector sparsity; "
        "the natural torus product has nonzero support-3 coefficient at "
        f"{hostile_vector}"
    )
    print(
        "full relation l1<=1902; the same exact rational-pi ledger is "
        "negative at core N=354"
    )
    print("ALL EXACT ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
