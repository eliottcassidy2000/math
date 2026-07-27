#!/usr/bin/env python3
"""Exact referee for THM-2464.

The companion uses only integer and Fraction arithmetic.  Every
truth-bearing executable check uses ``require`` so optimized Python
verifies the same claims.
"""

from fractions import Fraction as F


P = 13
ORDINARY_RADIUS = F(1, 14)
GUARD_RADIUS = F(1, 7)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_fraction(x):
    return x.numerator // x.denominator


def nearest_integer(x):
    return floor_fraction(x + F(1, 2))


def fractional_part(x):
    return x - floor_fraction(x)


def circle_distance(x):
    return abs(x - nearest_integer(x))


def danger(speed, x, radius=ORDINARY_RADIUS):
    return circle_distance(F(speed) * x) < radius


def phase_root_mask(rho, phase, radius):
    return frozenset(
        r
        for r in range(P)
        if circle_distance((phase + rho * r) / P) < radius
    )


def physical_root_mask(speed, y, radius):
    return frozenset(
        r
        for r in range(P)
        if circle_distance(F(speed) * (y + r) / P) < radius
    )


def interval_width(interval):
    return interval[1] - interval[0]


def interval_midpoint(interval):
    return (interval[0] + interval[1]) / 2


def fit_full_inverse_branch(source, target, multiplier):
    """Return one full pullback of target compactly inside source."""
    source_lo, source_hi = source
    target_lo, target_hi = target
    lower = floor_fraction(multiplier * source_lo - target_lo) - 2
    upper = floor_fraction(multiplier * source_hi - target_hi) + 3
    for k in range(lower, upper + 1):
        pullback = (
            (target_lo + k) / multiplier,
            (target_hi + k) / multiplier,
        )
        if source_lo < pullback[0] and pullback[1] < source_hi:
            return pullback, k
    return None


def equal_centered_subinterval(interval, target_width):
    require(interval_width(interval) >= target_width, "equal-width shrink")
    middle = interval_midpoint(interval)
    result = (middle - target_width / 2, middle + target_width / 2)
    require(
        interval[0] <= result[0] and result[1] <= interval[1],
        "centered subinterval containment",
    )
    return result


# ---------------------------------------------------------------------------
# 1. The literal shallow-depth multiplicity law, including endpoints.
# ---------------------------------------------------------------------------

ordinary_multiplicity_checks = 0
ordinary_endpoint_exceptions = 0
guard_multiplicity_checks = 0

ordinary_cases = (
    (F(0), 1, True),
    (F(-1, 14), 1, False),
    (F(1, 14), 1, False),
    (F(-1, 4), 2, False),
    (F(1, 4), 2, False),
)

guard_cases = (
    (F(0), 3, True),
    (F(-1, 28), 3, True),
    (F(1, 28), 3, True),
    (F(-1, 14), 3, False),
    (F(1, 14), 3, False),
    (F(-1, 7), 3, False),
    (F(1, 7), 3, False),
    (F(-1, 4), 4, False),
    (F(1, 4), 4, False),
)

for rho in range(1, P):
    for start in range(P):
        base_nearest = -rho * start
        for delta, expected_count, expected_shallow in ordinary_cases:
            phase = F(base_nearest) + delta
            root_count = len(
                phase_root_mask(rho, phase, ORDINARY_RADIUS)
            )
            shallow = circle_distance(phase) < ORDINARY_RADIUS
            require(
                root_count == expected_count,
                f"ordinary multiplicity rho={rho}, start={start}, delta={delta}",
            )
            require(
                shallow == expected_shallow,
                f"ordinary shallow bit rho={rho}, start={start}, delta={delta}",
            )
            if abs(delta) == ORDINARY_RADIUS:
                require(
                    root_count == 1 and not shallow,
                    "ordinary endpoint exception",
                )
                ordinary_endpoint_exceptions += 1
            ordinary_multiplicity_checks += 1

        guard_nearest = -rho * start - 1
        for delta, expected_count, expected_shallow in guard_cases:
            phase = F(guard_nearest) + delta
            root_count = len(phase_root_mask(rho, phase, GUARD_RADIUS))
            shallow = circle_distance(phase) < ORDINARY_RADIUS
            require(
                root_count == expected_count,
                f"guard multiplicity rho={rho}, start={start}, delta={delta}",
            )
            require(
                shallow == expected_shallow,
                f"guard shallow bit rho={rho}, start={start}, delta={delta}",
            )
            if root_count == 4:
                require(not shallow, "four-root guard forces shallow safe")
            guard_multiplicity_checks += 1


# ---------------------------------------------------------------------------
# 2. The sharp phase-before-speed lambda=2 threshold.
# ---------------------------------------------------------------------------

DANGER_TARGET = (F(-1, 14), F(1, 14))
SAFE_TARGET = (F(1, 7), F(3, 14))


def ordinary_phase_intervals(rho, start):
    negative_nearest = -rho * start
    positive_nearest = -rho * start - 1
    return (
        (
            F(negative_nearest) - F(1, 2),
            F(negative_nearest) - F(1, 14),
        ),
        (
            F(positive_nearest) + F(1, 14),
            F(positive_nearest) + F(1, 2),
        ),
    )


def guard_phase_intervals(rho, start):
    negative_nearest = -rho * start - 1
    positive_nearest = -rho * start - 2
    return (
        (
            F(negative_nearest) - F(1, 2),
            F(negative_nearest) - F(1, 7),
        ),
        (
            F(positive_nearest) + F(1, 7),
            F(positive_nearest) + F(1, 2),
        ),
    )


lambda_two_phase_checks = 0
lambda_two_widths = set()

require(
    P * F(3, 7) - (1 + interval_width(DANGER_TARGET)) == F(31, 7),
    "ordinary lambda-two mesh margin",
)
require(
    P * F(5, 14) - (1 + interval_width(DANGER_TARGET)) == F(7, 2),
    "guard lambda-two mesh margin",
)

for rho in range(1, P):
    for start in range(P):
        inverse = pow(rho, -1, P)
        expected_ordinary = frozenset((start, (start + inverse) % P))
        expected_guard = frozenset(
            (start + j * inverse) % P for j in range(4)
        )
        for source in ordinary_phase_intervals(rho, start):
            require(interval_width(source) == F(3, 7), "ordinary phase width")
            for desired, target in ((True, DANGER_TARGET), (False, SAFE_TARGET)):
                fitted = fit_full_inverse_branch(source, target, P)
                require(fitted is not None, "ordinary lambda-two branch")
                refined, _ = fitted
                test = interval_midpoint(refined)
                require(
                    danger(P, test) == desired,
                    "ordinary linked blocker truth",
                )
                require(
                    phase_root_mask(rho, test, ORDINARY_RADIUS)
                    == expected_ordinary,
                    "ordinary old root mask retained",
                )
                lambda_two_widths.add(interval_width(refined))
                lambda_two_phase_checks += 1

        for source in guard_phase_intervals(rho, start):
            require(interval_width(source) == F(5, 14), "guard phase width")
            for desired, target in ((True, DANGER_TARGET), (False, SAFE_TARGET)):
                fitted = fit_full_inverse_branch(source, target, P)
                require(fitted is not None, "guard lambda-two branch")
                refined, _ = fitted
                test = interval_midpoint(refined)
                require(
                    danger(P, test) == desired,
                    "guard linked blocker truth",
                )
                require(
                    phase_root_mask(rho, test, GUARD_RADIUS)
                    == expected_guard,
                    "guard old root mask retained",
                )
                lambda_two_widths.add(interval_width(refined))
                lambda_two_phase_checks += 1

require(
    lambda_two_widths == {F(1, 91), F(1, 182)},
    "lambda-two refined widths",
)
require(lambda_two_phase_checks == 1248, "lambda-two phase census")


# ---------------------------------------------------------------------------
# 3. Exact fixed-tuple joint cells for the THM-2462 control.
# ---------------------------------------------------------------------------

ROWS = (
    (0, 9, 7, 3, 12),
    (1, 2, 9, 5, 7),
    (2, 3, 10, 6, 3),
    (3, 2, 10, 6, 7),
    (4, 10, 7, 2, 11),
    (5, 6, 3, 9, 11),
    (6, 2, 9, 2, 7),
    (7, 8, 2, 2, 6),
    (8, 2, 6, 9, 11),
    (9, 8, 3, 2, 7),
    (10, 9, 4, 3, 3),
    (11, 2, 6, 9, 5),
    (12, 8, 6, 2, 11),
)

WORD_BITS = (0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0)

CENTERS = (
    F(3034865843, 3332313600),
    F(9058389143, 9996940800),
    F(17899453493, 19993881600),
    F(18124187051, 19993881600),
    F(8996707229, 9996940800),
    F(501517651, 555385600),
    F(362322469, 399877632),
    F(3604728023, 3998776320),
    F(823790063, 908812800),
    F(360620951, 399877632),
    F(18217325957, 19993881600),
    F(9062455079, 9996940800),
    F(3003819443, 3332313600),
)

PARENT_WIDTH = F(1, 6664627200)
PARENTS = tuple(
    (center - PARENT_WIDTH / 2, center + PARENT_WIDTH / 2)
    for center in CENTERS
)

SPEEDS = (14, 644, 27048, 1190112, 57125376, 2856268800)
RADII = (
    ORDINARY_RADIUS,
    ORDINARY_RADIUS,
    GUARD_RADIUS,
    ORDINARY_RADIUS,
    ORDINARY_RADIUS,
    ORDINARY_RADIUS,
)

U_K = SPEEDS[0]
U_Q = SPEEDS[1]
require(U_Q == 46 * U_K, "linked unit ratio")

# Keys are (K-linked danger bit, q-linked danger bit).
JOINT_COMPONENTS = {
    (True, True): (F(-1, 9016), F(1, 9016)),
    (True, False): (F(1, 9016), F(13, 9016)),
    (False, True): (F(321, 9016), F(323, 9016)),
    (False, False): (F(309, 9016), F(321, 9016)),
}

expected_joint_widths = {
    (True, True): F(1, 4508),
    (True, False): F(3, 2254),
    (False, True): F(1, 4508),
    (False, False): F(3, 2254),
}

joint_component_checks = 0
for pattern, component in JOINT_COMPONENTS.items():
    require(
        interval_width(component) == expected_joint_widths[pattern],
        f"joint component width {pattern}",
    )
    for numerator in (1, 2, 3):
        z = component[0] + interval_width(component) * F(numerator, 4)
        require(
            (danger(U_K, z), danger(U_Q, z)) == pattern,
            f"joint component truth {pattern}",
        )
        joint_component_checks += 1

COMMON_N = 9
COMMON_LAMBDA = COMMON_N + 1
COMMON_MULTIPLIER = P**COMMON_N

require(
    COMMON_MULTIPLIER * PARENT_WIDTH
    - (1 + F(3, 2254))
    == F(3931001773, 6664627200),
    "large joint-cell mesh margin",
)
require(
    COMMON_MULTIPLIER * PARENT_WIDTH
    - (1 + F(1, 4508))
    == F(3938393773, 6664627200),
    "small joint-cell mesh margin",
)

linked_intervals = []
linked_chart_checks = 0
for row_index, (row, parent) in enumerate(zip(ROWS, PARENTS)):
    g, q, a, b, c = row
    desired = (g % 2 == 0, bool(WORD_BITS[row_index]))
    fitted = fit_full_inverse_branch(
        parent, JOINT_COMPONENTS[desired], COMMON_MULTIPLIER
    )
    require(fitted is not None, f"linked joint cell row {row_index}")
    interval, _ = fitted
    linked_intervals.append(interval)

    y = interval_midpoint(interval)
    require(
        (danger(COMMON_MULTIPLIER * U_K, y),
         danger(COMMON_MULTIPLIER * U_Q, y))
        == desired,
        f"linked midpoint bits row {row_index}",
    )

    expected_masks = (
        frozenset((0, 1)),
        frozenset((q, (q + 2) % P)),
        frozenset((g, (g + 5) % P, (g + 10) % P, (g + 15) % P)),
        frozenset((a, (a + 1) % P)),
        frozenset((b, (b + 3) % P)),
        frozenset((c, (c + 5) % P)),
    )
    require(
        tuple(
            physical_root_mask(speed, y, radius)
            for speed, radius in zip(SPEEDS, RADII)
        )
        == expected_masks,
        f"old root word row {row_index}",
    )

    physical_c_k = P**COMMON_LAMBDA * U_K
    physical_c_q = P**COMMON_LAMBDA * U_Q
    for desired_bit, physical_speed in zip(desired, (physical_c_k, physical_c_q)):
        root_values = {
            danger(physical_speed, (y + r) / P)
            for r in range(P)
        }
        require(
            root_values == {desired_bit},
            f"root-constant linked bit row {row_index}",
        )
        linked_chart_checks += 1

require(
    {interval_width(interval) for interval in linked_intervals}
    == {F(1, 47805083173484), F(3, 23902541586742)},
    "linked interval widths",
)

LINKED_EQUAL_WIDTH = F(1, 47805083173484)
linked_equal_intervals = tuple(
    equal_centered_subinterval(interval, LINKED_EQUAL_WIDTH)
    for interval in linked_intervals
)
require(
    all(
        left[1] < right[0] or right[1] < left[0]
        for i, left in enumerate(linked_equal_intervals)
        for right in linked_equal_intervals[i + 1 :]
    ),
    "linked intervals remain disjoint",
)

PHYSICAL_C_K = P**COMMON_LAMBDA * U_K
PHYSICAL_C_Q = P**COMMON_LAMBDA * U_Q
require(PHYSICAL_C_K == 1930018885886, "physical K-linked speed")
require(PHYSICAL_C_Q == 88780868750756, "physical q-linked speed")


# ---------------------------------------------------------------------------
# 4. A fixed positive delayed word after the linked blockers.
# ---------------------------------------------------------------------------

DELAYED_N = 13
DELAYED_K = DELAYED_N + 1
DELAYED_MULTIPLIER = P**DELAYED_N
Q_TARGET = (F(2, 7), F(3, 7))

require(
    DELAYED_MULTIPLIER * LINKED_EQUAL_WIDTH == F(28561, 4508),
    "delayed-word mesh ratio",
)
require(
    DELAYED_MULTIPLIER * LINKED_EQUAL_WIDTH
    > 1 + interval_width(Q_TARGET),
    "delayed-word universal fit",
)

delayed_word_checks = 0
delayed_intervals = []
for row_index, source in enumerate(linked_equal_intervals):
    fitted = fit_full_inverse_branch(source, Q_TARGET, DELAYED_MULTIPLIER)
    require(fitted is not None, f"delayed word row {row_index}")
    interval, _ = fitted
    delayed_intervals.append(interval)
    y = interval_midpoint(interval)
    q_value = fractional_part(DELAYED_MULTIPLIER * y)
    require(
        Q_TARGET[0] < q_value < Q_TARGET[1],
        f"delayed word truth row {row_index}",
    )
    desired = (ROWS[row_index][0] % 2 == 0, bool(WORD_BITS[row_index]))
    require(
        (danger(COMMON_MULTIPLIER * U_K, y),
         danger(COMMON_MULTIPLIER * U_Q, y))
        == desired,
        f"delayed shrink keeps linked bits row {row_index}",
    )
    delayed_word_checks += 1

require(
    {interval_width(interval) for interval in delayed_intervals}
    == {F(1, 2120125746145771)},
    "delayed-word interval width",
)


# ---------------------------------------------------------------------------
# 5. Exact hostile boundaries.
# ---------------------------------------------------------------------------

SHALLOW_HOSTILE_PARENT = (F(0), F(1, 28))
for numerator in (1, 2, 3):
    y = (
        SHALLOW_HOSTILE_PARENT[0]
        + interval_width(SHALLOW_HOSTILE_PARENT) * F(numerator, 4)
    )
    require(danger(1, y), "fixed shallow clock is wholly dangerous")

require(
    fit_full_inverse_branch(SHALLOW_HOSTILE_PARENT, SAFE_TARGET, 1) is None,
    "fixed shallow safe request has no branch",
)

require(128 // 2 == 64, "conditional complete-bank halving")
require(64**2 == 4096, "conditional service denominator")
require((2 * (64 - 2) - 1) ** 2 == 15129, "conditional drift denominator")
require(4096 * 2028 == 8306688, "conditional root denominator")
require(384 // 2 == 192, "conditional three-word bank halving")

for numerator in (1, 2, 3):
    z = F(numerator, 17)
    require(
        not (danger(1, z) and not danger(1, z)),
        "same-factor opposite-bit joint cell is empty",
    )


print("THM2464 LINKED-BLOCKER CLOCK-CELL AUDIT")
print(f"ordinary_multiplicity_checks={ordinary_multiplicity_checks}")
print(f"ordinary_endpoint_exceptions={ordinary_endpoint_exceptions}")
print(f"guard_multiplicity_checks={guard_multiplicity_checks}")
print(f"lambda_two_phase_checks={lambda_two_phase_checks}")
print(f"lambda_two_min_phase_width={min(lambda_two_widths)}")
print(f"joint_component_checks={joint_component_checks}")
print(
    "joint_component_widths_SS_SD_DS_DD="
    + ",".join(
        str(expected_joint_widths[key])
        for key in ((False, False), (False, True), (True, False), (True, True))
    )
)
print(f"common_clock_n={COMMON_N} lambda={COMMON_LAMBDA}")
print(f"linked_speeds_cK_cq={PHYSICAL_C_K},{PHYSICAL_C_Q}")
print(f"linked_chart_checks={linked_chart_checks}")
print(f"linked_equal_width={LINKED_EQUAL_WIDTH}")
print(f"delayed_clock_n={DELAYED_N} k={DELAYED_K}")
print(
    "delayed_mesh_ratio="
    f"{DELAYED_MULTIPLIER * LINKED_EQUAL_WIDTH}"
)
print(f"delayed_word_checks={delayed_word_checks}")
print(f"delayed_word_width={interval_width(delayed_intervals[0])}")
print("shallow_fixed_clock_hostile=PASS")
print("same_factor_joint_cell_hostile=PASS")
print("conditional_bank_128_to_64=4096,15129,8306688,192")
print("optimized_require_path=PASS")
