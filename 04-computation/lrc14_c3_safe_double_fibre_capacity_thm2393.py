#!/usr/bin/env python3
"""Exact companion for THM-2393.

The script reconstructs the ten compatible low-cage unions, their
seven-root empty-fibre laws, and the thirteen-root refinement which closes
the oriented (2,1) boundary.  All arithmetic is exact.
"""

from collections import defaultdict
from fractions import Fraction as F
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def circle_norm(x: F) -> F:
    r = frac_part(x)
    return min(r, 1 - r)


def danger(speed: int, x: F) -> bool:
    return circle_norm(speed * x) < F(1, 14)


def boundary_grid(speeds: tuple[int, ...], root_scale: int = 1) -> tuple[F, ...]:
    """Boundaries in the base variable for x=(base+r)/root_scale."""

    points = {F(0), F(1)}
    for speed in speeds:
        for tooth in range(speed):
            for sign in (-1, 1):
                endpoint = F(root_scale * (14 * tooth + sign), 14 * speed)
                points.add(frac_part(endpoint))
    return tuple(sorted(points))


def exact_histogram(
    boundaries: tuple[F, ...],
    evaluator,
) -> dict[int, F]:
    out: dict[int, F] = defaultdict(F)
    for left, right in zip(boundaries, boundaries[1:]):
        if left == right:
            continue
        midpoint = (left + right) / 2
        out[evaluator(midpoint)] += right - left
    require(sum(out.values(), F()) == 1, "cell histogram lost Haar mass")
    return dict(sorted(out.items()))


def low_cage(a: int, b: int, x: F) -> bool:
    return (
        danger(b, x) or danger(13 * a, x)
    ) and (
        danger(13 * b, x) or danger(169 * a, x)
    )


def low_cage_speeds(a: int, b: int) -> tuple[int, ...]:
    return b, 13 * a, 13 * b, 169 * a


def low_cage_measure(a: int, b: int) -> F:
    grid = boundary_grid(low_cage_speeds(a, b))
    hist = exact_histogram(grid, lambda x: int(low_cage(a, b, x)))
    return hist.get(1, F())


def seven_count(a: int, b: int, y: F) -> int:
    return sum(low_cage(a, b, (y + root) / 7) for root in range(7))


def seven_histogram(a: int, b: int) -> dict[int, F]:
    grid = boundary_grid(low_cage_speeds(a, b), root_scale=7)
    hist = exact_histogram(grid, lambda y: seven_count(a, b, y))
    mean = sum(F(count) * mass for count, mass in hist.items())
    require(
        mean == 7 * low_cage_measure(a, b),
        "seven-root disintegration mean failed",
    )
    return hist


def thirteen_empty_histogram(a: int, b: int) -> dict[int, F]:
    grid = boundary_grid(low_cage_speeds(a, b), root_scale=7 * 13)

    def occupancy(z: F) -> int:
        return sum(
            seven_count(a, b, (z + root) / 13) == 0
            for root in range(13)
        )

    return exact_histogram(grid, occupancy)


def pair_intersection_measure(a: int, b: int) -> F:
    grid = boundary_grid((a, b))
    hist = exact_histogram(
        grid,
        lambda x: int(danger(a, x) and danger(b, x)),
    )
    return hist.get(1, F())


TYPES = (
    (1, 1),
    (1, 2),
    (2, 1),
    (1, 3),
    (3, 1),
    (4, 1),
    (2, 3),
    (3, 2),
    (3, 4),
    (4, 3),
)

EXPECTED_MASS = {
    (1, 1): F(193, 1183),
    (1, 2): F(114, 1183),
    (2, 1): F(239, 2366),
    (1, 3): F(263, 3549),
    (3, 1): F(95, 1183),
    (4, 1): F(331, 4732),
    (2, 3): F(43, 546),
    (3, 2): F(95, 1183),
    (3, 4): F(491, 7098),
    (4, 3): F(331, 4732),
}

EXPECTED_SEVEN_HIST = {
    (1, 1): {1: F(145, 169), 2: F(24, 169)},
    (1, 2): {0: F(66, 169), 1: F(92, 169), 2: F(11, 169)},
    (2, 1): {0: F(62, 169), 1: F(189, 338), 2: F(25, 338)},
    (1, 3): {0: F(88, 169), 1: F(223, 507), 2: F(20, 507)},
    (3, 1): {0: F(248, 507), 1: F(233, 507), 2: F(2, 39)},
    (4, 1): {0: F(93, 169), 1: F(277, 676), 2: F(27, 676)},
    (2, 3): {0: F(84, 169), 1: F(461, 1014), 2: F(49, 1014)},
    (3, 2): {0: F(248, 507), 1: F(233, 507), 2: F(2, 39)},
    (3, 4): {0: F(281, 507), 1: F(413, 1014), 2: F(1, 26)},
    (4, 3): {0: F(93, 169), 1: F(277, 676), 2: F(27, 676)},
}


# ---------------------------------------------------------------------------
# 1. Whole compatible-union table and seven-root empty-fibre table.
# ---------------------------------------------------------------------------

observed_mass = {pair: low_cage_measure(*pair) for pair in TYPES}
require(observed_mass == EXPECTED_MASS, "compatible-union mass table")

observed_seven = {pair: seven_histogram(*pair) for pair in TYPES}
require(observed_seven == EXPECTED_SEVEN_HIST, "seven-root histogram table")


# ---------------------------------------------------------------------------
# 2. Exact high-safe mass.
# ---------------------------------------------------------------------------

same_line_overlap = pair_intersection_measure(1, 13)
require(same_line_overlap == F(1, 91), "same-line 1/91 overlap")
high_pair_safe = 1 - F(2, 7) + same_line_overlap
require(high_pair_safe == F(66, 91), "high-pair safe mass")
high_safe = F(6, 7) * high_pair_safe
require(high_safe == F(396, 637), "q/top high-safe factorization")


# ---------------------------------------------------------------------------
# 3. First clean-hole floors.
# ---------------------------------------------------------------------------

seven_floors = {}
for pair in TYPES:
    p_zero = observed_seven[pair].get(0, F())
    seven_floors[pair] = max(F(), high_safe + p_zero - 1) / 7

seven_survivors = tuple(pair for pair in TYPES if seven_floors[pair] == 0)
require(seven_survivors == ((1, 1), (2, 1)), "wrong seven-root residual")

positive_seven = {
    pair: floor for pair, floor in seven_floors.items() if floor > 0
}
minimum_seven_pair = min(positive_seven, key=positive_seven.get)
minimum_seven_floor = positive_seven[minimum_seven_pair]
require(minimum_seven_pair == (1, 2), "wrong smallest seven-root floor type")
require(minimum_seven_floor == F(101, 57967), "wrong seven-root floor")


# ---------------------------------------------------------------------------
# 4. Thirteen-root refinement of the oriented (2,1) boundary.
# ---------------------------------------------------------------------------

empty_13_21 = thirteen_empty_histogram(2, 1)
require(
    empty_13_21 == {0: F(7, 13), 10: F(5, 13), 12: F(1, 13)},
    "wrong (2,1) thirteen-root empty-word law",
)

# Positive/hostile control on the already first-stage-closed orientation.
empty_13_12 = thirteen_empty_histogram(1, 2)
require(
    empty_13_12 == {0: F(1, 2), 10: F(6, 13), 12: F(1, 26)},
    "wrong (1,2) thirteen-root empty-word law",
)

nonempty_base_21 = 1 - empty_13_21[0]
require(nonempty_base_21 == F(6, 13), "wrong nonempty base mass")
base_overlap_floor = high_pair_safe + nonempty_base_21 - 1
require(base_overlap_floor == F(17, 91), "wrong base overlap floor")

# An open length-1/7 arc contains one or two points of a generic 13-grid.
grid_13 = boundary_grid((1,), root_scale=13)
danger_13_hist = exact_histogram(
    grid_13,
    lambda z: sum(danger(1, (z + root) / 13) for root in range(13)),
)
require(danger_13_hist == {1: F(1, 7), 2: F(6, 7)}, "13-grid count")
require(max(danger_13_hist) == 2, "13-grid danger cap")

event_y_floor = F(8, 13) * base_overlap_floor
require(event_y_floor == F(136, 1183), "wrong y-event floor")
second_fibre_floor = event_y_floor / 7
require(second_fibre_floor == F(136, 8281), "wrong clean-hole floor")


# ---------------------------------------------------------------------------
# 5. THM-2391 congruence filter and quantitative descendants.
# ---------------------------------------------------------------------------

mod49 = tuple(
    pair
    for pair in TYPES
    if (13 * pair[0] - pair[1]) % 49 == 0
    or (13 * pair[0] + pair[1]) % 49 == 0
)
mod343 = tuple(
    pair
    for pair in TYPES
    if (13 * pair[0] - pair[1]) % 343 == 0
    or (13 * pair[0] + pair[1]) % 343 == 0
)
require(mod49 == ((4, 3),), "wrong M>=2 congruence survivor")
require(mod343 == (), "wrong M>=3 congruence survivor")
require(seven_floors[(4, 3)] == F(1424, 57967), "M=2 floor")

global_clean_floor = F(1, 26754)
top_cell_floor = global_clean_floor / 52
tensor_cell_floor = global_clean_floor / 338
require(top_cell_floor == F(1, 1391208), "top-labelled cell floor")
require(tensor_cell_floor == F(1, 9042852), "transverse tensor floor")


# ---------------------------------------------------------------------------
# 6. Common-core boundary and optimization-safe positive control.
# ---------------------------------------------------------------------------

require(
    observed_seven[(1, 1)] == {1: F(145, 169), 2: F(24, 169)},
    "common-core boundary acquired an empty fibre",
)

guard_roots = {1, 2}
lower_q_roots = ({3}, {4}, {5}, {6})
blocker_roots = ({0}, {0})
full_incidence = [0] * 7
for root in guard_roots:
    full_incidence[root] += 1
for word in lower_q_roots + blocker_roots:
    for root in word:
        full_incidence[root] += 1
require(full_incidence == [2, 1, 1, 1, 1, 1, 1], "boundary word")
require(sum(full_incidence) == 8, "boundary incidence weight")

# Common dilation does not change any table.  At the finite-word level,
# multiplication by every unit modulo 7 and 13 only permutes root labels.
for unit in range(1, 7):
    require(
        {unit * root % 7 for root in range(7)} == set(range(7)),
        "septimal unit failed to permute roots",
    )
for unit in range(1, 13):
    require(
        {unit * root % 13 for root in range(13)} == set(range(13)),
        "thirteen-unit failed to permute roots",
    )


def fmt_hist(hist: dict[int, F]) -> str:
    return ",".join(f"{key}:{value}" for key, value in hist.items())


print("theorem=THM-2393")
print("status=PROVED-CANDIDATE+VERIFIED-EXACT; independent-audit=PENDING")
print(f"compatible_oriented_types={len(TYPES)}")
print(
    "whole_union_masses="
    + ";".join(f"{a}:{b}={observed_mass[(a,b)]}" for a, b in TYPES)
)
print(f"high_pair_safe={high_pair_safe}; full_high_safe={high_safe}")
print(
    "seven_empty_masses="
    + ";".join(
        f"{a}:{b}={observed_seven[(a,b)].get(0,F())}" for a, b in TYPES
    )
)
print(
    f"seven_root_closed=8/10; minimum_floor={minimum_seven_floor}"
    f"@{minimum_seven_pair[0]}:{minimum_seven_pair[1]}"
)
print(f"seven_root_residual={seven_survivors}")
print(f"thirteen_empty_2:1={fmt_hist(empty_13_21)}")
print(f"thirteen_grid_danger={fmt_hist(danger_13_hist)}")
print(
    f"second_fibre_event={event_y_floor}; clean_floor_2:1={second_fibre_floor}"
)
print(f"M2_filter={mod49}; M3_filter={mod343}")
print(f"M2_4:3_clean_floor={seven_floors[(4,3)]}")
print(
    f"uniform_noncommon_clean>{global_clean_floor};"
    f" top_cell>{top_cell_floor}; tensor_cell>{tensor_cell_floor}"
)
print("final_no_clean_residual=M1;(a,b)=(1,1);chain=(1,13,13,169)")
print("common_core_seven_hist=" + fmt_hist(observed_seven[(1, 1)]))
print("boundary_word=2,1,1,1,1,1,1")
print("row_decrement=0; ledger=165; LRC(14)=OPEN")
print("all_checks=PASS")
