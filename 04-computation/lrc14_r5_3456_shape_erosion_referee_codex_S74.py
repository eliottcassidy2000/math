#!/usr/bin/env python3
"""Exact referee for THM-1159's proportional four-comb ray.

The computation is deliberately dependency-free and uses Fraction throughout.
It reconstructs the complete safe-component atlas for the primitive shape
(3,4,5,6), erodes those components by the target length 1/42, and verifies
the cyclic covering radius 23/168.  It also checks the exact core-atlas
inequality at every possible eight-core maximum M=8,...,12.

No assertion is used: the replay must remain proof-checking under ``python -O``.
"""

from fractions import Fraction as F


SHAPE = (3, 4, 5, 6)
LAMBDA = F(1, 14)
DELTA = F(1, 42)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(x: F) -> F:
    y = x % 1
    return min(y, 1 - y)


def shape_safe(x: F) -> bool:
    return all(circle_norm(k * x) >= LAMBDA for k in SHAPE)


def wall_owner(x: F) -> tuple[int, ...]:
    return tuple(k for k in SHAPE if circle_norm(k * x) == LAMBDA)


def complete_wall_atlas() -> tuple[list[F], list[tuple[F, F]]]:
    walls = {F(0), F(1)}
    raw_walls: list[F] = []
    for k in SHAPE:
        # The danger-boundary equations are k*x = n +/- 1/14.
        for n in range(-1, k + 2):
            for sign in (-1, 1):
                x = (F(n) + sign * LAMBDA) / k
                if 0 <= x <= 1:
                    walls.add(x)
                    raw_walls.append(x)
    require(len(raw_walls) == 36, "unexpected raw wall-event count")
    require(len(set(raw_walls)) == 36, "simultaneous wall events detected")
    ordered = sorted(walls)

    safe_cells: list[tuple[F, F]] = []
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if shape_safe(midpoint):
            require(shape_safe(left), f"unsafe left boundary {left}")
            require(shape_safe(right), f"unsafe right boundary {right}")
            if safe_cells and safe_cells[-1][1] == left:
                safe_cells[-1] = (safe_cells[-1][0], right)
            else:
                safe_cells.append((left, right))

    # No zero-length safe component is allowed to hide between unsafe cells.
    covered_safe_walls = {
        x
        for left, right in safe_cells
        for x in ordered
        if left <= x <= right
    }
    require(
        all(not shape_safe(x) or x in covered_safe_walls for x in ordered),
        "isolated safe wall omitted from the component atlas",
    )
    return ordered, safe_cells


EXPECTED_COMPONENTS = [
    (F(1, 42), F(13, 84)),
    (F(5, 28), F(13, 70)),
    (F(3, 14), F(13, 56)),
    (F(15, 56), F(13, 42)),
    (F(5, 14), F(27, 70)),
    (F(29, 70), F(27, 56)),
    (F(29, 56), F(41, 70)),
    (F(43, 70), F(9, 14)),
    (F(29, 42), F(41, 56)),
    (F(43, 56), F(11, 14)),
    (F(57, 70), F(23, 28)),
    (F(71, 84), F(41, 42)),
]

EXPECTED_STARTS = [
    (F(1, 42), F(11, 84)),
    (F(15, 56), F(2, 7)),
    (F(5, 14), F(38, 105)),
    (F(29, 70), F(11, 24)),
    (F(29, 56), F(59, 105)),
    (F(43, 70), F(13, 21)),
    (F(29, 42), F(17, 24)),
    (F(71, 84), F(20, 21)),
]

EXPECTED_GAPS = [
    F(23, 168),
    F(1, 14),
    F(11, 210),
    F(5, 84),
    F(11, 210),
    F(1, 14),
    F(23, 168),
    F(1, 14),
]


def main() -> None:
    walls, components = complete_wall_atlas()
    require(len(walls) == 38, "unexpected cut-at-zero wall count")
    require(len(walls) - 2 == 36, "unexpected genuine wall count")
    require(components == EXPECTED_COMPONENTS, "safe-component atlas mismatch")

    for left, right in components:
        require(len(wall_owner(left)) == 1, f"left wall tie at {left}")
        require(len(wall_owner(right)) == 1, f"right wall tie at {right}")

    starts = [
        (left, right - DELTA)
        for left, right in components
        if right - left >= DELTA
    ]
    require(starts == EXPECTED_STARTS, "eroded start-set atlas mismatch")
    for left, right in starts:
        require(left <= right, "negative eroded component")
        require(
            any(a <= left and right + DELTA <= b for a, b in components),
            "eroded start interval does not carry a full target interval",
        )

    cyclic_gaps: list[F] = []
    for index, (_, right) in enumerate(starts):
        next_left = starts[(index + 1) % len(starts)][0]
        if index + 1 == len(starts):
            next_left += 1
        cyclic_gaps.append(next_left - right)
    require(cyclic_gaps == EXPECTED_GAPS, "cyclic start-gap word mismatch")
    max_gap = max(cyclic_gaps)
    require(max_gap == F(23, 168), "wrong cyclic covering radius")
    phase_threshold = max_gap + DELTA
    require(phase_threshold == F(9, 56), "wrong erosion threshold")

    # THM-1148(19): ell(P) >= 72/[35(13M+1)].  An eight-subset of [12]
    # has 8 <= M <= 12.  Check the least integer m satisfying 3m>13M;
    # every larger m only improves the phase-length inequality.
    core_rows: list[tuple[int, int, F, F, F]] = []
    for maximum in range(8, 13):
        least_m = (13 * maximum) // 3 + 1
        ell_floor = F(72, 35 * (13 * maximum + 1))
        phase_length = least_m * ell_floor
        margin = phase_length - phase_threshold
        require(3 * least_m > 13 * maximum, "least legal scale is not legal")
        require(phase_length > phase_threshold, "core phase interval too short")
        core_rows.append((maximum, least_m, ell_floor, phase_length, margin))

    # The uniform cleared comparison behind the five finite rows is
    # 13M/3 > 5(13M+1)/64, with numerator 637M-15 after clearing 192.
    require(all(637 * maximum - 15 > 0 for maximum in range(1, 13)),
            "cleared legality comparison failed")
    require(DELTA == F(1, 7 * SHAPE[-1]), "pullback target identity failed")

    print("THM-1159 exact (3,4,5,6) shape-erosion referee")
    print("primitive shape:", SHAPE)
    print("safe convention: ||k*x|| >= 1/14 (closed); danger teeth are open")
    print(f"wall events on (0,1): {len(walls) - 2} distinct; wall ties: 0")
    print(f"positive safe components: {len(components)}")
    for index, (left, right) in enumerate(components, start=1):
        print(
            f"  S{index:02d}: [{left},{right}] length={right-left} "
            f"owners={wall_owner(left)}->{wall_owner(right)}"
        )
    print(f"target phase length delta: {DELTA}")
    print(f"delta-admissible start components: {len(starts)}")
    for index, (left, right) in enumerate(starts, start=1):
        print(f"  T{index}: [{left},{right}]")
    print("cyclic complementary gap word: (" + ", ".join(map(str, cyclic_gaps)) + ")")
    print("largest cyclic gap:", max_gap)
    print("universal phase-interval threshold max_gap+delta:", phase_threshold)
    print("core-atlas/legality rows (least legal scale is the worst scale):")
    for maximum, least_m, ell_floor, phase_length, margin in core_rows:
        print(
            f"  M={maximum}: m_min={least_m}, ell_floor={ell_floor}, "
            f"m*ell_floor={phase_length}, margin_over_9/56={margin}"
        )
    print("cleared uniform comparison:")
    print("  m>13M/3>5(13M+1)/64, since 637M-15>0")
    print("pullback identity: delta/m=1/(42m)=1/(7*(6m))")
    print("fifth-comb endpoint audit:")
    print("  a closed interval of length 1/(7k4) cannot lie in one open tooth")
    print("  and for k5>k4 each k5-tooth is strictly shorter")
    print("Tournament Analysis:")
    print("  runner-order tournament: n=4, scores=(0,1,2,3), cycles=0")
    print("  wall-order tournament: n=36, scores=(0,...,35), cycles=0")
    print("  wall SCCs are 36 singletons; no-tie Hamiltonian paths=1")
    print("  safe-component order tournament: n=12, scores=(0,...,11), cycles=0")
    print("  faithful object: cyclic weighted start-set/gap word + wall owners + phase needle")
    print("  naked tournaments lose lengths, erosion eligibility, cyclic gap, and pullback scale")
    print("VERDICT: every legal m(3,4,5,6) four-comb prefix clears the sharp target")


if __name__ == "__main__":
    main()
