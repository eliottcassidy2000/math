"""Numerical monodromy scout for the old L divisor in the fixed Keller tower.

This is deliberately a numerical discovery tool, not a proof companion.  It
continues all 3^n inverse points around a small transverse loop about the
generic point (a,b,c)=(2/27,1,1) of L=0 and reports the resulting cycle types.

The inverse chart is the exact one used by THM-3533.  Matching uses the
chordal metric on each affine coordinate, so the two branches tending to
infinity remain numerically visible.  Repeated radii and step counts are
hostile controls against a continuation artefact.  The depth-five run is a
hostile against the closed form suggested by depths one through four.
"""

from __future__ import annotations

import cmath
from dataclasses import dataclass
from math import lcm

import numpy as np
from scipy.optimize import linear_sum_assignment


Point = tuple[complex, complex, complex]


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def l_value(point: Point) -> complex:
    a, b, c = point
    return 27 * a * a * c * c - 18 * a * b * c + 16 * a + b**3 * c - b * b


def forward(point: Point) -> Point:
    x, y, z = point
    unit = 1 + x * y
    return (
        unit**3 * z + y * y * unit * (4 + 3 * x * y),
        y + 3 * x * unit * unit * z + 3 * x * y * y * (4 + 3 * x * y),
        2 * x - 3 * x * x * y - x**3 * z,
    )


def relative_error(left: Point, right: Point) -> float:
    return max(
        abs(x - y) / max(1.0, abs(x), abs(y)) for x, y in zip(left, right)
    )


def inverse_points(target: Point) -> tuple[list[Point], float]:
    a, b, c = target
    lead = l_value(target)
    linear = 4 - 3 * b * c
    scale = max(abs(lead), abs(linear), abs(2 * c), 1.0)
    roots = np.roots(
        np.asarray([lead / scale, 0.0, linear / scale, -2 * c / scale], dtype=complex)
    )
    answer: list[Point] = []
    residual = 0.0
    for raw_x in roots:
        x = complex(raw_x)
        denominator = 12 * a * x * x - b * b * x * x + b * x + 2
        require(abs(denominator) > 1.0e-13, ("inverse denominator", target, x))
        y = b - 3 * a * x * (9 * a * c * x - b * x + 2) / denominator
        require(abs(x) > 1.0e-13, ("inverse x", target, x))
        z = (2 * x - 3 * x * x * y - c) / x**3
        point = (x, y, z)
        residual = max(residual, relative_error(forward(point), target))
        answer.append(point)
    require(len(answer) == 3, ("inverse degree", len(answer)))
    return answer, residual


def inverse_levels(target: Point, depth: int) -> tuple[list[list[Point]], float]:
    levels: list[list[Point]] = []
    parents = [target]
    residual = 0.0
    for _ in range(depth):
        children: list[Point] = []
        for parent in parents:
            new_children, new_residual = inverse_points(parent)
            children.extend(new_children)
            residual = max(residual, new_residual)
        levels.append(children)
        parents = children
    return levels, residual


def continue_levels(
    target: Point, previous: list[list[Point]], depth: int
) -> tuple[list[list[Point]], float, float]:
    """Continue in ancestry blocks using only local three-by-three matches."""

    levels: list[list[Point]] = []
    residual = 0.0
    max_step_cost = 0.0
    first_children, first_residual = inverse_points(target)
    first_ordered, first_cost = reorder_by_continuation(previous[0], first_children)
    levels.append(first_ordered)
    residual = max(residual, first_residual)
    max_step_cost = max(max_step_cost, first_cost)
    for level_index in range(1, depth):
        children: list[Point] = []
        for parent_index, parent in enumerate(levels[level_index - 1]):
            raw_children, new_residual = inverse_points(parent)
            old_block = previous[level_index][
                3 * parent_index : 3 * parent_index + 3
            ]
            ordered, step_cost = reorder_by_continuation(old_block, raw_children)
            children.extend(ordered)
            residual = max(residual, new_residual)
            max_step_cost = max(max_step_cost, step_cost)
        levels.append(children)
    return levels, residual, max_step_cost


def projective_pair(value: complex) -> tuple[complex, float]:
    size = abs(value)
    normalizer = np.hypot(1.0, size)
    return value / normalizer, 1.0 / normalizer


def chordal(left: complex, right: complex) -> float:
    left_num, left_den = projective_pair(left)
    right_num, right_den = projective_pair(right)
    return abs(left_num * right_den - right_num * left_den)


def point_distance(left: Point, right: Point) -> float:
    return float(sum(chordal(x, y) ** 2 for x, y in zip(left, right)) ** 0.5)


def distance_matrix(left: list[Point], right: list[Point]) -> np.ndarray:
    return np.asarray(
        [[point_distance(x, y) for y in right] for x in left], dtype=float
    )


def reorder_by_continuation(
    previous: list[Point], current: list[Point]
) -> tuple[list[Point], float]:
    costs = distance_matrix(previous, current)
    rows, columns = linear_sum_assignment(costs)
    require(list(rows) == list(range(len(previous))), ("assignment rows", rows))
    column_for_row = np.empty(len(previous), dtype=int)
    column_for_row[rows] = columns
    selected = [current[int(column_for_row[row])] for row in range(len(previous))]
    return selected, float(max(costs[row, column_for_row[row]] for row in rows))


def endpoint_permutation(initial: list[Point], final: list[Point]) -> tuple[int, ...]:
    costs = distance_matrix(final, initial)
    rows, columns = linear_sum_assignment(costs)
    require(list(rows) == list(range(len(final))), ("endpoint rows", rows))
    permutation = np.empty(len(final), dtype=int)
    permutation[rows] = columns
    require(float(max(costs[row, permutation[row]] for row in rows)) < 2.0e-6,
            ("endpoint mismatch", float(np.max(costs))))
    return tuple(int(value) for value in permutation)


def endpoint_tree_permutations(
    initial: list[list[Point]], final: list[list[Point]]
) -> tuple[tuple[int, ...], ...]:
    permutations = [endpoint_permutation(initial[0], final[0])]
    for level_index in range(1, len(initial)):
        parent_permutation = permutations[-1]
        child_permutation = [-1] * len(initial[level_index])
        for parent_index, parent_image in enumerate(parent_permutation):
            final_block = final[level_index][
                3 * parent_index : 3 * parent_index + 3
            ]
            initial_block = initial[level_index][
                3 * parent_image : 3 * parent_image + 3
            ]
            local = endpoint_permutation(initial_block, final_block)
            for child_index, child_image in enumerate(local):
                child_permutation[3 * parent_index + child_index] = (
                    3 * parent_image + child_image
                )
        require(-1 not in child_permutation, ("incomplete endpoint tree", level_index))
        permutations.append(tuple(child_permutation))
    return tuple(permutations)


def cycle_lengths(permutation: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((len(cycle) for cycle in permutation_cycles(permutation)), reverse=True))


def permutation_cycles(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    seen = [False] * len(permutation)
    cycles: list[tuple[int, ...]] = []
    for start in range(len(permutation)):
        if seen[start]:
            continue
        cursor = start
        cycle: list[int] = []
        while not seen[cursor]:
            seen[cursor] = True
            cycle.append(cursor)
            cursor = permutation[cursor]
        cycles.append(tuple(cycle))
    return tuple(cycles)


def compose(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[int, ...]:
    """Return left after right, in the source-to-image convention."""

    require(len(left) == len(right), ("composition degrees", len(left), len(right)))
    return tuple(left[right[index]] for index in range(len(right)))


def permutation_order(permutation: tuple[int, ...]) -> int:
    answer = 1
    for cycle in permutation_cycles(permutation):
        answer = lcm(answer, len(cycle))
    return answer


def cycle_histogram(permutation: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    lengths = cycle_lengths(permutation)
    return tuple((length, lengths.count(length)) for length in sorted(set(lengths)))


def root_action_and_sections(
    permutation: tuple[int, ...], depth: int
) -> tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    """Decompose an ancestry-preserving ternary permutation at the root."""

    require(len(permutation) == 3**depth, ("root section degree", depth, len(permutation)))
    block_size = 3 ** (depth - 1)
    root_action: list[int] = []
    sections: list[tuple[int, ...]] = []
    for source_root in range(3):
        source_start = source_root * block_size
        images = permutation[source_start : source_start + block_size]
        target_roots = {image // block_size for image in images}
        require(len(target_roots) == 1, ("root block image", source_root, target_roots))
        target_root = next(iter(target_roots))
        root_action.append(target_root)
        section = tuple(image % block_size for image in images)
        require(sorted(section) == list(range(block_size)),
                ("root section permutation", source_root, section))
        sections.append(section)
    require(sorted(root_action) == [0, 1, 2], ("root action", root_action))
    return tuple(root_action), tuple(sections)


def reflection_rotation_row(
    permutation: tuple[int, ...], depth: int
) -> tuple[
    int,
    tuple[int, ...],
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int], ...],
    int,
]:
    """Expose the two subtree involutions and their rotation product.

    The observed root action swaps blocks zero and one and fixes block two.
    If A and B are the corresponding source sections, then g^2 on block zero
    is C=B after A.  Every C-cycle therefore lifts to a g-cycle of twice its
    length, while the third block is fixed pointwise.
    """

    root_action, sections = root_action_and_sections(permutation, depth)
    require(root_action == (1, 0, 2), ("old-L root action", depth, root_action))
    reflection_a, reflection_b, fixed_section = sections
    identity = tuple(range(len(fixed_section)))
    require(fixed_section == identity, ("third root section", depth, fixed_section))
    require(compose(reflection_a, reflection_a) == identity,
            ("A is not an involution", depth))
    require(compose(reflection_b, reflection_b) == identity,
            ("B is not an involution", depth))
    rotation = compose(reflection_b, reflection_a)
    lifted_lengths = tuple(
        sorted((2 * length for length in cycle_lengths(rotation)), reverse=True)
    )
    global_nonfixed = tuple(length for length in cycle_lengths(permutation) if length > 1)
    require(lifted_lengths == global_nonfixed,
            ("reflection-rotation lift", depth, lifted_lengths, global_nonfixed))
    expected_exponent = 2 * len(rotation) - len(permutation_cycles(rotation))
    actual_exponent = len(permutation) - len(permutation_cycles(permutation))
    require(expected_exponent == actual_exponent,
            ("reflection-rotation exponent", depth, expected_exponent, actual_exponent))
    return (
        depth - 1,
        root_action,
        cycle_histogram(reflection_a),
        cycle_histogram(reflection_b),
        cycle_histogram(rotation),
        permutation_order(rotation),
    )


def lift_profile(
    parent: tuple[int, ...], child: tuple[int, ...]
) -> tuple[tuple[int, tuple[int, ...], int], ...]:
    require(len(child) == 3 * len(parent), ("lift degrees", len(parent), len(child)))
    profile: dict[tuple[int, tuple[int, ...]], int] = {}
    predicted_child_cycles: list[int] = []
    for cycle in permutation_cycles(parent):
        action = tuple(range(3))
        for parent_index in cycle:
            child_action = []
            for child_index in range(3):
                image = child[3 * parent_index + child_index]
                require(
                    image // 3 == parent[parent_index],
                    ("block projection", parent_index, child_index, image),
                )
                child_action.append(image % 3)
            require(sorted(child_action) == [0, 1, 2],
                    ("local child permutation", child_action))
            action = tuple(child_action[action[index]] for index in range(3))
        action_type = cycle_lengths(action)
        key = (len(cycle), action_type)
        profile[key] = profile.get(key, 0) + 1
        predicted_child_cycles.extend(len(cycle) * length for length in action_type)
    require(
        tuple(sorted(predicted_child_cycles, reverse=True)) == cycle_lengths(child),
        ("lift cycle reconstruction", profile, cycle_lengths(child)),
    )
    return tuple(
        (parent_length, action_type, count)
        for (parent_length, action_type), count in sorted(profile.items())
    )


@dataclass(frozen=True)
class Run:
    radius: float
    steps: int
    cycle_rows: tuple[tuple[int, tuple[int, ...], int], ...]
    lift_rows: tuple[tuple[int, tuple[tuple[int, tuple[int, ...], int], ...]], ...]
    permutations: tuple[tuple[int, ...], ...]
    max_step_cost: float
    max_forward_residual: float


def run(radius: float, steps: int, depth: int = 4) -> Run:
    tracked: list[list[Point]] | None = None
    initial: list[list[Point]] | None = None
    max_step_cost = 0.0
    max_forward_residual = 0.0
    for step in range(steps + 1):
        parameter = radius * cmath.exp(2j * np.pi * step / steps)
        target = (complex(2 / 27) + parameter, 1.0 + 0j, 1.0 + 0j)
        if tracked is None:
            levels, residual = inverse_levels(target, depth)
            # Generation order is ancestry-compatible: child i has parent
            # i//3.  Retaining it lets the endpoint permutations expose their
            # exact wreath lift products rather than only global cycle types.
            tracked = [list(level) for level in levels]
            initial = [list(level) for level in tracked]
            max_forward_residual = max(max_forward_residual, residual)
            continue
        tracked, residual, step_cost = continue_levels(target, tracked, depth)
        max_forward_residual = max(max_forward_residual, residual)
        max_step_cost = max(max_step_cost, step_cost)

    require(tracked is not None and initial is not None, "empty continuation")
    rows = []
    permutations = endpoint_tree_permutations(initial, tracked)
    for level_index, permutation in enumerate(permutations, start=1):
        cycles = cycle_lengths(permutation)
        exponent = len(permutation) - len(cycles)
        rows.append((level_index, cycles, exponent))
    lift_rows = tuple(
        (depth, lift_profile(parent, child))
        for depth, (parent, child) in enumerate(
            zip(permutations, permutations[1:]), start=1
        )
    )
    return Run(
        radius=radius,
        steps=steps,
        cycle_rows=tuple(rows),
        lift_rows=lift_rows,
        permutations=permutations,
        max_step_cost=max_step_cost,
        max_forward_residual=max_forward_residual,
    )


def main() -> None:
    runs = (
        run(1.0e-3, 240),
        run(3.0e-4, 480),
        run(1.0e-4, 720),
    )
    baseline = runs[0].cycle_rows
    baseline_lifts = runs[0].lift_rows
    for item in runs[1:]:
        require(item.cycle_rows == baseline, ("inconsistent cycle rows", runs))
        require(item.lift_rows == baseline_lifts, ("inconsistent lift rows", runs))
    require(baseline[0][1] == (2, 1), ("level-one positive control", baseline[0]))
    require(all(row[2] % 2 == (1 if row[0] == 1 else 0) for row in baseline),
            ("THM-3531 old-L parity", baseline))
    deep_runs = (
        run(1.0e-3, 360, depth=5),
        run(3.0e-4, 720, depth=5),
    )
    require(deep_runs[0].cycle_rows == deep_runs[1].cycle_rows,
            ("inconsistent depth-five cycle rows", deep_runs))
    require(deep_runs[0].lift_rows == deep_runs[1].lift_rows,
            ("inconsistent depth-five lift rows", deep_runs))
    require(deep_runs[0].cycle_rows[:4] == baseline,
            ("depth-five prefix mismatch", deep_runs[0].cycle_rows, baseline))
    reflection_rows = tuple(
        reflection_rotation_row(permutation, depth)
        for depth, permutation in enumerate(deep_runs[0].permutations, start=1)
        if depth >= 2
    )
    hostile_reflection_rows = tuple(
        reflection_rotation_row(permutation, depth)
        for depth, permutation in enumerate(deep_runs[1].permutations, start=1)
        if depth >= 2
    )
    require(reflection_rows == hostile_reflection_rows,
            ("inconsistent reflection-rotation rows", reflection_rows,
             hostile_reflection_rows))
    rotation_orbit_counts = tuple(
        sum(count for _, count in rotation_histogram)
        for _, _, _, _, rotation_histogram, _ in reflection_rows
    )
    require(rotation_orbit_counts == (2, 4, 8, 20), rotation_orbit_counts)
    a_fixed_counts = tuple(
        dict(a_histogram).get(1, 0)
        for _, _, a_histogram, _, _, _ in reflection_rows
    )
    b_fixed_counts = tuple(
        dict(b_histogram).get(1, 0)
        for _, _, _, b_histogram, _, _ in reflection_rows
    )
    require(a_fixed_counts == (1, 1, 1, 1), a_fixed_counts)
    require(b_fixed_counts == (3, 7, 7, 7), b_fixed_counts)
    depth_five = deep_runs[0].cycle_rows[4]
    candidate_exponent = 2 * 3**4 - 2**4
    require(depth_five[2] == 142, depth_five)
    require(depth_five[2] != candidate_exponent,
            ("closed-form hostile failed", depth_five, candidate_exponent))
    print("== fixed Keller old-L numerical inertia scout ==")
    for item in runs:
        print(
            f"radius={item.radius:.1e};steps={item.steps};"
            f"max_step_cost={item.max_step_cost:.3e};"
            f"max_forward_residual={item.max_forward_residual:.3e}"
        )
        for depth, cycles, exponent in item.cycle_rows:
            histogram = {length: cycles.count(length) for length in sorted(set(cycles))}
            print(
                f"depth={depth};degree={3**depth};cycles={histogram};"
                f"tame_permutation_exponent={exponent};parity={exponent % 2}"
            )
        for depth, profile in item.lift_rows:
            print(f"lift={depth}->{depth + 1};parent_cycle_child_products={profile}")
    for item in deep_runs:
        depth, cycles, exponent = item.cycle_rows[4]
        histogram = {length: cycles.count(length) for length in sorted(set(cycles))}
        print(
            f"depth5_hostile_radius={item.radius:.1e};steps={item.steps};"
            f"max_step_cost={item.max_step_cost:.3e};"
            f"max_forward_residual={item.max_forward_residual:.3e}"
        )
        print(
            f"depth={depth};degree={3**depth};cycles={histogram};"
            f"tame_permutation_exponent={exponent};parity={exponent % 2}"
        )
        print(
            "lift=4->5;parent_cycle_child_products="
            + str(item.lift_rows[3][1])
        )
    for (
        subtree_depth,
        root_action,
        a_histogram,
        b_histogram,
        rotation_histogram,
        rotation_order,
    ) in reflection_rows:
        rotation_orbits = sum(count for _, count in rotation_histogram)
        exponent = 2 * 3**subtree_depth - rotation_orbits
        print(
            f"reflection_rotation_subtree_depth={subtree_depth};"
            f"root_action={root_action};A_cycles={a_histogram};"
            f"B_cycles={b_histogram};C_equals_B_after_A_cycles={rotation_histogram};"
            f"C_orbits={rotation_orbits};C_order={rotation_order};"
            f"global_exponent=2*3^{subtree_depth}-C_orbits={exponent}"
        )
    print(
        f"depths1to4_closed_form_candidate_at5={candidate_exponent};"
        f"observed_at5={depth_five[2]};verdict=REFUTED_AT_DEPTH_5"
    )
    print(
        "closed_form_hidden_assumption=C_orbits_equals_2^r;"
        f"observed_C_orbits={rotation_orbit_counts};"
        "first_failure_at_r=4_is_20_not_16"
    )
    print("status=VERIFIED-NUMERICAL DISCOVERY ONLY;no exact inertia theorem or index claim")


if __name__ == "__main__":
    main()
