#!/usr/bin/env python3
"""Exact finite controls for THM-4085.

The universal theorem is proved by a triangular left-characteristic argument.
This companion exhausts bounded instances of that argument, the rational-ray
entropy corollary, the marked-block independence threshold, its first hostile,
and the finite-line cemetery class.  Every calculation is over finite Boolean
sets; there are no probabilistic samples.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations


def fail(label):
    raise RuntimeError(label)


def require(condition, label):
    if not condition:
        fail(label)


def local_rule(rule, left, center, right):
    index = (left << 2) | (center << 1) | right
    return (rule >> index) & 1


def is_left_permutive(rule):
    return all(
        local_rule(rule, 0, center, right)
        != local_rule(rule, 1, center, right)
        for center in (0, 1)
        for right in (0, 1)
    )


def dependency_interval(points):
    return (
        min(site - time for time, site in points),
        max(site + time for time, site in points),
    )


def outputs_from_assignment(rule, points, lower, upper, assignment):
    row = [(assignment >> offset) & 1 for offset in range(upper - lower + 1)]
    base = lower
    by_time = {}
    for index, (time, site) in enumerate(points):
        by_time.setdefault(time, []).append((index, site))
    values = [0] * len(points)
    maximum_time = max(time for time, site in points)
    for time in range(maximum_time + 1):
        for index, site in by_time.get(time, ()):
            values[index] = row[site - base]
        if time != maximum_time:
            row = [
                local_rule(rule, row[index], row[index + 1], row[index + 2])
                for index in range(len(row) - 2)
            ]
            base += 1
    return tuple(values)


def distribution(rule, points):
    lower, upper = dependency_interval(points)
    width = upper - lower + 1
    counts = [0] * (1 << len(points))
    for assignment in range(1 << width):
        values = outputs_from_assignment(rule, points, lower, upper, assignment)
        code = sum(value << index for index, value in enumerate(values))
        counts[code] += 1
    return lower, upper, counts


def check_uniform(rule, points, label):
    lower, upper, counts = distribution(rule, points)
    width = upper - lower + 1
    require(width >= len(points), label + ": width")
    expected = 1 << (width - len(points))
    require(all(count == expected for count in counts), label + ": uniform")
    return 1 << width, len(counts)


def marked_points(radius, time_index):
    return tuple(
        (time_index - radius, site)
        for site in range(-radius, radius - 1)
    )


def marked_target(radius, block_count):
    length = 2 * radius - 1
    return sum(1 << (block * length) for block in range(block_count))


def rule30_formula(left, center, right):
    return left ^ center ^ right ^ (center & right)


def direct_rule30_rows(initial, lower, upper, steps):
    row = {site: initial.get(site, 0) for site in range(lower, upper + 1)}
    rows = [row]
    for step in range(steps):
        previous = rows[-1]
        next_row = {
            site: rule30_formula(
                previous.get(site - 1, 0),
                previous.get(site, 0),
                previous.get(site + 1, 0),
            )
            for site in range(lower + step + 1, upper - step)
        }
        rows.append(next_row)
    return rows


def main():
    semantic = []
    gates = 0
    assignments_checked = 0

    left_permutive_rules = tuple(rule for rule in range(256) if is_left_permutive(rule))
    require(len(left_permutive_rules) == 16, "sixteen left-permutive ECA rules")
    require(30 in left_permutive_rules, "Rule 30 left permutive")

    cells = tuple(
        (time, site)
        for time in range(3)
        for site in range(-time, time + 1)
    )
    point_sets = []
    for size in range(1, 5):
        for points in combinations(cells, size):
            addresses = tuple(site - time for time, site in points)
            if len(set(addresses)) == size:
                point_sets.append(points)
    for rule in left_permutive_rules:
        for index, points in enumerate(point_sets):
            assignments, word_gates = check_uniform(
                rule, points, "left-permutive atlas rule=%d set=%d" % (rule, index)
            )
            assignments_checked += assignments
            gates += word_gates
    print(
        "left_permutive_atlas=rules:%d;cells:%d;distinct_address_sets:%d;"
        "assignments:%d;uniform_word_gates:%d"
        % (
            len(left_permutive_rules),
            len(cells),
            len(point_sets),
            assignments_checked,
            gates,
        )
    )
    semantic.append(("atlas_sets", len(point_sets)))
    semantic.append(("atlas_gates", gates))

    collision_points = ((0, 0), (1, 1))
    lower, upper, collision_counts = distribution(30, collision_points)
    require((lower, upper) == (0, 2), "collision dependency interval")
    require(collision_counts == [1, 3, 3, 1], "collision histogram")
    walsh_numerator = (
        collision_counts[0]
        + collision_counts[3]
        - collision_counts[1]
        - collision_counts[2]
    )
    require(walsh_numerator == -4, "collision Walsh numerator")
    print(
        "same_address_hostile=cells:(0,0),(1,1);address:0;"
        "counts_00_10_01_11:%s;walsh_character:%d/8"
        % (tuple(collision_counts), walsh_numerator)
    )
    semantic.append(("collision", tuple(collision_counts)))

    ray_records = []
    for numerator, denominator in ((-1, 2), (0, 1), (1, 3), (1, 2), (2, 3), (1, 1)):
        horizon = 6
        points = tuple(
            (time, (numerator * time) // denominator)
            for time in range(horizon)
        )
        addresses = tuple(site - time for time, site in points)
        distinct = len(set(addresses))
        chosen = []
        seen_addresses = set()
        for point, address in zip(points, addresses):
            if address not in seen_addresses:
                seen_addresses.add(address)
                chosen.append(point)
        assignments, word_gates = check_uniform(
            30, tuple(chosen), "rational ray selected iid"
        )
        assignments_checked += assignments
        gates += word_gates
        lower, upper, counts = distribution(30, points)
        total = 1 << (upper - lower + 1)
        maximum = max(counts)
        support = sum(count != 0 for count in counts)
        require(maximum * (1 << distinct) <= total, "rational ray min entropy")
        if numerator <= 0:
            require(distinct == horizon, "nonpositive ray distinct addresses")
            require(all(count == total // (1 << horizon) for count in counts), "nonpositive ray iid")
        if 0 <= numerator <= denominator:
            expected_distinct = horizon - (numerator * (horizon - 1)) // denominator
            require(distinct == expected_distinct, "rational ray address count")
        record = (
            numerator,
            denominator,
            tuple(site for time, site in points),
            addresses,
            distinct,
            support,
            maximum,
            total,
        )
        ray_records.append(record)
        print(
            "rational_ray=p/q:%d/%d;sites:%s;addresses:%s;distinct:%d;"
            "support:%d;max_word_mass:%d/%d;bound_denominator:%d"
            % (
                numerator,
                denominator,
                record[2],
                addresses,
                distinct,
                support,
                maximum,
                total,
                1 << distinct,
            )
        )
    semantic.append(("rays", tuple(ray_records)))

    marked_cases = ((1, (1, 2, 3)), (2, (2, 5)), (3, (3, 8)))
    marked_records = []
    for radius, time_indices in marked_cases:
        points = tuple(
            point
            for time_index in time_indices
            for point in marked_points(radius, time_index)
        )
        assignments, word_gates = check_uniform(
            30, points, "marked separated blocks radius=%d" % radius
        )
        assignments_checked += assignments
        gates += word_gates
        lower, upper, counts = distribution(30, points)
        target = marked_target(radius, len(time_indices))
        expected_target_count = 1 << (
            upper - lower + 1 - len(points)
        )
        require(counts[target] == expected_target_count, "marked target count")
        length = 2 * radius - 1
        intervals = tuple(
            (-time_index, 2 * radius - 2 - time_index)
            for time_index in time_indices
        )
        for left_index in range(len(intervals)):
            for right_index in range(left_index + 1, len(intervals)):
                require(
                    intervals[left_index][1] < intervals[right_index][0]
                    or intervals[right_index][1] < intervals[left_index][0],
                    "marked address intervals disjoint",
                )
        record = (
            radius,
            time_indices,
            length,
            intervals,
            counts[target],
            1 << (upper - lower + 1),
        )
        marked_records.append(record)
        print(
            "marked_blocks=r:%d;k:%s;length:%d;address_intervals:%s;"
            "joint_target_mass:%d/%d;independent_denominator:%d"
            % (
                radius,
                time_indices,
                length,
                intervals,
                counts[target],
                1 << (upper - lower + 1),
                1 << (length * len(time_indices)),
            )
        )
    semantic.append(("marked", tuple(marked_records)))

    hostile_points = marked_points(2, 2) + marked_points(2, 4)
    lower, upper, hostile_counts = distribution(30, hostile_points)
    hostile_code = marked_target(2, 2)
    require((lower, upper) == (-4, 2), "gap-two hostile dependency interval")
    require(hostile_counts[hostile_code] == 3, "gap-two hostile target count")
    free_solutions = 0
    formula_gates = 0
    for free_word in range(16):
        a = (free_word >> 0) & 1
        b = (free_word >> 1) & 1
        c = (free_word >> 2) & 1
        d = (free_word >> 3) & 1
        initial = {
            -4: a,
            -3: b,
            -2: 1,
            -1: 0,
            0: 0,
            1: c,
            2: d,
        }
        rows = direct_rule30_rows(initial, -4, 2, 2)
        second_block = tuple(rows[2][site] for site in (-2, -1, 0))
        expected_block = (a, b, 1 ^ c ^ d ^ (c & d))
        require(second_block == expected_block, "gap-two hostile Boolean formula")
        formula_gates += 3
        if second_block == (1, 0, 0):
            free_solutions += 1
    require(free_solutions == 3, "gap-two hostile free solutions")
    gap_three_points = marked_points(2, 2) + marked_points(2, 5)
    lower_three, upper_three, gap_three_counts = distribution(30, gap_three_points)
    require(gap_three_counts[hostile_code] == 8, "gap-three independent target count")
    require(1 << (upper_three - lower_three + 1) == 512, "gap-three denominator")
    print(
        "marked_threshold_hostile=r:2;gap:2;joint_mass:3/128;"
        "independent_mass:2/128;conditioned_second_block:(a,b,1+c+d+cd);"
        "free_solutions:3/16;gap3_mass:8/512"
    )
    gates += formula_gates + len(hostile_counts) + len(gap_three_counts)
    assignments_checked += (1 << (upper - lower + 1)) + (
        1 << (upper_three - lower_three + 1)
    )
    semantic.append(("gap2", hostile_counts[hostile_code], 128, free_solutions))
    semantic.append(("gap3", gap_three_counts[hostile_code], 512))

    cemetery_records = []
    for depth in range(1, 7):
        terminal_points = tuple(
            (depth - address, -address)
            for address in range(1, depth + 1)
        )
        lower, upper, counts = distribution(30, terminal_points)
        all_one_code = (1 << depth) - 1
        require((lower, upper) == (-depth, depth - 2), "cemetery interval")
        require(counts[all_one_code] == 1, "cemetery unique cylinder")
        joint_points = terminal_points + ((depth, 0),)
        joint_lower, joint_upper, joint_counts = distribution(30, joint_points)
        joint_code = (1 << (depth + 1)) - 1
        require((joint_lower, joint_upper) == (-depth, depth), "cemetery center interval")
        require(joint_counts[joint_code] == 1, "cemetery center unique cylinder")
        record = (
            depth,
            1 << (2 * depth - 1),
            1 << (2 * depth + 1),
        )
        cemetery_records.append(record)
        assignments_checked += (1 << (2 * depth - 1)) + (1 << (2 * depth + 1))
        gates += len(counts) + len(joint_counts)
    print(
        "cemetery_classes=k1_to_k6:%s;"
        "mass_forms:P(Z_k=infinity)=2^-(2k-1),P(center1_and_infinity)=2^-(2k+1)"
        % (tuple(cemetery_records),)
    )
    semantic.append(("cemetery", tuple(cemetery_records)))

    digest = sha256(repr(tuple(semantic)).encode("ascii")).hexdigest()
    print("semantic_sha256=" + digest)
    print("assignments_checked=%d" % assignments_checked)
    print("exact_gates=%d" % gates)
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
