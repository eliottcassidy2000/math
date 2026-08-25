#!/usr/bin/env python3
"""Independent ANF/triangular-inverse audit for THM-4085.

Unlike the primary direct-orbit census, this script builds algebraic normal
forms in the Boolean quotient x_i^2=x_i.  It checks the unique extreme-left
pivot monomial, replays bounded triangular inverses, derives the gap-two
hostile after symbolic substitution, and audits the cemetery cylinder by
solving its ANF equations exactly.
"""

from hashlib import sha256
from itertools import combinations


ZERO = frozenset()
ONE = frozenset((0,))


def fail(label):
    raise RuntimeError(label)


def require(condition, label):
    if not condition:
        fail(label)


def poly_xor(*polynomials):
    answer = set()
    for polynomial in polynomials:
        for monomial in polynomial:
            if monomial in answer:
                answer.remove(monomial)
            else:
                answer.add(monomial)
    return frozenset(answer)


def poly_mul(left, right):
    answer = set()
    for left_monomial in left:
        for right_monomial in right:
            monomial = left_monomial | right_monomial
            if monomial in answer:
                answer.remove(monomial)
            else:
                answer.add(monomial)
    return frozenset(answer)


def rule30_anf(left, center, right):
    return poly_xor(left, center, right, poly_mul(center, right))


def dependency_interval(points):
    return (
        min(site - time for time, site in points),
        max(site + time for time, site in points),
    )


def point_polynomials(points):
    lower, upper = dependency_interval(points)
    row = [frozenset((1 << offset,)) for offset in range(upper - lower + 1)]
    base = lower
    by_time = {}
    for index, (time, site) in enumerate(points):
        by_time.setdefault(time, []).append((index, site))
    values = [ZERO] * len(points)
    maximum_time = max(time for time, site in points)
    for time in range(maximum_time + 1):
        for index, site in by_time.get(time, ()):
            values[index] = row[site - base]
        if time != maximum_time:
            row = [
                rule30_anf(row[index], row[index + 1], row[index + 2])
                for index in range(len(row) - 2)
            ]
            base += 1
    return lower, upper, tuple(values)


def evaluate(polynomial, assignment):
    value = 0
    for monomial in polynomial:
        if assignment & monomial == monomial:
            value ^= 1
    return value


def substitute(polynomial, fixed_values):
    answer = set()
    for monomial in polynomial:
        survives = True
        reduced = monomial
        for variable, value in fixed_values.items():
            if monomial & variable:
                if value == 0:
                    survives = False
                    break
                reduced &= ~variable
        if survives:
            if reduced in answer:
                answer.remove(reduced)
            else:
                answer.add(reduced)
    return frozenset(answer)


def check_pivot_form(polynomial, pivot, label):
    require(pivot in polynomial, label + ": pivot singleton")
    containing = tuple(monomial for monomial in polynomial if monomial & pivot)
    require(containing == (pivot,), label + ": no nonlinear pivot term")
    lower_variables = pivot - 1
    require(
        all((monomial & lower_variables) == 0 for monomial in polynomial),
        label + ": no coordinate below left address",
    )
    return len(polynomial)


def assignment_from_code(variable_positions, code):
    assignment = 0
    for index, position in enumerate(variable_positions):
        if (code >> index) & 1:
            assignment |= 1 << position
    return assignment


def verify_triangular_family(points, label, exhaustive_inverse):
    lower, upper, polynomials = point_polynomials(points)
    width = upper - lower + 1
    addresses = tuple(site - time for time, site in points)
    require(len(set(addresses)) == len(points), label + ": distinct addresses")
    pivots = tuple(1 << (address - lower) for address in addresses)
    monomial_gates = 0
    for index, polynomial in enumerate(polynomials):
        monomial_gates += check_pivot_form(
            polynomial, pivots[index], label + ": cell=%d" % index
        )
    inverse_replays = 0
    evaluation_gates = 0
    if exhaustive_inverse:
        pivot_positions = {address - lower for address in addresses}
        free_positions = tuple(
            position for position in range(width) if position not in pivot_positions
        )
        solve_order = tuple(
            sorted(range(len(points)), key=lambda index: addresses[index], reverse=True)
        )
        for free_code in range(1 << len(free_positions)):
            free_assignment = assignment_from_code(free_positions, free_code)
            for target in range(1 << len(points)):
                assignment = free_assignment
                for index in solve_order:
                    pivot = pivots[index]
                    value_zero = evaluate(polynomials[index], assignment & ~pivot)
                    value_one = evaluate(polynomials[index], assignment | pivot)
                    evaluation_gates += 2
                    require(value_zero != value_one, label + ": pivot permutation")
                    target_bit = (target >> index) & 1
                    if value_one == target_bit:
                        assignment |= pivot
                    else:
                        assignment &= ~pivot
                    require(
                        evaluate(polynomials[index], assignment) == target_bit,
                        label + ": selected pivot",
                    )
                    evaluation_gates += 1
                for index, polynomial in enumerate(polynomials):
                    require(
                        evaluate(polynomial, assignment) == ((target >> index) & 1),
                        label + ": inverse replay",
                    )
                    evaluation_gates += 1
                inverse_replays += 1
    return width, monomial_gates, inverse_replays, evaluation_gates


def marked_points(radius, time_index):
    return tuple(
        (time_index - radius, site)
        for site in range(-radius, radius - 1)
    )


def main():
    semantic = []
    monomial_gates = 0
    inverse_replays = 0
    evaluation_gates = 0

    cell_records = []
    for time in range(7):
        for site in range(-time, time + 1):
            lower, upper, polynomials = point_polynomials(((time, site),))
            pivot = 1 << (site - time - lower)
            monomial_gates += check_pivot_form(
                polynomials[0], pivot, "single-cell pivot t=%d j=%d" % (time, site)
            )
            cell_records.append((time, site, len(polynomials[0]), upper - lower + 1))
    require(len(cell_records) == 49, "single-cell atlas cardinality")
    maximum_monomials = max(record[2] for record in cell_records)
    require(maximum_monomials == 1360, "depth-six ANF monomial count")
    print(
        "anf_cell_atlas=cells:%d;depth:0..6;max_monomials:%d;"
        "pivot_form:singleton_only"
        % (len(cell_records), maximum_monomials)
    )
    semantic.append(("cell_records", tuple(cell_records)))

    cells = tuple(
        (time, site)
        for time in range(3)
        for site in range(-time, time + 1)
    )
    family_count = 0
    for size in range(1, 4):
        for points in combinations(cells, size):
            if len({site - time for time, site in points}) != size:
                continue
            width, monomials, inverses, evaluations = verify_triangular_family(
                points,
                "ANF triangular atlas family=%d" % family_count,
                True,
            )
            require(width <= 9, "ANF triangular atlas width")
            monomial_gates += monomials
            inverse_replays += inverses
            evaluation_gates += evaluations
            family_count += 1
    require(family_count == 91, "ANF triangular family count")
    print(
        "anf_triangular_inverse_atlas=families:%d;inverse_replays:%d;"
        "evaluation_gates:%d"
        % (family_count, inverse_replays, evaluation_gates)
    )
    semantic.append(("triangular", family_count, inverse_replays, evaluation_gates))

    marked_records = []
    for radius, time_indices, exhaustive in (
        (1, (1, 2, 3), True),
        (2, (2, 5), True),
        (3, (3, 8), False),
    ):
        points = tuple(
            point
            for time_index in time_indices
            for point in marked_points(radius, time_index)
        )
        width, monomials, inverses, evaluations = verify_triangular_family(
            points, "marked ANF radius=%d" % radius, exhaustive
        )
        monomial_gates += monomials
        inverse_replays += inverses
        evaluation_gates += evaluations
        length = 2 * radius - 1
        intervals = tuple(
            (-time_index, 2 * radius - 2 - time_index)
            for time_index in time_indices
        )
        record = (radius, time_indices, length, intervals, width, inverses)
        marked_records.append(record)
        print(
            "anf_marked_blocks=r:%d;k:%s;address_intervals:%s;width:%d;"
            "full_inverse_replays:%d"
            % (radius, time_indices, intervals, width, inverses)
        )
    semantic.append(("marked", tuple(marked_records)))

    collision_points = ((0, 0), (1, 1))
    lower, upper, collision_polynomials = point_polynomials(collision_points)
    require((lower, upper) == (0, 2), "ANF collision interval")
    x0 = 1 << 0
    x1 = 1 << 1
    x2 = 1 << 2
    require(collision_polynomials[0] == frozenset((x0,)), "ANF collision first")
    require(
        collision_polynomials[1] == frozenset((x0, x1, x2, x1 | x2)),
        "ANF collision second",
    )
    collision_counts = [0, 0, 0, 0]
    for assignment in range(8):
        first = evaluate(collision_polynomials[0], assignment)
        second = evaluate(collision_polynomials[1], assignment)
        collision_counts[first | (second << 1)] += 1
        evaluation_gates += 2
    require(collision_counts == [1, 3, 3, 1], "ANF collision histogram")
    print(
        "anf_same_address_hostile=Y0:x0;Y1:x0+x1+x2+x1x2;"
        "counts:%s;nonzero_collision_character:-1/2"
        % (tuple(collision_counts),)
    )
    semantic.append(("collision", tuple(collision_counts)))

    hostile_block = tuple((2, site) for site in (-2, -1, 0))
    lower, upper, hostile_polynomials = point_polynomials(hostile_block)
    require((lower, upper) == (-4, 2), "ANF hostile interval")
    variables = {coordinate: 1 << (coordinate - lower) for coordinate in range(lower, upper + 1)}
    fixed = {variables[-2]: 1, variables[-1]: 0, variables[0]: 0}
    reduced = tuple(substitute(polynomial, fixed) for polynomial in hostile_polynomials)
    expected_reduced = (
        frozenset((variables[-4],)),
        frozenset((variables[-3],)),
        frozenset((0, variables[1], variables[2], variables[1] | variables[2])),
    )
    require(reduced == expected_reduced, "ANF gap-two reduced block")
    free_positions = tuple(coordinate - lower for coordinate in (-4, -3, 1, 2))
    hostile_solutions = 0
    for code in range(16):
        assignment = assignment_from_code(free_positions, code)
        values = tuple(evaluate(polynomial, assignment) for polynomial in reduced)
        evaluation_gates += 3
        if values == (1, 0, 0):
            hostile_solutions += 1
    require(hostile_solutions == 3, "ANF gap-two hostile solutions")
    print(
        "anf_gap_two_hostile=condition:first_block_100;"
        "second_block:(x_-4,x_-3,1+x_1+x_2+x_1x_2);solutions:3/16;"
        "joint_mass:3/128"
    )
    semantic.append(("gap2_solutions", hostile_solutions))

    cemetery_records = []
    for depth in range(1, 7):
        terminal_points = tuple(
            (depth - address, -address)
            for address in range(1, depth + 1)
        )
        lower, upper, polynomials = point_polynomials(terminal_points)
        require((lower, upper) == (-depth, depth - 2), "ANF cemetery interval")
        solutions = []
        for assignment in range(1 << (upper - lower + 1)):
            values = tuple(evaluate(polynomial, assignment) for polynomial in polynomials)
            evaluation_gates += len(polynomials)
            if all(value == 1 for value in values):
                solutions.append(assignment)
        require(solutions == [1], "ANF cemetery marked word")

        joint_points = terminal_points + ((depth, 0),)
        joint_lower, joint_upper, joint_polynomials = point_polynomials(joint_points)
        require((joint_lower, joint_upper) == (-depth, depth), "ANF cemetery joint interval")
        joint_solutions = []
        for assignment in range(1 << (joint_upper - joint_lower + 1)):
            values = tuple(
                evaluate(polynomial, assignment) for polynomial in joint_polynomials
            )
            evaluation_gates += len(joint_polynomials)
            if all(value == 1 for value in values):
                joint_solutions.append(assignment)
        require(joint_solutions == [1], "ANF cemetery center marked word")
        cemetery_records.append(
            (depth, upper - lower + 1, joint_upper - joint_lower + 1)
        )
    print(
        "anf_cemetery=k1_to_k6:%s;unique_words:1_then_all_zero;"
        "infinity_class_retained:true"
        % (tuple(cemetery_records),)
    )
    semantic.append(("cemetery", tuple(cemetery_records)))

    digest = sha256(repr(tuple(semantic)).encode("ascii")).hexdigest()
    print("semantic_sha256=" + digest)
    print("monomial_gates=%d" % monomial_gates)
    print("inverse_replays=%d" % inverse_replays)
    print("evaluation_gates=%d" % evaluation_gates)
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
