#!/usr/bin/env python3
"""Independent exact verifier for the provisional THM-4369 candidate.

This is a clean-room implementation.  It imports only Python's standard
library and imports no repository module or primary checker.  Packet
polynomials are built by convolution rather than binomial-coefficient lookup;
source-cone kernels are certified by integral pivot/elimination structure;
ballot formulas are compared with both a dynamic reducer and literal path
enumeration.
"""

from __future__ import annotations

from hashlib import sha256
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")


CHECKS = 0
LINEAR_CACHE: dict[tuple[int, int], tuple[int, ...]] = {}
G_CACHE: dict[tuple[int, int], tuple[int, ...]] = {}
F_CACHE: dict[tuple[int, int], tuple[int, ...]] = {}
CONSTANT_CACHE: dict[int, tuple[int, int, int]] = {}


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"clean-room check failed: {label}")


def tri(n: int) -> int:
    check(n >= 0, "triangular argument")
    return n * (n + 1) // 2


def choose(n: int, k: int) -> int:
    """Multiplicative exact binomial, with zero outside 0..n."""
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    answer = 1
    for j in range(1, k + 1):
        answer = answer * (n - k + j) // j
    return answer


def trimmed(values: list[int] | tuple[int, ...]) -> tuple[int, ...]:
    values = list(values)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def multiply(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return trimmed(answer)


def linear_power(sign: int, exponent: int) -> tuple[int, ...]:
    check(sign in (-1, 1) and exponent >= 0, "linear power input")
    key = (sign, exponent)
    if key in LINEAR_CACHE:
        return LINEAR_CACHE[key]
    answer = (1,)
    factor = (1, sign)
    for _ in range(exponent):
        answer = multiply(answer, factor)
    LINEAR_CACHE[key] = answer
    return answer


def shifted(poly: tuple[int, ...], amount: int, scale: int = 1) -> tuple[int, ...]:
    check(amount >= 0, "polynomial shift")
    return trimmed([0] * amount + [scale * value for value in poly])


def add_scaled(
    accumulator: tuple[int, ...], poly: tuple[int, ...], scale: int
) -> tuple[int, ...]:
    length = max(len(accumulator), len(poly))
    answer = [0] * length
    for index in range(length):
        answer[index] = (
            (accumulator[index] if index < len(accumulator) else 0)
            + scale * (poly[index] if index < len(poly) else 0)
        )
    return trimmed(answer)


def packet_g(point: tuple[int, int]) -> tuple[int, ...]:
    if point in G_CACHE:
        return G_CACHE[point]
    u, v = point
    check(u >= 1 and v >= 1, "G positive coordinates")
    answer = shifted(linear_power(1, v), u - 1)
    G_CACHE[point] = answer
    return answer


def trace_f(point: tuple[int, int]) -> tuple[int, ...]:
    if point in F_CACHE:
        return F_CACHE[point]
    u, v = point
    check(u >= 1 and v >= 1, "F positive coordinates")
    sign = -1 if (u - 1) % 2 else 1
    answer = shifted(linear_power(-1, v - 1), u - 1, sign)
    F_CACHE[point] = answer
    return answer


def relation_image(
    terms: dict[tuple[int, int], int], maker
) -> tuple[int, ...]:
    answer = (0,)
    for point in sorted(terms):
        answer = add_scaled(answer, maker(point), terms[point])
    return answer


def divide_by_one_plus_w(poly: tuple[int, ...]) -> tuple[int, ...]:
    """Exact division by 1+w in ascending coefficient order."""
    check(len(poly) >= 2, "linear division degree")
    quotient = [poly[0]]
    for index in range(1, len(poly) - 1):
        quotient.append(poly[index] - quotient[-1])
    check(quotient[-1] == poly[-1], "exact (1+w) division")
    return trimmed(quotient)


def remove_one_plus_w(poly: tuple[int, ...], times: int) -> tuple[int, ...]:
    answer = poly
    for _ in range(times):
        answer = divide_by_one_plus_w(answer)
    return answer


def order_at_zero(poly: tuple[int, ...]) -> int:
    for index, coefficient in enumerate(poly):
        if coefficient:
            return index
    raise RuntimeError("zero polynomial has no finite valuation")


def order_at_minus_one(poly: tuple[int, ...]) -> int:
    order = 0
    current = poly
    while len(current) > 1 and sum(
        coefficient * ((-1) ** index)
        for index, coefficient in enumerate(current)
    ) == 0:
        current = divide_by_one_plus_w(current)
        order += 1
    return order


def constants(ell: int) -> tuple[int, int, int]:
    if ell not in CONSTANT_CACHE:
        check(ell >= 2, "ell lower bound")
        CONSTANT_CACHE[ell] = (ell + 1) // 2, (ell + 2) // 3, ell // 2
    return CONSTANT_CACHE[ell]


def valid_type(ell: int, point: tuple[int, int]) -> bool:
    u, v = point
    if u < 1 or v < 1:
        return False
    s, rho, floor_half = constants(ell)
    return v >= rho and u + s - 1 >= v and u + v - 1 >= floor_half


def slack(ell: int, point: tuple[int, int]) -> int:
    s, _, _ = constants(ell)
    u, v = point
    return u + s - 1 - v


def source_types_bruteforce(ell: int, cutoff: int) -> list[tuple[int, int]]:
    answer = []
    for tau in range(1, cutoff + 1):
        for u in range(1, tau + 1):
            point = (u, tau + 1 - u)
            if valid_type(ell, point):
                answer.append(point)
    return answer


def source_layer_formula(ell: int, tau: int) -> list[tuple[int, int]]:
    s, rho, floor_half = constants(ell)
    if tau < floor_half:
        return []
    lower = (tau + 2 - s + 1) // 2
    upper = tau + 1 - rho
    return [(u, tau + 1 - u) for u in range(lower, upper + 1)]


def is_boundary(ell: int, cutoff: int, point: tuple[int, int]) -> bool:
    u, v = point
    return u + v - 1 == cutoff or slack(ell, point) == 0


def boundary_formula(ell: int, cutoff: int) -> list[tuple[int, int]]:
    s, rho, _ = constants(ell)
    rank = cutoff - rho + 1
    return [(u, min(u + s - 1, cutoff + 1 - u)) for u in range(1, rank + 1)]


def boundary_split_formula(ell: int, cutoff: int) -> list[tuple[int, int]]:
    s, rho, _ = constants(ell)
    h = (cutoff + 1 - s) // 2
    wall = [(u, u + s - 1) for u in range(1, h + 1)]
    top = [(u, cutoff + 1 - u) for u in range(h + 1, cutoff + 2 - rho)]
    return wall + top


def local_terms(point: tuple[int, int]) -> dict[tuple[int, int], int]:
    u, v = point
    return {point: 1, (u, v + 1): -1, (u + 1, v): 1}


def address(point: tuple[int, int]) -> int:
    u, v = point
    tau = u + v - 1
    return tri(tau - 1) + u


def decode_address(number: int) -> tuple[int, int]:
    check(number >= 1, "positive natural address")
    tau = 1
    while tri(tau) < number:
        tau += 1
    u = number - tri(tau - 1)
    return u, tau + 1 - u


def bareiss_det(matrix: list[list[int]]) -> int:
    """Exact fraction-free determinant, used only as an independent minor test."""
    n = len(matrix)
    check(n >= 1 and all(len(row) == n for row in matrix), "square determinant")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(n - 1):
        pivot_row = next(
            (row for row in range(column, n) if work[row][column] != 0), None
        )
        check(pivot_row is not None, "Bareiss nonzero pivot")
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, n):
            for index in range(column + 1, n):
                numerator = (
                    work[row][index] * pivot
                    - work[row][column] * work[column][index]
                )
                check(numerator % previous == 0, "Bareiss exact division")
                work[row][index] = numerator // previous
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def first_contact_count(depth: int, time: int) -> int:
    check(depth >= 1 and time >= depth and (time - depth) % 2 == 0,
          "first-contact domain")
    rights = (time - depth) // 2
    numerator = depth * choose(time, rights)
    check(numerator % time == 0, "ballot divisibility")
    return numerator // time


def closed_normal_form(
    ell: int, cutoff: int, point: tuple[int, int]
) -> dict[tuple[int, int], int]:
    u, v = point
    depth = slack(ell, point)
    horizon = cutoff - (u + v - 1)
    if depth == 0 or horizon == 0:
        return {point: 1}
    check(depth >= 1 and horizon >= 1, "internal ballot parameters")
    answer: dict[tuple[int, int], int] = {}
    for time in range(depth, horizon):
        if (time - depth) % 2:
            continue
        rights = (time - depth) // 2
        endpoint = (u + rights, v + time - rights)
        answer[endpoint] = ((-1) ** rights) * first_contact_count(depth, time)
    for rights in range(horizon + 1):
        final_depth = depth + 2 * rights - horizon
        if final_depth < 0:
            continue
        if final_depth == 0:
            count = first_contact_count(depth, horizon)
        else:
            count = choose(horizon, rights) - choose(horizon, rights + depth)
        if count:
            endpoint = (u + rights, v + horizon - rights)
            answer[endpoint] = answer.get(endpoint, 0) + ((-1) ** rights) * count
    return {point_: value for point_, value in answer.items() if value}


def dynamic_normal_form(
    ell: int, cutoff: int, point: tuple[int, int]
) -> dict[tuple[int, int], int]:
    active = {point: 1}
    finished: dict[tuple[int, int], int] = {}
    while active:
        next_active: dict[tuple[int, int], int] = {}
        for current in sorted(active):
            coefficient = active[current]
            if is_boundary(ell, cutoff, current):
                finished[current] = finished.get(current, 0) + coefficient
                continue
            u, v = current
            for child, multiplier in (((u, v + 1), 1), ((u + 1, v), -1)):
                next_active[child] = next_active.get(child, 0) + multiplier * coefficient
        active = {key: value for key, value in next_active.items() if value}
    return {key: value for key, value in finished.items() if value}


def literal_absorption_counts(depth: int, horizon: int) -> dict[tuple, int]:
    """Enumerate actual stopped paths; no recurrence or reflection principle."""
    stack = [(0, depth, 0)]
    counts: dict[tuple, int] = {}
    while stack:
        time, current_depth, rights = stack.pop()
        if current_depth == 0:
            key = ("wall", time, rights, 0)
            counts[key] = counts.get(key, 0) + 1
        elif time == horizon:
            key = ("top", time, rights, current_depth)
            counts[key] = counts.get(key, 0) + 1
        else:
            stack.append((time + 1, current_depth - 1, rights))
            stack.append((time + 1, current_depth + 1, rights + 1))
    return counts


def formula_absorption_counts(depth: int, horizon: int) -> dict[tuple, int]:
    counts: dict[tuple, int] = {}
    for time in range(depth, horizon):
        if (time - depth) % 2 == 0:
            rights = (time - depth) // 2
            counts[("wall", time, rights, 0)] = first_contact_count(depth, time)
    for rights in range(horizon + 1):
        final_depth = depth + 2 * rights - horizon
        if final_depth < 0:
            continue
        if final_depth == 0:
            count = first_contact_count(depth, horizon)
            key = ("wall", horizon, rights, 0)
        else:
            count = choose(horizon, rights) - choose(horizon, rights + depth)
            key = ("top", horizon, rights, final_depth)
        if count:
            counts[key] = count
    return counts


def realizer_formula(
    ell: int, point: tuple[int, int]
) -> list[tuple[int, int, int, int]]:
    s, _, _ = constants(ell)
    u, v = point
    n0 = u + s - 1
    lower = max(0, ell - 2 * v)
    upper = min(v, n0 - v)
    return [
        (2 * v + exponent_e - ell,
         n0 - v - exponent_e,
         v - exponent_e,
         exponent_e)
        for exponent_e in range(lower, upper + 1)
    ]


def realizer_bruteforce(
    ell: int, point: tuple[int, int]
) -> list[tuple[int, int, int, int]]:
    """Search every possible e; the defining equations determine a,b,c."""
    s, _, _ = constants(ell)
    u, v = point
    n0 = u + s - 1
    answer = []
    for exponent_e in range(v + 1):
        exponent_c = v - exponent_e
        exponent_b = n0 - exponent_c - 2 * exponent_e
        exponent_a = 2 * exponent_c + 3 * exponent_e - ell
        candidate = (exponent_a, exponent_b, exponent_c, exponent_e)
        if min(candidate) >= 0:
            answer.append(candidate)
    return answer


def type_from_realizer(
    ell: int, realizer: tuple[int, int, int, int]
) -> tuple[int, int]:
    a, b, c, exponent_e = realizer
    check(min(realizer) >= 0, "realizer nonnegative")
    check(a == 2 * c + 3 * exponent_e - ell, "realizer diagonal equation")
    n0 = b + c + 2 * exponent_e
    v = c + exponent_e
    s, _, _ = constants(ell)
    return n0 - s + 1, v


def expanded_packet_from_realizer(
    ell: int, realizer: tuple[int, int, int, int]
) -> dict[tuple[int, int], int]:
    """Return coefficients keyed by (x exponent,t exponent)."""
    a, b, c, exponent_e = realizer
    check(a == 2 * c + 3 * exponent_e - ell, "expanded packet diagonal")
    n0 = b + c + 2 * exponent_e
    total = c + exponent_e
    coefficients = linear_power(1, total)
    return {
        (2 * n0 - ell + 2 * offset, n0 + offset): coefficient
        for offset, coefficient in enumerate(coefficients)
    }


def selected_realisers(values: list[tuple[int, int, int, int]]):
    indices = {0, len(values) // 2, len(values) - 1}
    return [values[index] for index in sorted(indices)]


def ambient_audit() -> tuple[int, int, int]:
    cases = 0
    local_cells = 0
    coefficient_checks = 0
    for cutoff in range(1, 29):
        nodes = [
            (u, tau + 1 - u)
            for tau in range(1, cutoff + 1)
            for u in range(1, tau + 1)
        ]
        parents = [
            (u, tau + 1 - u)
            for tau in range(1, cutoff)
            for u in range(1, tau + 1)
        ]
        top = [(u, cutoff + 1 - u) for u in range(1, cutoff + 1)]
        check(len(nodes) == tri(cutoff), "ambient node count")
        check(len(parents) == tri(cutoff - 1), "ambient parent count")
        check([address(point) for point in nodes] == list(range(1, tri(cutoff) + 1)),
              "ambient address interval")
        check(all(decode_address(address(point)) == point for point in nodes),
              "ambient address decoder")
        for point in parents:
            terms = local_terms(point)
            check(relation_image(terms, packet_g) == (0,), "ambient G circuit")
            check(relation_image(terms, trace_f) == (0,), "ambient F circuit")
            u, v = point
            check(address((u + 1, v)) == address((u, v + 1)) + 1,
                  "ambient consecutive children")
            local_cells += 1
        matrix = []
        for index, point in enumerate(top):
            reduced = remove_one_plus_w(packet_g(point), 1)
            row = list(reduced) + [0] * (cutoff - len(reduced))
            check(all(value == 0 for value in row[:index]) and row[index] == 1,
                  "ambient unit pivot")
            matrix.append(row)
            coefficient_checks += len(row)
        if cutoff <= 12:
            check(bareiss_det(matrix) == 1, "ambient exact unit minor")
        cases += 1
    return cases, local_cells, coefficient_checks


def source_audit() -> tuple[int, int, int, int, int, int, int, int]:
    cases = 0
    node_total = 0
    circuit_total = 0
    boundary_coordinates = 0
    normal_coordinates = 0
    realizer_total = 0
    fibre_difference_total = 0
    lifted_choice_total = 0
    for ell in range(2, 43):
        s, rho, floor_half = constants(ell)
        check(floor_half >= rho, "nonempty minimum rank")
        for cutoff in range(floor_half, floor_half + 15):
            points_list = source_types_bruteforce(ell, cutoff)
            points = set(points_list)
            layered = [
                point
                for tau in range(floor_half, cutoff + 1)
                for point in source_layer_formula(ell, tau)
            ]
            check(points_list == layered, "source layer quantifiers")
            check(all(valid_type(ell, point) for point in points),
                  "source inequalities")
            expected_count = sum(
                tau + 2 - rho - ((tau + 2 - s + 1) // 2)
                for tau in range(floor_half, cutoff + 1)
            )
            check(len(points_list) == expected_count, "source cardinality formula")

            boundary = [point for point in points_list if is_boundary(ell, cutoff, point)]
            internal = [point for point in points_list if not is_boundary(ell, cutoff, point)]
            rank = cutoff - rho + 1
            check(boundary == boundary_formula(ell, cutoff), "min-wall/top boundary")
            check(boundary == boundary_split_formula(ell, cutoff),
                  "strict-wall plus top split")
            check(len(boundary) == rank, "source boundary rank")
            check(len(internal) == len(points_list) - rank, "source nullity")
            check(set(boundary).isdisjoint(internal), "boundary/internal disjoint")

            boundary_matrix: list[list[int]] = []
            for index, point in enumerate(boundary):
                reduced = remove_one_plus_w(packet_g(point), rho)
                row = list(reduced) + [0] * (rank - len(reduced))
                check(len(row) == rank, "boundary polynomial degree bound")
                check(all(value == 0 for value in row[:index]) and row[index] == 1,
                      "boundary integral pivot")
                boundary_matrix.append(row)
                boundary_coordinates += rank
            if rank <= 11:
                check(bareiss_det(boundary_matrix) == 1,
                      "source exact unit determinant")

            internal_index = {point: index for index, point in enumerate(internal)}
            for parent in internal:
                u, v = parent
                children = ((u, v + 1), (u + 1, v))
                check(slack(ell, parent) >= 1 and u + v - 1 < cutoff,
                      "source internal quantifiers")
                check(all(child in points for child in children),
                      "source children retained")
                terms = local_terms(parent)
                check(relation_image(terms, packet_g) == (0,), "source G circuit")
                check(relation_image(terms, trace_f) == (0,), "source F circuit")
                parent_index = internal_index[parent]
                for child in children:
                    if child in internal_index:
                        check(internal_index[child] > parent_index,
                              "integral circuit unit-pivot order")
                check(address(children[1]) == address(children[0]) + 1,
                      "source address adjacency")

            total_realisers = 0
            total_fibre_differences = 0
            fibre_bank: dict[tuple[int, int], list[tuple[int, int, int, int]]] = {}
            for point in points_list:
                formula = realizer_formula(ell, point)
                brute = realizer_bruteforce(ell, point)
                check(formula == brute and bool(brute), "exact nonempty monomial fibre")
                fibre_bank[point] = brute
                first_packet = expanded_packet_from_realizer(ell, brute[0])
                for realizer in brute:
                    check(type_from_realizer(ell, realizer) == point,
                          "realizer reconstructs type")
                    check(expanded_packet_from_realizer(ell, realizer) == first_packet,
                          "whole individual fibre has one expanded packet")
                total_realisers += len(brute)
                total_fibre_differences += len(brute) - 1

            for parent in internal:
                u, v = parent
                triple = (parent, (u, v + 1), (u + 1, v))
                choice_banks = [selected_realisers(fibre_bank[point]) for point in triple]
                for left in choice_banks[0]:
                    for middle in choice_banks[1]:
                        for right in choice_banks[2]:
                            check(
                                add_scaled(
                                    add_scaled(packet_g(type_from_realizer(ell, left)),
                                               packet_g(type_from_realizer(ell, middle)), -1),
                                    packet_g(type_from_realizer(ell, right)), 1,
                                ) == (0,),
                                "independent monomial-fibre circuit lift",
                            )
                            lifted_choice_total += 1

            check(total_fibre_differences + len(internal)
                  == total_realisers - rank,
                  "full individual-monomial kernel rank")

            for point in points_list:
                dynamic = dynamic_normal_form(ell, cutoff, point)
                closed = closed_normal_form(ell, cutoff, point)
                check(dynamic == closed, "dynamic versus ballot normal form")
                check(set(closed).issubset(boundary), "normal form boundary support")
                check(relation_image(closed, packet_g) == packet_g(point),
                      "normal form G equality")
                check(relation_image(closed, trace_f) == trace_f(point),
                      "normal form F equality")
                normal_coordinates += len(closed)

            cases += 1
            node_total += len(points_list)
            circuit_total += len(internal)
            realizer_total += total_realisers
            fibre_difference_total += total_fibre_differences
    return (
        cases,
        node_total,
        circuit_total,
        boundary_coordinates,
        normal_coordinates,
        realizer_total,
        fibre_difference_total,
        lifted_choice_total,
    )


def ballot_edge_audit() -> tuple[int, int]:
    cases = 0
    stopped_paths = 0
    for depth in range(1, 13):
        for horizon in range(1, 13):
            literal = literal_absorption_counts(depth, horizon)
            formula = formula_absorption_counts(depth, horizon)
            check(literal == formula, "literal paths versus reflection/first-contact formula")
            check(all(count > 0 for count in formula.values()), "positive ballot counts")
            stopped_paths += sum(literal.values())
            cases += 1
    # Explicit theorem-boundary branches: already on the wall, already on top,
    # first contact at H=1, d=H, inaccessible parity, and out-of-range binomial.
    check(closed_normal_form(2, 1, (1, 1)) == {(1, 1): 1}, "ell=2 minimal prefix")
    check(closed_normal_form(10, 8, (1, 5)) == {(1, 5): 1}, "wall delta case")
    check(closed_normal_form(10, 8, (3, 6)) == {(3, 6): 1}, "top delta case")
    check(formula_absorption_counts(1, 1)
          == {("wall", 1, 0, 0): 1, ("top", 1, 1, 2): 1},
          "H=1 top-intersection split")
    check(choose(4, -1) == 0 and choose(4, 5) == 0, "out-of-range binomial zero")
    check(("wall", 2, 0, 0) not in formula_absorption_counts(1, 3),
          "parity-inaccessible wall exit absent")
    return cases, stopped_paths


def named_controls() -> list[str]:
    ell = 10
    s, rho, floor_half = constants(ell)
    check((s, rho, floor_half) == (5, 4, 5), "named ell constants")

    parent, vertical, right = (2, 4), (2, 5), (3, 4)
    check(all(valid_type(ell, point) for point in (parent, vertical, right)),
          "named circuit source-valid")
    check([address(point) for point in (parent, vertical, right)] == [12, 17, 18],
          "named circuit addresses")
    check(relation_image(local_terms(parent), packet_g) == (0,),
          "named packet circuit")
    representatives = [realizer_formula(ell, point)[0]
                       for point in (parent, vertical, right)]
    check(representatives == [(0, 0, 2, 2), (0, 1, 5, 0), (0, 1, 2, 2)],
          "named monomial representatives")

    wall_parent, invalid_vertical, valid_right = (1, 5), (1, 6), (2, 5)
    check(valid_type(ell, wall_parent), "wall hostile parent valid")
    check(not valid_type(ell, invalid_vertical), "wall hostile child invalid")
    check(valid_type(ell, valid_right), "wall hostile other child valid")
    check([address(point) for point in (wall_parent, invalid_vertical, valid_right)]
          == [11, 16, 17], "wall hostile addresses")
    check(relation_image(local_terms(wall_parent), packet_g) == (0,),
          "wall hostile remains ambient identity")

    cutoff = 8
    points = source_types_bruteforce(ell, cutoff)
    boundary = [point for point in points if is_boundary(ell, cutoff, point)]
    top_only = [point for point in points if sum(point) - 1 == cutoff]
    internal = [point for point in points if not is_boundary(ell, cutoff, point)]
    check((len(points), len(boundary), len(internal)) == (10, 5, 5),
          "named source dimensions")
    check([address(point) for point in boundary] == [11, 23, 31, 32, 33],
          "named boundary addresses")
    check(len(top_only) == 3 < len(boundary), "top-only boundary hostile count")
    check(all(order_at_zero(packet_g(point)) >= 2 for point in top_only)
          and order_at_zero(packet_g(wall_parent)) == 0,
          "top-only boundary hostile missing valuation")
    probe_normal = closed_normal_form(ell, cutoff, parent)
    expected_normal = {(2, 6): 1, (3, 6): -2, (4, 5): 3, (5, 4): -1}
    check(probe_normal == expected_normal, "named ballot normal form")

    # Valuations recover every type and rule out support-one/two relations.
    valuation_map = {}
    for point in source_types_bruteforce(ell, 13):
        poly = packet_g(point)
        pair = (order_at_zero(poly), order_at_minus_one(poly))
        check(pair == (point[0] - 1, point[1]), "computed two-boundary valuations")
        check(pair not in valuation_map, "valuation injectivity")
        valuation_map[pair] = point

    reflected_hostile = (3, 6)
    reflected = (6, 3)
    check(valid_type(ell, reflected_hostile) and not valid_type(ell, reflected),
          "reciprocal-address source-closure hostile")
    check(address(reflected_hostile) + address(reflected) == 8 * 8 + 1,
          "reciprocal-address block reflection")

    return [
        "universe:ambient_R=1..28;source_ell=2..42;R=f..f+14;"
        "literal_ballot_d,H=1..12",
        "kernel:ambient=all_local_cells;source=source_internal_local_cells;"
        "integral_unit_pivots=yes",
        "image:G=(1+w)^rho*Z[w]_{<=R-rho};"
        "F=(1-z)^(rho-1)*Z[z]_{<=R-rho}",
        "named_circuit:ell=10;types=(2,4)-(2,5)+(3,4);"
        "addresses=12-17+18;monomials=p^2y^2-u*p^5+u*p^2y^2",
        "wall_hostile:ell=10;(1,5)-(1,6)+(2,5);"
        "addresses=11-16+17;(1,6)_not_source",
        "top_only_hostile:ell=10;R=8;top_nodes=3;full_boundary=5;"
        "missing_wall_valuation=0",
        "ballot_probe:ell=10;R=8;(2,4)->"
        + ",".join(f"{point}:{value}" for point, value in sorted(probe_normal.items())),
        "reflection_hostile:ell=10;(3,6)_source;(6,3)_not_source;"
        "address_sum=65",
        "scope:finite_packet_kernel_only;JC(2)=OPEN",
    ]


def main() -> None:
    ambient_cases, ambient_cells, ambient_coefficients = ambient_audit()
    (
        source_cases,
        source_nodes,
        source_circuits,
        boundary_coordinates,
        normal_coordinates,
        realizers,
        fibre_differences,
        lifted_choices,
    ) = source_audit()
    ballot_cases, literal_stopped_paths = ballot_edge_audit()
    ledger = named_controls()
    digest = sha256(("\n".join(ledger) + "\n").encode("ascii")).hexdigest()

    print("THM-4369 clean-room independent verifier: PASS")
    print(f"checks={CHECKS}")
    print(
        f"ambient_cases={ambient_cases} ambient_local_cells={ambient_cells} "
        f"ambient_boundary_coefficients={ambient_coefficients}"
    )
    print(
        f"source_cases={source_cases} source_nodes={source_nodes} "
        f"source_local_circuits={source_circuits}"
    )
    print(
        f"boundary_matrix_coordinates={boundary_coordinates} "
        f"boundary_normal_coordinates={normal_coordinates}"
    )
    print(
        f"source_monomial_realizers={realizers} "
        f"within_fibre_differences={fibre_differences} "
        f"independent_lift_choices={lifted_choices}"
    )
    print(
        f"literal_ballot_cases={ballot_cases} "
        f"literal_stopped_paths={literal_stopped_paths}"
    )
    for line in ledger:
        print(line)
    print(f"semantic_digest={digest}")


if __name__ == "__main__":
    main()
