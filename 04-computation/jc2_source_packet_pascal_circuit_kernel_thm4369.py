#!/usr/bin/env python3
"""Primary exact certificate for proposed THM-4369.

Starting from THM-4368's boundary coordinates, this script checks the local
Pascal circuit, the full integral kernel on ambient triangular prefixes, and
the corresponding kernel on the source-realizable cone.  It also verifies the
canonical boundary basis, the signed ballot normal form, natural-address
rewrites, and named source/interior/boundary hostiles.

The implementation imports no repository computation.  Every audit uses
``require`` rather than ``assert`` and therefore remains active under
``python -O``.
"""

from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
from math import comb
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("check failed: " + label)


def ceil_div(a: int, b: int) -> int:
    if b <= 0:
        raise ValueError("positive divisor required")
    return -((-a) // b)


def triangular(n: int) -> int:
    if n < 0:
        raise ValueError("nonnegative triangular index required")
    return n * (n + 1) // 2


def trim(poly: list[int]) -> list[int]:
    result = poly[:]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_add(left: list[int], right: list[int]) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        result[index] = (left[index] if index < len(left) else 0) + (
            right[index] if index < len(right) else 0
        )
    return trim(result)


def poly_scale(value: int, poly: list[int]) -> list[int]:
    return trim([value * coefficient for coefficient in poly])


def poly_linear_combination(
    terms: dict[tuple[int, int], int],
    packet_function,
) -> list[int]:
    result = [0]
    for point, coefficient in terms.items():
        result = poly_add(result, poly_scale(coefficient, packet_function(*point)))
    return result


def g_packet(u: int, v: int) -> list[int]:
    """Coefficients of w^(u-1)(1+w)^v."""
    require(u >= 1 and v >= 1, "G packet positive coordinates")
    result = [0] * (u + v)
    for offset in range(v + 1):
        result[u - 1 + offset] = comb(v, offset)
    return result


def f_trace(u: int, v: int) -> list[int]:
    """Coefficients of (-1)^(u-1) z^(u-1)(1-z)^(v-1)."""
    require(u >= 1 and v >= 1, "F trace positive coordinates")
    result = [0] * (u + v - 1)
    sign = -1 if (u - 1) % 2 else 1
    for offset in range(v):
        result[u - 1 + offset] = sign * ((-1) ** offset) * comb(v - 1, offset)
    return result


def reduced_g_packet(u: int, v: int, rho: int) -> list[int]:
    """G packet after removing the common factor (1+w)^rho."""
    require(v >= rho >= 1, "reduced G factor range")
    result = [0] * (u + v - rho)
    for offset in range(v - rho + 1):
        result[u - 1 + offset] = comb(v - rho, offset)
    return result


def address(u: int, v: int) -> int:
    tau = u + v - 1
    return triangular(tau - 1) + u


def source_constants(ell: int) -> tuple[int, int, int]:
    require(ell >= 2, "source intercept lower bound")
    s = ceil_div(ell, 2)
    rho = ceil_div(ell, 3)
    floor_half = ell - s
    return s, rho, floor_half


def source_valid(ell: int, u: int, v: int) -> bool:
    if u < 1 or v < 1:
        return False
    s, rho, floor_half = source_constants(ell)
    n0 = u + s - 1
    return v >= rho and n0 >= v and u + v - 1 >= floor_half


def source_points(ell: int, R: int) -> list[tuple[int, int]]:
    _, _, floor_half = source_constants(ell)
    result = []
    for tau in range(1, R + 1):
        for u in range(1, tau + 1):
            v = tau + 1 - u
            if source_valid(ell, u, v):
                result.append((u, v))
    require(all(u + v - 1 >= floor_half for u, v in result),
            "source minimum layer")
    return result


def layer_interval(ell: int, tau: int) -> list[tuple[int, int]]:
    s, rho, floor_half = source_constants(ell)
    if tau < floor_half:
        return []
    lower = ceil_div(tau + 2 - s, 2)
    upper = tau + 1 - rho
    if lower > upper:
        return []
    return [(u, tau + 1 - u) for u in range(lower, upper + 1)]


def depth_from_wall(ell: int, u: int, v: int) -> int:
    s, _, _ = source_constants(ell)
    return u + s - 1 - v


def source_boundary(ell: int, R: int, point: tuple[int, int]) -> bool:
    u, v = point
    return u + v - 1 == R or depth_from_wall(ell, u, v) == 0


def source_internal(ell: int, R: int, point: tuple[int, int]) -> bool:
    u, v = point
    return u + v - 1 < R and depth_from_wall(ell, u, v) >= 1


def realizing_monomials(ell: int, u: int, v: int) -> list[tuple[int, int, int, int]]:
    """All (a,b,c,e) producing the source packet type (u,v)."""
    s, _, _ = source_constants(ell)
    n0 = u + s - 1
    lower = max(0, ell - 2 * v)
    upper = min(v, n0 - v)
    result = []
    for exponent_e in range(lower, upper + 1):
        exponent_c = v - exponent_e
        exponent_a = 2 * v + exponent_e - ell
        exponent_b = n0 - v - exponent_e
        result.append((exponent_a, exponent_b, exponent_c, exponent_e))
    return result


def dict_add(
    left: dict[tuple[int, int], int],
    right: dict[tuple[int, int], int],
    right_scale: int = 1,
) -> dict[tuple[int, int], int]:
    result = dict(left)
    for key, value in right.items():
        result[key] = result.get(key, 0) + right_scale * value
        if result[key] == 0:
            del result[key]
    return result


def recursive_normal_forms(
    ell: int,
    R: int,
    points: set[tuple[int, int]],
    boundary: set[tuple[int, int]],
):
    @lru_cache(maxsize=None)
    def normal_form(u: int, v: int) -> tuple[tuple[tuple[int, int], int], ...]:
        point = (u, v)
        require(point in points, "normal-form point in source cone")
        if point in boundary:
            return ((point, 1),)
        down = dict(normal_form(u, v + 1))
        up = dict(normal_form(u + 1, v))
        result = dict_add(down, up, -1)
        return tuple(sorted(result.items()))

    return normal_form


def first_hit_count(depth: int, time: int) -> int:
    require(depth >= 1 and time >= depth and (time - depth) % 2 == 0,
            "first-hit parameter range")
    up_steps = (time - depth) // 2
    numerator = depth * comb(time, up_steps)
    require(numerator % time == 0, "first-hit count integral")
    return numerator // time


def ballot_normal_form(
    ell: int,
    R: int,
    point: tuple[int, int],
) -> dict[tuple[int, int], int]:
    u, v = point
    tau = u + v - 1
    depth = depth_from_wall(ell, u, v)
    horizon = R - tau
    if horizon == 0 or depth == 0:
        return {point: 1}

    result: dict[tuple[int, int], int] = {}

    # Paths first absorbed at the n0=N wall before reaching the top layer.
    for time in range(depth, horizon):
        if (time - depth) % 2:
            continue
        up_steps = (time - depth) // 2
        down_steps = time - up_steps
        endpoint = (u + up_steps, v + down_steps)
        coefficient = ((-1) ** up_steps) * first_hit_count(depth, time)
        result[endpoint] = result.get(endpoint, 0) + coefficient

    # Paths reaching the top without earlier contact with the wall.
    for up_steps in range(horizon + 1):
        final_depth = depth + 2 * up_steps - horizon
        if final_depth < 0:
            continue
        if final_depth == 0:
            count = first_hit_count(depth, horizon)
        else:
            reflected_index = up_steps + depth
            reflected = comb(horizon, reflected_index) if reflected_index <= horizon else 0
            count = comb(horizon, up_steps) - reflected
        if count == 0:
            continue
        endpoint = (u + up_steps, v + horizon - up_steps)
        coefficient = ((-1) ** up_steps) * count
        result[endpoint] = result.get(endpoint, 0) + coefficient

    return {key: value for key, value in result.items() if value}


def rank_mod_prime(matrix: list[list[int]], prime: int) -> int:
    if not matrix:
        return 0
    work = [[entry % prime for entry in row] for row in matrix]
    rows = len(work)
    columns = len(work[0])
    require(all(len(row) == columns for row in work), "rank matrix rectangle")
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], prime - 2, prime)
        work[pivot_row] = [(entry * inverse) % prime
                           for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or work[row][column] == 0:
                continue
            multiplier = work[row][column]
            work[row] = [
                (work[row][index] - multiplier * work[pivot_row][index]) % prime
                for index in range(columns)
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def ambient_triangle_checks() -> tuple[int, int, int]:
    triangle_cases = 0
    relation_coordinates = 0
    total_kernel_rank = 0
    for R in range(1, 41):
        points = [(u, tau + 1 - u)
                  for tau in range(1, R + 1)
                  for u in range(1, tau + 1)]
        parents = [(u, tau + 1 - u)
                   for tau in range(1, R)
                   for u in range(1, tau + 1)]
        top = [(u, R + 1 - u) for u in range(1, R + 1)]
        require(len(points) == triangular(R), "ambient triangular node count")
        require(len(parents) == triangular(R - 1), "ambient circuit count")
        require([address(*point) for point in points]
                == list(range(1, triangular(R) + 1)),
                "ambient address prefix")

        for u, v in parents:
            g_relation = {
                (u, v): 1,
                (u, v + 1): -1,
                (u + 1, v): 1,
            }
            f_relation = dict(g_relation)
            require(poly_linear_combination(g_relation, g_packet) == [0],
                    "ambient G Pascal circuit")
            require(poly_linear_combination(f_relation, f_trace) == [0],
                    "ambient F Pascal circuit")
            tau = u + v - 1
            parent_address = address(u, v)
            child_address = address(u, v + 1)
            require(parent_address == triangular(tau - 1) + u,
                    "ambient parent block address")
            require(child_address == triangular(tau) + u
                    and address(u + 1, v) == child_address + 1,
                    "ambient consecutive child addresses")
            relation_coordinates += len(g_packet(u, v)) + len(f_trace(u, v))

        # After removing (1+w), the top packets are a unitriangular
        # Bernstein basis of Z[w]_{<=R-1}.
        for index, (u, v) in enumerate(top):
            reduced = reduced_g_packet(u, v, 1)
            require(all(coefficient == 0 for coefficient in reduced[:index]),
                    "ambient top leading zeros")
            require(reduced[index] == 1, "ambient top unit diagonal")
            require(len(reduced) <= R, "ambient top degree bound")
        triangle_cases += 1
        total_kernel_rank += len(parents)
    return triangle_cases, relation_coordinates, total_kernel_rank


def source_cone_checks() -> tuple[int, int, int, int, int, int, int, int]:
    cone_cases = 0
    source_nodes = 0
    source_relations = 0
    normal_form_coordinates = 0
    ballot_abs_sum = 0
    rank_cases = 0
    total_realizers = 0
    within_fibre_relations = 0

    for ell in range(2, 61):
        s, rho, floor_half = source_constants(ell)
        for R in range(floor_half, floor_half + 25):
            points_list = source_points(ell, R)
            points = set(points_list)
            expected = [point
                        for tau in range(floor_half, R + 1)
                        for point in layer_interval(ell, tau)]
            require(points_list == expected, "source layer interval formula")
            require(all(source_valid(ell, *point) for point in points),
                    "source point feasibility")

            boundary_list = [point for point in points_list
                             if source_boundary(ell, R, point)]
            internal_list = [point for point in points_list
                             if source_internal(ell, R, point)]
            boundary = set(boundary_list)
            require(boundary.isdisjoint(internal_list),
                    "source boundary/internal disjoint")
            require(boundary | set(internal_list) == points,
                    "source boundary/internal exhaustive")
            require(len(boundary_list) == R - rho + 1,
                    "source boundary rank count")
            require(len(internal_list) == len(points_list) - (R - rho + 1),
                    "source circuit nullity count")

            # The canonical boundary has exactly one point for every initial
            # w-valuation, and its reduced G matrix is upper unitriangular.
            boundary_by_u = sorted(boundary_list)
            require([u for u, _ in boundary_by_u]
                    == list(range(1, R - rho + 2)),
                    "source boundary consecutive valuations")
            for index, (u, v) in enumerate(boundary_by_u):
                reduced = reduced_g_packet(u, v, rho)
                require(index == u - 1, "source boundary valuation index")
                require(all(coefficient == 0 for coefficient in reduced[:index]),
                        "source boundary leading zeros")
                require(reduced[index] == 1, "source boundary unit diagonal")
                require(len(reduced) <= R - rho + 1,
                        "source boundary quotient degree")

            normal_form = recursive_normal_forms(ell, R, points, boundary)
            case_realizers = 0
            case_fibre_relations = 0
            for point in points_list:
                recursive = dict(normal_form(*point))
                ballot = ballot_normal_form(ell, R, point)
                require(recursive == ballot, "recursive equals ballot normal form")
                require(set(recursive).issubset(boundary),
                        "normal form supported on boundary")
                require(poly_linear_combination(recursive, g_packet)
                        == g_packet(*point), "G boundary normal form")
                require(poly_linear_combination(recursive, f_trace)
                        == f_trace(*point), "F boundary normal form")
                normal_form_coordinates += len(recursive)
                ballot_abs_sum += sum(abs(value) for value in recursive.values())

                # THM-4368's complete source exponent fibre.
                u, v = point
                n0 = u + s - 1
                realizers = realizing_monomials(ell, u, v)
                lower = max(0, ell - 2 * v)
                upper = min(v, n0 - v)
                require(len(realizers) == upper - lower + 1 and realizers,
                        "complete source fibre multiplicity")
                for exponent_a, exponent_b, exponent_c, exponent_e in realizers:
                    require(min(exponent_a, exponent_b,
                                exponent_c, exponent_e) >= 0,
                            "source exponents nonnegative")
                    require(exponent_c + exponent_e == v
                            and exponent_b + exponent_c + 2 * exponent_e == n0,
                            "source type reconstruction")
                    require(exponent_a
                            == 2 * exponent_c + 3 * exponent_e - ell,
                            "source diagonal reconstruction")
                case_realizers += len(realizers)
                case_fibre_relations += len(realizers) - 1

            for u, v in internal_list:
                children = ((u, v + 1), (u + 1, v))
                require(all(child in points for child in children),
                        "source internal circuit closes")
                relation = {(u, v): 1, children[0]: -1, children[1]: 1}
                require(poly_linear_combination(relation, g_packet) == [0],
                        "source G local circuit")
                require(poly_linear_combination(relation, f_trace) == [0],
                        "source F local circuit")
                require(address(*children[1]) == address(*children[0]) + 1,
                        "source circuit consecutive children")

            # Direct modular ranks are a separate bounded check of the image
            # rank, independent of the rewrite recursion.
            if ell <= 24 and R <= floor_half + 10:
                coefficient_rows = []
                for degree in range(R + 1):
                    coefficient_rows.append([
                        g_packet(*point)[degree]
                        if degree < len(g_packet(*point)) else 0
                        for point in points_list
                    ])
                for prime in (101, 1009):
                    require(rank_mod_prime(coefficient_rows, prime)
                            == R - rho + 1,
                            "source direct modular image rank")
                rank_cases += 1

            require(case_fibre_relations + len(internal_list)
                    == case_realizers - (R - rho + 1),
                    "full monomial-kernel rank splits by fibre and circuit")

            cone_cases += 1
            source_nodes += len(points_list)
            source_relations += len(internal_list)
            total_realizers += case_realizers
            within_fibre_relations += case_fibre_relations

    return (cone_cases, source_nodes, source_relations,
            normal_form_coordinates, ballot_abs_sum, rank_cases,
            total_realizers, within_fibre_relations)


def named_controls() -> list[str]:
    ell = 10
    s, rho, floor_half = source_constants(ell)
    require((s, rho, floor_half) == (5, 4, 5), "named ell=10 constants")

    # First small source-internal circuit used in the theorem.
    parent = (2, 4)
    down = (2, 5)
    up = (3, 4)
    require(all(source_valid(ell, *point) for point in (parent, down, up)),
            "named source circuit feasibility")
    require([address(*point) for point in (parent, down, up)] == [12, 17, 18],
            "named source circuit addresses")
    relation = {parent: 1, down: -1, up: 1}
    require(poly_linear_combination(relation, g_packet) == [0],
            "named source packet identity")
    require(poly_linear_combination(relation, f_trace) == [0],
            "named trace identity")

    named_types = []
    for u, v in (parent, down, up):
        n0 = u + s - 1
        exponent_e = max(0, ell - 2 * v)
        exponent_c = v - exponent_e
        exponent_a = 2 * v + exponent_e - ell
        exponent_b = n0 - v - exponent_e
        named_types.append((u, v, n0, (exponent_a, exponent_b,
                                      exponent_c, exponent_e)))
    require(named_types == [
        (2, 4, 6, (0, 0, 2, 2)),
        (2, 5, 6, (0, 1, 5, 0)),
        (3, 4, 7, (0, 1, 2, 2)),
    ], "named source monomial representatives")

    # The adjacent ambient circuit based at the source wall is not a source
    # circuit, because its v-child violates n0>=N.
    wall_parent = (1, 5)
    invalid_child = (1, 6)
    valid_child = (2, 5)
    require(source_valid(ell, *wall_parent), "named wall parent valid")
    require(not source_valid(ell, *invalid_child), "named wall child invalid")
    require(source_valid(ell, *valid_child), "named wall other child valid")
    require(depth_from_wall(ell, *wall_parent) == 0,
            "named source wall depth zero")
    require([address(*point) for point in
             (wall_parent, invalid_child, valid_child)] == [11, 16, 17],
            "named wall-crossing addresses")

    # No two distinct packet types can be proportional: their two boundary
    # valuations recover (u,v).  The displayed three-term circuit is therefore
    # support-minimal.
    sample = source_points(ell, 12)
    valuation_pairs = {(u - 1, v): (u, v) for u, v in sample}
    require(len(valuation_pairs) == len(sample),
            "two valuations separate source packet types")

    R = 8
    points = source_points(ell, R)
    boundary = [point for point in points if source_boundary(ell, R, point)]
    internal = [point for point in points if source_internal(ell, R, point)]
    require(len(points) == 10 and len(boundary) == 5 and len(internal) == 5,
            "named ell=10 R=8 dimensions")
    require([address(*point) for point in sorted(boundary)]
            == [11, 23, 31, 32, 33], "named source boundary addresses")

    normal_form = recursive_normal_forms(ell, R, set(points), set(boundary))
    probe = (2, 4)
    probe_normal = dict(normal_form(*probe))
    require(probe_normal == ballot_normal_form(ell, R, probe),
            "named ballot decomposition")

    return [
        "local_source_circuit:ell=10;addresses=12-17+18;"
        "types=(N,n0):(4,6)-(5,6)+(4,7);"
        "monomials=p^2y^2-u p^5+u p^2y^2",
        "wall_crossing_hostile:ell=10;addresses=11-16+17;"
        "pair=(1,6)_invalid_since_n0=5<N=6",
        "source_prefix:ell=10,R=8;nodes=10;rank=5;nullity=5;"
        "boundary_addresses=11,23,31,32,33",
        "ballot_probe:ell=10,R=8,parent=(2,4);normal="
        + ",".join(f"{point}:{value}" for point, value in sorted(probe_normal.items())),
    ]


def main() -> None:
    ambient_cases, ambient_coordinates, ambient_kernel_rank = ambient_triangle_checks()
    (cone_cases, source_nodes, source_relations, normal_form_coordinates,
     ballot_abs_sum, rank_cases, total_realizers,
     within_fibre_relations) = source_cone_checks()
    ledger = named_controls()
    semantic_digest = sha256(("\n".join(ledger) + "\n").encode("ascii")).hexdigest()

    print("THM-4369 source-packet Pascal circuit primary: PASS")
    print(f"checks={CHECKS}")
    print(f"ambient_triangle_cases={ambient_cases} "
          f"relation_coordinates={ambient_coordinates} "
          f"summed_kernel_rank={ambient_kernel_rank}")
    print(f"source_cone_cases={cone_cases} source_nodes={source_nodes} "
          f"source_local_relations={source_relations}")
    print(f"boundary_normal_form_coordinates={normal_form_coordinates} "
          f"ballot_abs_coefficient_sum={ballot_abs_sum}")
    print(f"direct_modular_rank_cases={rank_cases} primes=101,1009")
    print(f"source_monomial_realizers={total_realizers} "
          f"within_type_fibre_relations={within_fibre_relations}")
    print("ambient_kernel=local_cells;source_kernel=source-internal_cells")
    print("source_image_G=(1+w)^rho*Z[w]_{<=R-rho}")
    print("source_image_F=(1-z)^(rho-1)*Z[z]_{<=R-rho}")
    for line in ledger:
        print(line)
    print(f"semantic_digest={semantic_digest}")


if __name__ == "__main__":
    main()
