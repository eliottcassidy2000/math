#!/usr/bin/env python3
"""Exact referee for THM-4115's ear-cut Walsh variance refinement.

The proof is all-order.  This referee supplies finite hostile controls by an
independent subset-DP boundary implementation.  It exhausts every labelled
tournament through order five, derives the integral cut weights and
orientation field, and checks the complete Boolean fibre.  It also audits two
large selected-bank parents from THM-4102/4104.

Every executable gate uses ``require`` and survives ``python -O``.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations


EXPECTED_SEMANTIC_SHA256 = "7580aaf0fd320afad386d6401b3e95dea31c91fc2cb3fc67b085773dcdf0836e"


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def fraction_row(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def odd_ceiling(value: Fraction) -> int:
    ceiling = (value.numerator + value.denominator - 1) // value.denominator
    return ceiling if ceiling % 2 else ceiling + 1


def canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def tournament_from_code(order: int, code: int) -> list[int]:
    adjacency = [0] * order
    bit_index = 0
    for left in range(order):
        for right in range(left + 1, order):
            if (code >> bit_index) & 1:
                adjacency[left] |= 1 << right
            else:
                adjacency[right] |= 1 << left
            bit_index += 1
    return adjacency


def add_cut_ear(adjacency: list[int], signature: int) -> list[int]:
    order = len(adjacency)
    child = list(adjacency) + [0]
    for vertex in range(order):
        if (signature >> vertex) & 1:
            child[order] |= 1 << vertex
        else:
            child[vertex] |= 1 << order
    return child


def reachable(adjacency: list[int], start: int) -> int:
    seen = 1 << start
    frontier = seen
    while frontier:
        bit = frontier & -frontier
        frontier -= bit
        vertex = bit.bit_length() - 1
        fresh = adjacency[vertex] & ~seen
        seen |= fresh
        frontier |= fresh
    return seen


def converse(adjacency: list[int]) -> list[int]:
    order = len(adjacency)
    reversed_adjacency = [0] * order
    for source in range(order):
        available = adjacency[source]
        while available:
            bit = available & -available
            available -= bit
            target = bit.bit_length() - 1
            reversed_adjacency[target] |= 1 << source
    return reversed_adjacency


def is_strong(adjacency: list[int]) -> bool:
    full = (1 << len(adjacency)) - 1
    return reachable(adjacency, 0) == full and reachable(converse(adjacency), 0) == full


def all_end_counts(adjacency: list[int]) -> list[list[int]]:
    order = len(adjacency)
    dp = [[0] * order for _ in range(1 << order)]
    for vertex in range(order):
        dp[1 << vertex][vertex] = 1
    for mask in range(1, 1 << order):
        for vertex, count in enumerate(dp[mask]):
            if count == 0:
                continue
            available = adjacency[vertex] & ~mask
            while available:
                bit = available & -available
                available -= bit
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += count
    return dp


def boundary_state(adjacency: list[int]) -> tuple[list[int], list[int], list[list[int]]]:
    order = len(adjacency)
    full = (1 << order) - 1
    end_dp = all_end_counts(adjacency)
    reverse_end_dp = all_end_counts(converse(adjacency))
    starts = list(reverse_end_dp[full])
    ends = list(end_dp[full])
    exposed = [[0] * order for _ in range(order)]
    for left_end in range(order):
        for right_start in range(order):
            if left_end == right_start:
                continue
            free = full ^ (1 << left_end) ^ (1 << right_start)
            subset = free
            total = 0
            while True:
                left = subset | (1 << left_end)
                right = full ^ left
                total += end_dp[left][left_end] * reverse_end_dp[right][right_start]
                if subset == 0:
                    break
                subset = (subset - 1) & free
            exposed[left_end][right_start] = total
    return starts, ends, exposed


def cut_field(
    starts: list[int], ends: list[int], exposed: list[list[int]]
) -> tuple[list[list[int]], list[int]]:
    order = len(starts)
    weights = [[0] * order for _ in range(order)]
    for left, right in combinations(range(order), 2):
        symmetric = exposed[left][right] + exposed[right][left]
        require(symmetric % 2 == 0, "integral symmetric cut weight")
        weights[left][right] = weights[right][left] = symmetric // 2
    field: list[int] = []
    for vertex in range(order):
        column = sum(exposed[other][vertex] for other in range(order))
        row = sum(exposed[vertex])
        require((column - row) % 2 == 0, "integral orientation field")
        field.append(starts[vertex] - ends[vertex] + (column - row) // 2)
    require(sum(field) == 0, "zero-sum orientation field")
    return weights, field


def insertion_count(
    state: tuple[list[int], list[int], list[list[int]]], signature: int
) -> int:
    starts, ends, exposed = state
    order = len(starts)
    total = 0
    for vertex in range(order):
        total += starts[vertex] if (signature >> vertex) & 1 else ends[vertex]
    for left_end in range(order):
        if (signature >> left_end) & 1:
            continue
        for right_start in range(order):
            if (signature >> right_start) & 1:
                total += exposed[left_end][right_start]
    return total


def defect_layers(adjacency: list[int], words: tuple[tuple[int, ...], ...]) -> tuple[int, int]:
    zero = 0
    one = 0
    for word in words:
        defects = sum(
            1
            for index in range(len(word) - 1)
            if not ((adjacency[word[index]] >> word[index + 1]) & 1)
        )
        zero += int(defects == 0)
        one += int(defects == 1)
    return zero, one


def profile(
    adjacency: list[int],
    parent_words: tuple[tuple[int, ...], ...] | None = None,
    check_children: bool = False,
) -> dict[str, object]:
    order = len(adjacency)
    require(order >= 2, "profile order at least two")
    state = boundary_state(adjacency)
    parent_h = sum(state[1])
    require(parent_h == sum(state[0]), "start/end Hamilton totals")
    weights, field = cut_field(*state)
    edge_pairs = tuple(combinations(range(order), 2))
    total_weight = sum(weights[left][right] for left, right in edge_pairs)
    defect_one = 2 * total_weight - (order - 1) * parent_h
    require(defect_one >= 0, "nonnegative one-defect count")
    if parent_words is not None:
        literal_h, literal_f1 = defect_layers(adjacency, parent_words)
        require((literal_h, literal_f1) == (parent_h, defect_one),
                "literal zero/one-defect layers")

    values: list[int] = []
    walsh_values: list[Fraction] = []
    for signature in range(1 << order):
        value = insertion_count(state, signature)
        values.append(value)
        cut = sum(
            weights[left][right]
            for left, right in edge_pairs
            if ((signature >> left) & 1) != ((signature >> right) & 1)
        )
        charge = sum(
            field[vertex]
            for vertex in range(order)
            if (signature >> vertex) & 1
        )
        require(value == parent_h + cut + charge, "cut-field identity")
        epsilon = [1 if (signature >> vertex) & 1 else -1 for vertex in range(order)]
        walsh = (
            Fraction(parent_h)
            + Fraction(total_weight, 2)
            + Fraction(sum(field[vertex] * epsilon[vertex] for vertex in range(order)), 2)
            - Fraction(
                sum(
                    weights[left][right] * epsilon[left] * epsilon[right]
                    for left, right in edge_pairs
                ),
                2,
            )
        )
        require(walsh.denominator == 1 and walsh.numerator == value,
                "degree-two Walsh identity")
        walsh_values.append(walsh)

    require(values[0] == values[-1] == parent_h, "constant source/sink cuts")
    require(all(value >= parent_h for value in values), "ear insertion monotonicity")
    require(all(value % 2 for value in values), "Redei parity")
    if check_children and is_strong(adjacency):
        for signature in range(1, (1 << order) - 1):
            require(is_strong(add_cut_ear(adjacency, signature)),
                    "strong nonconstant ear inheritance")

    mean = Fraction(sum(values), len(values))
    expected_mean = Fraction((order + 3) * parent_h + defect_one, 4)
    require(mean == Fraction(parent_h) + Fraction(total_weight, 2) == expected_mean,
            "uniform mean identity")
    variance = sum((Fraction(value) - mean) ** 2 for value in values) / len(values)
    walsh_variance = Fraction(
        sum(entry * entry for entry in field)
        + sum(weights[left][right] ** 2 for left, right in edge_pairs),
        4,
    )
    require(variance == walsh_variance, "Parseval variance identity")
    require(mean > parent_h, "positive mean increment")

    maximum = max(values)
    nonconstant_maximum = max(values[1:-1])
    require(maximum == nonconstant_maximum, "maximum occurs on a nonconstant cut")
    maximum_neighbor_drops = sorted({
        maximum - values[signature ^ (1 << vertex)]
        for signature, value in enumerate(values)
        if value == maximum
        for vertex in range(order)
        if values[signature ^ (1 << vertex)] < maximum
    })
    require(maximum_neighbor_drops, "maximum has a strictly lower cube neighbor")
    exact_floor = mean + variance / (mean - parent_h)
    require(maximum >= exact_floor, "support-sensitive second-moment floor")
    edge_count = order * (order - 1) // 2
    cauchy_floor = mean + (mean - parent_h) / edge_count
    require(variance >= (mean - parent_h) ** 2 / edge_count,
            "cut-weight Cauchy floor")
    require(exact_floor >= cauchy_floor, "exact floor dominates Cauchy floor")
    full_floor = (
        Fraction((order + 1) * (order + 2), 4 * order) * parent_h
        + Fraction(edge_count + 1, 4 * edge_count) * defect_one
    )
    universal_floor = Fraction((order + 1) * (order + 2), 4 * order) * parent_h
    require(cauchy_floor == full_floor >= universal_floor,
            "closed full and universal floors")
    require(maximum >= odd_ceiling(exact_floor) >= odd_ceiling(universal_floor),
            "odd-ceiling refinement")

    return {
        "order": order,
        "H": parent_h,
        "F1": defect_one,
        "W": total_weight,
        "field_square": sum(entry * entry for entry in field),
        "weight_square": sum(weights[left][right] ** 2 for left, right in edge_pairs),
        "mean": fraction_row(mean),
        "variance": fraction_row(variance),
        "exact_floor": fraction_row(exact_floor),
        "exact_oddceil": odd_ceiling(exact_floor),
        "universal_floor": fraction_row(universal_floor),
        "universal_oddceil": odd_ceiling(universal_floor),
        "maximum": maximum,
        "maximum_neighbor_drops": maximum_neighbor_drops,
        "minimum_maximum_neighbor_drop": maximum_neighbor_drops[0],
        "image": sorted(set(values[1:-1])),
        "strong": is_strong(adjacency),
    }


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    # The n=1 boundary has zero denominator and is deliberately outside the theorem.
    one_vertex_values = [1, 1]
    require(one_vertex_values[0] == one_vertex_values[1], "order-one zero-variance boundary")

    words = {
        order: tuple(permutations(range(order)))
        for order in range(2, 6)
    }
    order_rows: list[dict[str, object]] = []
    total_tournaments = 0
    total_cuts = 0
    total_strong_parents = 0
    total_strong_ears = 0
    triple_images: dict[tuple[int, int, Fraction], tuple[int, ...]] = {}
    triple_image_conflicts = 0

    selected_profiles: dict[int, dict[str, object]] = {}
    for order in range(2, 6):
        tournament_count = 1 << (order * (order - 1) // 2)
        strong_parents = 0
        strong_ears = 0
        strong_without_unit_maximum_descent = 0
        minimum_exact_slack: int | None = None
        minimum_universal_slack: int | None = None
        for code in range(tournament_count):
            adjacency = tournament_from_code(order, code)
            row = profile(adjacency, words[order], check_children=True)
            total_tournaments += 1
            total_cuts += 1 << order
            if row["strong"]:
                strong_parents += 1
                strong_ears += (1 << order) - 2
                if order >= 4:
                    strong_f1_floor = order - 1 if order % 2 == 0 else 2
                    require(
                        row["F1"] >= strong_f1_floor,
                        "strong cyclic-order additive F1 floor",
                    )
                if row["minimum_maximum_neighbor_drop"] != 2:
                    strong_without_unit_maximum_descent += 1
                key = (int(row["H"]), int(row["F1"]), Fraction(*row["variance"]))
                image = tuple(int(value) for value in row["image"])
                previous = triple_images.get(key)
                if previous is not None and previous != image:
                    triple_image_conflicts += 1
                triple_images[key] = image
                exact_slack = int(row["maximum"]) - int(row["exact_oddceil"])
                universal_slack = int(row["maximum"]) - int(row["universal_oddceil"])
                minimum_exact_slack = (
                    exact_slack if minimum_exact_slack is None
                    else min(minimum_exact_slack, exact_slack)
                )
                minimum_universal_slack = (
                    universal_slack if minimum_universal_slack is None
                    else min(minimum_universal_slack, universal_slack)
                )
            if (order, code) in {(2, 0), (3, 2), (5, 1015), (5, 759)}:
                selected_profiles[code if order != 2 else -1] = row
        total_strong_parents += strong_parents
        total_strong_ears += strong_ears
        order_rows.append({
            "order": order,
            "tournaments": tournament_count,
            "strong_parents": strong_parents,
            "cuts": tournament_count * (1 << order),
            "strong_ears": strong_ears,
            "strong_without_maximum_minus_two_neighbor": (
                strong_without_unit_maximum_descent
            ),
            "minimum_exact_oddceil_slack": minimum_exact_slack,
            "minimum_universal_oddceil_slack": minimum_universal_slack,
        })

    require(triple_image_conflicts == 0, "no same-triple image conflict through order five")
    order_two = selected_profiles[-1]
    require(
        (order_two["H"], order_two["F1"], order_two["mean"],
         order_two["variance"], order_two["maximum"])
        == (1, 1, [3, 2], [3, 4], 3),
        "order-two factor control",
    )
    cycle = selected_profiles[2]
    require(
        (cycle["H"], cycle["F1"], cycle["W"], cycle["mean"],
         cycle["variance"], cycle["exact_floor"], cycle["maximum"])
        == (3, 0, 3, [9, 2], [3, 4], [5, 1], 5),
        "C3 sharp equality",
    )
    hostile_a = selected_profiles[1015]
    hostile_b = selected_profiles[759]
    require(
        (hostile_a["H"], hostile_a["F1"], hostile_a["mean"])
        == (hostile_b["H"], hostile_b["F1"], hostile_b["mean"])
        == (9, 30, [51, 2]),
        "equal mean hostile",
    )
    require(
        hostile_a["variance"] == [305, 4]
        and hostile_b["variance"] == [315, 4]
        and hostile_a["image"] != hostile_b["image"],
        "variance separates the known mean hostile",
    )
    require(
        hostile_a["maximum_neighbor_drops"] == [4, 8, 18]
        and hostile_b["maximum_neighbor_drops"] == [6, 10, 12, 14],
        "frozen maximum-neighbor local-step hostiles",
    )
    require(
        [row["strong_without_maximum_minus_two_neighbor"] for row in order_rows]
        == [0, 0, 0, 400],
        "local maximum-minus-two hostile census",
    )

    large_rows = []
    for order, code in (
        (9, 68_164_491_031),
        (10, 10_852_575_419_951),
    ):
        row = profile(tournament_from_code(order, code))
        large_rows.append({"code": code, **row})

    expected_large = [
        {
            "order": 9,
            "code": 68_164_491_031,
            "H": 3_357,
            "F1": 16_434,
            "W": 21_645,
            "variance": [13_030_011, 4],
            "exact_floor": [34_825_587, 2_405],
            "exact_oddceil": 14_481,
            "universal_oddceil": 10_259,
            "maximum": 15_487,
        },
        {
            "order": 10,
            "code": 10_852_575_419_951,
            "H": 15_621,
            "F1": 95_625,
            "W": 118_107,
            "variance": [366_134_213, 4],
            "exact_floor": [9_002_648_278, 118_107],
            "exact_oddceil": 76_225,
            "universal_oddceil": 51_551,
            "maximum": 93_751,
        },
    ]
    for observed, expected in zip(large_rows, expected_large):
        require(
            all(observed[key] == value for key, value in expected.items()),
            f"selected-bank exact row order {expected['order']}",
        )

    # Algebraic iteration of the universal coefficient.
    iteration_rows = []
    for start, stop in ((3, 10), (9, 11), (10, 12), (5, 30)):
        product = Fraction(1)
        for order in range(start, stop):
            product *= Fraction((order + 1) * (order + 2), 4 * order)
        closed = Fraction(
            stop,
            start,
        )
        for value in range(start + 2, stop + 2):
            closed *= value
        closed /= 4 ** (stop - start)
        require(product == closed, "recursive product closed form")
        iteration_rows.append([start, stop, product.numerator, product.denominator])

    ledger = {
        "order_rows": order_rows,
        "totals": [total_tournaments, total_cuts, total_strong_parents, total_strong_ears],
        "order_one_boundary": one_vertex_values,
        "order_two_control": order_two,
        "c3_equality": cycle,
        "equal_mean_variance_hostiles": [hostile_a, hostile_b],
        "same_triple_image_conflicts_through_order_five": triple_image_conflicts,
        "large_rows": large_rows,
        "iteration_rows": iteration_rows,
        "formula": "Var=1/4*(sum_i h_i^2+sum_ij w_ij^2)",
        "floor": "M>=mu+Var/(mu-H)>=((n+1)(n+2)/(4n))*H",
        "scope": "maxima and variance only; interval overlap remains OPEN",
    }
    semantic = canonical_digest(ledger)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic ledger digest")

    print("THM-4115 EAR-CUT WALSH VARIANCE PRIMARY REFEREE")
    print("status=PRIMARY_PASS")
    print("walsh=X_S=H+W/2+sum_i(h_i eps_i)/2-sum_ij(w_ij eps_i eps_j)/2")
    print("variance=1/4*(sum_i h_i^2+sum_ij w_ij^2)")
    print("exact_floor=M>=mu+Var/(mu-H)")
    print("universal_floor=M>=((n+1)(n+2)/(4n))*H")
    print("order_rows=", order_rows)
    print("totals_tournaments_cuts_strong_parents_strong_ears=", ledger["totals"])
    print("order_one_denominator_hostile=", one_vertex_values)
    print("order_two_H_F1_mean_var_max=", [
        order_two["H"], order_two["F1"], order_two["mean"],
        order_two["variance"], order_two["maximum"],
    ])
    print("C3_H_F1_W_mean_var_floor_max=", [
        cycle["H"], cycle["F1"], cycle["W"], cycle["mean"],
        cycle["variance"], cycle["exact_floor"], cycle["maximum"],
    ])
    print("equal_mean_codes_H_F1_mean_variances=", [
        [1015, hostile_a["H"], hostile_a["F1"], hostile_a["mean"],
         hostile_a["variance"], hostile_a["maximum_neighbor_drops"]],
        [759, hostile_b["H"], hostile_b["F1"], hostile_b["mean"],
         hostile_b["variance"], hostile_b["maximum_neighbor_drops"]],
    ])
    print("same_H_F1_variance_different_image_conflicts_through_order5=", triple_image_conflicts)
    print("selected_large_rows=", [
        {
            key: row[key]
            for key in (
                "order", "code", "H", "F1", "W", "variance", "exact_floor",
                "exact_oddceil", "universal_oddceil", "maximum",
            )
        }
        for row in large_rows
    ])
    print("recursive_product_rows=", iteration_rows)
    print("semantic_sha256=", semantic)
    print("semantic_preserved=full cut fibre, Walsh degrees 0/1/2, and lower support H")
    print("semantic_destroyed=individual cut image after moments")
    print("semantic_scope=maxima and variance only; interval overlap remains OPEN")
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
