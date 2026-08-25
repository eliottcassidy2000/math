#!/usr/bin/env python3
"""Independent literal audit of THM-4111's uniform cut-ear average.

This implementation imports no tournament compiler and does not use the
primary subset DP.  It decodes every labelled tournament through order five,
enumerates its vertex orderings to obtain the zero- and one-defect layers,
then enumerates every ordering of every cut-ear child literally.

The audit also freezes a distributional hostile: two strong order-five
tournaments have the same ``(H,F1)`` and hence the same uniform ear mean, but
different nonconstant ear images.  Thus the average proves growth of maxima,
not interval propagation.

Every executable gate uses ``require`` and survives ``python -O``.
"""

from __future__ import annotations

import json
from fractions import Fraction
from hashlib import sha256
from itertools import permutations


EXPECTED_SEMANTIC_SHA256 = "1580575c67bc5b8cdea70db9c7efef256cad50b8a5ed3ce18e7a73303dfad9f0"


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def tournament_from_code(order: int, code: int) -> tuple[int, ...]:
    adjacency = [0] * order
    bit_index = 0
    for left in range(order):
        for right in range(left + 1, order):
            if (code >> bit_index) & 1:
                adjacency[left] |= 1 << right
            else:
                adjacency[right] |= 1 << left
            bit_index += 1
    return tuple(adjacency)


def add_cut_ear(adjacency: tuple[int, ...], cut: int) -> tuple[int, ...]:
    order = len(adjacency)
    child = list(adjacency) + [0]
    for vertex in range(order):
        if (cut >> vertex) & 1:
            child[order] |= 1 << vertex
        else:
            child[vertex] |= 1 << order
    return tuple(child)


def defect_count(adjacency: tuple[int, ...], word: tuple[int, ...]) -> int:
    defects = 0
    for left, right in zip(word, word[1:]):
        defects += not ((adjacency[left] >> right) & 1)
    return defects


def zero_one_layers(
    adjacency: tuple[int, ...], words: tuple[tuple[int, ...], ...]
) -> tuple[int, int]:
    zero = 0
    one = 0
    for word in words:
        defects = defect_count(adjacency, word)
        zero += defects == 0
        one += defects == 1
    return zero, one


def literal_hamilton_count(
    adjacency: tuple[int, ...], words: tuple[tuple[int, ...], ...]
) -> int:
    total = 0
    for word in words:
        total += defect_count(adjacency, word) == 0
    return total


def is_strong(adjacency: tuple[int, ...]) -> bool:
    order = len(adjacency)
    full = (1 << order) - 1
    for root in range(order):
        reached = 1 << root
        while True:
            enlarged = reached
            for vertex in range(order):
                if (reached >> vertex) & 1:
                    enlarged |= adjacency[vertex]
            if enlarged == reached:
                break
            reached = enlarged
        if reached != full:
            return False
    return True


def odd_ceiling(value: Fraction) -> int:
    ceiling = (value.numerator + value.denominator - 1) // value.denominator
    return ceiling if ceiling % 2 else ceiling + 1


def fraction_row(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def ear_profile(
    adjacency: tuple[int, ...], child_words: tuple[tuple[int, ...], ...]
) -> list[int]:
    return [
        literal_hamilton_count(add_cut_ear(adjacency, cut), child_words)
        for cut in range(1 << len(adjacency))
    ]


def main() -> None:
    words = {
        order: tuple(permutations(range(order)))
        for order in range(2, 7)
    }
    order_rows: list[dict[str, object]] = []
    total_parents = 0
    total_ears = 0
    total_parent_words = 0
    total_child_words = 0
    total_strong_ears = 0

    expected_rows = {
        2: (2, 0, [7, 4], [0, 1, 1, 3, 2, 1]),
        3: (8, 2, [1, 6], [2, 3, 0, 5, 5, 1]),
        4: (64, 24, [5, 4], [4, 5, 7, 15, 79, 7]),
        5: (1_024, 544, [13, 15], [40, 15, 18, 43, 179, 5]),
    }

    for order in range(2, 6):
        tournament_count = 1 << (order * (order - 1) // 2)
        parent_words = words[order]
        child_words = words[order + 1]
        strong_parents = 0
        strong_ear_checks = 0
        worst: tuple[Fraction, int, int, int, int, Fraction] | None = None

        for code in range(tournament_count):
            adjacency = tournament_from_code(order, code)
            parent_h, defect_one = zero_one_layers(adjacency, parent_words)
            ear_values = ear_profile(adjacency, child_words)
            nonconstant = ear_values[1:-1]

            predicted_sum = (1 << (order - 2)) * (
                (order + 3) * parent_h + defect_one
            )
            require(sum(ear_values) == predicted_sum,
                    "literal all-cut average identity")
            require(ear_values[0] == ear_values[-1] == parent_h,
                    "literal constant source/sink ears")
            require(all(value % 2 for value in ear_values),
                    "literal Redei parity")

            nonconstant_mean = Fraction(
                predicted_sum - 2 * parent_h, (1 << order) - 2
            )
            require(Fraction(sum(nonconstant), len(nonconstant)) == nonconstant_mean,
                    "literal nonconstant mean")
            require(
                max(nonconstant)
                >= odd_ceiling(Fraction(order + 3, 4) * parent_h),
                "literal selected-growth lower bound",
            )

            candidate = (
                Fraction(max(nonconstant), parent_h) - Fraction(order + 3, 4),
                code,
                parent_h,
                defect_one,
                max(nonconstant),
                nonconstant_mean,
            )
            if worst is None or candidate < worst:
                worst = candidate

            if is_strong(adjacency):
                strong_parents += 1
                for cut in range(1, (1 << order) - 1):
                    require(is_strong(add_cut_ear(adjacency, cut)),
                            "literal strong-ear inheritance")
                    strong_ear_checks += 1

        require(worst is not None, "nonempty labelled census")
        slack, code, parent_h, defect_one, maximum, mean = worst
        row = {
            "order": order,
            "tournaments": tournament_count,
            "strong_parents": strong_parents,
            "identity_checks": tournament_count,
            "ear_checks": tournament_count * (1 << order),
            "parent_words_scanned": tournament_count * len(parent_words),
            "child_words_scanned": tournament_count * (1 << order) * len(child_words),
            "strong_ear_checks": strong_ear_checks,
            "minimum_max_ratio_slack": fraction_row(slack),
            "minimum_row": [
                code,
                parent_h,
                defect_one,
                maximum,
                *fraction_row(mean),
            ],
        }
        expected_count, expected_strong, expected_slack, expected_minimum = (
            expected_rows[order]
        )
        require(tournament_count == expected_count, "labelled tournament count")
        require(strong_parents == expected_strong, "labelled strong-parent count")
        require(row["minimum_max_ratio_slack"] == expected_slack,
                "minimum growth-slack row")
        require(row["minimum_row"] == expected_minimum,
                "minimum growth witness row")
        order_rows.append(row)

        total_parents += tournament_count
        total_ears += tournament_count * (1 << order)
        total_parent_words += tournament_count * len(parent_words)
        total_child_words += tournament_count * (1 << order) * len(child_words)
        total_strong_ears += strong_ear_checks

    require(
        [total_parents, total_ears, total_parent_words, total_child_words,
         total_strong_ears]
        == [1_098, 33_864, 124_468, 23_717_424, 16_668],
        "aggregate literal census counts",
    )

    cycle = tournament_from_code(3, 2)
    cycle_h, cycle_f1 = zero_one_layers(cycle, words[3])
    cycle_ears = ear_profile(cycle, words[4])
    require(is_strong(cycle), "C3 strong control")
    require((cycle_h, cycle_f1) == (3, 0), "C3 zero/one-defect layers")
    require(cycle_ears == [3, 5, 5, 5, 5, 5, 5, 3],
            "C3 literal ear profile")
    require(max(cycle_ears[1:-1]) < Fraction(7, 4) * cycle_h,
            "C3 defeats the n+4 coefficient")

    hostile_rows: list[dict[str, object]] = []
    expected_images = {
        1_015: [15, 17, 19, 23, 25, 27, 29, 33, 37, 41],
        759: [15, 17, 19, 23, 29, 31, 33, 37, 43],
    }
    for code in (1_015, 759):
        adjacency = tournament_from_code(5, code)
        parent_h, defect_one = zero_one_layers(adjacency, words[5])
        ears = ear_profile(adjacency, words[6])
        image = sorted(set(ears[1:-1]))
        require(is_strong(adjacency), "equal-mean hostile parent is strong")
        require((parent_h, defect_one, sum(ears), sum(ears[1:-1]))
                == (9, 30, 816, 798),
                "equal-mean hostile scalar data")
        require(image == expected_images[code], "equal-mean hostile ear image")
        hostile_rows.append(
            {
                "code": code,
                "H": parent_h,
                "F1": defect_one,
                "all_cut_sum": sum(ears),
                "nonconstant_sum": sum(ears[1:-1]),
                "nonconstant_mean": fraction_row(Fraction(798, 30)),
                "image": image,
            }
        )
    first_image = set(hostile_rows[0]["image"])
    second_image = set(hostile_rows[1]["image"])
    require(first_image != second_image, "equal mean does not determine ear image")
    require(
        [sorted(first_image - second_image), sorted(second_image - first_image)]
        == [[25, 27, 41], [31, 43]],
        "equal-mean hostile symmetric difference",
    )

    selected_rows = []
    for order, parent_maximum, observed_child_maximum in (
        (9, 3_357, 15_621),
        (10, 15_621, 93_751),
    ):
        forced = odd_ceiling(Fraction(order + 3, 4) * parent_maximum)
        require(observed_child_maximum >= forced,
                "inherited selected-bank maximum control")
        selected_rows.append(
            [order, parent_maximum, forced, observed_child_maximum]
        )

    for initial_order in range(3, 9):
        product = Fraction(1)
        for order in range(initial_order, initial_order + 6):
            product *= Fraction(order + 3, 4)
        closed = Fraction(
            1,
            4**6,
        )
        for factor in range(initial_order + 3, initial_order + 9):
            closed *= factor
        require(product == closed, "six-step factorial product identity")

    ledger = {
        "orders": order_rows,
        "totals": [
            total_parents,
            total_ears,
            total_parent_words,
            total_child_words,
            total_strong_ears,
        ],
        "c3_hostile": [cycle_h, cycle_f1, cycle_ears],
        "equal_mean_hostiles": hostile_rows,
        "hostile_image_differences": [
            sorted(first_image - second_image),
            sorted(second_image - first_image),
        ],
        "selected_bank_rows": selected_rows,
        "formula": "sum_cut_H=2^(n-2)*((n+3)*H+F1)",
        "scope": "image maxima only; (H,F1) does not determine an ear image",
    }
    semantic = canonical_digest(ledger)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                "independent semantic ledger digest")

    print("THM-4111 UNIFORM CUT-EAR AVERAGE INDEPENDENT AUDIT")
    print("status=INDEPENDENT_AUDIT_PASS")
    print("method=literal parent permutations plus literal child permutations; no primary import")
    print("formula=sum_S H(T+x_S)=2^(n-2)*((n+3)H(T)+F1(T))")
    print("exhaustive_order_rows=", order_rows)
    print("totals_parents_ears_parentwords_childwords_strongears=", ledger["totals"])
    print("C3_hostile_H_F1_ears=", ledger["c3_hostile"])
    print("equal_H_F1_mean_hostiles=", hostile_rows)
    print("hostile_image_differences=", ledger["hostile_image_differences"])
    print("selected_bank_order_parent_forced_observed=", selected_rows)
    print("semantic_sha256=", semantic)
    print("semantic_preserved=literal cuts, path words, zero/one-defect layers, strongness")
    print("semantic_destroyed=individual cut distribution after taking the mean")
    print("semantic_scope=unbounded image maxima only; interval overlap remains OPEN")
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
