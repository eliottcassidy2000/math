#!/usr/bin/env python3
"""Independent literal referee for the reserved THM-4115 candidate.

The program imports no tournament compiler and no THM-4097/THM-4111 helper.
For every labelled tournament of orders two through five it literally scans
all parent permutations and all permutations of every labelled cut-ear child.
It then reconstructs the quadratic Boolean cut polynomial solely from its
singleton values and second Boolean differences.

The exact gates cover the degree-two Walsh identity, its mean and variance,
the lower-support second-moment maximum bound, the resulting universal and
F1-sensitive coefficients, Redei parity/odd ceilings, strong-ear inheritance,
and the closed recursive product.  Frozen hostile controls delimit the scope:
order one has a zero denominator, and equal (H,F1,mean) need not determine the
variance or cut image.  Only maxima are controlled; interval overlap is open.

Every executable gate uses ``require`` and survives ``python -O``.  Output is
forced to LF so the normal and optimized streams are byte-stable on Windows.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
from math import factorial


EXPECTED_SEMANTIC_SHA256 = "d805f598ce72128f22bd46bf4c8ba696b471be4934455af00b02837fa508cc0c"


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def fraction_row(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def tournament_from_code(order: int, code: int) -> tuple[int, ...]:
    """Decode LSB-first upper-pair orientations into outgoing bit masks."""
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
    """Adjoin x, with x->v exactly when the bit of v belongs to cut."""
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
        defects += 1 - ((adjacency[left] >> right) & 1)
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


def ear_profile(
    adjacency: tuple[int, ...], child_words: tuple[tuple[int, ...], ...]
) -> list[int]:
    return [
        literal_hamilton_count(add_cut_ear(adjacency, cut), child_words)
        for cut in range(1 << len(adjacency))
    ]


def maximum_neighbor_drops(values: list[int], order: int) -> list[int]:
    maximum = max(values)
    drops = sorted({
        maximum - values[cut ^ (1 << vertex)]
        for cut, value in enumerate(values)
        if value == maximum
        for vertex in range(order)
        if values[cut ^ (1 << vertex)] < maximum
    })
    require(bool(drops), "maximum has a strictly lower cube neighbor")
    return drops


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


def reconstruct_boolean_cut_field(
    values: list[int], order: int
) -> tuple[int, list[Fraction], dict[tuple[int, int], Fraction]]:
    """Recover H, h_i, w_ij from f(0), f(e_i), and second differences."""
    base = values[0]
    weights: dict[tuple[int, int], Fraction] = {}
    for left in range(order):
        for right in range(left + 1, order):
            second_difference = (
                values[(1 << left) | (1 << right)]
                - values[1 << left]
                - values[1 << right]
                + base
            )
            weights[left, right] = Fraction(-second_difference, 2)
            require(
                weights[left, right].denominator == 1,
                "integral pair weight from Boolean second difference",
            )

    linear: list[Fraction] = []
    for vertex in range(order):
        incident = sum(
            (
                weight
                for pair, weight in weights.items()
                if vertex in pair
            ),
            Fraction(0),
        )
        linear.append(Fraction(values[1 << vertex] - base) - incident)
        require(
            linear[-1].denominator == 1,
            "integral linear field from singleton difference",
        )
    require(sum(linear, Fraction(0)) == 0, "zero-sum linear cut field")
    return base, linear, weights


def cut_polynomial_value(
    base: int,
    linear: list[Fraction],
    weights: dict[tuple[int, int], Fraction],
    cut: int,
) -> Fraction:
    order = len(linear)
    value = Fraction(base)
    for vertex in range(order):
        if (cut >> vertex) & 1:
            value += linear[vertex]
    for (left, right), weight in weights.items():
        if ((cut >> left) & 1) != ((cut >> right) & 1):
            value += weight
    return value


def walsh_polynomial_value(
    mean: Fraction,
    linear: list[Fraction],
    weights: dict[tuple[int, int], Fraction],
    cut: int,
) -> Fraction:
    signs = [1 - 2 * ((cut >> vertex) & 1) for vertex in range(len(linear))]
    value = mean
    value -= Fraction(1, 2) * sum(
        (linear[vertex] * signs[vertex] for vertex in range(len(linear))),
        Fraction(0),
    )
    value -= Fraction(1, 2) * sum(
        (
            weight * signs[left] * signs[right]
            for (left, right), weight in weights.items()
        ),
        Fraction(0),
    )
    return value


def audit_profile(
    order: int,
    parent_h: int,
    defect_one: int,
    values: list[int],
    label: str,
) -> dict[str, object]:
    base, linear, weights = reconstruct_boolean_cut_field(values, order)
    require(base == parent_h, f"{label}: base field equals parent H")

    pair_count = order * (order - 1) // 2
    require(len(weights) == pair_count, f"{label}: complete pair field")
    weight_sum = sum(weights.values(), Fraction(0))
    mean = Fraction(sum(values), 1 << order)
    predicted_mean = Fraction((order + 3) * parent_h + defect_one, 4)
    require(mean == predicted_mean, f"{label}: exact all-cut mean")
    require(
        mean == Fraction(parent_h) + weight_sum / 2,
        f"{label}: Walsh constant coefficient",
    )
    require(
        weight_sum == Fraction((order - 1) * parent_h + defect_one, 2),
        f"{label}: full F1 pair-weight coefficient",
    )

    for cut, observed in enumerate(values):
        require(
            cut_polynomial_value(base, linear, weights, cut) == observed,
            f"{label}: exact degree-two Boolean cut identity",
        )
        require(
            walsh_polynomial_value(mean, linear, weights, cut) == observed,
            f"{label}: exact degree-two Walsh identity",
        )

    variance = sum(
        ((Fraction(value) - mean) ** 2 for value in values), Fraction(0)
    ) / (1 << order)
    walsh_variance = Fraction(1, 4) * (
        sum((entry * entry for entry in linear), Fraction(0))
        + sum((entry * entry for entry in weights.values()), Fraction(0))
    )
    require(variance == walsh_variance, f"{label}: exact Walsh variance")
    require(all(value >= parent_h for value in values), f"{label}: X at least H")
    require(all(value % 2 for value in values), f"{label}: Redei odd parity")

    delta = mean - parent_h
    require(delta > 0, f"{label}: positive second-moment denominator")
    maximum = max(values)
    local_drops = maximum_neighbor_drops(values, order)
    support_rhs = mean + variance / delta
    require(maximum >= support_rhs, f"{label}: lower-support variance bound")
    require(
        maximum >= odd_ceiling(support_rhs),
        f"{label}: odd ceiling of variance bound",
    )

    require(
        sum((entry * entry for entry in weights.values()), Fraction(0))
        >= weight_sum * weight_sum / pair_count,
        f"{label}: pair-weight Cauchy bound",
    )
    cauchy_rhs = mean + delta / pair_count
    full_f1_rhs = (
        Fraction((order + 1) * (order + 2), 4 * order) * parent_h
        + Fraction(order * (order - 1) + 2, 4 * order * (order - 1))
        * defect_one
    )
    universal_rhs = (
        Fraction((order + 1) * (order + 2), 4 * order) * parent_h
    )
    require(cauchy_rhs == full_f1_rhs, f"{label}: full F1 coefficient algebra")
    require(support_rhs >= cauchy_rhs, f"{label}: variance refines Cauchy bound")
    require(full_f1_rhs >= universal_rhs, f"{label}: discard nonnegative F1")
    require(maximum >= full_f1_rhs, f"{label}: full F1 maximum bound")
    require(maximum >= universal_rhs, f"{label}: universal maximum bound")
    require(
        maximum >= odd_ceiling(full_f1_rhs),
        f"{label}: full F1 odd ceiling",
    )
    require(
        maximum >= odd_ceiling(universal_rhs),
        f"{label}: universal odd ceiling",
    )

    return {
        "mean": fraction_row(mean),
        "variance": fraction_row(variance),
        "maximum": maximum,
        "maximum_neighbor_drops": local_drops,
        "minimum_maximum_neighbor_drop": local_drops[0],
        "support_rhs": fraction_row(support_rhs),
        "full_f1_rhs": fraction_row(full_f1_rhs),
        "universal_rhs": fraction_row(universal_rhs),
        "weight_sum": fraction_row(weight_sum),
        "weight_square_sum": fraction_row(
            sum((entry * entry for entry in weights.values()), Fraction(0))
        ),
        "linear_square_sum": fraction_row(
            sum((entry * entry for entry in linear), Fraction(0))
        ),
    }


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    words = {
        order: tuple(permutations(range(order)))
        for order in range(1, 7)
    }

    # The support-sensitive quotient is deliberately not extended to n=1:
    # both cuts give H=1, so mu-H=0 and Var=0.
    singleton = tournament_from_code(1, 0)
    singleton_h, singleton_f1 = zero_one_layers(singleton, words[1])
    singleton_values = ear_profile(singleton, words[2])
    singleton_mean = Fraction(sum(singleton_values), 2)
    singleton_variance = sum(
        ((Fraction(value) - singleton_mean) ** 2 for value in singleton_values),
        Fraction(0),
    ) / 2
    require(
        [singleton_h, singleton_f1, singleton_values]
        == [1, 0, [1, 1]],
        "n=1 denominator hostile profile",
    )
    require(
        singleton_mean - singleton_h == 0 and singleton_variance == 0,
        "n=1 zero-over-zero denominator boundary",
    )
    denominator_hostile = {
        "order": 1,
        "H": singleton_h,
        "F1": singleton_f1,
        "values": singleton_values,
        "mean": fraction_row(singleton_mean),
        "variance": fraction_row(singleton_variance),
        "mu_minus_H": [0, 1],
    }

    order_rows: list[dict[str, object]] = []
    total_parents = 0
    total_ears = 0
    total_parent_words = 0
    total_child_words = 0
    total_strong_ears = 0
    expected_counts = {
        2: (2, 0),
        3: (8, 2),
        4: (64, 24),
        5: (1_024, 544),
    }

    factor_control: dict[str, object] | None = None
    cycle_control: dict[str, object] | None = None
    hostile_profiles: dict[int, dict[str, object]] = {}

    for order in range(2, 6):
        tournament_count = 1 << (order * (order - 1) // 2)
        parent_words = words[order]
        child_words = words[order + 1]
        strong_parents = 0
        strong_ear_checks = 0
        strong_without_unit_maximum_descent = 0
        minimum_support_slack: tuple[Fraction, int] | None = None
        minimum_universal_slack: tuple[Fraction, int] | None = None
        variance_range: tuple[Fraction, Fraction] | None = None

        for code in range(tournament_count):
            adjacency = tournament_from_code(order, code)
            parent_h, defect_one = zero_one_layers(adjacency, parent_words)
            values = ear_profile(adjacency, child_words)
            report = audit_profile(
                order,
                parent_h,
                defect_one,
                values,
                f"n={order},code={code}",
            )

            variance = Fraction(*report["variance"])
            support_rhs = Fraction(*report["support_rhs"])
            universal_rhs = Fraction(*report["universal_rhs"])
            support_candidate = (Fraction(max(values)) - support_rhs, code)
            universal_candidate = (Fraction(max(values)) - universal_rhs, code)
            if minimum_support_slack is None or support_candidate < minimum_support_slack:
                minimum_support_slack = support_candidate
            if minimum_universal_slack is None or universal_candidate < minimum_universal_slack:
                minimum_universal_slack = universal_candidate
            if variance_range is None:
                variance_range = (variance, variance)
            else:
                variance_range = (
                    min(variance_range[0], variance),
                    max(variance_range[1], variance),
                )

            if order == 2 and code == 0:
                factor = Fraction((order + 1) * (order + 2), 4 * order)
                require(factor == Fraction(3, 2), "n=2 recursive factor")
                factor_control = {
                    "code": code,
                    "H": parent_h,
                    "F1": defect_one,
                    "values": values,
                    "factor": fraction_row(factor),
                    "audit": report,
                }

            if order == 3 and code == 2:
                require(is_strong(adjacency), "C3 strong equality parent")
                require(
                    [parent_h, defect_one, values]
                    == [3, 0, [3, 5, 5, 5, 5, 5, 5, 3]],
                    "C3 literal sharp profile",
                )
                require(
                    report["variance"] == [3, 4]
                    and report["support_rhs"] == [5, 1]
                    and report["universal_rhs"] == [5, 1]
                    and report["maximum"] == 5,
                    "C3 sharp equality in both bounds",
                )
                cycle_control = {
                    "code": code,
                    "H": parent_h,
                    "F1": defect_one,
                    "values": values,
                    "audit": report,
                }

            if order == 5 and code in (1_015, 759):
                hostile_profiles[code] = {
                    "code": code,
                    "H": parent_h,
                    "F1": defect_one,
                    "mean": report["mean"],
                    "variance": report["variance"],
                    "image": sorted(set(values)),
                    "nonconstant_image": sorted(set(values[1:-1])),
                    "maximum": max(values),
                    "maximum_neighbor_drops": report["maximum_neighbor_drops"],
                }

            if is_strong(adjacency):
                strong_parents += 1
                if order >= 4:
                    strong_f1_floor = order - 1 if order % 2 == 0 else 2
                    require(
                        defect_one >= strong_f1_floor,
                        "literal strong cyclic-order additive F1 floor",
                    )
                if report["minimum_maximum_neighbor_drop"] != 2:
                    strong_without_unit_maximum_descent += 1
                for cut in range(1, (1 << order) - 1):
                    require(
                        is_strong(add_cut_ear(adjacency, cut)),
                        "literal strong-ear inheritance",
                    )
                    strong_ear_checks += 1

        require(minimum_support_slack is not None, "nonempty support-slack row")
        require(minimum_universal_slack is not None, "nonempty universal-slack row")
        require(variance_range is not None, "nonempty variance range")
        require(
            (tournament_count, strong_parents) == expected_counts[order],
            "labelled tournament and strong-parent counts",
        )

        order_rows.append(
            {
                "order": order,
                "tournaments": tournament_count,
                "strong_parents": strong_parents,
                "ear_checks": tournament_count * (1 << order),
                "parent_words_scanned": tournament_count * len(parent_words),
                "child_words_scanned": tournament_count * (1 << order) * len(child_words),
                "strong_ear_checks": strong_ear_checks,
                "strong_without_maximum_minus_two_neighbor": (
                    strong_without_unit_maximum_descent
                ),
                "variance_range": [
                    fraction_row(variance_range[0]),
                    fraction_row(variance_range[1]),
                ],
                "minimum_support_slack_code": [
                    *fraction_row(minimum_support_slack[0]),
                    minimum_support_slack[1],
                ],
                "minimum_universal_slack_code": [
                    *fraction_row(minimum_universal_slack[0]),
                    minimum_universal_slack[1],
                ],
            }
        )
        total_parents += tournament_count
        total_ears += tournament_count * (1 << order)
        total_parent_words += tournament_count * len(parent_words)
        total_child_words += tournament_count * (1 << order) * len(child_words)
        total_strong_ears += strong_ear_checks

    totals = [
        total_parents,
        total_ears,
        total_parent_words,
        total_child_words,
        total_strong_ears,
    ]
    require(
        totals == [1_098, 33_864, 124_468, 23_717_424, 16_668],
        "aggregate literal census counts",
    )
    require(factor_control is not None, "n=2 factor control captured")
    require(cycle_control is not None, "C3 sharp equality control captured")

    require(set(hostile_profiles) == {759, 1_015}, "variance hostiles captured")
    first = hostile_profiles[1_015]
    second = hostile_profiles[759]
    require(
        [first["H"], first["F1"], first["mean"]]
        == [second["H"], second["F1"], second["mean"]]
        == [9, 30, [51, 2]],
        "hostiles share H, F1, and all-cut mean",
    )
    require(
        first["variance"] == [305, 4]
        and second["variance"] == [315, 4],
        "hostiles have frozen distinct variances",
    )
    expected_images = {
        1_015: [15, 17, 19, 23, 25, 27, 29, 33, 37, 41],
        759: [15, 17, 19, 23, 29, 31, 33, 37, 43],
    }
    require(
        first["nonconstant_image"] == expected_images[1_015]
        and second["nonconstant_image"] == expected_images[759],
        "hostiles have frozen distinct nonconstant images",
    )
    require(
        first["nonconstant_image"] != second["nonconstant_image"],
        "equal mean does not determine cut image",
    )
    require(
        first["maximum_neighbor_drops"] == [4, 8, 18]
        and second["maximum_neighbor_drops"] == [6, 10, 12, 14],
        "literal maximum-neighbor local-step hostiles",
    )
    require(
        [row["strong_without_maximum_minus_two_neighbor"] for row in order_rows]
        == [0, 0, 0, 400],
        "literal local maximum-minus-two hostile census",
    )

    recursive_rows: list[dict[str, object]] = []
    for initial_order in range(2, 10):
        for final_order in range(initial_order + 1, initial_order + 9):
            product = Fraction(1)
            for order in range(initial_order, final_order):
                product *= Fraction((order + 1) * (order + 2), 4 * order)
            closed = Fraction(
                final_order * factorial(final_order + 1),
                initial_order
                * factorial(initial_order + 1)
                * 4 ** (final_order - initial_order),
            )
            require(product == closed, "recursive sharp-factor product algebra")
        terminal = initial_order + 8
        terminal_closed = Fraction(
            terminal * factorial(terminal + 1),
            initial_order * factorial(initial_order + 1) * 4**8,
        )
        recursive_rows.append(
            {
                "initial_order": initial_order,
                "final_order": terminal,
                "eight_step_product": fraction_row(terminal_closed),
            }
        )

    ledger = {
        "denominator_hostile": denominator_hostile,
        "orders": order_rows,
        "totals": totals,
        "n2_factor_control": factor_control,
        "c3_sharp_control": cycle_control,
        "equal_mean_variance_hostiles": [first, second],
        "recursive_rows": recursive_rows,
        "walsh_formula": "f=mu-(1/2)sum_i h_i r_i-(1/2)sum_ij w_ij r_i*r_j",
        "variance_formula": "Var=(sum_i h_i^2+sum_ij w_ij^2)/4",
        "support_bound": "max>=mu+Var/(mu-H)",
        "full_f1_bound": (
            "max>=((n+1)(n+2)/(4n))H"
            "+((n(n-1)+2)/(4n(n-1)))F1"
        ),
        "scope": "image maxima and recursive growth only; interval overlap remains OPEN",
    }
    semantic = canonical_digest(ledger)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic ledger digest")

    print("THM-4115 EAR-CUT WALSH VARIANCE INDEPENDENT AUDIT")
    print("status=INDEPENDENT_AUDIT_PASS")
    print("method=literal parent and child permutations; Boolean singleton/second differences")
    print("degree2_walsh=f=mu-(1/2)sum_i h_i*r_i-(1/2)sum_ij w_ij*r_i*r_j")
    print("variance=Var(X)=(sum_i h_i^2+sum_ij w_ij^2)/4")
    print("support_bound=max(X)>=mu+Var(X)/(mu-H)")
    print("full_F1_bound=max>=((n+1)(n+2)/(4n))H+((n(n-1)+2)/(4n(n-1)))F1")
    print("universal_bound=max>=((n+1)(n+2)/(4n))H")
    print("n1_denominator_hostile=", denominator_hostile)
    print("exhaustive_order_rows=", order_rows)
    print("totals_parents_ears_parentwords_childwords_strongears=", totals)
    print("n2_factor_control=", factor_control)
    print("C3_sharp_equality=", cycle_control)
    print("equal_mean_distinct_variance_hostiles=", [first, second])
    print("recursive_product_rows=", recursive_rows)
    print("recursive_product=M_m/M_n>=m*(m+1)!/(n*(n+1)!*4^(m-n))")
    print("semantic_sha256=", semantic)
    print("semantic_preserved=literal path counts, full cut field, Walsh norm, lower support")
    print("semantic_destroyed=interval placement and overlap after taking a maximum")
    print("semantic_scope=maxima/variance only; interval overlap remains OPEN")
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
