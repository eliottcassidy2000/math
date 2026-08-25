#!/usr/bin/env python3
"""Exact controls for THM-4111's uniform cut-ear average.

The proof is a double count, not a finite extrapolation.  This referee checks
every labelled tournament through order five by a subset Hamiltonian-path DP,
enumerates the zero/one-defect word layers independently, and tests every cut
ear.  It also freezes the directed 3-cycle boundary that rules out replacing
the coefficient ``n+3`` by ``n+4`` even for strong parents.

Every executable gate uses ``require`` and survives ``python -O``.
"""

from __future__ import annotations

import json
from fractions import Fraction
from hashlib import sha256
from itertools import permutations


EXPECTED_SEMANTIC_SHA256 = "0bb800d8caa1c1fd449657fb9f68a33842062ec4f3866f34488d9c8ea3251915"


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def tournament_from_code(order: int, code: int) -> list[int]:
    adjacency = [0] * order
    bit = 0
    for left in range(order):
        for right in range(left + 1, order):
            if (code >> bit) & 1:
                adjacency[left] |= 1 << right
            else:
                adjacency[right] |= 1 << left
            bit += 1
    return adjacency


def hamilton_count(adjacency: list[int]) -> int:
    order = len(adjacency)
    full = (1 << order) - 1
    paths = [[0] * order for _ in range(1 << order)]
    for vertex in range(order):
        paths[1 << vertex][vertex] = 1
    for mask in range(1, full + 1):
        for last in range(order):
            count = paths[mask][last]
            if not count:
                continue
            available = adjacency[last] & (full ^ mask)
            while available:
                next_bit = available & -available
                paths[mask | next_bit][next_bit.bit_length() - 1] += count
                available ^= next_bit
    return sum(paths[full])


def one_defect_count(adjacency: list[int]) -> int:
    total = 0
    for word in permutations(range(len(adjacency))):
        defects = sum(
            not ((adjacency[word[index]] >> word[index + 1]) & 1)
            for index in range(len(word) - 1)
        )
        total += defects == 1
    return total


def add_ear(adjacency: list[int], cut: int) -> list[int]:
    order = len(adjacency)
    child = adjacency[:] + [0]
    for vertex in range(order):
        if (cut >> vertex) & 1:
            child[order] |= 1 << vertex
        else:
            child[vertex] |= 1 << order
    return child


def is_strong(adjacency: list[int]) -> bool:
    order = len(adjacency)
    full = (1 << order) - 1
    for root in range(order):
        reached = 1 << root
        old = -1
        while reached != old:
            old = reached
            for vertex in range(order):
                if (reached >> vertex) & 1:
                    reached |= adjacency[vertex]
        if reached != full:
            return False
    return True


def odd_ceiling(value: Fraction) -> int:
    integer_ceiling = (value.numerator + value.denominator - 1) // value.denominator
    return integer_ceiling if integer_ceiling % 2 else integer_ceiling + 1


def fraction_row(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def main() -> None:
    order_rows: list[dict[str, object]] = []
    total_tournaments = 0
    total_ears = 0
    total_strong_ears = 0

    for order in range(2, 6):
        tournament_count = 1 << (order * (order - 1) // 2)
        strong_parents = 0
        identity_checks = 0
        ear_checks = 0
        strong_ear_checks = 0
        worst: tuple[Fraction, int, int, int, int, Fraction] | None = None

        for code in range(tournament_count):
            adjacency = tournament_from_code(order, code)
            parent_h = hamilton_count(adjacency)
            defect_one = one_defect_count(adjacency)
            ear_values = [
                hamilton_count(add_ear(adjacency, cut))
                for cut in range(1 << order)
            ]
            nonconstant = ear_values[1:-1]
            total_sum = sum(ear_values)
            predicted_sum = (1 << (order - 2)) * (
                (order + 3) * parent_h + defect_one
            )
            require(total_sum == predicted_sum, "all-cut exact average identity")
            require(ear_values[0] == ear_values[-1] == parent_h,
                    "constant source/sink ears")

            predicted_nonconstant_mean = Fraction(
                predicted_sum - 2 * parent_h, (1 << order) - 2
            )
            actual_nonconstant_mean = Fraction(sum(nonconstant), len(nonconstant))
            require(actual_nonconstant_mean == predicted_nonconstant_mean,
                    "nonconstant-ear exact mean")
            require(
                max(nonconstant) >= odd_ceiling(Fraction(order + 3, 4) * parent_h),
                "coarse odd-ceiling growth bound",
            )
            require(all(value % 2 for value in ear_values), "Redei parity control")

            slack = Fraction(max(nonconstant), parent_h) - Fraction(order + 3, 4)
            candidate = (
                slack,
                code,
                parent_h,
                defect_one,
                max(nonconstant),
                actual_nonconstant_mean,
            )
            if worst is None or candidate < worst:
                worst = candidate

            if is_strong(adjacency):
                strong_parents += 1
                for cut in range(1, (1 << order) - 1):
                    require(is_strong(add_ear(adjacency, cut)), "strong-ear inheritance")
                    strong_ear_checks += 1

            identity_checks += 1
            ear_checks += 1 << order

        require(worst is not None, "nonempty order census")
        slack, code, parent_h, defect_one, maximum, nonconstant_mean = worst
        order_rows.append(
            {
                "order": order,
                "tournaments": tournament_count,
                "strong_parents": strong_parents,
                "identity_checks": identity_checks,
                "ear_checks": ear_checks,
                "strong_ear_checks": strong_ear_checks,
                "minimum_max_ratio_slack": fraction_row(slack),
                "minimum_row": [
                    code,
                    parent_h,
                    defect_one,
                    maximum,
                    *fraction_row(nonconstant_mean),
                ],
            }
        )
        total_tournaments += tournament_count
        total_ears += ear_checks
        total_strong_ears += strong_ear_checks

    # The directed 3-cycle is the sharp boundary for a tempting stronger
    # coefficient.  In the LSB-first encoding, code 2 is 0->2->1->0.
    cycle = tournament_from_code(3, 2)
    cycle_h = hamilton_count(cycle)
    cycle_f1 = one_defect_count(cycle)
    cycle_ears = [hamilton_count(add_ear(cycle, cut)) for cut in range(8)]
    require(is_strong(cycle), "C3 is strong")
    require((cycle_h, cycle_f1) == (3, 0), "C3 zero/one-defect layers")
    require(cycle_ears == [3, 5, 5, 5, 5, 5, 5, 3], "C3 ear multiset")
    require(max(cycle_ears[1:-1]) == 5 < Fraction(3 + 4, 4) * cycle_h,
            "C3 refutes n+4 coefficient")

    # Frozen observed maxima from THM-4097/4102/4104.  Only the lower bounds
    # are consequences of THM-4111; the larger observed values remain finite
    # exact data from those theorems.
    selected_rows = []
    for order, parent_maximum, child_maximum in (
        (9, 3_357, 15_621),
        (10, 15_621, 93_751),
    ):
        forced = odd_ceiling(Fraction(order + 3, 4) * parent_maximum)
        require(child_maximum >= forced, "observed selected-bank growth control")
        selected_rows.append([order, parent_maximum, forced, child_maximum])

    ledger = {
        "orders": order_rows,
        "totals": [total_tournaments, total_ears, total_strong_ears],
        "c3_hostile": [cycle_h, cycle_f1, cycle_ears],
        "selected_bank_rows": selected_rows,
        "formula": "sum_cut_H=2^(n-2)*((n+3)*H+F1)",
        "scope": "image maxima, not solid-interval endpoints or overlap",
    }
    semantic = canonical_digest(ledger)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic ledger digest")

    print("THM-4111 UNIFORM CUT-EAR AVERAGE PRIMARY REFEREE")
    print("status=PROVISIONAL_PRIMARY_PASS")
    print("formula=sum_S H(T+x_S)=2^(n-2)*((n+3)H(T)+F1(T))")
    print("exhaustive_order_rows=", order_rows)
    print("totals_tournaments_ears_strong_ears=", ledger["totals"])
    print("C3_hostile_H_F1_ears=", ledger["c3_hostile"])
    print("C3_refutes_stronger_coefficient=(n+4)/4")
    print("selected_bank_order_parent_forced_observed=", selected_rows)
    print("semantic_sha256=", semantic)
    print("semantic_source=all labelled tournaments through order five")
    print("semantic_map=uniformly average the exact cut-ear boundary")
    print("semantic_preserved=Hamiltonian count, cut signature, and one-defect layer")
    print("semantic_destroyed=individual cut values after averaging")
    print("semantic_scope=unbounded recursive image maxima only; interval overlap remains OPEN")
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
