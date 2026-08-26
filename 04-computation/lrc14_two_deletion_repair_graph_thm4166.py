#!/usr/bin/env python3
"""Exact Fraction audit for THM-4166.

The script reconstructs all two-deletion Haar banks for the fixed THM-4156
pool, builds every repair graph Gamma_q for 1 <= q <= 200 with q outside the
pool, and computes its exact vertex-cover number through a maximum independent
set search.  Every threshold comparison is made with Fraction; equality is an
edge by definition.
"""

from fractions import Fraction as F
from itertools import combinations
from math import comb, floor
import sys


sys.stdout.reconfigure(newline="\n")


DELTA = F(1, 14)
THRESHOLD = F(4, 63)
ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(value for value in POOL if value not in ANCHORS)
SEARCH_LIMIT = 200


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def intersect(intervals, speed):
    """Intersect a disjoint closed interval bank with the safe comb G_speed."""
    answer = []
    for left, right in intervals:
        low = max(0, floor(speed * left) - 1)
        high = min(speed - 1, floor(speed * right) + 1)
        for tooth in range(low, high + 1):
            a = max(left, F(14 * tooth + 1, 14 * speed))
            b = min(right, F(14 * tooth + 13, 14 * speed))
            if a <= b:
                if answer and answer[-1][1] == a:
                    answer[-1] = (answer[-1][0], b)
                else:
                    answer.append((a, b))
    return tuple(answer)


def components(speeds):
    bank = ((F(0), F(1)),)
    for speed in speeds:
        bank = intersect(bank, speed)
    return bank


def mass(bank):
    return sum((right - left for left, right in bank), F(0))


def maximum_independent_set(adjacency):
    """Exact maximum clique search in the complement graph."""
    order = len(adjacency)
    universe = (1 << order) - 1
    complement = tuple(
        universe & ~(1 << vertex) & ~adjacency[vertex]
        for vertex in range(order)
    )
    best_size = 0
    best_mask = 0

    def search(chosen, candidates, excluded):
        nonlocal best_size, best_mask
        if chosen.bit_count() + candidates.bit_count() <= best_size:
            return
        if not candidates and not excluded:
            size = chosen.bit_count()
            if size > best_size:
                best_size = size
                best_mask = chosen
            return
        pivot_pool = candidates | excluded
        if pivot_pool:
            pivot = max(
                (v for v in range(order) if pivot_pool >> v & 1),
                key=lambda v: (candidates & complement[v]).bit_count(),
            )
            extensions = candidates & ~complement[pivot]
        else:
            extensions = candidates
        while extensions:
            bit = extensions & -extensions
            vertex = bit.bit_length() - 1
            search(
                chosen | bit,
                candidates & complement[vertex],
                excluded & complement[vertex],
            )
            candidates &= ~bit
            excluded |= bit
            extensions &= ~bit
            if chosen.bit_count() + candidates.bit_count() <= best_size:
                return

    search(0, universe, 0)
    return best_size, best_mask


def graph_for(newcomer, leave_two):
    adjacency = [0] * len(OPTIONAL)
    equality_edges = []
    edges = []
    positive_margins = []
    negative_deficits = []
    for (first, second), bank in leave_two.items():
        value = mass(intersect(bank, newcomer))
        difference = value - THRESHOLD
        if difference >= 0:
            i = OPTIONAL.index(first)
            j = OPTIONAL.index(second)
            adjacency[i] |= 1 << j
            adjacency[j] |= 1 << i
            edges.append((first, second))
            if difference == 0:
                equality_edges.append((first, second))
            else:
                positive_margins.append((difference, first, second))
        else:
            negative_deficits.append((-difference, first, second))
    return (
        tuple(adjacency),
        tuple(edges),
        tuple(equality_edges),
        tuple(positive_margins),
        tuple(negative_deficits),
    )


def fraction_text(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    require(len(POOL) == 30, "pool size changed")
    require(len(OPTIONAL) == 27, "optional-pool size changed")
    require(set(ANCHORS) <= set(POOL), "anchors left the pool")
    require(set(ANCHORS).isdisjoint(OPTIONAL), "anchor/optional collision")

    leave_two = {}
    for first, second in combinations(OPTIONAL, 2):
        speeds = tuple(v for v in POOL if v not in (first, second))
        leave_two[(first, second)] = components(speeds)
    require(len(leave_two) == comb(27, 2) == 351, "leave-two count changed")

    labels = tuple(q for q in range(1, SEARCH_LIMIT + 1) if q not in POOL)
    require(len(labels) == 175, "q universe changed")

    rows = []
    global_positive = []
    global_negative = []
    equality_total = 0
    for newcomer in labels:
        adjacency, edges, equalities, positive, negative = graph_for(
            newcomer, leave_two
        )
        alpha, witness_mask = maximum_independent_set(adjacency)
        tau = len(OPTIONAL) - alpha
        require(alpha == witness_mask.bit_count(), "MIS witness size mismatch")
        for vertex in range(len(OPTIONAL)):
            if witness_mask >> vertex & 1:
                require(
                    adjacency[vertex] & witness_mask == 0,
                    "reported MIS witness is not independent",
                )
        rows.append(
            {
                "q": newcomer,
                "adjacency": adjacency,
                "edge_count": len(edges),
                "equality_count": len(equalities),
                "alpha": alpha,
                "tau": tau,
                "witness": tuple(
                    OPTIONAL[index]
                    for index in range(len(OPTIONAL))
                    if witness_mask >> index & 1
                ),
            }
        )
        equality_total += len(equalities)
        global_positive.extend((d, newcomer, r, s) for d, r, s in positive)
        global_negative.extend((d, newcomer, r, s) for d, r, s in negative)

    admitted = tuple(row for row in rows if row["tau"] > 7)
    rejected = tuple(row for row in rows if row["tau"] <= 7)
    admitted_labels = tuple(row["q"] for row in admitted)
    expected_admitted = (
        4, 5, 9, 13, 18, 21, 27, 32, 33, 34, 39, 41, 43, 49, 54, 56,
        66, 68, 77, 78, 81, 82, 86, 91, 97, 98, 101, 107, 115, 117,
        118, 119, 121, 131, 133, 138, 139, 142, 149, 152, 156, 164,
        167, 171, 174, 178, 181, 182, 187, 191, 194, 196, 199,
    )
    require(admitted_labels == expected_admitted, "admitted q census changed")
    require(len(admitted) == 53 and len(rejected) == 122, "census split changed")
    require(equality_total == 0, "unexpected threshold equality")

    tau_histogram = {}
    for row in rows:
        tau_histogram[row["tau"]] = tau_histogram.get(row["tau"], 0) + 1
    expected_histogram = {
        0: 25, 1: 25, 2: 25, 3: 7, 4: 10, 5: 11, 6: 11, 7: 8,
        8: 6, 9: 12, 10: 15, 11: 3, 12: 8, 13: 3, 15: 2, 16: 1,
        18: 2, 19: 1,
    }
    require(tau_histogram == expected_histogram, "tau histogram changed")

    new_body_count = len(admitted) * comb(27, 7)
    old_body_count = comb(27, 8)
    require(new_body_count == 47_065_590, "new body count changed")
    require(old_body_count == 2_220_075, "old body count changed")

    print("THM4166_PRIMARY_EXACT_FRACTION_AUDIT")
    print("pool_size", len(POOL), "optional_size", len(OPTIONAL))
    print("leave_two_banks", len(leave_two))
    print("q_universe", f"1..{SEARCH_LIMIT} outside P", len(labels))
    print("threshold", fraction_text(THRESHOLD), "comparison", ">=")
    print("threshold_equalities", equality_total)
    print("tau_histogram", tuple(sorted(tau_histogram.items())))
    print("admitted_count", len(admitted))
    print("admitted_q", admitted_labels)
    print("admitted_q_tau", tuple((row["q"], row["tau"]) for row in admitted))
    smallest_positive = min(global_positive)
    smallest_negative = min(global_negative)
    print(
        "smallest_positive_margin",
        fraction_text(smallest_positive[0]),
        "q_r_s",
        smallest_positive[1:],
    )
    print(
        "smallest_negative_deficit",
        fraction_text(smallest_negative[0]),
        "q_r_s",
        smallest_negative[1:],
    )
    print("positive_control_q27", next(
        (row["edge_count"], row["tau"], row["alpha"])
        for row in rows if row["q"] == 27
    ))
    print("hostile_control_q200", next(
        (row["edge_count"], row["tau"], row["alpha"], row["witness"])
        for row in rows if row["q"] == 200
    ))
    print("bodies_per_q", comb(27, 7))
    print("new_two_deletion_bodies", new_body_count)
    print("with_old_thm4156_bodies", new_body_count + old_body_count)
    print("SEMANTIC_LEDGER_BEGIN")
    for row in rows:
        masks = ",".join(f"{mask:07x}" for mask in row["adjacency"])
        print(
            f"q={row['q']:03d};edges={row['edge_count']:03d};"
            f"eq={row['equality_count']};alpha={row['alpha']:02d};"
            f"tau={row['tau']:02d};adj={masks}"
        )
    print("SEMANTIC_LEDGER_END")
    print("PRIMARY_AUDIT_PASS")


if __name__ == "__main__":
    main()
