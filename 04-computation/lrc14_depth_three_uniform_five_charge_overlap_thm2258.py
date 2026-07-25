#!/usr/bin/env python3
"""Exact two-path verifier for THM-2258's five-charge overlap exclusion.

For a scalar (3,4,5) survivor, THM-2243 bounds the literal five-unit
residual

    p = measure(C_H minus union_i D_(q_i))

by B=8907541/62377224.  The union bound then forces the sum of the five
individual charges measure(C_H intersect D_(q_i)) to be at least 5/7-B.
THM-2080's exact odd-guard charge law and its 1/(4ab) tail reduce all five
distinct reduced ratios to a 13-row bank.  Only 15 of its 1287 five-subsets
meet the forced charge threshold.

Path A evaluates the THM-2080 formula exactly on the complete small-product
bank.  Path B independently reconstructs every charge and every literal
five-comb residual from rational interval atoms.  After common-scale
normalization, the 15 hostile marginal rows all have residual at least
317/1155, strictly above B.  Thus near-extremal marginal charge forces
catastrophic overlap and the profile (3,4,5) is empty.

All arithmetic is Fraction-exact; no floating-point value enters a proof
assertion.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm


Q = Fraction
P = 13
TARGET = Q(961, 6930)
PROFILE = (3, 4, 5)
UNIFORM_BELLMAN_BOUND = Q(8907541, 62377224)
FORCED_CHARGE_SUM = Q(5, 7) - UNIFORM_BELLMAN_BOUND


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle_distance(value):
    value %= 1
    return min(value, 1 - value)


def thm2080_eligibility_overlap(a, b):
    """Exact measure of E_a intersect D_b for coprime reduced (a,b)."""
    require(a > 0 and b > 0 and gcd(a, b) == 1, "ratio is not reduced")
    x = Q(a % 14, 14)
    y = Q(b % 7, 7)
    correction = min(x, y) + max(x + y - 1, Q()) - 2 * x * y
    return Q(2, 49) + Q(2, a * b) * correction


def guard_danger_capacity(a, b):
    """Exact measure of C_a intersect D_b."""
    return Q(1, 7) - thm2080_eligibility_overlap(a, b)


def comb_endpoints(speed, numerator):
    """Endpoints for {||speed*x||<numerator/14} on the circle."""
    require(speed > 0 and numerator in (1, 2), "endpoint request drift")
    endpoints = set()
    for tooth in range(speed):
        endpoints.add(Q(14 * tooth - numerator, 14 * speed) % 1)
        endpoints.add(Q(14 * tooth + numerator, 14 * speed) % 1)
    return endpoints


def atom_sweep(guard, dangers, *, capacity_mode=False):
    """Direct exact interval measure, independent of the charge formula."""
    endpoints = {Q(), Q(1)}
    endpoints.update(comb_endpoints(guard, 2))
    for danger in dangers:
        endpoints.update(comb_endpoints(danger, 1))
    endpoints = tuple(sorted(endpoints))
    measure = Q()
    active_atoms = 0
    for left, right in zip(endpoints, endpoints[1:]):
        midpoint = (left + right) / 2
        in_guard = circle_distance(guard * midpoint) > Q(1, 7)
        danger_bits = tuple(
            circle_distance(danger * midpoint) < Q(1, 14)
            for danger in dangers
        )
        active = (
            in_guard and any(danger_bits)
            if capacity_mode
            else in_guard and not any(danger_bits)
        )
        if active:
            measure += right - left
            active_atoms += 1
    return measure, active_atoms, len(endpoints) - 1


def direct_capacity_audit():
    """Check every admissible reduced row with ab<=50 by interval atoms."""
    checks = 0
    for a in range(1, 51, 2):
        for b in range(1, 50 // a + 1):
            if gcd(a, b) != 1 or a % P == 0 or b % P == 0:
                continue
            direct, _, _ = atom_sweep(a, (b,), capacity_mode=True)
            require(
                direct == guard_danger_capacity(a, b),
                f"direct capacity formula failed at {(a, b)}",
            )
            checks += 1
    return checks


TOP_CAPACITIES = (
    Q(5, 42),
    Q(9, 77),
    Q(4, 35),
    Q(4, 35),
)
FOURTH_CAPACITY = TOP_CAPACITIES[-1]
TOP_FOUR_PRODUCT_CUTOFF = Q(1, 4) / (FOURTH_CAPACITY - Q(5, 49))
FIFTH_CUTOFF = FORCED_CHARGE_SUM - sum(TOP_CAPACITIES, Q())
PRODUCT_CUTOFF = Q(1, 4) / (FIFTH_CUTOFF - Q(5, 49))


def top_four_ratios():
    """Certify the four global maxima using the THM-2080 tail."""
    # A capacity at least FOURTH_CAPACITY has ab<=20 because
    # TOP_FOUR_PRODUCT_CUTOFF=245/12 lies strictly between 20 and 21.
    rows = []
    for a in range(1, 21, 2):
        for b in range(1, 20 // a + 1):
            if gcd(a, b) != 1 or a % P == 0 or b % P == 0:
                continue
            rows.append((guard_danger_capacity(a, b), a, b))
    return tuple(
        sorted(rows, key=lambda row: (-row[0], row[1], row[2]))[:4]
    )


TOP_FOUR_RATIOS = top_four_ratios()
EXPECTED_TOP_FOUR_RATIOS = (
    (Q(5, 42), 1, 6),
    (Q(9, 77), 11, 1),
    (Q(4, 35), 1, 5),
    (Q(4, 35), 3, 5),
)
require(
    TOP_FOUR_RATIOS == EXPECTED_TOP_FOUR_RATIOS,
    "global top-four charge spectrum drift",
)


DIRECT_CAPACITY_CHECKS = direct_capacity_audit()


def candidate_ratios():
    """Complete thirteen-unit ratio bank above the forced fifth cutoff."""
    # THM-2080 gives capacity <=5/49+1/(4ab).  PRODUCT_CUTOFF is
    # 50.591..., so a row at least FIFTH_CUTOFF has ab<=50.
    rows = []
    for a in range(1, 51, 2):
        for b in range(1, 50 // a + 1):
            if gcd(a, b) != 1 or a % P == 0 or b % P == 0:
                continue
            capacity = guard_danger_capacity(a, b)
            if capacity >= FIFTH_CUTOFF:
                rows.append((capacity, a, b))
    return tuple(sorted(rows, key=lambda row: (-row[0], row[1], row[2])))


CANDIDATES = candidate_ratios()
EXPECTED_CANDIDATES = (
    (Q(5, 42), 1, 6),
    (Q(9, 77), 11, 1),
    (Q(4, 35), 1, 5),
    (Q(4, 35), 3, 5),
    (Q(1, 9), 9, 1),
    (Q(1, 9), 9, 2),
    (Q(17, 154), 11, 2),
    (Q(19, 175), 25, 1),
    (Q(3, 28), 1, 4),
    (Q(3, 28), 1, 12),
    (Q(3, 28), 1, 20),
    (Q(3, 28), 3, 4),
    (Q(3, 28), 5, 4),
)
require(CANDIDATES == EXPECTED_CANDIDATES, "candidate ratio bank drift")


def normalized_row(subset):
    """Normalize fixed-H reduced ratios by their common integer scale."""
    guard = 1
    for _, a, _ in subset:
        guard = lcm(guard, a)
    dangers = tuple(sorted(guard * b // a for _, a, b in subset))
    require(len(set(dangers)) == 5, "distinct reduced ratios collided")
    require(
        guard % P != 0 and all(danger % P != 0 for danger in dangers),
        "thirteen-unit normalization drift",
    )
    residual, active_atoms, atom_count = atom_sweep(guard, dangers)
    ratios = tuple(sorted((a, b) for _, a, b in subset))
    charge_sum = sum((capacity for capacity, _, _ in subset), Q())
    direct_charge_sum = sum(
        (atom_sweep(a, (b,), capacity_mode=True)[0] for _, a, b in subset),
        Q(),
    )
    require(charge_sum == direct_charge_sum, "two-path charge sum drift")
    overlap_rebate = residual + charge_sum - Q(5, 7)
    return (
        ratios,
        guard,
        dangers,
        charge_sum,
        residual,
        overlap_rebate,
        active_atoms,
        atom_count,
    )


QUALIFYING_SUBSETS = tuple(
    subset
    for subset in combinations(CANDIDATES, 5)
    if sum((row[0] for row in subset), Q()) >= FORCED_CHARGE_SUM
)
ROWS = tuple(sorted(normalized_row(subset) for subset in QUALIFYING_SUBSETS))

EXPECTED_ROWS = (
    (
        ((1, 4), (1, 5), (1, 6), (3, 5), (11, 1)),
        33,
        (3, 55, 132, 165, 198),
        FORCED_CHARGE_SUM + Q(183527, 1143582440),
        Q(131, 462),
        Q(31, 220),
        298,
        1167,
    ),
    (
        ((1, 5), (1, 6), (1, 12), (3, 5), (11, 1)),
        33,
        (3, 55, 165, 198, 396),
        FORCED_CHARGE_SUM + Q(183527, 1143582440),
        Q(733, 2310),
        Q(269, 1540),
        316,
        1695,
    ),
    (
        ((1, 5), (1, 6), (1, 20), (3, 5), (11, 1)),
        33,
        (3, 55, 165, 198, 660),
        FORCED_CHARGE_SUM + Q(183527, 1143582440),
        Q(661, 2310),
        Q(221, 1540),
        538,
        2223,
    ),
    (
        ((1, 5), (1, 6), (3, 4), (3, 5), (11, 1)),
        33,
        (3, 44, 55, 165, 198),
        FORCED_CHARGE_SUM + Q(183527, 1143582440),
        Q(317, 1155),
        Q(29, 220),
        262,
        991,
    ),
    (
        ((1, 5), (1, 6), (3, 5), (5, 4), (11, 1)),
        165,
        (15, 132, 275, 825, 990),
        FORCED_CHARGE_SUM + Q(183527, 1143582440),
        Q(1121, 3850),
        Q(3431, 23100),
        1162,
        4775,
    ),
    (
        ((1, 5), (1, 6), (3, 5), (9, 1), (11, 1)),
        99,
        (9, 11, 165, 495, 594),
        FORCED_CHARGE_SUM + Q(42493973, 10292241960),
        Q(6103, 20790),
        Q(46, 297),
        674,
        2707,
    ),
    (
        ((1, 5), (1, 6), (3, 5), (9, 2), (11, 1)),
        99,
        (9, 22, 165, 495, 594),
        FORCED_CHARGE_SUM + Q(42493973, 10292241960),
        Q(5977, 20790),
        Q(221, 1485),
        674,
        2707,
    ),
    (
        ((1, 5), (1, 6), (3, 5), (11, 1), (11, 2)),
        33,
        (3, 6, 55, 165, 198),
        FORCED_CHARGE_SUM + Q(3896457, 1143582440),
        Q(237, 770),
        Q(389, 2310),
        234,
        915,
    ),
    (
        ((1, 5), (1, 6), (3, 5), (11, 1), (25, 1)),
        825,
        (33, 75, 1375, 4125, 4950),
        FORCED_CHARGE_SUM + Q(9086081, 5717912200),
        Q(2407, 8250),
        Q(207, 1375),
        5602,
        22551,
    ),
    (
        ((1, 5), (1, 6), (9, 1), (9, 2), (11, 1)),
        99,
        (9, 11, 22, 495, 594),
        FORCED_CHARGE_SUM + Q(9820189, 10292241960),
        Q(3041, 10395),
        Q(3133, 20790),
        660,
        2399,
    ),
    (
        ((1, 5), (1, 6), (9, 1), (11, 1), (11, 2)),
        99,
        (9, 11, 18, 495, 594),
        FORCED_CHARGE_SUM + Q(342047, 1470320280),
        Q(3092, 10395),
        Q(46, 297),
        656,
        2435,
    ),
    (
        ((1, 5), (1, 6), (9, 2), (11, 1), (11, 2)),
        99,
        (9, 18, 22, 495, 594),
        FORCED_CHARGE_SUM + Q(342047, 1470320280),
        Q(2047, 6930),
        Q(353, 2310),
        638,
        2413,
    ),
    (
        ((1, 6), (3, 5), (9, 1), (9, 2), (11, 1)),
        99,
        (9, 11, 22, 165, 594),
        FORCED_CHARGE_SUM + Q(9820189, 10292241960),
        Q(3203, 10395),
        Q(3457, 20790),
        458,
        1735,
    ),
    (
        ((1, 6), (3, 5), (9, 1), (11, 1), (11, 2)),
        99,
        (9, 11, 18, 165, 594),
        FORCED_CHARGE_SUM + Q(342047, 1470320280),
        Q(6451, 20790),
        Q(317, 1890),
        442,
        1771,
    ),
    (
        ((1, 6), (3, 5), (9, 2), (11, 1), (11, 2)),
        99,
        (9, 18, 22, 165, 594),
        FORCED_CHARGE_SUM + Q(342047, 1470320280),
        Q(2113, 6930),
        Q(25, 154),
        442,
        1771,
    ),
)
require(
    ROWS == EXPECTED_ROWS,
    "qualifying normalized residual atlas drift",
)


def main():
    require(P == 13 and PROFILE == (3, 4, 5), "frozen universe drift")
    require(
        UNIFORM_BELLMAN_BOUND == Q(8907541, 62377224),
        "uniform Bellman bound drift",
    )
    require(
        UNIFORM_BELLMAN_BOUND - TARGET == Q(42493973, 10292241960),
        "uniform hostile gap drift",
    )
    require(
        FORCED_CHARGE_SUM == Q(5092517, 8911032),
        "forced five-charge sum drift",
    )
    require(
        TOP_FOUR_PRODUCT_CUTOFF == Q(245, 12)
        and 20 < TOP_FOUR_PRODUCT_CUTOFF < 21,
        "top-four small-product reduction drift",
    )
    require(
        FIFTH_CUTOFF == Q(122343163, 1143582440),
        "fifth-ratio cutoff drift",
    )
    require(
        PRODUCT_CUTOFF == Q(2001269270, 39557541)
        and 50 < PRODUCT_CUTOFF < 51,
        "small-product reduction drift",
    )
    require(DIRECT_CAPACITY_CHECKS == 102, "direct capacity audit drift")
    require(len(CANDIDATES) == 13, "candidate ratio count drift")
    require(len(tuple(combinations(CANDIDATES, 5))) == 1287,
            "candidate subset universe drift")
    require(len(QUALIFYING_SUBSETS) == 15, "qualifying subset count drift")
    require(all(row[4] > UNIFORM_BELLMAN_BOUND for row in ROWS),
            "uniform residual atlas failed to close")

    worst = min(ROWS, key=lambda row: row[4])
    require(
        worst[0] == ((1, 5), (1, 6), (3, 4), (3, 5), (11, 1)),
        "worst ratio row drift",
    )
    require(worst[4] == Q(317, 1155), "worst residual drift")
    require(
        worst[4] - UNIFORM_BELLMAN_BOUND
        == Q(150561431, 1143582440),
        "worst residual margin drift",
    )
    require(
        min(row[5] for row in ROWS) == Q(29, 220),
        "minimum overlap rebate drift",
    )
    transcript = "\n".join(
        (
            f"{row[0]}|H={row[1]}|q={row[2]}|"
            f"sum={row[3]}|p={row[4]}|rebate={row[5]}|"
            f"active={row[6]}|atoms={row[7]}"
        )
        for row in ROWS
    )
    digest = sha256(transcript.encode("ascii")).hexdigest()

    print("THM-2258 LRC14 DEPTH-THREE UNIFORM FIVE-CHARGE OVERLAP")
    print(f"profile={PROFILE} target={TARGET}")
    print(
        f"uniform_bellman_bound={UNIFORM_BELLMAN_BOUND} "
        f"bound_minus_target={UNIFORM_BELLMAN_BOUND-TARGET}"
    )
    print(f"forced_five_charge_sum={FORCED_CHARGE_SUM}")
    print(
        f"global_top_four={TOP_FOUR_RATIOS} "
        f"top_four_product_cutoff={TOP_FOUR_PRODUCT_CUTOFF}"
    )
    print(
        f"fifth_ratio_cutoff={FIFTH_CUTOFF} "
        f"product_cutoff={PRODUCT_CUTOFF} hence_ab<=50"
    )
    print(f"direct_interval_capacity_checks={DIRECT_CAPACITY_CHECKS}")
    print(f"candidate_ratios={CANDIDATES}")
    print(
        "candidate_five_subsets=1287 "
        f"qualifying_capacity_subsets={len(QUALIFYING_SUBSETS)}"
    )
    for (
        ratios,
        guard,
        dangers,
        charge_sum,
        residual,
        rebate,
        active,
        atoms,
    ) in ROWS:
        print(
            f"ratios={ratios} H={guard} q={dangers} "
            f"charge_sum={charge_sum} "
            f"residual={residual} residual_minus_bound="
            f"{residual-UNIFORM_BELLMAN_BOUND} "
            f"overlap_rebate={rebate} "
            f"active_atoms={active}/{atoms}"
        )
    print(
        f"worst_residual={worst[4]} "
        f"worst_minus_uniform_bound={worst[4]-UNIFORM_BELLMAN_BOUND}"
    )
    print(f"minimum_overlap_rebate={min(row[5] for row in ROWS)}")
    print(f"row_digest={digest}")
    print("consequence=scalar_depth_three_profile_(3,4,5)_EMPTY")
    print("independent_ledger_effect=166_to_165")
    print("canonical_ledger=165_already_by_THM-2257_no_additional_decrement")
    print("scope=LRC14_OPEN_165_first_depth_one_rows_and_non_scalar_frontiers")


if __name__ == "__main__":
    main()
