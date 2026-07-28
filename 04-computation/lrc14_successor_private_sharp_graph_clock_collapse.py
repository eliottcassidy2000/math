#!/usr/bin/env python3
"""Exact successor-private incidence refinement of the sharp graph.

This companion combines three already constructed pieces of the canonical
LRC(14) common-x carrier without assigning them a theorem number:

* THM-2616's globally typed raw carrier and all delayed digits;
* THM-2629's sharp coefficient graph ``q=h, r=-h-1``; and
* THM-2630's later-tooth half bit and predecessor/carry identity.

For every positive middle rail, every nonzero equivariant clock step, every
retained future digit, and both values of the tooth half and binary carry, we
retain the exact predecessor digit before integration.  These 42,768 atoms
partition the corresponding sharp-graph events exactly.  Their positive
support defines a common-x incidence matrix in the two numerical digit labels
``(j,h)`` at the imposed clock pairing ``(ell4,ell5=ell4+d)``.

The test is deliberately performed before any Bockstein unit test, endpoint
restriction, marginalization over clocks, or choice of a hidden Perron sheet.
The matrix is **not** a proved transition: no theorem identifies its printed
``h`` label with the incoming ``j`` label of another incidence cell, or glues
the common-x witnesses across clocks.  The zero products below concern only
fixed-source, constant-offset candidate families.  Exact varying-offset and
clock-marginal hostiles are included to prevent a broader no-go reading.
"""

from collections import Counter
from itertools import permutations

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = 13
Q7 = 7
INV2 = 7
PREDECESSOR_SPEED = P**5
SUCCESSOR_SPEED = P**6
FUTURE_DIGITS = tuple(range(1, 12))
CLOCK_STEPS = tuple(range(1, Q7))

EXPECTED_CANDIDATES = 162 * len(CLOCK_STEPS) * len(FUTURE_DIGITS) * 2 * 2
EXPECTED_POSITIVE_ATOMS = 7_436
EXPECTED_PARTITION_CHECKS = 162 * len(CLOCK_STEPS) * len(FUTURE_DIGITS)
EXPECTED_ZERO_CELL_COUNTS = {1: 18, 2: 39, 3: 16, 4: 16, 5: 39, 6: 18}
EXPECTED_UNIVERSALLY_EMPTY_ELL4 = {
    1: (5,),
    2: (0, 4, 6),
    3: (2,),
    4: (5,),
    5: (0, 1, 3),
    6: (2,),
}
EXPECTED_D1_ZERO_CELLS = (
    (1, 5), (2, 5), (3, 5), (4, 5), (5, 5),
    (6, 0), (6, 1), (6, 2), (6, 5), (6, 6),
    (7, 3), (7, 4), (7, 5),
    (8, 5), (9, 5), (10, 5), (11, 5), (12, 5),
)
GENERIC_S = tuple(s for s in range(1, P) if s not in (6, 7))
VARYING_OFFSET_ITINERARY = (0, 1, 2, 3, 6, 4, 5)
BEST_NORMALIZED_ITINERARY = (0, 6, 2, 4, 3, 5, 1)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def zero_matrix():
    return [[0] * P for _ in range(P)]


def identity_matrix():
    return [[int(i == j) for j in range(P)] for i in range(P)]


def boolean_product(left, right):
    """Formal Boolean product after imposing equality of printed labels."""
    return [
        [
            int(any(left[i][k] and right[k][j] for k in range(P)))
            for j in range(P)
        ]
        for i in range(P)
    ]


def support_size(matrix):
    return sum(sum(row) for row in matrix)


def boolean_union(matrices):
    """Forget a matrix-family label by Boolean union (a hostile quotient)."""
    result = zero_matrix()
    for matrix in matrices:
        for i in range(P):
            for j in range(P):
                result[i][j] |= matrix[i][j]
    return result


def boolean_power(matrix, exponent):
    result = identity_matrix()
    for _ in range(exponent):
        result = boolean_product(result, matrix)
    return result


def itinerary_product(matrices, s, itinerary):
    """Formal product for a displayed varying-offset clock itinerary."""
    require(len(itinerary) == Q7 and set(itinerary) == set(range(Q7)),
            "itinerary is not a permutation of the seven clock labels")
    result = identity_matrix()
    for index, ell4 in enumerate(itinerary):
        ell5 = itinerary[(index + 1) % Q7]
        d = (ell5 - ell4) % Q7
        require(d in CLOCK_STEPS, "itinerary repeats a clock label")
        result = boolean_product(result, matrices[d, s, ell4])
    return result


def ambient_relation_controls():
    """Check that the local affine/digit grammar is not itself obstructed.

    This is the carrier-free relaxation retaining only

        h = -2*j-kappa-epsilon-1 mod 13

    and the condition that the h-digit interval meets the chosen kappa half.
    It intentionally forgets every rail, word, target gate, and owner clock.
    """
    relation = zero_matrix()
    for j in range(P):
        for kappa in (0, 1):
            for epsilon in (0, 1):
                h = (-2 * j - kappa - epsilon - 1) % P
                half_is_possible = (
                    (kappa == 0 and h <= 6)
                    or (kappa == 1 and h >= 6)
                )
                if half_is_possible:
                    relation[j][h] = 1

    row_histograms = []
    support_sizes = []
    power = identity_matrix()
    for exponent in range(1, 5):
        power = boolean_product(power, relation)
        support_sizes.append(support_size(power))
        row_histograms.append(Counter(sum(row) for row in power))

    require(
        row_histograms[0] == Counter({2: 10, 3: 2, 1: 1}),
        "ambient one-edge row histogram changed",
    )
    require(
        support_sizes == [27, 56, 113, 169],
        "ambient Boolean-power support sequence changed",
    )
    require(
        row_histograms[3] == Counter({13: 13}),
        "ambient fourth power is no longer the full relation",
    )
    twisted_returns = sum(
        power[j][(j - 7 * a) % P]
        for a in range(1, P)
        for j in range(P)
    )
    require(twisted_returns == 12 * 13,
            "ambient full relation lost a twisted return")
    return tuple(support_sizes), tuple(row_histograms), twisted_returns


def rebuild_private_incidence_bank():
    """Build every exact sharp-graph predecessor/successor incidence atom."""
    (module, prefixes, _, _, rails, present, starts) = cross.build_carrier_data()
    require(len(rails) == 162, "THM-2616 middle-rail census changed")
    require(module.C3 == 2 * PREDECESSOR_SPEED,
            "later probe no longer has speed 2*13^5")
    require(old.R == SUCCESSOR_SPEED,
            "delayed word no longer has clock 13^6")

    matrices = {
        (d, s, ell): zero_matrix()
        for d in CLOCK_STEPS
        for s in range(1, P)
        for ell in range(Q7)
    }
    candidates = 0
    positive_atoms = 0
    partition_checks = 0

    for s, ell, rail_digit, pieces in rails:
        require(rail_digit in (0, 12), "middle rail left its two-digit chart")
        for d in CLOCK_STEPS:
            future_clock = (ell + d) % Q7
            for h in FUTURE_DIGITS:
                # THM-2629's sharp graph.  Its forbidden r=0 value would occur
                # only at h=12, already absent from this delayed-word carrier.
                r = (-h - 1) % P
                require(r != 0, "sharp graph entered the absent deep sheet")
                overlap = old.intersect_weighted_union(
                    pieces,
                    present[future_clock, (-h) % P],
                    starts[future_clock, (-h) % P],
                )
                whole = old.intersect_weighted_comb(
                    overlap, module.C3, 182, 14 * r - 13, 14 * r + 13
                )
                whole_value = old.delayed_weighted_numerator(
                    whole, prefixes[future_clock][h]
                )

                split_total = 0
                for epsilon, left, right in (
                    (0, 14 * r, 14 * r + 13),
                    (1, 14 * r - 13, 14 * r),
                ):
                    half_tooth = old.intersect_weighted_comb(
                        overlap, module.C3, 182, left, right
                    )
                    for kappa in (0, 1):
                        candidates += 1
                        # r=2*j+kappa+epsilon, with 2^(-1)=7 mod 13.
                        j = INV2 * (r - epsilon - kappa) % P
                        require(
                            h == (-2 * j - kappa - epsilon - 1) % P,
                            "sharp graph/predecessor incidence equation failed",
                        )
                        predecessor_piece = old.intersect_weighted_comb(
                            half_tooth,
                            PREDECESSOR_SPEED,
                            P,
                            j,
                            j + 1,
                        )
                        carry_piece = old.intersect_weighted_comb(
                            predecessor_piece,
                            SUCCESSOR_SPEED,
                            2,
                            kappa,
                            kappa + 1,
                        )
                        value = old.delayed_weighted_numerator(
                            carry_piece, prefixes[future_clock][h]
                        )
                        require(value >= 0, "negative exact atom")
                        split_total += value
                        if value:
                            positive_atoms += 1
                            matrices[d, s, ell][j][h] = 1

                require(
                    split_total == whole_value,
                    "successor-private atoms do not partition the graph event",
                )
                partition_checks += 1

    require(candidates == EXPECTED_CANDIDATES,
            "successor-private candidate universe changed")
    require(positive_atoms == EXPECTED_POSITIVE_ATOMS,
            "successor-private positive census changed")
    require(partition_checks == EXPECTED_PARTITION_CHECKS,
            "sharp-graph partition-check universe changed")
    return matrices, candidates, positive_atoms, partition_checks


def constant_offset_candidate_controls(matrices):
    """Test fixed-source, constant-offset numerical candidate products."""
    zero_cells = {}
    universally_empty = {}
    zero_products = 0
    product_support_histograms = {}

    for d in CLOCK_STEPS:
        zero_cells[d] = tuple(
            (s, ell)
            for s in range(1, P)
            for ell in range(Q7)
            if support_size(matrices[d, s, ell]) == 0
        )
        require(
            len(zero_cells[d]) == EXPECTED_ZERO_CELL_COUNTS[d],
            f"clock-step {d} zero-cell census changed",
        )
        universally_empty[d] = tuple(
            ell4 for ell4 in range(Q7)
            if all(support_size(matrices[d, s, ell4]) == 0
                   for s in range(1, P))
        )
        require(
            universally_empty[d] == EXPECTED_UNIVERSALLY_EMPTY_ELL4[d],
            f"clock-step {d} universal empty-ell4 set changed",
        )

        product_hist = Counter()
        clock_order = tuple((i * d) % Q7 for i in range(Q7))
        require(len(set(clock_order)) == Q7,
                "nonzero clock step failed to traverse F_7")
        for s in range(1, P):
            product = identity_matrix()
            for ell in clock_order:
                product = boolean_product(product, matrices[d, s, ell])
            product_support = support_size(product)
            product_hist[product_support] += 1
            zero_products += int(product_support == 0)
        require(product_hist == Counter({0: 12}),
                f"constant-offset candidate {d} acquired formal support")
        product_support_histograms[d] = product_hist

    require(zero_cells[1] == EXPECTED_D1_ZERO_CELLS,
            "canonical adjacent-clock zero-cell atlas changed")
    require(zero_products == 12 * len(CLOCK_STEPS) == 72,
            "not all fixed-source constant-offset candidates vanish")
    return (
        zero_cells, universally_empty,
        product_support_histograms, zero_products,
    )


def broader_no_go_hostiles(matrices):
    """Exhibit why the constant-offset zero is not a chronology no-go."""
    varying_support = {
        s: support_size(itinerary_product(
            matrices, s, VARYING_OFFSET_ITINERARY
        ))
        for s in range(1, P)
    }
    require(
        all(varying_support[s] == 70 for s in GENERIC_S)
        and varying_support[6] == varying_support[7] == 0,
        "displayed varying-offset hostile changed",
    )

    # Normalize a seven-clock itinerary to begin at labelled clock zero.
    # Exhaust all 6!=720 possibilities for one generic source and both
    # exceptional sources.  The displayed maximizing itinerary is then
    # checked on every other generic source.
    best = {}
    for s in (1, 6, 7):
        candidates = []
        for tail in permutations(range(1, Q7)):
            itinerary = (0,) + tail
            candidates.append((
                support_size(itinerary_product(matrices, s, itinerary)),
                itinerary,
            ))
        maximum = max(value for value, _ in candidates)
        first_maximizer = min(
            itinerary for value, itinerary in candidates if value == maximum
        )
        best[s] = (maximum, first_maximizer)
    require(best[1] == (88, BEST_NORMALIZED_ITINERARY),
            "generic normalized-itinerary maximum changed")
    require(best[6][0] == best[7][0] == 0,
            "an exceptional source acquired a varying-offset candidate")
    require(
        all(support_size(itinerary_product(
            matrices, s, BEST_NORMALIZED_ITINERARY
        )) == 88 for s in GENERIC_S),
        "displayed best itinerary is not uniformly generic",
    )

    # This intentionally forgets ell4 before taking a seventh power.  Its
    # positive result is a hostile to any claim that the zero survives clock
    # marginalization; it is not a physically composable transition.
    marginal_support = {}
    for d in CLOCK_STEPS:
        for s in range(1, P):
            marginal = boolean_union(
                matrices[d, s, ell4] for ell4 in range(Q7)
            )
            marginal_support[d, s] = support_size(
                boolean_power(marginal, Q7)
            )
    require(
        Counter(marginal_support.values()) == Counter({110: 71, 75: 1}),
        "clock-marginal seventh-power hostile changed",
    )
    require(marginal_support[5, 6] == 75,
            "clock-marginal exceptional cell changed")
    require(
        all(value == 110 for key, value in marginal_support.items()
            if key != (5, 6)),
        "a generic clock-marginal seventh power changed",
    )
    return varying_support, best, marginal_support


def histogram_string(counter):
    return tuple(sorted(counter.items()))


def main():
    ambient_sizes, ambient_rows, ambient_twists = ambient_relation_controls()
    matrices, candidates, positives, partitions = rebuild_private_incidence_bank()
    (
        zero_cells, universally_empty, product_hists, zero_products,
    ) = constant_offset_candidate_controls(matrices)
    varying_support, best, marginal_support = broader_no_go_hostiles(matrices)

    print("LRC14 exact successor-private sharp-graph incidence audit")
    print("scope=canonical THM-2616 raw nonnegative carrier; no unit or endpoint restriction")
    print("typing=common-x incidence at imposed (ell4,ell5=ell4+d); not a proved transition")
    print("graph=q=h, r=-h-1; incidence law h=-2*j-kappa-epsilon-1 mod13")
    print(f"rails=162 clock_steps=6 future_digits=11 candidates={candidates}")
    print(
        f"positive_labelled_atom_occurrences={positives} "
        f"partition_checks={partitions}"
    )
    print(
        "zero_incidence_cells_by_offset="
        + str(tuple((d, len(zero_cells[d])) for d in CLOCK_STEPS))
    )
    print(
        "universally_empty_ell4_by_offset="
        + str(tuple((d, universally_empty[d]) for d in CLOCK_STEPS))
    )
    print(f"d1_zero_incidence_cells={zero_cells[1]}")
    print(
        "fixed_source_constant_offset_candidate_support_hist="
        + str(tuple((d, histogram_string(product_hists[d])) for d in CLOCK_STEPS))
    )
    print(f"fixed_source_constant_offset_zero_candidates={zero_products}/72")
    print(
        f"varying_offset_itinerary={VARYING_OFFSET_ITINERARY} "
        f"support_by_s={tuple(sorted(varying_support.items()))}"
    )
    print(
        f"displayed_itinerary_uniform_generic_support=88 "
        f"itinerary={BEST_NORMALIZED_ITINERARY} exceptional_s6_s7_best=0"
    )
    print(
        "clock_marginal_seventh_power_support_hist="
        + str(histogram_string(Counter(marginal_support.values())))
        + " exceptional=(d,s,support)=(5,6,75)"
    )
    print(f"ambient_relation_power_supports={ambient_sizes}")
    print(
        "ambient_relation_row_histograms="
        + str(tuple(histogram_string(hist) for hist in ambient_rows))
    )
    print(f"ambient_power4_twisted_returns={ambient_twists}/156")
    print(
        "verdict=PASS: fixed-source constant-offset candidate families meet "
        "universal empty incidence cells; varying offsets and clock marginalization survive"
    )
    print(
        "semantics=no adjacent-clock gluing or chronology follows; only "
        "label-preserving restrictions of a fixed (s,d,ell4) zero stay zero"
    )


if __name__ == "__main__":
    main()
