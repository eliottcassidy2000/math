#!/usr/bin/env python3
"""Exact delayed-guard cospan test for successor-private incidences.

The delayed word in THM-2616 contains the guard-safe sector.  This script
rebuilds the lawful disjoint Boolean cospan

    guard-safe  disjoint-union  guard-danger = guard-free word core,

then refines every sector *before integration* by the sharp graph
``q=h, r=-h-1`` and the exact successor-private labels
``(j,h,epsilon,kappa)``.  The guard-danger leg restores the future digit
``h=0``, whose sharp deep label is the admissible ``r=12``.  The other new
future digit ``h=12`` would require the structurally absent ``r=0`` sheet;
that zero is checked separately rather than silently discarded.

The primary question is whether the danger leg fills any of the eighteen
zero common-x incidence cells for the imposed pairing ``ell5=ell4+1``.  The
computation also exhausts all six nonzero equivariant steps, all twelve
source shifts, and all three labelled sectors (safe, danger, and their
guard-free union).  All arithmetic is exact integer interval arithmetic
inherited from THM-2616.
No Bockstein, endpoint restriction, or hidden-sheet unit claim is used.

The arrays are not proved transitions.  Their ``j`` and ``h`` labels are
numerically paired only inside one common-x incidence cell, and no physical
gluing identifies the output of one cell with the input of another.  The
fixed-source constant-offset products are therefore candidate products.  A
separate, intentionally lossy clock-marginal hostile is computed to show that
the zero statement is not preserved when ``ell4`` is forgotten.
"""

from collections import Counter

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old
import lrc14_successor_private_sharp_graph_clock_collapse as sharp


P = 13
Q7 = 7
INV2 = 7
PREDECESSOR_SPEED = P**5
SUCCESSOR_SPEED = P**6
ADMISSIBLE_GRAPH_DIGITS = tuple(range(12))
SECTORS = ("safe", "danger", "guard_free")

EXPECTED_CANDIDATES = 162 * 6 * 12 * 2 * 2
EXPECTED_GRAPH_PARTITIONS = 162 * 6 * 12
EXPECTED_R0_CHECKS = 162 * 6
EXPECTED_POSITIVE_ATOMS = {
    "safe": 7_436,
    "danger": 1_636,
    "guard_free": 8_360,
}
EXPECTED_POSITIVE_BY_H = {
    "safe": {
        1: 466, 2: 668, 3: 712, 4: 690, 5: 734, 6: 1_182,
        7: 602, 8: 734, 9: 690, 10: 712, 11: 246,
    },
    "danger": {0: 660, 1: 642, 11: 334},
    "guard_free": {
        0: 660, 1: 642, 2: 668, 3: 712, 4: 690, 5: 734,
        6: 1_182, 7: 602, 8: 734, 9: 690, 10: 712, 11: 334,
    },
}
EXPECTED_WORD_GRID_MASSES = {
    "guard_free": 15_722_865_944_474,
    "safe": 10_786_297_480_556,
    "danger": 4_936_568_463_918,
}
EXPECTED_DELAYED_TOTALS = {
    "safe": 64_717_784_883_336,
    "danger": 29_619_410_783_508,
    "guard_free": 94_337_195_666_844,
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def interval_mass(intervals):
    return sum(right - left for left, right in intervals)


def build_guard_cospan(module, inherited_safe_prefixes):
    """Build guard-safe, guard-danger, and guard-free delayed prefixes."""
    guard_free = module.make_comb(module.C2, 182, -13, 13)
    for index in module.UNIT_IDX:
        guard_free = module.subtract_comb(
            guard_free, module.W[index], 182, -13, 13
        )
    guard_free = module.subtract_comb(
        guard_free, module.C3, 182, -13, 13
    )
    safe = module.subtract_comb(
        guard_free, module.W[module.GUARD], 91, -13, 13
    )
    danger = module.intersect_comb(
        guard_free, module.W[module.GUARD], 91, -13, 13
    )

    require(safe == module.build_word_Ta(),
            "rebuilt guard-safe word disagrees with THM-2616")
    require(not old.intersect_sorted(safe, danger),
            "guard-safe and guard-danger sectors overlap")

    word_masses = {
        "safe": interval_mass(safe),
        "danger": interval_mass(danger),
        "guard_free": interval_mass(guard_free),
    }
    require(word_masses == EXPECTED_WORD_GRID_MASSES,
            "guard cospan word masses changed")
    require(word_masses["safe"] + word_masses["danger"]
            == word_masses["guard_free"],
            "guard sectors do not partition the word core")

    words = {"safe": safe, "danger": danger, "guard_free": guard_free}
    prefixes = {
        sector: [[None] * P for _ in range(Q7)] for sector in SECTORS
    }
    digit_supports = {sector: [] for sector in SECTORS}
    delayed_totals = Counter()
    for ell in range(Q7):
        clock_words = {
            sector: module.subtract_comb(
                words[sector], module.C1, 182,
                26 * ell - 13, 26 * ell + 13,
            )
            for sector in SECTORS
        }
        for h in range(P):
            digit_intervals = {
                sector: old.sat.intersect_interval(
                    clock_words[sector],
                    h * old.T // P,
                    (h + 1) * old.T // P,
                )
                for sector in SECTORS
            }
            require(not old.intersect_sorted(
                digit_intervals["safe"], digit_intervals["danger"]
            ), "delayed guard sectors overlap inside a digit")
            digit_masses = {
                sector: interval_mass(digit_intervals[sector])
                for sector in SECTORS
            }
            require(
                digit_masses["safe"] + digit_masses["danger"]
                == digit_masses["guard_free"],
                "delayed digit does not split into safe plus danger",
            )
            for sector in SECTORS:
                prefixes[sector][ell][h] = module.make_prefix(
                    digit_intervals[sector]
                )
                delayed_totals[sector] += digit_masses[sector]

        for sector in SECTORS:
            digit_supports[sector].append(tuple(
                h for h in range(P)
                if prefixes[sector][ell][h][2][-1] > 0
            ))

    require(prefixes["safe"] == inherited_safe_prefixes,
            "safe prefixes disagree with THM-2616's inherited bank")
    require(dict(delayed_totals) == EXPECTED_DELAYED_TOTALS,
            "delayed guard-sector total masses changed")
    require(
        digit_supports["safe"][0] == tuple(range(2, 11))
        and all(support == tuple(range(1, 12))
                for support in digit_supports["safe"][1:]),
        "safe delayed-digit support changed",
    )
    require(
        all(support == (0, 1, 11, 12)
            for support in digit_supports["danger"]),
        "danger delayed-digit support changed",
    )
    require(
        all(support == tuple(range(P))
            for support in digit_supports["guard_free"]),
        "guard-free delayed word lost a digit",
    )
    return prefixes, word_masses, digit_supports, dict(delayed_totals)


def rebuild_cospan_incidence_bank():
    """Refine all admissible sharp-graph events by both guard sectors."""
    (module, inherited_safe, _, _, rails,
     present, starts) = cross.build_carrier_data()
    require(len(rails) == 162, "THM-2616 middle-rail census changed")
    prefixes, word_masses, digit_supports, delayed_totals = (
        build_guard_cospan(module, inherited_safe)
    )

    matrices = {
        sector: {
            (d, s, ell): sharp.zero_matrix()
            for d in sharp.CLOCK_STEPS
            for s in range(1, P)
            for ell in range(Q7)
        }
        for sector in SECTORS
    }
    positives = Counter()
    positives_by_h = {sector: Counter() for sector in SECTORS}
    candidates = 0
    atom_partition_checks = 0
    graph_partition_checks = 0
    hole_candidates = 0
    hole_danger_positives = 0

    for s, ell, rail_digit, pieces in rails:
        require(rail_digit in (0, 12), "middle rail left its two-digit chart")
        for d in sharp.CLOCK_STEPS:
            future_clock = (ell + d) % Q7
            for h in ADMISSIBLE_GRAPH_DIGITS:
                r = (-h - 1) % P
                require(r != 0, "admissible graph entered r=0")
                overlap = old.intersect_weighted_union(
                    pieces,
                    present[future_clock, (-h) % P],
                    starts[future_clock, (-h) % P],
                )
                whole = old.intersect_weighted_comb(
                    overlap, module.C3, 182, 14 * r - 13, 14 * r + 13
                )
                whole_values = {
                    sector: old.delayed_weighted_numerator(
                        whole, prefixes[sector][future_clock][h]
                    )
                    for sector in SECTORS
                }
                require(
                    whole_values["safe"] + whole_values["danger"]
                    == whole_values["guard_free"],
                    "whole graph event violates the guard partition",
                )
                split_totals = Counter()

                for epsilon, left, right in (
                    (0, 14 * r, 14 * r + 13),
                    (1, 14 * r - 13, 14 * r),
                ):
                    half_tooth = old.intersect_weighted_comb(
                        overlap, module.C3, 182, left, right
                    )
                    for kappa in (0, 1):
                        candidates += 1
                        in_d1_hole = (
                            d == 1 and (s, ell) in sharp.EXPECTED_D1_ZERO_CELLS
                        )
                        hole_candidates += int(in_d1_hole)
                        j = INV2 * (r - epsilon - kappa) % P
                        require(
                            h == (-2 * j - kappa - epsilon - 1) % P,
                            "sharp graph/predecessor equation failed",
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
                        values = {
                            sector: old.delayed_weighted_numerator(
                                carry_piece,
                                prefixes[sector][future_clock][h],
                            )
                            for sector in SECTORS
                        }
                        require(
                            values["safe"] + values["danger"]
                            == values["guard_free"],
                            "private atom violates the guard partition",
                        )
                        atom_partition_checks += 1
                        for sector in SECTORS:
                            split_totals[sector] += values[sector]
                            if values[sector]:
                                positives[sector] += 1
                                positives_by_h[sector][h] += 1
                                matrices[sector][d, s, ell][j][h] = 1
                        if in_d1_hole and values["danger"]:
                            hole_danger_positives += 1

                require(
                    all(split_totals[sector] == whole_values[sector]
                        for sector in SECTORS),
                    "private atoms do not partition a graph event",
                )
                graph_partition_checks += 1

    require(candidates == EXPECTED_CANDIDATES,
            "guard-cospan candidate universe changed")
    require(atom_partition_checks == EXPECTED_CANDIDATES,
            "atom partition-check universe changed")
    require(graph_partition_checks == EXPECTED_GRAPH_PARTITIONS,
            "graph partition-check universe changed")
    require(dict(positives) == EXPECTED_POSITIVE_ATOMS,
            "guard-sector positive census changed")
    require(
        all(dict(positives_by_h[sector]) == EXPECTED_POSITIVE_BY_H[sector]
            for sector in SECTORS),
        "guard-sector positive-by-digit census changed",
    )
    require(hole_candidates == 1_584,
            "canonical d=1 hole candidate universe changed")
    require(hole_danger_positives == 0,
            "danger sector unexpectedly entered a canonical d=1 hole")

    # The danger leg restores h=12 as a delayed digit, but its graph partner
    # r=0 is disjoint from the present c3-safe factor before any delayed word
    # is consulted.  Check that structural puncture on every rail and step.
    r0_window = module.make_comb(module.C3, 182, -13, 13)
    r0_checks = 0
    for _, ell, _, pieces in rails:
        for d in sharp.CLOCK_STEPS:
            future_clock = (ell + d) % Q7
            overlap = old.intersect_weighted_union(
                pieces,
                present[future_clock, (-12) % P],
                starts[future_clock, (-12) % P],
            )
            require(not old.intersect_weighted_union(overlap, r0_window),
                    "present c3-safe factor met the r=0 danger sheet")
            r0_checks += 1
    require(r0_checks == EXPECTED_R0_CHECKS,
            "structural r=0 check universe changed")

    # The guard-free matrix must be the Boolean union of the two labelled
    # sector matrices, while the labels themselves remain available above.
    for d in sharp.CLOCK_STEPS:
        for s in range(1, P):
            for ell in range(Q7):
                for j in range(P):
                    for h in range(P):
                        require(
                            matrices["guard_free"][d, s, ell][j][h]
                            == int(
                                matrices["safe"][d, s, ell][j][h]
                                or matrices["danger"][d, s, ell][j][h]
                            ),
                            "guard-free support is not safe union danger",
                        )

    return (
        matrices, candidates, atom_partition_checks, graph_partition_checks,
        r0_checks, positives, positives_by_h, hole_candidates,
        hole_danger_positives, word_masses, digit_supports, delayed_totals,
    )


def constant_offset_candidate_controls(matrices):
    """Compute zero incidences and fixed-source constant-offset candidates."""
    zero_cells = {sector: {} for sector in SECTORS}
    universally_empty = {sector: {} for sector in SECTORS}
    product_histograms = {sector: {} for sector in SECTORS}
    zero_products = Counter()
    twisted_returns = Counter()

    for sector in SECTORS:
        for d in sharp.CLOCK_STEPS:
            zero_cells[sector][d] = tuple(
                (s, ell)
                for s in range(1, P)
                for ell in range(Q7)
                if sharp.support_size(matrices[sector][d, s, ell]) == 0
            )
            require(
                len(zero_cells[sector][d])
                == sharp.EXPECTED_ZERO_CELL_COUNTS[d],
                f"{sector} offset-{d} zero-incidence census changed",
            )
            universally_empty[sector][d] = tuple(
                ell4 for ell4 in range(Q7)
                if all(sharp.support_size(
                    matrices[sector][d, s, ell4]
                ) == 0 for s in range(1, P))
            )
            require(
                universally_empty[sector][d]
                == sharp.EXPECTED_UNIVERSALLY_EMPTY_ELL4[d],
                f"{sector} offset-{d} universal empty-ell4 set changed",
            )
            order = tuple((i * d) % Q7 for i in range(Q7))
            histogram = Counter()
            for s in range(1, P):
                product = sharp.identity_matrix()
                for ell in order:
                    product = sharp.boolean_product(
                        product, matrices[sector][d, s, ell]
                    )
                product_support = sharp.support_size(product)
                histogram[product_support] += 1
                zero_products[sector] += int(product_support == 0)
                twisted_returns[sector] += sum(
                    product[j][(j - 7 * a) % P]
                    for a in range(1, P)
                    for j in range(P)
                )
            require(histogram == Counter({0: 12}),
                    f"{sector} offset-{d} candidate acquired support")
            product_histograms[sector][d] = histogram

        require(zero_cells[sector][1] == sharp.EXPECTED_D1_ZERO_CELLS,
                f"{sector} canonical zero atlas changed")
        require(zero_products[sector] == 72,
                f"{sector} does not have 72 zero constant-offset candidates")
        require(twisted_returns[sector] == 0,
                f"{sector} acquired a twisted return")

    require(
        all(zero_cells["safe"][d] == zero_cells["danger"][d]
            == zero_cells["guard_free"][d]
            for d in sharp.CLOCK_STEPS),
        "the delayed guard sector changed the guard-invariant zero atlas",
    )
    return (
        zero_cells, universally_empty,
        product_histograms, zero_products, twisted_returns,
    )


def quotient_hostiles(matrices):
    """Show that varying offsets or forgetting ell4 destroys a broad no-go."""
    safe_varying, safe_best, safe_marginal = sharp.broader_no_go_hostiles(
        matrices["safe"]
    )
    require(
        all(safe_varying[s] == 70 for s in sharp.GENERIC_S)
        and safe_best[1][0] == 88,
        "safe-sector hostiles disagree with the base incidence audit",
    )

    varying_support = {sector: {} for sector in SECTORS}
    displayed_best_support = {sector: {} for sector in SECTORS}
    marginal_support = {sector: {} for sector in SECTORS}
    for sector in SECTORS:
        for s in range(1, P):
            varying_support[sector][s] = sharp.support_size(
                sharp.itinerary_product(
                    matrices[sector], s, sharp.VARYING_OFFSET_ITINERARY
                )
            )
            displayed_best_support[sector][s] = sharp.support_size(
                sharp.itinerary_product(
                    matrices[sector], s, sharp.BEST_NORMALIZED_ITINERARY
                )
            )
        for d in sharp.CLOCK_STEPS:
            for s in range(1, P):
                marginal = sharp.boolean_union(
                    matrices[sector][d, s, ell4]
                    for ell4 in range(Q7)
                )
                marginal_support[sector][d, s] = sharp.support_size(
                    sharp.boolean_power(marginal, Q7)
                )

    require(
        all(varying_support["guard_free"][s] == 104
            for s in sharp.GENERIC_S)
        and varying_support["guard_free"][6] == 0
        and varying_support["guard_free"][7] == 0,
        "guard-free varying-offset hostile changed",
    )
    require(
        all(displayed_best_support["guard_free"][s] == 143
            for s in sharp.GENERIC_S)
        and displayed_best_support["guard_free"][6] == 0
        and displayed_best_support["guard_free"][7] == 0,
        "guard-free displayed-itinerary hostile changed",
    )
    require(
        all(value == 0 for value in varying_support["danger"].values())
        and all(value == 0
                for value in displayed_best_support["danger"].values()),
        "danger-only itinerary hostile changed",
    )
    require(
        Counter(safe_marginal.values()) == Counter({110: 71, 75: 1})
        and safe_marginal[5, 6] == 75,
        "safe clock-marginal hostile changed",
    )
    require(
        Counter(marginal_support["danger"].values()) == Counter({0: 72}),
        "danger-only clock marginal acquired support",
    )
    require(
        Counter(marginal_support["guard_free"].values())
        == Counter({156: 72}),
        "guard-free clock-marginal hostile changed",
    )
    for d in sharp.CLOCK_STEPS:
        for s in range(1, P):
            marginal = sharp.boolean_union(
                matrices["guard_free"][d, s, ell4]
                for ell4 in range(Q7)
            )
            seventh = sharp.boolean_power(marginal, Q7)
            require(
                all(
                    seventh[j][h] == int(h != 12)
                    for j in range(P) for h in range(P)
                ),
                "guard-free marginal seventh power lost its exact rectangle",
            )
    return varying_support, displayed_best_support, marginal_support


def sorted_counter(counter):
    return tuple(sorted(counter.items()))


def main():
    ambient_sizes, _, ambient_twists = sharp.ambient_relation_controls()
    (
        matrices, candidates, atom_partitions, graph_partitions, r0_checks,
        positives, positives_by_h, hole_candidates, hole_danger_positives,
        word_masses, digit_supports, delayed_totals,
    ) = rebuild_cospan_incidence_bank()
    (
        zero_cells, universally_empty,
        product_hists, zero_products, twisted_returns,
    ) = constant_offset_candidate_controls(matrices)
    varying_support, displayed_best_support, marginal_support = (
        quotient_hostiles(matrices)
    )

    print("LRC14 exact guard-cospan successor-private incidence audit")
    print("scope=canonical raw nonnegative carrier; delayed guard split retained before integration")
    print("typing=common-x incidence at imposed (ell4,ell5=ell4+d); not a proved transition")
    print("graph=q=h, r=-h-1; admissible h=0,...,11; h=12->r=0 checked empty")
    print(
        "guard_word_grid_masses="
        + str(tuple((sector, word_masses[sector]) for sector in SECTORS))
    )
    print(
        "delayed_prefix_total_masses="
        + str(tuple((sector, delayed_totals[sector]) for sector in SECTORS))
    )
    print(f"safe_digit_support_clock0={digit_supports['safe'][0]}")
    print("safe_digit_support_clocks1to6=(1,...,11)")
    print(f"danger_digit_support_all_clocks={digit_supports['danger'][0]}")
    print("guard_free_digit_support_all_clocks=(0,...,12)")
    print(
        f"candidates={candidates} atom_partition_checks={atom_partitions} "
        f"graph_partition_checks={graph_partitions} structural_r0_checks={r0_checks}"
    )
    print(
        "positive_labelled_atom_occurrences="
        + str(tuple((sector, positives[sector]) for sector in SECTORS))
    )
    print(
        "positive_atoms_by_h="
        + str(tuple(
            (sector, sorted_counter(positives_by_h[sector]))
            for sector in SECTORS
        ))
    )
    print(
        f"canonical_d1_hole_candidates={hole_candidates} "
        f"danger_positive={hole_danger_positives} filled_cells=0/18"
    )
    print(
        "common_zero_incidence_cells_by_offset="
        + str(tuple(
            (d, len(zero_cells["guard_free"][d]))
            for d in sharp.CLOCK_STEPS
        ))
    )
    print(
        "common_universally_empty_ell4_by_offset="
        + str(tuple(
            (d, universally_empty["guard_free"][d])
            for d in sharp.CLOCK_STEPS
        ))
    )
    print(f"common_d1_zero_incidence_cells={zero_cells['guard_free'][1]}")
    print(
        "fixed_source_constant_offset_zero_candidates="
        + str(tuple((sector, zero_products[sector]) for sector in SECTORS))
    )
    print(
        "fixed_source_constant_offset_support_hist_by_sector="
        + str(tuple(
            (sector, tuple(
                (d, sorted_counter(product_hists[sector][d]))
                for d in sharp.CLOCK_STEPS
            ))
            for sector in SECTORS
        ))
    )
    print(
        "fixed_source_constant_offset_twisted_entries="
        + str(tuple((sector, twisted_returns[sector]) for sector in SECTORS))
    )
    print(
        f"varying_offset_itinerary={sharp.VARYING_OFFSET_ITINERARY} "
        + "generic_supports="
        + str(tuple(
            (sector, varying_support[sector][1]) for sector in SECTORS
        ))
        + " exceptional_s6_s7=0"
    )
    print(
        f"displayed_normalized_itinerary={sharp.BEST_NORMALIZED_ITINERARY} "
        + "generic_supports="
        + str(tuple(
            (sector, displayed_best_support[sector][1]) for sector in SECTORS
        ))
        + " exceptional_s6_s7=0"
    )
    print(
        "clock_marginal_seventh_power_support_hist="
        + str(tuple(
            (sector, sorted_counter(Counter(
                marginal_support[sector].values()
            )))
            for sector in SECTORS
        ))
        + " safe_exception=(d,s,support)=(5,6,75)"
    )
    print(
        "guard_free_clock_marginal_seventh_power_shape="
        "F13x(F13\\{12}) for all 72 (d,s)"
    )
    print(
        f"ambient_relation_power_supports={ambient_sizes} "
        f"power4_twisted_returns={ambient_twists}/156"
    )
    print(
        "verdict=PASS: danger restores h=0 but fills 0/18 fixed-offset "
        "incidence holes; only fixed-source constant-offset candidates vanish"
    )
    print(
        "semantics=the incidence zeros are invariant under the delayed guard split, "
        "but varying offsets and ell4-marginalization survive; no chronology follows"
    )


if __name__ == "__main__":
    main()
