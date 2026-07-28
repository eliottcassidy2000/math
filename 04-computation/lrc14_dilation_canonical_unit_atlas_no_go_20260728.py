#!/usr/bin/env python3
"""Exact scout for the canonical-unit atlas under the numerical D handoff.

THM-2635 gives the state-independent canonical fixed-half unit labels

    U^0={9},  U^1={3,8,10}.

For the dilation D(x)={13x}, THM-2670 identifies the next predecessor digit
with the current future digit: j'=h.  This scout applies only that numerical
digit handoff, then rebuilds every possible child half/carry atom in the
complete THM-2616 carrier.  It retains the safe/danger/guard-free cospan and
counts both rail-labelled occurrences and distinct (d,s,ell) base cells.

The result is deliberately not called a fibre product.  No D-map on the
current rail, source, base cell, or component labels is constructed here, and
cell-dependent unit sets are not tested.  The exact conclusion is only that
the four-label *state-independent* canonical unit atlas is not closed under
the numerical D digit handoff.
"""

from collections import Counter

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old
import lrc14_guard_cospan_successor_private_clock_collapse as guard
import lrc14_successor_private_sharp_graph_clock_collapse as sharp


P = 13
Q7 = 7
CANONICAL_UNITS = {0: (9,), 1: (3, 8, 10)}
SOURCE_UNITS = tuple(
    (epsilon, h)
    for epsilon in (0, 1)
    for h in CANONICAL_UNITS[epsilon]
)
BITS = tuple((epsilon, kappa) for epsilon in (0, 1) for kappa in (0, 1))
SECTORS = guard.SECTORS
EXPECTED_RAILS = 162
EXPECTED_OCCURRENCE_CANDIDATES = EXPECTED_RAILS * len(sharp.CLOCK_STEPS)
EXPECTED_DISTINCT_CELLS = 12 * Q7 * len(sharp.CLOCK_STEPS)

# (source_epsilon, source_h, child_epsilon, child_kappa):
#     (child_r, child_h, positive rail occurrences, positive base cells)
EXPECTED = {
    (0, 9, 0, 0): (5, 7, 0, 0),
    (0, 9, 0, 1): (6, 6, 301, 301),
    (0, 9, 1, 0): (6, 6, 301, 301),
    (0, 9, 1, 1): (7, 5, 0, 0),
    (1, 3, 0, 0): (6, 6, 290, 290),
    (1, 3, 0, 1): (7, 5, 0, 0),
    (1, 3, 1, 0): (7, 5, 367, 347),
    (1, 3, 1, 1): (8, 4, 0, 0),
    (1, 8, 0, 0): (3, 9, 0, 0),
    (1, 8, 0, 1): (4, 8, 367, 347),
    (1, 8, 1, 0): (4, 8, 0, 0),
    (1, 8, 1, 1): (5, 7, 301, 301),
    (1, 10, 0, 0): (7, 5, 367, 347),
    (1, 10, 0, 1): (8, 4, 0, 0),
    (1, 10, 1, 0): (8, 4, 345, 325),
    (1, 10, 1, 1): (9, 3, 0, 0),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def child_labels(predecessor, epsilon, kappa):
    """THM-2670's successor-private affine law on printed digit labels."""
    r = (2 * predecessor + epsilon + kappa) % P
    h = (-r - 1) % P
    return r, h


def half_compatible(h, kappa):
    """Elementary intersection condition for the successor half-bit."""
    return (kappa == 0 and h <= 6) or (kappa == 1 and h >= 6)


def rebuild():
    module, inherited, _, _, rails, present, starts = cross.build_carrier_data()
    prefixes, _, _, _ = guard.build_guard_cospan(module, inherited)
    require(len(rails) == EXPECTED_RAILS, "THM-2616 rail census changed")
    require(module.C3 == 2 * sharp.PREDECESSOR_SPEED,
            "later tooth speed changed")
    require(old.R == sharp.SUCCESSOR_SPEED,
            "successor digit speed changed")

    occurrences = {sector: Counter() for sector in SECTORS}
    cells = {sector: {key: set() for key in EXPECTED} for sector in SECTORS}
    numerators = {sector: Counter() for sector in SECTORS}
    partition_checks = 0

    for s, ell, rail_digit, pieces in rails:
        require(rail_digit in (0, 12), "middle rail left its two-digit chart")
        for d in sharp.CLOCK_STEPS:
            future = (ell + d) % Q7
            cell = (d, s, ell)
            for source_epsilon, source_h in SOURCE_UNITS:
                predecessor = source_h  # numerical D handoff j'=h
                for child_epsilon, child_kappa in BITS:
                    key = (source_epsilon, source_h,
                           child_epsilon, child_kappa)
                    r, h = child_labels(
                        predecessor, child_epsilon, child_kappa
                    )
                    require(r != 0, "scout unexpectedly entered the r=0 chart")

                    overlap = old.intersect_weighted_union(
                        pieces,
                        present[future, (-h) % P],
                        starts[future, (-h) % P],
                    )
                    if child_epsilon == 0:
                        left, right = 14 * r, 14 * r + 13
                    else:
                        left, right = 14 * r - 13, 14 * r
                    atom = old.intersect_weighted_comb(
                        overlap, module.C3, 182, left, right
                    )
                    atom = old.intersect_weighted_comb(
                        atom, sharp.PREDECESSOR_SPEED, P,
                        predecessor, predecessor + 1,
                    )
                    atom = old.intersect_weighted_comb(
                        atom, sharp.SUCCESSOR_SPEED, 2,
                        child_kappa, child_kappa + 1,
                    )
                    values = {
                        sector: old.delayed_weighted_numerator(
                            atom, prefixes[sector][future][h]
                        )
                        for sector in SECTORS
                    }
                    require(
                        values["safe"] + values["danger"]
                        == values["guard_free"],
                        "guard cospan failed on a D-handoff atom",
                    )
                    partition_checks += 1
                    for sector, value in values.items():
                        numerators[sector][key] += value
                        if value > 0:
                            occurrences[sector][key] += 1
                            cells[sector][key].add(cell)

    expected_checks = (
        EXPECTED_OCCURRENCE_CANDIDATES * len(SOURCE_UNITS) * len(BITS)
    )
    require(partition_checks == expected_checks,
            "D-handoff atom universe changed")
    require(len({(d, s, ell) for d in sharp.CLOCK_STEPS
                 for s in range(1, P) for ell in range(Q7)})
            == EXPECTED_DISTINCT_CELLS,
            "base-cell universe changed")

    observed = {}
    for key in EXPECTED:
        source_epsilon, source_h, child_epsilon, child_kappa = key
        r, h = child_labels(source_h, child_epsilon, child_kappa)
        observed[key] = (
            r, h,
            occurrences["guard_free"][key],
            len(cells["guard_free"][key]),
        )
    require(observed == EXPECTED, "canonical D-handoff table changed")

    # The danger leg contributes nothing, and safe equals guard-free on every
    # row at both the positivity and exact-numerator levels.
    require(all(occurrences["danger"][key] == 0 for key in EXPECTED),
            "danger sector acquired a positive handoff atom")
    require(all(numerators["danger"][key] == 0 for key in EXPECTED),
            "danger sector acquired nonzero exact mass")
    require(all(occurrences["safe"][key]
                == occurrences["guard_free"][key] for key in EXPECTED),
            "safe and guard-free occurrence counts differ")
    require(all(cells["safe"][key] == cells["guard_free"][key]
                for key in EXPECTED),
            "safe and guard-free positive cells differ")

    # Every positive exact atom obeys the elementary half/carry constraint.
    for key, (_, h, count, _) in observed.items():
        child_kappa = key[3]
        require(count == 0 or half_compatible(h, child_kappa),
                "positive atom violates half compatibility")

    unit_children = tuple(
        key for key, (_, h, _, _) in observed.items()
        if h in CANONICAL_UNITS[key[2]]
    )
    require(unit_children == (
        (1, 8, 0, 0),
        (1, 8, 1, 0),
        (1, 10, 1, 1),
    ), "formula-compatible canonical-unit child list changed")
    require(all(observed[key][2:] == (0, 0) for key in unit_children),
            "a canonical uniform-unit child became positive")
    require(all(not half_compatible(observed[key][1], key[3])
                for key in unit_children),
            "unit-child zero lost its half-compatibility explanation")

    positive_nonunits = tuple(
        key for key, (_, h, count, _) in observed.items()
        if count > 0 and h not in CANONICAL_UNITS[key[2]]
    )
    require(len(positive_nonunits) == 8,
            "positive nonunit row count changed")
    return observed, unit_children, positive_nonunits, partition_checks


def main():
    observed, unit_children, positive_nonunits, partition_checks = rebuild()
    print("LRC14 exact canonical-unit atlas / numerical D-handoff scout")
    print("scope=state-independent THM-2635 unit labels under j'=h only")
    print("not_tested=physical D fibre product; rail/source/cell/component transport; cell-dependent units")
    print(f"canonical_units=U0{CANONICAL_UNITS[0]} U1{CANONICAL_UNITS[1]}")
    print(
        f"candidate_universe_per_source_and_bit={EXPECTED_OCCURRENCE_CANDIDATES} "
        f"distinct_base_cells={EXPECTED_DISTINCT_CELLS} "
        f"partition_checks={partition_checks}"
    )
    print("table_columns=(source_epsilon,source_h,child_epsilon,child_kappa,r,h,positive_occurrences,positive_cells)")
    for key in EXPECTED:
        print(f"row={key + observed[key]}")
    print("danger_positive_rows=0; safe_equals_guard_free=exact")
    print(f"formula_compatible_uniform_unit_children={unit_children}")
    print(f"positive_nonunit_rows={positive_nonunits}")
    print("mechanism=every possible uniform-unit child violates kappa=0=>h<=6, kappa=1=>h>=6")
    print("verdict=PASS: the canonical uniform unit atlas is not closed under the numerical D digit handoff")
    print("boundary=this is not an emptiness theorem for the physical two-edge D fibre product")


if __name__ == "__main__":
    main()
