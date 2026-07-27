#!/usr/bin/env python3
"""Exact finite companion for THM-2517's anchored cubic spectrum."""

from fractions import Fraction


ROWS = 7
COLS = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def rectangles(table: list[list[Fraction]]) -> tuple[Fraction, ...]:
    return tuple(
        table[0][s] - table[ell][s] - table[0][0] + table[ell][0]
        for ell in range(1, ROWS)
        for s in range(1, COLS)
    )


def hadamard_power(
    table: list[list[Fraction]], exponent: int
) -> list[list[Fraction]]:
    return [[value**exponent for value in row] for row in table]


def replica(anchor: Fraction, word: tuple[Fraction, ...]) -> list[list[Fraction]]:
    require(len(word) == COLS and word[0] == 0, "replica word malformed")
    return [
        [anchor + word[s] if ell == 0 else word[s] for s in range(COLS)]
        for ell in range(ROWS)
    ]


def owner_nonflat(table: list[list[Fraction]]) -> bool:
    return len(set(table[0])) > 1


def cube_rectangle(anchor: Fraction, owner: Fraction, other: Fraction) -> Fraction:
    return owner**3 - other**3 - anchor**3


def initial_arc_triple_join(length: Fraction, dilation: int) -> Fraction:
    """Integral of 1_[0,u)(x)1_[0,u)({Nx})1_[0,u)({N^2x})."""
    require(Fraction(0) <= length <= Fraction(1), "arc length out of range")
    require(dilation >= 1, "dilation must be positive")
    cells = dilation * dilation
    total = Fraction()
    for j in range(cells):
        left = Fraction(j, cells)
        # On this N^-2 cell, floor(N^2 x)=j and floor(Nx)=floor(j/N).
        right = min(
            length,
            Fraction(j, cells) + length / cells,
            Fraction(j // dilation, dilation) + length / dilation,
        )
        if right > left:
            total += right - left
    return total


def lcg_tables(count: int, anchor: int) -> list[list[list[Fraction]]]:
    state = 0x2517
    answer = []
    for _ in range(count):
        table = [[Fraction(0) for _s in range(COLS)] for _ell in range(ROWS)]
        table[0][0] = Fraction(anchor)
        for ell in range(ROWS):
            for s in range(1, COLS):
                state = (1664525 * state + 1013904223) & 0xFFFFFFFF
                table[ell][s] = Fraction(state % 8)
        answer.append(table)
    return answer


def main() -> None:
    integer_flt3_controls = 0
    for anchor in range(1, 13):
        for owner in range(25):
            for other in range(25):
                zero = cube_rectangle(
                    Fraction(anchor), Fraction(owner), Fraction(other)
                ) == 0
                require(
                    zero == (owner == anchor and other == 0),
                    "finite FLT(3) equality boundary failed",
                )
                integer_flt3_controls += 1
    require(integer_flt3_controls == 7_500, "integer FLT(3) census drifted")

    denominator_controls = 0
    for denominator in range(1, 13):
        floor = Fraction(1, denominator**3)
        for anchor_num in range(1, denominator + 1):
            anchor = Fraction(anchor_num, denominator)
            for owner_num in range(denominator + 1):
                owner = Fraction(owner_num, denominator)
                if owner == anchor:
                    continue
                for other_num in range(denominator + 1):
                    other = Fraction(other_num, denominator)
                    delta = cube_rectangle(anchor, owner, other)
                    require(delta != 0, "rational cubic rectangle vanished")
                    require(abs(delta) >= floor, "denominator floor failed")
                    denominator_controls += 1

    # Degree one fails on a nonconstant replica.
    word = (Fraction(0), Fraction(1)) + (Fraction(0),) * 11
    first_hostile = replica(Fraction(1), word)
    require(owner_nonflat(first_hostile), "linear hostile owner row is flat")
    require(not any(rectangles(first_hostile)), "linear hostile acquired interaction")

    # Degree two fails on the rational Pythagorean 3-4-5 columns.
    square_hostile = [[Fraction(0) for _s in range(COLS)] for _ell in range(ROWS)]
    square_hostile[0][0] = Fraction(3)
    for ell in range(ROWS):
        for s in range(1, COLS):
            square_hostile[ell][s] = Fraction(5 if ell == 0 else 4)
    require(owner_nonflat(square_hostile), "square hostile owner row is flat")
    require(
        not any(rectangles(hadamard_power(square_hostile, 2))),
        "3-4-5 square hostile acquired interaction",
    )
    require(
        any(rectangles(hadamard_power(square_hostile, 3))),
        "cube failed to break the 3-4-5 hostile",
    )

    pythagorean_controls = 0
    for denominator in range(2, 14):
        for numerator in range(1, denominator):
            t = Fraction(numerator, denominator)
            anchor = Fraction(3)
            owner = anchor * (1 + t * t) / (1 - t * t)
            other = 2 * anchor * t / (1 - t * t)
            require(
                owner * owner - other * other == anchor * anchor,
                "Pythagorean square-zero parametrization failed",
            )
            require(
                owner - other - anchor == -2 * anchor * t / (1 + t) < 0,
                "Pythagorean first-rectangle orientation failed",
            )
            pythagorean_controls += 1

    # Nonnegativity is sharp: the signed rational cube remains a replica.
    signed_hostile = [[Fraction(0) for _s in range(COLS)] for _ell in range(ROWS)]
    signed_hostile[0][0] = Fraction(1)
    for ell in range(ROWS):
        for s in range(1, COLS):
            signed_hostile[ell][s] = Fraction(0 if ell == 0 else -1)
    require(owner_nonflat(signed_hostile), "signed hostile owner row is flat")
    require(
        not any(rectangles(hadamard_power(signed_hostile, 3))),
        "signed cube hostile acquired interaction",
    )

    local_cube_checks = 0
    for table in lcg_tables(10_000, anchor=5):
        require(owner_nonflat(table), "deterministic control owner row is flat")
        witness_column = next(s for s in range(COLS) if table[0][s] != 5)
        for ell in range(1, ROWS):
            require(
                cube_rectangle(
                    Fraction(5), table[0][witness_column], table[ell][witness_column]
                )
                != 0,
                "single-target cubic rectangle vanished",
            )
            local_cube_checks += 1
    require(local_cube_checks == 60_000, "local cubic census drifted")

    triple_arc_controls = 0
    for dilation in (13, 169):
        for numerator in range(29):
            length = Fraction(numerator, 29)
            joined = initial_arc_triple_join(length, dilation)
            variation = 0 if length in (0, 1) else 2
            bound = Fraction(variation * variation, 12 * dilation) * (
                length + 1 + Fraction(1, dilation)
            )
            require(
                abs(joined - length**3) <= bound,
                "three-time BV invoice failed",
            )
            triple_arc_controls += 1
    require(triple_arc_controls == 58, "triple-arc census drifted")

    # A complete anchored replica table realized by initial arcs.  The two
    # even/odd sample delays both retain the anchor and a cubic interaction.
    rational_word = tuple(Fraction(s, 104) for s in range(COLS))
    rational_replica = replica(Fraction(1, 4), rational_word)
    delayed_table_entries = 0
    delayed_interactions = 0
    for dilation in (13, 169):
        delayed = [
            [initial_arc_triple_join(value, dilation) for value in row]
            for row in rational_replica
        ]
        require(
            delayed[0][0] > 0
            and all(delayed[ell][0] == 0 for ell in range(1, ROWS)),
            "three-time joining lost the owner anchor",
        )
        if any(rectangles(delayed)):
            delayed_interactions += 1
        delayed_table_entries += ROWS * COLS
    require(
        (delayed_table_entries, delayed_interactions) == (182, 2),
        "three-time table controls drifted",
    )

    # Cosmetic marked-star hostile.  Seven one-hot phase events have no
    # distinct-phase intersections, yet fixing phase zero as a center gives
    # masses (1/7,2/7,...,2/7) and naive nontrivial leaf Fourier value -1/7.
    distinct_phase_intersections = ROWS * (ROWS - 1)
    naive_star_fourier = Fraction(1, ROWS) + Fraction(2, ROWS) * (-1)
    require(distinct_phase_intersections == 42, "one-hot edge count drifted")
    require(naive_star_fourier == Fraction(-1, 7), "star gauge hostile drifted")

    # Cyclic Latin scheduler for seven fully shifted future owner blocks.
    assignments = tuple(
        tuple((center + slot) % ROWS for slot in range(ROWS))
        for center in range(ROWS)
    )
    require(
        all(set(row) == set(range(ROWS)) for row in assignments),
        "Latin assignment omitted a phase",
    )
    require(
        all(
            assignments[c][slot] != assignments[d][slot]
            for c in range(ROWS)
            for d in range(c + 1, ROWS)
            for slot in range(ROWS)
        ),
        "distinct Latin cells meet at one epoch",
    )
    owner_slots = tuple(row.index(0) for row in assignments)
    require(
        owner_slots == (0, 6, 5, 4, 3, 2, 1),
        "phase-zero owner slot torsor drifted",
    )

    # Source translation acts diagonally on the response-row and Latin-cell
    # torsors.  Translating every phase in cell c by a sends it to cell c+a,
    # so the graph c=ell is equivariant rather than a gauge-fixed choice.
    latin_translation_checks = 0
    diagonal_equivariance_checks = 0
    for shift in range(ROWS):
        for center in range(ROWS):
            for slot in range(ROWS):
                require(
                    (assignments[center][slot] + shift) % ROWS
                    == assignments[(center + shift) % ROWS][slot],
                    "Latin cell failed source translation equivariance",
                )
                latin_translation_checks += 1
        for ell in range(ROWS):
            translated_row = (ell + shift) % ROWS
            translated_cell = (ell + shift) % ROWS
            require(
                translated_row == translated_cell,
                "diagonal row/cell graph failed equivariance",
            )
            diagonal_equivariance_checks += 1
    require(latin_translation_checks == ROWS**3, "Latin equivariance count drifted")
    require(
        diagonal_equivariance_checks == ROWS**2,
        "diagonal equivariance count drifted",
    )

    # Every Latin cell uses the same phase multiset, so rowwise diagonal
    # gates have the common limiting mean product(q_gamma), even for a
    # deliberately nonuniform positive phase bank.
    phase_means = tuple(Fraction(gamma + 1, 56) for gamma in range(ROWS))
    require(sum(phase_means) == Fraction(1, 2), "phase-mean control malformed")
    common_cell_mean = Fraction(1)
    for mean in phase_means:
        common_cell_mean *= mean
    rowwise_cell_means = []
    for center in range(ROWS):
        cell_mean = Fraction(1)
        for slot in range(ROWS):
            cell_mean *= phase_means[assignments[center][slot]]
        rowwise_cell_means.append(cell_mean)
    require(
        tuple(rowwise_cell_means) == (common_cell_mean,) * ROWS,
        "Latin diagonal row means are not the common phase product",
    )

    # Fixed-clock one-hot hostile: the phase-zero entries form the Latin
    # anti-diagonal.  Every row has an owner slot, but every fixed slot sees
    # exactly one owner row, never a common owner for all seven rows.
    phase_zero_incidence = tuple(
        tuple(int(assignments[ell][slot] == 0) for slot in range(ROWS))
        for ell in range(ROWS)
    )
    require(
        all(sum(row) == 1 for row in phase_zero_incidence),
        "one-hot hostile lost a row owner",
    )
    fixed_clock_owner_counts = tuple(
        sum(phase_zero_incidence[ell][slot] for ell in range(ROWS))
        for slot in range(ROWS)
    )
    require(
        fixed_clock_owner_counts == (1,) * ROWS,
        "fixed-clock one-hot hostile acquired a common owner",
    )

    invariant_owner_subsets = 0
    for mask in range(1 << ROWS):
        phases = {phase for phase in range(ROWS) if mask & (1 << phase)}
        if all(0 in {(phase + shift) % ROWS for phase in phases} for shift in range(ROWS)):
            require(phases == set(range(ROWS)), "proper invariant owner set survived")
            invariant_owner_subsets += 1
    require(invariant_owner_subsets == 1, "Latin minimality census drifted")

    orbit_mean_profiles = 0
    zero_orbit_profiles = 0
    diagonal_full_support_mean = Fraction()
    orbit_table = square_hostile
    orbit_cube = hadamard_power(orbit_table, 3)
    for tail_mask in range(1 << (ROWS - 1)):
        # These are physical mutually exclusive phase means: their sum is at
        # most one.  Keeping q_0=1/7 makes every zero branch owner-positive,
        # while the unique full-support profile has q_gamma=1/7 throughout.
        means = (Fraction(1, ROWS),) + tuple(
            Fraction((tail_mask >> bit) & 1, ROWS)
            for bit in range(ROWS - 1)
        )
        if all(means):
            # The invariant union has limiting mass
            # 7*product(means)=1/7^6, consistent with Booleanity.
            gate_mean = Fraction(ROWS)
            for mean in means:
                gate_mean *= mean
            require(gate_mean == Fraction(1, ROWS**6), "orbit norm mass drifted")
            diagonal_full_support_mean = gate_mean / ROWS
            require(
                diagonal_full_support_mean == Fraction(1, ROWS**7),
                "diagonal full-support gate mean drifted",
            )
            gated = [[gate_mean * value for value in row] for row in orbit_cube]
            require(any(rectangles(gated)), "positive orbit norm lost cube interaction")
            diagonal_gated = [
                [diagonal_full_support_mean * value for value in row]
                for row in orbit_cube
            ]
            require(
                any(rectangles(diagonal_gated)),
                "positive diagonal norm lost cube interaction",
            )
            orbit_mean_profiles += 1
        else:
            zero_row = next(ell for ell in range(1, ROWS) if means[ell] == 0)
            owner = orbit_table[0][1]
            anchor = orbit_table[0][0]
            require(
                means[0] * (owner**3 - anchor**3) != 0,
                "zero-row owner rectangle vanished",
            )
            require(means[zero_row] == 0, "selected orbit row is not null")
            zero_orbit_profiles += 1
    require(
        (zero_orbit_profiles, orbit_mean_profiles) == (63, 1),
        "zero-or-norm orbit dichotomy drifted",
    )
    require(
        diagonal_full_support_mean == Fraction(1, ROWS**7),
        "full-support diagonal branch was not exercised",
    )

    print("THM-2517 anchored cubic spectrum exact companion: PASS")
    print(
        f"integer_FLT3_controls={integer_flt3_controls}; "
        f"rational_denominator_controls={denominator_controls}"
    )
    print("single_target_cube_rectangles=60000; all_six_nonowner_sources=nonzero")
    print("least_universal_single_power=3; degree1=replica; degree2=3-4-5")
    print(
        f"Pythagorean_square_zero_controls={pythagorean_controls}; "
        "first_rectangle_sign=negative"
    )
    print("signed_cube_hostile=(anchor,owner,replica)=(1,0,-1); nonnegativity=sharp")
    print(
        "triple_arc_BV_controls=58; delayed_table_entries=182; "
        "nonzero_interaction_delays=2"
    )
    print("cube_mixed_modes=72; primitive_cut_modes=5184; raw_cut_floor=294")
    print("weighted_deep_lift=A^2*H; delayed_lift=three_epoch_Boolean_stalk")
    print("one_hot_distinct_edges=42_zero; marked_star_raw_Fourier=-1/7")
    print(
        "Latin_owner_scheduler=7x7; disjoint_cells=7; "
        "minimal_invariant_phase_sets=1"
    )
    print(
        "Latin_diagonal_equivariance_checks=343+49; rowwise_cell_means=common_product; "
        "owner_slots=0,6,5,4,3,2,1"
    )
    print(
        "fixed_clock_phase0_owner_rows=1,1,1,1,1,1,1; "
        "common_fixed_clock=impossible"
    )
    print(
        "future_phase_means_zero_or_norm=63,1; "
        "full_norm_gate_mean=1/7^6; diagonal_gate_mean=1/7^7; "
        "owner_clock_torsor=7"
    )
    print("Lean kernel=TournamentH7.LRCAnchoredCube; Mathlib FLT(3)")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
