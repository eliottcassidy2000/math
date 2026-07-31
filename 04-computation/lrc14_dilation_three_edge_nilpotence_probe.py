#!/usr/bin/env python3
"""Exact three-edge nilpotence of the D-compatible LRC(14) carrier.

Let ``D(x)={13*x}``.  The typed two-edge companion proves that a candidate
handoff from

    E_i = E_(s_i; owner=ell_(i+1), shallow=ell_i; j_i,h_i)

to ``E_(i+1)`` is formed by evaluating the latter at ``D x``.  This transports
the current owner clock to the following shallow clock and gives
``h_i=j_(i+1)``.  Here we construct the genuine three-event object

    E_0(x) AND E_1(Dx) AND E_2(D^2 x).                    (1)

Every THM-2616 middle rail is contained in the raw arrival interval

    I_6=[6/13,7/13).

Put ``z=Dx`` and ``w=Dz=D^2x``.  The second and third rail factors in (1)
force ``z,w in I_6``.  The identity ``C1=13`` is load-bearing: the shallow
clock of ``E_i`` is the phase clock of ``D^(i+1)x``.  Thus the shallow clock
``ell_0`` of ``E_0`` labels ``z``, while the shallow clock ``ell_1`` of
``E_1`` labels ``w``.  Owner/shallow covariance also makes ``ell_1`` the
owner clock of ``E_0``.  Exact interval arithmetic gives

    I_6 intersect D^-1(I_6) = [84/169,85/169),

whose lower and upper halves have clock pairs only ``(3,3)`` and ``(4,4)``.
Thus ``ell_0=ell_1`` is necessary.  Every incidence edge under study has a
nonzero F_7 clock step, hence requires ``ell_0!=ell_1``.  The three-edge
fibre product is therefore empty before any target, sharp-root,
predecessor/carry-label,
guard, Bockstein, endpoint, or Perron restriction is consulted.

The script classifies all 7*6^3 admissible clock quadruples and all 12^3
source triples.  It also counts every inherited rail-labelled and nominal
predecessor/carry-atom triple eliminated by the same source-independent phase
gate.  A full exact positive two-edge numerator certificate is replayed as
the sharp control,
so the D-compatible support transfer has exact nilpotence index three, not
two, within this raw carrier model.

This is a support theorem for the stated D handoff and inherited arrival-six
middle rails.  It does not assert a Bockstein unit, endpoint section, scalar
cover, chronology outside this handoff, or LRC(14).
"""

from collections import Counter

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_dilation_reversed_clock_fibre_product_probe as two
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = 13
Q7 = 7
ARRIVAL_DIGIT = 6
LABELLED_ATOMS_PER_FIXED_RAIL_EDGE = 11 * 2 * 2

EXPECTED_RAIL_TRIPLE_MULTIPLICITY_HISTOGRAM = Counter({
    1: 900,
    2: 36_612,
    4: 483_948,
    8: 2_091_276,
})
EXPECTED_RAIL_LABELLED_TRIPLES = 18_740_124
EXPECTED_NOMINAL_LABELLED_ATOM_TRIPLES = 1_596_358_722_816
EXPECTED_TWO_EDGE_WITNESS = 126_816_337_986_097_204_341_478_787_325_120


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def interval_mass(intervals):
    return sum(right - left for left, right in intervals)


def exact_two_edge_positive_control(
        module, prefixes, rails, present, present_starts, grid):
    """Replay the first exact nonzero D-fibre witness from the companion."""
    # Current event: source 1, owner 1, shallow 3, rail index 2,
    # (j,h,epsilon,kappa)=(5,2,0,0).
    current = two.build_atom(
        module, prefixes, present, present_starts,
        rails[2], 3, 2, 0, 0,
    )
    # Following event at Dx: source 1, owner 0, shallow 1, rail index 0,
    # (j,h,epsilon,kappa)=(2,6,1,1).
    following = two.build_atom(
        module, prefixes, present, present_starts,
        rails[0], 1, 6, 1, 1,
    )
    require(rails[2][:3] == (1, 1, 12)
            and rails[0][:3] == (1, 0, 12),
            "positive-control rail metadata changed")
    require(
        tuple(current[name] for name in
              ("j", "h", "epsilon", "kappa")) == (5, 2, 0, 0)
        and tuple(following[name] for name in
                  ("j", "h", "epsilon", "kappa")) == (2, 6, 1, 1),
            "positive-control predecessor/carry labels changed",
    )
    require(current["h"] == following["j"],
            "D-compatible digit handoff failed in the positive control")

    base = two.intersect_weighted(
        two.scale_weighted(current["pieces"], P),
        two.pullback_dilation_weighted(following["pieces"], grid),
    )
    current_q = two.prefix_intervals(
        prefixes[current["future"]][current["h"]]
    )
    following_q_pullback = two.preimage_times_13(
        two.prefix_intervals(
            prefixes[following["future"]][following["h"]]
        ),
        grid,
    )
    joint_q = old.intersect_sorted(current_q, following_q_pullback)
    joint_prefix = module.make_prefix([
        (P * left, P * right) for left, right in joint_q
    ])
    value = two.delayed_weighted_numerator(
        base, joint_prefix, P**6, P * grid
    )
    require(len(base) == 23 and len(joint_q) == 66,
            "positive-control interval presentation changed")
    require(value == EXPECTED_TWO_EDGE_WITNESS,
            "exact two-edge positive numerator certificate changed")
    return value, len(base), len(joint_q)


def main():
    (module, prefixes, _, _, rails,
     present, present_starts) = cross.build_carrier_data()
    grid = old.T
    require(module.C1 == P,
            "the shallow clock is no longer the base-thirteen phase clock")
    require(grid % (2 * P * P) == 0,
            "base grid no longer resolves the nilpotence boundary")
    require(len(rails) == 162,
            "THM-2616 middle-rail census changed")

    arrival = [(
        ARRIVAL_DIGIT * grid // P,
        (ARRIVAL_DIGIT + 1) * grid // P,
    )]
    arrival_left, arrival_right = arrival[0]
    rail_piece_count = 0
    for _, _, _, pieces in rails:
        for left, right, weight in pieces:
            require(arrival_left <= left < right <= arrival_right,
                    "a middle-rail piece left the arrival-six interval")
            require(weight > 0, "a retained rail piece lost positivity")
            rail_piece_count += 1
    require(rail_piece_count == 1_050,
            "middle-rail weighted-piece census changed")

    # Check the owner/shallow covariance used to type every D handoff.
    owner_cells = old.base.clock_cells(P, Q7, grid, P * P)
    dilation_clock_covariance_checks = 0
    for ell in range(Q7):
        shallow = module.make_comb(
            module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        scaled_owner = two.merge_adjacent(
            (P * left, P * right) for left, right in owner_cells[ell]
        )
        pulled_shallow = two.merge_adjacent(
            (left + branch * grid, right + branch * grid)
            for left, right in shallow
            for branch in range(P)
        )
        require(scaled_owner == pulled_shallow,
                "D owner-to-shallow covariance changed")
        dilation_clock_covariance_checks += 1

    # The exact two-arrival phase cylinder for z=Dx and Dz=D^2x.
    arrival_return = old.intersect_sorted(
        arrival, two.preimage_times_13(arrival, grid)
    )
    expected_arrival_return = [(
        84 * grid // (P * P),
        85 * grid // (P * P),
    )]
    require(arrival_return == expected_arrival_return,
            "I_6 intersect D^-1(I_6) changed")
    require(interval_mass(arrival_return) == grid // (P * P),
            "two-arrival phase-cylinder mass changed")

    # Classify all shallow-clock pairs (ell_0,ell_1) on that cylinder.
    # Because C1=13, the shallow clock of E0 at x is the depth-one phase
    # clock ell_0 of z=Dx, and the shallow clock of E1 at z is the phase
    # clock ell_1 of w=Dz.  The latter is also E0's owner label by the
    # checked owner/shallow covariance.  Verify the first identification as
    # an exact interval equality before using it.  Exact support on (z,w) is
    # diagonal and split at z=1/2.
    phase_clock_cells = old.base.clock_cells(P, Q7, grid, 1)
    shallow_phase_identification_checks = 0
    for ell in range(Q7):
        shallow = module.make_comb(
            module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        phase_pullback = two.preimage_times_13(
            phase_clock_cells[ell], grid
        )
        require(
            two.merge_adjacent(shallow)
            == two.merge_adjacent(phase_pullback),
            "shallow clock did not become the next-point phase clock",
        )
        shallow_phase_identification_checks += 1
    phase_pair_masses = {}
    phase_pair_checks = 0
    for ell0 in range(Q7):
        first = old.intersect_sorted(
            arrival_return, phase_clock_cells[ell0]
        )
        for ell1 in range(Q7):
            support = old.intersect_sorted(
                first,
                two.preimage_times_13(phase_clock_cells[ell1], grid),
            )
            mass = interval_mass(support)
            if mass:
                phase_pair_masses[ell0, ell1] = mass
            phase_pair_checks += 1
    require(phase_pair_checks == Q7 * Q7,
            "phase clock-pair universe changed")
    require(phase_pair_masses == {
        (3, 3): grid // 338,
        (4, 4): grid // 338,
    }, "three-edge phase trap lost its exact diagonal support")
    require(sum(phase_pair_masses.values()) == grid // (P * P),
            "diagonal phase pieces do not partition the arrival return")

    # Classify the complete nonrepeating clock/source/rail universe.  A
    # three-event chain has four clock vertices.  Its first edge already has
    # ell_0 != ell_1, while the phase gate above requires equality.
    rails_per_cell = Counter((s, ell) for s, ell, _, _ in rails)
    require(Counter(rails_per_cell.values()) == Counter({2: 78, 1: 6}),
            "one/two-rail cell census changed")
    clock_quadruples = 0
    clock_source_candidates = 0
    rail_labelled_triples = 0
    rail_multiplicity_histogram = Counter()
    surviving_clock_source_candidates = 0
    for ell0 in range(Q7):
        for ell1 in range(Q7):
            if ell0 == ell1:
                continue
            require((ell0, ell1) not in phase_pair_masses,
                    "a nonzero first clock edge survived the phase trap")
            for ell2 in range(Q7):
                if ell1 == ell2:
                    continue
                for ell3 in range(Q7):
                    if ell2 == ell3:
                        continue
                    clock_quadruples += 1
                    for s0 in range(1, P):
                        for s1 in range(1, P):
                            for s2 in range(1, P):
                                multiplicity = (
                                    rails_per_cell[s0, ell1]
                                    * rails_per_cell[s1, ell2]
                                    * rails_per_cell[s2, ell3]
                                )
                                require(multiplicity in (1, 2, 4, 8),
                                        "rail-triple multiplicity changed")
                                clock_source_candidates += 1
                                rail_labelled_triples += multiplicity
                                rail_multiplicity_histogram[multiplicity] += 1
                                # The source-independent phase support is
                                # already empty, so no candidate survives.
                                surviving_clock_source_candidates += 0

    require(clock_quadruples == Q7 * 6**3 == 1_512,
            "nonrepeating clock-quadruple census changed")
    require(clock_source_candidates == clock_quadruples * 12**3
            == 2_612_736,
            "clock/source triple census changed")
    require(rail_multiplicity_histogram
            == EXPECTED_RAIL_TRIPLE_MULTIPLICITY_HISTOGRAM,
            "rail-triple multiplicity histogram changed")
    require(rail_labelled_triples == EXPECTED_RAIL_LABELLED_TRIPLES,
            "rail-labelled three-event universe changed")
    nominal_labelled_atom_triples = (
        rail_labelled_triples * LABELLED_ATOMS_PER_FIXED_RAIL_EDGE**3
    )
    require(nominal_labelled_atom_triples
            == EXPECTED_NOMINAL_LABELLED_ATOM_TRIPLES,
            "nominal predecessor/carry-atom triple universe changed")
    require(surviving_clock_source_candidates == 0,
            "a three-event D-compatible candidate survived")

    witness, witness_base_pieces, witness_q_pieces = (
        exact_two_edge_positive_control(
            module, prefixes, rails, present, present_starts, grid
        )
    )

    print("LRC14 exact D-compatible three-edge nilpotence probe")
    print(
        "scope=THM-2616 arrival-six raw nonnegative carrier; all "
        "sharp/predecessor-carry restrictions downstream"
    )
    print("object=E0(x) AND E1(Dx) AND E2(D^2x), D(x)={13x}")
    print("typing=owner_i=shallow_(i+1), h_i=j_(i+1); sources may drift independently")
    print(
        f"rails={len(rails)} rail_weighted_pieces={rail_piece_count} "
        f"dilation_clock_covariance_checks={dilation_clock_covariance_checks} "
        "shallow_next_phase_checks="
        f"{shallow_phase_identification_checks}"
    )
    print("two_arrival_phase_interval=[84/169,85/169) mass=1/169")
    print(
        "phase_clock_pair_support="
        + str(tuple(sorted(
            (pair, mass // (grid // 338))
            for pair, mass in phase_pair_masses.items()
        )))
        + " units_of_1/338"
    )
    print("nonzero_clock_edge_phase_support=empty (diagonal repeat clocks only)")
    print(
        f"clock_quadruples={clock_quadruples} source_triples=1728 "
        f"clock_source_candidates={clock_source_candidates} survivors=0"
    )
    print(
        "rail_triple_multiplicity_hist="
        + str(tuple(sorted(rail_multiplicity_histogram.items())))
    )
    print(
        f"rail_labelled_three_event_candidates={rail_labelled_triples} "
        "nominal_predecessor_carry_atom_triples="
        f"{nominal_labelled_atom_triples}"
    )
    print(
        f"positive_two_edge_numerator_certificate={witness} "
        f"base_pieces={witness_base_pieces} delayed_joint_pieces={witness_q_pieces}"
    )
    print(
        "verdict=PASS: D-compatible support transfer has exact nilpotence "
        "index three on the inherited nonzero-clock-edge carrier"
    )
    print(
        "mechanism=arrival-six return preserves the lower/upper half, forcing "
        "ell0=ell1 against the first stored edge (shallow ell0, owner ell1)"
    )
    print(
        "semantics=all longer D-chains vanish already at the arrival/clock "
        "skeleton; no Bockstein, endpoint, Perron, scalar-cover, or LRC claim"
    )


if __name__ == "__main__":
    main()
