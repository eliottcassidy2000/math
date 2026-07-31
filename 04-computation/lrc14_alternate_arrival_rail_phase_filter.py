#!/usr/bin/env python3
"""Exact LRC(14) rail filter for alternate arrival-clock evasions.

The general p=13,q=7 return classification says that a fixed arrival digit
can support a nonzero shallow-clock return exactly for

    m in {0,2,4,8,10,12}.

This script asks which of those evasion *arrival labels* also occur in the
upstream rail tensor used by THM-2616.  It rebuilds THM-2584's exact
``K[s][m][t][ell]`` tensor on the canonical typed row and restricts to the
same inherited deep rail labels ``t in {0,12}``.  Every one of the

    13 * 2 * 12 * 7 = 2,184

arrival/deep/source/clock cells is independently reconstructed from the
nonnegative step profile.  The arrival-six bank is then compared piecewise
with ``cross.build_carrier_data()``, proving that this is the actual parent
of THM-2616's 162 rails rather than a lookalike tensor.

The exact arrival support is only ``{0,6,12}``:

* m=0: 81 positive cells, all on t=0;
* m=6: 162 positive cells, 81 on each of t=0,12;
* m=12: 81 positive cells, all on t=12.

All other arrival digits and all other deep labels are zero.  Intersecting
with the phase-evasion set leaves precisely ``{0,12}``.  Each endpoint sheet
has off-diagonal return mass ``1/2366``: clock pair ``0->1`` for m=0 and
``0->6`` for m=12.  Thus the phase-nearest digits 4,8 and varying arrivals
5->6,7->6 do not exist on this inherited rail tensor; the reflected endpoint
sheets are the only already-present alternate-arrival candidates.

The two endpoint banks form a sharp base-cell cospan: each covers 81 of the
84 ``(s,ell)`` cells, their overlap has size 78, and their union is all 84.
Canonical half-open reflection, equal to literal reflection away from null
endpoints,

    (m,t,s,ell,[a,b),w)
      -> (12-m,12-t,13-s,-ell,[T-b,T-a),w)

identifies all 81 weighted cell profiles (537 constituent intervals on each
leg).  On the supported arrival
alphabet ``{0,6,12}``, every ordered phase transition has positive
off-diagonal clock mass except ``6->6``.  Endpoint loops have off-diagonal
mass ``1/2366``; all six varying-arrival transitions have mass ``1/169``.

This is a sharp support filter, not an endpoint-carrier construction.  A
positive rail cell times a positive phase interval does not yet prove their
same-x intersection, delayed-word positivity, a Bockstein unit, endpoint
section, coherent three-event chain, scalar exclusion, or LRC(14).
"""

from collections import Counter
from fractions import Fraction

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old
import lrc_central_arrival_clock_return_classification as phase


P = 13
Q7 = 7
DEEP_RAIL_LABELS = (0, 12)

EXPECTED_SUPPORT_COUNTS = {
    0: (81, 0),
    1: (0, 0),
    2: (0, 0),
    3: (0, 0),
    4: (0, 0),
    5: (0, 0),
    6: (81, 81),
    7: (0, 0),
    8: (0, 0),
    9: (0, 0),
    10: (0, 0),
    11: (0, 0),
    12: (0, 81),
}

COMMON_DENOMINATOR = 6_939_029_398_456_584_394_954_868_880

EXPECTED_NUMERATORS = {
    (0, 0): 42_233_787_706_033_899_602_125_200,
    (0, 12): 0,
    (6, 0): 42_220_444_026_338_661_554_256_240,
    (6, 12): 42_220_444_026_338_661_554_256_240,
    (12, 0): 0,
    (12, 12): 42_233_787_706_033_899_602_125_200,
}

EXPECTED_MISSING = {
    (0, 0): ((7, 4), (7, 5), (7, 6)),
    (6, 0): ((7, 1), (7, 2), (7, 3)),
    (6, 12): ((6, 4), (6, 5), (6, 6)),
    (12, 12): ((6, 1), (6, 2), (6, 3)),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    module = old.load_present_module()
    require(module.W == old.base.W,
            "canonical typed row changed before the arrival audit")
    require(old.P == P and old.Q7 == Q7,
            "inherited 13/7 scale constants changed")

    # This is exactly the THM-2584 tensor construction imported by THM-2616.
    e4 = old.base.build_set(old.base.PAT_E, old.base.ZELL)
    qb = old.base.build_set(old.host.PAT_QB, old.base.ZELL)
    ust, uv, vst, vv = old.rail.packet_profiles(e4, qb)
    _, _, k_tensor = old.rail.exact_tensor(e4, qb)
    owner = old.base.clock_cells(P, Q7, old.T, P * P)
    deep = old.rail.deep_cells()

    reconstructed = {}
    cells_checked = 0
    support_counts = {}
    numerator_totals = {}
    missing = {}
    for arrival_digit in range(P):
        arrival = [(
            arrival_digit * old.T // P,
            (arrival_digit + 1) * old.T // P,
        )]
        for deep_label in DEEP_RAIL_LABELS:
            positives = []
            numerator_total = 0
            for s in range(1, P):
                rvst, rvv = old.base.rotate_profile(
                    vst, vv, s * old.T // P, old.T
                )
                profile_starts, profile_values, _, _ = old.base.product_cum(
                    ust, uv, rvst, rvv, old.T
                )
                for ell in range(Q7):
                    cell = old.intersect_sorted(
                        old.intersect_sorted(owner[ell], deep[deep_label]),
                        arrival,
                    )
                    pieces = old.profile_on_intervals(
                        cell, profile_starts, profile_values
                    )
                    numerator = P * sum(
                        weight * (right - left)
                        for left, right, weight in pieces
                    )
                    require(
                        numerator
                        == k_tensor[s][arrival_digit][deep_label][ell],
                        "step-profile reconstruction disagrees with K tensor",
                    )
                    cells_checked += 1
                    numerator_total += numerator
                    if numerator:
                        positives.append((s, ell))
                        reconstructed[
                            arrival_digit, s, ell, deep_label
                        ] = pieces
            support_counts[arrival_digit, deep_label] = len(positives)
            numerator_totals[arrival_digit, deep_label] = numerator_total
            missing[arrival_digit, deep_label] = tuple(
                (s, ell)
                for s in range(1, P)
                for ell in range(Q7)
                if (s, ell) not in positives
            )

    require(cells_checked == P * len(DEEP_RAIL_LABELS) * 12 * Q7
            == 2_184,
            "arrival/deep/source/clock universe changed")
    require(
        {
            m: tuple(support_counts[m, t] for t in DEEP_RAIL_LABELS)
            for m in range(P)
        } == EXPECTED_SUPPORT_COUNTS,
        "arrival/deep support census changed",
    )
    for key, value in EXPECTED_NUMERATORS.items():
        require(numerator_totals[key] == value,
                f"arrival/deep total numerator changed at {key}")
    for key, value in EXPECTED_MISSING.items():
        require(missing[key] == value,
                f"exceptional source/clock holes changed at {key}")

    endpoint_support_0 = {
        (s, ell) for s in range(1, P) for ell in range(Q7)
        if (0, s, ell, 0) in reconstructed
    }
    endpoint_support_12 = {
        (s, ell) for s in range(1, P) for ell in range(Q7)
        if (12, s, ell, 12) in reconstructed
    }
    endpoint_overlap = endpoint_support_0 & endpoint_support_12
    endpoint_union = endpoint_support_0 | endpoint_support_12
    endpoint_private_0 = tuple(sorted(
        endpoint_support_0 - endpoint_support_12
    ))
    endpoint_private_12 = tuple(sorted(
        endpoint_support_12 - endpoint_support_0
    ))
    require(
        (len(endpoint_support_0), len(endpoint_support_12),
         len(endpoint_overlap), len(endpoint_union)) == (81, 81, 78, 84),
        "endpoint base-cell cospan census changed",
    )
    require(endpoint_private_0 == ((6, 1), (6, 2), (6, 3))
            and endpoint_private_12 == ((7, 4), (7, 5), (7, 6)),
            "endpoint cospan private cells changed")

    reflection_profile_checks = 0
    reflection_interval_checks = 0
    for s, ell in sorted(endpoint_support_0):
        pieces = reconstructed[0, s, ell, 0]
        reflected = sorted(
            (old.T - right, old.T - left, weight)
            for left, right, weight in pieces
        )
        target_key = (12, P - s, (-ell) % Q7, 12)
        require(target_key in reconstructed,
                "endpoint reflection missed its target rail cell")
        require(reflected == reconstructed[target_key],
                "endpoint weighted pieces violate exact reflection")
        require(
            k_tensor[s][0][0][ell]
            == k_tensor[P - s][12][12][(-ell) % Q7],
            "endpoint tensor masses violate exact reflection",
        )
        reflection_profile_checks += 1
        reflection_interval_checks += len(pieces)
    require(reflection_profile_checks == 81,
            "endpoint piecewise reflection census changed")
    require(reflection_interval_checks == 537,
            "endpoint constituent-interval census changed")

    # No unlisted deep root carries any mass anywhere in the full tensor.
    other_deep_positive_cells = sum(
        k_tensor[s][m][t][ell] > 0
        for m in range(P)
        for t in range(P)
        if t not in DEEP_RAIL_LABELS
        for s in range(1, P)
        for ell in range(Q7)
    )
    require(other_deep_positive_cells == 0,
            "an unlisted deep rail label acquired positive mass")

    # Piecewise identity with THM-2616's actual arrival-six rail bank.
    (_, _, _, _, inherited_rails, _, _) = cross.build_carrier_data()
    inherited_bank = {
        (s, ell, deep_label): pieces
        for s, ell, deep_label, pieces in inherited_rails
    }
    rebuilt_central_bank = {
        (s, ell, deep_label): pieces
        for (m, s, ell, deep_label), pieces in reconstructed.items()
        if m == 6
    }
    require(len(inherited_bank) == len(rebuilt_central_bank) == 162,
            "arrival-six rail bank census changed")
    require(inherited_bank == rebuilt_central_bank,
            "audited tensor is not piecewise identical to THM-2616 rails")

    rail_arrivals = tuple(
        m for m in range(P)
        if any(support_counts[m, t] for t in DEEP_RAIL_LABELS)
    )
    require(rail_arrivals == (0, 6, 12),
            "inherited rail arrival support changed")
    phase_evasions = tuple(
        m for m in range(P)
        if phase.off_diagonal_mass(
            phase.forward_cross_arrival_support(Q7, m, m)
        )
    )
    require(phase_evasions == (0, 2, 4, 8, 10, 12),
            "13/7 phase-evasion set changed")
    candidate_evasion_arrivals = tuple(
        sorted(set(rail_arrivals) & set(phase_evasions))
    )
    require(candidate_evasion_arrivals == (0, 12),
            "rail/phase filter no longer leaves only endpoint arrivals")

    phase_unit = Fraction(1, 2 * Q7 * P * P)
    endpoint_phase_support = {
        m: phase.normalized_support(
            phase.forward_cross_arrival_support(Q7, m, m),
            phase_unit,
        )
        for m in candidate_evasion_arrivals
    }
    require(endpoint_phase_support == {
        0: (((0, 0), 13), ((0, 1), 1)),
        12: (((0, 0), 13), ((0, 6), 1)),
    }, "endpoint arrival phase support changed")
    require(all(
        phase.off_diagonal_mass(
            phase.forward_cross_arrival_support(Q7, m, m)
        ) == phase_unit
        for m in candidate_evasion_arrivals
    ), "endpoint off-diagonal return mass changed")

    # Full 3x3 phase transition atlas on the arrival alphabet that actually
    # exists in the inherited rail tensor.  This remains a phase marginal:
    # it is not multiplied into independently positive rail cells.
    phase_transition_atlas = {
        (first, second): phase.normalized_support(
            phase.forward_cross_arrival_support(Q7, first, second),
            phase_unit,
        )
        for first in rail_arrivals
        for second in rail_arrivals
    }
    expected_phase_transition_atlas = {
        (0, 0): (((0, 0), 13), ((0, 1), 1)),
        (0, 6): (((0, 3), 7), ((0, 4), 7)),
        (0, 12): (((0, 6), 1), ((1, 0), 13)),
        (6, 0): (((3, 0), 13), ((3, 1), 1)),
        (6, 6): (((3, 3), 7), ((4, 4), 7)),
        (6, 12): (((4, 0), 13), ((4, 6), 1)),
        (12, 0): (((0, 1), 1), ((6, 0), 13)),
        (12, 6): (((0, 3), 7), ((0, 4), 7)),
        (12, 12): (((0, 0), 13), ((0, 6), 1)),
    }
    require(phase_transition_atlas == expected_phase_transition_atlas,
            "supported-arrival 3x3 phase transition atlas changed")
    off_diagonal_mass_matrix = tuple(
        tuple(int(
            phase.off_diagonal_mass(
                phase.forward_cross_arrival_support(Q7, first, second)
            ) / phase_unit
        ) for second in rail_arrivals)
        for first in rail_arrivals
    )
    require(off_diagonal_mass_matrix == (
        (1, 14, 14),
        (14, 0, 14),
        (14, 14, 1),
    ), "supported-arrival off-diagonal mass matrix changed")

    print("LRC14 exact alternate-arrival rail/phase filter")
    print("tensor=THM-2584 K[s][m][t][ell] on the THM-2616 canonical typed row")
    print("deep_rail_labels=(0,12); all 13 arrivals, 12 sources, 7 clocks checked")
    print(
        f"arrival_deep_cells_checked={cells_checked} "
        f"other_deep_positive_cells={other_deep_positive_cells}"
    )
    print(
        "support_counts_by_arrival_(t0,t12)="
        + str(tuple(
            (m, tuple(support_counts[m, t] for t in DEEP_RAIL_LABELS))
            for m in range(P)
        ))
    )
    print(
        f"positive_numerator_totals_(denominator={COMMON_DENOMINATOR})="
        + str(tuple(
            (key, numerator_totals[key])
            for key in ((0, 0), (6, 0), (6, 12), (12, 12))
        ))
    )
    print(
        "exceptional_missing_cells="
        + str(tuple(sorted(EXPECTED_MISSING.items())))
    )
    print(
        "endpoint_projected_label_cospan="
        "(left=81,right=81,overlap=78,union=84,"
        f"left_private={endpoint_private_0},"
        f"right_private={endpoint_private_12})"
    )
    print(
        "endpoint_reflection=(m,t,s,ell)->"
        f"(12-m,12-t,13-s,-ell); a.e._half_open_profile_checks="
        f"{reflection_profile_checks}/81; constituent_intervals_each_leg="
        f"{reflection_interval_checks}"
    )
    print("arrival6_piecewise_identity_with_THM2616=162/162")
    print(f"rail_arrival_support={rail_arrivals}")
    print(f"phase_same_arrival_evasions={phase_evasions}")
    print(
        "rail_label_intersect_phase_evasion_labels="
        f"{candidate_evasion_arrivals}"
    )
    print(
        "endpoint_phase_support_units_1_over_2366="
        + str(tuple(sorted(endpoint_phase_support.items())))
    )
    print("endpoint_offdiagonal_mass=m0:m12=1/2366:1/2366")
    print(
        "supported_arrival_phase_atlas_units_1_over_2366="
        + str(tuple(sorted(phase_transition_atlas.items())))
    )
    print(
        "supported_arrival_offdiagonal_mass_matrix_rows_0_6_12="
        + str(off_diagonal_mass_matrix)
    )
    print(
        "verdict=PASS: only the reflected endpoint arrival sheets 0 and 12 "
        "are already-present rail-level candidates to evade the central trap"
    )
    print(
        "semantics=rail positivity and phase support are separate marginals; "
        "no same-x endpoint chain, delayed unit, endpoint section, scalar "
        "exclusion, or LRC(14) follows"
    )
    print(
        "next_test=intersect each phase transition with the exact pulled-back "
        "source/owner rail pieces; do not lift the 3x3 phase graph formally"
    )


if __name__ == "__main__":
    main()
