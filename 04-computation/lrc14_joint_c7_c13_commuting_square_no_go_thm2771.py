#!/usr/bin/env python3
"""Exact THM-2771 fixed-clock chart x target-label square audit.

The target-label operator acts on the labelled raw-section family before the
chart pullback.  The chart pullback is the fixed tau translation followed by
the common rail cut.  Their two object-level composites are compared as
literal interval unions.  The nonconstant right-cofiber target profile is
then separated from this (zero) commutator.
"""

from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_fully_marked_root_zero_target_profile_thm2749 as marked
import lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751 as wing


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def signed_profile(channels):
    """Canonical signed step profile from (interval-union, coefficient)."""
    events = {}
    for pieces, coefficient in channels:
        for left, right in pieces:
            events[left] = events.get(left, 0) + coefficient
            events[right] = events.get(right, 0) - coefficient
    result = []
    value = 0
    previous = None
    for point in sorted(events):
        if previous is not None and previous < point and value:
            if result and result[-1][1] == previous and result[-1][2] == value:
                result[-1] = (result[-1][0], point, value)
            else:
                result.append((previous, point, value))
        value += events[point]
        previous = point
    require(value == 0, "signed interval profile did not close")
    return tuple(result)


def l1_mass(profile):
    return sum((right - left) * abs(value) for left, right, value in profile)


def rank_q(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = 1 / work[rank][column]
        work[rank] = [inverse * x for x in work[rank]]
        for i in range(len(work)):
            if i == rank or not work[i][column]:
                continue
            value = work[i][column]
            work[i] = [a - value * b for a, b in zip(work[i], work[rank])]
        rank += 1
    return rank


def main():
    module, _prefixes, _a, _b, rails, present, _starts = (
        marked.lift.m.core.build_carrier_data()
    )
    source = marked.two.exclusive_source(module, 3)
    clock_one = tuple(module.make_comb(module.C1, 182, 13, 39))
    _source_weight, _target_weight, rail_common = marked.rail_data(
        rails, marked.RAIL
    )
    static = tuple(marked.two.intersect_sorted(source, clock_one))
    static = tuple(module.subtract_comb(static, module.W[1], 182, -13, 13))
    static = tuple(module.subtract_comb(static, module.C2, 182, -13, 13))

    raw = []
    A = []
    B = []
    for t in range(P):
        row = tuple(module.subtract_comb(
            static, module.W[2], 182, -14 * t - 13, -14 * t + 13
        ))
        row = tuple(module.subtract_comb(
            row, module.C3, 182, 14 * t - 13, 14 * t + 13
        ))
        raw.append(row)
        A.append(marked.intersect(rail_common, row))
        B.append(marked.intersect(
            rail_common, marked.shift_union(row, marked.SHIFT)
        ))
    raw, A, B = tuple(raw), tuple(A), tuple(B)

    # A labelled target shift sigma_delta sends the state (chart,t) to
    # (chart,t+delta); tau sends (source,t) to (target,t).  Both composites
    # therefore realize the same literal B_(t+delta) interval object.
    commutator_failures = []
    for delta in range(P):
        for t in range(P):
            chart_after_target = marked.intersect(
                rail_common,
                marked.shift_union(raw[(t + delta) % P], marked.SHIFT),
            )
            target_after_chart = B[(t + delta) % P]
            if chart_after_target != target_after_chart:
                commutator_failures.append((delta, t))
    require(not commutator_failures,
            "chart pullback and lawful target co-shift stopped commuting")

    # The fixed-clock normalized right cofiber from THM-2751, after removing
    # its common scalar G.  It is nonconstant and every primitive C13
    # character survives, but its target differences are boundary data.
    r = tuple(2 if 3 <= t <= 11 else 121 if t == 12 else 0
              for t in range(P))
    difference_matrix = tuple(
        tuple(r[(t + delta) % P] - r[t] for t in range(P))
        for delta in range(P)
    )
    require(rank_q(difference_matrix) == 12,
            "right-cofiber target boundary lost its full augmentation rank")

    # Record the literal alternating four-corner boundary after the two
    # relative-present cuts and source seam.  It may be nonzero even though
    # the operation commutator is zero: it is d_target d_chart(A), not
    # [d_target,d_chart](A).
    cut_a = tuple(tuple(
        wing.cut_carrier(module, present, A[t], ell) for ell in range(7)
    ) for t in range(P))
    cut_b = tuple(tuple(
        wing.cut_carrier(module, present, B[t], ell) for ell in range(7)
    ) for t in range(P))
    boundary_census = []
    for delta in range(P):
        nonzero_cells = 0
        total_l1 = 0
        for t in range(P):
            td = (t + delta) % P
            for ell in range(7):
                a_t = cut_a[t][ell]
                a_td = cut_a[td][ell]
                b_t = cut_b[t][ell]
                b_td = cut_b[td][ell]
                boundary = signed_profile((
                    (b_td, 1), (b_t, -1), (a_td, -1), (a_t, 1),
                ))
                if boundary:
                    nonzero_cells += 1
                    total_l1 += l1_mass(boundary)
        boundary_census.append((delta, nonzero_cells, total_l1))
    expected_boundary_census = (
        (0, 0, 0),
        (1, 21, 89595253440),
        (2, 35, 90705938400),
        (3, 49, 91816623360),
        (4, 56, 92186851680),
        (5, 56, 92186851680),
        (6, 56, 92186851680),
        (7, 56, 92186851680),
        (8, 56, 92186851680),
        (9, 56, 92186851680),
        (10, 49, 91816623360),
        (11, 35, 90705938400),
        (12, 21, 89595253440),
    )
    require(tuple(boundary_census) == expected_boundary_census,
            "exact four-corner boundary census changed")

    print("FIXED-CLOCK CHART x TARGET-LABEL FOUR-CORNER AUDIT")
    print("status=THM-2771 FINITE-EXACT candidate; object commutator versus boundary profile")
    print("object_commutator_failures=() over 13*13 labelled squares")
    print(f"primitive_right_cofiber={r}")
    print(f"right_cofiber_difference_rank_Q={rank_q(difference_matrix)}")
    print(f"four_corner_boundary_census_(delta,nonzero_t_ell_cells,L1_grid_mass)={tuple(boundary_census)}")
    print("verdict=target and chart actions commute exactly; nonzero Delta_target(r) and the alternating four-corner support are mixed coboundaries, not a transition commutator/holonomy 2-cocycle")
    print("THM2591_SCOPE: does not pay the mixed-square invoice without an independently noncommuting semantic/common-ancestry transition or reference")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
