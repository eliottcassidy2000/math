#!/usr/bin/env python3
"""Exact companion for THM-2829.

This companion starts after the fixed-cell THM-2820 covariance
audit.  It asks whether the q=11 source-chart commutator boundary can repair
the fixed t=4 q2 failure by choosing another lawful cell in THM-2806's full
9*9*7 bank.  It separates:

* native fixed-factor and weighted carrier support;
* endpoint/address/delayed-semantic survival; and
* the outer THM-2584 Boolean ancestry sheet.

It proves only the finite statements asserted in THM-2829 and claims no row
decrement or LRC(14) consequence.
"""

from __future__ import annotations

from hashlib import sha256
import importlib.util
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    COMP / "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
    COMP / "lrc14_semantic_arm_right_wing_central_digit_thm2782.py":
        "7fbc6bb1ec303ded98eaad6e5d8205eb3d247258ada32b6f9904fc439ebb11fb",
    COMP / "lrc14_replica_dichotomy_typed_row_opus_20260727.py":
        "6ba64a68a9fd008d2e06949b1f1cf75012f1f4e734f75f55ce0af58ae20ad7b9",
}
for path, expected in PINNED.items():
    actual = sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == expected, f"pinned dependency changed: {path.name}")


import lrc14_literal_fixed_sheet_allocation_thm2806 as allocation
import lrc14_semantic_arm_right_wing_central_digit_thm2782 as arm


# One inherited loader pins raw LF bytes.  Load the already LF-verified
# dependency directly on Windows without weakening the check above.
def crlf_safe_present_loader():
    path = COMP / "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
    spec = importlib.util.spec_from_file_location(
        "q11_present_base", path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


allocation.physical.target.load_present_module = crlf_safe_present_loader

P = 13
Q = 11


def circular_shift_interval(interval, shift, period):
    left = (interval[0] + shift) % period
    right = (interval[1] + shift) % period
    require(left < right, "test interval crossed the circle seam")
    return left, right


def safe_comb_contains(interval, module, speed, denominator, lo, hi):
    """Test containment in subtract_comb without materializing its union."""
    require(
        module.T % (denominator * speed) == 0,
        "comb grid ceased to resolve",
    )
    unit = module.T // (denominator * speed)
    step = denominator * unit
    length = (hi - lo) * unit
    base = (lo % denominator) * unit
    left, right = interval
    k0 = (left - base) // step
    for k in range(k0 - 1, k0 + 3):
        forbidden_left = base + k * step
        forbidden_right = forbidden_left + length
        if max(left, forbidden_left) < min(right, forbidden_right):
            return False
    return True


def weighted_hit(interval, pieces):
    hits = tuple(
        (left, right, weight)
        for left, right, weight in pieces
        if left <= interval[0] and interval[1] <= right
    )
    require(len(hits) <= 1, "test atom met two weighted pieces")
    return hits


def source_center(address):
    return arm.arm_source_center(address, 0) * arm.T


def target_center(address):
    return arm.arm_target_center(address, 0) * arm.T


def cylinder(center, half_width):
    return center - half_width, center + half_width


def interval_contained(interval, container):
    return container[0] <= interval[0] and interval[1] <= container[1]


def main():
    (
        _module, full_module, details, e3, clocks, q_pairs
    ) = allocation.build_geometry()
    period = full_module.T
    allocation_unit = period // P
    I = allocation.ATOM_INTERVAL
    J0 = tuple(endpoint + allocation.physical.SHIFT for endpoint in I)
    J11 = circular_shift_interval(
        J0, Q * allocation_unit, period
    )

    def signature(interval, s, t, clock):
        return (
            allocation.contained(interval, e3),
            allocation.contained(interval, clocks[clock]),
            safe_comb_contains(
                interval, full_module, full_module.W[1], 182,
                -14 * s - 13, -14 * s + 13,
            ),
            safe_comb_contains(
                interval, full_module, full_module.W[2], 182,
                -14 * t - 13, -14 * t + 13,
            ),
            safe_comb_contains(
                interval, full_module, full_module.C2, 182,
                14 * s - 13, 14 * s + 13,
            ),
            safe_comb_contains(
                interval, full_module, full_module.C3, 182,
                14 * t - 13, 14 * t + 13,
            ),
        )

    full_cells = []
    for s in allocation.COMMON_S:
        for t in allocation.COMMON_T:
            for clock in range(7):
                if all(signature(I, s, t, clock)):
                    if (
                        all(signature(J0, s, t, clock))
                        and all(signature(J11, s, t, clock))
                    ):
                        full_cells.append((s, t, clock))
    expected_cells = tuple(
        (s, t, 1)
        for s in allocation.COMMON_S
        for t in range(5, 12)
    )
    require(
        tuple(full_cells) == expected_cells and len(full_cells) == 63,
        "q=11 full-factor reselection bank changed",
    )

    # Spot-check the arithmetic containment shortcut against the canonical
    # interval builders on the cheapest cell.
    witness_cell = (0, 5, 1)
    witness_factors = allocation.section_factors(
        full_module, e3, clocks, *witness_cell
    )
    for interval in (I, J0, J11):
        require(
            tuple(
                allocation.contained(interval, factor)
                for factor in witness_factors
            ) == signature(interval, *witness_cell) == (True,) * 6,
            "direct and materialized factor signatures disagree",
        )

    # The base clock-one carrier has equal literal weights at I and J0.
    # Translating the target carrier and endpoint together by q*T/13, then
    # pulling back by SHIFT+q*T/13, returns exactly the same weighted atom.
    source_base, target_base = details[1]
    source_hit = weighted_hit(I, source_base)
    target_hit = weighted_hit(J0, target_base)
    target_h11 = allocation.physical.overlap.shift_weighted(
        target_base, Q * allocation_unit
    )
    target_h11_hit = weighted_hit(J11, target_h11)
    require(
        source_hit == ((*I, allocation.ATOM_WEIGHT),)
        and target_hit == ((*J0, allocation.ATOM_WEIGHT),)
        and target_h11_hit == ((*J11, allocation.ATOM_WEIGHT),),
        "q=11 weighted carrier atom changed",
    )
    pulled_h11 = allocation.physical.overlap.shift_weighted(
        target_h11,
        -(allocation.physical.SHIFT + Q * allocation_unit),
    )
    pulled_hit = weighted_hit(I, pulled_h11)
    require(
        pulled_hit == source_hit,
        "q=11 source/target pulled weights ceased to agree",
    )

    # The physical address/root cylinder is still of low source residue 8
    # (equivalently target residue 7), but its high address has moved.
    n_source = (
        6716 + 9 * P**5
    ) % P**6
    n_target = (n_source - 1) % P**6
    half_width = arm.relative.Q_RADIUS * period
    require(
        (n_source, n_target) == (3348353, 3348352)
        and (n_source % P, n_target % P) == (8, 7)
        and interval_contained(
            cylinder(source_center(n_source), half_width), J11
        )
        and interval_contained(
            cylinder(target_center(n_target), half_width), J11
        ),
        "q=11 address/root cylinder changed",
    )

    # Delayed semantic/root data survive the physical q=11 translate.
    source_semantic = allocation.physical.relative.private.delayed_carry_pair(
        ((*I, 1),), q_pairs[1], {}
    )
    target_semantic = allocation.physical.relative.private.delayed_carry_pair(
        ((*J11, 1),), q_pairs[1], {}
    )
    require(
        source_semantic[12] == (0, allocation.SEMANTIC_VALUE)
        and target_semantic[6] == (0, allocation.SEMANTIC_VALUE),
        "q=11 delayed semantic coefficient changed",
    )

    # Both right endpoint origins remain live with exactly the base scalar.
    endpoint_base = allocation.endpoint_base
    phase_exponent = (
        endpoint_base.Y0
        * endpoint_base.endpoint.RDIL
        * allocation_unit
    )
    require(
        phase_exponent % endpoint_base.endpoint.NN == 0,
        "allocation-unit endpoint phase ceased to be trivial",
    )
    endpoint_value = (
        231164267889491750,
        630230755085920022,
    )
    endpoint_rows = {}
    target_atom = ((*J0, 1),)
    for address in (
        allocation.RIGHT_ORIGIN,
        allocation.add(
            allocation.RIGHT_ORIGIN, allocation.TARGET_STEP
        ),
    ):
        ell = endpoint_base.REPS[address]
        present = tuple(
            endpoint_base.endpoint.build_set(endpoint_base.PAT_E3, ell)
        )
        starts = tuple(left for left, _right in present)
        rows = []
        for q in (0, Q):
            shifted = allocation.physical.overlap.shift_weighted(
                target_atom, q * allocation_unit
            )
            restricted = allocation.indexed_weighted_intersection(
                shifted, present, starts
            )
            values = tuple(
                allocation.endpoint_sum(
                    restricted, -endpoint_base.Y0, embedding
                )
                for embedding in endpoint_base.endpoint.MODS
            )
            mass = sum(
                (right - left) * weight
                for left, right, weight in restricted
            )
            rows.append((q, len(restricted), mass, values))
        endpoint_rows[address] = tuple(rows)
    expected_endpoint_row = (
        (0, 1, 26444880, endpoint_value),
        (Q, 1, 26444880, endpoint_value),
    )
    require(
        all(
            rows == expected_endpoint_row
            for rows in endpoint_rows.values()
        ),
        "q=11 two-right-origin endpoint edge changed",
    )

    # The outer THM-2584 U_Q ancestry gate is the first exact obstruction.
    # Enumerate all 13^5 ancestry labels for every allocation residue.  The
    # q=11 target has no active U_Q label at all, so changing s,t,clock or
    # the delayed carry/root label cannot repair the missing common sheet.
    old = allocation.physical.relative.lift.m.core.old
    q_intervals = tuple(
        old.base.build_set(old.host.PAT_QB, old.base.ZELL)
    )
    q_starts = tuple(left for left, _right in q_intervals)
    depth = old.rail.DEPTH
    midpoint = sum(I) // 2
    source_flags = tuple(
        allocation.interval_hit(
            midpoint + ancestry * period,
            depth, q_intervals, q_starts,
        ) is not None
        for ancestry in range(depth)
    )
    source_count = sum(source_flags)
    target_counts = []
    common_counts = []
    first_common = []
    target_flags_by_q = []
    for q in range(P):
        offset = allocation.physical.SHIFT + q * allocation_unit
        target_flags = tuple(
            allocation.interval_hit(
                midpoint + offset + ancestry * period,
                depth, q_intervals, q_starts,
            ) is not None
            for ancestry in range(depth)
        )
        target_flags_by_q.append(target_flags)
        target_counts.append(sum(target_flags))
        common = tuple(
            source and target
            for source, target in zip(source_flags, target_flags)
        )
        common_counts.append(sum(common))
        first_common.append(next(
            (index for index, value in enumerate(common) if value),
            None,
        ))
    expected_target_counts = (
        66099, 0, 0, 0, 0, 0, 0, 65652, 0, 0, 0, 0, 0,
    )
    expected_common_counts = (
        66099, 0, 0, 0, 0, 0, 0, 65612, 0, 0, 0, 0, 0,
    )
    require(
        depth == P**5
        and source_count == 66099
        and tuple(target_counts) == expected_target_counts
        and tuple(common_counts) == expected_common_counts
        and first_common[0] == first_common[7] == 59162
        and all(
            first_common[q] is None
            for q in range(P) if q not in (0, 7)
        ),
        "THM-2584 U_Q ancestry spectrum changed",
    )

    # A residue arrow q -> q+h has a carry on the natural-extension sheet.
    # In this numerator convention 13*a+q -> 13*a+q+h, so the sheet label
    # changes by floor((q+h)/13).  The residue-only common count therefore
    # ceases to be the lifted common count when q+h wraps.
    arrow_pairs = (
        (0, 0), (0, 7), (3, 10), (3, 4), (11, 2), (11, 9),
    )

    def lifted_common_weights(base_point):
        base_source = tuple(
            allocation.interval_hit(
                base_point + ancestry * period,
                depth, q_intervals, q_starts,
            ) is not None
            for ancestry in range(depth)
        )
        base_targets = {}
        for residue in (0, 7):
            offset = allocation.physical.SHIFT + residue * allocation_unit
            base_targets[residue] = tuple(
                allocation.interval_hit(
                    base_point + offset + ancestry * period,
                    depth, q_intervals, q_starts,
                ) is not None
                for ancestry in range(depth)
            )
        return tuple(
            (
                q,
                h,
                (q + h) // P,
                sum(
                    base_source[ancestry]
                    and base_targets[(q + h) % P][
                        (ancestry + (q + h) // P) % depth
                    ]
                    for ancestry in range(depth)
                ),
            )
            for q, h in arrow_pairs
        )

    natural_lifted_common = tuple(
        (
            q,
            h,
            (q + h) // P,
            sum(
                source_flags[ancestry]
                and target_flags_by_q[(q + h) % P][
                    (ancestry + (q + h) // P) % depth
                ]
                for ancestry in range(depth)
            ),
        )
        for q, h in arrow_pairs
    )
    expected_natural_lifted_common = (
        (0, 0, 0, 66099),
        (0, 7, 0, 65612),
        (3, 10, 1, 65579),
        (3, 4, 0, 65612),
        (11, 2, 1, 65579),
        (11, 9, 1, 65098),
    )
    require(
        natural_lifted_common == expected_natural_lifted_common
        and lifted_common_weights(I[0])
        == expected_natural_lifted_common
        and lifted_common_weights(I[1] - 1)
        == expected_natural_lifted_common,
        "borrow-aware whole-cylinder ancestry weights changed",
    )
    same_section_common = (
        66099, 65612, 66099, 65612, 66099, 65612,
    )
    carry_defect = tuple(
        lifted[-1] - same
        for lifted, same in zip(
            natural_lifted_common, same_section_common
        )
    )
    require(
        carry_defect == (0, 0, -520, 0, -520, -514)
        and tuple(value % P for value in carry_defect)
        == (0, 0, 0, 0, 0, 6)
        and gcd(abs(carry_defect[-1]), 91) == 1,
        "borrow-response modular selector changed",
    )
    q7_flags = target_flags_by_q[7]
    reverse_q11_h9 = sum(
        source_flags[ancestry]
        and q7_flags[(ancestry - 1) % depth]
        for ancestry in range(depth)
    )
    require(
        reverse_q11_h9 == 65619,
        "reverse-orientation q11/h9 control changed",
    )
    e3_support = tuple(
        q for q in range(P)
        if allocation.contained(
            circular_shift_interval(
                J0, q * allocation_unit, period
            ),
            e3,
        )
    )
    ancestry_support = tuple(
        q for q, count in enumerate(target_counts) if count
    )
    require(
        e3_support == (0, 3, 11)
        and ancestry_support == (0, 7)
        and not (
            (set(e3_support) - {0})
            & (set(ancestry_support) - {0})
        ),
        "E3/ancestry support obstruction changed",
    )

    # Exhaust the other two canonical outer words at q=11.  No fixed word
    # has a same-word common label.  The nearest positive cospan changes
    # the outer word from QB={b} on the source to QA={a} on the target.
    qa_pattern = old.base.PAT_QA
    qab_pattern = dict(qa_pattern)
    qab_pattern[8] = "in"
    require(
        (qa_pattern[7], qa_pattern[8])
        == ("in", "out")
        and (old.host.PAT_QB[7], old.host.PAT_QB[8])
        == ("out", "in")
        and (qab_pattern[7], qab_pattern[8])
        == ("in", "in"),
        "canonical outer-word signatures changed",
    )

    def word_flags(pattern, offset):
        pieces = tuple(old.base.build_set(pattern, old.base.ZELL))
        starts = tuple(left for left, _right in pieces)
        flags = tuple(
            allocation.interval_hit(
                midpoint + offset + ancestry * period,
                depth, pieces, starts,
            ) is not None
            for ancestry in range(depth)
        )
        return pieces, starts, flags

    qa_pieces, qa_starts, qa_source = word_flags(qa_pattern, 0)
    _qa_pieces_11, _qa_starts_11, qa_target_11 = word_flags(
        qa_pattern, allocation.physical.SHIFT + Q * allocation_unit
    )
    _qab_pieces, _qab_starts, qab_source = word_flags(qab_pattern, 0)
    _qab_pieces_11, _qab_starts_11, qab_target_11 = word_flags(
        qab_pattern, allocation.physical.SHIFT + Q * allocation_unit
    )
    cross_qb_to_qa_11 = tuple(
        ancestry
        for ancestry in range(depth)
        if source_flags[ancestry] and qa_target_11[ancestry]
    )
    require(
        sum(qa_source) == 0
        and sum(qa_target_11) == 10787
        and sum(qab_source) == 10786
        and sum(qab_target_11) == 0
        and target_counts[Q] == 0
        and all(
            not (left[ancestry] and right[ancestry])
            for left, right in (
                (qa_source, qa_target_11),
                (source_flags, target_flags_by_q[Q]),
                (qab_source, qab_target_11),
            )
            for ancestry in range(depth)
        )
        and len(cross_qb_to_qa_11) == 449
        and cross_qb_to_qa_11[0] == 59306
        and cross_qb_to_qa_11[-1] == 311961
        and {
            ancestry % (P * P)
            for ancestry in cross_qb_to_qa_11
        } == {156},
        "q11 canonical outer-word/cross-word census changed",
    )
    q11_target_offset = allocation.physical.SHIFT + Q * allocation_unit
    for ancestry in cross_qb_to_qa_11:
        source_indices = tuple(
            allocation.interval_hit(
                endpoint + ancestry * period,
                depth, q_intervals, q_starts,
            )
            for endpoint in (I[0], I[1] - 1)
        )
        target_indices = tuple(
            allocation.interval_hit(
                endpoint + q11_target_offset + ancestry * period,
                depth, qa_pieces, qa_starts,
            )
            for endpoint in (I[0], I[1] - 1)
        )
        require(
            source_indices[0] == source_indices[1] is not None
            and target_indices[0] == target_indices[1] is not None,
            "q11 QB-to-QA label ceased to be whole-cylinder",
        )

    print("Q11 SEMANTIC-RESELECTION AND ANCESTRY PROBE")
    print("status=THM-2829 VERIFIED-EXACT companion; no row decrement or LRC14")
    print(
        f"I={I}; J0={J0}; J11={J11}; "
        f"allocation_unit={allocation_unit}"
    )
    print(
        f"q11_full_factor_cells={len(full_cells)}; "
        "characterization=COMMON_S_x_{5,...,11}_x_{clock1}; "
        f"cheapest_witness={witness_cell}"
    )
    print(
        "weighted_carrier="
        f"(source_I={source_hit},target_J0={target_hit},"
        f"target_h11_J11={target_h11_hit},"
        f"pullback_to_I={pulled_hit})"
    )
    print(
        f"address_cylinders=(source={n_source}_mod13_{n_source % P},"
        f"target={n_target}_mod13_{n_target % P}); "
        "delayed_semantic=(source_carry12,target_carry6)="
        f"({source_semantic[12]},{target_semantic[6]})"
    )
    print(
        f"right_endpoint_rows={endpoint_rows}; "
        f"phase_quotient={phase_exponent // endpoint_base.endpoint.NN}"
    )
    print(
        f"ancestry_source_count={source_count}; "
        f"target_counts_by_q={tuple(target_counts)}; "
        f"common_counts_by_q={tuple(common_counts)}; "
        f"first_common_by_q={tuple(first_common)}"
    )
    print(
        "natural_lifted_common_weights_(q,h,carry,count)="
        f"{natural_lifted_common}; "
        f"carry_defect_lift_minus_same={carry_defect}; "
        f"carry_defect_mod13={tuple(value % P for value in carry_defect)}; "
        f"beta_arrow_abs_defect_gcd91={gcd(abs(carry_defect[-1]), 91)}; "
        f"reverse_q11_h9_count={reverse_q11_h9}; "
        "whole_cylinder=left_mid_right_agree"
    )
    print(
        f"E3_q_support={e3_support}; "
        f"THM2584_UQ_target_q_support={ancestry_support}; "
        "nonzero_intersection=empty"
    )
    print(
        "q11_outer_words="
        f"(QA_source={sum(qa_source)},QA_target={sum(qa_target_11)},"
        f"QB_source={source_count},QB_target={target_counts[Q]},"
        f"QAB_source={sum(qab_source)},"
        f"QAB_target={sum(qab_target_11)}); "
        "same_word_common=(0,0,0); "
        "directed_QB_to_QA="
        f"(count={len(cross_qb_to_qa_11)},"
        f"first={cross_qb_to_qa_11[0]},"
        f"last={cross_qb_to_qa_11[-1]},"
        "residue_mod169=156,whole_cylinder=1)"
    )
    print(
        "CONCLUSION: semantic reselection repairs the fixed-cell q2 loss.  "
        "At q=11, 63 lawful cells retain all six native factors; the "
        "cheapest cell (0,5,1) also retains equal weighted source/target "
        "carrier atoms, the residue-8/root cylinder, delayed carry/root "
        "data, and both right endpoint scalars.  The first failure is "
        "strictly earlier than endpoint transport: the target THM-2584 "
        "U_Q ancestry gate has zero active labels at q=11.  Across the "
        "whole allocation orbit, nonzero ancestry occurs only at q=7, "
        "where E3 fails, while E3 survives only at q=3,11.  Thus the "
        "fixed-atom 567-cell bank has an ancestry/normal transgression "
        "sidecar debt, not semantic-cell, carrier-weight, root, or endpoint "
        "repair.  The canonical lifted q11,h9 arrow remains positive with "
        "65098 common sheets after its +1 ancestry borrow, but canon does "
        "not supply that future action on the physical packet.  Static "
        "outer-word replacement also fails: the nearest positive q11 "
        "cospan has 449 whole-cylinder labels but changes QB={b} to QA={a}."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
