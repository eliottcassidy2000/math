#!/usr/bin/env python3
"""Finite-exact first-normal-lift test on the THM-3285 canonical horn.

The fixed horn source addresses are all congruent to seven modulo thirteen.
THM-2813's relative lift A_t therefore fixes each source address, while the
successor target address n+1 is on the adjacent residue-eight sheet and moves
by t*13^5.  In physical allocation coordinates this is the target-only shift

    q_alloc = 7*t mod 13.

For all thirteen labels this scout constructs the shifted target unit atom by
that inherited successor/transvection map.  It keeps the source atom fixed,
retains all six semantic gates, and enumerates the raw Boolean
bare/source/target/both K4 before any endpoint Fourier transform.  A candidate
may enter the 169-origin/two-field endpoint covariance scan only if all three
horn vertices have typed co-support and the raw K4 is nonflat at that common
support.

No candidate reaches that expensive scan.  The outer source semantic gates
are invariantly absent, and the strongest middle candidate t=6 (q=3) has all
six source/target gates but raw K4 (1,1,1,1) and mixed face zero at its unique
co-support.  The normal lift faithfully recovers t; it does not create the
missing endpoint difference q, determinant Delta, or a two-axis mixed face.

This is a bounded finite-exact partial scout, not a no-go for semantic-cell
reselection or arbitrary endpoint couplings, and not an LRC(14) result.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


P = 13
ADDRESS_MODULUS = P**6
TOP = P**5
N0 = 3_454_614
N_PLUS = N0 + P
N_A = N0 + 689_364
ADDRESSES = (N0, N_PLUS, N_A)
CENTRAL_INDICES = (0, 1, 53_028)
CLOCK = 1
SIGMA = 0
TAU = 3
EXPECTED_WEIGHT = 27_581_135_604
EXPECTED_MASS = Fraction(60_781_651_775_958_960, 371_293)
EXPECTED_COEFFICIENT = 790_161_473_087_466_480


PINNED = {
    "04-computation/lrc14_endpoint_origin_k4_horn_partial_scout_20260803.py":
        "5f1258cd19561941de93135b857bc4187435d10bfa0c1ca0be186c3d8e039e0c",
    "05-knowledge/results/lrc14_endpoint_origin_k4_horn_partial_scout_20260803.out":
        "af9c946d7832c8dfc6442846e7616fadb153cda33f4ef5386272671fe856cbe7",
    "04-computation/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py":
        "c42d66498f460f2142ea375fe9d4047b82c62b872b35d5a1634d2bb4c80a68ee",
    "05-knowledge/results/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.out":
        "e89dce3307e5d374e8583f92e1b2da1214e44929e52fdd42c6532d61adb3e246",
    "04-computation/lrc14_affine_lift_transvection_horn_decoder_thm2813.py":
        "255ce911a18f33d6b700213d6362886970b12170c324d39bed82a69a51b63e83",
    "05-knowledge/results/lrc14_affine_lift_transvection_horn_decoder_thm2813.out":
        "32f61740a0e7e73384c3d1ff1a83ba30d4cd3ebeca5d228b57e0bf2510928d58",
    "04-computation/lrc14_boolean_norm_cotangent_boundary_thm2820.py":
        "8f9a51e0fd616cd616eef0080bcd054b2a0c191704e1f289e78ea28364476376",
    "05-knowledge/results/lrc14_boolean_norm_cotangent_boundary_thm2820.out":
        "d9f5cf1e0999f414e0c9af838070441def6baf1c65df266fbdb44b9f578cce58",
    "04-computation/lrc14_residue8_common_allocation_covariance_thm2820.py":
        "779e2fab9b6aa80097b4d3756c32cdb040d4c2a2e9dd31ec9c7effcf11b780ae",
    "05-knowledge/results/lrc14_residue8_common_allocation_covariance_thm2820.out":
        "e4809c77178a3b66901e4b68eca517cfcb74e2516e16ece719cfca1efe4cbe7e",
    "04-computation/lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806.py":
        "7fea18161046f8de35b2e6ef04c88a13485de61045bc363f89c5ebfec8f76480",
    "05-knowledge/results/lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806.out":
        "cba545e7beff4fe889c76ae681c47c806969ccd33aa79def54527709a6ffafc6",
    "04-computation/lrc14_carrier_allocation_pullback_k4_segre_thm2772.py":
        "10b681f575fb51eb4b1af4bc909fba89846b85bd5da36fc069dac97ae2ebe409",
    "05-knowledge/results/lrc14_carrier_allocation_pullback_k4_segre_thm2772.out":
        "9cd1e82634822e4997d9966a8252b566d099050e42593a47c49b8a0f387190e5",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


ACTUAL_HASHES = tuple((path, lf_hash(ROOT / path)) for path in PINNED)
require(
    all(digest == PINNED[path] for path, digest in ACTUAL_HASHES),
    "a pinned inheritance artifact changed",
)


import lrc14_endpoint_origin_k4_horn_partial_scout_20260803 as frozen
import lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806 as tau12


physical = frozen.physical
allocation = frozen.allocation


def relative_lift(label, address):
    """THM-2813 relative lift A_label on Z/13^6."""
    return (
        address
        + label * TOP * ((address - 7) % P)
    ) % ADDRESS_MODULUS


def raw_k4(source_mask, target_mask):
    """Enumerate the inherited Boolean K4 before endpoint/carrier DFT.

    The corner order is (bare, source, target, both).  The first corner is
    the unprojected unit atom; the next two apply the source and target
    carrier selectors; the last applies both.  This is the B/P/Q/H square of
    THM-2820, not an endpoint-current coefficient square.
    """
    rows = []
    co_support = []
    mixed_support = []
    for source_harmonic in range(P):
        for target_harmonic in range(P):
            source_present = int(source_mask[source_harmonic])
            target_present = int(target_mask[target_harmonic])
            corners = (
                1,
                source_present,
                target_present,
                source_present * target_present,
            )
            mixed = corners[0] - corners[1] - corners[2] + corners[3]
            row = (
                (source_harmonic, target_harmonic),
                corners,
                mixed,
            )
            rows.append(row)
            if source_present and target_present:
                co_support.append(row)
            if mixed:
                mixed_support.append(row)
    return tuple(rows), tuple(co_support), tuple(mixed_support)


def endpoint_scan(source_atoms, target_atoms):
    """Full inherited scan, called only after every cheap gate survives."""
    endpoint_base = allocation.endpoint_base
    present_sets = endpoint_base.present_cache()
    terminal = tuple(
        endpoint_base.endpoint.build_set(
            endpoint_base.PAT_Q12, endpoint_base.ZERO
        )
    )
    q_starts = tuple(left for left, _right in terminal)
    tabs = tuple(
        endpoint_base.endpoint.make_tabs(
            terminal, endpoint_base.X0, (embedding,)
        )
        for embedding in endpoint_base.endpoint.MODS
    )
    context = (endpoint_base, present_sets, terminal, q_starts, tabs)
    fields = tuple(
        frozen.inverse_endpoint_banks(source_atom, target_atom, context)
        for source_atom, target_atom in zip(source_atoms, target_atoms)
    )
    require(
        len(fields) == 3
        and all(len(vertex) == 2 for vertex in fields)
        and all(
            len(left) == len(right) == P**2
            for vertex in fields
            for _prime, _zeta, left, right in vertex
        ),
        "a surviving endpoint scan did not cover 169 origins in two fields",
    )
    covariance = []
    for field_index, (prime, zeta, _left0, right0) in enumerate(fields[0]):
        _prime_plus, _zeta_plus, _left_plus, right_plus = fields[1][field_index]
        _prime_a, _zeta_a, _left_a, right_a = fields[2][field_index]
        require(
            all(right0.values())
            and all(right_plus.values())
            and all(right_a.values()),
            "a surviving target endpoint bank lost full support",
        )
        covariance.append((
            prime,
            frozen.affine_covariances(
                right0, right_plus, zeta, prime, endpoint_base.KEYS
            ),
            frozen.affine_covariances(
                right_plus, right_a, zeta, prime, endpoint_base.KEYS
            ),
            frozen.affine_covariances(
                right0, right_a, zeta, prime, endpoint_base.KEYS
            ),
        ))
    return fields, tuple(covariance)


def build_horn_geometry():
    (
        _module,
        _rails,
        _present,
        overlap_details,
        full_module,
        e3,
        clocks,
        q_pairs,
    ) = physical.build_physical_geometry()
    source_base, target_base = overlap_details[CLOCK]
    section = tuple(
        physical.physical.target.source_present_section(
            full_module, e3, CLOCK, SIGMA, TAU, clocks
        )
    )
    shifted_section = physical.shift_intervals(section, -physical.SHIFT)
    source_one = physical.intersect_weighted(source_base, section)
    target_one = physical.intersect_weighted(target_base, section)
    target_pullback = physical.shift_weighted(target_one, -physical.SHIFT)
    common_pullback = physical.intersect_weighted(
        target_pullback, physical.support_of(source_one)
    )
    right_pullback = physical.intersect_weighted(
        target_pullback,
        physical.complement(physical.support_of(source_one)),
    )
    target_base_pullback = physical.shift_weighted(
        target_base, -physical.SHIFT
    )
    require(
        target_pullback
        == physical.intersect_weighted(target_base_pullback, shifted_section),
        "target section stopped commuting with pullback",
    )
    common = physical.shift_weighted(common_pullback, physical.SHIFT)
    right = physical.shift_weighted(right_pullback, physical.SHIFT)

    cylinders = tuple(
        frozen.target_cylinder(index) for index in CENTRAL_INDICES
    )
    common_cells = tuple(
        physical.intersect_weighted(common, cylinder)
        for cylinder in cylinders
    )
    right_cells = tuple(
        physical.intersect_weighted(right, cylinder)
        for cylinder in cylinders
    )
    positive_cells = tuple(
        common_cell or right_cell
        for common_cell, right_cell in zip(common_cells, right_cells)
    )
    require(
        tuple(bool(cell) for cell in common_cells) == (False, True, False)
        and tuple(bool(cell) for cell in right_cells) == (True, False, True),
        "canonical R-M-R word changed",
    )
    require(
        all(
            frozen.weighted_stats(cell, q_pairs)
            == (1, EXPECTED_MASS, EXPECTED_COEFFICIENT)
            and cell[0][2] == EXPECTED_WEIGHT
            for cell in positive_cells
        ),
        "a canonical horn cell changed",
    )

    target_atoms = []
    for cell in positive_cells:
        left, right_endpoint, _weight = cell[0]
        integer_left = frozen.ceil_fraction(left)
        require(
            left < integer_left < integer_left + 1 < right_endpoint,
            "a horn cell lost its canonical integral unit atom",
        )
        target_atoms.append(((integer_left, integer_left + 1, 1),))
    target_atoms = tuple(target_atoms)
    source_atoms = tuple(
        physical.shift_weighted(atom, -physical.SHIFT)
        for atom in target_atoms
    )
    require(
        target_atoms == (
            ((142_082_432_180_573, 142_082_432_180_574, 1),),
            ((142_088_047_310_093, 142_088_047_310_094, 1),),
            ((142_004_622_528_653, 142_004_622_528_654, 1),),
        ),
        "canonical target unit atoms changed",
    )
    factors = allocation.section_factors(
        full_module, e3, clocks, SIGMA, TAU, CLOCK
    )
    return {
        "full_module": full_module,
        "e3": e3,
        "clocks": clocks,
        "source_base": source_base,
        "target_base": target_base,
        "cylinders": cylinders,
        "source_atoms": source_atoms,
        "target_atoms": target_atoms,
        "factors": factors,
    }


def tau12_control():
    """Rebuild the inherited target-only hostile before any normal label."""
    context = tau12.arm.build_context()
    module, _delayed, _present, source, clocks, _weight, rail_common = context
    objects = tau12.arm.carrier_objects(
        module, source, clocks, rail_common, CLOCK, SIGMA, 12
    )
    counts = tuple(len(objects[key]) for key in ("A", "B", "M", "L", "R"))
    require(
        counts == (240, 241, 0, 240, 241)
        and not objects["M"]
        and objects["L"] == objects["A"]
        and objects["R"] == objects["B"],
        "tau-twelve target-only control changed",
    )
    return counts


def main():
    geometry = build_horn_geometry()
    full_module = geometry["full_module"]
    source_base = geometry["source_base"]
    target_base = geometry["target_base"]
    cylinders = geometry["cylinders"]
    source_atoms = geometry["source_atoms"]
    target_atoms = geometry["target_atoms"]
    factors = geometry["factors"]
    unit = full_module.T // P
    endpoint_origins = tuple(allocation.endpoint_base.KEYS)
    require(
        full_module.T % P == 0
        and pow(7, -1, P) == 2
        and len(endpoint_origins) == P**2
        and set(endpoint_origins)
        == {(left, right) for left in range(P) for right in range(P)},
        "allocation-unit, normal decoder, or endpoint-origin plane changed",
    )

    factor_names = ("E3", "clock", "q1", "q2", "c2", "c3")
    source_signatures = tuple(
        tuple(
            allocation.contained(atom[0][:2], factor)
            for factor in factors
        )
        for atom in source_atoms
    )
    require(
        source_signatures == (
            (False, True, True, True, False, True),
            (True, True, True, True, True, True),
            (False, True, True, True, True, True),
        ),
        "base source signatures stopped reproducing the frozen scout",
    )

    normal_rows = []
    target_signatures_by_vertex = [[], [], []]
    middle_local_survivors = []
    endpoint_entrants = []
    endpoint_scan_rows = []
    raw_table_checks = 0
    raw_origin_copies_checked = 0

    for label in range(P):
        q_alloc = 7 * label % P
        require(2 * q_alloc % P == label, "normal label decoder changed")
        shifted_targets = []
        label_rows = []
        for vertex, (address, source_atom, target_atom, cylinder) in enumerate(
            zip(ADDRESSES, source_atoms, target_atoms, cylinders)
        ):
            source_address = address
            target_address = (address + 1) % ADDRESS_MODULUS
            source_image = relative_lift(label, source_address)
            target_image = relative_lift(label, target_address)
            require(
                source_address % P == 7
                and target_address % P == 8
                and source_image == source_address
                and target_image
                == (target_address + label * TOP) % ADDRESS_MODULUS,
                "the successor/transvection address square changed",
            )
            require(
                (
                    target_image // TOP
                    - target_address // TOP
                ) % P == label,
                "the adjacent-sheet top-digit cocycle stopped decoding t",
            )
            physical_shift = q_alloc * unit
            require(
                (7 * label * unit) % full_module.T == physical_shift,
                "address normal displacement stopped matching q=7t",
            )
            shifted_target = physical.shift_weighted(
                target_atom, physical_shift
            )
            shifted_cylinder = physical.shift_intervals(
                cylinder, physical_shift
            )
            require(
                physical.intersect_weighted(
                    shifted_target, shifted_cylinder
                ) == shifted_target
                and physical.shift_weighted(
                    shifted_target, -physical_shift
                ) == target_atom,
                "normal target atom was not built by the inherited shift",
            )
            shifted_targets.append(shifted_target)

            target_signature = tuple(
                allocation.contained(shifted_target[0][:2], factor)
                for factor in factors
            )
            target_signatures_by_vertex[vertex].append(target_signature)
            source_mask = frozen.carrier_mask(
                source_atom[0][:2], source_base, -1, unit
            )
            target_mask = frozen.carrier_mask(
                shifted_target[0][:2], target_base, +1, unit
            )
            require(
                tuple(index for index, bit in enumerate(source_mask) if bit)
                == (0,)
                and tuple(
                    index for index, bit in enumerate(target_mask) if bit
                ) == (q_alloc,),
                "normal carrier masks stopped being delta_0/delta_q",
            )

            table, co_support, mixed_support = raw_k4(
                source_mask, target_mask
            )
            raw_table_checks += len(table)
            raw_origin_copies_checked += len(table) * len(endpoint_origins)
            require(
                len(table) == P**2
                and len(co_support) == 1
                and co_support[0]
                == ((0, q_alloc), (1, 1, 1, 1), 0)
                and len(mixed_support) == (P - 1) ** 2
                and all(row[0] != (0, q_alloc) for row in mixed_support),
                "raw K4 co-support/mixed-face boundary changed",
            )

            source_gate = all(source_signatures[vertex])
            target_gate = all(target_signature)
            semantic_cosupport = source_gate and target_gate
            nonflat_common_face = co_support[0][2] != 0
            label_rows.append((
                vertex,
                source_address,
                target_address,
                target_image,
                tuple(int(bit) for bit in source_signatures[vertex]),
                tuple(int(bit) for bit in target_signature),
                co_support[0][0],
                co_support[0][1],
                co_support[0][2],
                semantic_cosupport,
                nonflat_common_face,
            ))

        shifted_targets = tuple(shifted_targets)
        outer_semantic_cosupport = label_rows[0][9] and label_rows[2][9]
        middle_semantic_cosupport = label_rows[1][9]
        all_raw_common_faces_nonflat = all(row[10] for row in label_rows)
        if middle_semantic_cosupport:
            middle_local_survivors.append((label, q_alloc))
        if (
            outer_semantic_cosupport
            and middle_semantic_cosupport
            and all_raw_common_faces_nonflat
        ):
            endpoint_entrants.append((label, q_alloc))
            fields, covariance = endpoint_scan(source_atoms, shifted_targets)
            endpoint_scan_rows.append((label, q_alloc, fields, covariance))
        normal_rows.append((
            label,
            q_alloc,
            tuple(row[3] for row in label_rows),
            tuple(row[5] for row in label_rows),
            outer_semantic_cosupport,
            middle_semantic_cosupport,
            all_raw_common_faces_nonflat,
        ))

    target_signatures_by_vertex = tuple(
        tuple(rows) for rows in target_signatures_by_vertex
    )
    require(
        all(
            tuple(
                label for label, signature in enumerate(vertex_rows)
                if all(signature)
            ) == (0, 6)
            for vertex_rows in target_signatures_by_vertex
        )
        and all(
            vertex_rows[6] == (True,) * 6
            for vertex_rows in target_signatures_by_vertex
        )
        and target_signatures_by_vertex[0][9]
        == (True, True, True, False, True, True)
        and target_signatures_by_vertex[1][9]
        == (False, True, True, False, True, True)
        and target_signatures_by_vertex[2][9]
        == (True, True, True, False, True, True),
        "q=3/q=11 semantic controls changed",
    )
    require(
        middle_local_survivors == [(0, 0), (6, 3)]
        and not endpoint_entrants
        and not endpoint_scan_rows,
        "the staged endpoint-scan boundary changed",
    )

    # Label zero must reproduce the pre-DFT portion of the frozen scout.
    label_zero = normal_rows[0]
    require(
        label_zero[1] == 0
        and label_zero[3]
        == tuple(tuple(int(bit) for bit in signature)
                 for signature in ((True,) * 6,) * 3)
        and not label_zero[4]
        and label_zero[5]
        and not label_zero[6],
        "normal label zero stopped reproducing the current horn failure",
    )
    frozen_output = (
        RESULTS
        / "lrc14_endpoint_origin_k4_horn_partial_scout_20260803.out"
    ).read_text(encoding="utf-8")
    require(
        "middle_endpoint_fibre=all_169_left_and_right_origins_nonzero_in_both_fields;"
        in frozen_output
        and "raw_twist00_vector_is_flat_and_has_Mobius_face_zero"
        in frozen_output
        and "VERDICT=OUTCOME_2" in frozen_output,
        "the pinned t=0 downstream control lost its stated failure",
    )

    # THM-2820's q=3/q=11 either-or belongs to a distinct fixed semantic
    # cell.  It is retained only as an external hostile/positive control.
    thm2820_output = (
        RESULTS
        / "lrc14_residue8_common_allocation_covariance_thm2820.out"
    ).read_text(encoding="utf-8")
    require(
        "unique_nonzero_full_fixed_boundary=(t=6,q=3)" in thm2820_output
        and "factor_full_q=(0, 3); two_point_endpoint_edge_q=(0, 11);"
        in thm2820_output
        and "nonzero_intersection=empty" in thm2820_output,
        "THM-2820 exclusive-cell control changed",
    )

    tau12_counts = tau12_control()
    tau12_output = (
        RESULTS
        / "lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806.out"
    ).read_text(encoding="utf-8")
    require(
        "M_empty_before_address_restriction=yes" in tau12_output,
        "tau-twelve stored hostile changed",
    )

    source_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(
        isinstance(node, ast.Assert) for node in ast.walk(source_tree)
    )
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(source_tree)
    )
    require(
        assert_nodes == 0 and float_literals == 0,
        "source validity gate found assert or floating-point literals",
    )

    row_digest = sha256(repr(tuple(normal_rows)).encode("utf-8")).hexdigest()

    print("LRC14 CANONICAL HORN FIRST-NORMAL-LIFT PARTIAL SCOUT")
    print("status=FINITE-EXACT PARTIAL;not_a_canonical_theorem")
    print("dependency_hashes=BEGIN")
    for path, digest in ACTUAL_HASHES:
        print(f"{digest}  {path}")
    print("dependency_hashes=END")
    print(
        "universe=rail8;(physical_clock,sigma,tau)=(1,0,3);"
        f"source_addresses={ADDRESSES};source_residue=7;"
        "successor_target_residue=8;normal_labels=F13;"
        "endpoint_origin_sidecar=F13^2_but_not_DFT_evaluated_without_gate_survivor;"
        "certified_fields_available=2"
    )
    print(
        "lawful_normal_map=A_t(n)=n+t*13^5*(n-7_mod13);"
        "source_addresses_fixed;target_successors_move_by_t*13^5;"
        "physical_target_shift=q_alloc*T/13;q_alloc=7t_mod13;"
        "normal_decoder=t=2*q_alloc_mod13"
    )
    print(
        f"factor_order={factor_names};source_signatures={source_signatures};"
        "target_signatures_by_vertex_then_t="
        f"{tuple(tuple(tuple(int(bit) for bit in signature) for signature in vertex_rows) for vertex_rows in target_signatures_by_vertex)}"
    )
    print(
        "raw_K4_order=(bare,source,target,both);"
        "formula_at_harmonics_(k,l)=(1,delta0(k),delta_q(l),delta0(k)*delta_q(l));"
        "unique_raw_cosupport=(0,q)_with_K4_(1,1,1,1)_and_mixed_face_0;"
        "raw_mixed_aggregate_support=144_joint_absent_cells_disjoint_from_cosupport"
    )
    print(
        f"raw_table_checks={raw_table_checks};"
        f"raw_table_x_endpoint_origin_cross_product_size={raw_origin_copies_checked};"
        "raw_flatness_is_uniform_in_all_169_endpoint_origins_before_DFT"
    )
    print(
        "normal_rows_(t,q,target_images,target_factor_signatures_at_3_vertices,"
        "outer_semantic_cosupport,middle_semantic_cosupport,"
        f"all_common_faces_nonflat)={tuple(normal_rows)}"
    )
    print(f"normal_row_digest={row_digest}")
    print(
        f"middle_local_six_gate_survivors={tuple(middle_local_survivors)};"
        "strongest=(t=6,q_alloc=3):all_six_source_and_target_gates_at_middle;"
        "raw_cosupport=(0,3);raw_K4=(1,1,1,1);mixed_face=0"
    )
    print(
        "outer_semantic_cosupport_survivors=0_of_13;"
        "n0_source_invariantly_loses_(E3,c2);"
        "na_source_invariantly_loses_E3;"
        "A_t_fixes_the_residue7_source_sheet"
    )
    print(
        "t0_control=reproduces_frozen_pre_DFT_failure;"
        "inherited_downstream_control=middle_169_origins_present_but_raw_common_face_flat_"
        "and_vertical_affine_covariance_absent"
    )
    print(
        "THM2820_external_exclusive_cell_control=q3_keeps_all_six_fixed_factors_"
        "but_loses_translated_endpoint_companion;q11_keeps_two_point_endpoint_edge_"
        "but_loses_q2;different_(sigma,tau)_cell_not_folded_into_this_horn"
    )
    print(
        f"tau12_target_only_control_(A,B,M,L,R)_piece_counts={tau12_counts};"
        "M_empty_before_address_or_normal_label"
    )
    print(
        f"endpoint_scan_entrants={tuple(endpoint_entrants)};"
        f"endpoint_169x2_covariance_scan_calls={len(endpoint_scan_rows)};"
        "cheap_gate_correctly_stopped_all_candidates"
    )
    print(
        "FIRST_HORN_GLOBAL_FAILURE=the_normal_lift_moves_only_the_target_successor_"
        "and_cannot_restore_the_outer_source_E3_semantic_gates"
    )
    print(
        "FIRST_STRONG_LOCAL_FAILURE=one_faithful_target_axis_normal_translation_"
        "only_reindexes_delta_q_and_preserves_raw_common_K4_flatness;"
        "it_does_not_supply_two_independent_endpoint_translations_or_Delta"
    )
    print(
        "STRONGEST_SURVIVOR=t_is_recovered_exactly_from_q_alloc_and_t6_has_a_"
        "fully_typed_middle_successor_pair_on_the_transported_raw_carrier_copy;"
        "the_endpoint_difference_q_E=(L-R),endpoint_origin_(L,R),and_Delta=det(L,R)_"
        "remain_unconstructed_and_must_not_be_identified_with_q_alloc"
    )
    print(
        "RESELECTION_BOUNDARY=neighboring_semantic_cells_may_change_the_q3_q11_"
        "factor_endpoint_tradeoff_but_are_a_separate_typed_branch_requiring_fresh_"
        "carrier_weight_ancestry_address_and_endpoint_audits"
    )
    print(
        "scope=canonical_integral_unit_subatoms_and_all_13_inherited_first_normal_"
        "labels;no_claim_against_semantic_cell_reselection_or_arbitrary_endpoint_"
        "couplings;no_root_Cech_map;no_row_exclusion;no_LRC14"
    )
    print(
        f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
