#!/usr/bin/env python3
"""Independent hostile audit of the THM-2806 tau-twelve addendum.

This script does not import the primary addendum or its ``carrier_objects``
constructor.  It starts from the independently audited one-sided source and
target carriers, inserts the full semantic section, pulls the target object
back to the source chart, and forms A∩B, A\\B, B\\A directly.  It then
restricts the target right wing to the three exact THM-2807 cylinders.

The separate 81-label bank is rebuilt with the same direct set operations.
No endpoint-current or Fourier code is used.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

P = 13
Q = 7
ADDRESS_MODULUS = P**6
TAU = Fraction(7, ADDRESS_MODULUS)
N0 = 3_454_614
N_PLUS = N0 + P
N_A = N0 + 689_364
CENTRAL_INDICES = (0, 1, 53_028)
EXPECTED_WEIGHT = 27_581_135_604
EXPECTED_MASS = Fraction(60_781_651_775_958_960, 371_293)
EXPECTED_COEFFICIENT = 790_161_473_087_466_480
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)
EXPECTED_COMMON_CLOCKS = (0, 81, 56, 56, 0, 0, 0)

PINNED = {
    "lrc14_opposite_wing_common_atom_virtual_sector_audit_thm2782.py":
        "f696867c5242e86593de64f3ee7944bf44a3887e7818291ed059c7358f959c1a",
    "lrc14_positive_graded_address_two_simplex_thm2807.py":
        "11cdbe3c6cc7f9d5b6b24863ced71eb91cc84adc67fe38a3f8a3e637362453fb",
    "lrc14_literal_fixed_sheet_allocation_independent_audit_thm2806.py":
        "90386206fa441edaa121b54b501347f14005e997b4c0e8f95947f2fac14050b4",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for dependency, expected in PINNED.items():
    require(
        lf_hash(COMP / dependency) == expected,
        f"pinned dependency changed: {dependency}",
    )


import lrc14_opposite_wing_common_atom_virtual_sector_audit_thm2782 as audit


def frac(value):
    return value - value.numerator // value.denominator


def target_cylinder(central_index):
    center = frac(
        audit.physical.relative.Z
        + (N0 + P * central_index + 1) * TAU
    )
    coordinate = center * audit.PERIOD
    half = audit.physical.relative.Q_RADIUS * audit.PERIOD
    require(
        0 < coordinate - half < coordinate + half < audit.PERIOD,
        "independent target cylinder crossed the circle seam",
    )
    return ((coordinate - half, coordinate + half),)


def main():
    (
        _module,
        _rails,
        _present,
        overlap_details,
        full_module,
        e3,
        clock_combs,
        q_pairs,
    ) = audit.build_physical_geometry()

    def build_face(clock, sigma, tau):
        source_base, target_base = overlap_details[clock]
        target_pullback = audit.shift_weighted(
            target_base, -audit.SHIFT
        )
        section = tuple(
            audit.physical.target.source_present_section(
                full_module, e3, clock, sigma, tau, clock_combs
            )
        )
        shifted_section = audit.shift_intervals(section, -audit.SHIFT)

        source_one = audit.intersect_weighted(source_base, section)
        target_one = audit.intersect_weighted(target_base, section)
        target_one_pullback = audit.shift_weighted(
            target_one, -audit.SHIFT
        )
        target_direct_pullback = audit.intersect_weighted(
            target_pullback, shifted_section
        )
        require(
            target_one_pullback == target_direct_pullback,
            "target pullback stopped commuting with section restriction",
        )

        common_source = audit.intersect_weighted(
            source_one, audit.support_of(target_one_pullback)
        )
        common_target_pullback = audit.intersect_weighted(
            target_one_pullback, audit.support_of(source_one)
        )
        source_wing = audit.intersect_weighted(
            source_one, audit.complement(audit.support_of(target_one_pullback))
        )
        target_wing_pullback = audit.intersect_weighted(
            target_one_pullback,
            audit.complement(audit.support_of(source_one)),
        )
        target_wing = audit.shift_weighted(
            target_wing_pullback, audit.SHIFT
        )
        require(
            audit.disjoint_union(common_source, source_wing) == source_one
            and audit.disjoint_union(
                common_target_pullback, target_wing_pullback
            ) == target_one_pullback
            and not audit.intersect_intervals(
                audit.support_of(common_source),
                audit.support_of(target_wing_pullback),
            ),
            "independent one-sided decomposition failed",
        )
        return {
            "A": source_one,
            "B": target_one_pullback,
            "M_source": common_source,
            "M_target": common_target_pullback,
            "L": source_wing,
            "R": target_wing_pullback,
            "R_target": target_wing,
        }

    face12 = build_face(1, 0, 12)
    full_section_counts = tuple(
        len(face12[key]) for key in ("A", "B", "M_source", "L", "R")
    )
    require(
        full_section_counts == (0, 241, 0, 0, 241)
        and not face12["M_source"]
        and not face12["M_target"]
        and not face12["A"]
        and face12["B"] == face12["R"],
        "independent tau-twelve target-only face changed",
    )

    selected_rows = []
    target_cells = []
    for address, central_index in zip(
        (N0, N_PLUS, N_A), CENTRAL_INDICES
    ):
        cylinder = target_cylinder(central_index)
        cell = audit.intersect_weighted(
            face12["R_target"], cylinder
        )
        mass = audit.weighted_mass(cell)
        coefficient = audit.physical.relative.private.delayed_carry_pair(
            cell, q_pairs[1], {}
        )[6][1]
        require(
            len(cell) == 1
            and cell[0][:2] == cylinder[0]
            and cell[0][2] == EXPECTED_WEIGHT
            and mass == EXPECTED_MASS
            and coefficient == EXPECTED_COEFFICIENT,
            f"independent tau-twelve positive cell changed at {address}",
        )
        target_cells.append(tuple(cell))
        selected_rows.append((
            address, central_index, len(cell), mass, coefficient
        ))

    for index in (1, 2):
        address_gap = (N_PLUS, N_A)[index - 1] - N0
        require(
            audit.shift_weighted(
                target_cells[0], address_gap * audit.SHIFT
            )
            == target_cells[index],
            "independent whole-cylinder translation failed",
        )

    face3 = build_face(1, 0, 3)
    require(
        face3["M_source"] and not face12["M_source"],
        "independent common-positive/tau-twelve hostile changed",
    )
    table3 = []
    table12 = []
    for central_index in CENTRAL_INDICES:
        cylinder = target_cylinder(central_index)
        for table, face in ((table3, face3), (table12, face12)):
            cell = audit.intersect_weighted(face["R_target"], cylinder)
            table.append((
                audit.weighted_mass(cell),
                audit.physical.relative.private.delayed_carry_pair(
                    cell, q_pairs[1], {}
                )[6][1],
            ))
    positive = (EXPECTED_MASS, EXPECTED_COEFFICIENT)
    require(
        tuple(table3) == (positive, (0, 0), positive)
        and tuple(table12) == (positive, positive, positive),
        "independent tau-three/tau-twelve table changed",
    )

    common_clock_counts = []
    right_clock_counts = []
    paired_support_checks = 0
    disjoint_checks = 0
    pair_type_histogram = Counter()
    for clock in range(Q):
        common_count = 0
        right_count = 0
        for sigma in COMMON_S:
            for tau in COMMON_T:
                face = build_face(clock, sigma, tau)
                common_present = bool(face["M_source"])
                right_present = bool(face["R"])
                common_count += common_present
                right_count += right_present
                paired_support_checks += common_present == right_present
                pair_type_histogram[
                    (common_present, right_present)
                ] += 1
                require(
                    not audit.intersect_intervals(
                        audit.support_of(face["M_source"]),
                        audit.support_of(face["R"]),
                    ),
                    "independent common/right supports met",
                )
                disjoint_checks += 1
        common_clock_counts.append(common_count)
        right_clock_counts.append(right_count)

    common_clock_counts = tuple(common_clock_counts)
    right_clock_counts = tuple(right_clock_counts)
    require(
        common_clock_counts == EXPECTED_COMMON_CLOCKS
        and right_clock_counts == EXPECTED_COMMON_CLOCKS
        and paired_support_checks == disjoint_checks == 567
        and pair_type_histogram
        == Counter({(False, False): 374, (True, True): 193})
        and 12 not in COMMON_T,
        "independent separate-bank census changed",
    )

    source_text = Path(__file__).read_text()
    tree = ast.parse(source_text)
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "optimized-mode validity gate contains assert",
    )

    print("THM-2806 TAU-TWELVE SIMPLEX INDEPENDENT HOSTILE AUDIT")
    print("status=FINITE-EXACT independent audit")
    print(f"dependency_sha256={tuple(PINNED.items())}")
    print(
        "implementation=full_one_sided_source_and_target_carriers;"
        "insert_full_semantic_section;pull_target_back;direct_set_differences;"
        "no_primary_addendum_or_carrier_objects_import"
    )
    print(
        "universe=(clock,sigma,tau)=(1,0,12);"
        f"addresses=({N0},{N_PLUS},{N_A});rail8;target_root1;"
        "delayed_carry6"
    )
    print(
        f"tau12_full_section_piece_counts_(A,B,M,L,R)="
        f"{full_section_counts};"
        "source_empty_target_only_before_cylinder=yes"
    )
    for row in selected_rows:
        print(
            f"tau12_vertex=address:{row[0]};central_index:{row[1]};"
            f"pieces:{row[2]};mass:{row[3]};coefficient:{row[4]}"
        )
    print(
        f"independent_selected_tables=tau3:{tuple(table3)};"
        f"tau12:{tuple(table12)};whole_cylinder_translation=3/3"
    )
    print(
        f"separate_common_bank_clock_counts={common_clock_counts};"
        f"right_cofiber_clock_counts={right_clock_counts};"
        f"pair_type_histogram={tuple(sorted(pair_type_histogram.items()))}"
    )
    print(
        f"cellwise_nonempty_support_agreement={paired_support_checks}/567;"
        f"common_vs_right_disjoint={disjoint_checks}/567;"
        "tau12_not_in_common_target_bank=yes"
    )
    print(
        "VERDICT=ACCEPT:the_THM2807_tau12_simplex_has_no_common_"
        "allocation_atom_and_is_literally_target_only_in_the_independent_"
        "full_section_constructor"
    )
    print(
        "SCOPE=item1_for_this_fixed_simplex_only;no_endpoint_origin;"
        "no_Fourier_or_root_deck_map;no_row_exclusion;no_LRC14"
    )
    print("ALL INDEPENDENT CHECKS PASSED")


if __name__ == "__main__":
    main()
