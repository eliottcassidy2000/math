#!/usr/bin/env python3
"""Exact scout for a same-ancestry allocation-switching LRC horn.

The universe is the genuine THM-2806 common-label bank at physical clock one,
on the fixed rail-eight sheet:

    sigma in {0,1,2,3,8,9,10,11,12},
    tau   in {3,4,5,6,7,8,9,10,11}.

For each label, form the full one-sided source carrier A and pulled target
carrier B before address restriction.  In the source chart put

    M=A intersect B,                 R=B minus A.

Push the B-side copies of M and R to the target chart, impose the inherited
root-one/carry-six address cylinders, and inspect the three THM-2807 vertices
n_0, n_+=n_0+13, n_a=n_0+689364.

The proof path is deliberately redundant.  Constructor A starts from the
independently audited full one-sided weighted carriers and uses direct set
operations.  Constructor B starts from the canonical unweighted
``carrier_objects`` decomposition and applies the physical restrictions in
the inherited order.  All M/R address cells must agree exactly.

The literal rail ancestry is audited separately by rebuilding all contributor
walls from THM-2791's path-sheet definitions.  The three cylinders lie in one
wall-free chamber, and the complete middle contributor set is enumerated.

This is a finite-exact carrier result.  It constructs no endpoint origin,
endpoint current, K4 allocation square, factorwise covariance, row exclusion,
or LRC(14) proof.
"""

from __future__ import annotations

import ast
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

P = 13
ADDRESS_MODULUS = P**6
N0 = 3_454_614
N_PLUS = N0 + P
N_A = N0 + 689_364
ADDRESSES = (N0, N_PLUS, N_A)
CENTRAL_INDICES = (0, 1, 53_028)
SIGMAS = (0, 1, 2, 3, 8, 9, 10, 11, 12)
TAUS = (3, 4, 5, 6, 7, 8, 9, 10, 11)
EXPECTED_HORN_SIGMAS = (0, 3, 8, 9, 10, 11, 12)
EXPECTED_WEIGHT = 27_581_135_604
EXPECTED_MASS = Fraction(60_781_651_775_958_960, 371_293)
EXPECTED_COEFFICIENT = 790_161_473_087_466_480
EXPECTED_PATH_DIGEST = (
    "15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd"
)

PINNED = {
    "04-computation/"
    "lrc14_opposite_wing_common_atom_virtual_sector_audit_thm2782.py":
        "f696867c5242e86593de64f3ee7944bf44a3887e7818291ed059c7358f959c1a",
    "04-computation/"
    "lrc14_semantic_arm_right_wing_central_digit_thm2782.py":
        "7fbc6bb1ec303ded98eaad6e5d8205eb3d247258ada32b6f9904fc439ebb11fb",
    "04-computation/lrc14_full_arm_orbit_path_sheet_audit_thm2791.py":
        "1e00b6711db0d878fca70047f5f1532518084dbf6928551cd28fe51283fde543",
    "05-knowledge/results/lrc14_full_arm_orbit_path_sheet_audit_thm2791.out":
        "2325aa126dc4a97af4ba4bde0348389bfdfe5720dec6f88792b9a06baa40afd3",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


ACTUAL_HASHES = tuple((path, lf_hash(ROOT / path)) for path in PINNED)
require(
    all(digest == PINNED[path] for path, digest in ACTUAL_HASHES),
    "a pinned carrier or ancestry dependency changed",
)


import lrc14_full_arm_orbit_path_sheet_audit_thm2791 as sheet
import lrc14_opposite_wing_common_atom_virtual_sector_audit_thm2782 as direct
import lrc14_semantic_arm_right_wing_central_digit_thm2782 as canonical


def frac(value):
    return value - value.numerator // value.denominator


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


TARGET_STEP = Fraction(7, ADDRESS_MODULUS)


def target_cylinder(central_index):
    center = frac(
        direct.physical.relative.Z
        + (N0 + P * central_index + 1) * TARGET_STEP
    )
    coordinate = center * direct.PERIOD
    half = direct.physical.relative.Q_RADIUS * direct.PERIOD
    require(
        0 < coordinate - half < coordinate + half < direct.PERIOD,
        "selected target cylinder crossed the circle seam",
    )
    return ((coordinate - half, coordinate + half),)


(
    _direct_module,
    _direct_rails,
    _direct_present,
    overlap_details,
    full_module,
    e3,
    clock_combs,
    q_pairs,
) = direct.build_physical_geometry()


def build_direct_face(clock, sigma, tau):
    """Build A, B, M and R from full one-sided weighted carriers."""
    source_base, target_base = overlap_details[clock]
    target_base_pullback = direct.shift_weighted(
        target_base, -direct.SHIFT
    )
    section = tuple(
        direct.physical.target.source_present_section(
            full_module, e3, clock, sigma, tau, clock_combs
        )
    )
    shifted_section = direct.shift_intervals(section, -direct.SHIFT)

    source_one = direct.intersect_weighted(source_base, section)
    target_one = direct.intersect_weighted(target_base, section)
    target_one_pullback = direct.shift_weighted(
        target_one, -direct.SHIFT
    )
    target_direct_pullback = direct.intersect_weighted(
        target_base_pullback, shifted_section
    )
    require(
        target_one_pullback == target_direct_pullback,
        "target pullback stopped commuting with section restriction",
    )

    common_source = direct.intersect_weighted(
        source_one, direct.support_of(target_one_pullback)
    )
    common_target_pullback = direct.intersect_weighted(
        target_one_pullback, direct.support_of(source_one)
    )
    source_wing = direct.intersect_weighted(
        source_one,
        direct.complement(direct.support_of(target_one_pullback)),
    )
    right_pullback = direct.intersect_weighted(
        target_one_pullback,
        direct.complement(direct.support_of(source_one)),
    )
    require(
        direct.disjoint_union(common_source, source_wing) == source_one
        and direct.disjoint_union(
            common_target_pullback, right_pullback
        ) == target_one_pullback
        and direct.support_of(common_source)
        == direct.support_of(common_target_pullback)
        and not direct.intersect_intervals(
            direct.support_of(common_target_pullback),
            direct.support_of(right_pullback),
        ),
        f"direct allocation decomposition failed at {(clock, sigma, tau)}",
    )
    return {
        "B": target_one,
        "M": direct.shift_weighted(
            common_target_pullback, direct.SHIFT
        ),
        "R": direct.shift_weighted(right_pullback, direct.SHIFT),
    }


canonical_context = canonical.build_context()
(
    canonical_module,
    canonical_delayed,
    canonical_present,
    canonical_source,
    canonical_clocks,
    canonical_source_weight,
    canonical_rail_common,
) = canonical_context
canonical_source_safe = canonical.marked.complement(
    canonical_present[1, 7]
)
canonical_target_safe = canonical.marked.complement(
    canonical.marked.shift_union(
        canonical_present[1, 7], canonical.marked.SHIFT
    )
)


def canonical_cell(objects, key, central_index):
    """Apply the inherited filters to canonical M or R."""
    carrier = canonical.marked.intersect(
        objects[key], canonical_source_safe
    )
    carrier = canonical.marked.intersect(
        carrier, canonical_target_safe
    )
    carrier = canonical.marked.shift_union(
        carrier, -canonical.marked.SHIFT
    )
    weighted = canonical.marked.restrict_weighted(
        canonical_source_weight, carrier
    )
    weighted = canonical.marked.private.old.intersect_weighted_comb(
        weighted, canonical_module.C3, 182, 1, 13
    )
    weighted = canonical.marked.restrict_weighted(
        weighted, canonical.target_cylinder(N0, central_index)
    )
    return tuple(weighted)


def cell_stats(cell):
    coefficient = direct.physical.relative.private.delayed_carry_pair(
        cell, q_pairs[1], {}
    )[6][1]
    return len(cell), direct.weighted_mass(cell), coefficient


records = {}
pattern_labels = defaultdict(list)
cross_implementation_checks = 0
positive_stats = set()
pattern_digest = sha256()

for sigma in SIGMAS:
    for tau in TAUS:
        direct_face = build_direct_face(1, sigma, tau)
        canonical_objects = canonical.carrier_objects(
            canonical_module,
            canonical_source,
            canonical_clocks,
            canonical_rail_common,
            1,
            sigma,
            tau,
        )
        cells = {key: [] for key in ("B", "M", "R")}
        for central_index in CENTRAL_INDICES:
            cylinder = target_cylinder(central_index)
            require(
                cylinder
                == canonical.target_cylinder(N0, central_index),
                "the two target-cylinder definitions diverged",
            )
            for key in ("B", "M", "R"):
                cells[key].append(
                    direct.intersect_weighted(direct_face[key], cylinder)
                )
            for key in ("M", "R"):
                comparison = canonical_cell(
                    canonical_objects, key, central_index
                )
                require(
                    cells[key][-1] == comparison,
                    f"constructors disagree at {(sigma, tau, key, central_index)}",
                )
                cross_implementation_checks += 1

        m_bits = tuple(bool(cell) for cell in cells["M"])
        r_bits = tuple(bool(cell) for cell in cells["R"])
        b_bits = tuple(bool(cell) for cell in cells["B"])
        pattern = (m_bits, r_bits)
        pattern_labels[pattern].append((sigma, tau))
        pattern_digest.update(
            (
                f"{sigma},{tau}:"
                f"M{''.join(str(int(bit)) for bit in m_bits)}:"
                f"R{''.join(str(int(bit)) for bit in r_bits)};"
            ).encode("ascii")
        )

        for index in range(3):
            require(
                not (
                    cells["M"][index] and cells["R"][index]
                ),
                "M and R met after address restriction",
            )
            if cells["M"][index] or cells["R"][index]:
                chosen = cells["M"][index] or cells["R"][index]
                positive_stats.add(cell_stats(chosen))
                require(
                    cells["B"][index] == chosen,
                    "B stopped decomposing as the selected M/R cell",
                )
            else:
                require(
                    not cells["B"][index],
                    "B survived where both M and R vanished",
                )
        records[sigma, tau] = {
            "B": tuple(cells["B"]),
            "M": tuple(cells["M"]),
            "R": tuple(cells["R"]),
            "bits": (b_bits, m_bits, r_bits),
        }


pattern_counts = Counter(
    {pattern: len(labels) for pattern, labels in pattern_labels.items()}
)
expected_patterns = Counter({
    ((False, True, False), (True, False, True)): 63,
    ((False, False, False), (False, False, True)): 9,
    ((False, False, False), (True, False, True)): 9,
})
require(
    pattern_counts == expected_patterns,
    "the complete 81-label allocation-pattern census changed",
)

horn_pattern = ((False, True, False), (True, False, True))
horn_labels = tuple(pattern_labels[horn_pattern])
expected_horn_labels = tuple(
    (sigma, tau)
    for sigma in EXPECTED_HORN_SIGMAS
    for tau in TAUS
)
require(
    horn_labels == expected_horn_labels
    and tuple(pattern_labels[((False, False, False), (False, False, True))])
    == tuple((1, tau) for tau in TAUS)
    and tuple(pattern_labels[((False, False, False), (True, False, True))])
    == tuple((2, tau) for tau in TAUS),
    "the horn/nonhorn label factorization changed",
)
require(
    cross_implementation_checks == 81 * 3 * 2,
    "not every M/R address cell was cross-checked",
)
require(
    positive_stats
    == {(1, EXPECTED_MASS, EXPECTED_COEFFICIENT)},
    "positive cells stopped being one common weighted cylinder",
)


gaps = (N_PLUS - N0, N_A - N_PLUS, N_A - N0)
translation_fractions = tuple(
    Fraction(gap * direct.SHIFT, direct.PERIOD) % 1
    for gap in gaps
)
require(
    translation_fractions
    == (
        Fraction(7, 371_293),
        Fraction(28_553, 28_561),
        Fraction(371_196, 371_293),
    ),
    "physical translation fractions changed",
)
translation_checks = 0
for label in horn_labels:
    row = records[label]
    r0, _r1, ra = row["R"]
    _m0, mplus, _ma = row["M"]
    require(
        direct.shift_weighted(r0, gaps[0] * direct.SHIFT) == mplus,
        f"R(n0) stopped translating to M(n+) at {label}",
    )
    require(
        direct.shift_weighted(mplus, gaps[1] * direct.SHIFT) == ra,
        f"M(n+) stopped translating to R(na) at {label}",
    )
    require(
        direct.shift_weighted(r0, gaps[2] * direct.SHIFT) == ra,
        f"R(n0) stopped translating directly to R(na) at {label}",
    )
    translation_checks += 3
require(
    translation_checks == 3 * len(horn_labels),
    "not every horn translation was checked",
)


# Rebuild THM-2791's literal contributor walls.  A single open chamber
# containing all three cylinders proves identity of every (a,b,e') label set;
# no equinumerous but nonidentity correspondence is being substituted.
e_intervals = tuple(sheet.base.build_set(sheet.base.PAT_E, sheet.base.ZELL))
q_intervals = tuple(sheet.base.build_set(sheet.host.PAT_QB, sheet.base.ZELL))
q_depth, q_starts = sheet.scaled_intervals(q_intervals, sheet.DEPTH)
e_depth_pack, e_depth_pack_starts = sheet.scaled_intervals(
    e_intervals, sheet.DEPTH * sheet.PACK
)
e_depth, e_depth_starts = sheet.scaled_intervals(
    e_intervals, sheet.DEPTH
)

cylinders = tuple(target_cylinder(index)[0] for index in CENTRAL_INDICES)
hull_left = min(left for left, _right in cylinders)
hull_right = max(right for _left, right in cylinders)
u_events = tuple(sorted(set(
    sheet.mapped_events(q_intervals, sheet.DEPTH)
    + sheet.mapped_events(e_intervals, sheet.DEPTH * sheet.PACK)
)))
v_events = sheet.mapped_events(
    e_intervals,
    sheet.DEPTH,
    sheet.RAIL_DISPLACEMENT * sheet.T // P,
)
u_chamber = sheet.chamber(u_events, hull_left, hull_right)
v_chamber = sheet.chamber(v_events, hull_left, hull_right)
common_chamber = (
    max(u_chamber[0], v_chamber[0]),
    min(u_chamber[1], v_chamber[1]),
)
require(
    common_chamber[0] < hull_left
    and hull_right < common_chamber[1],
    "one of the horn cylinders meets a literal ancestry wall",
)

probes = tuple(ceil_fraction(left) for left, _right in cylinders)
require(
    all(
        left < probe < right
        for probe, (left, right) in zip(probes, cylinders)
    ),
    "an ancestry probe missed its open cylinder",
)

arguments = (
    q_depth,
    q_starts,
    e_depth_pack,
    e_depth_pack_starts,
    e_depth,
    e_depth_starts,
)
middle_u, middle_v = sheet.contributor_sets(probes[1], *arguments)
middle_digest = sheet.path_digest(middle_u, middle_v)
require(
    len(middle_u) == 966_606
    and len(middle_v) == 28_534
    and len(middle_u) * len(middle_v) == EXPECTED_WEIGHT
    and middle_digest == EXPECTED_PATH_DIGEST,
    "the newly enumerated middle ancestry set differs from THM-2791",
)


def supplied_path_active(coordinate):
    a_label, b_label, e_label = sheet.SUPPLIED_PATH
    y_numerator = coordinate + a_label * sheet.T
    z_numerator = y_numerator + b_label * sheet.DEPTH * sheet.T
    rotated = (
        coordinate
        - sheet.RAIL_DISPLACEMENT * sheet.T // P
    ) % sheet.T
    e_numerator = rotated + e_label * sheet.T
    return (
        sheet.contains_scaled(q_depth, q_starts, y_numerator)
        and sheet.contains_scaled(
            e_depth_pack, e_depth_pack_starts, z_numerator
        )
        and sheet.contains_scaled(e_depth, e_depth_starts, e_numerator)
    )


require(
    all(supplied_path_active(probe) for probe in probes),
    "the supplied literal Boolean path is absent at a horn vertex",
)
collision_labels = tuple(
    (
        (P * probe // sheet.T - sheet.RAIL_DISPLACEMENT) % P,
        P * probe // sheet.T,
        (2 * P * probe // sheet.T) % P,
    )
    for probe in probes
)
require(
    collision_labels == ((5, 6, 12),) * 3,
    "collision/root/deep labels changed along the horn",
)


source_tree = ast.parse(Path(__file__).read_text())
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


def bits(pattern):
    return "".join(str(int(value)) for value in pattern)


print("LRC14 SAME-ANCESTRY ALLOCATION-SWITCHING HORN SCOUT")
print("status=FINITE-EXACT candidate;not-yet-canonical")
print("dependency_hashes=BEGIN")
for path, digest in ACTUAL_HASHES:
    print(f"{digest}  {path}")
print("dependency_hashes=END")
print(
    "universe=rail8;physical_clock=1;"
    f"sigma={SIGMAS};tau={TAUS};"
    f"addresses={ADDRESSES};central_indices={CENTRAL_INDICES};"
    "both_relative_safeties;target_root1;delayed_carry6"
)
print(
    "typing=A=source_one_sided;B=pulled_target_one_sided;"
    "M=A_intersect_B;R=B_minus_A;"
    "M_and_R_are_pushed_to_the_common_target_coordinate_before_cylinder_cut"
)
print(
    "implementation_A=full_weighted_one_sided_carriers_plus_direct_set_operations;"
    "implementation_B=canonical_carrier_objects_plus_inherited_filter_order;"
    f"exact_MR_cell_agreement={cross_implementation_checks}/486"
)
print(
    "complete_81_label_pattern_histogram="
    f"M010_R101:{pattern_counts[horn_pattern]},"
    f"M000_R001:{pattern_counts[((False, False, False), (False, False, True))]},"
    f"M000_R101:{pattern_counts[((False, False, False), (True, False, True))]}"
)
print(
    f"horn_labels={len(horn_labels)}=clock1_x_sigma{EXPECTED_HORN_SIGMAS}_x_tau{TAUS};"
    "allocation_word=R-M-R;B_word=1-1-1"
)
print(
    "nonhorn_boundary=sigma1:M000_R001;sigma2:M000_R101;"
    f"pattern_digest={pattern_digest.hexdigest()}"
)
print(
    f"positive_whole_cylinder=(pieces:1,weight:{EXPECTED_WEIGHT},"
    f"mass:{EXPECTED_MASS},coefficient:{EXPECTED_COEFFICIENT})"
)
print(
    f"translation_fractions=(R_n0_to_M_nplus:{translation_fractions[0]},"
    f"M_nplus_to_R_na:{translation_fractions[1]},"
    f"R_n0_to_R_na:{translation_fractions[2]});"
    f"exact_translation_checks={translation_checks}/189"
)
print(
    f"ancestry_wall_counts=(U:{len(u_events)},V:{len(v_events)});"
    f"common_chamber={common_chamber};probes={probes}"
)
print(
    f"middle_literal_factor_counts=({len(middle_u)},{len(middle_v)});"
    f"product={len(middle_u) * len(middle_v)};"
    f"path_set_digest={middle_digest};"
    "same_chamber_identity_for_all_three_vertices=yes"
)
print(
    f"supplied_path={sheet.SUPPLIED_PATH};active_at_all_three=yes;"
    f"collision_labels={collision_labels}"
)
print(
    "MECHANISM=the_full_target_carrier_B_is_one_positive_three_vertex_"
    "address_horn_while_its_disjoint_allocation_decomposition_switches_"
    "right_cofiber_to_common_carrier_to_right_cofiber_on_one_literal_rail_"
    "ancestry_chamber"
)
print(
    "CONNECTION=the_THM2791_tau3_endpoint_chord_acquires_a_genuine_common_"
    "middle_vertex_inside_the_THM2806_common_label_bank;the_tau12_"
    "target_only_hostile_remains_outside_this_bank"
)
print(
    "LOSS=no_endpoint_origin_or_endpoint_current;no_four_state_K4_face;"
    "no_factorwise_translation_covariance;no_global_packet_action;"
    "no_row_exclusion;no_LRC14"
)
print(
    f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})"
)
print("ALL EXACT CHECKS PASSED")
