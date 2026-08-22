#!/usr/bin/env python3
"""Independent hostile audit of the THM-3285 allocation-switching horn.

This audit does not import the candidate companion.  It reconstructs each
clock-one labelled one-sided carrier from the independently audited THM-2782
geometry, pushes the source carrier into the target chart, cuts by the address
cylinder first, and only then decomposes the target cell into its intersection
with the pushed source and its right cofiber.  This reverses the candidate's
decompose-then-cut order.

Literal ancestry is checked by separately enumerating the complete U and V
contributor label sets at all three vertices, not only at the new middle one.
The audit deliberately computes no endpoint-origin or current observable.
"""

from __future__ import annotations

import ast
from bisect import bisect_right
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path
import sys


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

CANDIDATE_SCRIPT = (
    "04-computation/"
    "lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py"
)
CANDIDATE_OUTPUT = (
    "05-knowledge/results/"
    "lrc14_allocation_switching_same_ancestry_horn_scout_20260803.out"
)
DIRECT_SCRIPT = (
    "04-computation/"
    "lrc14_opposite_wing_common_atom_virtual_sector_audit_thm2782.py"
)
SHEET_SCRIPT = "04-computation/lrc14_full_arm_orbit_path_sheet_audit_thm2791.py"
SHEET_OUTPUT = (
    "05-knowledge/results/lrc14_full_arm_orbit_path_sheet_audit_thm2791.out"
)
TAU12_SCRIPT = (
    "04-computation/"
    "lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806.py"
)
TAU12_OUTPUT = (
    "05-knowledge/results/"
    "lrc14_tau12_simplex_allocation_support_no_go_addendum_thm2806.out"
)
PINS = {
    CANDIDATE_SCRIPT:
        "c42d66498f460f2142ea375fe9d4047b82c62b872b35d5a1634d2bb4c80a68ee",
    CANDIDATE_OUTPUT:
        "e89dce3307e5d374e8583f92e1b2da1214e44929e52fdd42c6532d61adb3e246",
    DIRECT_SCRIPT:
        "f696867c5242e86593de64f3ee7944bf44a3887e7818291ed059c7358f959c1a",
    SHEET_SCRIPT:
        "1e00b6711db0d878fca70047f5f1532518084dbf6928551cd28fe51283fde543",
    SHEET_OUTPUT:
        "2325aa126dc4a97af4ba4bde0348389bfdfe5720dec6f88792b9a06baa40afd3",
    TAU12_SCRIPT:
        "7fea18161046f8de35b2e6ef04c88a13485de61045bc363f89c5ebfec8f76480",
    TAU12_OUTPUT:
        "cba545e7beff4fe889c76ae681c47c806969ccd33aa79def54527709a6ffafc6",
}


def lf_hash(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(payload).hexdigest()


ACTUAL_HASHES = tuple((path, lf_hash(ROOT / path)) for path in PINS)
require(all(digest == PINS[path] for path, digest in ACTUAL_HASHES),
        "candidate or hostile-control dependency drift")

source_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(source_tree)),
        "optimization-sensitive assert found")
require(not any(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(source_tree)
), "floating literal found")


import lrc14_full_arm_orbit_path_sheet_audit_thm2791 as sheet
import lrc14_opposite_wing_common_atom_virtual_sector_audit_thm2782 as direct


P = 13
ADDRESS_MODULUS = P**6
N0 = 3_454_614
N_PLUS = 3_454_627
N_A = 4_143_978
ADDRESSES = (N0, N_PLUS, N_A)
CENTRAL_INDICES = (0, 1, 53_028)
SIGMAS = (0, 1, 2, 3, 8, 9, 10, 11, 12)
TAUS = (3, 4, 5, 6, 7, 8, 9, 10, 11)
HORN_SIGMAS = (0, 3, 8, 9, 10, 11, 12)
TARGET_STEP = Fraction(7, ADDRESS_MODULUS)
EXPECTED_WEIGHT = 27_581_135_604
EXPECTED_MASS = Fraction(60_781_651_775_958_960, 371_293)
EXPECTED_COEFFICIENT = 790_161_473_087_466_480
EXPECTED_PATH_DIGEST = (
    "15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd"
)

require(direct.COMMON_S == SIGMAS and direct.COMMON_T == TAUS,
        "the inherited 9-by-9 common-label bank changed")
require(ADDRESS_MODULUS == 4_826_809
        and N_PLUS - N0 == P
        and N_A - N0 == 689_364,
        "address arithmetic drift")


def frac(value):
    return value - value.numerator // value.denominator


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def target_cylinder(central_index):
    """Rebuild the selected target cylinder without the candidate helper."""
    address = N0 + P * central_index
    center = frac(
        direct.physical.relative.Z + (address + 1) * TARGET_STEP
    )
    coordinate = center * direct.PERIOD
    half = direct.physical.relative.Q_RADIUS * direct.PERIOD
    require(0 < coordinate - half < coordinate + half < direct.PERIOD,
            "audit cylinder crossed the circle seam")
    return ((coordinate - half, coordinate + half),)


def translate_weighted(pieces, shift):
    """Translate half-open weighted intervals modulo the period independently."""
    period = direct.PERIOD
    shift %= period
    translated = []
    for left, right, weight in pieces:
        length = right - left
        new_left = (left + shift) % period
        new_right = new_left + length
        if new_right <= period:
            translated.append((new_left, new_right, weight))
        else:
            translated.append((new_left, period, weight))
            translated.append((0, new_right - period, weight))
    return direct.normalize_weighted(translated)


def manual_present_section(module, e3, clock_combs, clock, sigma, tau):
    """Spell out every clock and graft factor to guard MISTAKE-313."""
    intervals = direct.physical.target.intersect_sorted(
        e3, clock_combs[clock]
    )
    intervals = module.subtract_comb(
        intervals, module.W[1], 182, -14 * sigma - 13, -14 * sigma + 13
    )
    intervals = module.subtract_comb(
        intervals, module.W[2], 182, -14 * tau - 13, -14 * tau + 13
    )
    intervals = module.subtract_comb(
        intervals, module.C2, 182, 14 * sigma - 13, 14 * sigma + 13
    )
    intervals = module.subtract_comb(
        intervals, module.C3, 182, 14 * tau - 13, 14 * tau + 13
    )
    return tuple(intervals)


(
    _module,
    _rails,
    _present,
    overlap_details,
    full_module,
    e3,
    clock_combs,
    q_pairs,
) = direct.build_physical_geometry()

source_base, target_base = overlap_details[1]
cylinders = tuple(target_cylinder(index) for index in CENTRAL_INDICES)

records = {}
pattern_labels = defaultdict(list)
pattern_digest = sha256()
positive_stats = Counter()
section_clock_checks = 0
decomposition_checks = 0
full_objects = {}

for sigma, tau in product(SIGMAS, TAUS):
    section = manual_present_section(
        full_module, e3, clock_combs, 1, sigma, tau
    )
    inherited_section = tuple(
        direct.physical.target.source_present_section(
            full_module, e3, 1, sigma, tau, clock_combs
        )
    )
    require(section == inherited_section,
            ("manual fixed-clock section differs", sigma, tau))
    section_clock_checks += 1

    source_one = direct.intersect_weighted(source_base, section)
    target_one = direct.intersect_weighted(target_base, section)
    source_in_target = translate_weighted(source_one, direct.SHIFT)
    full_m = direct.intersect_weighted(
        target_one, direct.support_of(source_in_target)
    )
    full_r = direct.intersect_weighted(
        target_one, direct.complement(direct.support_of(source_in_target))
    )
    require(direct.disjoint_union(full_m, full_r) == target_one,
            ("target-chart M/R decomposition failed", sigma, tau))
    full_objects[sigma, tau] = (full_m, full_r, section)

    cells = {"B": [], "M": [], "R": []}
    for cylinder in cylinders:
        # Address-first route: cut B first, then allocate relative to pushed A.
        b_cell = direct.intersect_weighted(target_one, cylinder)
        m_cell = direct.intersect_weighted(
            b_cell, direct.support_of(source_in_target)
        )
        r_cell = direct.intersect_weighted(
            b_cell, direct.complement(direct.support_of(source_in_target))
        )
        require(direct.disjoint_union(m_cell, r_cell) == b_cell,
                ("address-first decomposition failed", sigma, tau, cylinder))
        require(
            m_cell == direct.intersect_weighted(full_m, cylinder)
            and r_cell == direct.intersect_weighted(full_r, cylinder),
            ("cut/decomposition square failed", sigma, tau, cylinder),
        )
        decomposition_checks += 2
        cells["B"].append(b_cell)
        cells["M"].append(m_cell)
        cells["R"].append(r_cell)

        for allocated in (m_cell, r_cell):
            if not allocated:
                continue
            require(len(allocated) == 1, "positive cell split into pieces")
            left, right, weight = allocated[0]
            mass = direct.weighted_mass(allocated)
            coefficient = direct.physical.relative.private.delayed_carry_pair(
                allocated, q_pairs[1], {}
            )[6][1]
            require(mass == (right - left) * weight,
                    "weighted mass stopped being literal cylinder mass")
            require(coefficient == ADDRESS_MODULUS * mass,
                    "carry-six coefficient lost whole-cylinder identity")
            positive_stats[(len(allocated), weight, mass, coefficient)] += 1

    m_bits = tuple(bool(cell) for cell in cells["M"])
    r_bits = tuple(bool(cell) for cell in cells["R"])
    b_bits = tuple(bool(cell) for cell in cells["B"])
    pattern_labels[m_bits, r_bits].append((sigma, tau))
    pattern_digest.update(
        (
            f"{sigma},{tau}:"
            f"M{''.join(str(int(bit)) for bit in m_bits)}:"
            f"R{''.join(str(int(bit)) for bit in r_bits)};"
        ).encode("ascii")
    )
    records[sigma, tau] = {
        "B": tuple(cells["B"]),
        "M": tuple(cells["M"]),
        "R": tuple(cells["R"]),
        "bits": (b_bits, m_bits, r_bits),
    }

horn_pattern = ((False, True, False), (True, False, True))
sigma_one_pattern = ((False, False, False), (False, False, True))
sigma_two_pattern = ((False, False, False), (True, False, True))
pattern_counts = Counter(
    {pattern: len(labels) for pattern, labels in pattern_labels.items()}
)
require(pattern_counts == Counter({
    horn_pattern: 63,
    sigma_one_pattern: 9,
    sigma_two_pattern: 9,
}), "independent 81-label pattern census disagrees")

horn_labels = tuple(pattern_labels[horn_pattern])
require(horn_labels == tuple(product(HORN_SIGMAS, TAUS)),
        "63-label Cartesian factorization failed")
require(tuple(pattern_labels[sigma_one_pattern])
        == tuple((1, tau) for tau in TAUS),
        "sigma-one hostile row changed")
require(tuple(pattern_labels[sigma_two_pattern])
        == tuple((2, tau) for tau in TAUS),
        "sigma-two hostile row changed")
require(all(records[label]["bits"][0] == (True, True, True)
            for label in horn_labels),
        "undecomposed horn target carrier stopped being 111")
require(section_clock_checks == 81 and decomposition_checks == 486,
        "fixed-clock or address-first census incomplete")
require(positive_stats == Counter({
    (1, EXPECTED_WEIGHT, EXPECTED_MASS, EXPECTED_COEFFICIENT): 216,
}), "positive cells lost their common whole-cylinder statistic")
require(pattern_digest.hexdigest()
        == "9fd5090fb89493d07115ee753a0eadeba82d70d0e52663aa11d60c121f4f9d4f",
        "full pattern digest disagrees with candidate")


gaps = (N_PLUS - N0, N_A - N_PLUS, N_A - N0)
translation_fractions = tuple(
    Fraction(gap * direct.SHIFT, direct.PERIOD) % 1 for gap in gaps
)
require(translation_fractions == (
    Fraction(7, 371_293),
    Fraction(28_553, 28_561),
    Fraction(371_196, 371_293),
), "translation fractions changed")
require((translation_fractions[0] + translation_fractions[1]) % 1
        == translation_fractions[2],
        "horn translations stopped composing")

translation_checks = 0
for label in horn_labels:
    row = records[label]
    r0, _rplus, ra = row["R"]
    _m0, mplus, _ma = row["M"]
    require(translate_weighted(r0, gaps[0] * direct.SHIFT) == mplus,
            ("independent R0-to-Mplus translation failed", label))
    require(translate_weighted(mplus, gaps[1] * direct.SHIFT) == ra,
            ("independent Mplus-to-Ra translation failed", label))
    require(translate_weighted(r0, gaps[2] * direct.SHIFT) == ra,
            ("independent R0-to-Ra translation failed", label))
    translation_checks += 3
require(translation_checks == 189, "translation audit incomplete")

# The sigma-two hostile retains only the endpoint chord.  Sigma one retains
# only the far endpoint.  These are literal controls against filling the horn.
sigma_two_endpoint_checks = 0
for tau in TAUS:
    r0, _rplus, ra = records[2, tau]["R"]
    require(translate_weighted(r0, gaps[2] * direct.SHIFT) == ra,
            ("sigma-two endpoint hostile stopped translating", tau))
    sigma_two_endpoint_checks += 1

# A local cylinder translation is not a global allocation action.  At the
# canonical label, neither full allocated piece translates to the next one.
full_m_03, full_r_03, section_03 = full_objects[0, 3]
require(translate_weighted(full_r_03, gaps[0] * direct.SHIFT) != full_m_03,
        "unexpected global R-to-M allocation action appeared")
require(translate_weighted(full_r_03, gaps[2] * direct.SHIFT) != full_r_03,
        "unexpected global R-periodicity appeared")
require(direct.shift_intervals(section_03, gaps[0] * direct.SHIFT)
        != section_03,
        "unexpected global factor-section covariance appeared")


# Independently enumerate the literal U(a,b) and V(e') ancestry labels at all
# three vertices.  This directly tests identity of label sets, not only counts.
e_intervals = tuple(sheet.base.build_set(sheet.base.PAT_E, sheet.base.ZELL))
q_intervals = tuple(sheet.base.build_set(sheet.host.PAT_QB, sheet.base.ZELL))
q_depth, q_starts = sheet.scaled_intervals(q_intervals, sheet.DEPTH)
e_depth_pack, e_depth_pack_starts = sheet.scaled_intervals(
    e_intervals, sheet.DEPTH * sheet.PACK
)
e_depth, e_depth_starts = sheet.scaled_intervals(e_intervals, sheet.DEPTH)

probe_coordinates = tuple(ceil_fraction(cylinder[0][0]) for cylinder in cylinders)
require(all(cylinder[0][0] < probe < cylinder[0][1]
            for probe, cylinder in zip(probe_coordinates, cylinders)),
        "integral ancestry probe missed a cylinder")

contributor_arguments = (
    q_depth,
    q_starts,
    e_depth_pack,
    e_depth_pack_starts,
    e_depth,
    e_depth_starts,
)
reference_u = None
reference_v = None
ancestry_records = []
for probe in probe_coordinates:
    u_labels, v_labels = sheet.contributor_sets(
        probe, *contributor_arguments
    )
    current_digest = sheet.path_digest(u_labels, v_labels)
    ancestry_records.append(
        (len(u_labels), len(v_labels), current_digest)
    )
    if reference_u is None:
        reference_u, reference_v = u_labels, v_labels
    else:
        require(u_labels == reference_u and v_labels == reference_v,
                "literal ancestry label sets differ across horn vertices")

require(tuple(ancestry_records) == (
    (966_606, 28_534, EXPECTED_PATH_DIGEST),
    (966_606, 28_534, EXPECTED_PATH_DIGEST),
    (966_606, 28_534, EXPECTED_PATH_DIGEST),
), "three-vertex ancestry enumeration drift")
require(len(reference_u) * len(reference_v) == EXPECTED_WEIGHT,
        "literal ancestry factorization stopped giving the rail weight")
require(sheet.SUPPLIED_PATH[0] * sheet.PACK + sheet.SUPPLIED_PATH[1]
        in reference_u and sheet.SUPPLIED_PATH[2] in reference_v,
        "supplied Boolean ancestry path disappeared")


def mapped_events(intervals, multiplier, shift=0):
    return tuple(sorted(set(
        (multiplier * endpoint + shift) % sheet.T
        for interval in intervals
        for endpoint in interval
    )))


def bracketing_chamber(events, left, right):
    index = bisect_right(events, left)
    require(0 < index < len(events) and events[index] > right,
            "a raw contributor wall crosses the horn hull")
    return events[index - 1], events[index]


hull_left = min(cylinder[0][0] for cylinder in cylinders)
hull_right = max(cylinder[0][1] for cylinder in cylinders)
u_events = tuple(sorted(set(
    mapped_events(q_intervals, sheet.DEPTH)
    + mapped_events(e_intervals, sheet.DEPTH * sheet.PACK)
)))
v_events = mapped_events(
    e_intervals,
    sheet.DEPTH,
    sheet.RAIL_DISPLACEMENT * sheet.T // P,
)
u_chamber = bracketing_chamber(u_events, hull_left, hull_right)
v_chamber = bracketing_chamber(v_events, hull_left, hull_right)
common_chamber = (
    max(u_chamber[0], v_chamber[0]),
    min(u_chamber[1], v_chamber[1]),
)
require(common_chamber == (140_890_500_190_440, 144_190_879_112_280)
        and common_chamber[0] < hull_left < hull_right < common_chamber[1],
        "literal common no-wall chamber changed")
require(len(u_events) == 34 and len(v_events) == 14,
        "raw ancestry-wall census changed")

collision_labels = tuple(
    (
        (P * probe // sheet.T - sheet.RAIL_DISPLACEMENT) % P,
        P * probe // sheet.T,
        (2 * P * probe // sheet.T) % P,
    )
    for probe in probe_coordinates
)
require(collision_labels == ((5, 6, 12),) * 3,
        "collision/root/deep labels changed")

# Scope audit against MISTAKE-281/300/310/313.  The candidate source contains
# no endpoint/current computation, and this hostile path likewise stops at the
# literal rail labels.  The imported direct module may contain unrelated
# endpoint audit utilities; none is called or exported into the horn records.
candidate_tree = ast.parse((ROOT / CANDIDATE_SCRIPT).read_text(encoding="utf-8"))
candidate_import_roots = tuple(
    alias.name
    for node in ast.walk(candidate_tree)
    if isinstance(node, (ast.Import, ast.ImportFrom))
    for alias in node.names
)
require(not any("endpoint" in name or "current" in name
                for name in candidate_import_roots),
        "candidate directly imported an endpoint/current constructor")
require(set(records[(0, 3)]) == {"B", "M", "R", "bits"},
        "candidate carrier record silently gained endpoint/current data")
tau12_text = (ROOT / TAU12_OUTPUT).read_text(encoding="utf-8")
require("M_empty_before_address_restriction=yes" in tau12_text
        and "ALL EXACT CHECKS PASSED" in tau12_text,
        "tau-twelve target-only hostile transcript changed")


print("THM-3285 SAME-ANCESTRY ALLOCATION HORN INDEPENDENT HOSTILE AUDIT")
print("status=FINITE-EXACT independent audit;candidate_not_promoted_here")
print("dependency_hashes=BEGIN")
for path, digest in ACTUAL_HASHES:
    print(f"{digest}  {path}")
print("dependency_hashes=END")
print(
    "universe=clock1_x_sigma%s_x_tau%s=81;manual_clocked_sections=81/81"
    % (SIGMAS, TAUS)
)
print(
    "independent_route=push_source_to_target_chart_then_address_cut_then_M=B_intersect_A_and_R=B_minus_A;"
    f"cut_decomposition_checks={decomposition_checks}/486"
)
print(
    "pattern_histogram=M010_R101:63,M000_R001:9,M000_R101:9;"
    f"pattern_digest={pattern_digest.hexdigest()}"
)
print(
    f"horn_factorization=clock1_x_sigma{HORN_SIGMAS}_x_tau{TAUS}=63:PASS"
)
print("sigma1_hostile=M000_R001_on_all_9_tau:PASS")
print("sigma2_hostile=M000_R101_on_all_9_tau:PASS")
print(
    f"positive_cells=216;whole_cylinder=(weight:{EXPECTED_WEIGHT},"
    f"mass:{EXPECTED_MASS},coefficient:{EXPECTED_COEFFICIENT});"
    "coefficient=13^6*mass:PASS"
)
print(
    f"translation_fractions={translation_fractions};"
    f"horn_translation_checks={translation_checks}/189;"
    f"sigma2_endpoint_controls={sigma_two_endpoint_checks}/9"
)
print(
    "global_hostiles=full_R_to_M_step_fails;full_R_diagonal_periodicity_fails;"
    "full_factor_section_step_covariance_fails"
)
print(
    f"ancestry_probes={probe_coordinates};"
    f"three_literal_enumerations={tuple(ancestry_records)}"
)
print(
    f"ancestry_wall_counts=(U:{len(u_events)},V:{len(v_events)});"
    f"common_no_wall_chamber={common_chamber};"
    "literal_UV_identity_at_all_three=PASS"
)
print(
    f"supplied_path={sheet.SUPPLIED_PATH}:active;"
    f"collision_labels={collision_labels}"
)
print("tau12_external_hostile=common_carrier_empty_before_address:PASS")
print(
    "scope_audit=no_endpoint_origin_enumerated;no_current_observable_evaluated;"
    "no_K4_allocation;no_factorwise_or_global_action;no_LRC14"
)
print("mistake_controls=281_same_atom_only;300_truth_sets_not_vocabulary;310_physical_support_literal;313_clock_factor_explicit")
print("source_ast=(assert_nodes=0,float_literals=0)")
print("verdict=ACCEPT_THM3285_FINITE_RESPONSE_CARRIER_CLAIMS_WITH_STATED_SCOPE")
