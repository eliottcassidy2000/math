#!/usr/bin/env python3
"""Exact partial endpoint-origin test for the THM-3285 canonical horn.

The tested universe is the single label

    (physical clock, sigma, tau) = (1, 0, 3)

at the three THM-3285 addresses.  The script independently reconstructs the
R--M--R whole-cylinder word, then chooses the least integral unit interval
strictly inside each rational address cylinder.  These three unit selectors
translate with the whole cylinders and are the cheapest inputs accepted by
the inherited integer endpoint-current constructor.

For every selector the script keeps separate:

* the raw source and target one-sided carrier masks;
* the six typed semantic factors;
* all 169 endpoint origins in both certified exact-order fields; and
* the four fixed-sheet bare/source/target/both central coefficients.

It also exhausts the natural affine endpoint covariance class

    F'(x) = c chi_u(x) F(x+s)

on the three target endpoint banks.  This does not invent a map from the
address horn to THM-2772's endpoint-origin pullback.  It tests one explicit
origin and a canonical subatom, and records the first coordinate that blocks
extension to the whole horn.

The result is a finite-exact partial obstruction, not a theorem about every
possible endpoint coupling, a row exclusion, or LRC(14).
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
    "04-computation/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py":
        "c42d66498f460f2142ea375fe9d4047b82c62b872b35d5a1634d2bb4c80a68ee",
    "05-knowledge/results/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.out":
        "e89dce3307e5d374e8583f92e1b2da1214e44929e52fdd42c6532d61adb3e246",
    "04-computation/lrc14_carrier_allocation_pullback_k4_segre_thm2772.py":
        "10b681f575fb51eb4b1af4bc909fba89846b85bd5da36fc069dac97ae2ebe409",
    "05-knowledge/results/lrc14_carrier_allocation_pullback_k4_segre_thm2772.out":
        "9cd1e82634822e4997d9966a8252b566d099050e42593a47c49b8a0f387190e5",
    "04-computation/lrc14_positive_graded_address_two_simplex_thm2807.py":
        "11cdbe3c6cc7f9d5b6b24863ced71eb91cc84adc67fe38a3f8a3e637362453fb",
    "05-knowledge/results/lrc14_positive_graded_address_two_simplex_thm2807.out":
        "a6fbc42c5a9fa7b84fa42a6fb625228851385777155e22653e360e190ea765d9",
    "04-computation/lrc14_affine_lift_transvection_horn_decoder_thm2813.py":
        "255ce911a18f33d6b700213d6362886970b12170c324d39bed82a69a51b63e83",
    "05-knowledge/results/lrc14_affine_lift_transvection_horn_decoder_thm2813.out":
        "32f61740a0e7e73384c3d1ff1a83ba30d4cd3ebeca5d228b57e0bf2510928d58",
    "04-computation/lrc14_boolean_norm_cotangent_boundary_thm2820.py":
        "8f9a51e0fd616cd616eef0080bcd054b2a0c191704e1f289e78ea28364476376",
    "05-knowledge/results/lrc14_boolean_norm_cotangent_boundary_thm2820.out":
        "d9f5cf1e0999f414e0c9af838070441def6baf1c65df266fbdb44b9f578cce58",
    "04-computation/lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
    "05-knowledge/results/lrc14_literal_fixed_sheet_allocation_thm2806.out":
        "a970776ed95128b5745c1fd370af768778b409d931a57b68006e04271e00813f",
    "04-computation/lrc14_opposite_wing_common_atom_virtual_sector_audit_thm2782.py":
        "f696867c5242e86593de64f3ee7944bf44a3887e7818291ed059c7358f959c1a",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


ACTUAL_HASHES = tuple(
    (path, lf_hash(ROOT / path)) for path in PINNED
)
require(
    all(digest == PINNED[path] for path, digest in ACTUAL_HASHES),
    "a pinned inheritance artifact changed",
)


import lrc14_literal_fixed_sheet_allocation_thm2806 as allocation
import lrc14_opposite_wing_common_atom_virtual_sector_audit_thm2782 as physical


def fractional_part(value):
    return value - value.numerator // value.denominator


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def target_cylinder(central_index):
    center = fractional_part(
        physical.physical.relative.Z
        + (N0 + P * central_index + 1)
        * Fraction(7, ADDRESS_MODULUS)
    )
    coordinate = center * physical.PERIOD
    half_width = physical.physical.relative.Q_RADIUS * physical.PERIOD
    require(
        0 < coordinate - half_width < coordinate + half_width < physical.PERIOD,
        "a selected target cylinder crossed the circle seam",
    )
    return ((coordinate - half_width, coordinate + half_width),)


def weighted_stats(cell, q_pairs):
    coefficient = physical.physical.relative.private.delayed_carry_pair(
        cell, q_pairs[CLOCK], {}
    )[6][1]
    return len(cell), physical.weighted_mass(cell), coefficient


def carrier_mask(interval, carrier, sign, unit):
    mask = []
    for harmonic in range(P):
        translated = allocation.physical.overlap.shift_weighted(
            carrier, sign * harmonic * unit
        )
        overlap = allocation.intersect_sorted(
            (interval,), allocation.support_of(translated)
        )
        require(
            not overlap or overlap == (interval,),
            "a carrier twist cut through the canonical unit atom",
        )
        mask.append(bool(overlap))
    return tuple(mask)


def inverse_endpoint_banks(source_atom, target_atom, endpoint_context):
    endpoint_base, present_sets, terminal, q_starts, tabs = endpoint_context
    mods = endpoint_base.endpoint.MODS
    left_cache = {}
    right_cache = {}
    left_dual = {}
    right_dual = {}

    for address, ell in endpoint_base.REPS.items():
        present = present_sets[address]
        starts = tuple(left for left, _right in present)
        source_restricted = allocation.indexed_weighted_intersection(
            source_atom, present, starts
        )
        target_restricted = allocation.indexed_weighted_intersection(
            target_atom, present, starts
        )

        if source_restricted not in left_cache:
            values = []
            for field_index, embedding in enumerate(mods):
                value = allocation.x_sweep(
                    source_restricted,
                    terminal,
                    q_starts,
                    embedding,
                    tabs[field_index],
                )[0]
                values.append(value)
            left_cache[source_restricted] = tuple(values)
        if target_restricted not in right_cache:
            right_cache[target_restricted] = tuple(
                allocation.endpoint_sum(
                    target_restricted, -endpoint_base.Y0, embedding
                )
                for embedding in mods
            )

        left_values = []
        for field_index, (prime, root) in enumerate(mods):
            zeta = pow(root, endpoint_base.NRED // P, prime)
            phase = pow(zeta, ell[endpoint_base.DEEP], prime)
            left_values.append(
                left_cache[source_restricted][field_index] * phase % prime
            )
        left_dual[address] = tuple(left_values)
        right_dual[address] = right_cache[target_restricted]

    fields = []
    for field_index, (prime, root) in enumerate(mods):
        zeta = pow(root, endpoint_base.NRED // P, prime)
        powers = tuple(pow(zeta, exponent, prime) for exponent in range(P))
        left = {}
        right = {}
        for point in endpoint_base.KEYS:
            left[point] = sum(
                left_dual[address][field_index]
                * powers[
                    -(
                        address[0] * point[0]
                        + address[1] * point[1]
                    ) % P
                ]
                for address in endpoint_base.KEYS
            ) % prime
            right[point] = sum(
                right_dual[address][field_index]
                * powers[
                    (
                        address[0] * point[0]
                        + address[1] * point[1]
                    ) % P
                ]
                for address in endpoint_base.KEYS
            ) % prime
        fields.append((prime, zeta, left, right))
    return tuple(fields)


def point_add(left, right):
    return (
        (left[0] + right[0]) % P,
        (left[1] + right[1]) % P,
    )


def affine_covariances(bare, carried, zeta, prime, points):
    """Exhaust F'(x)=c*zeta^(u.x)*F(x+s) for full-support banks."""
    require(
        all(bare.values()) and all(carried.values()),
        "affine covariance shortcut requires full support",
    )
    logarithm = {
        pow(zeta, exponent, prime): exponent for exponent in range(P)
    }
    zero = (0, 0)
    basis = ((1, 0), (0, 1))
    answers = []
    for shift in points:
        scalar = (
            carried[zero] * pow(bare[shift], -1, prime)
        ) % prime
        character = []
        valid = True
        for vector in basis:
            ratio = (
                carried[vector]
                * pow(bare[point_add(vector, shift)], -1, prime)
                * pow(scalar, -1, prime)
            ) % prime
            if ratio not in logarithm:
                valid = False
                break
            character.append(logarithm[ratio])
        if not valid:
            continue
        character = tuple(character)
        for point in points:
            phase = pow(
                zeta,
                (
                    character[0] * point[0]
                    + character[1] * point[1]
                ) % P,
                prime,
            )
            predicted = (
                scalar * phase * bare[point_add(point, shift)]
            ) % prime
            if predicted != carried[point]:
                valid = False
                break
        if valid:
            answers.append((shift, character, scalar))
    return tuple(answers)


def hadamard(products, prime):
    bare, source, target, both = products
    return (
        sum(products) % prime,
        (bare + source - target - both) % prime,
        (bare - source + target - both) % prime,
        (bare - source - target + both) % prime,
    )


def main():
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

    cells = {"B": [], "M": [], "R": []}
    cylinders = []
    for central_index in CENTRAL_INDICES:
        cylinder = target_cylinder(central_index)
        cylinders.append(cylinder)
        cells["B"].append(physical.intersect_weighted(target_one, cylinder))
        cells["M"].append(physical.intersect_weighted(common, cylinder))
        cells["R"].append(physical.intersect_weighted(right, cylinder))
    cells = {key: tuple(value) for key, value in cells.items()}
    bits = {
        key: tuple(bool(cell) for cell in value)
        for key, value in cells.items()
    }
    require(
        bits == {
            "B": (True, True, True),
            "M": (False, True, False),
            "R": (True, False, True),
        },
        "canonical horn stopped having the R-M-R allocation word",
    )
    positive_cells = tuple(
        cells["M"][index] or cells["R"][index]
        for index in range(3)
    )
    require(
        all(
            weighted_stats(cell, q_pairs)
            == (1, EXPECTED_MASS, EXPECTED_COEFFICIENT)
            and cell[0][2] == EXPECTED_WEIGHT
            for cell in positive_cells
        ),
        "a horn cell stopped being the expected whole weighted cylinder",
    )

    gaps = (N_PLUS - N0, N_A - N_PLUS, N_A - N0)
    require(
        tuple(address % P**2 for address in ADDRESSES) == (85, 98, 98)
        and N_A - N_PLUS == P**2 * 4_079,
        "the quotient horn or full-depth vertical edge changed",
    )
    require(
        physical.shift_weighted(
            positive_cells[0], gaps[0] * physical.SHIFT
        ) == positive_cells[1]
        and physical.shift_weighted(
            positive_cells[1], gaps[1] * physical.SHIFT
        ) == positive_cells[2]
        and physical.shift_weighted(
            positive_cells[0], gaps[2] * physical.SHIFT
        ) == positive_cells[2],
        "whole-cylinder horn translation changed",
    )

    target_atoms = []
    source_atoms = []
    for cell in positive_cells:
        left, right_endpoint, _weight = cell[0]
        integer_left = ceil_fraction(left)
        require(
            left < integer_left < integer_left + 1 < right_endpoint,
            "a horn cylinder lost its canonical internal unit atom",
        )
        target_atom = ((integer_left, integer_left + 1, 1),)
        source_atom = physical.shift_weighted(
            target_atom, -physical.SHIFT
        )
        target_atoms.append(target_atom)
        source_atoms.append(source_atom)
    target_atoms = tuple(target_atoms)
    source_atoms = tuple(source_atoms)
    require(
        target_atoms == (
            ((142_082_432_180_573, 142_082_432_180_574, 1),),
            ((142_088_047_310_093, 142_088_047_310_094, 1),),
            ((142_004_622_528_653, 142_004_622_528_654, 1),),
        )
        and source_atoms == (
            ((142_082_000_247_533, 142_082_000_247_534, 1),),
            ((142_087_615_377_053, 142_087_615_377_054, 1),),
            ((142_004_190_595_613, 142_004_190_595_614, 1),),
        ),
        "canonical integer endpoint atoms changed",
    )
    require(
        physical.shift_weighted(
            target_atoms[0], gaps[0] * physical.SHIFT
        ) == target_atoms[1]
        and physical.shift_weighted(
            target_atoms[1], gaps[1] * physical.SHIFT
        ) == target_atoms[2]
        and physical.shift_weighted(
            target_atoms[0], gaps[2] * physical.SHIFT
        ) == target_atoms[2],
        "canonical target endpoint atoms stopped following the horn",
    )

    factor_names = ("E3", "clock", "q1", "q2", "c2", "c3")
    factors = allocation.section_factors(
        full_module, e3, clocks, SIGMA, TAU, CLOCK
    )
    source_signatures = tuple(
        tuple(
            allocation.contained(atom[0][:2], factor)
            for factor in factors
        )
        for atom in source_atoms
    )
    target_signatures = tuple(
        tuple(
            allocation.contained(atom[0][:2], factor)
            for factor in factors
        )
        for atom in target_atoms
    )
    require(
        source_signatures == (
            (False, True, True, True, False, True),
            (True, True, True, True, True, True),
            (False, True, True, True, True, True),
        )
        and target_signatures == ((True,) * 6,) * 3,
        "the exact six-factor horn signatures changed",
    )

    unit = full_module.T // P
    raw_carrier_masks = tuple(
        (
            carrier_mask(atom_s[0][:2], source_base, -1, unit),
            carrier_mask(atom_t[0][:2], target_base, +1, unit),
        )
        for atom_s, atom_t in zip(source_atoms, target_atoms)
    )
    delta_zero = (True,) + (False,) * (P - 1)
    require(
        raw_carrier_masks == ((delta_zero, delta_zero),) * 3,
        "a raw carrier mask stopped being the fixed-sheet delta at zero",
    )

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
    endpoint_context = (
        endpoint_base, present_sets, terminal, q_starts, tabs
    )
    endpoint_fields = tuple(
        inverse_endpoint_banks(source_atom, target_atom, endpoint_context)
        for source_atom, target_atom in zip(source_atoms, target_atoms)
    )

    expected_origin_values = {
        352_341_050_142_921_841: (
            (93_139_741_941_715_911, 55_846_846_208_753_388),
            (217_164_758_098_570_837, 9_366_449_861_706_094),
            (0, 199_849_486_328_444_172),
        ),
        956_354_278_959_359_281: (
            (752_032_664_919_432_295, 889_340_730_416_083_895),
            (240_770_923_428_182_654, 330_542_147_584_671_122),
            (0, 444_326_171_470_978_060),
        ),
    }
    support_rows = []
    origin_rows = []
    typed_k4_rows = []
    for vertex_index, address in enumerate(ADDRESSES):
        source_gate = all(source_signatures[vertex_index])
        target_gate = all(target_signatures[vertex_index])
        for prime, _zeta, left, right_bank in endpoint_fields[vertex_index]:
            support_row = (
                address,
                prime,
                sum(bool(value) for value in left.values()),
                sum(bool(value) for value in right_bank.values()),
            )
            support_rows.append(support_row)
            origin_values = (left[(0, 0)], right_bank[(0, 0)])
            require(
                origin_values == expected_origin_values[prime][vertex_index],
                "an explicit-origin endpoint value changed",
            )
            origin_rows.append((address, prime, origin_values))

            source_present = left[(0, 0)] if source_gate else 0
            target_present = right_bank[(0, 0)] if target_gate else 0
            source_bare = P * source_present % prime
            target_bare = P * target_present % prime
            factors_at_origin = (
                source_bare,
                source_present,
                target_bare,
                target_present,
            )
            products = (
                source_bare * target_bare % prime,
                source_present * target_bare % prime,
                source_bare * target_present % prime,
                source_present * target_present % prime,
            )
            transform = hadamard(products, prime)
            typed_k4_rows.append(
                (address, prime, factors_at_origin, products, transform)
            )
            if vertex_index == 1:
                require(
                    all(factors_at_origin)
                    and all(products)
                    and all(transform)
                    and transform[3]
                    == (P - 1) ** 2 * products[3] % prime,
                    "the middle central K4 certificate lost nonvanishing",
                )
                require(
                    products[3] - products[3]
                    - products[3] + products[3] == 0,
                    "the raw common-twist K4 stopped being flat",
                )
            else:
                require(
                    factors_at_origin[:2] == (0, 0)
                    and products == (0, 0, 0, 0),
                    "an outer typed source fibre unexpectedly appeared",
                )

    require(
        tuple(
            (row[0], row[2], row[3]) for row in support_rows[::2]
        )
        == (
            (N0, 169, 169),
            (N_PLUS, 169, 169),
            (N_A, 0, 169),
        )
        and tuple(
            (row[0], row[2], row[3]) for row in support_rows[1::2]
        )
        == (
            (N0, 169, 169),
            (N_PLUS, 169, 169),
            (N_A, 0, 169),
        ),
        "the two-field endpoint support census changed",
    )

    covariance_rows = []
    expected_positive_covariances = {
        352_341_050_142_921_841:
            (((0, 0), (0, 0), 85_672_738_390_772_541),),
        956_354_278_959_359_281:
            (((0, 0), (0, 0), 503_393_058_368_290_190),),
    }
    for field_index, (prime, zeta, _left0, right0) in enumerate(
        endpoint_fields[0]
    ):
        _prime_plus, _zeta_plus, _left_plus, right_plus = (
            endpoint_fields[1][field_index]
        )
        _prime_a, _zeta_a, _left_a, right_a = endpoint_fields[2][field_index]
        covariance_0_plus = affine_covariances(
            right0, right_plus, zeta, prime, endpoint_base.KEYS
        )
        covariance_plus_a = affine_covariances(
            right_plus, right_a, zeta, prime, endpoint_base.KEYS
        )
        covariance_0_a = affine_covariances(
            right0, right_a, zeta, prime, endpoint_base.KEYS
        )
        require(
            covariance_0_plus == expected_positive_covariances[prime]
            and not covariance_plus_a
            and not covariance_0_a,
            "the target endpoint affine-covariance boundary changed",
        )
        covariance_rows.append(
            (
                prime,
                covariance_0_plus,
                covariance_plus_a,
                covariance_0_a,
            )
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

    print("LRC14 ENDPOINT-ORIGIN K4 HORN PARTIAL SCOUT")
    print("status=FINITE-EXACT PARTIAL;not_a_canonical_theorem")
    print("dependency_hashes=BEGIN")
    for path, digest in ACTUAL_HASHES:
        print(f"{digest}  {path}")
    print("dependency_hashes=END")
    print(
        "universe=rail8;(physical_clock,sigma,tau)=(1,0,3);"
        f"addresses={ADDRESSES};central_indices={CENTRAL_INDICES};"
        "both_relative_safeties;target_root1;delayed_carry6;"
        "endpoint_origins=F13^2=169;certified_fields=2"
    )
    print(
        "inheritance=THM2772_K4_requires_one_common_endpoint_atom;"
        "THM2807_supplies_the_address_triangle;"
        "THM2813_says_the_residue7_horn_is_fixed_and_needs_a_normal_jet;"
        "THM2820_separates_raw_common_support_from_target_only_endpoint_motion"
    )
    print(
        f"whole_cylinder_word=R-M-R;B=1-1-1;weight={EXPECTED_WEIGHT};"
        f"mass={EXPECTED_MASS};coefficient={EXPECTED_COEFFICIENT};"
        "whole_cylinder_translations=3/3"
    )
    print(
        f"canonical_target_unit_atoms={target_atoms};"
        f"canonical_source_unit_atoms={source_atoms};unit_translations=3/3"
    )
    print(
        "raw_carrier_twist_supports="
        f"{tuple((tuple(i for i, bit in enumerate(source) if bit), tuple(i for i, bit in enumerate(target) if bit)) for source, target in raw_carrier_masks)}"
    )
    print(
        "six_factor_names="
        f"{factor_names};source_signatures={source_signatures};"
        f"target_signatures={target_signatures}"
    )
    print(
        "first_typed_failure=n0_loses_(E3,c2);"
        "na_loses_E3_only;raw_source_carrier_is_present_at_all_three"
    )
    print(f"endpoint_support_rows_(address,field,left,right)={tuple(support_rows)}")
    print(f"explicit_origin00_rows_(address,field,(left,right))={tuple(origin_rows)}")
    print(
        "typed_origin00_K4_rows_(address,field,"
        f"(Pbare,Ppresent,Qbare,Qpresent),(B,S,T,H),Hadamard)={tuple(typed_k4_rows)}"
    )
    print(
        "middle_endpoint_fibre=all_169_left_and_right_origins_nonzero_in_both_fields;"
        "origin00_has_four_nonzero_central_allocation_states;"
        "raw_twist00_vector_is_flat_and_has_Mobius_face_zero"
    )
    print(
        "outer_typed_current_fibres=empty;"
        "na_underlying_left_endpoint_bank=0_of_169_before_semantic_gating;"
        "target_R_cylinders_remain_positive"
    )
    print(
        "target_endpoint_affine_covariances_"
        f"(field,n0_to_nplus,nplus_to_na,n0_to_na)={tuple(covariance_rows)}"
    )
    print(
        "full_depth_edge=nplus_to_na=Z2^4079;"
        "quotient_mod169_collapses_nplus_and_na;"
        "target_endpoint_affine_companion_on_vertical_edge=none_in_both_fields"
    )
    print(
        "FIRST_MISSING_COORDINATE=a_typed_map_from_the_address_and_outer_rail_"
        "ancestry_sheet_to_one_common_(L,R,q,Delta)_endpoint_atom_that_"
        "transports_E3_and_the_endpoint_present_selector;R_is_not_a_source_"
        "allocation_absence_bit"
    )
    print(
        "VERDICT=OUTCOME_2:middle_endpoint_fibre_present_but_outer_typed_"
        "endpoint_companions_and_the_Z2^4079_affine_lift_fail"
    )
    print(
        "scope=canonical_unit_subatom_inside_each_THM3285_whole_cylinder;"
        "fixed_origin00_plus_complete_169_origin_support_and_affine_covariance_"
        "census;no_claim_against_arbitrary_new_couplings;no_root_Cech_map;"
        "no_row_exclusion;no_LRC14"
    )
    print(
        f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
