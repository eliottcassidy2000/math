#!/usr/bin/env python3
"""Exact C13 reflection/Gabor hostile sidecar for the THM-3285 endpoint scout.

This script reconstructs the canonical ``(clock,sigma,tau)=(1,0,3)`` unit
subatoms and all three 169-origin endpoint banks through the pinned endpoint
scout.  It then keeps the address and six-factor semantic sidecars while
testing the owner-phase reflection

    r(a,b) = (a, 2*a-b)  on F_13^2.

Two exact tests are separated deliberately:

* the integer Hermitian support compression ``direct_sum_v D_v J_r D_v``;
* covariance of the actual middle endpoint amplitudes under every allowed
  finite Heisenberg gauge ``f(x) -> c*zeta^(u.x)*f(x+s)``.

The result is a hostile finite-exact negative control.  It is not a row
exclusion, a physical owner/current construction, or an LRC(14) theorem.

Reproduce from the repository root with:

    python 04-computation/lrc14_c13_reflection_gabor_hostile_sidecar_20260812.py
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

P = 13
SCOUT_PATH = (
    "04-computation/"
    "lrc14_endpoint_origin_k4_horn_partial_scout_20260803.py"
)
SCOUT_OUTPUT_PATH = (
    "05-knowledge/results/"
    "lrc14_endpoint_origin_k4_horn_partial_scout_20260803.out"
)
REPRODUCTION_COMMAND = (
    "python 04-computation/"
    "lrc14_c13_reflection_gabor_hostile_sidecar_20260812.py"
)
STORED_OUTPUT_PATH = (
    "05-knowledge/results/"
    "lrc14_c13_reflection_gabor_hostile_sidecar_20260812.out"
)
PINNED = {
    SCOUT_PATH:
        "5f1258cd19561941de93135b857bc4187435d10bfa0c1ca0be186c3d8e039e0c",
    SCOUT_OUTPUT_PATH:
        "af9c946d7832c8dfc6442846e7616fadb153cda33f4ef5386272671fe856cbe7",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


ACTUAL_HASHES = tuple((path, lf_hash(ROOT / path)) for path in PINNED)
require(
    all(digest == PINNED[path] for path, digest in ACTUAL_HASHES),
    "the canonical endpoint scout or its stored output changed",
)


import lrc14_endpoint_origin_k4_horn_partial_scout_20260803 as scout


EXPECTED_ATOMS = (
    (
        ((142_082_000_247_533, 142_082_000_247_534, 1),),
        ((142_087_615_377_053, 142_087_615_377_054, 1),),
        ((142_004_190_595_613, 142_004_190_595_614, 1),),
    ),
    (
        ((142_082_432_180_573, 142_082_432_180_574, 1),),
        ((142_088_047_310_093, 142_088_047_310_094, 1),),
        ((142_004_622_528_653, 142_004_622_528_654, 1),),
    ),
)
EXPECTED_SOURCE_SIGNATURES = (
    (False, True, True, True, False, True),
    (True, True, True, True, True, True),
    (False, True, True, True, True, True),
)
EXPECTED_TARGET_SIGNATURES = ((True,) * 6,) * 3
EXPECTED_SUPPORT_ROWS = (
    (3_454_614, 352_341_050_142_921_841, 169, 169),
    (3_454_614, 956_354_278_959_359_281, 169, 169),
    (3_454_627, 352_341_050_142_921_841, 169, 169),
    (3_454_627, 956_354_278_959_359_281, 169, 169),
    (4_143_978, 352_341_050_142_921_841, 0, 169),
    (4_143_978, 956_354_278_959_359_281, 0, 169),
)
EXPECTED_ORIGIN00_ROWS = (
    (3_454_614, 352_341_050_142_921_841,
     (93_139_741_941_715_911, 55_846_846_208_753_388)),
    (3_454_614, 956_354_278_959_359_281,
     (752_032_664_919_432_295, 889_340_730_416_083_895)),
    (3_454_627, 352_341_050_142_921_841,
     (217_164_758_098_570_837, 9_366_449_861_706_094)),
    (3_454_627, 956_354_278_959_359_281,
     (240_770_923_428_182_654, 330_542_147_584_671_122)),
    (4_143_978, 352_341_050_142_921_841,
     (0, 199_849_486_328_444_172)),
    (4_143_978, 956_354_278_959_359_281,
     (0, 444_326_171_470_978_060)),
)


def reconstruct_packet():
    """Rebuild canonical atoms, signatures, and endpoint banks from source."""
    (
        _module,
        _rails,
        _present,
        overlap_details,
        full_module,
        e3,
        clocks,
        _q_pairs,
    ) = scout.physical.build_physical_geometry()
    source_base, target_base = overlap_details[scout.CLOCK]
    section = tuple(
        scout.physical.physical.target.source_present_section(
            full_module,
            e3,
            scout.CLOCK,
            scout.SIGMA,
            scout.TAU,
            clocks,
        )
    )
    target_one = scout.physical.intersect_weighted(target_base, section)
    source_one = scout.physical.intersect_weighted(source_base, section)
    target_pullback = scout.physical.shift_weighted(
        target_one, -scout.physical.SHIFT
    )
    common_pullback = scout.physical.intersect_weighted(
        target_pullback, scout.physical.support_of(source_one)
    )
    right_pullback = scout.physical.intersect_weighted(
        target_pullback,
        scout.physical.complement(scout.physical.support_of(source_one)),
    )
    common = scout.physical.shift_weighted(
        common_pullback, scout.physical.SHIFT
    )
    right = scout.physical.shift_weighted(
        right_pullback, scout.physical.SHIFT
    )

    cells = {"B": [], "M": [], "R": []}
    for central_index in scout.CENTRAL_INDICES:
        cylinder = scout.target_cylinder(central_index)
        cells["B"].append(
            scout.physical.intersect_weighted(target_one, cylinder)
        )
        cells["M"].append(
            scout.physical.intersect_weighted(common, cylinder)
        )
        cells["R"].append(
            scout.physical.intersect_weighted(right, cylinder)
        )
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
        "the reconstructed horn stopped having the R-M-R word",
    )

    positive_cells = tuple(
        cells["M"][index] or cells["R"][index]
        for index in range(3)
    )
    target_atoms = []
    source_atoms = []
    for cell in positive_cells:
        left, right_endpoint, _weight = cell[0]
        integer_left = scout.ceil_fraction(left)
        require(
            left < integer_left < integer_left + 1 < right_endpoint,
            "a reconstructed horn cylinder lost its internal unit atom",
        )
        target_atom = ((integer_left, integer_left + 1, 1),)
        source_atom = scout.physical.shift_weighted(
            target_atom, -scout.physical.SHIFT
        )
        target_atoms.append(target_atom)
        source_atoms.append(source_atom)
    source_atoms = tuple(source_atoms)
    target_atoms = tuple(target_atoms)
    require(
        (source_atoms, target_atoms) == EXPECTED_ATOMS,
        "the reconstructed canonical endpoint atoms changed",
    )

    factors = scout.allocation.section_factors(
        full_module,
        e3,
        clocks,
        scout.SIGMA,
        scout.TAU,
        scout.CLOCK,
    )
    source_signatures = tuple(
        tuple(
            scout.allocation.contained(atom[0][:2], factor)
            for factor in factors
        )
        for atom in source_atoms
    )
    target_signatures = tuple(
        tuple(
            scout.allocation.contained(atom[0][:2], factor)
            for factor in factors
        )
        for atom in target_atoms
    )
    require(
        source_signatures == EXPECTED_SOURCE_SIGNATURES
        and target_signatures == EXPECTED_TARGET_SIGNATURES,
        "the reconstructed six-factor signatures changed",
    )

    endpoint_base = scout.allocation.endpoint_base
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
        scout.inverse_endpoint_banks(
            source_atom, target_atom, endpoint_context
        )
        for source_atom, target_atom in zip(source_atoms, target_atoms)
    )
    return (
        endpoint_base,
        source_signatures,
        target_signatures,
        endpoint_fields,
    )


def reflection(point):
    a, b = point
    return (a, (2 * a - b) % P)


def bank_digest(points, left, right):
    rows = (
        f"{a},{b}:{left[(a, b)]}:{right[(a, b)]}\n"
        for a, b in points
    )
    return sha256("".join(rows).encode("ascii")).hexdigest()


def compressed_inertia(points, mask):
    """Inertia of D J_r D by exact reflection-orbit decomposition."""
    seen = set()
    positive = 0
    negative = 0
    active_fixed = 0
    active_pairs = 0
    for point in points:
        if point in seen:
            continue
        image = reflection(point)
        if image == point:
            seen.add(point)
            if mask[point]:
                positive += 1
                active_fixed += 1
        else:
            seen.update((point, image))
            if mask[point] and mask[image]:
                positive += 1
                negative += 1
                active_pairs += 1
    zero = len(points) - positive - negative
    frobenius_square = sum(
        bool(mask[point] and mask[reflection(point)])
        for point in points
    )
    return (
        positive,
        negative,
        zero,
        active_fixed,
        active_pairs,
        frobenius_square,
    )


def main():
    (
        endpoint_base,
        source_signatures,
        target_signatures,
        endpoint_fields,
    ) = reconstruct_packet()
    points = tuple(endpoint_base.KEYS)
    require(
        len(points) == P**2
        and len(set(points)) == P**2
        and all(reflection(reflection(point)) == point for point in points),
        "the endpoint universe or reflection involution changed",
    )

    support_rows = []
    origin00_rows = []
    bank_digests = []
    for vertex_index, address in enumerate(scout.ADDRESSES):
        for prime, _zeta, left, right in endpoint_fields[vertex_index]:
            support_rows.append(
                (
                    address,
                    prime,
                    sum(bool(value) for value in left.values()),
                    sum(bool(value) for value in right.values()),
                )
            )
            origin00_rows.append(
                (address, prime, (left[(0, 0)], right[(0, 0)]))
            )
            bank_digests.append(
                (address, prime, bank_digest(points, left, right))
            )
    support_rows = tuple(support_rows)
    origin00_rows = tuple(origin00_rows)
    bank_digests = tuple(bank_digests)
    require(
        support_rows == EXPECTED_SUPPORT_ROWS,
        "the canonical two-field support census changed",
    )
    require(
        origin00_rows == EXPECTED_ORIGIN00_ROWS,
        "the canonical explicit-origin values changed",
    )

    typed_masks = []
    typed_support_rows = []
    inertia_rows = []
    for vertex_index, address in enumerate(scout.ADDRESSES):
        source_gate = all(source_signatures[vertex_index])
        target_gate = all(target_signatures[vertex_index])
        mask = {
            point: bool(
                source_gate
                and target_gate
                and all(
                    left[point] and right[point]
                    for _prime, _zeta, left, right
                    in endpoint_fields[vertex_index]
                )
            )
            for point in points
        }
        require(
            all(mask[point] == mask[reflection(point)] for point in points),
            "a typed mask is not reflection invariant",
        )
        typed_masks.append(mask)
        typed_support_rows.append(
            (address, source_gate, target_gate, sum(mask.values()))
        )
        inertia_rows.append((address,) + compressed_inertia(points, mask))
    typed_support_rows = tuple(typed_support_rows)
    inertia_rows = tuple(inertia_rows)
    require(
        tuple(row[-1] for row in typed_support_rows) == (0, 169, 0),
        "the empty/full/empty typed co-support hostile pattern changed",
    )

    fixed_points = tuple(
        point for point in points if reflection(point) == point
    )
    two_cycles = (len(points) - len(fixed_points)) // 2
    require(
        len(fixed_points) == 13 and two_cycles == 78,
        "the reflection orbit census changed",
    )
    total_inertia = (
        sum(row[1] for row in inertia_rows),
        sum(row[2] for row in inertia_rows),
        sum(row[3] for row in inertia_rows),
    )
    total_trace = sum(row[4] for row in inertia_rows)
    total_frobenius_square = sum(row[6] for row in inertia_rows)
    require(
        total_inertia == (91, 78, 338)
        and total_trace == 13
        and total_frobenius_square == 169,
        "the address-preserving Hermitian packet changed",
    )

    middle_fields = endpoint_fields[1]
    covariance_rows = []
    identity_controls = []
    identity = ((0, 0), (0, 0), 1)
    expected_matches = (13, 13, 0, 0)
    for prime, zeta, left, right in middle_fields:
        reflected_left = {
            point: left[reflection(point)] for point in points
        }
        reflected_right = {
            point: right[reflection(point)] for point in points
        }
        left_control = scout.affine_covariances(
            left, left, zeta, prime, points
        )
        right_control = scout.affine_covariances(
            right, right, zeta, prime, points
        )
        require(
            identity in left_control and identity in right_control,
            "the exact Heisenberg search lost its identity control",
        )
        identity_controls.append(
            (prime, identity in left_control, identity in right_control)
        )
        cases = (
            ("L_to_reflected_L", left, reflected_left),
            ("R_to_reflected_R", right, reflected_right),
            ("L_to_reflected_R", left, reflected_right),
            ("R_to_reflected_L", right, reflected_left),
        )
        field_rows = []
        for case_name, bare, carried in cases:
            matches = sum(
                bare[point] == carried[point] for point in points
            )
            covariances = scout.affine_covariances(
                bare, carried, zeta, prime, points
            )
            field_rows.append((case_name, matches, covariances))
        require(
            tuple(row[1] for row in field_rows) == expected_matches
            and all(not row[2] for row in field_rows),
            "reflection entered a middle finite-Heisenberg covariance class",
        )
        covariance_rows.append((prime, tuple(field_rows)))
    identity_controls = tuple(identity_controls)
    covariance_rows = tuple(covariance_rows)

    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_nodes = sum(
        isinstance(node, ast.Assert) for node in ast.walk(tree)
    )
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(
        assert_nodes == 0 and float_literals == 0,
        "the sidecar source stopped being exact under optimized Python",
    )

    print("LRC14 C13 REFLECTION/GABOR HOSTILE SIDECAR")
    print("status=FINITE-EXACT HOSTILE SIDECAR;not_a_canonical_theorem")
    print(f"reproduction_command={REPRODUCTION_COMMAND}")
    print(f"stored_output={STORED_OUTPUT_PATH}")
    print("dependency_hashes=BEGIN")
    for path, digest in ACTUAL_HASHES:
        print(f"{digest}  {path}")
    print("dependency_hashes=END")
    print(
        "universe=THM3285_unit_subatom;"
        "(clock,sigma,tau)=(1,0,3);"
        f"addresses={scout.ADDRESSES};endpoint_origins=F13^2={len(points)};"
        "certified_fields=2"
    )
    print(
        "reflection=r(a,b)=(a,2a-b)_mod13;"
        f"fixed_points={len(fixed_points)};two_cycles={two_cycles};"
        "fixed_locus=b=a"
    )
    print(
        "matrix=direct_sum_over_addresses(D_v*J_r*D_v);"
        "J_r[e_(a,b)]=e_(a,2a-b);"
        "D_v=typed_common_origin_mask_across_both_certified_fields"
    )
    print(
        "factor_sidecars=(E3,clock,q1,q2,c2,c3);"
        f"source_signatures={source_signatures};"
        f"target_signatures={target_signatures}"
    )
    print(f"raw_support_rows_(address,field,left,right)={support_rows}")
    print(f"origin00_provenance_(address,field,(left,right))={origin00_rows}")
    print(f"bank_sha256_rows_(address,field,digest)={bank_digests}")
    print(
        "typed_support_rows_(address,source_gate,target_gate,common_origins)="
        f"{typed_support_rows}"
    )
    print(
        "block_inertia_rows_(address,nplus,nminus,nzero,"
        "active_fixed,active_pairs,frobenius_sq)="
        f"{inertia_rows}"
    )
    print(
        f"total_inertia_(nplus,nminus,nzero)={total_inertia};"
        f"trace={total_trace};frobenius_sq={total_frobenius_square}"
    )
    print(
        "heisenberg_gauge=f'(x)=c*zeta^(u.x)*f(x+s);"
        "exhausts_all_s,u_in_F13^2_and_nonzero_scalar_for_full_support"
    )
    print(
        "identity_positive_controls_(field,left,right)="
        f"{identity_controls}"
    )
    print(
        "reflection_covariance_rows_"
        "(field,((case,strict_matches,covariances),...))="
        f"{covariance_rows}"
    )
    print(
        "PASS_NEGATIVE=middle_support_is_reflection_invariant_but_neither_"
        "middle_bank_nor_left_right_swap_lies_in_its_reflected_exact_"
        "finite_Heisenberg_orbit_in_either_certified_field"
    )
    print(
        "VERDICT=the_(91,78,338)_support_inertia_is_an_involution_orbit_"
        "constant_that_fires_on_the_hostile_middle_fibre_while_both_outer_"
        "typed_blocks_are_zero;no_native_Anthropic_inertia_to_LRC_transfer"
    )
    print(
        "scope=one_canonical_THM3285_unit_subatom;address_and_outer_"
        "co-support_sidecars_retained;no_row_exclusion;no_owner_current_"
        "or_exit;no_LRC14"
    )
    print(
        f"source_ast=(assert_nodes={assert_nodes},"
        f"float_literals={float_literals})"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
