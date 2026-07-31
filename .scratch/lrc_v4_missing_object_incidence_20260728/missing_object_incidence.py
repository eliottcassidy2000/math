#!/usr/bin/env python3
"""Exact incidence audit for the rank-587 missing THM-2859 V4 object.

Scratch status only.  The script distinguishes:

* a literal existing collar vertex;
* an integral signed combination of the 587 whole rooted paths;
* the virtual marginal identity R+M1-M2; and
* the full joint semantic/carrier/factor mask required by X=Pi(R).

No executable Python assertion statement is used.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import ast
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge" / "results"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_digest(path):
    return sha256(
        path.read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()


PINS = {
    COMP / "lrc14_nearest_half_step_common_right_collar_thm2825.py":
        "bd9ffe7f6815b5c563bd483c300118fbdd683f3d9303babbab7912e031747c9a",
    RESULTS / "lrc14_nearest_half_step_common_right_collar_thm2825.out":
        "c4a31e5ee0aa5af69faa3efbe315d0968a85cba49c2d77c0ca93a229bc39fa0c",
    COMP / "lrc14_nearest_half_step_common_right_collar_path_operator_thm2825.py":
        "7f0780e70161cbfafc4499d02a7d1c5aa40366b6dfa9935b00dc652e3d54c8e0",
    RESULTS / "lrc14_nearest_half_step_common_right_collar_path_operator_thm2825.out":
        "6cff846944c59d8e243e74afe8829046feb288f0cf23da1c998037dcf70411b2",
    COMP / "lrc14_horn_collar_v4_globalization_thm2859.py":
        "f7b7ee404f7a55eac10af0f2da19069197cad67c5a943326c070175233d289c6",
    RESULTS / "lrc14_horn_collar_v4_globalization_thm2859.out":
        "1719513ea3e009c5ae22c4bba618cd7c16285bd41f881233d8a5f68ca58e3bdc",
}
for path, expected in PINS.items():
    require(
        lf_digest(path) == expected,
        f"pinned dependency changed: {path.name}",
    )


import lrc14_nearest_half_step_common_right_collar_thm2825 as collar


FACTOR_NAMES = ("E3", "clock", "q1", "q2", "c2", "c3")
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)
H = collar.H


def safe_comb_contains(interval, module, speed, denominator, lo, hi):
    require(
        module.T % (denominator * speed) == 0,
        "comb grid ceased to resolve",
    )
    unit = module.T // (denominator * speed)
    period = denominator * unit
    width = (hi - lo) * unit
    base = (lo % denominator) * unit
    left, right = interval
    first = (left - base) // period
    return not any(
        max(left, base + number * period)
        < min(right, base + number * period + width)
        for number in range(first - 1, first + 3)
    )


def root_factor_signature(
    interval, module, e3, clocks, clock, sigma, target
):
    copies = collar.copies
    return (
        copies.contained(interval, e3),
        copies.contained(interval, clocks[clock]),
        safe_comb_contains(
            interval, module, module.W[1], 182,
            -14 * sigma - 13, -14 * sigma + 13,
        ),
        safe_comb_contains(
            interval, module, module.W[2], 182,
            -14 * target - 13, -14 * target + 13,
        ),
        safe_comb_contains(
            interval, module, module.C2, 182,
            14 * sigma - 13, 14 * sigma + 13,
        ),
        safe_comb_contains(
            interval, module, module.C3, 182,
            14 * target - 13, 14 * target + 13,
        ),
    )


def gf2_incidence_certificate(rows):
    """Return rank and one exact augmented inconsistency certificate."""
    pivots = {}
    first_certificate = None
    ordered = sorted(
        enumerate(rows),
        key=lambda pair: (
            pair[1][0].bit_count(),
            repr(pair[1][2]),
        ),
    )
    for original_index, (bits_initial, rhs_initial, _key) in ordered:
        bits = bits_initial
        rhs = rhs_initial
        certificate = 1 << original_index
        while bits:
            pivot = (bits & -bits).bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (bits, rhs, certificate)
                break
            old_bits, old_rhs, old_certificate = pivots[pivot]
            bits ^= old_bits
            rhs ^= old_rhs
            certificate ^= old_certificate
        if not bits and rhs and first_certificate is None:
            first_certificate = certificate

    if first_certificate is not None:
        check_bits = 0
        check_rhs = 0
        support = []
        for index, (bits, rhs, key) in enumerate(rows):
            if (first_certificate >> index) & 1:
                check_bits ^= bits
                check_rhs ^= rhs
                support.append(key)
        require(
            check_bits == 0 and check_rhs == 1,
            "GF(2) inconsistency certificate failed replay",
        )
    else:
        support = []
    return len(pivots), tuple(support)


def corner_features(semantic, carrier):
    """Dimension and the two separate sign characters."""
    return (
        1,
        (-1) ** semantic,
        (-1) ** carrier,
    )


def joint_character(semantic, carrier):
    return (-1) ** (semantic + carrier)


def determinant3(columns):
    matrix = tuple(zip(*columns))
    return (
        matrix[0][0] * (
            matrix[1][1] * matrix[2][2]
            - matrix[1][2] * matrix[2][1]
        )
        - matrix[0][1] * (
            matrix[1][0] * matrix[2][2]
            - matrix[1][2] * matrix[2][0]
        )
        + matrix[0][2] * (
            matrix[1][0] * matrix[2][1]
            - matrix[1][1] * matrix[2][0]
        )
    )


def determinant4(columns):
    matrix = tuple(zip(*columns))
    total = 0
    for column in range(4):
        minor_columns = []
        for kept_column in range(4):
            if kept_column == column:
                continue
            minor_columns.append(
                tuple(matrix[row][kept_column] for row in range(1, 4))
            )
        total += (
            (-1) ** column
            * matrix[0][column]
            * determinant3(tuple(minor_columns))
        )
    return total


def main():
    (
        _module,
        _rails,
        _present,
        details,
        full_module,
        e3,
        clocks,
        q_pairs,
        _delayed,
        _source_weight,
        _target_weight,
        _rail_common,
    ) = collar.copies.physical_setup()

    semantic_caches = tuple({} for _clock in range(7))
    records = []
    labelled_vertex_owner = {}
    physical_rows = defaultdict(int)
    physical_occurrences = defaultdict(list)
    physical_root_rhs = Counter()
    rooted_common_count = 0

    for clock in range(7):
        for sigma in COMMON_S:
            for target in COMMON_T:
                (
                    _source,
                    _target_physical,
                    common,
                    right,
                ) = collar.cell_objects(
                    details,
                    full_module,
                    e3,
                    clocks,
                    clock,
                    sigma,
                    target,
                )
                if not right:
                    continue
                cell = (clock, sigma, target)
                common_by_left = {piece[0]: piece for piece in common}
                require(
                    len(common_by_left) == len(common),
                    "common interval left endpoints collided",
                )
                for root_index, root in enumerate(right):
                    path = [root]
                    cursor = root[0] + H
                    while cursor in common_by_left:
                        path.append(common_by_left[cursor])
                        cursor += H
                    require(
                        len(path) >= 3,
                        "rooted path lost its first two common vertices",
                    )

                    path_index = len(records)
                    for position, piece in enumerate(path):
                        labelled_key = (
                            cell,
                            piece[0],
                            piece[1],
                            piece[2],
                        )
                        require(
                            labelled_key not in labelled_vertex_owner,
                            "two rooted paths met in the labelled forest",
                        )
                        labelled_vertex_owner[labelled_key] = (
                            path_index,
                            position,
                        )
                        bit = 1 << path_index
                        require(
                            not physical_rows[piece] & bit,
                            "one path revisited a physical interval",
                        )
                        physical_rows[piece] |= bit
                        physical_occurrences[piece].append(
                            (
                                path_index,
                                position,
                                cell,
                                root_index,
                            )
                        )
                    rooted_common_count += len(path) - 1
                    physical_root_rhs[root] += 1

                    root_live = (
                        collar.semantic_value(
                            root,
                            q_pairs[clock],
                            semantic_caches[clock],
                        )
                        != (0, 0)
                    )
                    factor_signature = root_factor_signature(
                        root[:2],
                        full_module,
                        e3,
                        clocks,
                        clock,
                        sigma,
                        target,
                    )
                    holes = tuple(
                        name
                        for name, present in zip(
                            FACTOR_NAMES, factor_signature
                        )
                        if not present
                    )
                    m1_live = (
                        collar.semantic_value(
                            path[1],
                            q_pairs[clock],
                            semantic_caches[clock],
                        )
                        != (0, 0)
                    )
                    m2_live = (
                        collar.semantic_value(
                            path[2],
                            q_pairs[clock],
                            semantic_caches[clock],
                        )
                        != (0, 0)
                    )
                    require(
                        (m1_live, m2_live)
                        == (not root_live, root_live),
                        "initial ladder semantic pattern changed",
                    )
                    records.append({
                        "cell": cell,
                        "root_index": root_index,
                        "root": root,
                        "path": tuple(path),
                        "root_live": root_live,
                        "holes": holes,
                    })

    require(
        len(records) == 587
        and rooted_common_count == 54754
        and len(labelled_vertex_owner) == 587 + 54754,
        "rooted forest census changed",
    )

    # Literal vertex test.  A desired X has R's carrier/factor type but
    # M1's semantic value.  Common vertices have the wrong carrier type, so
    # only roots can be candidates.  The root type census is already enough
    # to rule every candidate out.
    existing_root_types = Counter(
        (record["holes"], record["root_live"])
        for record in records
    )
    desired_root_types = Counter(
        (record["holes"], not record["root_live"])
        for record in records
    )
    require(
        existing_root_types
        == Counter({
            (("E3",), True): 319,
            (("E3", "c2"), True): 37,
            (("c2",), True): 217,
            (("q1",), False): 14,
        })
        and desired_root_types
        == Counter({
            (("E3",), False): 319,
            (("E3", "c2"), False): 37,
            (("c2",), False): 217,
            (("q1",), True): 14,
        })
        and set(existing_root_types).isdisjoint(desired_root_types),
        "root/desired decorated type boundary changed",
    )

    physical_root_semantics = defaultdict(set)
    physical_root_hole_semantics = set()
    desired_physical_root_hole_semantics = set()
    for record in records:
        root = record["root"]
        live = record["root_live"]
        holes = record["holes"]
        physical_root_semantics[root].add(live)
        physical_root_hole_semantics.add((root, holes, live))
        desired_physical_root_hole_semantics.add(
            (root, holes, not live)
        )
    require(
        len(physical_root_semantics) == 37
        and all(
            len(values) == 1
            for values in physical_root_semantics.values()
        )
        and len(physical_root_hole_semantics) == 42
        and physical_root_hole_semantics.isdisjoint(
            desired_physical_root_hole_semantics
        ),
        "physical root mask-class refinement changed",
    )

    # Whole-path incidence over the labelled forest.  Every column contains
    # its unique root row and at least one unique common row.  Asking for the
    # root selector therefore forces z_i=1 and z_i=0 on every path.
    local_path_contradictions = tuple(
        (
            index,
            record["cell"],
            record["root_index"],
            record["root"],
            record["path"][1],
        )
        for index, record in enumerate(records)
    )
    require(
        len(local_path_contradictions) == 587,
        "labelled whole-path contradiction census changed",
    )
    labelled_path_rank = 587
    labelled_augmented_rank = 588

    # Even after forgetting cell labels and allowing arbitrary signed
    # multiplicities, reduce the physical incidence equations modulo two.
    # A GF(2) inconsistency is an integral no-go.
    physical_equations = []
    all_physical_keys = set(physical_rows) | set(physical_root_rhs)
    for key in sorted(all_physical_keys):
        physical_equations.append(
            (
                physical_rows.get(key, 0),
                physical_root_rhs.get(key, 0) & 1,
                key,
            )
        )
    physical_rank_mod2, certificate_support = gf2_incidence_certificate(
        physical_equations
    )
    require(
        certificate_support,
        "physical label-forgotten whole-path system became GF(2)-consistent",
    )
    certificate_payload = repr(certificate_support).encode()
    certificate_digest = sha256(certificate_payload).hexdigest()
    certificate_details = tuple(
        (
            key,
            physical_rows[key].bit_count(),
            physical_root_rhs.get(key, 0),
            tuple(physical_occurrences[key]),
        )
        for key in certificate_support
    )
    certificate_integer_rhs = tuple(
        physical_root_rhs.get(key, 0)
        for key in certificate_support
    )
    certificate_support_rhs = tuple(
        int(physical_root_rhs.get(key, 0) != 0)
        for key in certificate_support
    )
    require(
        len(certificate_details) == 2
        and certificate_details[0][1]
        == certificate_details[1][1]
        and physical_rows[certificate_support[0]]
        == physical_rows[certificate_support[1]]
        and certificate_integer_rhs == (7, 0)
        and certificate_support_rhs == (1, 0)
        and tuple(
            occurrence[0]
            for occurrence in certificate_details[0][3]
        )
        == tuple(
            occurrence[0]
            for occurrence in certificate_details[1][3]
        )
        and all(
            occurrence[1] == 0
            for occurrence in certificate_details[0][3]
        )
        and all(
            occurrence[1] == 1
            for occurrence in certificate_details[1][3]
        ),
        "minimal two-row physical certificate changed",
    )

    # The exact virtual survivor.  If the joint semantic/carrier incidence
    # is quotiented to dimension plus its two separate marginals, then
    # X=R+M1-M2.  The mixed character detects the discarded interaction.
    corners = {
        "R": (0, 0),
        "M1": (1, 1),
        "M2": (0, 1),
        "X": (1, 0),
    }
    virtual_coefficients = {"R": 1, "M1": 1, "M2": -1}
    marginal_sum = tuple(
        sum(
            virtual_coefficients[name]
            * corner_features(*corners[name])[coordinate]
            for name in virtual_coefficients
        )
        for coordinate in range(3)
    )
    target_marginals = corner_features(*corners["X"])
    mixed_virtual = sum(
        virtual_coefficients[name]
        * joint_character(*corners[name])
        for name in virtual_coefficients
    )
    mixed_target = joint_character(*corners["X"])
    mixed_defect_per_root = mixed_virtual - mixed_target
    old_marginal_determinant = determinant3(tuple(
        corner_features(*corners[name])
        for name in ("R", "M1", "M2")
    ))
    full_joint_determinant = determinant4(tuple(
        corner_features(*corners[name])
        + (joint_character(*corners[name]),)
        for name in ("R", "M1", "M2", "X")
    ))
    require(
        marginal_sum == target_marginals
        and virtual_coefficients == {"R": 1, "M1": 1, "M2": -1}
        and abs(old_marginal_determinant) == 4
        and abs(full_joint_determinant) == 16
        and mixed_defect_per_root == 4,
        "virtual marginal identity changed",
    )

    # The M1 corner has exactly the desired semantic population, but not the
    # desired carrier/factor type.  This names X as the missing fibre product.
    r_semantics = Counter(
        record["root_live"] for record in records
    )
    m1_semantics = Counter(
        not record["root_live"] for record in records
    )
    desired_semantics = Counter(
        not record["root_live"] for record in records
    )
    require(
        r_semantics == Counter({True: 573, False: 14})
        and m1_semantics == desired_semantics
        == Counter({False: 573, True: 14}),
        "semantic fibre-product counts changed",
    )

    undecorated_carrier_masks = {
        "source_empty_target_delta0"
        for _record in records
    }
    factor_refined_carrier_masks = {
        ("source_empty_target_delta0", record["holes"])
        for record in records
    }
    desired_decorated_masks = set(desired_root_types)
    require(
        len(undecorated_carrier_masks) == 1
        and len(factor_refined_carrier_masks) == 4
        and len(desired_decorated_masks) == 4,
        "mask refinement count changed",
    )

    source_tree = ast.parse(Path(__file__).read_text())
    require(
        not any(
            isinstance(node, ast.Assert)
            for node in ast.walk(source_tree)
        ),
        "executable assertion statement found",
    )

    print("THM-2859 RANK-587 MISSING-OBJECT INCIDENCE AUDIT")
    print(
        "pinned="
        f"{tuple((path.name, digest) for path, digest in PINS.items())}"
    )
    print(
        f"universe=paths:{len(records)},rooted_common:{rooted_common_count},"
        f"labelled_vertices:{len(labelled_vertex_owner)},"
        f"physical_vertices:{len(all_physical_keys)}"
    )
    print(
        f"existing_root_types={tuple(sorted(existing_root_types.items(), key=str))}"
    )
    print(
        f"desired_X_types={tuple(sorted(desired_root_types.items(), key=str))};"
        "intersection=empty"
    )
    print(
        "literal_vertex_no_go="
        "X needs R carrier/factor type and M1 semantic;"
        " roots have the first but wrong semantic, common vertices have the "
        "second on alternating levels but wrong carrier"
    )
    print(
        "physical_root_refinement="
        f"intervals:{len(physical_root_semantics)},"
        f"interval_hole_semantic_classes:{len(physical_root_hole_semantics)},"
        "semantic_conflicts:0;desired_intersection:0"
    )
    print(
        "labelled_whole_path_system="
        "A*z=root_selector;each path gives equations z_i=1 at its root "
        "and z_i=0 at its first common vertex;"
        f"contradictions:{len(local_path_contradictions)},"
        f"rank_A:{labelled_path_rank},"
        f"rank_augmented:{labelled_augmented_rank}"
    )
    print(
        "first_labelled_witness="
        f"{local_path_contradictions[0]}"
    )
    print(
        "physical_quotient_signed_system="
        "A_phys*z=root_multiplicity;"
        f"equations:{len(physical_equations)},"
        f"rank_mod2:{physical_rank_mod2},"
        f"rank_augmented_mod2:{physical_rank_mod2 + 1},"
        f"certificate_rows:{len(certificate_support)},"
        f"certificate_digest:{certificate_digest},"
        f"identical_integer_rows_rhs:{certificate_integer_rhs},"
        f"identical_support_rows_rhs:{certificate_support_rhs},"
        f"certificate_details:{certificate_details}"
    )
    print(
        "virtual_marginal_survivor="
        f"X=R+M1-M2;coefficients:{virtual_coefficients};"
        f"old_marginal_determinant:{old_marginal_determinant};"
        f"full_joint_determinant:{full_joint_determinant};"
        f"marginals:{marginal_sum};"
        f"mixed_character_virtual_target_defect:"
        f"{(mixed_virtual, mixed_target, mixed_defect_per_root)};"
        f"all_roots_defect:{587 * mixed_defect_per_root}"
    )
    print(
        f"semantic_counts=R:{tuple(sorted(r_semantics.items()))},"
        f"M1=X:{tuple(sorted(m1_semantics.items()))}"
    )
    print(
        "mask_class_answer="
        "no new undecorated carrier support is needed, but one faithful "
        "joint parity tag splits into 4 absent factor-refined classes, "
        "42 absent physical interval/hole/semantic classes, and 587 "
        "labelled atoms"
    )
    print(
        "cheapest_surviving_enlargement="
        "adjoin the parity-tagged root copy X=Pi(R), equivalently the "
        "rootwise fibre product (R carrier/factors) x_(root index) "
        "(M1 semantic);dimension=587;one new joint tag;no smaller actual "
        "filler can factor rank(S)=587"
    )
    print(
        "scope=exact labelled and physical-incidence obstruction;"
        "the virtual marginal identity is not an object, positive packet, "
        "or typed V4 globalization"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
