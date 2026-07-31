#!/usr/bin/env python3
"""Exact q3/q11/q7-horn versus half-step-collar hinge audit.

The computation composes promoted THM-2825, THM-2847, and THM-2851 without
identifying their differently typed actions.  It finds the exact q=0 physical
hinge and checks the split-versus-nonsplit 169-state obstruction.  It proves no
row exclusion and no LRC(14) conclusion.

No Python ``assert`` statement is used.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge" / "results"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    COMP / "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
    COMP / "lrc14_nearest_half_step_common_right_collar_thm2825.py":
        "bd9ffe7f6815b5c563bd483c300118fbdd683f3d9303babbab7912e031747c9a",
    COMP / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    COMP / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    RESULTS / "lrc14_nearest_half_step_common_right_collar_thm2825.out":
        "c4a31e5ee0aa5af69faa3efbe315d0968a85cba49c2d77c0ca93a229bc39fa0c",
    RESULTS / "lrc14_right_cofiber_positive_copy_stratification_thm2818.out":
        "225bd77c27d5972e7dad663e46be3c4c20e2b9449615018773e6822360356a33",
    RESULTS / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.out":
        "155fce129c750a9505fdda3c71a250ff3a57fcd4044bb1df941da83c08baee1d",
    RESULTS / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out":
        "424fd2e83a618f862a5ee1b5f073a93fe236d10cdc5412eab1b54dec5e537eac",
}
for path, expected in PINNED.items():
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"pinned dependency changed: {path.name}",
    )


import lrc14_literal_fixed_sheet_allocation_thm2806 as allocation
import lrc14_nearest_half_step_common_right_collar_thm2825 as collar


P = 13
H = collar.H
I = allocation.ATOM_INTERVAL
FACTOR_NAMES = ("E3", "clock", "q1", "q2", "c2", "c3")


def section_factors(full_module, e3, clocks, clock, sigma, target):
    universe = ((0, full_module.T),)
    return (
        e3,
        clocks[clock],
        full_module.subtract_comb(
            universe, full_module.W[1], 182,
            -14 * sigma - 13, -14 * sigma + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.W[2], 182,
            -14 * target - 13, -14 * target + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C2, 182,
            14 * sigma - 13, 14 * sigma + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C3, 182,
            14 * target - 13, 14 * target + 13,
        ),
    )


def carrier_supports(source, target):
    unit = collar.copies.T // P
    return (
        tuple(
            collar.copies.support_of(
                collar.copies.physical.overlap.shift_weighted(
                    source, -twist * unit
                )
            )
            for twist in range(P)
        ),
        tuple(
            collar.copies.support_of(
                collar.copies.physical.overlap.shift_weighted(
                    target, twist * unit
                )
            )
            for twist in range(P)
        ),
    )


def carrier_mask(piece, supports):
    interval = piece[:2]
    target_interval = tuple(endpoint + collar.copies.SHIFT for endpoint in interval)
    return (
        tuple(
            bool(collar.copies.intersect_sorted((interval,), support))
            for support in supports[0]
        ),
        tuple(
            bool(collar.copies.intersect_sorted((target_interval,), support))
            for support in supports[1]
        ),
    )


def mat_mul(left, right, p):
    return (
        (
            (left[0][0] * right[0][0] + left[0][1] * right[1][0]) % p,
            (left[0][0] * right[0][1] + left[0][1] * right[1][1]) % p,
        ),
        (
            (left[1][0] * right[0][0] + left[1][1] * right[1][0]) % p,
            (left[1][0] * right[0][1] + left[1][1] * right[1][1]) % p,
        ),
    )


def mat_pow(matrix, exponent, p):
    answer = ((1, 0), (0, 1))
    base = matrix
    while exponent:
        if exponent & 1:
            answer = mat_mul(answer, base, p)
        base = mat_mul(base, base, p)
        exponent //= 2
    return answer


def mat_vec(matrix, vector, p):
    return (
        (matrix[0][0] * vector[0] + matrix[0][1] * vector[1]) % p,
        (matrix[1][0] * vector[0] + matrix[1][1] * vector[1]) % p,
    )


def affine_pow(linear, shift, exponent, p):
    out_linear = ((1, 0), (0, 1))
    out_shift = (0, 0)
    base_linear = linear
    base_shift = shift
    while exponent:
        if exponent & 1:
            out_shift = tuple(
                (x + y) % p
                for x, y in zip(
                    mat_vec(out_linear, base_shift, p),
                    out_shift,
                )
            )
            out_linear = mat_mul(out_linear, base_linear, p)
        base_shift = tuple(
            (x + y) % p
            for x, y in zip(
                mat_vec(base_linear, base_shift, p),
                base_shift,
            )
        )
        base_linear = mat_mul(base_linear, base_linear, p)
        exponent //= 2
    return out_linear, out_shift


def affine_p_primary_audit(p):
    identity = ((1, 0), (0, 1))
    linear_candidates = []
    affine_candidates = 0
    nontrivial_after_p = 0
    order_histogram = Counter()
    for a in range(p):
        for b in range(p):
            for c in range(p):
                for d in range(p):
                    matrix = ((a, b), (c, d))
                    determinant = (a * d - b * c) % p
                    if determinant == 0 or mat_pow(matrix, p, p) != identity:
                        continue
                    linear_candidates.append(matrix)
                    for x in range(p):
                        for y in range(p):
                            affine_candidates += 1
                            power_p = affine_pow(matrix, (x, y), p, p)
                            require(
                                power_p == (identity, (0, 0)),
                                "an affine p-primary candidate survived p steps",
                            )
                            if power_p != (identity, (0, 0)):
                                nontrivial_after_p += 1
                            if matrix == identity and (x, y) == (0, 0):
                                order_histogram[1] += 1
                            else:
                                order_histogram[p] += 1
    return (
        len(linear_candidates),
        affine_candidates,
        nontrivial_after_p,
        tuple(sorted(order_histogram.items())),
    )


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

    # THM-2847 proves that these, and only these, are its E3-only horn cells.
    horn = tuple(
        (1, sigma, target)
        for sigma in (0, 3, 8, 9, 12)
        for target in (5, 6, 9, 10)
    )

    delta_zero = (True,) + (False,) * 12
    zero = (False,) * 13
    interval_q_membership = Counter()
    root_count = 0
    common_count = 0
    distinguished_path_lengths = Counter()
    distinguished_semantics = Counter()
    factor_triples = Counter()
    carrier_triples = Counter()
    physical_roots = Counter()
    physical_m1 = Counter()
    physical_m2 = Counter()
    distinguished_rows = []

    unit = full_module.T // P
    for cell in horn:
        clock, sigma, target = cell
        source, target_physical, common, right = collar.cell_objects(
            details,
            full_module,
            e3,
            clocks,
            clock,
            sigma,
            target,
        )
        require(right and common, "a horn cell left the collar bank")
        root_count += len(right)
        common_count += len(common)
        common_by_left = {piece[0]: piece for piece in common}
        right_by_left = {piece[0]: piece for piece in right}

        rpiece = right_by_left.get(I[0] - 2 * H)
        m1 = common_by_left.get(I[0] - H)
        m2 = common_by_left.get(I[0])
        require(
            rpiece is not None and m1 is not None and m2 is not None,
            "distinguished q0 two-step hinge disappeared",
        )
        require(
            rpiece[:2] == collar.copies.SELECTED_MINUS
            and m2[:2] == I,
            "distinguished physical triple changed",
        )

        length = 0
        left = m1[0]
        while left in common_by_left:
            length += 1
            left += H
        distinguished_path_lengths[length] += 1
        semantic = tuple(
            collar.semantic_value(piece, q_pairs[clock], {}) != (0, 0)
            for piece in (rpiece, m1, m2)
        )
        require(
            semantic == (True, False, True),
            "distinguished hinge semantic grading changed",
        )
        distinguished_semantics[semantic] += 1

        factors = section_factors(
            full_module, e3, clocks, clock, sigma, target
        )
        factor_rows = tuple(
            collar.copies.factor_masks(piece[:2], factors)
            for piece in (rpiece, m1, m2)
        )
        factor_triples[factor_rows] += 1
        require(
            factor_rows[0]
            == (
                (False, True, True, True, True, True),
                (True,) * 6,
                (True,) * 6,
                (False, True, True, True, True, True),
            )
            and factor_rows[1] == ((True,) * 6,) * 4
            and factor_rows[2] == ((True,) * 6,) * 4,
            "uniform E3-only collar-root defect changed",
        )

        supports = carrier_supports(source, target_physical)
        carrier_rows = tuple(
            carrier_mask(piece, supports) for piece in (rpiece, m1, m2)
        )
        carrier_triples[carrier_rows] += 1
        require(
            carrier_rows
            == (
                (zero, delta_zero),
                (delta_zero, delta_zero),
                (delta_zero, delta_zero),
            ),
            "uniform carrier hinge changed",
        )

        for q in range(P):
            q_interval = tuple(endpoint + q * unit for endpoint in I)
            if q_interval in {piece[:2] for piece in right}:
                membership = "right"
            elif q_interval in {piece[:2] for piece in common}:
                membership = (
                    "M2"
                    if q_interval == m2[:2]
                    else "other_common"
                )
            else:
                membership = "outside"
            interval_q_membership[q, membership] += 1

        physical_roots[rpiece[:2]] += 1
        physical_m1[m1[:2]] += 1
        physical_m2[m2[:2]] += 1
        distinguished_rows.append(
            (cell, len(common), len(right), length)
        )

    require(
        len(horn) == 20 and root_count == 52 and common_count == 4076,
        "horn/collar cell census changed",
    )
    require(
        distinguished_path_lengths
        == Counter({144: 8, 14: 4, 40: 4, 118: 4}),
        "distinguished path lengths changed",
    )
    require(
        interval_q_membership
        == Counter(
            {(0, "M2"): 20}
            | {(q, "outside"): 20 for q in range(1, P)}
        ),
        "allocation-orbit/collar fibre product changed",
    )
    require(
        len(physical_roots) == len(physical_m1) == len(physical_m2) == 1
        and tuple(physical_roots.values()) == (20,)
        and tuple(physical_m1.values()) == (20,)
        and tuple(physical_m2.values()) == (20,),
        "label forgetting multiplicity changed",
    )

    # THM-2818 independently proves the selected source endpoint masks:
    # R=SELECTED_MINUS has size 0, while M2=I has size 81.  Their target
    # masks agree and have size 81.  The output pin above protects this
    # inherited statement; unequal cardinality rules out every permutation,
    # not only split endpoint translations.
    endpoint_boundary = (0, 81, 81, 81)

    # Formal grade bookkeeping only.  THM-2851's algebraic horn attachment
    # changes (outer-word,E3) by (1,1).  The even collar changes
    # (semantic,carrier) by (0,1).  If the two second coordinates could be
    # physically identified, their diagonal sum would be the missing (1,0).
    horn_attachment_degree = (1, 1)
    even_collar_degree = (0, 1)
    formal_diagonal_degree = tuple(
        (x + y) % 2
        for x, y in zip(horn_attachment_degree, even_collar_degree)
    )
    require(
        formal_diagonal_degree == (1, 0),
        "formal punctured-V4 degree cancellation changed",
    )

    # Exact p=13 check of the general affine-exponent lemma.  The p=2 run is
    # the sharp hostile: AGL(2,2)=S4 does contain order-four elements.
    affine_13 = affine_p_primary_audit(13)
    require(
        affine_13[2] == 0
        and affine_13[3] == ((1, 1), (13, affine_13[1] - 1)),
        "AGL(2,13) acquired a 169-cycle",
    )

    # Find the p=2 boundary directly rather than reusing the odd-prime proof.
    identity2 = ((1, 0), (0, 1))
    order4_witnesses = []
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    matrix = ((a, b), (c, d))
                    if (a * d - b * c) % 2 == 0:
                        continue
                    for x in range(2):
                        for y in range(2):
                            g2 = affine_pow(matrix, (x, y), 2, 2)
                            g4 = affine_pow(matrix, (x, y), 4, 2)
                            if (
                                g2 != (identity2, (0, 0))
                                and g4 == (identity2, (0, 0))
                            ):
                                order4_witnesses.append((matrix, (x, y)))
    require(order4_witnesses, "p=2 affine order-four boundary disappeared")

    profile_histogram = Counter(
        (cell[1], common, roots, length)
        for cell, common, roots, length in distinguished_rows
    )
    print("THM-2859 Q3/Q11/Q7 HORN -- HALF-STEP COLLAR -- WITT HINGE")
    print("status=PROVED + FINITE-EXACT; no row or LRC14 conclusion")
    print(
        f"horn_cells={len(horn)}; collar_roots={root_count}; "
        f"collar_common={common_count}; "
        f"distinguished_hinges={len(distinguished_rows)}"
    )
    print(
        "distinguished_profiles="
        f"{tuple(sorted(profile_histogram.items()))};"
        "target_labels=(5,6,9,10)"
    )
    print(
        "physical_hinge="
        f"R={next(iter(physical_roots))};"
        f"M1={next(iter(physical_m1))};"
        f"M2={next(iter(physical_m2))};"
        "labels_per_physical_hinge=20"
    )
    print(
        f"distinguished_path_lengths="
        f"{tuple(sorted(distinguished_path_lengths.items()))};"
        f"semantic={tuple(sorted(distinguished_semantics.items()))}"
    )
    print(
        "factor_hinge=R:011111/111111/111111/011111;"
        "M1:111111^4;M2:111111^4;"
        f"uniform_labels={next(iter(factor_triples.values()))}"
    )
    print(
        "carrier_hinge=R:zero/delta0;"
        "M1:delta0/delta0;M2:delta0/delta0;"
        f"uniform_labels={next(iter(carrier_triples.values()))}"
    )
    print(
        "allocation_interval_membership="
        f"{tuple(sorted(interval_q_membership.items()))};"
        "q3_q7_q11_physical_fibre_product=empty"
    )
    print(
        "endpoint_boundary="
        f"(R_source={endpoint_boundary[0]},M2_source={endpoint_boundary[1]},"
        f"R_target={endpoint_boundary[2]},M2_target={endpoint_boundary[3]});"
        "source_cardinality_obstruction=81;target_identity_match=1"
    )
    print(
        f"formal_grades=horn{horn_attachment_degree}+"
        f"collar{even_collar_degree}={formal_diagonal_degree};"
        "status=grade_quotient_only_not_physical_composition"
    )
    print(
        f"AGL2_F13_p_primary="
        f"(linear_candidates={affine_13[0]},"
        f"affine_candidates={affine_13[1]},"
        f"order_histogram={affine_13[3]});"
        "no_order169=1"
    )
    print(
        f"p2_boundary_order4_witness_count={len(order4_witnesses)};"
        f"first={order4_witnesses[0]}"
    )
    print(
        "state_invoice=endpoint_F13^2_has_169_states_but_split_exponent13;"
        "THM2851_requires_nonsplit_C169=W2(F13)_additive;"
        "cardinality_is_paid_but_extension_class_is_not"
    )
    print(
        "next_test=replace_split_endpoint_translation_by_a_physically_typed_"
        "Witt_carry_law_or_prove_every_lawful_endpoint_rechart_is_affine"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
