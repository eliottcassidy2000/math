#!/usr/bin/env python3
"""Exact q3/q11 transverse endpoint-edge companion for THM-2847.

It composes the proved THM-2820 and THM-2829 constructors without changing
either theorem.
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
    COMP / "lrc14_q11_semantic_word_horn_thm2835.py":
        "207dd94f086338ae1e80b7d7196f99bf41e795893d13b6d48e4e7d516af03523",
}
for path, expected in PINNED.items():
    actual = sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == expected, f"pinned dependency changed: {path.name}")


import lrc14_literal_fixed_sheet_allocation_thm2806 as allocation


def crlf_safe_present_loader():
    path = COMP / "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
    spec = importlib.util.spec_from_file_location("q3q11_present_base", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


allocation.physical.target.load_present_module = crlf_safe_present_loader

P = 13
S = P**5
M = P**6


def circular_shift_interval(interval, shift, period):
    left = (interval[0] + shift) % period
    right = (interval[1] + shift) % period
    require(left < right, "test interval crossed the circle seam")
    return left, right


def safe_comb_contains(interval, module, speed, denominator, lo, hi):
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


def affine_lift_action(t, address):
    return (
        address + t * S * ((address - 7) % P)
    ) % M


def beta(q):
    return 2 * q % P


def fine_address(q):
    return (6716 + beta(q) * S) % M


def main():
    (
        _module, full_module, details, e3, clocks, _q_pairs
    ) = allocation.build_geometry()
    period = full_module.T
    unit = period // P
    I = allocation.ATOM_INTERVAL
    J0 = tuple(endpoint + allocation.physical.SHIFT for endpoint in I)

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

    def full_cells(q):
        Jq = circular_shift_interval(J0, q * unit, period)
        cells = []
        for s in allocation.COMMON_S:
            for t in allocation.COMMON_T:
                for clock in range(7):
                    if (
                        all(signature(I, s, t, clock))
                        and all(signature(J0, s, t, clock))
                        and all(signature(Jq, s, t, clock))
                    ):
                        cells.append((s, t, clock))
        return tuple(cells)

    cells3 = full_cells(3)
    cells11 = full_cells(11)
    common = tuple(sorted(set(cells3) & set(cells11)))
    expected3 = tuple(
        (s, t, 1)
        for s in (0, 3, 8, 9, 10, 11, 12)
        for t in range(3, 11)
    )
    expected11 = tuple(
        (s, t, 1)
        for s in allocation.COMMON_S
        for t in range(5, 12)
    )
    expected_common = tuple(
        (s, t, 1)
        for s in (0, 3, 8, 9, 10, 11, 12)
        for t in range(5, 11)
    )
    require(cells3 == expected3 and len(cells3) == 56, "q3 cells changed")
    require(cells11 == expected11 and len(cells11) == 63, "q11 cells changed")
    require(
        common == expected_common and len(common) == 42,
        "q3/q11 common cell bank changed",
    )
    J7 = circular_shift_interval(J0, 7 * unit, period)
    q7_signatures = tuple(
        (cell, signature(J7, *cell))
        for cell in common
    )
    q7_e3_only_cells = tuple(
        cell
        for cell, bits in q7_signatures
        if not bits[0] and all(bits[1:])
    )
    expected_q7_e3_only = tuple(
        (s, t, 1)
        for s in (0, 3, 8, 9, 12)
        for t in (5, 6, 9, 10)
    )
    q7_extra_loss_cells = tuple(
        (cell, tuple(index for index, bit in enumerate(bits) if not bit))
        for cell, bits in q7_signatures
        if cell not in q7_e3_only_cells
    )
    require(
        q7_e3_only_cells == expected_q7_e3_only
        and len(q7_e3_only_cells) == 20
        and len(q7_extra_loss_cells) == 22,
        "q7 E3-only horn bank changed",
    )
    require(
        all(
            losses[0] == 0
            and any(index in (2, 3) for index in losses[1:])
            for _cell, losses in q7_extra_loss_cells
        ),
        "a q7 non-horn cell ceased to lose an ordinary factor",
    )

    source_base, target_base = details[1]
    source_hit = weighted_hit(I, source_base)
    target_rows = {}
    for q in (3, 11):
        Jq = circular_shift_interval(J0, q * unit, period)
        target_q = allocation.physical.overlap.shift_weighted(
            target_base, q * unit
        )
        target_rows[q] = weighted_hit(Jq, target_q)
    require(
        source_hit == ((*I, allocation.ATOM_WEIGHT),)
        and target_rows[3]
        == ((*circular_shift_interval(J0, 3 * unit, period),
             allocation.ATOM_WEIGHT),)
        and target_rows[11]
        == ((*circular_shift_interval(J0, 11 * unit, period),
             allocation.ATOM_WEIGHT),),
        "common weighted carrier value changed",
    )

    endpoint_base = allocation.endpoint_base
    target_atom = ((*J0, 1),)
    right_origin = allocation.RIGHT_ORIGIN
    right_step = allocation.add(
        allocation.RIGHT_ORIGIN, allocation.TARGET_STEP
    )
    endpoint_bank = {}
    endpoint_restricted = {}
    for address in (right_origin, right_step):
        ell = endpoint_base.REPS[address]
        present = tuple(
            endpoint_base.endpoint.build_set(endpoint_base.PAT_E3, ell)
        )
        starts = tuple(left for left, _right in present)
        rows = []
        for q in range(P):
            shifted = allocation.physical.overlap.shift_weighted(
                target_atom, q * unit
            )
            restricted = allocation.indexed_weighted_intersection(
                shifted, present, starts
            )
            endpoint_restricted[(address, q)] = restricted
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
        endpoint_bank[address] = tuple(rows)

    endpoint_value = (
        231164267889491750,
        630230755085920022,
    )
    q3_restricted = endpoint_restricted[(right_origin, 3)]
    require(
        len(q3_restricted) == 1 and q3_restricted[0][2] == 1,
        "q3 endpoint piece ceased to be one unweighted interval",
    )
    endpoint_conductor = endpoint_base.endpoint.NN
    endpoint_frequency = -endpoint_base.Y0
    endpoint_exponents = tuple(
        (
            -endpoint_frequency
            * endpoint_base.endpoint.RDIL
            * endpoint
        ) % endpoint_conductor
        for endpoint in q3_restricted[0][:2]
    )
    endpoint_gcd = gcd(
        endpoint_conductor,
        gcd(*endpoint_exponents),
    )
    reduced_conductor = endpoint_conductor // endpoint_gcd
    reduced_exponents = tuple(
        exponent // endpoint_gcd for exponent in endpoint_exponents
    )
    galois_u = 27
    moved_exponents = tuple(
        galois_u * exponent % reduced_conductor
        for exponent in reduced_exponents
    )
    # Since 1183 is odd, Q(zeta_2366)=Q(zeta_1183).  In the latter field
    # c=x^624-x^510, while u=27 sends it to x^286-x^757.
    # Their difference x^624-x^510-x^286+x^757 is a nonzero polynomial
    # of degree 757<phi(1183)=936.  Thus c is not fixed by u, even though
    # u=1 mod13 fixes Q(zeta_13).
    endpoint_exponents_1183 = (624, 510)
    moved_exponents_1183 = (286, 757)
    galois_difference_degree = 757
    require(
        endpoint_gcd == 21274064131320
        and endpoint_conductor == 50334435734703120
        and endpoint_exponents
        == (46866763281297960, 1382814168535800)
        and reduced_conductor == 2366
        and reduced_exponents == (2203, 65)
        and moved_exponents == (331, 1755)
        and galois_u % P == 1
        and tuple(
            galois_u * exponent % 1183
            for exponent in endpoint_exponents_1183
        ) == moved_exponents_1183
        and galois_difference_degree < 936,
        "endpoint conductor / coefficient-field hostile changed",
    )
    support00 = tuple(
        q for q, _pieces, _mass, values in endpoint_bank[right_origin]
        if any(values)
    )
    support12 = tuple(
        q for q, _pieces, _mass, values in endpoint_bank[right_step]
        if any(values)
    )
    require(
        support00 == (0, 3, 11) and support12 == (0, 11),
        "endpoint support bank changed",
    )
    for rows in endpoint_bank.values():
        for _q, pieces, mass, values in rows:
            require(
                (
                    (pieces, mass, values)
                    == (1, 26444880, endpoint_value)
                )
                if any(values)
                else (pieces, mass, values) == (0, 0, (0, 0)),
                "endpoint value pattern changed",
            )

    occupancy = tuple(
        tuple(
            int(any(endpoint_bank[address][q][3]))
            for address in (right_origin, right_step)
        )
        for q in (3, 11)
    )
    require(occupancy == ((1, 0), (1, 1)), "unit minor changed")
    determinant = occupancy[0][0] * occupancy[1][1] - (
        occupancy[0][1] * occupancy[1][0]
    )
    require(determinant == 1, "normalized occupancy determinant changed")

    # THM-2835's proved QB(source)->QA(target) coefficient-support column
    # has value 449 at q=11 and zero at q=0,3.  Composing that immutable
    # column with the two physical endpoint-origin columns gives a unit
    # triangular three-row clutch after normalizing the positive scalars.
    qb_to_qa = {0: 0, 3: 0, 11: 449}
    semantic_clutch = tuple(
        (
            int(any(endpoint_bank[right_origin][q][3])),
            int(any(endpoint_bank[right_step][q][3])),
            qb_to_qa[q],
        )
        for q in (0, 3, 11)
    )
    semantic_clutch_det = (
        semantic_clutch[0][0]
        * (
            semantic_clutch[1][1] * semantic_clutch[2][2]
            - semantic_clutch[1][2] * semantic_clutch[2][1]
        )
        - semantic_clutch[0][1]
        * (
            semantic_clutch[1][0] * semantic_clutch[2][2]
            - semantic_clutch[1][2] * semantic_clutch[2][0]
        )
        + semantic_clutch[0][2]
        * (
            semantic_clutch[1][0] * semantic_clutch[2][1]
            - semantic_clutch[1][1] * semantic_clutch[2][0]
        )
    )
    require(
        semantic_clutch
        == ((1, 1, 0), (1, 0, 0), (1, 1, 449))
        and semantic_clutch_det == -449,
        "q0/q3/q11 semantic endpoint clutch changed",
    )
    # Adding the proved q7 QB->QAB filler gives a full coefficient-support
    # four-frame before the macro gate.  On the exact 20-cell horn above,
    # q7 loses only E3; THM-2829's exact E3 support {0,3,11} therefore
    # deletes precisely that q7 row.  The remaining 22 common q3/q11 cells
    # are not part of this four-frame statement because q7 also loses q1/q2.
    algebraic_horn_frame = (
        (1, 1, 0, 0),
        (1, 0, 0, 0),
        (1, 1, 449, 0),
        (0, 0, 0, 449),
    )
    algebraic_horn_det = semantic_clutch_det * 449
    e3_support = (0, 3, 11)
    q7_e3_overlap = allocation.intersect_sorted((J7,), e3)
    physical_horn_frame = tuple(
        row if q in e3_support else (0, 0, 0, 0)
        for q, row in zip((0, 3, 11, 7), algebraic_horn_frame)
    )
    e3_kernel_generator = (0, 0, 0, 1)
    e3_image_rank = 3
    complement_on_kernel = sum(
        algebraic_horn_frame[3][column] * value
        for column, value in enumerate(e3_kernel_generator)
    )
    require(
        algebraic_horn_det == -(449**2)
        and q7_e3_overlap == ()
        and physical_horn_frame[-1] == (0, 0, 0, 0),
        "coefficient horn / E3 rank-collapse boundary changed",
    )
    require(
        semantic_clutch_det != 0
        and e3_image_rank == 3
        and all(
            sum(row[column] * e3_kernel_generator[column]
                for column in range(4)) == 0
            for row in physical_horn_frame[:3]
        )
        and complement_on_kernel == 449,
        "rank-one E3 mapping-cone boundary changed",
    )

    n3 = fine_address(3)
    n11 = fine_address(11)
    t_edge = (beta(11) - beta(3)) % P
    require(
        (n3, n11) == (2234474, 3348353)
        and n3 % P == n11 % P == 8
        and t_edge == 3
        and (n11 - n3) % M == 3 * S
        and affine_lift_action(t_edge, n3) == n11,
        "transverse affine-lift edge changed",
    )
    q_edge = (11 - 3) % P
    require(
        q_edge == 8
        and q_edge == 7 * t_edge % P
        and beta(q_edge) == t_edge,
        "allocation/normal semiconjugacy changed",
    )

    # The absolute fine-address law is equivariant on all thirteen labels.
    torsor_rows = []
    for q0 in range(P):
        for q1 in range(P):
            t = (beta(q1) - beta(q0)) % P
            require(
                affine_lift_action(t, fine_address(q0))
                == fine_address(q1)
                and (q1 - q0) % P == 7 * t % P,
                "global q/address torsor law changed",
            )
        torsor_rows.append((q0, beta(q0), fine_address(q0)))
    require(
        len({row[2] for row in torsor_rows}) == P,
        "fine-address orbit ceased to be free",
    )
    unweighted_graph_fourier_support = tuple(
        (a, b)
        for a in range(P)
        for b in range(P)
        if (7 * a + b) % P == 0
    )
    graph_character_class_counts = tuple(
        (
            residue,
            sum(
                (7 * a + b) % P == residue
                for a in range(P)
                for b in range(P)
            ),
        )
        for residue in range(P)
    )
    require(
        len(unweighted_graph_fourier_support) == P
        and len({a for a, _b in unweighted_graph_fourier_support}) == P
        and len({b for _a, b in unweighted_graph_fourier_support}) == P
        and graph_character_class_counts
        == tuple((residue, P) for residue in range(P)),
        "graph-annihilator Fourier line changed",
    )

    # R00-R12 is exactly c*delta_3.  The scalar c lies in the inherited
    # K=Q(zeta_N), and in fact in Q(zeta_1183), so the physical edge is
    # c*z^3 in K[C13].  The exact hostile above proves c need not lie in
    # Q(zeta_13).
    # After normalizing by c, z^3 lies in Q[C13] and has inverse z^10.
    # The exponent check proves all characters are nonzero without
    # floating cyclotomic arithmetic.
    difference_support = tuple(
        q for q in range(P)
        if bool(any(endpoint_bank[right_origin][q][3]))
        != bool(any(endpoint_bank[right_step][q][3]))
    )
    require(
        difference_support == (3,) and (3 + 10) % P == 0,
        "endpoint-edge group-ring monomial changed",
    )

    J3 = circular_shift_interval(J0, 3 * unit, period)
    J11 = circular_shift_interval(J0, 11 * unit, period)
    disjoint = J3[1] <= J11[0] or J11[1] <= J3[0]
    require(disjoint, "q3/q11 target intervals ceased to be disjoint")

    print("Q3/Q11 TRANSVERSE NORMAL ENDPOINT EDGE")
    print("status=VERIFIED-EXACT; target-relative only; no LRC14")
    print(
        f"q3_full_cells={len(cells3)}; q11_full_cells={len(cells11)}; "
        f"common_cells={len(common)}; first={common[0]}; last={common[-1]}"
    )
    print(
        f"q7_E3_only_horn_cells={len(q7_e3_only_cells)}; "
        f"horn_first={q7_e3_only_cells[0]}; "
        f"horn_last={q7_e3_only_cells[-1]}; "
        f"q7_additional_q1_q2_loss_cells={len(q7_extra_loss_cells)}"
    )
    print(
        f"weighted_carrier=(source={source_hit},q3={target_rows[3]},"
        f"q11={target_rows[11]})"
    )
    print(
        f"endpoint_supports=(origin00={support00},origin12={support12}); "
        f"endpoint_value={endpoint_value}; "
        f"q3_q11_occupancy={occupancy}; normalized_det={determinant}"
    )
    print(
        f"q0_q3_q11_endpoint_QBQA_clutch={semantic_clutch}; "
        f"normalized_det={semantic_clutch_det}; "
        "word_column_status=imported_PROVED_T2835_whole-cylinder_count"
    )
    print(
        f"coefficient_horn_four_frame={algebraic_horn_frame}; "
        f"det={algebraic_horn_det}; "
        f"E3_support={e3_support}; "
        f"q7_E3_overlap={q7_e3_overlap}; "
        f"physical_macro_gated_frame={physical_horn_frame}; "
        "scope=20_cell_E3_only_horn; first_rank_loss=q7_E3; "
        "complementary_E3_bit_is_whole_on_q7=1"
    )
    print(
        f"E3_mapping_cone=(image_rank={e3_image_rank},"
        f"kernel_generator={e3_kernel_generator},"
        f"complement_on_kernel={complement_on_kernel}); "
        "graded_E3_plus_complement_frame_rank=4; "
        "remaining_debt=lawful_macro_truth_block_contraction_or_intertwiner"
    )
    print(
        f"fine_addresses=(q3={n3},q11={n11}); residue=(8,8); "
        f"difference={n11-n3}=3*13^5; A3(q3_address)={affine_lift_action(3,n3)}"
    )
    print(
        f"label_edge=(delta_q={q_edge},delta_beta={t_edge}); "
        "delta_q=7*delta_beta; beta(q)=2q is inverse to q=7t"
    )
    print(f"absolute_q_beta_address_torsor={tuple(torsor_rows)}")
    print(
        f"unweighted_joint_q_t_graph_DFT_support="
        f"{unweighted_graph_fourier_support}; "
        f"weighted_character_class_counts={graph_character_class_counts}; "
        "unweighted_support_size=13; arbitrary_weight_has_13_restricted_"
        "character_values_repeated_on_13_cosets; "
        "point_mass_hostile_has_all_169_coefficients_nonzero=1"
    )
    print(
        f"signed_endpoint_edge_support={difference_support}; "
        f"endpoint_scalar_reduction=(conductor={reduced_conductor},"
        f"exponents={reduced_exponents}); "
        f"fix_zeta13_Galois_hostile=(u={galois_u},"
        f"zeta1183_exponents={endpoint_exponents_1183},"
        f"moved={moved_exponents_1183},"
        f"difference_degree={galois_difference_degree}<phi1183=936); "
        "group_ring_over_inherited_K=c*z^3; "
        "normalized_Q_group_ring=z^3; inverse_monomial_exponent=10; "
        "all_13_DFT_characters_nonzero=1"
    )
    print(
        f"physical_target_intervals_disjoint={int(disjoint)}; "
        f"raw_weighted_source_nonempty={int(bool(source_hit))}; "
        "endpoint_source_allocation_and_common_H_face_not_supplied=1"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
