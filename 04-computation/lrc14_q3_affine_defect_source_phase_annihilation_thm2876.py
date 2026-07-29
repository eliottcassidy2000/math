#!/usr/bin/env python3
"""Exact companion for THM-2876 q3 affine-defect phase annihilation.

The q3 step-2 and step-68 endpoint masks are nearest under
g(a,b)=(a,2b), with signed defect A3 tensor (delta_12-delta_3).
This probe applies the 26 lawful THM-2868 multiplier samples to that
two-ended defect.  It keeps endpoint-only and full-current pairings
separate, and compares the result with THM-2874's macro-E3 Bockstein seam.

No executable Python assert is used.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import ast
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_digest(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINS = {
    COMP / "lrc14_two_origin_endpoint_projective_kummer_thm2868.py":
        "3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5",
    RESULTS / "lrc14_two_origin_endpoint_projective_kummer_thm2868.out":
        "ce4b8961b3dfa468d51b884b91002872440b05bf70d7e5eecfc3318d4100e2f9",
    COMP / "lrc14_endpoint_kummer_galois_bockstein_thm2874.py":
        "3f15c44dc5f66c660ac3605cc25814adc39594bf193aa130a0f5353d6a6178b0",
    RESULTS / "lrc14_endpoint_kummer_galois_bockstein_thm2874.out":
        "90b993b56508ef3603f94104596b899ed9ec7084a2b58ead1604882873ef5455",
    COMP / "lrc14_endpoint_kummer_hermitian_complement_addendum_thm2874.py":
        "44dfdefbf5392e7840f74e63d190a96a484af71ba9bd31df3ce62a22b827d67e",
    RESULTS / "lrc14_endpoint_kummer_hermitian_complement_addendum_thm2874.out":
        "f914d934c40ef58ea5df0f0df0c61c357c9ab9073db3ba7cbb044d8564886cab",
    COMP / "lrc14_horn_collar_endpoint_carry_thm2859.py":
        "6e062f3cc57c80fcff372c272bc138e280205bb953e484f1cc267340774260f0",
    COMP / "lrc14_horn_collar_prony_typed_descent_gate_thm2859.py":
        "ff9a954e65209d0b96de7d9215ccc6a38dfdbb16245414564a63237924efad28",
}
for path, expected in PINS.items():
    require(
        lf_digest(path) == expected,
        f"pinned dependency changed: {path.name}",
    )


import lrc14_two_origin_endpoint_projective_kummer_thm2868 as atlas
import lrc14_horn_collar_endpoint_carry_thm2859 as hinge


P = 13
A3 = tuple(range(10))
PLUS_TARGET = tuple((a, 12) for a in A3)
MINUS_TARGET = tuple((a, 3) for a in A3)
MINUS_PREIMAGE = tuple((a, 8) for a in A3)
SEAM_SECTION = {3: 0, 11: 8, 7: 4}
SEAM_EDGES = ((3, 8, 11), (11, 9, 7), (3, 4, 7))


def dft(values, omega, prime):
    return tuple(
        sum(
            value * pow(omega, (-frequency * index) % P, prime)
            for index, value in enumerate(values)
        ) % prime
        for frequency in range(P)
    )


def mask_points(mask):
    return frozenset(
        point for point, value in zip(hinge.KEYS, mask) if value
    )


def affine_image(points):
    return frozenset((a, 2 * b % P) for a, b in points)


def split_pair(left_value, right_value, left_node, right_node, prime):
    inverse_difference = pow((left_node - right_node) % prime, -1, prime)
    split_left = (
        right_value - right_node * left_value
    ) * inverse_difference % prime
    split_minus_right = (
        left_node * left_value - right_value
    ) * inverse_difference % prime
    return split_left, split_minus_right


def carried(q, h):
    return (q + h) // P


def main():
    allocation = atlas.allocation
    _module, full, _details, e3, clocks, _q_pairs = (
        allocation.build_geometry()
    )
    period = full.T
    unit = period // P
    delta = 66 * hinge.H
    source_start_base = allocation.ATOM_INTERVAL
    source_end_base = (
        source_start_base[0] + delta,
        source_start_base[1] + delta,
    )
    target_start_base = tuple(
        value + allocation.physical.SHIFT
        for value in source_start_base
    )
    target_end_base = tuple(
        value + allocation.physical.SHIFT
        for value in source_end_base
    )
    source_start = atlas.horn.circular_shift_interval(
        source_start_base, 3 * unit, period
    )
    source_end = atlas.horn.circular_shift_interval(
        source_end_base, 3 * unit, period
    )
    target_start = atlas.horn.circular_shift_interval(
        target_start_base, 3 * unit, period
    )
    target_end = atlas.horn.circular_shift_interval(
        target_end_base, 3 * unit, period
    )
    target_q7 = atlas.horn.circular_shift_interval(
        target_start_base, 7 * unit, period
    )
    require(
        source_end
        == (source_start[0] + delta, source_start[1] + delta)
        and target_end
        == (target_start[0] + delta, target_start[1] + delta),
        "q3 physical intervals stopped differing by 66h",
    )
    require(
        all(
            allocation.contained(interval, e3)
            for interval in (
                source_start,
                source_end,
                target_start,
                target_end,
            )
        )
        and not allocation.contained(target_q7, e3),
        "macro-E3 truth boundary changed",
    )

    start_mask = hinge.endpoint_mask(target_start)
    end_mask = hinge.endpoint_mask(target_end)
    start_description = hinge.cartesian_description(start_mask)
    end_description = hinge.cartesian_description(end_mask)
    require(
        start_description == (
            A3,
            (0, 3, 4, 5, 8, 9, 10, 11, 12),
            90,
            True,
        )
        and end_description == (
            A3,
            (0, 5, 6, 7, 8, 9, 10, 11, 12),
            90,
            True,
        ),
        "q3 endpoint rectangles changed",
    )
    start_points = mask_points(start_mask)
    end_points = mask_points(end_mask)
    transported_start = affine_image(start_points)
    require(
        end_points - transported_start == frozenset(PLUS_TARGET)
        and transported_start - end_points == frozenset(MINUS_TARGET)
        and all((a, 8) in start_points for a in A3)
        and all((a, 3) in transported_start for a in A3),
        "rank-one nearest-affine q3 boundary changed",
    )

    # The raw endpoint restriction is constant on every address whose mask
    # is one, and empty on every address in the Boolean complement.  This
    # makes the same-carrier balanced boundary zero and keeps the global
    # complement as a hostile coarse object.
    start_present = {
        point: start_mask[hinge.KEY_INDEX[point]] for point in hinge.KEYS
    }
    end_present = {
        point: end_mask[hinge.KEY_INDEX[point]] for point in hinge.KEYS
    }
    require(
        all(start_present[point] for point in PLUS_TARGET)
        and all(start_present[point] for point in MINUS_TARGET)
        and all(end_present[point] for point in PLUS_TARGET)
        and not any(end_present[point] for point in MINUS_TARGET)
        and sum(not value for value in start_present.values()) == 79
        and sum(not value for value in end_present.values()) == 79,
        "boundary-address occupancy changed",
    )

    target_start_carrier = ((*target_start, 1),)
    target_end_carrier = ((*target_end, 1),)
    start_positive_restrictions = tuple(
        atlas.indexed_restriction(target_start_carrier, (a, 12))
        for a in A3
    )
    start_negative_preimage_restrictions = tuple(
        atlas.indexed_restriction(target_start_carrier, (a, 8))
        for a in A3
    )
    end_positive_restrictions = tuple(
        atlas.indexed_restriction(target_end_carrier, (a, 12))
        for a in A3
    )
    end_negative_restrictions = tuple(
        atlas.indexed_restriction(target_end_carrier, (a, 3))
        for a in A3
    )
    require(
        all(
            restriction == target_start_carrier
            for restriction in start_positive_restrictions
        )
        and all(
            restriction == target_start_carrier
            for restriction in start_negative_preimage_restrictions
        )
        and all(
            restriction == target_end_carrier
            for restriction in end_positive_restrictions
        )
        and not any(end_negative_restrictions)
        and atlas.indexed_restriction(target_start_carrier, (0, 3))
        == target_start_carrier
        and not atlas.indexed_restriction(
            target_start_carrier, atlas.STEPPED
        ),
        "literal ten-copy endpoint restriction controls changed",
    )
    samples = tuple(
        sample
        for n in atlas.UNIT_SECTIONS
        for sample in (n, n + 1)
    )
    require(
        len(samples) == 26
        and len(set(samples)) == 26
        and all(sample in atlas.FREQUENCY_MEASUREMENTS for sample in samples),
        "THM-2868 lawful 26-sample bank changed",
    )

    phase_audit = hinge.semilinear_character_audit()
    require(
        phase_audit["delta"] == delta
        and phase_audit["characters"] == (7, 6)
        and phase_audit["exact_powers"] == (396, 1254),
        "THM-2859 paired phase audit changed",
    )
    source_phase_rows = {
        row[0]: row for row in phase_audit["specializations"]
    }

    field_rows = []
    for prime, root in atlas.endpoint.MODS:
        xi = pow(root, atlas.endpoint.NN // 2366, prime)
        omega = pow(xi, 182, prime)
        omega3 = pow(omega, 3, prime)
        endpoint_phase = pow(omega, 6, prime)
        source_phase = pow(omega, 7, prime)
        boundary_scalar = 10 * (endpoint_phase - 1) % prime
        conjugate_boundary_scalar = (
            10 * (pow(omega, -6, prime) - 1)
        ) % prime
        boundary_norm = (
            boundary_scalar * conjugate_boundary_scalar % prime
        )
        require(
            endpoint_phase * source_phase % prime == 1
            and boundary_scalar
            == 10 * (omega3 - 1) * (omega3 + 1) % prime
            and boundary_scalar != 0
            and boundary_norm != 0,
            "endpoint/source phase or boundary scalar changed",
        )

        raw_start = {}
        raw_end = {}
        same_start_boundary = {}
        same_end_boundary = {}
        two_ended_boundary = {}
        global_start_complement = {}
        global_end_complement = {}
        for sample in samples:
            frequency = -(12 + 26 * sample)
            start_value = atlas.weighted_endpoint_sum(
                target_start_carrier, frequency, (prime, root)
            )
            end_value = atlas.weighted_endpoint_sum(
                target_end_carrier, frequency, (prime, root)
            )
            require(
                start_value != 0
                and end_value == endpoint_phase * start_value % prime,
                "66h endpoint phase failed on a lawful raw sample",
            )
            raw_start[sample] = start_value
            raw_end[sample] = end_value

            # At the start, both b=12 and b=3 are present and all ten
            # address restrictions equal the same carrier coefficient.
            same_start_boundary[sample] = (
                len(PLUS_TARGET) - len(MINUS_TARGET)
            ) * start_value % prime
            # At the end, b=12 is present and b=3 is absent.
            same_end_boundary[sample] = (
                len(PLUS_TARGET) * end_value
            ) % prime
            # The genuine affine edge pairs end b=12 with the transported
            # start preimage b=8.  Its endpoint-only coefficient is the
            # difference of the two physical interval values.
            two_ended_boundary[sample] = (
                len(PLUS_TARGET) * end_value
                - len(MINUS_PREIMAGE) * start_value
            ) % prime
            global_start_complement[sample] = 0
            global_end_complement[sample] = 0
            require(
                same_start_boundary[sample] == 0
                and same_end_boundary[sample]
                == 10 * end_value % prime
                and two_ended_boundary[sample]
                == boundary_scalar * start_value % prime,
                "address-boundary raw coefficient changed",
            )

        start_left, start_right, start_weight = target_start_carrier[0]
        end_left, end_right, end_weight = target_end_carrier[0]
        require(start_weight == end_weight == 1, "endpoint weight changed")
        left_node = pow(
            root,
            (26 * atlas.endpoint.RDIL * start_left) % atlas.endpoint.NN,
            prime,
        )
        right_node = pow(
            root,
            (26 * atlas.endpoint.RDIL * start_right) % atlas.endpoint.NN,
            prime,
        )
        end_left_node = pow(
            root,
            (26 * atlas.endpoint.RDIL * end_left) % atlas.endpoint.NN,
            prime,
        )
        end_right_node = pow(
            root,
            (26 * atlas.endpoint.RDIL * end_right) % atlas.endpoint.NN,
            prime,
        )
        left_alpha = pow(
            root,
            (12 * atlas.endpoint.RDIL * start_left) % atlas.endpoint.NN,
            prime,
        )
        right_alpha = pow(
            root,
            (12 * atlas.endpoint.RDIL * start_right) % atlas.endpoint.NN,
            prime,
        )
        end_left_alpha = pow(
            root,
            (12 * atlas.endpoint.RDIL * end_left) % atlas.endpoint.NN,
            prime,
        )
        end_right_alpha = pow(
            root,
            (12 * atlas.endpoint.RDIL * end_right) % atlas.endpoint.NN,
            prime,
        )
        require(
            left_node == end_left_node
            and right_node == end_right_node
            and end_left_alpha == endpoint_phase * left_alpha % prime
            and end_right_alpha == endpoint_phase * right_alpha % prime
            and left_node != right_node,
            "66h shift changed a Prony node rather than only the weights",
        )

        start_left_branch = []
        start_right_branch = []
        defect_left_branch = []
        defect_right_branch = []
        for r, (formal, offset, actual) in enumerate(zip(
            atlas.FREQUENCY_SECTIONS,
            atlas.SECTION_OFFSETS,
            atlas.UNIT_SECTIONS,
        )):
            start_split = split_pair(
                raw_start[actual],
                raw_start[actual + 1],
                left_node,
                right_node,
                prime,
            )
            defect_split = split_pair(
                two_ended_boundary[actual],
                two_ended_boundary[actual + 1],
                left_node,
                right_node,
                prime,
            )
            left_transport = pow(
                pow(left_node, offset, prime), -1, prime
            )
            right_transport = pow(
                pow(right_node, offset, prime), -1, prime
            )
            start_left_value = start_split[0] * left_transport % prime
            start_right_value = start_split[1] * right_transport % prime
            defect_left_value = defect_split[0] * left_transport % prime
            defect_right_value = defect_split[1] * right_transport % prime
            require(
                defect_left_value
                == boundary_scalar * start_left_value % prime
                and defect_right_value
                == boundary_scalar * start_right_value % prime
                and start_left_value
                == left_alpha * pow(left_node, formal, prime) % prime
                and start_right_value
                == -right_alpha * pow(right_node, formal, prime) % prime,
                f"local Prony transport failed at section {r}",
            )
            start_left_branch.append(start_left_value)
            start_right_branch.append(start_right_value)
            defect_left_branch.append(defect_left_value)
            defect_right_branch.append(defect_right_value)

        start_left_branch = tuple(start_left_branch)
        start_right_branch = tuple(start_right_branch)
        defect_left_branch = tuple(defect_left_branch)
        defect_right_branch = tuple(defect_right_branch)
        require(
            all(
                defect_left_branch[(r + 1) % P]
                == omega3 * defect_left_branch[r] % prime
                and defect_right_branch[(r + 1) % P]
                == defect_right_branch[r]
                for r in range(P)
            ),
            "defect branch characters changed",
        )

        start_ratio = tuple(
            left * pow(right, -1, prime) % prime
            for left, right in zip(start_left_branch, start_right_branch)
        )
        defect_ratio = tuple(
            left * pow(right, -1, prime) % prime
            for left, right in zip(defect_left_branch, defect_right_branch)
        )
        require(
            defect_ratio == start_ratio
            and all(
                start_ratio[(r + 1) % P]
                == omega3 * start_ratio[r] % prime
                for r in range(P)
            ),
            "projective defect atlas changed",
        )

        normalized_sum = tuple((1 + value) % prime for value in start_ratio)
        normalized_edge = tuple(
            normalized_sum[(r + 1) % P]
            * (1 + pow(start_ratio[r], -1, prime))
            % prime
            for r in range(P)
        )
        conjugate_edge = tuple(
            (1 + pow(start_ratio[(r + 1) % P], -1, prime))
            * normalized_sum[r]
            % prime
            for r in range(P)
        )
        require(
            all(
                normalized_edge[r]
                == omega3 * conjugate_edge[r] % prime
                for r in range(P)
            )
            and all(normalized_edge)
            and len(set(normalized_edge)) == P,
            "defect Hermitian phase or separation failed",
        )
        edge_transform = dft(normalized_edge, omega, prime)
        edge_support = tuple(
            index for index, value in enumerate(edge_transform) if value
        )
        require(
            edge_support == (0, 3, 10),
            "defect Hermitian Fourier support changed",
        )

        # Restoring the actual source coefficient at both physical ends
        # annihilates the affine boundary: source gains omega^7 and target
        # gains omega^6, so every paired full-current sample is invariant.
        phase_row = source_phase_rows[prime]
        (
            _row_prime,
            source_before,
            source_after,
            source_character,
            _target_before,
            _target_after,
            target_character,
        ) = phase_row
        require(
            source_before == atlas.COMMON_SOURCE[prime]
            and source_character == 7
            and target_character == 6
            and source_after
            == source_phase * source_before % prime,
            "actual source coefficient phase changed",
        )
        source_copies_before = tuple(source_before for _a in A3)
        source_copies_after = tuple(source_after for _a in A3)
        require(
            len(set(source_copies_before)) == 1
            and len(set(source_copies_after)) == 1,
            "source P stopped being common over the ten target addresses",
        )
        full_boundary_terms = {}
        full_boundary = []
        for sample in samples:
            frequency = -(12 + 26 * sample)
            paired_terms = tuple(
                (
                    source_copies_after[a]
                    * atlas.weighted_endpoint_sum(
                        end_positive_restrictions[a],
                        frequency,
                        (prime, root),
                    )
                    - source_copies_before[a]
                    * atlas.weighted_endpoint_sum(
                        start_negative_preimage_restrictions[a],
                        frequency,
                        (prime, root),
                    )
                ) % prime
                for a in A3
            )
            full_boundary_terms[sample] = paired_terms
            full_boundary.append(sum(paired_terms) % prime)
        require(
            all(
                not any(paired_terms)
                for paired_terms in full_boundary_terms.values()
            )
            and not any(full_boundary),
            "ten-copy paired full-current affine boundary became nonzero",
        )

        # After localization at 13, the constant endpoint scalar lies on
        # the augmentation line of omega^3-1.  Integrally, the factor 10 is
        # not a cyclotomic unit.  The scalar does not alter any seam
        # transition, so the collapsed frequency seam stays flat.
        scalar_over_bockstein = (
            boundary_scalar * pow((omega3 - 1) % prime, -1, prime)
        ) % prime
        require(
            scalar_over_bockstein == 10 * (omega3 + 1) % prime,
            "augmentation-line comparison changed",
        )
        require(20 % P == 7, "cotangent scalar ratio changed")
        flat_edge = {}
        ancestry_edge = {}
        for source_q, h, target_q in SEAM_EDGES:
            displacement = (
                SEAM_SECTION[target_q] - SEAM_SECTION[source_q]
            ) % P
            flat_edge[source_q, target_q] = pow(
                omega, 3 * displacement, prime
            )
            ancestry_edge[source_q, target_q] = pow(
                omega, 3 * carried(source_q, h), prime
            )
        flat_direct = flat_edge[3, 7]
        flat_via = flat_edge[3, 11] * flat_edge[11, 7] % prime
        ancestry_direct = ancestry_edge[3, 7]
        ancestry_via = (
            ancestry_edge[3, 11] * ancestry_edge[11, 7] % prime
        )
        require(
            flat_via == flat_direct
            and ancestry_via == omega3 * ancestry_direct % prime,
            "THM-2874 flat/Bockstein boundary changed",
        )

        field_rows.append({
            "prime": prime,
            "endpoint_phase": endpoint_phase,
            "source_phase": source_phase,
            "boundary_scalar": boundary_scalar,
            "boundary_norm": boundary_norm,
            "scalar_over_bockstein": scalar_over_bockstein,
            "edge_support": edge_support,
            "edge_distinct": len(set(normalized_edge)),
            "full_boundary_nonzero": sum(value != 0 for value in full_boundary),
            "flat_holonomy": flat_via * pow(flat_direct, -1, prime) % prime,
            "ancestry_holonomy":
                ancestry_via * pow(ancestry_direct, -1, prime) % prime,
        })

    source_tree = ast.parse(Path(__file__).read_text())
    require(
        not any(
            isinstance(node, ast.Assert) for node in ast.walk(source_tree)
        ),
        "executable Python assert statement found",
    )

    print("THM-2876 Q3 AFFINE DEFECT / SOURCE-PHASE ANNIHILATION")
    print(f"pinned={tuple((path.name, digest) for path, digest in PINS.items())}")
    print(
        f"q3_intervals=start={target_start};end={target_end};"
        f"delta={delta}=66h"
    )
    print(
        f"q3_rectangles={start_description}->{end_description};"
        "nearest_affine=(a,b)->(a,2b);"
        f"positive_defect={PLUS_TARGET};negative_defect={MINUS_TARGET};"
        f"negative_preimage={MINUS_PREIMAGE}"
    )
    print(
        "raw_address_boundary=same_start:0,"
        "same_end:10*C_end,"
        "two_ended_endpoint_only:10*(omega^6-1)*C_start;"
        "global_endpoint_complement:0"
    )
    print(
        "lawful_samples=26;local_offsets="
        f"{atlas.SECTION_OFFSETS};Prony_nodes_unchanged_by_66h=1"
    )
    for row in field_rows:
        print(
            f"field={row['prime']};"
            f"endpoint_phase={row['endpoint_phase']};"
            f"source_phase={row['source_phase']};"
            f"K=10*(omega^6-1)={row['boundary_scalar']};"
            f"KbarK={row['boundary_norm']};"
            f"K/(omega^3-1)={row['scalar_over_bockstein']};"
            f"Hermitian_support={row['edge_support']};"
            f"Hermitian_distinct={row['edge_distinct']};"
            f"full_boundary_nonzero={row['full_boundary_nonzero']}/26;"
            f"flat_holonomy={row['flat_holonomy']};"
            f"ancestry_holonomy={row['ancestry_holonomy']}"
        )
    print(
        "ENDPOINT_ONLY_VERDICT=the entire transported Prony atlas is scaled "
        "by the F-rational trivial-character scalar K=10*(omega^6-1); "
        "U remains chi3, V remains trivial, the projective ratio and "
        "normalized THM-2861 edge are unchanged, and the source-retaining "
        "signed coefficient edge scale is multiplied by K*bar(K)"
    )
    print(
        "FULL_CURRENT_VERDICT=zero: the actual source phase omega^7 is "
        "inverse to the endpoint phase omega^6, so the two-ended paired "
        "full-current boundary cancels termwise on all 10 address copies "
        "and all 26 samples; P is common over the 10 target addresses"
    )
    print(
        "BOCKSTEIN_VERDICT=K=10*(omega^3-1)*(omega^3+1) lies on the same "
        "13-localized/cotangent augmentation line, with cotangent ratio 7 "
        "mod13; 10 is not an integral cyclotomic unit, and K is a constant "
        "vertex scalar that leaves frequency seam holonomy one; true "
        "macro-E3 ancestry holonomy remains omega^3"
    )
    print(
        "HOSTILE=the 79-address global endpoint-support complement has "
        "empty restriction on its own carrier and is not the macro-E3 "
        "complement; q3 start/end stay macro-E3 while q7 is not-E3, and "
        "no q7 coefficient or common-support transporter follows from "
        "this affine defect"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
