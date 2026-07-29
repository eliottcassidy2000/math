#!/usr/bin/env python3
"""Exact hostile audit of the THM-2878 event-twisted coefficient atlas.

The positive result is a coefficient-level lift on 169 decorated states.
The negative result is equally important: restoring the two endpoint
origins under an origin-independent common truth-label rule leaves the
signed fixed-source coefficient supported only at q3.

This is exact evidence for THM-2882.  It proves no row exclusion or
LRC(14) result.
"""

from __future__ import annotations

from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    COMP / "lrc14_two_origin_endpoint_projective_kummer_thm2868.py":
        "3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5",
    RESULTS / "lrc14_two_origin_endpoint_projective_kummer_thm2868.out":
        "ce4b8961b3dfa468d51b884b91002872440b05bf70d7e5eecfc3318d4100e2f9",
    COMP / "lrc14_endpoint_kummer_hermitian_complement_addendum_thm2874.py":
        "44dfdefbf5392e7840f74e63d190a96a484af71ba9bd31df3ce62a22b827d67e",
    RESULTS / "lrc14_endpoint_kummer_hermitian_complement_addendum_thm2874.out":
        "f914d934c40ef58ea5df0f0df0c61c357c9ab9073db3ba7cbb044d8564886cab",
    COMP / "lrc14_endpoint_factor_exit_carry_transducer_thm2878.py":
        "b379b9278f6c0d0864908bbc2da2123f4d208eb83c35738d12f651119e7a3366",
    RESULTS / "lrc14_endpoint_factor_exit_carry_transducer_thm2878.out":
        "35bdec6bc5b2bb3c0287bd5aee26c66e8485876e066bf423e2fadb3a94727224",
    COMP / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    RESULTS / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out":
        "424fd2e83a618f862a5ee1b5f073a93fe236d10cdc5412eab1b54dec5e537eac",
    RESULTS / "lrc14_q11_semantic_word_horn_thm2835.out":
        "1ebe0cbaf7d4ef13defed0bdb5b37df1364880acdbfc6139b243ab9df65f6bf6",
    RESULTS / "lrc14_prime_power_unit_mass_q11_response_thm2839.out":
        "495829603ea0c3944f83d7ae269cbbc5cbdec9fdc452395e96c78de8b2e7e24b",
    COMP / "lrc14_q0_q3_one_fibre_selector_provenance_obstruction_thm2880.py":
        "7d379279f08f4df8f16d1fc699f4c6f7a9b657fa151e54d2794ca883b2ceee24",
    RESULTS / "lrc14_q0_q3_one_fibre_selector_provenance_obstruction_thm2880.out":
        "7c7a74d38eac97155b55178d08b317b4658f52cff9651e93113f6b556d19f2e1",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_two_origin_endpoint_projective_kummer_thm2868 as atlas
import lrc14_endpoint_factor_exit_carry_transducer_thm2878 as transducer


allocation = atlas.allocation
endpoint_base = atlas.endpoint_base
endpoint = atlas.endpoint
P = 13
PATTERN_E3 = "SSSSSSSSD"
ZERO_E3_RESIDUES = (0, 3, 11)


def complement(intervals, modulus):
    result = []
    cursor = 0
    for left, right in intervals:
        require(0 <= left < right <= modulus, "bad endpoint interval")
        require(cursor <= left, "endpoint intervals overlap")
        if cursor < left:
            result.append((cursor, left))
        cursor = right
    if cursor < modulus:
        result.append((cursor, modulus))
    return tuple(result)


def restricted_piece(carrier, block):
    return allocation.indexed_weighted_intersection(
        carrier,
        block,
        tuple(left for left, _right in block),
    )


def piece_parameters(piece, embedding):
    prime, root = embedding
    require(len(piece) == 1 and piece[0][2] == 1, "piece is not one atom")
    left, right, _weight = piece[0]
    return (
        pow(root, 12 * endpoint.RDIL * left % endpoint.NN, prime),
        pow(root, 12 * endpoint.RDIL * right % endpoint.NN, prime),
        pow(root, 26 * endpoint.RDIL * left % endpoint.NN, prime),
        pow(root, 26 * endpoint.RDIL * right % endpoint.NN, prime),
    )


def marker_exits(patterns, q, h):
    total = 0
    for offset in range(h):
        source = patterns[(q + offset) % P][3]
        target = patterns[(q + offset + 1) % P][3]
        total += source == "D" and target == "S"
    return total


def main() -> None:
    (
        _module, full_module, _details, _e3, _clocks, _q_pairs
    ) = allocation.build_geometry()
    period = full_module.T
    unit = period // P
    interval = allocation.ATOM_INTERVAL
    target = tuple(value + allocation.physical.SHIFT for value in interval)
    target_atom = ((*target, 1),)
    origins = (atlas.ORIGIN, atlas.STEPPED)

    require(
        period == endpoint_base.T
        and (12 * endpoint.RDIL * unit) % endpoint.NN == 0
        and (26 * endpoint.RDIL * unit) % endpoint.NN == 0,
        "target translation acquired an endpoint phase",
    )

    target_carriers = tuple(
        allocation.physical.overlap.shift_weighted(
            target_atom, q * unit
        )
        for q in range(P)
    )
    require(
        all(
            len(carrier) == 1
            and carrier[0][2] == 1
            and carrier[0][1] - carrier[0][0] == 26444880
            for carrier in target_carriers
        )
        and len({carrier for carrier in target_carriers}) == P,
        "the thirteen translated target atoms changed",
    )

    # The zero endpoint chart has PAT_E3 precisely at q=0,3,11.
    zero_ell = endpoint_base.REPS[atlas.ORIGIN]
    zero_e3 = tuple(endpoint.build_set(endpoint_base.PAT_E3, zero_ell))
    zero_not_e3 = complement(zero_e3, endpoint_base.T)
    zero_pieces = {}
    zero_truth = []
    zero_patterns = []
    for q, carrier in enumerate(target_carriers):
        e3_piece = restricted_piece(carrier, zero_e3)
        complement_piece = restricted_piece(carrier, zero_not_e3)
        pattern = transducer.literal_pattern(
            transducer.TARGET_ATOM, q, zero_ell
        )
        zero_patterns.append(pattern)
        if e3_piece == carrier and not complement_piece:
            truth = "E3"
            chosen = e3_piece
            require(pattern == PATTERN_E3, "PAT_E3 literal type changed")
        else:
            truth = "not-E3"
            chosen = complement_piece
            require(
                chosen == carrier and not e3_piece and pattern != PATTERN_E3,
                "full-pattern complement typing changed",
            )
        zero_truth.append(truth)
        zero_pieces[q] = chosen
    zero_truth = tuple(zero_truth)
    zero_patterns = tuple(zero_patterns)
    require(
        zero_ell == (0,) * 9
        and tuple(q for q, truth in enumerate(zero_truth) if truth == "E3")
        == ZERO_E3_RESIDUES,
        "zero-chart truth atlas changed",
    )

    # Independently derive kappa from the unique positive u3 D->S event.
    event_candidates = []
    for factor in range(9):
        for source_state, target_state in (("S", "D"), ("D", "S")):
            good = True
            for q in range(P):
                for h in range(P):
                    count = 0
                    for offset in range(h):
                        source = zero_patterns[(q + offset) % P][factor]
                        target_state_seen = zero_patterns[
                            (q + offset + 1) % P
                        ][factor]
                        count += (
                            source == source_state
                            and target_state_seen == target_state
                        )
                    if count != (q + h) // P:
                        good = False
            if good:
                event_candidates.append(
                    (factor, f"{source_state}_to_{target_state}")
                )
    require(
        event_candidates == [(3, "D_to_S")]
        and all(
            marker_exits(zero_patterns, q, h) == (q + h) // P
            for q in range(P) for h in range(P)
        ),
        "THM-2878 event carry changed",
    )

    formal_sections = atlas.FREQUENCY_SECTIONS
    offsets = atlas.SECTION_OFFSETS
    actual_sections = atlas.UNIT_SECTIONS
    sample_pairs = tuple((value, value + 1) for value in actual_sections)
    all_actual_samples = tuple(
        value for pair in sample_pairs for value in pair
    )
    frequency_pairs = tuple(
        (12 + 26 * left, 12 + 26 * right)
        for left, right in sample_pairs
    )
    require(
        formal_sections == tuple(1 + 42 * r for r in range(P))
        and offsets == (0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0)
        and all(
            gcd(value, 91) == 1
            for pair in sample_pairs for value in pair
        )
        and all(
            gcd(value, 91) == 1
            for pair in frequency_pairs for value in pair
        ),
        "the 26-sample 91-unit atlas changed",
    )

    field_rows = []
    for embedding in endpoint.MODS:
        prime, root = embedding
        xi = pow(root, endpoint.NN // 2366, prime)
        omega = pow(xi, 182, prime)
        source_value = atlas.COMMON_SOURCE[prime]
        require(
            pow(omega, P, prime) == 1
            and omega != 1
            and source_value != 0,
            "field normalization changed",
        )

        parameters = tuple(
            piece_parameters(zero_pieces[q], embedding) for q in range(P)
        )
        require(
            len(set(parameters)) == 1,
            "truth-polarized target atoms lost their common Prony nodes",
        )
        alpha_left, alpha_right, lambda_left, lambda_right = parameters[0]
        require(
            lambda_left == pow(xi, 13, prime)
            and lambda_right == pow(xi, 169, prime)
            and lambda_left != lambda_right
            and pow(lambda_left, 42, prime) == pow(omega, 3, prime)
            and pow(lambda_right, 42, prime) == 1,
            "Prony node character types changed",
        )
        inverse_difference = pow(
            (lambda_left - lambda_right) % prime, -1, prime
        )

        split_u = {}
        split_v = {}
        raw_pairs = {}
        for q in range(P):
            piece = zero_pieces[q]
            for r, (formal, offset, actual) in enumerate(zip(
                formal_sections, offsets, actual_sections
            )):
                raw = []
                for sample in (actual, actual + 1):
                    frequency = -(12 + 26 * sample)
                    value = allocation.endpoint_sum(
                        piece, frequency, embedding
                    )
                    require(value != 0, "truth-polarized sample vanished")
                    raw.append(source_value * value % prime)
                current, current_next = raw
                left_at_actual = (
                    current_next - lambda_right * current
                ) * inverse_difference % prime
                right_at_actual = (
                    lambda_left * current - current_next
                ) * inverse_difference % prime
                left = (
                    left_at_actual
                    * pow(pow(lambda_left, offset, prime), -1, prime)
                ) % prime
                right = (
                    right_at_actual
                    * pow(pow(lambda_right, offset, prime), -1, prime)
                ) % prime
                require(
                    left
                    == source_value * alpha_left
                    * pow(lambda_left, formal, prime) % prime
                    and right
                    == -source_value * alpha_right
                    * pow(lambda_right, formal, prime) % prime,
                    "Prony transport changed",
                )
                split_u[q, r] = left
                split_v[q, r] = right
                raw_pairs[q, r] = tuple(raw)

        require(
            all(
                split_u[q, r] == split_u[0, r]
                and split_v[q, r] == split_v[0, r]
                and raw_pairs[q, r] == raw_pairs[0, r]
                for q in range(P) for r in range(P)
            )
            and len({split_u[0, r] for r in range(P)}) == P
            and len({split_v[0, r] for r in range(P)}) == 1
            and all(
                split_u[0, (r + 1) % P]
                == pow(omega, 3, prime) * split_u[0, r] % prime
                for r in range(P)
            ),
            "the all-q split atlas is not q-flat/chi3-by-section",
        )

        # Put the collapsed section r=q-3 over q, then restore the missing
        # ancestry coordinate a.  The lifted state (a,q) uses section
        #
        #                 r(a,q)=a+q-3 mod13.
        #
        # It is distinguishable by (target support, split-left value).
        section = {
            (a, q): (a + q - 3) % P
            for a in range(P) for q in range(P)
        }
        state_u = {
            (a, q): split_u[q, section[a, q]]
            for a in range(P) for q in range(P)
        }
        state_v = {
            (a, q): split_v[q, section[a, q]]
            for a in range(P) for q in range(P)
        }
        state_s = {
            state: (state_u[state] + state_v[state]) % prime
            for state in state_u
        }
        decorated_states = {
            (zero_pieces[q], state_u[a, q])
            for a in range(P) for q in range(P)
        }
        require(
            len(section) == P**2
            and len(decorated_states) == P**2,
            "support plus chi3 coefficient stopped distinguishing C169",
        )

        edge_checks = 0
        lifted_state_checks = 0
        gauge_checks = 0
        recombined_scalar_failures = 0
        for q in range(P):
            for h in range(P):
                carry = marker_exits(zero_patterns, q, h)
                flat = pow(omega, 3 * h, prime)
                event = pow(omega, 3 * carry, prime)
                twisted = event * flat % prime
                q_next = (q + h) % P
                a_next = carry % P
                ratio = (
                    state_u[a_next, q_next]
                    * pow(state_u[0, q], -1, prime)
                ) % prime
                require(
                    ratio == twisted,
                    "event-twisted coefficient edge changed",
                )
                require(
                    state_v[a_next, q_next] == state_v[0, q],
                    "trivial Prony channel acquired event phase",
                )
                recombined_scalar_failures += (
                    state_s[a_next, q_next]
                    != twisted * state_s[0, q] % prime
                )
                edge_checks += 1

                # Uniformity over all initial ancestry states.
                for a in range(P):
                    destination = ((a + carry) % P, q_next)
                    state_ratio = (
                        state_u[destination]
                        * pow(state_u[a, q], -1, prime)
                    ) % prime
                    require(
                        state_ratio == twisted,
                        "lifted C169 state transport changed",
                    )
                    lifted_state_checks += 1

                # Gauge away the flat frequency section.  The remainder is
                # exactly the THM-2878 event/ancestry character.
                phi_source = pow(omega, -3 * (q - 3), prime)
                phi_target = pow(omega, -3 * (q_next - 3), prime)
                require(
                    event * phi_source % prime
                    == phi_target * twisted % prime,
                    "flat-section gauge did not leave the event character",
                )
                gauge_checks += 1

        composition_checks = 0
        full_state_composition_checks = 0
        for q in range(P):
            for h in range(P):
                for k in range(P):
                    q_after_h = (q + h) % P
                    q_reduced = (q + h + k) % P
                    carry_h = marker_exits(zero_patterns, q, h)
                    carry_k = marker_exits(
                        zero_patterns, q_after_h, k
                    )
                    reduced_step = (h + k) % P
                    carry_reduced = marker_exits(
                        zero_patterns, q, reduced_step
                    )
                    central = (h + k) // P
                    require(
                        carry_h + carry_k == carry_reduced + central,
                        "event cocycle identity changed",
                    )
                    edge_h = (
                        pow(omega, 3 * carry_h, prime)
                        * pow(omega, 3 * h, prime)
                    ) % prime
                    edge_k = (
                        pow(omega, 3 * carry_k, prime)
                        * pow(omega, 3 * k, prime)
                    ) % prime
                    edge_reduced = (
                        pow(omega, 3 * carry_reduced, prime)
                        * pow(omega, 3 * reduced_step, prime)
                    ) % prime
                    require(
                        edge_k * edge_h % prime
                        == pow(omega, 3 * central, prime)
                        * edge_reduced % prime,
                        "twisted edge composition changed",
                    )
                    composition_checks += 1

                    # One representative a suffices algebraically, but audit
                    # all thirteen total-state lifts as a hostile control.
                    for a in range(P):
                        after_h = (
                            (a + carry_h) % P, q_after_h
                        )
                        after_hk = (
                            (after_h[0] + carry_k) % P,
                            q_reduced,
                        )
                        direct = (
                            (a + carry_reduced + central) % P,
                            q_reduced,
                        )
                        require(
                            after_hk == direct,
                            "C169 lifted-state composition changed",
                        )
                        full_state_composition_checks += 1

        require(
            edge_checks == P**2
            and lifted_state_checks == P**3
            and gauge_checks == P**2
            and recombined_scalar_failures == 144
            and composition_checks == P**3
            and full_state_composition_checks == P**4,
            "transport/composition census changed",
        )

        def edge_data(q, h):
            carry = marker_exits(zero_patterns, q, h)
            flat = pow(omega, 3 * h, prime)
            event = pow(omega, 3 * carry, prime)
            return carry, flat, event, flat * event % prime

        edge_3_11 = edge_data(3, 8)
        edge_11_7 = edge_data(11, 9)
        edge_3_7 = edge_data(3, 4)
        flat_holonomy = (
            edge_11_7[1] * edge_3_11[1]
            * pow(edge_3_7[1], -1, prime)
        ) % prime
        event_holonomy = (
            edge_11_7[2] * edge_3_11[2]
            * pow(edge_3_7[2], -1, prime)
        ) % prime
        twisted_holonomy = (
            edge_11_7[3] * edge_3_11[3]
            * pow(edge_3_7[3], -1, prime)
        ) % prime
        require(
            tuple(edge[0] for edge in (
                edge_3_11, edge_11_7, edge_3_7
            )) == (0, 1, 0)
            and flat_holonomy == 1
            and event_holonomy == twisted_holonomy
            == pow(omega, 3, prime)
            and twisted_holonomy * pow(omega, -3, prime) % prime == 1,
            "q3/q11/q7 holonomy changed",
        )

        # Restore the actual two endpoint-origin selector.  Use the same
        # truth label E3/not-E3 at both origins.  The zero chart dictates
        # that label.  The signed origin difference remains delta_q3.
        origin_blocks = {}
        origin_piece = {}
        origin_truth_piece = {}
        addresswise_full_truth = {}
        for origin_index, address in enumerate(origins):
            ell = endpoint_base.REPS[address]
            e3_block = tuple(
                endpoint.build_set(endpoint_base.PAT_E3, ell)
            )
            not_e3_block = complement(e3_block, endpoint_base.T)
            origin_blocks[origin_index] = {
                "E3": e3_block, "not-E3": not_e3_block
            }
            full_truth = []
            for q, carrier in enumerate(target_carriers):
                pieces = {
                    truth: restricted_piece(carrier, block)
                    for truth, block in origin_blocks[origin_index].items()
                }
                for truth, piece in pieces.items():
                    origin_truth_piece[origin_index, q, truth] = piece
                origin_piece[origin_index, q] = pieces[zero_truth[q]]
                containing = tuple(
                    truth for truth, piece in pieces.items()
                    if piece == carrier
                )
                require(
                    len(containing) == 1,
                    "originwise Boolean selector ceased to be total",
                )
                full_truth.append(containing[0])
            addresswise_full_truth[origin_index] = tuple(full_truth)

        require(
            addresswise_full_truth[0] == zero_truth
            and tuple(
                q for q, truth in enumerate(addresswise_full_truth[1])
                if truth == "E3"
            ) == (0, 11)
            and addresswise_full_truth[0][3] == "E3"
            and addresswise_full_truth[1][3] == "not-E3",
            "originwise selector-typing dichotomy changed",
        )

        signed_support_by_section = []
        q11_origin_values = []
        q7_origin_values = []
        for actual in all_actual_samples:
            signed_row = []
            for q in range(P):
                values = []
                for origin_index in range(2):
                    value = allocation.endpoint_sum(
                        origin_piece[origin_index, q],
                        -(12 + 26 * actual),
                        embedding,
                    )
                    values.append(value)
                signed_row.append((values[0] - values[1]) % prime)
                if q == 11:
                    q11_origin_values.append(tuple(values))
                if q == 7:
                    q7_origin_values.append(tuple(values))
            signed_support_by_section.append(
                tuple(q for q, value in enumerate(signed_row) if value)
            )
        require(
            tuple(signed_support_by_section)
            == ((3,),) * len(all_actual_samples)
            and all(left == right != 0 for left, right in q11_origin_values)
            and all(left == right != 0 for left, right in q7_origin_values),
            "origin-restored truth current is no longer delta_q3",
        )

        # If each origin instead chooses whichever truth block contains its
        # atom, both origins are full and equal at every q.  This produces
        # the 13-q coefficient atlas only by making the selector
        # origin-dependent at q3; the signed current then vanishes globally.
        addresswise_signed_support = []
        for actual in all_actual_samples:
            signed = []
            for q in range(P):
                values = []
                for origin_index in range(2):
                    truth = addresswise_full_truth[origin_index][q]
                    piece = restricted_piece(
                        target_carriers[q],
                        origin_blocks[origin_index][truth],
                    )
                    values.append(allocation.endpoint_sum(
                        piece,
                        -(12 + 26 * actual),
                        embedding,
                    ))
                signed.append((values[0] - values[1]) % prime)
            addresswise_signed_support.append(
                tuple(q for q, value in enumerate(signed) if value)
            )
        require(
            tuple(addresswise_signed_support)
            == ((),) * len(all_actual_samples),
            "origin-dependent total selector stopped cancelling identically",
        )

        # Classify all four Boolean selector pairs at every q.  Encode
        # E3=1 and complement=0.  If t_o is the containing truth and s_o
        # the selected truth, the origin contributes iff s_o=t_o.  Thus
        # signed support is nonzero exactly when
        #
        #       s_0 XOR s_1 XOR t_0 XOR t_1 = 1.
        #
        # The truth-origin parity t_0 XOR t_1 is delta_q3.  In particular,
        # q7/q11 support requires an origin-odd selector.
        bit_truth = {0: "not-E3", 1: "E3"}
        truth_bit = {"not-E3": 0, "E3": 1}
        truth_origin_parity = tuple(
            truth_bit[addresswise_full_truth[0][q]]
            ^ truth_bit[addresswise_full_truth[1][q]]
            for q in range(P)
        )
        require(
            truth_origin_parity
            == tuple(int(q == 3) for q in range(P)),
            "truth-origin parity stopped being delta_q3",
        )
        selector_classification_checks = 0
        for actual in all_actual_samples:
            for q in range(P):
                base_value = allocation.endpoint_sum(
                    zero_pieces[q], -(12 + 26 * actual), embedding
                )
                require(base_value != 0, "base endpoint value vanished")
                t0 = truth_bit[addresswise_full_truth[0][q]]
                t1 = truth_bit[addresswise_full_truth[1][q]]
                for s0 in range(2):
                    for s1 in range(2):
                        values = tuple(
                            allocation.endpoint_sum(
                                origin_truth_piece[
                                    origin_index, q,
                                    bit_truth[(s0, s1)[origin_index]]
                                ],
                                -(12 + 26 * actual),
                                embedding,
                            )
                            for origin_index in range(2)
                        )
                        f0 = int(s0 == t0)
                        f1 = int(s1 == t1)
                        signed = (values[0] - values[1]) % prime
                        expected_signed = (f0 - f1) * base_value % prime
                        support_parity = s0 ^ s1 ^ t0 ^ t1
                        require(
                            signed == expected_signed
                            and bool(signed) == bool(support_parity),
                            "Boolean selector parity law changed",
                        )
                        selector_classification_checks += 1
        require(
            selector_classification_checks == 2 * P**2 * 4,
            "selector classification census changed",
        )

        # There is a unique selector producing the positive zero-origin
        # coefficient at each q: retain the containing zero-origin block
        # and select the empty block at the stepped origin.  This makes the
        # all-q atlas a precise origin-polarized coefficient difference,
        # but it is not the inherited common physical selector.
        positive_selector = tuple(
            (
                addresswise_full_truth[0][q],
                bit_truth[
                    1 ^ truth_bit[addresswise_full_truth[1][q]]
                ],
            )
            for q in range(P)
        )
        require(
            (positive_selector[3], positive_selector[11],
             positive_selector[7])
            == (
                ("E3", "E3"),
                ("E3", "not-E3"),
                ("not-E3", "E3"),
            ),
            "positive seam selector changed",
        )
        positive_selector_parity = tuple(
            truth_bit[left] ^ truth_bit[right]
            for left, right in positive_selector
        )
        selector_toggle_edges = tuple(
            q for q in range(P)
            if positive_selector_parity[q]
            != positive_selector_parity[(q + 1) % P]
        )
        require(
            positive_selector_parity
            == tuple(1 ^ int(q == 3) for q in range(P))
            and selector_toggle_edges == (2, 3)
            and positive_selector_parity[12]
            == positive_selector_parity[0] == 1
            and marker_exits(zero_patterns, 12, 1) == 1
            and positive_selector_parity[0]
            != positive_selector_parity[3]
            and marker_exits(zero_patterns, 0, 3) == 0,
            "selector parity and carry event ceased to be independent",
        )
        positive_selector_checks = 0
        for actual in all_actual_samples:
            for q in range(P):
                selected = positive_selector[q]
                values = tuple(
                    allocation.endpoint_sum(
                        origin_truth_piece[
                            origin_index, q, selected[origin_index]
                        ],
                        -(12 + 26 * actual),
                        embedding,
                    )
                    for origin_index in range(2)
                )
                base_value = allocation.endpoint_sum(
                    zero_pieces[q], -(12 + 26 * actual), embedding
                )
                require(
                    values[0] == base_value
                    and values[1] == 0,
                    "positive origin polarization changed",
                )
                positive_selector_checks += 1
        require(
            positive_selector_checks == 2 * P**2,
            "positive selector census changed",
        )

        # A fixed-section target translation is coefficient-flat.  The
        # nontrivial edge E comes from changing frequency section according
        # to r(a,q), not from translating the source-retaining current.
        fixed_section_ratios = tuple(
            split_u[(q + h) % P, 0]
            * pow(split_u[q, 0], -1, prime) % prime
            for q in range(P) for h in range(P)
        )
        require(
            fixed_section_ratios == (1,) * P**2,
            "raw fixed-section q-translation acquired a phase",
        )

        field_rows.append((
            prime,
            parameters[0][2:],
            len(decorated_states),
            edge_3_11,
            edge_11_7,
            edge_3_7,
            flat_holonomy,
            event_holonomy,
            twisted_holonomy,
            split_u[0, 0],
            split_v[0, 0],
            tuple(signed_support_by_section),
            recombined_scalar_failures,
            selector_classification_checks,
        ))

    word_output = (
        RESULTS / "lrc14_q11_semantic_word_horn_thm2835.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    spectrum_output = (
        RESULTS / "lrc14_prime_power_unit_mass_q11_response_thm2839.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    bockstein_output = (
        RESULTS / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    require(
        "natural_beta_lift=(a,11)->(a+1,7); "
        "word_state=QAB_on_449/449; QA=QB=0" in word_output
        and "status=VERIFIED-EXACT COMPANION; coefficient/support theorem "
        "only; no current, row decrement, or LRC14" in word_output
        and "marked_gauge_vertical_clutch=" in spectrum_output
        and "verdict=word_gauge_retains_H_to_C_but_endpoint_gauge_moves_"
        "source_to_all_safe_no_word" in spectrum_output
        and "semantic_leg=THM2835_PINNED:"
        "449_QA(q11,a)_to_QAB(q7,a+1)" in bockstein_output,
        "pinned QA/QAB support/current boundary changed",
    )

    print("THM-2882 EVENT-TWISTED ALL-q COEFFICIENT C169 LIFT AUDIT")
    print(
        f"zero_truth={zero_truth}; E3_residues={ZERO_E3_RESIDUES}; "
        f"literal_patterns={zero_patterns}"
    )
    print(
        "all_13_zero_chart_atoms=one_full_translated_weight1_atom; "
        "PAT_E3_iff_literal_pattern_SSSSSSSSD; "
        "full_complement_otherwise"
    )
    print(
        f"formal_sections={formal_sections}; offsets={offsets}; "
        f"sample_pairs={sample_pairs}; frequency_pairs={frequency_pairs}; "
        "all_52_entries_are_91_units=1"
    )
    print(
        "event=u3_D_to_S_at_q12_to_q0; "
        f"unique_oriented_event={event_candidates}; "
        "kappa(q,h)=floor((q+h)/13)"
    )
    print(
        "C169_states=n=13a+q_with_base13_set_chart_(a,q)_not_product_group; "
        "section_r=a+q-3_mod13; "
        "vertex=(truth_atom_Jq,U_r); 169_distinct_support_coefficient_pairs"
    )
    print(
        "transport=L_h(a,q)=(a+kappa(q,h),q+h_mod13); "
        "flat_F=omega^(3h); event_tau=omega^(3kappa); "
        "twisted_E=tau*F"
    )
    print(
        "transport_checks=(base_edges=169,lifted_states=2197,"
        "gauge=169); composition_checks=(projective=2197,"
        "all_C169_states=28561); "
        "law=E(q+h_mod13,k)E(q,h)=omega^(3floor((h+k)/13))"
        "E(q,h+k_mod13)"
    )
    for row in field_rows:
        print(
            f"field={row[0]}; common_nodes={row[1]}; "
            f"decorated_states={row[2]}; "
            f"edge_q3_q11={row[3]}; edge_q11_q7={row[4]}; "
            f"edge_q3_q7={row[5]}; flat_seam_curvature={row[6]}; "
            f"event_seam_curvature={row[7]}; "
            f"twisted_seam_curvature={row[8]}; "
            "lifted_loop_closed_by_T^-1_has_holonomy=1; "
            f"U0={row[9]}; V0={row[10]}; "
            f"recombined_scalar_transport_failures={row[12]}/169"
        )
    print(
        "gauge=phi(q)=omega^(-3(q-3)); "
        "tau(q,h)phi(q)=phi(q+h_mod13)E(q,h); "
        "event-twisted system_is_gauge_equivalent_to_ancestry_character"
    )
    print(
        "selector_parity_law=nonzero_iff_"
        "s0_XOR_s1_XOR_t0_XOR_t1=1; "
        "truth_origin_parity=t0_XOR_t1=delta_q3; "
        "q7/q11_nonzero_requires_origin-odd_selector"
    )
    print(
        "origin-independent_same-truth-label_selector: "
        "origin_supports=(zero_E3={0,3,11},stepped_E3={0,11}); "
        "signed_current_support={3}_for_all_26_samples; "
        "q11_E3_origins_equal_nonzero_then_cancel; "
        "q7_complement_origins_equal_nonzero_then_cancel"
    )
    print(
        "addresswise_full_atom_selector: q3_requires_E3_at_origin0_but_"
        "not-E3_at_origin12; signed_current_support=empty; "
        "selector_is_origin_dependent"
    )
    print(
        "unique_zero-origin-positive_selector_within_Boolean_"
        "E3/complement_class_on_seam="
        "(q3:(E3,E3),q11:(E3,not-E3),q7:(not-E3,E3)); "
        "all_q_signed_coefficient_equals_origin0_atom; "
        "this_closes_coefficient_support_only_by_projecting_away_origin12"
    )
    print(
        "selector_vs_carry_independence=selector_parity="
        "1_XOR_delta_q3,toggles_only_at_q2_to_q3_and_q3_to_q4; "
        "u3_carry_q12_to_q0_has_selector_parity_constant_1; "
        "q0_to_q3_changes_selector_parity_but_kappa(0,3)=0; "
        "event_twist_cannot_supply_origin-odd_selector"
    )
    print(
        "physicality_hostile=fixed_frequency_section_target_translation_"
        "ratio_is_1_on_all_169_edges; nontrivial_E_comes_from_external_"
        "frequency-section_reindex_plus_event_decoration"
    )
    print(
        "channel_scope=rank-one_charged_split_U_or_centered_c-A_only; "
        "two-node_pair_transport=diag(E,1); recombined_S=U+V_is_not_"
        "scalar-E-transported_on_144/169_edges"
    )
    print(
        "QA_QAB_gate=449-sheet_(a,11)->(a+1,7)_support_is_pinned, "
        "but no coefficient-to-Lambda basepoint/alignment and marked "
        "endpoint gauge moves the fixed source to all-safe/no-word; "
        "inherited_common-selector_q11/q7_signed_coefficients_are_zero"
    )
    print(
        "VERDICT=exact_coefficient-level_C169_lift_and_exact_omega3_"
        "reduced_Bockstein_seam_curvature/vertical_T_clutch; "
        "honest_lifted_loop_is_flat; not_a_lawful_operation_on_the_actual_"
        "source-retaining_current; no_row_or_LRC14_conclusion"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
