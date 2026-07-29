#!/usr/bin/env python3
"""Exact q0-to-q3 one-fibre completion and provenance audit.

This scratch companion starts from the promoted THM-2877 injection

    g(a,b) = (a+1,b+8),
    E_q3 = g(E_q0) disjoint-union ({9} x B3),

and keeps three operations rigorously separate:

1. endpoint-address support/reindexing;
2. cyclotomic coefficient and Prony transport; and
3. physical interval, source, and ancestry transport.

It tests the canonical two-origin selector, its pushforward and pullback,
the positive and character-three polarizations of the missing fibre, all
26 lawful THM-2868 raw samples, the THM-2874 Hermitian edge, macro E3,
the QA/QAB columns, and the literal U/V ancestry gate.

No executable Python assert statement is used.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path
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
    COMP / "lrc14_semilinear_endpoint_rectangle_classification_thm2877.py":
        "10a4c965d02f9fab60f135d0bf10184096eeb70d30d0d9ef7ad4fcf5fc1aa447",
    RESULTS / "lrc14_semilinear_endpoint_rectangle_classification_thm2877.out":
        "9ab7526b1e268176a98f0cf21238bc4cbabad0ce30c28b6565d0ca7b9d9a751e",
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
    RESULTS / "lrc14_horn_collar_prony_typed_descent_gate_thm2859.out":
        "487363c8d8d34cf703dd83fa6d3867b5932796454792a6394da44354bd59278b",
    COMP / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    RESULTS / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.out":
        "155fce129c750a9505fdda3c71a250ff3a57fcd4044bb1df941da83c08baee1d",
}
for path, expected in PINS.items():
    require(
        lf_digest(path) == expected,
        f"pinned dependency changed: {path.name}",
    )


import lrc14_two_origin_endpoint_projective_kummer_thm2868 as atlas
import lrc14_horn_collar_endpoint_carry_thm2859 as endpoint
import lrc14_horn_collar_prony_typed_descent_gate_thm2859 as typed_gate


P = 13
A0 = frozenset((0, 1, 2, 3, 4, 5, 6, 7, 12))
B0 = frozenset((0, 1, 2, 3, 4, 5, 8, 9, 10))
A3 = frozenset(range(10))
B3 = frozenset((0, 3, 4, 5, 8, 9, 10, 11, 12))
ORIGINS = (atlas.ORIGIN, atlas.STEPPED)
# These are marked right/target-endpoint selector addresses.  They are not
# the physical left/source origin carrying COMMON_SOURCE.


def g(point):
    a, b = point
    return ((a + 1) % P, (b + 8) % P)


def g_inverse(point):
    a, b = point
    return ((a - 1) % P, (b - 8) % P)


def points(mask):
    return frozenset(
        point for point, value in zip(endpoint.KEYS, mask) if value
    )


def moment(support, frequency, omega, prime):
    first, second = frequency
    return sum(
        pow(
            omega,
            (-first * a - second * b) % P,
            prime,
        )
        for a, b in support
    ) % prime


def split_pair(left_value, right_value, left_node, right_node, prime):
    inverse_difference = pow(
        (left_node - right_node) % prime, -1, prime
    )
    return (
        (right_value - right_node * left_value)
        * inverse_difference % prime,
        (left_node * left_value - right_value)
        * inverse_difference % prime,
    )


def dft(values, omega, prime):
    return tuple(
        sum(
            value * pow(omega, (-frequency * index) % P, prime)
            for index, value in enumerate(values)
        ) % prime
        for frequency in range(P)
    )


def sequence_digest(values):
    payload = ",".join(str(value) for value in values).encode()
    return sha256(payload).hexdigest()


def main():
    allocation = atlas.allocation
    _module, full, details, e3, clocks, _q_pairs = (
        allocation.build_geometry()
    )
    period = full.T
    unit = period // P
    source_q0 = allocation.ATOM_INTERVAL
    source_q3 = atlas.horn.circular_shift_interval(
        source_q0, 3 * unit, period
    )
    target_q0 = tuple(
        value + allocation.physical.SHIFT for value in source_q0
    )
    target_q3 = atlas.horn.circular_shift_interval(
        target_q0, 3 * unit, period
    )
    source_carrier = ((*source_q0, 1),)
    target_q0_carrier = ((*target_q0, 1),)
    target_q3_carrier = ((*target_q3, 1),)

    source_mask = endpoint.endpoint_mask(source_q0)
    source_q3_mask = endpoint.endpoint_mask(source_q3)
    target_q0_mask = endpoint.endpoint_mask(target_q0)
    target_mask = endpoint.endpoint_mask(target_q3)
    source_support = points(source_mask)
    target_support = points(target_mask)
    image_support = frozenset(g(point) for point in source_support)
    fibre = target_support - image_support
    preimage_fibre = frozenset(g_inverse(point) for point in fibre)
    expected_source = frozenset(
        (a, b) for a in A0 for b in B0
    )
    expected_target = frozenset(
        (a, b) for a in A3 for b in B3
    )
    expected_fibre = frozenset((9, b) for b in B3)
    expected_preimage_fibre = frozenset((8, b) for b in B0)
    require(
        source_support == expected_source
        and points(target_q0_mask) == expected_source
        and points(source_q3_mask) == expected_target
        and target_support == expected_target
        and image_support.isdisjoint(fibre)
        and image_support | fibre == target_support
        and fibre == expected_fibre
        and preimage_fibre == expected_preimage_fibre
        and len(image_support) == 81
        and len(fibre) == 9,
        "q0-to-q3 one-fibre decomposition changed",
    )

    moved_origins = tuple(g(origin) for origin in ORIGINS)
    pulled_origins = tuple(g_inverse(origin) for origin in ORIGINS)
    require(
        ORIGINS == ((0, 0), (12, 0))
        and moved_origins == ((1, 8), (0, 8))
        and pulled_origins == ((12, 5), (11, 5))
        and tuple(point in source_support for point in ORIGINS) == (True, True)
        and tuple(point in target_support for point in ORIGINS) == (True, False)
        and tuple(point in target_support for point in moved_origins)
        == (True, True)
        and tuple(point in source_support for point in pulled_origins)
        == (True, False)
        and fibre.isdisjoint(ORIGINS)
        and fibre.isdisjoint(moved_origins),
        "named-origin pushforward/pullback dichotomy changed",
    )

    # The direct evaluator has already proved that every address mask value
    # is constant on the entire tiny carrier.  Replay the selected charts
    # with the interval constructor itself, including all nine missing
    # target points and all nine empty source preimages.
    require(
        atlas.indexed_restriction(source_carrier, atlas.ORIGIN)
        == source_carrier
        and atlas.indexed_restriction(source_carrier, pulled_origins[0])
        == source_carrier
        and not atlas.indexed_restriction(
            source_carrier, pulled_origins[1]
        )
        and all(
            not atlas.indexed_restriction(source_carrier, address)
            for address in preimage_fibre
        )
        and all(
            atlas.indexed_restriction(target_q3_carrier, address)
            == target_q3_carrier
            for address in fibre
        )
        and all(
            atlas.indexed_restriction(target_q3_carrier, address)
            == target_q3_carrier
            for address in moved_origins
        )
        and atlas.indexed_restriction(
            target_q3_carrier, ORIGINS[0]
        ) == target_q3_carrier
        and not atlas.indexed_restriction(
            target_q3_carrier, ORIGINS[1]
        ),
        "literal source/target chart restrictions changed",
    )

    # The physical allocation carrier exists and all native q3 factors
    # survive, but this is interval data rather than an address action.
    q3_cells = atlas.full_cells(
        3,
        source_q0,
        target_q0,
        unit,
        period,
        full,
        e3,
        clocks,
    )
    source_base, target_base = details[1]
    target_q3_weighted = allocation.physical.overlap.shift_weighted(
        target_base, 3 * unit
    )
    source_hit = atlas.horn.weighted_hit(source_q0, source_base)
    target_hit = atlas.horn.weighted_hit(target_q3, target_q3_weighted)
    require(
        len(q3_cells) == 56
        and source_hit
        == ((*source_q0, allocation.ATOM_WEIGHT),)
        and target_hit
        == ((*target_q3, allocation.ATOM_WEIGHT),)
        and all(
            allocation.contained(interval, e3)
            for interval in (
                source_q0,
                target_q0,
                source_q3,
                target_q3,
            )
        ),
        "physical carrier, q3 cell bank, or macro-E3 truth changed",
    )

    # Cheap cross-test against the directed factor-u3 exit marker.  At the
    # zero representative, u3 is dangerous only at target residue q12, so
    # its unique D->S unit exit is q12->q0.  The present q0->q3 step does
    # not wrap and sees no event.  Moreover every endpoint representative
    # has ell_u3=0, so neither g nor any selector/fibre choice can alter it.
    u3_index = 3
    u3_speed = atlas.endpoint.W[u3_index]
    u3_danger = full.make_comb(u3_speed, 182, -13, 13)
    u3_states = []
    for residue in range(P):
        interval = atlas.horn.circular_shift_interval(
            target_q0, residue * unit, period
        )
        hit = allocation.intersect_sorted((interval,), u3_danger)
        if allocation.contained(interval, u3_danger):
            u3_states.append("D")
        elif not hit:
            u3_states.append("S")
        else:
            raise RuntimeError("u3 is partial on a selected target atom")
    u3_states = tuple(u3_states)

    def u3_exits(start, step):
        return sum(
            u3_states[(start + offset) % P] == "D"
            and u3_states[(start + offset + 1) % P] == "S"
            for offset in range(step)
        )

    require(
        tuple(
            residue for residue, state in enumerate(u3_states)
            if state == "D"
        ) == (12,)
        and u3_exits(0, 3) == 0
        and u3_exits(3, 10) == 1
        and u3_exits(11, 9) == 1
        and all(
            representative[u3_index] == 0
            for representative in atlas.endpoint_base.REPS.values()
        ),
        "factor-u3 directed exit marker changed",
    )

    # Recompute the literal ancestry at q0 and q3 in both physical frames.
    # The q0 root has a nonempty U factor; q3 has none.  Address relabelling
    # cannot affect these interval-midpoint contributor sets.
    ancestry = endpoint.independent.ancestry
    ancestry_args = typed_gate.ancestry_arguments()
    ancestry_bank = {}
    for q, source_interval, target_interval in (
        (0, source_q0, target_q0),
        (3, source_q3, target_q3),
    ):
        for frame, interval in (
            ("source", source_interval),
            ("target", target_interval),
        ):
            ancestry_bank[q, frame] = ancestry.contributor_sets(
                sum(interval) // 2, *ancestry_args
            )
    require(
        ancestry_bank[0, "source"] == ancestry_bank[0, "target"]
        and ancestry_bank[3, "source"] == ancestry_bank[3, "target"],
        "source/target ancestry equality within a q-frame changed",
    )
    q0_ancestry = ancestry_bank[0, "source"]
    q3_ancestry = ancestry_bank[3, "source"]
    q0_ancestry_row = (
        tuple(map(len, q0_ancestry)),
        ancestry.path_digest(*q0_ancestry),
    )
    q3_ancestry_row = (
        tuple(map(len, q3_ancestry)),
        ancestry.path_digest(*q3_ancestry),
    )
    require(
        q0_ancestry_row
        == (
            (966606, 28534),
            "15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd",
        )
        and q3_ancestry_row
        == (
            (0, 28535),
            "bdf95bb641f8d0ab0d8ebe764fecda782b9a7fe4995031c7f23e2f02868de931",
        )
        and not q3_ancestry[0],
        "q0/q3 ancestry gate changed",
    )

    # Pin the separately typed semantic columns.  The endpoint-origin
    # boundary q0=(1,1), q3=(1,0) is present, while QA and QAB are zero in
    # both rows.  Thus the one-fibre address completion has no semantic
    # ancestry leg to q11/q7.
    horn_output = (
        RESULTS / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.out"
    ).read_text()
    typed_output = (
        RESULTS / "lrc14_horn_collar_prony_typed_descent_gate_thm2859.out"
    ).read_text()
    require(
        "q0_q3_q11_endpoint_QBQA_clutch="
        "((1, 1, 0), (1, 0, 0), (1, 1, 449))" in horn_output
        and "coefficient_horn_four_frame="
        "((1, 1, 0, 0), (1, 0, 0, 0), "
        "(1, 1, 449, 0), (0, 0, 0, 449))" in horn_output
        and "FIRST_STRICT_LOSS=physical same-labelled-cell forest vertex;"
        in typed_output
        and "q3/q11 pulled intervals are not forest vertices and have empty "
        "U ancestry" in typed_output,
        "QA/QAB or strict physical type-gate certificate changed",
    )

    # A q-shift by 3 units is literally invisible to both endpoint Prony
    # weights and nodes.  Hence q0 and q3 use the same raw bank; no Galois
    # or cyclotomic coefficient rechart is required.
    q_displacement = 3 * unit
    require(
        (12 * atlas.endpoint.RDIL * q_displacement)
        % atlas.endpoint.NN == 0
        and (26 * atlas.endpoint.RDIL * q_displacement)
        % atlas.endpoint.NN == 0,
        "q0-to-q3 displacement acquired an endpoint coefficient phase",
    )
    samples = tuple(
        sample
        for actual in atlas.UNIT_SECTIONS
        for sample in (actual, actual + 1)
    )
    require(
        len(samples) == 26
        and len(set(samples)) == 26
        and all(sample in atlas.FREQUENCY_MEASUREMENTS for sample in samples),
        "lawful 26-sample bank changed",
    )

    field_rows = []
    strict_missing_fibre_term_checks = 0
    for prime, root in atlas.endpoint.MODS:
        xi = pow(root, atlas.endpoint.NN // 2366, prime)
        omega = pow(xi, 182, prime)
        omega3 = pow(omega, 3, prime)
        address_phase = pow(omega, 2, prime)
        address_phase_inverse = pow(address_phase, -1, prime)
        common_source = atlas.COMMON_SOURCE[prime]
        require(common_source != 0, "common source coefficient vanished")

        source_character = moment(
            source_support, (0, 3), omega, prime
        )
        image_character = moment(
            image_support, (0, 3), omega, prime
        )
        fibre_character = moment(fibre, (0, 3), omega, prime)
        target_character = moment(
            target_support, (0, 3), omega, prime
        )
        require(
            source_character != 0
            and fibre_character != 0
            and image_character
            == address_phase * source_character % prime
            and image_character == 9 * fibre_character % prime
            and target_character == 10 * fibre_character % prime
            and target_character
            == (image_character + fibre_character) % prime,
            "address character-three phase or fibre decomposition changed",
        )

        raw_q0 = {}
        raw_q3 = {}
        canonical_full = {}
        pushed_full = {}
        pulled_full = {}
        strict_positive_fibre = {}
        strict_chi3_fibre = {}
        free_positive_fibre = {}
        free_chi3_fibre = {}
        for sample in samples:
            frequency = -(12 + 26 * sample)
            value_q0 = atlas.weighted_endpoint_sum(
                target_q0_carrier, frequency, (prime, root)
            )
            value_q3 = atlas.weighted_endpoint_sum(
                target_q3_carrier, frequency, (prime, root)
            )
            require(
                value_q0 == value_q3 != 0,
                "q0/q3 endpoint numerator bank diverged",
            )
            raw_q0[sample] = value_q0
            raw_q3[sample] = value_q3

            # Expand the source-induced missing-fibre expression pointwise,
            # rather than relying only on its already proved empty preimage.
            for target_address in fibre:
                source_address = g_inverse(target_address)
                source_coefficient = (
                    common_source
                    if source_address in source_support
                    else 0
                )
                require(
                    source_coefficient * value_q3 % prime == 0,
                    "missing-fibre source term became nonzero",
                )
                strict_missing_fibre_term_checks += 1

            # The target selector pulls back to one present and one absent
            # q0 chart, so it reproduces the canonical q3 current.  The
            # source selector pushes forward to two present q3 charts and
            # cancels.  The missing fibre meets neither selector.
            canonical_full[sample] = common_source * value_q3 % prime
            pulled_full[sample] = common_source * value_q0 % prime
            pushed_full[sample] = 0

            # On the actual missing fibre every source preimage is absent.
            # Therefore every paired full-current term is zero, independently
            # of the chosen target weights.  Attaching a free duplicate of P
            # would give the displayed nonzero endpoint-only alternatives.
            strict_positive_fibre[sample] = 0
            strict_chi3_fibre[sample] = 0
            free_positive_fibre[sample] = (
                9 * common_source * value_q3 % prime
            )
            free_chi3_fibre[sample] = (
                fibre_character * common_source * value_q3 % prime
            )

        require(
            canonical_full == pulled_full
            and all(canonical_full.values())
            and not any(pushed_full.values())
            and not any(strict_positive_fibre.values())
            and not any(strict_chi3_fibre.values())
            and all(free_positive_fibre.values())
            and all(free_chi3_fibre.values()),
            "named-selector or missing-fibre current dichotomy changed",
        )

        q0_left, q0_right, q0_weight = target_q0_carrier[0]
        q3_left, q3_right, q3_weight = target_q3_carrier[0]
        q0_left_node = pow(
            root,
            (26 * atlas.endpoint.RDIL * q0_left) % atlas.endpoint.NN,
            prime,
        )
        q0_right_node = pow(
            root,
            (26 * atlas.endpoint.RDIL * q0_right) % atlas.endpoint.NN,
            prime,
        )
        q3_left_node = pow(
            root,
            (26 * atlas.endpoint.RDIL * q3_left) % atlas.endpoint.NN,
            prime,
        )
        q3_right_node = pow(
            root,
            (26 * atlas.endpoint.RDIL * q3_right) % atlas.endpoint.NN,
            prime,
        )
        q0_left_alpha = pow(
            root,
            (12 * atlas.endpoint.RDIL * q0_left) % atlas.endpoint.NN,
            prime,
        )
        q0_right_alpha = pow(
            root,
            (12 * atlas.endpoint.RDIL * q0_right) % atlas.endpoint.NN,
            prime,
        )
        q3_left_alpha = pow(
            root,
            (12 * atlas.endpoint.RDIL * q3_left) % atlas.endpoint.NN,
            prime,
        )
        q3_right_alpha = pow(
            root,
            (12 * atlas.endpoint.RDIL * q3_right) % atlas.endpoint.NN,
            prime,
        )
        require(
            q0_weight == q3_weight == 1
            and q0_left_node == q3_left_node
            and q0_right_node == q3_right_node
            and q0_left_alpha == q3_left_alpha
            and q0_right_alpha == q3_right_alpha
            and q0_left_node != q0_right_node,
            "q0/q3 Prony nodes or endpoint weights changed",
        )

        pulled_left = []
        pulled_right = []
        free_fibre_left = []
        free_fibre_right = []
        for r, (formal, offset, actual) in enumerate(zip(
            atlas.FREQUENCY_SECTIONS,
            atlas.SECTION_OFFSETS,
            atlas.UNIT_SECTIONS,
        )):
            split = split_pair(
                pulled_full[actual],
                pulled_full[actual + 1],
                q0_left_node,
                q0_right_node,
                prime,
            )
            free_split = split_pair(
                free_chi3_fibre[actual],
                free_chi3_fibre[actual + 1],
                q0_left_node,
                q0_right_node,
                prime,
            )
            left_transport = pow(
                pow(q0_left_node, offset, prime), -1, prime
            )
            right_transport = pow(
                pow(q0_right_node, offset, prime), -1, prime
            )
            left = split[0] * left_transport % prime
            right = split[1] * right_transport % prime
            fibre_left = free_split[0] * left_transport % prime
            fibre_right = free_split[1] * right_transport % prime
            require(
                left
                == common_source
                * q0_left_alpha
                * pow(q0_left_node, formal, prime)
                % prime
                and right
                == -common_source
                * q0_right_alpha
                * pow(q0_right_node, formal, prime)
                % prime
                and fibre_left == fibre_character * left % prime
                and fibre_right == fibre_character * right % prime,
                f"pulled/free-fibre Prony split failed at section {r}",
            )
            pulled_left.append(left)
            pulled_right.append(right)
            free_fibre_left.append(fibre_left)
            free_fibre_right.append(fibre_right)

        pulled_left = tuple(pulled_left)
        pulled_right = tuple(pulled_right)
        free_fibre_left = tuple(free_fibre_left)
        free_fibre_right = tuple(free_fibre_right)
        require(
            all(
                pulled_left[(r + 1) % P]
                == omega3 * pulled_left[r] % prime
                and pulled_right[(r + 1) % P] == pulled_right[r]
                and free_fibre_left[(r + 1) % P]
                == omega3 * free_fibre_left[r] % prime
                and free_fibre_right[(r + 1) % P]
                == free_fibre_right[r]
                for r in range(P)
            ),
            "pulled/free-fibre character-three atlas changed",
        )

        ratio = tuple(
            left * pow(right, -1, prime) % prime
            for left, right in zip(pulled_left, pulled_right)
        )
        free_fibre_ratio = tuple(
            left * pow(right, -1, prime) % prime
            for left, right in zip(free_fibre_left, free_fibre_right)
        )
        require(
            ratio == free_fibre_ratio
            and all(
                ratio[(r + 1) % P] == omega3 * ratio[r] % prime
                for r in range(P)
            ),
            "projective Kummer ratio changed",
        )
        normalized_sum = tuple((1 + value) % prime for value in ratio)
        hermitian_edge = tuple(
            normalized_sum[(r + 1) % P]
            * (1 + pow(ratio[r], -1, prime))
            % prime
            for r in range(P)
        )
        reverse_edge = tuple(
            (1 + pow(ratio[(r + 1) % P], -1, prime))
            * normalized_sum[r]
            % prime
            for r in range(P)
        )
        hermitian_transform = dft(hermitian_edge, omega, prime)
        hermitian_support = tuple(
            index
            for index, value in enumerate(hermitian_transform)
            if value
        )
        require(
            all(
                hermitian_edge[r] == omega3 * reverse_edge[r] % prime
                for r in range(P)
            )
            and all(hermitian_edge)
            and len(set(hermitian_edge)) == P
            and hermitian_support == (0, 3, 10),
            "THM-2874 Hermitian edge replay changed",
        )

        # Normalizing the address-character phase is an optional scalar
        # gauge, not a coefficient automorphism.  It scales U and V equally,
        # leaves their characters and ratio unchanged, and cannot turn a
        # zero pushed selector or zero source fibre into a nonzero current.
        gauged_left = tuple(
            address_phase_inverse * value % prime
            for value in pulled_left
        )
        gauged_right = tuple(
            address_phase_inverse * value % prime
            for value in pulled_right
        )
        gauged_ratio = tuple(
            left * pow(right, -1, prime) % prime
            for left, right in zip(gauged_left, gauged_right)
        )
        require(
            gauged_ratio == ratio
            and all(
                gauged_left[(r + 1) % P]
                == omega3 * gauged_left[r] % prime
                and gauged_right[(r + 1) % P] == gauged_right[r]
                for r in range(P)
            ),
            "address-phase scalar gauge changed the coefficient characters",
        )

        field_rows.append({
            "prime": prime,
            "raw_digest": sequence_digest(
                tuple(raw_q3[sample] for sample in samples)
            ),
            "C_first": raw_q3[samples[0]],
            "P": common_source,
            "address_phase": address_phase,
            "fibre_character": fibre_character,
            "U0": pulled_left[0],
            "V0": pulled_right[0],
            "hermitian_support": hermitian_support,
            "free_chi3_fibre_U0": free_fibre_left[0],
        })

    require(
        strict_missing_fibre_term_checks
        == len(atlas.endpoint.MODS) * len(fibre) * len(samples),
        "missing-fibre pointwise term census changed",
    )

    source_tree = ast.parse(Path(__file__).read_text())
    require(
        not any(
            isinstance(node, ast.Assert)
            for node in ast.walk(source_tree)
        ),
        "executable Python assert statement found",
    )

    print("Q0->Q3 ONE-FIBRE COMPLETION / PROVENANCE AUDIT")
    print(
        f"physical_intervals=q0_source={source_q0};q0_target={target_q0};"
        f"q3_source={source_q3};q3_target={target_q3};"
        f"q3_cells={len(q3_cells)};weight={allocation.ATOM_WEIGHT};"
        "macro_E3=(1,1,1,1)"
    )
    print(
        "address_exact=E_q3=g(E_q0) disjoint_union F;"
        f"g=(+1,+8);sizes=(81,9,90);F={tuple(sorted(fibre))};"
        f"g_inverse(F)={tuple(sorted(preimage_fibre))}"
    )
    print(
        f"named_right_endpoint_selector=q0={ORIGINS};"
        f"pushforward={moved_origins};q3_pullback={pulled_origins};"
        "occupancy=source(1,1),target(1,0),pushforward_target(1,1),"
        "pullback_source(1,0);F_hits_neither_selector=1"
    )
    print(
        f"ancestry=q0{q0_ancestry_row};q3{q3_ancestry_row};"
        "common_U=0;QA_QAB_rows=q0(0,0),q3(0,0);"
        "strict_same_labelled_forest_q3=0"
    )
    print(
        f"u3_exit_marker=states={''.join(u3_states)};"
        "danger_residues={12};unique_D_to_S=q12_to_q0;"
        f"q0_to_q3_events={u3_exits(0, 3)};"
        "uniform_over_169_addresses=1;"
        "relation=orthogonal_to_all_three_routes"
    )
    print(
        f"lawful_samples={len(samples)};unit_pairs="
        f"{tuple((n, n + 1) for n in atlas.UNIT_SECTIONS)};"
        "q0_to_q3_endpoint_weight_phase=1;Prony_nodes_identical=1"
    )
    for row in field_rows:
        print(
            f"field={row['prime']};raw_digest={row['raw_digest']};"
            f"C_first={row['C_first']};P={row['P']};"
            f"address_chi3_phase={row['address_phase']};"
            f"fibre_chi3_scalar={row['fibre_character']};"
            f"U0={row['U0']};V0={row['V0']};"
            f"Hermitian_support={row['hermitian_support']};"
            f"free_chi3_fibre_U0={row['free_chi3_fibre_U0']}"
        )
    print(
        "SELECTOR_DICHOTOMY=q3 right-endpoint selector pulled back through "
        "g gives the exact nonzero 26-sample P*C bank but changes the "
        "marked right-endpoint selector labels to ((12,5),(11,5)); the "
        "physical source atom and P are unchanged.  The q0 selector pushed "
        "forward through g lands on two identical present q3 charts and "
        "cancels at every endpoint frequency"
    )
    print(
        "FIBRE_DICHOTOMY=positive and address-chi3 target polarizations of "
        "F are nonzero endpoint-only Prony atlases if a free duplicate of P "
        "is attached; under the actual source-induced transport every one "
        "of the 9*26 terms in each field is zero because g_inverse(F) is "
        f"source-absent; pointwise_zero_checks="
        f"{strict_missing_fibre_term_checks}"
    )
    print(
        "COEFFICIENT_RECHART=no automorphism is needed: q0 and q3 raw "
        "banks and both Prony nodes are identical.  The address chi3 phase "
        "omega^2 may be removed by a common F-rational scalar gauge; this "
        "preserves U=chi3,V=trivial,the Kummer ratio,and the normalized "
        "Hermitian edge, but cannot repair either zero route"
    )
    print(
        "PHYSICAL_VERDICT=no lawful named-right-origin-preserving "
        "one-fibre completion within g plus F.  The nonzero pulled-selector "
        "atlas changes right-endpoint selector provenance while retaining "
        "the same physical source atom/P, then stops at q3's empty literal "
        "U ancestry; the missing fibre has no source coefficient termwise, "
        "and q0/q3 carry no QA/QAB sidecar.  The address completion neither "
        "contracts the E3/complement horn nor supplies a physical ancestry "
        "action"
    )
    print(
        "U3_CROSS_VERDICT=the directed exit is a genuine positive carry "
        "marker for wrapped residue paths, but q0->q3 is the nonwrapping "
        "h=3 path and has zero exits.  The address coordinate 12 in the "
        "pulled selector is not target residue q12; u3 is address-invariant, "
        "so the marker changes none of the pulled, pushed, or missing-fibre "
        "coefficient verdicts and supplies no absent U or QA/QAB label"
    )
    print(
        "scope=finite exact address/coefficient/physical typing theorem "
        "for g, F, and pointwise source-induced coefficients; broader "
        "signed correspondences or new ancestry sheets are not excluded; "
        "no q7 transporter, row exclusion, LRC14 decrement, or proof"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
