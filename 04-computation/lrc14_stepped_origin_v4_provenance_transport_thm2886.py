#!/usr/bin/env python3
"""Exact companion for the THM-2886 origin-odd provenance transport.

The finite question is whether the proved semantic arrow

    QA(q11,a) -> QAB(q7,a+1)

on the 449 THM-2835 sheets can carry an origin-odd E3/complement selector
while preserving the fixed source and the raw recombined current.

The same-block Boolean answer is negative, but the coupled carrier

    H=ker(e_E3 XOR semantic_slot7 XOR semantic_slot8)

selects a canonical stepped-only q3 subcopy.  Transporting that subcopy,
including its empty zero-origin component, gives an exact negative-oriented
q3/q11/q7 current and a positive answer at the finite endpoint/semantic
level.  The next implication fails: the event lift acts by diag(E,1) on
the Prony pair and therefore not scalarly on U+V.

The audit distinguishes three levels:

1. target-only Boolean coefficient transport;
2. a lawful diagonal source/target truth-block current;
3. the event-charged Prony channel versus its recombination U+V.

It does not promote the sheetwise H-product subcopy to the canonical current
interface, and proves no LRC row exclusion.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import product
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
    COMP / "lrc14_q11_semantic_word_horn_thm2835.py":
        "207dd94f086338ae1e80b7d7196f99bf41e795893d13b6d48e4e7d516af03523",
    RESULTS / "lrc14_q11_semantic_word_horn_thm2835.out":
        "1ebe0cbaf7d4ef13defed0bdb5b37df1364880acdbfc6139b243ab9df65f6bf6",
    COMP / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    COMP / "lrc14_two_origin_endpoint_projective_kummer_thm2868.py":
        "3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5",
    RESULTS / "lrc14_two_origin_endpoint_projective_kummer_thm2868.out":
        "ce4b8961b3dfa468d51b884b91002872440b05bf70d7e5eecfc3318d4100e2f9",
    COMP / "lrc14_endpoint_factor_exit_carry_transducer_thm2878.py":
        "b379b9278f6c0d0864908bbc2da2123f4d208eb83c35738d12f651119e7a3366",
    RESULTS / "lrc14_endpoint_factor_exit_carry_transducer_thm2878.out":
        "35bdec6bc5b2bb3c0287bd5aee26c66e8485876e066bf423e2fadb3a94727224",
    COMP / "lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.py":
        "3ed346e0c631b34bd61f0c4d27d7f161e8d35b70decfb95f5207c5f57893d005",
    RESULTS / "lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.out":
        "0faa0a24f6ba8b6c88b6bbfc4f225e38667097b1a937d977741453499884901c",
    COMP / "lrc14_macro_semantic_diagonal_horn_carrier_thm2884.py":
        "b739be20e741d5c061e0febcc8aef9b0f58f4ae8a648aa803610e0dad991929f",
    RESULTS / "lrc14_macro_semantic_diagonal_horn_carrier_thm2884.out":
        "8c3829b1052a641ca08a5e5bda86d9d5e8bd1584f5b2911c57c9fad9da41d4b6",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_q11_semantic_word_horn_thm2835 as word_horn
import lrc14_two_origin_endpoint_projective_kummer_thm2868 as atlas
import lrc14_endpoint_factor_exit_carry_transducer_thm2878 as transducer


allocation = atlas.allocation
endpoint_base = atlas.endpoint_base
endpoint = atlas.endpoint
P = 13
E3 = 1
NOT_E3 = 0
ORIGINS = (atlas.ORIGIN, atlas.STEPPED)


def complement(intervals, modulus):
    result = []
    cursor = 0
    for left, right in intervals:
        require(0 <= left < right <= modulus, "bad truth-block interval")
        require(cursor <= left, "truth-block intervals overlap")
        if cursor < left:
            result.append((cursor, left))
        cursor = right
    if cursor < modulus:
        result.append((cursor, modulus))
    return tuple(result)


def restrict(carrier, block):
    return allocation.indexed_weighted_intersection(
        carrier,
        block,
        tuple(left for left, _right in block),
    )


def marker_exits(patterns, q, h):
    return sum(
        patterns[(q + offset) % P][3] == "D"
        and patterns[(q + offset + 1) % P][3] == "S"
        for offset in range(h)
    )


def transform_pair(pair, permutation, functions):
    return tuple(
        functions[output_origin][pair[permutation[output_origin]]]
        for output_origin in range(2)
    )


def build_all_out(base, pattern, ell):
    """Reconstruct the all-out semantic word from the ambient circle."""
    require(
        set(pattern.values()) <= {"out", "gout"}
        and "gout" in pattern.values(),
        "all-out word received a positive or guardless pattern",
    )
    intervals = [(0, base.T_DEN)]
    for index, mode in pattern.items():
        if mode == "gout":
            intervals = base.subtract_comb(
                intervals,
                base.W[index],
                91,
                -13 - 7 * ell[index],
                13 - 7 * ell[index],
            )
    for _speed, index in sorted(
        (base.W[index], index)
        for index, mode in pattern.items()
        if mode == "out"
    ):
        intervals = base.subtract_comb(
            intervals,
            base.W[index],
            182,
            -13 - 14 * ell[index],
            13 - 14 * ell[index],
        )
    return tuple(intervals)


def chi(state):
    return state[0] ^ state[1] ^ state[2]


def main() -> None:
    (
        _module, full_module, _details, macro_e3, clocks, _q_pairs
    ) = allocation.build_geometry()
    period = full_module.T
    unit = period // P
    source_interval = allocation.ATOM_INTERVAL
    target0 = tuple(
        value + allocation.physical.SHIFT for value in source_interval
    )
    source_atom = ((*source_interval, 1),)
    target_atom = ((*target0, 1),)
    target_carriers = tuple(
        allocation.physical.overlap.shift_weighted(target_atom, q * unit)
        for q in range(P)
    )
    require(
        period == endpoint_base.T
        and tuple(macro_e3)
        == tuple(endpoint.build_set(
            endpoint_base.PAT_E3,
            endpoint_base.REPS[atlas.ORIGIN],
        ))
        and (12 * endpoint.RDIL * unit) % endpoint.NN == 0
        and (26 * endpoint.RDIL * unit) % endpoint.NN == 0,
        "shared geometry or target-phase normalization changed",
    )

    # Reconstruct both endpoint-origin truth partitions.  The fixed source
    # lies in E3 at both origins.  The target truth differs between origins
    # only at q3.
    truth_blocks = {}
    source_pieces = {}
    target_pieces = {}
    source_truth = []
    target_truth = [[], []]
    for origin_index, address in enumerate(ORIGINS):
        ell = endpoint_base.REPS[address]
        e3_block = tuple(endpoint.build_set(endpoint_base.PAT_E3, ell))
        not_e3_block = complement(e3_block, period)
        truth_blocks[origin_index] = {
            E3: e3_block,
            NOT_E3: not_e3_block,
        }
        source_candidates = {}
        for truth in (NOT_E3, E3):
            source_candidates[truth] = restrict(
                source_atom, truth_blocks[origin_index][truth]
            )
            source_pieces[origin_index, truth] = source_candidates[truth]
        containing_source = tuple(
            truth for truth, piece in source_candidates.items()
            if piece == source_atom
        )
        require(
            len(containing_source) == 1
            and not source_candidates[1 - containing_source[0]],
            "fixed source is not Boolean-total at an endpoint origin",
        )
        source_truth.append(containing_source[0])

        for q, carrier in enumerate(target_carriers):
            candidates = {}
            for truth in (NOT_E3, E3):
                candidates[truth] = restrict(
                    carrier, truth_blocks[origin_index][truth]
                )
                target_pieces[origin_index, q, truth] = candidates[truth]
            containing_target = tuple(
                truth for truth, piece in candidates.items()
                if piece == carrier
            )
            require(
                len(containing_target) == 1
                and not candidates[1 - containing_target[0]],
                "target atom is not Boolean-total at an endpoint origin",
            )
            target_truth[origin_index].append(containing_target[0])

    source_truth = tuple(source_truth)
    target_truth = tuple(tuple(row) for row in target_truth)
    truth_parity = tuple(
        target_truth[0][q] ^ target_truth[1][q] for q in range(P)
    )
    require(
        source_truth == (E3, E3)
        and tuple(q for q in range(P) if target_truth[0][q] == E3)
        == (0, 3, 11)
        and tuple(q for q in range(P) if target_truth[1][q] == E3)
        == (0, 11)
        and truth_parity == tuple(int(q == 3) for q in range(P)),
        "source/target truth atlas changed",
    )

    # Exhaust all four selector pairs at every residue.  A target-only
    # coefficient contributes at origin o iff s_o=t_o.  A lawful diagonal
    # current additionally requires the fixed source to lie in that same
    # selected truth block.
    selector_table = {}
    native_selector_support_union = set()
    for q in range(P):
        for selector in product((NOT_E3, E3), repeat=2):
            target_flags = tuple(
                int(selector[o] == target_truth[o][q]) for o in range(2)
            )
            diagonal_flags = tuple(
                int(
                    selector[o] == source_truth[o]
                    and selector[o] == target_truth[o][q]
                )
                for o in range(2)
            )
            target_signed = target_flags[0] - target_flags[1]
            diagonal_signed = diagonal_flags[0] - diagonal_flags[1]
            selector_table[q, selector] = (
                target_flags,
                diagonal_flags,
                target_signed,
                diagonal_signed,
            )
            if diagonal_signed:
                native_selector_support_union.add(q)
    require(
        native_selector_support_union == {0, 3, 11}
        and all(
            selector_table[7, selector][1] == (0, 0)
            for selector in product((NOT_E3, E3), repeat=2)
        ),
        "a Boolean selector produced a native diagonal q7 current",
    )

    # Unique selector retaining the full zero-origin target atom and the
    # empty stepped-origin block.  Its parity is 1 XOR delta_q3.
    positive_selector = tuple(
        (target_truth[0][q], 1 ^ target_truth[1][q])
        for q in range(P)
    )
    positive_parity = tuple(left ^ right for left, right in positive_selector)
    target_only_support = tuple(
        q for q in range(P)
        if selector_table[q, positive_selector[q]][2]
    )
    native_current_support = tuple(
        q for q in range(P)
        if selector_table[q, positive_selector[q]][3]
    )
    require(
        positive_parity == tuple(1 ^ int(q == 3) for q in range(P))
        and positive_selector[3] == (E3, E3)
        and positive_selector[11] == (E3, NOT_E3)
        and positive_selector[7] == (NOT_E3, E3)
        and target_only_support == tuple(range(P))
        and native_current_support == (0, 3, 11),
        "positive-selector coefficient/current split changed",
    )

    # Classify every componentwise Boolean truth operation, with or without
    # swapping the two marked origins.  The map must act on a block label
    # independently of whether the block is being used by the fixed source,
    # target, or selector.  This is the naturality content of "one
    # operation".  All four functions {0,1}->{0,1} are included.
    boolean_functions = tuple(product((NOT_E3, E3), repeat=2))
    permutations = ((0, 1), (1, 0))
    q11_selector = positive_selector[11]
    q7_selector = positive_selector[7]
    q11_truth = tuple(target_truth[o][11] for o in range(2))
    q7_truth = tuple(target_truth[o][7] for o in range(2))
    target_selector_candidates = []
    source_preserving_candidates = []
    bijective_target_selector_candidates = []
    for permutation in permutations:
        for functions in product(boolean_functions, repeat=2):
            transformed_selector = transform_pair(
                q11_selector, permutation, functions
            )
            transformed_target = transform_pair(
                q11_truth, permutation, functions
            )
            transformed_source = transform_pair(
                source_truth, permutation, functions
            )
            if (
                transformed_selector == q7_selector
                and transformed_target == q7_truth
            ):
                candidate = (
                    permutation,
                    functions,
                    transformed_source,
                )
                target_selector_candidates.append(candidate)
                if all(function in ((0, 1), (1, 0)) for function in functions):
                    bijective_target_selector_candidates.append(candidate)
                if transformed_source == source_truth:
                    source_preserving_candidates.append(candidate)
    require(
        len(target_selector_candidates) == 2
        and len(bijective_target_selector_candidates) == 1
        and bijective_target_selector_candidates[0]
        == ((0, 1), ((1, 0), (1, 0)), (NOT_E3, NOT_E3))
        and not source_preserving_candidates
        and all(
            candidate[2] == (NOT_E3, NOT_E3)
            for candidate in target_selector_candidates
        ),
        "one-operation Boolean rechart classification changed",
    )

    # Recover the 449 semantic sheets independently from the literal word
    # cylinders, including the previously omitted all-out Qempty word.
    old = word_horn.allocation.physical.relative.lift.m.core.old
    patterns = {
        "QA": dict(word_horn.atlas.PAT_QA),
        "QB": dict(word_horn.atlas.PAT_QB),
        "QAB": dict(word_horn.atlas.PAT_QAB),
    }
    qempty_pattern = {
        0: "gout",
        1: "out",
        2: "out",
        3: "out",
        4: "out",
        5: "out",
        6: "out",
        7: "out",
        8: "out",
    }
    require(
        (
            (patterns["QB"][7], patterns["QB"][8]),
            (patterns["QA"][7], patterns["QA"][8]),
            (patterns["QAB"][7], patterns["QAB"][8]),
        )
        == (("out", "in"), ("in", "out"), ("in", "in")),
        "semantic outer-word bits changed",
    )
    word_pieces = {
        name: tuple(old.base.build_set(pattern, old.base.ZELL))
        for name, pattern in patterns.items()
    }
    word_pieces["Qempty"] = build_all_out(
        old.base, qempty_pattern, old.base.ZELL
    )
    word_starts = {
        name: tuple(left for left, _right in pieces)
        for name, pieces in word_pieces.items()
    }
    source_qb = word_horn.word_flags(
        source_interval, 0, period, word_pieces["QB"], word_starts["QB"]
    )
    target_qa11 = word_horn.word_flags(
        source_interval,
        allocation.physical.SHIFT + 11 * unit,
        period,
        word_pieces["QA"],
        word_starts["QA"],
    )
    target_qab7 = word_horn.word_flags(
        source_interval,
        allocation.physical.SHIFT + 7 * unit,
        period,
        word_pieces["QAB"],
        word_starts["QAB"],
    )
    target_qempty3 = word_horn.word_flags(
        source_interval,
        allocation.physical.SHIFT + 3 * unit,
        period,
        word_pieces["Qempty"],
        word_starts["Qempty"],
    )
    labels = tuple(
        ancestry
        for ancestry in range(P**5)
        if source_qb[ancestry] and target_qa11[ancestry]
    )
    semantic_successors = tuple(
        ancestry
        for ancestry in labels
        if target_qab7[(ancestry + 1) % (P**5)]
    )
    q3_qempty_labels = tuple(
        ancestry for ancestry in labels if target_qempty3[ancestry]
    )
    q11_sheet_gate = selector_table[
        11, positive_selector[11]
    ][3]
    q7_sheet_gate = selector_table[
        7, positive_selector[7]
    ][3]
    q7_target_only_gate = selector_table[
        7, positive_selector[7]
    ][2]
    require(
        len(labels) == len(semantic_successors) == len(q3_qempty_labels) == 449
        and labels[0] == 59306
        and labels[-1] == 311961
        and {ancestry % 169 for ancestry in labels} == {156}
        and q11_sheet_gate == 1
        and q7_sheet_gate == 0
        and q7_target_only_gate == 1,
        "449-sheet selector/current gate changed",
    )

    # The coupled physical Boolean carrier
    #
    #       H=ker(e_E3 XOR semantic_slot7 XOR semantic_slot8)
    #
    # contains the fixed source QB, q11 QA, q7 QAB, and the stepped-origin
    # q3 Qempty state on every sheet.  A local H projector is nevertheless
    # origin-even at q11/q7.  Its genuinely useful datum is the canonical
    # stepped-only q3 image, which can be transported as a marked subcopy.
    outer_bits = {
        "Qempty": (0, 0),
        "QB": (0, 1),
        "QA": (1, 0),
        "QAB": (1, 1),
    }
    semantic_states = {
        word: (outer[0] ^ outer[1], *outer)
        for word, outer in outer_bits.items()
    }
    h_states = tuple(
        state for state in product((0, 1), repeat=3) if chi(state) == 0
    )
    require(
        h_states == ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0))
        and set(semantic_states.values()) == set(h_states)
        and semantic_states
        == {
            "Qempty": (0, 0, 0),
            "QB": (1, 0, 1),
            "QA": (1, 1, 0),
            "QAB": (0, 1, 1),
        },
        "coupled H carrier changed",
    )
    source_h_membership = tuple(
        int(chi((source_truth[o], *outer_bits["QB"])) == 0)
        for o in range(2)
    )
    checkpoint_words = {3: "Qempty", 11: "QA", 7: "QAB"}
    checkpoint_h_membership = {
        q: tuple(
            int(chi((target_truth[o][q], *outer_bits[word])) == 0)
            for o in range(2)
        )
        for q, word in checkpoint_words.items()
    }
    local_h_signed_amplitude = {
        q: pair[0] - pair[1]
        for q, pair in checkpoint_h_membership.items()
    }
    require(
        source_h_membership == (1, 1)
        and checkpoint_h_membership
        == {3: (0, 1), 11: (1, 1), 7: (1, 1)}
        and local_h_signed_amplitude == {3: -1, 11: 0, 7: 0},
        "local H origin-incidence boundary changed",
    )

    # Start from H's canonical stepped-only q3 piece and require target
    # translation to carry empty to empty and full to full at each marked
    # origin.  This uniquely forces the negative-orientation selector
    #
    #       q3:(notE3,notE3)
    #       q11:(notE3,E3)
    #       q7:(E3,notE3).
    #
    # It has the same required parity 1 XOR delta_q3 as the positive
    # selector, and q11->q7 is one global truth complement.
    negative_selector = tuple(
        (1 ^ target_truth[0][q], target_truth[1][q])
        for q in range(P)
    )
    transported_selector = {
        3: negative_selector[3],
        11: negative_selector[11],
        7: negative_selector[7],
    }
    transported_pieces = {
        (q, origin): target_pieces[
            origin, q, transported_selector[q][origin]
        ]
        for q in (3, 11, 7)
        for origin in range(2)
    }
    require(
        transported_selector
        == {3: (NOT_E3, NOT_E3), 11: (NOT_E3, E3), 7: (E3, NOT_E3)}
        and tuple(
            transported_selector[q][0] ^ transported_selector[q][1]
            for q in (3, 11, 7)
        ) == (0, 1, 1)
        and all(
            allocation.physical.overlap.shift_weighted(
                transported_pieces[3, origin], 8 * unit
            ) == transported_pieces[11, origin]
            for origin in range(2)
        )
        and all(
            allocation.physical.overlap.shift_weighted(
                transported_pieces[11, origin], 9 * unit
            ) == transported_pieces[7, origin]
            for origin in range(2)
        )
        and tuple(bool(transported_pieces[3, origin]) for origin in range(2))
        == (False, True)
        and tuple(bool(transported_pieces[11, origin]) for origin in range(2))
        == (False, True)
        and tuple(bool(transported_pieces[7, origin]) for origin in range(2))
        == (False, True),
        "H-seeded origin-provenance transport changed",
    )

    # Verify the target-cylinder arithmetic before quotienting by the
    # physical period.  The wrapped q11->q7 edge is exactly the ancestry
    # increment a->a+1 on all 449 labels.
    target_line = tuple(
        value + allocation.physical.SHIFT for value in source_interval
    )
    target_line3 = tuple(value + 3 * unit for value in target_line)
    target_line11 = tuple(value + 11 * unit for value in target_line)
    target_line7 = tuple(value + 7 * unit for value in target_line)
    line_transport_checks = 0
    for ancestry in labels:
        cylinder3 = tuple(
            value + ancestry * period for value in target_line3
        )
        cylinder11 = tuple(
            value + ancestry * period for value in target_line11
        )
        cylinder7_successor = tuple(
            value + (ancestry + 1) * period for value in target_line7
        )
        require(
            tuple(value + 8 * unit for value in cylinder3) == cylinder11
            and tuple(value + 9 * unit for value in cylinder11)
            == cylinder7_successor,
            "literal target-cylinder transport changed",
        )
        line_transport_checks += 2

    # The fixed horn cell loses only macro E3 at q7; every ordinary factor
    # survives.  Thus the coupled H candidate pays exactly one cross-grade
    # invoice rather than hiding an additional q1/q2 defect.
    horn_cell = (0, 5, 1)
    target3_interval = (
        (target0[0] + 3 * unit) % period,
        (target0[1] + 3 * unit) % period,
    )
    target11_interval = (
        (target0[0] + 11 * unit) % period,
        (target0[1] + 11 * unit) % period,
    )
    target7_interval = (
        (target0[0] + 7 * unit) % period,
        (target0[1] + 7 * unit) % period,
    )
    source_signature = atlas.signature(
        source_interval, *horn_cell, full_module, macro_e3, clocks,
    )
    target3_signature = atlas.signature(
        target3_interval, *horn_cell, full_module, macro_e3, clocks,
    )
    target11_signature = atlas.signature(
        target11_interval, *horn_cell, full_module, macro_e3, clocks,
    )
    target7_signature = atlas.signature(
        target7_interval, *horn_cell, full_module, macro_e3, clocks,
    )
    require(
        all(source_signature)
        and all(target3_signature)
        and all(target11_signature)
        and target7_signature == (False, True, True, True, True, True),
        "fixed horn cell acquired an extra q7 factor loss",
    )

    # Carry and selector parity are independent.  In particular the 449
    # q11->q7 arrows carry once but preserve odd selector parity.
    zero_ell = endpoint_base.REPS[atlas.ORIGIN]
    literal_patterns = tuple(
        transducer.literal_pattern(transducer.TARGET_ATOM, q, zero_ell)
        for q in range(P)
    )
    selector_toggle_edges = tuple(
        q for q in range(P)
        if positive_parity[q] != positive_parity[(q + 1) % P]
    )
    require(
        marker_exits(literal_patterns, 11, 9) == 1
        and positive_parity[11] == positive_parity[7] == 1
        and marker_exits(literal_patterns, 12, 1) == 1
        and positive_parity[12] == positive_parity[0] == 1
        and marker_exits(literal_patterns, 0, 3) == 0
        and positive_parity[0] != positive_parity[3]
        and selector_toggle_edges == (2, 3),
        "selector/carry independence changed",
    )

    # Exact finite-field recombined-current audit over all 26 samples.
    # "Cross-graded" means that the fixed E3 source scalar is reused even
    # when the selected target lies in not-E3.  "Native diagonal" requires
    # source and target to occupy the same selected truth block.
    all_samples = tuple(
        sample
        for actual in atlas.UNIT_SECTIONS
        for sample in (actual, actual + 1)
    )
    field_rows = []
    recombined_current_checks = 0
    for prime, root in endpoint.MODS:
        xi = pow(root, endpoint.NN // 2366, prime)
        omega = pow(xi, 182, prime)
        source_value = atlas.COMMON_SOURCE[prime]
        require(
            source_value != 0
            and pow(omega, P, prime) == 1
            and omega != 1,
            "endpoint-field normalization changed",
        )
        cross_support_rows = []
        diagonal_support_rows = []
        q11_cross_values = []
        q7_cross_values = []
        q11_diagonal_values = []
        q7_diagonal_values = []
        transported_raw = {3: [], 11: [], 7: []}
        transported_sheet_sums = {3: [], 11: [], 7: []}
        for sample in all_samples:
            cross_values = []
            diagonal_values = []
            for q in range(P):
                selector = positive_selector[q]
                target_values = tuple(
                    allocation.endpoint_sum(
                        target_pieces[o, q, selector[o]],
                        -(12 + 26 * sample),
                        (prime, root),
                    )
                    for o in range(2)
                )
                cross_current = (
                    source_value * (target_values[0] - target_values[1])
                ) % prime
                diagonal_terms = tuple(
                    (
                        source_value
                        if source_pieces[o, selector[o]] == source_atom
                        else 0
                    )
                    * target_values[o] % prime
                    for o in range(2)
                )
                diagonal_current = (
                    diagonal_terms[0] - diagonal_terms[1]
                ) % prime
                cross_values.append(cross_current)
                diagonal_values.append(diagonal_current)
                recombined_current_checks += 1
                if q == 11:
                    q11_cross_values.append(cross_current)
                    q11_diagonal_values.append(diagonal_current)
                if q == 7:
                    q7_cross_values.append(cross_current)
                    q7_diagonal_values.append(diagonal_current)
            for q in (3, 11, 7):
                selector = transported_selector[q]
                origin_values = tuple(
                    allocation.endpoint_sum(
                        transported_pieces[q, origin],
                        -(12 + 26 * sample),
                        (prime, root),
                    )
                    for origin in range(2)
                )
                require(
                    origin_values[0] == 0
                    and origin_values[1] != 0,
                    "transported selector lost zero-minus-stepped typing",
                )
                current = (
                    source_value * (origin_values[0] - origin_values[1])
                ) % prime
                transported_raw[q].append(current)
                transported_sheet_sums[q].append(449 * current % prime)
            cross_support_rows.append(tuple(
                q for q, value in enumerate(cross_values) if value
            ))
            diagonal_support_rows.append(tuple(
                q for q, value in enumerate(diagonal_values) if value
            ))
        require(
            tuple(cross_support_rows) == (tuple(range(P)),) * len(all_samples)
            and tuple(diagonal_support_rows)
            == ((0, 3, 11),) * len(all_samples)
            and all(
                q11 == q7 != 0
                for q11, q7 in zip(q11_cross_values, q7_cross_values)
            )
            and all(value != 0 for value in q11_diagonal_values)
            and all(value == 0 for value in q7_diagonal_values),
            "cross-graded/native recombined-current boundary changed",
        )
        require(
            all(
                transported_raw[3][index]
                == transported_raw[11][index]
                == transported_raw[7][index]
                != 0
                for index in range(len(all_samples))
            )
            and all(
                transported_sheet_sums[3][index]
                == transported_sheet_sums[11][index]
                == transported_sheet_sums[7][index]
                != 0
                for index in range(len(all_samples))
            )
            and all(
                transported_raw[11][index]
                == -q11_cross_values[index] % prime
                and transported_raw[7][index]
                == -q7_cross_values[index] % prime
                for index in range(len(all_samples))
            ),
            "H-seeded transported raw current changed",
        )

        # Downstream hostile after freely adjoining the missing cross-grade:
        # the event lift acts on (U,V) as diag(E,1).  On the q11->q7 carry
        # leg E=omega^4, so U transports but U+V does not transform by E.
        q0_piece = target_pieces[0, 0, positive_selector[0][0]]
        left, right, weight = q0_piece[0]
        require(weight == 1, "truth-polarized atom weight changed")
        alpha_left = pow(
            root, 12 * endpoint.RDIL * left % endpoint.NN, prime
        )
        alpha_right = pow(
            root, 12 * endpoint.RDIL * right % endpoint.NN, prime
        )
        lambda_left = pow(
            root, 26 * endpoint.RDIL * left % endpoint.NN, prime
        )
        lambda_right = pow(
            root, 26 * endpoint.RDIL * right % endpoint.NN, prime
        )
        inverse_difference = pow(
            (lambda_left - lambda_right) % prime, -1, prime
        )
        split_u = []
        split_v = []
        for formal, offset, actual in zip(
            atlas.FREQUENCY_SECTIONS,
            atlas.SECTION_OFFSETS,
            atlas.UNIT_SECTIONS,
        ):
            raw = []
            for sample in (actual, actual + 1):
                raw.append(
                    source_value
                    * allocation.endpoint_sum(
                        q0_piece,
                        -(12 + 26 * sample),
                        (prime, root),
                    )
                    % prime
                )
            current, current_next = raw
            left_at_actual = (
                current_next - lambda_right * current
            ) * inverse_difference % prime
            right_at_actual = (
                lambda_left * current - current_next
            ) * inverse_difference % prime
            split_u.append(
                left_at_actual
                * pow(pow(lambda_left, offset, prime), -1, prime)
                % prime
            )
            split_v.append(
                right_at_actual
                * pow(pow(lambda_right, offset, prime), -1, prime)
                % prime
            )
            require(
                split_u[-1]
                == source_value * alpha_left
                * pow(lambda_left, formal, prime) % prime
                and split_v[-1]
                == -source_value * alpha_right
                * pow(lambda_right, formal, prime) % prime,
                "Prony split formula changed",
            )
        require(
            all(value != 0 for value in split_u)
            and all(value != 0 for value in split_v)
            and all(
                split_u[(r + 1) % P]
                == pow(omega, 3, prime) * split_u[r] % prime
                for r in range(P)
            )
            and len(set(split_v)) == 1,
            "charged/trivial Prony channels changed",
        )
        edge_scalar = pow(omega, 4, prime)
        scalar_recombination_failures = 0
        for ancestry_residue in range(P):
            source_section = (ancestry_residue + 11 - 3) % P
            target_section = (ancestry_residue + 1 + 7 - 3) % P
            require(
                split_u[target_section]
                == edge_scalar * split_u[source_section] % prime
                and split_v[target_section] == split_v[source_section],
                "q11->q7 split-channel transport changed",
            )
            source_recombined = (
                split_u[source_section] + split_v[source_section]
            ) % prime
            target_recombined = (
                split_u[target_section] + split_v[target_section]
            ) % prime
            failure = (
                target_recombined
                - edge_scalar * source_recombined
            ) % prime
            require(
                failure
                == (1 - edge_scalar) * split_v[source_section] % prime
                and failure != 0,
                "recombined scalar hostile disappeared",
            )
            scalar_recombination_failures += 1
        field_rows.append((
            prime,
            source_value,
            q11_cross_values[0],
            q7_cross_values[0],
            q11_diagonal_values[0],
            q7_diagonal_values[0],
            transported_raw[3][0],
            transported_raw[11][0],
            transported_raw[7][0],
            transported_sheet_sums[7][0],
            edge_scalar,
            split_u[8],
            split_v[8],
            scalar_recombination_failures,
        ))

    require(
        len(all_samples) == 26
        and recombined_current_checks == 2 * 26 * P
        and all(row[-1] == P for row in field_rows),
        "finite-field audit census changed",
    )

    print("LRC14 ORIGIN-ODD SELECTOR / 449-SHEET PHYSICALITY AUDIT")
    print(
        "status=PROVED+VERIFIED-EXACT_CANDIDATE; "
        "no row exclusion or LRC14 conclusion"
    )
    print(
        "universe=13 residues x 2 marked origins x 2 E3-truth blocks x "
        "4 selector pairs; 32 object-independent componentwise Boolean "
        "recharts (including origin swap); 449 whole semantic sheets; "
        "2 certified endpoint fields x 26 raw samples; "
        "2 fields x 13 split-channel carry states"
    )
    print(
        f"fixed_source_truth={source_truth}; "
        f"target_E3_supports=(origin0={tuple(q for q in range(P) if target_truth[0][q])},"
        f"origin1={tuple(q for q in range(P) if target_truth[1][q])}); "
        f"truth_parity={truth_parity}=delta_q3"
    )
    print(
        f"positive_selector_bits={positive_selector}; "
        f"selector_parity={positive_parity}=1_XOR_delta_q3; "
        f"toggle_edges={selector_toggle_edges}"
    )
    print(
        "seam_selector="
        f"(q11={positive_selector[11]},q7={positive_selector[7]}); "
        "target_truth=(q11=(1,1),q7=(0,0)); "
        "the_unique_bijective_target/selector_transport_is_global_"
        "E3_complement"
    )
    print(
        f"all_object-independent_target/selector_recharts="
        f"{target_selector_candidates}; "
        f"source_preserving_recharts={source_preserving_candidates}; "
        "same-block_baseline_no_go=every_target/selector_rechart_sends_"
        "source_(1,1)_to_(0,0)"
    )
    print(
        f"selector_target_only_support={target_only_support}; "
        f"lawful_same-block_current_support={native_current_support}; "
        f"union_over_all_selectors={tuple(sorted(native_selector_support_union))}; "
        "q7_same-block_origin_flags=(0,0)_for_all_four_selectors"
    )
    print(
        f"H=ker(e_E3_XOR_u7_XOR_v8)={h_states}; "
        f"states={semantic_states}; fixed_source_QB_membership="
        f"{source_h_membership}"
    )
    print(
        f"local_H_membership_by_origin={checkpoint_h_membership}; "
        f"local_H_signed_zero-minus-stepped={local_h_signed_amplitude}; "
        "local_H_restriction_is_exactly_q3_(zero_empty,stepped_full)_"
        "but_re-evaluation_is_origin-even_and_cancels_at_q11/q7"
    )
    print(
        f"H_seeded_transported_selector={transported_selector}; "
        "ordered_pairs_are_(zero_origin,stepped_origin); "
        "occupancy=(0,1)_at_q3,q11,q7; signed_orientation=zero-minus-"
        "stepped=-1; parity=(0,1,1)=1_XOR_delta_q3_on_the_seam"
    )
    print(
        "typed_objects=source_(I,a,QB,E3,fixed)_to_"
        "targets_(J3,a,Qempty,stepped-notE3),"
        "(J11,a,QA,stepped-E3),"
        "(J7,a+1,QAB,stepped-notE3); "
        "zero-origin_component_is_transported_empty"
    )
    print(
        f"endpoint_piece_maps=q3_plus8U_to_q11_and_q11_plus9U_to_q7_"
        "hold_at_both_origins; literal_cylinder_checks="
        f"{line_transport_checks}; q11_plus9U_wrap_is_a_to_a+1"
    )
    print(
        f"fixed_horn_cell={horn_cell}; factor_order=(E3,clock,q1,q2,c2,c3); "
        f"signatures=(source={source_signature},q3={target3_signature},"
        f"q11={target11_signature},q7={target7_signature}); "
        "q7_loses_only_E3_and_preserves_the_other_five_factors"
    )
    print(
        f"semantic_sheets=(count={len(labels)},first={labels[0]},"
        f"last={labels[-1]},source_residue_mod169=156,"
        "successor_residue_mod169=157); "
        "Qempty(q3,a)=449/449; QA(q11,a)->QAB(q7,a+1)=449/449"
    )
    print(
        f"same-block_sheetwise_gate=(q11={q11_sheet_gate}*449,"
        f"q7={q7_sheet_gate}*449); "
        f"q7_target_only_coefficient_gate={q7_target_only_gate}*449; "
        "minimal_witness=a=59306"
    )
    print(
        "selector_carry_independence=q11_to_q7_has_carry1_but_parity1_to1; "
        "q12_to_q0_has_carry1_but_parity1_to1; "
        "q0_to_q3_has_carry0_but_parity1_to0"
    )
    for row in field_rows:
        print(
            f"field={row[0]}; fixed_source={row[1]}; "
            f"crossgraded_raw_(q11,q7)=({row[2]},{row[3]}); "
            f"native_diagonal_raw_(q11,q7)=({row[4]},{row[5]}); "
            f"H_seeded_raw_(q3,q11,q7)=({row[6]},{row[7]},{row[8]}); "
            f"449_sheet_sum_q7={row[9]}; "
            f"charged_edge_E=omega^4={row[10]}; "
            f"sample_(U8,V8)=({row[11]},{row[12]}); "
            f"recombined_scalar_failures={row[13]}/13"
        )
    print(
        "recombined_current_checks=676; crossgraded_support=all13_on_all_"
        "52_field/sample_rows; native_diagonal_support={0,3,11}; "
        "q11_crossgraded_equals_q7_crossgraded_nonzero_but_q7_native=0"
    )
    print(
        "H_seeded_origin-resolved_raw_current=q3=q11=q7_nonzero_on_all_"
        "52_field/sample_rows; 449-fold_sums_nonzero; this_closes_the_"
        "selector_parity_and_fixed-source_coefficient_transport_on_the_"
        "typed_sheetwise_horn"
    )
    print(
        "downstream_hostile_after_freely_adjoining_crossgrade="
        "(U,V)->(E*U,V), E=omega^4 on q11->q7; "
        "S'=E*U+V differs from E*S by (1-E)*V !=0 on 26/26 "
        "field/state rows"
    )
    print(
        "POSITIVE_RESULT=the_H-selected_stepped-only_q3_subcopy_has_a_"
        "unique_empty/full_preserving_originwise_transport_to_q11_and_q7; "
        "the_single_q11_QA(a)_to_q7_QAB(a+1)_operation_is_source-fixed,"
        "origin-odd,truth-complementing,and_raw-current-preserving"
    )
    print(
        "FIRST_FAILED_IMPLICATION=the_proved_raw_recombined_current_does_"
        "not_become_a_scalar_charged_carry_current: diag(E,1)_on_(U,V)_"
        "does_not_descend_through_U+V; additionally_canon_still_needs_the_"
        "sheetwise_H/product_subcopy_registered_as_a_full_current_functor"
    )
    print(
        "VERDICT=YES_at_finite_exact_origin-resolved_endpoint/semantic_"
        "raw-current_level; NO_scalar_Bockstein_transport_on_recombined_"
        "U+V_and_no_LRC14_conclusion"
    )
    print(
        "NEXT_TEST=formalize_the_H-graded_transported_subcopy_in_the_"
        "canonical_current_interface_then_seek_a_lawful_two-channel_or_"
        "Hermitian_observable_instead_of_a_scalar_action_on_U+V"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
